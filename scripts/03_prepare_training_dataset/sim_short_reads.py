#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Short-read community simulator (ISS)
====================================================

This script generates paired-end Illumina-like reads with InSilicoSeq (iss):
- Builds per-class merged FASTAs with normalized headers: >Class|filename|orig_id
- Filters sequences to length >= 500 bp (faster simulation)
- Derives per-record weights from per-genome class TSVs (generic scenario)
- Runs ISS in CHUNKS with bounded retries to reach exact requested pairs
- Concatenates chunk pieces into final class outputs
- Merges class outputs into 'generic' and 'filtered' mixes (fixed % per class)
- Verifies pair counts; optionally shuffles pairs reproducibly (streaming, no RAM blow-up)

Outputs are organized under:
  <COMM>/reads/short/<split>/
    generic_R1.fastq.gz  generic_R2.fastq.gz
    filtered_R1.fastq.gz filtered_R2.fastq.gz
    (optional) generic_shuf_R1.fastq.gz generic_shuf_R2.fastq.gz

Assumes per-class selection files exist:
  <COMM>/genomes/{prokaryote|microeukaryote|virus|plasmid}.sel

Assumes per-class *generic* abundance TSVs exist:
  <COMM>/abundance/{prokaryote|microeukaryote|virus|plasmid}_generic.tsv
Each line: <genome_path>\t<weight>

Requires coreutils `paste`, `shuf` if --shuffle is used (same as the Bash version).
"""

from __future__ import annotations
import argparse
import gzip
import math
import os
import re
import shutil
import subprocess as sp
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

# ----------------------------- logging & utils -----------------------------

def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def need_exe(path_or_name: str, desc: str) -> None:
    exe = shutil.which(path_or_name) if not os.access(path_or_name, os.X_OK) else path_or_name
    if not exe:
        raise RuntimeError(f"Missing tool: {desc} ({path_or_name})")
    # no return; just validate existence

# Count FASTQ pairs in a (possibly gz) file
def count_pairs_fastq(path: Path) -> int:
    opener = gzip.open if str(path).endswith('.gz') else open
    n_lines = 0
    with opener(path, 'rt', errors='ignore') as f:
        for _ in f:
            n_lines += 1
    return n_lines // 4

def atomic_move(src: Path, dst: Path) -> None:
    tmp = dst.with_suffix(dst.suffix + ".part")
    shutil.copy2(src, tmp)
    tmp.replace(dst)

def gz_concat(out_gz: Path, files: List[Path]) -> None:
    """Repack concat into a fresh gzip (like bash gz_repack_concat)."""
    ensure_dir(out_gz.parent)
    tmp = out_gz.with_suffix(out_gz.suffix + ".part")
    with gzip.open(tmp, 'wb', compresslevel=1, mtime=0) as fo:
        for f in files:
            if str(f).endswith('.gz'):
                with gzip.open(f, 'rb') as fi:
                    shutil.copyfileobj(fi, fo)
            else:
                with open(f, 'rb') as fi:
                    shutil.copyfileobj(fi, fo)
    # quick integrity test
    with gzip.open(tmp, 'rb') as fi:
        for _ in fi:
            break
    tmp.replace(out_gz)

# ------------------------------ domain logic -------------------------------

CLASS_LC_TO_CAP = {
    "prokaryote": "Prokaryote",
    "microeukaryote": "Microeukaryote",
    "virus": "Virus",
    "plasmid": "Plasmid",
}
CLASSES = ["prokaryote", "microeukaryote", "virus", "plasmid"]

def build_class_fasta(comm: Path, cls_lc: str, refs_dir: Path) -> Path:
    """Merge FASTAs listed in genomes/<cls>.sel with normalized headers."""
    cls_cap = CLASS_LC_TO_CAP[cls_lc]
    sel = comm / "genomes" / f"{cls_lc}.sel"
    out_fa = refs_dir / f"{cls_cap}.fasta"
    if not sel.exists() or sel.stat().st_size == 0:
        raise RuntimeError(f"Missing selection list: {sel}")
    ensure_dir(refs_dir)
    with out_fa.open('wt') as out:
        with sel.open('rt') as fsel:
            for p in (Path(line.strip()) for line in fsel if line.strip()):
                if not p.exists() or p.stat().st_size == 0:
                    log(f"[WARN] empty/missing {p}")
                    continue
                base = p.name
                opener = gzip.open if str(p).endswith('.gz') else open
                with opener(p, 'rt', errors='ignore') as fi:
                    for line in fi:
                        if line.startswith('>'):
                            # take first token after '>'
                            tok = line[1:].strip().split()[0]
                            out.write(f">{cls_cap}|{base}|{tok}\n")
                        else:
                            out.write(line.strip() + "\n")
    if out_fa.stat().st_size == 0:
        raise RuntimeError(f"Built empty FASTA: {out_fa}")
    return out_fa

def filter_fasta_minlen(in_fa: Path, min_len: int = 500) -> None:
    """In-place filter records shorter than min_len."""
    tmp = in_fa.with_suffix(in_fa.suffix + ".filtered")
    with in_fa.open('rt', errors='ignore') as fi, tmp.open('wt') as fo:
        hdr = None
        buf = []
        def flush():
            if hdr and sum(len(s) for s in buf) >= min_len:
                fo.write(hdr + "\n")
                for s in buf:
                    fo.write(s + "\n")
        for line in fi:
            if line.startswith('>'):
                flush()
                hdr = line.strip()
                buf = []
            else:
                buf.append(line.strip())
        flush()
    tmp.replace(in_fa)

def read_pergen_tsv(tsv: Path) -> Dict[str, float]:
    """Read <genome_path>\t<weight> → map[basename]=weight_sum."""
    m: Dict[str, float] = {}
    with tsv.open('rt') as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            toks = re.split(r'\s+', ln)
            if len(toks) < 2:
                continue
            base = Path(toks[0]).name
            w = float(toks[1])
            m[base] = m.get(base, 0.0) + w
    return m

def write_perrecord_weights(refs_fa: Path, pergen_tsv: Path, out_tsv: Path) -> None:
    """Distribute genome-level weights to each record proportional to record length.
       Record IDs must match the header produced in build_class_fasta (Class|filename|orig_id)."""
    g_weight = read_pergen_tsv(pergen_tsv)
    # accumulate lengths per record (group by filename in header)
    lens: Dict[str, List[Tuple[str, int]]] = {}
    with refs_fa.open('rt', errors='ignore') as f:
        hdr = None
        fname = None
        L = 0
        for line in f:
            if line.startswith('>'):
                if hdr is not None and fname is not None:
                    lens.setdefault(fname, []).append((hdr, L))
                hdr = line[1:].strip().split()[0]   # Class|fname|orig
                parts = hdr.split('|')
                fname = parts[1] if len(parts) >= 2 else Path(hdr).name
                L = 0
            else:
                L += len(line.strip())
        if hdr is not None and fname is not None:
            lens.setdefault(fname, []).append((hdr, L))
    # weight to each record
    records: List[Tuple[str, float]] = []
    for fname, recs in lens.items():
        gw = g_weight.get(fname, 0.0)
        tot = sum(L for _, L in recs) or 1.0
        for hdr, L in recs:
            records.append((hdr, gw * (L / tot)))
    S = sum(w for _, w in records)
    ensure_dir(out_tsv.parent)
    with out_tsv.open('wt') as out:
        if S <= 0.0:
            n = len(records) or 1
            for hdr, _ in records:
                out.write(f"{hdr}\t{1.0/n:.10g}\n")
        else:
            inv = 1.0 / S
            for hdr, w in records:
                out.write(f"{hdr}\t{(w*inv):.10g}\n")

def pick_mate(prefix: Path, mate: int) -> Optional[Path]:
    """Find ISS output for a given prefix and mate."""
    candidates = [
        f"{prefix}_R{mate}.fastq.gz",
        f"{prefix}_{mate}.fastq.gz",
        f"{prefix}_R{mate}.fastq",
        f"{prefix}_{mate}.fastq",
    ]
    for c in candidates:
        p = Path(c)
        if p.exists() and p.stat().st_size > 0:
            return p
    return None

@dataclass
class IssChunkSettings:
    iss_bin: Path
    model: str
    base_seed: int
    chunks: int
    chunk_max_retries: int
    cpus_per_iss: int = 1

def run_iss_generate(settings: IssChunkSettings, refs_fa: Path, abund_rec: Path,
                     need_pairs: int, tmp_prefix: Path, seed_offset: int) -> Tuple[Path, Path, int]:
    """Run one ISS generate call for 'need_pairs' pairs, returns (R1, R2, produced_pairs)."""
    n_reads = need_pairs * 2  # pairs → reads
    seed_i = settings.base_seed + seed_offset
    cmd = [
        str(settings.iss_bin), "generate",
        "--genomes", str(refs_fa),
        "--abundance_file", str(abund_rec),
        "--model", settings.model,
        "--seed", str(seed_i),
        "--n_reads", str(n_reads),
        "--cpus", str(settings.cpus_per_iss),
        "--output", str(tmp_prefix),
        "--compress",
    ]
    log(f"[ISS] {' '.join(cmd)}")
    rc = sp.call(cmd)
    if rc != 0:
        raise RuntimeError(f"iss generate failed rc={rc}")

    r1 = pick_mate(tmp_prefix, 1)
    r2 = pick_mate(tmp_prefix, 2)
    if not r1 or not r2:
        raise RuntimeError("Missing ISS outputs (R1/R2)")
    p1 = count_pairs_fastq(r1)
    p2 = count_pairs_fastq(r2)
    if p1 != p2:
        raise RuntimeError(f"Mate mismatch: R1={p1}, R2={p2}")
    return r1, r2, p1

def gen_group_chunked(group_label: str, refs_fa: Path, abund_rec: Path,
                      n_pairs: int, prefbase: Path, tmpd: Path,
                      settings: IssChunkSettings) -> Tuple[Path, Path]:
    """Produce exactly n_pairs for a class by chunking & retrying. Returns (final_R1, final_R2)."""
    final_r1 = prefbase.with_name(prefbase.name + f"_{group_label}_R1.fastq.gz")
    final_r2 = prefbase.with_name(prefbase.name + f"_{group_label}_R2.fastq.gz")
    for p in (final_r1, final_r2):
        if p.exists():
            p.unlink()

    per = n_pairs // settings.chunks
    rem = n_pairs - per * settings.chunks
    log(f"[ISS:{group_label}] pairs={n_pairs} chunks={settings.chunks} per={per} (+1 first {rem})")

    pieces_r1: List[Path] = []
    pieces_r2: List[Path] = []

    for i in range(settings.chunks):
        want = per + (1 if i < rem else 0)
        if want <= 0:
            continue
        produced = 0
        attempts = 0
        while produced < want:
            attempts += 1
            if attempts > settings.chunk_max_retries:
                raise RuntimeError(f"{group_label} c{i}: retries exceeded")
            need = want - produced
            pref = tmpd / f"{group_label}.c{i}.t{attempts}"
            # clean any leftovers
            for sfx in ["_R1.fastq.gz", "_R2.fastq.gz", "_1.fastq.gz", "_2.fastq.gz",
                        "_R1.fastq", "_R2.fastq", "_1.fastq", "_2.fastq"]:
                try:
                    (pref.parent / (pref.name + sfx)).unlink()
                except FileNotFoundError:
                    pass
            r1, r2, got = run_iss_generate(settings, refs_fa, abund_rec, need, pref, seed_offset=(i*100 + attempts))
            produced += got
            pieces_r1.append(r1)
            pieces_r2.append(r2)

    gz_concat(final_r1, pieces_r1)
    gz_concat(final_r2, pieces_r2)
    # quick integrity check
    _ = count_pairs_fastq(final_r1)
    _ = count_pairs_fastq(final_r2)
    return final_r1, final_r2

def assert_pairs(label: str, fq1: Path, fq2: Path, want_pairs: int) -> None:
    got1 = count_pairs_fastq(fq1)
    got2 = count_pairs_fastq(fq2)
    if got1 != want_pairs or got2 != want_pairs:
        raise RuntimeError(f"{label}: wanted {want_pairs} pairs, got R1={got1} R2={got2}")

def shuffle_pairs_stream(seed: int, in_r1: Path, in_r2: Path, out_r1: Path, out_r2: Path) -> None:
    """Replicate the bash streaming shuffle (paste+shuf) without loading whole files in RAM.
       Requires `paste` and `shuf` in PATH."""
    need_exe("paste", "paste")
    need_exe("shuf", "shuf")
    # Build the pipeline equivalent to the bash helper.
    # paste - - - - groups 4 lines into a single tab-separated record.
    # We do two pastes and then paste them together column-wise.
    # Using /bin/sh -c keeps it readable; provide seed via --random-source.
    cmd = [
        "/bin/sh", "-c",
        (
            "paste - - - - <(gzip -dc \"$0\") | "
            "paste - - - - <(gzip -dc \"$1\") | "
            "shuf --random-source=<(yes \"$2\") | "
            "tee >(cut -f1-4 | tr '\\t' '\\n' | gzip -c > \"$3\") | "
            "cut -f5-8 | tr '\\t' '\\n' | gzip -c > \"$4\""
        )
    ]
    env = os.environ.copy()
    # Use bash for process substitution
    proc = sp.Popen(
        ["/bin/bash", "-c",
         f"paste <(gzip -dc {shlex(in_r1)}) - - - - <(gzip -dc {shlex(in_r2)}) - - - - | "
         f"shuf --random-source=<(yes {seed}) | "
         f"tee >(cut -f1-4 | tr '\\t' '\\n' | gzip -c > {shlex(out_r1)}) | "
         f"cut -f5-8 | tr '\\t' '\\n' | gzip -c > {shlex(out_r2)}"],
        env=env, stdout=sp.PIPE, stderr=sp.STDOUT)
    out, _ = proc.communicate()
    if proc.returncode != 0:
        sys.stdout.buffer.write(out or b"")
        raise RuntimeError("shuffle pipeline failed")

def shlex(p: Path) -> str:
    return "'" + str(p).replace("'", "'\"'\"'") + "'"

# ------------------------------ main routine -------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--root', type=Path, default=Path("/path/to/training"),
                    help="Training root that contains train_communities/*")
    ap.add_argument('--comm', type=Path, help="Full path to a community directory")
    ap.add_argument('--comm-name', type=str, help="Community name under <root>/train_communities (e.g., TrainC_short_01)")
    ap.add_argument('--split', choices=['train', 'val', 'test'], default='train',
                    help='Dataset split subfolder for outputs')
    ap.add_argument('--iss-bin', type=Path, default=Path("/path/to/conda_envs/iss-py310/bin/iss"))
    ap.add_argument('--model', type=str, default="HS25", help='ISS model (HS25, MSv1, …)')
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--chunks', type=int, default=10)
    ap.add_argument('--chunk-max-retries', type=int, default=10)
    ap.add_argument('--generic-pairs', type=int, default=30_000_000)
    ap.add_argument('--filtered-pairs', type=int, default=8_000_000)
    ap.add_argument('--shuffle', action='store_true', help='Produce shuffled copies of generic pairs')
    ap.add_argument('--cpus-per-iss', type=int, default=1, help='CPUs per ISS generate call')
    ap.add_argument('--comm-arg', type=str, default="ALL",
                    help='If "ALL", process all TrainC_short_*; else process only this community name.')
    args = ap.parse_args()

    # Resolve community path(s)
    if args.comm:
        base_comm_root = args.comm.parent.parent  # best-effort; unused
        comms = [args.comm]
    else:
        if not args.comm_name and args.comm_arg != "ALL":
            # Allow passing a single name via --comm-arg (like original bash positional)
            args.comm_name = args.comm_arg
        if args.comm_name:
            comms = [args.root / "train_communities" / args.comm_name]
        else:
            # ALL → glob TrainC_short_*
            comms = sorted((args.root / "train_communities").glob("TrainC_short_*"))

    # Tool checks
    need_exe(str(args.iss_bin), "iss")
    need_exe("gzip", "gzip")
    # shuffle uses paste+shuf; we'll check on demand

    for comm in comms:
        if not comm.exists():
            log(f"[WARN] Skip missing community: {comm}")
            continue

        log(f"====== [BEGIN] {comm.name} ======")
        refs_dir = comm / "refs"
        abund_dir = comm / "abundance"
        perrec_dir = comm / "per_record"
        out_dir = comm / "reads" / "short" / args.split
        ensure_dir(refs_dir); ensure_dir(perrec_dir); ensure_dir(out_dir)

        # Build class FASTAs with normalized headers
        class_fastas: Dict[str, Path] = {}
        for cls in CLASSES:
            fa = build_class_fasta(comm, cls, refs_dir)
            filter_fasta_minlen(fa, min_len=500)
            class_fastas[CLASS_LC_TO_CAP[cls]] = fa

        # Per-record weights from *generic* per-genome TSVs
        perrec_files: Dict[str, Path] = {}
        for cls in CLASSES:
            cls_cap = CLASS_LC_TO_CAP[cls]
            pergen_tsv = abund_dir / f"{cls}_generic.tsv"
            out_tsv = perrec_dir / f"{cls}_generic_perrecord.tsv"
            if not pergen_tsv.exists():
                raise RuntimeError(f"Missing per-genome TSV for {cls}: {pergen_tsv}")
            write_perrecord_weights(class_fastas[cls_cap], pergen_tsv, out_tsv)
            perrec_files[cls_cap] = out_tsv

        # Temporary work area
        tmp = out_dir / f".tmp_{os.environ.get('SLURM_JOB_ID') or os.getpid()}_{comm.name}"
        ensure_dir(tmp)

        try:
            # Class splits (pairs)
            G = args.generic_pairs
            F = args.filtered_pairs
            splits_generic = {
                "Prokaryote":     G * 56 // 100,
                "Microeukaryote": G * 24 // 100,
                "Virus":          G * 10 // 100,
                "Plasmid":        G * 10 // 100,
            }
            # distribute rounding remainder to Prokaryote
            splits_generic["Prokaryote"] += (G - sum(splits_generic.values()))

            splits_filtered = {
                "Prokaryote":     F * 14 // 100,
                "Microeukaryote": F *  6 // 100,
                "Virus":          F * 40 // 100,
                "Plasmid":        F * 40 // 100,
            }
            splits_filtered["Prokaryote"] += (F - sum(splits_filtered.values()))

            settings = IssChunkSettings(
                iss_bin=args.iss_bin, model=args.model, base_seed=args.seed,
                chunks=args.chunks, chunk_max_retries=args.chunk_max_retries,
                cpus_per_iss=args.cpus_per_iss,
            )

            # Generate GENERIC per class
            log(f"[GENERIC] target pairs={G}")
            prok_r1, prok_r2 = gen_group_chunked("Prokaryote", class_fastas["Prokaryote"],
                                                 perrec_files["Prokaryote"], splits_generic["Prokaryote"],
                                                 out_dir / "prok", tmp, settings)
            euk_r1, euk_r2   = gen_group_chunked("Microeukaryote", class_fastas["Microeukaryote"],
                                                 perrec_files["Microeukaryote"], splits_generic["Microeukaryote"],
                                                 out_dir / "euk", tmp, settings)
            vir_r1, vir_r2   = gen_group_chunked("Virus", class_fastas["Virus"],
                                                 perrec_files["Virus"], splits_generic["Virus"],
                                                 out_dir / "vir", tmp, settings)
            pla_r1, pla_r2   = gen_group_chunked("Plasmid", class_fastas["Plasmid"],
                                                 perrec_files["Plasmid"], splits_generic["Plasmid"],
                                                 out_dir / "pla", tmp, settings)

            gen_r1 = out_dir / "generic_R1.fastq.gz"
            gen_r2 = out_dir / "generic_R2.fastq.gz"
            gz_concat(gen_r1, [prok_r1, euk_r1, vir_r1, pla_r1])
            gz_concat(gen_r2, [prok_r2, euk_r2, vir_r2, pla_r2])
            assert_pairs("GENERIC", gen_r1, gen_r2, G)
            log(f"[GENERIC] {comm.name} pairs={count_pairs_fastq(gen_r1)}")

            # Optional shuffle
            if args.shuffle:
                need_exe("paste", "paste")
                need_exe("shuf", "shuf")
                shuf_r1 = out_dir / "generic_shuf_R1.fastq.gz"
                shuf_r2 = out_dir / "generic_shuf_R2.fastq.gz"
                shuffle_pairs_stream(args.seed, gen_r1, gen_r2, shuf_r1, shuf_r2)
                assert_pairs("GENERIC_SHUF", shuf_r1, shuf_r2, G)

            # FILTERED per class (re-using per-record weights from generic per-genome)
            log(f"[FILTERED] target pairs={F}")
            f_prok_r1, f_prok_r2 = gen_group_chunked("Prokaryote", class_fastas["Prokaryote"],
                                                     perrec_files["Prokaryote"], splits_filtered["Prokaryote"],
                                                     out_dir / "f_prok", tmp, settings)
            f_euk_r1, f_euk_r2   = gen_group_chunked("Microeukaryote", class_fastas["Microeukaryote"],
                                                     perrec_files["Microeukaryote"], splits_filtered["Microeukaryote"],
                                                     out_dir / "f_euk", tmp, settings)
            f_vir_r1, f_vir_r2   = gen_group_chunked("Virus", class_fastas["Virus"],
                                                     perrec_files["Virus"], splits_filtered["Virus"],
                                                     out_dir / "f_vir", tmp, settings)
            f_pla_r1, f_pla_r2   = gen_group_chunked("Plasmid", class_fastas["Plasmid"],
                                                     perrec_files["Plasmid"], splits_filtered["Plasmid"],
                                                     out_dir / "f_pla", tmp, settings)

            fil_r1 = out_dir / "filtered_R1.fastq.gz"
            fil_r2 = out_dir / "filtered_R2.fastq.gz"
            gz_concat(fil_r1, [f_prok_r1, f_euk_r1, f_vir_r1, f_pla_r1])
            gz_concat(fil_r2, [f_prok_r2, f_euk_r2, f_vir_r2, f_pla_r2])
            assert_pairs("FILTERED", fil_r1, fil_r2, F)
            log(f"[FILTERED] {comm.name} pairs={count_pairs_fastq(fil_r1)}")

            log(f"====== [END] {comm.name} ======")

        finally:
            # Clean temporary workspace
            try:
                shutil.rmtree(tmp)
            except Exception:
                pass

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log(f"[ERR] {e}")
        sys.exit(1)
