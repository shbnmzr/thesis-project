#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NanoSim community runner
========================================

This script:
- Builds/uses merged abundance files for two scenarios (generic & filtered)
- Prepares per-genome refset symlinks (or unzips .gz)
- Allocates read counts exactly to the requested total
- Probes the model with a 10-read test
- Runs per-genome NanoSim (optionally in parallel) with resume & per-task timeout
- Materializes outputs (handles .gz) and merges all parts into one FASTA/FASTQ
- Aligns final reads with minimap2 + samtools (map-ont)
"""
from __future__ import annotations
import argparse
import concurrent.futures as cf
import gzip
import os
import re
import shutil
import subprocess as sp
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

# ----------------------------- utils & logging -----------------------------

def log(msg: str) -> None:
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{now}] {msg}")
    sys.stdout.flush()


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


# ------------------------ abundance merge (if needed) -----------------------

def mk_merged_if_needed(comm: Path, scen: str, out: Path) -> None:
    """Merge class TSVs into one abundance file (genome \t weight)."""
    if out.exists() and out.stat().st_size > 0:
        return
    abr = comm / 'abundance'
    files = [
        abr / f'prokaryote_{scen}.tsv',
        abr / f'microeukaryote_{scen}.tsv',
        abr / f'virus_{scen}.tsv',
        abr / f'plasmid_{scen}.tsv',
    ]
    weights: defaultdict[str, float] = defaultdict(float)
    total = 0.0
    for fp in files:
        if not fp.exists():
            continue
        with fp.open('rt') as f:
            for ln in f:
                ln = ln.strip()
                if not ln:
                    continue
                toks = ln.split('\t')
                if len(toks) < 2:
                    continue
                genome, w = toks[0], float(toks[1])
                weights[genome] += w
                total += w
    if total <= 0:
        total = 1.0
    ensure_dir(out.parent)
    with out.open('wt') as fo:
        for g, w in weights.items():
            fo.write(f"{g}\t{w/total:.10g}\n")
    if not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"Merged abundance missing/empty: {out}")


# ---------------------------- model flag picker ----------------------------

def pick_model_flags(python_bin: Path, simulator: Path, model_dir: Path, scen_out: Path) -> Tuple[str, str, str]:
    """Detect NanoSim genome mode flags for model; support -md/--model_dir/-c prefix."""
    try:
        help_txt = sp.check_output([str(python_bin), str(simulator), 'genome', '-h'],
                                   stderr=sp.STDOUT, text=True)
    except sp.CalledProcessError as e:
        help_txt = e.output or ''
    if re.search(r'(^|\s)-md(\b|,)', help_txt):
        return ('-md', str(model_dir), 'dir')
    if '--model_dir' in help_txt:
        return ('--model_dir', str(model_dir), 'dir')
    # Fallback to prefix (-c): symlink matching prefix to scen_out/model_copy
    mp_candidates = list(model_dir.glob('*_model_profile'))
    if not mp_candidates:
        raise RuntimeError(f"No '*_model_profile' in {model_dir}")
    prefix = mp_candidates[0].name.replace('_model_profile', '')
    model_copy = scen_out / 'model_copy'
    if model_copy.exists():
        shutil.rmtree(model_copy)
    ensure_dir(model_copy)
    for f in model_dir.glob(prefix + '_*'):
        (model_copy / f.name).symlink_to(f)
    return ('-c', str(model_copy / prefix), 'prefix')


# ----------------------------- output helpers ------------------------------
READ_STEMS = ('', '_aligned_reads', '_unaligned_reads', '_chimeric_reads')
READ_EXTS = ('.fastq', '.fq', '.fasta', '.fa', '.fna')


def materialize_outputs(outprefix: Path) -> List[Path]:
    """Return list of existing (possibly gunzipped) read files for given outprefix."""
    found: List[Path] = []
    for stem in READ_STEMS:
        for ext in READ_EXTS:
            p = outprefix.with_name(outprefix.name + stem + ext)
            if p.exists() and p.stat().st_size > 0:
                found.append(p)
            gz = Path(str(p) + '.gz')
            if gz.exists() and gz.stat().st_size > 0:
                # materialize to plain ext (idempotent)
                with gzip.open(gz, 'rb') as fi, p.open('wb') as fo:
                    shutil.copyfileobj(fi, fo)
                found.append(p)
    if not found:
        # fallback glob
        for p in outprefix.parent.glob(outprefix.name + '*reads.*'):
            q = p
            if p.suffix == '.gz':
                q = p.with_suffix('')
                with gzip.open(p, 'rb') as fi, q.open('wb') as fo:
                    shutil.copyfileobj(fi, fo)
            if q.stat().st_size > 0:
                found.append(q)
    return found


def has_any_output(outprefix: Path) -> bool:
    try:
        return len(materialize_outputs(outprefix)) > 0
    except Exception:
        return False


# ----------------------------- core per-scenario ---------------------------

def sanitize_gid(path: Path) -> str:
    gid = path.name
    for suf in ('.gz', '.fasta', '.fa', '.fna'):
        if gid.endswith(suf):
            gid = gid[: -len(suf)]
    gid = re.sub(r'[^A-Za-z0-9]', '_', gid)
    return gid


def allocate_counts(pairs: List[Tuple[Path, float]], total_reads: int) -> List[Tuple[Path, int]]:
    # exact-sum allocation by sorted fractional parts
    weights = [w for _, w in pairs]
    tot = sum(weights) or 1.0
    exp = [total_reads * (w / tot) for w in weights]
    fl = [int(x) for x in exp]
    fr = [x - int(x) for x in exp]
    s = sum(fl)
    r = total_reads - s
    idx = list(range(len(pairs)))
    idx.sort(key=lambda i: fr[i], reverse=True)
    for k in range(r):
        fl[idx[k]] += 1
    return [(pairs[i][0], fl[i]) for i in range(len(pairs)) if fl[i] > 0]


def build_refset_and_pairs(abund_path: Path, scen_out: Path) -> Tuple[Path, Path, List[Tuple[Path, float]]]:
    refset = scen_out / 'refset'
    parts_dir = scen_out / 'reads' / 'parts'
    ensure_dir(refset); ensure_dir(parts_dir)
    pairs_tsv = parts_dir / 'pairs.tsv'
    pairs: List[Tuple[Path, float]] = []
    with abund_path.open('rt') as f, pairs_tsv.open('wt') as fo:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            genome, w = ln.split('\t')[:2]
            wv = float(w)
            base = Path(genome).name
            dest = refset / (base[:-3] if base.endswith('.gz') else base)
            if base.endswith('.gz'):
                if not dest.exists():
                    with gzip.open(genome, 'rb') as fi, dest.open('wb') as out:
                        shutil.copyfileobj(fi, out)
            else:
                if not dest.exists():
                    (refset / base).symlink_to(Path(genome))
            fo.write(f"{dest}\t{wv}\n")
            pairs.append((dest, wv))
    if not pairs:
        raise RuntimeError(f"No genome-weight pairs from {abund_path}")
    return refset, pairs_tsv, pairs


def concat_ref_for_mapping(refset: Path, scen_out: Path) -> Path:
    ref_all = scen_out / 'ref.all.fasta'
    if not ref_all.exists() or ref_all.stat().st_size == 0:
        with ref_all.open('wb') as fo:
            for f in sorted(refset.iterdir()):
                if f.is_file() and f.stat().st_size > 0:
                    with open(f, 'rb') as fi:
                        shutil.copyfileobj(fi, fo)
    if not ref_all.exists() or ref_all.stat().st_size == 0:
        raise RuntimeError('Empty concatenated reference')
    return ref_all


def run_cmd(cmd: Sequence[str], log_to: Optional[Path] = None, timeout_sec: Optional[int] = None) -> int:
    env = os.environ.copy()
    if log_to is None:
        return sp.call(cmd, env=env, timeout=timeout_sec)
    ensure_dir(log_to.parent)
    with log_to.open('wb') as f:
        try:
            p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.STDOUT, env=env)
            assert p.stdout is not None
            for line in iter(p.stdout.readline, b''):
                f.write(line)
                sys.stdout.buffer.write(line)
                sys.stdout.flush()
            return p.wait(timeout=timeout_sec)
        except sp.TimeoutExpired:
            p.kill()
            f.write(b'\n[TIMEOUT] Command exceeded timeout.\n')
            return 124


def run_probe(python_bin: Path, simulator: Path, ref: Path, model_flag: str, model_arg: str,
              outpfx: Path, seed: int, fastq_mode: bool) -> None:
    cmd = [str(python_bin), '-u', str(simulator), 'genome', '-rg', str(ref),
           model_flag, model_arg, '-o', str(outpfx), '--seed', str(seed), '-t', '1']
    if fastq_mode:
        cmd += ['--fastq', '-n', '10']
    else:
        cmd += ['-n', '10']
    log(f"Probe: {' '.join(cmd)}")
    rc = run_cmd(cmd)
    if rc != 0:
        raise RuntimeError(f"Probe failed (rc={rc})")
    outs = materialize_outputs(outpfx)
    if not outs:
        raise RuntimeError("Probe produced no files")


def simulate_one(args: Tuple[Path, int, int, int, int, Path, Path, str, str, Path, int, bool, Path, int]) -> Tuple[int, str, int]:
    (
        ref_path, nreads, idx, total_genomes, tpj,
        python_bin, simulator, model_flag, model_arg,
        outdir, seed, fastq_mode, scen_out, timeout_min
    ) = args
    gid = sanitize_gid(ref_path)
    outpfx = outdir / f"{idx}_{gid}"
    if has_any_output(outpfx):
        return (idx, gid, 0)
    cmd = [str(python_bin), '-u', str(simulator), 'genome', '-rg', str(ref_path),
           model_flag, model_arg, '-o', str(outpfx), '--seed', str(seed), '-t', str(tpj)]
    if fastq_mode:
        cmd += ['--fastq']
    cmd += ['-n', str(nreads)]
    log_to = outpfx.with_suffix('.log')
    log(f"(RUN) ({idx}/{total_genomes}) {gid} : {nreads} reads | t/job={tpj}")
    rc = run_cmd(cmd, log_to=log_to, timeout_sec=timeout_min*60 if timeout_min > 0 else None)
    if rc == 0:
        outs = materialize_outputs(outpfx)
        if not outs:
            rc = 6
    return (idx, gid, rc)


def merge_parts(parts_dir: Path, final_reads: Path, fastq_mode: bool) -> int:
    pattern_exts = ['*.fa', '*.fna', '*.fasta', '*.fastq', '*.fa.gz', '*.fna.gz', '*.fasta.gz', '*.fastq.gz']
    files: List[Path] = []
    for pat in pattern_exts:
        files += sorted(parts_dir.glob(pat))
    written = 0
    with final_reads.open('wb') as fo:
        for f in files:
            if not f.is_file() or f.stat().st_size == 0:
                continue
            if f.suffix == '.gz':
                with gzip.open(f, 'rb') as fi:
                    shutil.copyfileobj(fi, fo)
            else:
                with open(f, 'rb') as fi:
                    shutil.copyfileobj(fi, fo)
            written += 1
    if final_reads.stat().st_size == 0:
        raise RuntimeError(f"Final reads file empty: {final_reads}")
    return written


def align_all(minimap2: Path, samtools: Path, ref_all: Path, final_reads: Path, bam_out: Path,
              threads: int) -> Tuple[Path, Path]:
    paf = bam_out.with_suffix('').with_name(bam_out.stem.replace('.sorted', ''))
    paf = paf.parent / (paf.name.replace('.bam', '') + '.paf')
    ensure_dir(bam_out.parent)
    # PAF
    cmd_paf = [str(minimap2), '-t', str(threads), '-x', 'map-ont', str(ref_all), str(final_reads)]
    log(f"Align (PAF): {' '.join(cmd_paf)} -> {paf}")
    with paf.open('wb') as fo:
        rc = sp.call(cmd_paf, stdout=fo)
    if rc != 0:
        raise RuntimeError(f"minimap2 PAF rc={rc}")
    # BAM
    cmd_bam = [str(minimap2), '-t', str(threads), '-ax', 'map-ont', str(ref_all), str(final_reads)]
    log(f"Align (BAM): {' '.join(cmd_bam)} | sort -> {bam_out}")
    p1 = sp.Popen(cmd_bam, stdout=sp.PIPE)
    p2 = sp.Popen([str(samtools), 'sort', '-@', str(threads), '-o', str(bam_out)], stdin=p1.stdout)
    assert p1.stdout is not None
    p1.stdout.close()
    rc2 = p2.wait()
    if rc2 != 0:
        raise RuntimeError(f"samtools sort rc={rc2}")
    rc3 = sp.call([str(samtools), 'index', str(bam_out)])
    if rc3 != 0:
        raise RuntimeError(f"samtools index rc={rc3}")
    return paf, bam_out


# ----------------------------------- main ----------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--comm', required=True, type=Path, help='Community root directory')
    ap.add_argument('--model-dir', required=True, type=Path, help='NanoSim model directory or prefix source')
    ap.add_argument('--python-bin', type=Path, default=Path('/path/to/miniconda3/envs/nanosim/bin/python'))
    ap.add_argument('--simulator', type=Path, default=Path('/path/to/NanoSim/src/simulator.py'))
    ap.add_argument('--minimap2', type=Path, default=Path('/path/to/miniconda3/envs/minimap2_env/bin/minimap2'))
    ap.add_argument('--samtools', type=Path, default=Path('/path/to/miniconda3/envs/minimap2_env/bin/samtools'))

    ap.add_argument('--generic-total-reads', type=int, default=3_000_000)
    ap.add_argument('--filtered-total-reads', type=int, default=800_000)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--threads', type=int, default=16)
    ap.add_argument('--parallel-jobs', type=int, default=4)
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--dryrun-first', type=int, default=0, help='If >0, simulate only first N genomes')
    ap.add_argument('--fastq', action='store_true', help='Emit FASTQ instead of FASTA')
    ap.add_argument('--timeout-min', type=int, default=45, help='Per-genome timeout (minutes); 0 to disable')
    ap.add_argument('--split', choices=['train', 'val', 'test'], default='train',
                    help='Dataset split tag to segregate outputs (e.g., train/val/test)')

    args = ap.parse_args()

    # Resource hygiene
    for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS']:
        os.environ[k] = '1'
    os.environ['MALLOC_ARENA_MAX'] = '2'
    os.environ['PYTHONMALLOC'] = 'malloc'

    # Prefer node-local scratch for Python temp if set by user/SLURM
    tmpdir = os.environ.get('TMPDIR') or os.environ.get('TMP') or None
    if tmpdir:
        ensure_dir(Path(tmpdir))
        log(f"TMPDIR = {tmpdir}")

    comm = args.comm
    outroot_generic = comm / f"sim_nanosim_{args.split}_generic"
    outroot_filtered = comm / f"sim_nanosim_{args.split}_filtered"

    abund_generic = comm / 'abundance' / 'merged_generic.tsv'
    abund_filtered = comm / 'abundance' / 'merged_filtered.tsv'

    # Build merged abundance if missing
    mk_merged_if_needed(comm, 'generic', abund_generic)
    mk_merged_if_needed(comm, 'filtered', abund_filtered)

    def run_scenario(scen: str, abund: Path, total_reads: int, scen_out: Path) -> None:
        log(f"[INFO][{scen}] Abundance: {abund}")
        log(f"[INFO][{scen}] Outroot:   {scen_out}")
        ensure_dir(scen_out / 'reads' / 'parts')
        ensure_dir(scen_out / 'align')

        refset, _pairs_tsv, pairs = build_refset_and_pairs(abund, scen_out)
        ref_all = concat_ref_for_mapping(refset, scen_out)

        counts = allocate_counts(pairs, total_reads)
        total_genomes = len(counts)
        log(f"[INFO][{scen}] Genomes to simulate: {total_genomes}")

        # Pick model flags
        model_flag, model_arg, style = pick_model_flags(args.python_bin, args.simulator, args.model_dir, scen_out)
        log(f"[INFO][{scen}] Genome-mode model arg: {model_flag} {model_arg} (style={style})")

        # Probe 10 reads using first genome with nonzero count
        probe_ref = next((r for r, n in counts if n > 0), None)
        if probe_ref is None:
            raise RuntimeError('No genome with positive count for probe')
        probe_out = scen_out / 'reads' / 'parts' / '00_probe'
        log(f"[PRE][{scen}] 10-read probe (FASTQ_MODE={'1' if args.fastq else '0'})")
        run_probe(args.python_bin, args.simulator, probe_ref, model_flag, model_arg, probe_out, args.seed, args.fastq)
        log(f"[OK ][{scen}] Probe succeeded.")

        # Parallel setup
        par_jobs = max(1, int(args.parallel_jobs))
        tpj = max(1, args.threads // par_jobs)
        log(f"[INFO][{scen}] Parallel jobs: {par_jobs} | Threads/job: {tpj}")

        # Optionally limit to first N for dryrun
        work = counts[: args.dryrun_first] if args.dryrun_first and args.dryrun_first > 0 else counts

        parts_dir = scen_out / 'reads' / 'parts'
        # Run per-genome
        if par_jobs == 1:
            for idx, (ref_path, nreads) in enumerate(work, start=1):
                outpfx = parts_dir / f"{idx}_{sanitize_gid(ref_path)}"
                if args.resume and has_any_output(outpfx):
                    log(f"[RESUME][{scen}] ({idx}/{total_genomes}) {outpfx.name}")
                    continue
                _, gid, rc = simulate_one((ref_path, nreads, idx, total_genomes, tpj,
                                           args.python_bin, args.simulator, model_flag, model_arg,
                                           parts_dir, args.seed, args.fastq, scen_out, args.timeout_min))
                if rc != 0:
                    raise RuntimeError(f"Per-genome failed rc={rc} for {gid}")
        else:
            tasks = []
            with cf.ThreadPoolExecutor(max_workers=par_jobs) as ex:
                for idx, (ref_path, nreads) in enumerate(work, start=1):
                    outpfx = parts_dir / f"{idx}_{sanitize_gid(ref_path)}"
                    if args.resume and has_any_output(outpfx):
                        log(f"[RESUME][{scen}] ({idx}/{total_genomes}) {outpfx.name}")
                        continue
                    t = ex.submit(simulate_one, (ref_path, nreads, idx, total_genomes, tpj,
                                                 args.python_bin, args.simulator, model_flag, model_arg,
                                                 parts_dir, args.seed, args.fastq, scen_out, args.timeout_min))
                    tasks.append(t)
                for t in tasks:
                    idx, gid, rc = t.result()
                    if rc != 0:
                        raise RuntimeError(f"Per-genome failed rc={rc} for {gid}")
        log(f"[INFO][{scen}] Per-genome sims completed for {len(work)} genomes.")

        # Merge all parts
        final_reads = scen_out / 'reads' / ('sim.fastq' if args.fastq else 'sim.fasta')
        parts_merged = merge_parts(parts_dir, final_reads, args.fastq)
        log(f"[INFO][{scen}] Merged {parts_merged} files -> {final_reads}")

        # Align
        bam_out = scen_out / 'align' / 'sim_vs_ref.sorted.bam'
        paf, bam = align_all(args.minimap2, args.samtools, ref_all, final_reads, bam_out, args.threads)
        log(f"[DONE][{scen}] Reads: {final_reads}")
        log(f"[DONE][{scen}] PAF:   {paf}")
        log(f"[DONE][{scen}] BAM:   {bam}")

    # Run both scenarios
    run_scenario('generic', abund_generic, args.generic_total_reads, outroot_generic)
    run_scenario('filtered', abund_filtered, args.filtered_total_reads, outroot_filtered)


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        log(f"[ERR] {e}")
        sys.exit(1)
