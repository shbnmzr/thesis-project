#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generic metaSPAdes runner for communities (generic & filtered).

Features:
- Runs metaSPAdes for each community (generic + filtered short reads).
- Resumes if prior SPAdes state exists; otherwise starts fresh.
- Chooses tmpdir from SLURM env or sensible fallbacks.
- Derives memory budget from SLURM and passes -m GB to SPAdes.
- Filters contigs to ≥500 bp.
- Invokes gfa_extract_validate.py and writes a sentinel _DONE.
- Skips steps if outputs already exist.
- Supports single target community (by name) or ALL (glob TrainC_short_*).

Example:
  python assemble_short_reads.py \
    --train-base /path/to/training/train_communities \
    --spades-bin /path/to/SPAdes-3.15.5-Linux/bin/metaspades.py \
    --gfa-script /path/to/scripts/gfa_extract_validate.py \
    --threads 16 --keep-seq 0 --only-assembler 0 \
    --target TrainC_short_01

Or for all:
  python assemble_short_reads.py --target ALL
"""
from __future__ import annotations
import argparse
import os
import sys
import subprocess as sp
from pathlib import Path

def log(msg: str) -> None:
    print(msg, flush=True)

def run(cmd, cwd=None, log_file: Path | None = None) -> int:
    env = os.environ.copy()
    if log_file is None:
        return sp.call(cmd, cwd=cwd, env=env)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("wb") as f:
        p = sp.Popen(cmd, cwd=cwd, env=env, stdout=sp.PIPE, stderr=sp.STDOUT)
        assert p.stdout is not None
        for line in iter(p.stdout.readline, b""):
            f.write(line)
            sys.stdout.buffer.write(line)
            sys.stdout.flush()
        return p.wait()

def calc_spades_mem_gb() -> int:
    # Follow bash logic: derive from SLURM_MEM_PER_NODE or SLURM_MEM_PER_CPU (MB), minus 8GB safety
    mb = os.environ.get("SLURM_MEM_PER_NODE") or os.environ.get("SLURM_MEM_PER_CPU") or "400000"
    try:
        gb = max(16, int(int(mb) / 1024))
    except ValueError:
        gb = 400000 // 1024
    if gb > 16:
        gb = gb - 8
    if gb < 16:
        gb = 16
    return gb

def choose_tmpdir() -> str:
    candidates = [
        os.environ.get("SLURM_TMPDIR", ""),
        f"/scratch/{os.environ.get('USER','user')}/{os.environ.get('SLURM_JOB_ID','')}",
        f"/tmp/{os.environ.get('USER','user')}_{os.environ.get('SLURM_JOB_ID','')}",
    ]
    for cand in candidates:
        if not cand:
            continue
        p = Path(cand)
        try:
            p.mkdir(parents=True, exist_ok=True)
        except Exception:
            continue
        if os.access(str(p), os.W_OK):
            return str(p)
    return ""

def spades_ok(asm_dir: Path) -> bool:
    contigs = asm_dir / "contigs.fasta"
    scaff = asm_dir / "scaffolds.fasta"
    logp = asm_dir / "spades.log"
    if not (contigs.exists() and contigs.stat().st_size > 0 and
            scaff.exists() and scaff.stat().st_size > 0 and
            logp.exists()):
        return False
    try:
        txt = logp.read_text(errors="ignore")
    except Exception:
        return False
    return "======= SPAdes pipeline finished" in txt

def spades_can_continue(asm_dir: Path) -> bool:
    return any([
        (asm_dir / "params.txt").exists(),
        (asm_dir / "dataset.info").exists(),
        any(x.is_dir() and x.name.startswith("K") for x in asm_dir.glob("K*"))
    ])

def filter_ge500(inp: Path, outp: Path) -> None:
    if outp.exists() and outp.stat().st_size > 0:
        return
    kept = total = 0
    with inp.open() as f, outp.open("w") as w:
        name = None; seq = []
        def flush():
            nonlocal kept, total, name, seq
            if name is None:
                return
            s = "".join(seq); total += 1
            if len(s) >= 500:
                w.write(f">{name}\n")
                for i in range(0, len(s), 80):
                    w.write(s[i:i+80] + "\n")
                kept += 1
        for line in f:
            if line.startswith(">"):
                flush()
                name = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        flush()
    log(f"[filter {outp}] kept {kept}/{total} contigs ≥500 bp")

def graph_extract(gfa_script: Path, asm_dir: Path, out_dir: Path, keep_seq: int, python_bin: Path | None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    sentinel = out_dir / "_DONE"
    if sentinel.exists():
        log(f"[SKIP] graph extract already done: {out_dir}")
        return
    py = str(python_bin) if python_bin else sys.executable
    cmd = [py, str(gfa_script), str(asm_dir), str(out_dir), str(keep_seq)]
    rc = run(cmd)
    if rc != 0:
        raise SystemExit(f"[ERR] gfa_extract_validate.py failed (rc={rc})")
    sentinel.touch()

def assemble_one_kind(comm_dir: Path, kind: str, spades_bin: Path, threads: int, only_assembler: int) -> None:
    r1 = comm_dir / "reads" / "short" / f"{kind}_R1.fastq.gz"
    r2 = comm_dir / "reads" / "short" / f"{kind}_R2.fastq.gz"
    if not (r1.exists() and r1.stat().st_size > 0 and r2.exists() and r2.stat().st_size > 0):
        log(f"[WARN] {kind}: reads missing in {comm_dir.name}; skipping")
        return

    asm_dir = comm_dir / f"assembly/metaSPAdes_{kind}"
    eval_dir = comm_dir / f"eval_{kind}"
    asm_dir.mkdir(parents=True, exist_ok=True)
    eval_dir.mkdir(parents=True, exist_ok=True)

    log_file = asm_dir / "spades.log"
    contigs = asm_dir / "contigs.fasta"

    if spades_ok(asm_dir):
        log(f"[SKIP] {kind}: SPAdes already completed.")
    else:
        tmpdir = choose_tmpdir()
        mem_gb = calc_spades_mem_gb()
        t = threads
        if spades_can_continue(asm_dir):
            log(f"[CONT] {kind}: Resuming metaSPAdes (t={t} mem={mem_gb}GB)")
            cmd = [
                str(spades_bin),
                "--continue",
                "-o", str(asm_dir),
                "-t", str(t),
                "-m", str(mem_gb),
            ]
            if tmpdir:
                cmd[1:1] = ["--tmp-dir", tmpdir]
            rc = run(cmd, log_file=log_file)
        else:
            log(f"[RUN ] {kind}: Fresh metaSPAdes (t={t} mem={mem_gb}GB only_asm={only_assembler})")
            # Clean previous partial outputs if any
            for child in asm_dir.iterdir():
                try:
                    if child.is_dir():
                        for c in child.rglob("*"):
                            pass
                except Exception:
                    pass
            # We let SPAdes overwrite in place; metaSPAdes expects clean dir
            # Build command
            cmd = [
                str(spades_bin),
                "-1", str(r1),
                "-2", str(r2),
                "--meta",
                "-o", str(asm_dir),
                "-t", str(t),
                "-m", str(mem_gb),
            ]
            if tmpdir:
                cmd[1:1] = ["--tmp-dir", tmpdir]
            if only_assembler == 1:
                cmd.append("--only-assembler")
            rc = run(cmd, log_file=log_file)

        if rc != 0 or not spades_ok(asm_dir):
            raise SystemExit(f"[ERR] {kind}: SPAdes did not finish cleanly. See {log_file}")

    # ≥500 bp contigs (idempotent)
    ge500 = eval_dir / "contigs_ge500.fasta"
    if ge500.exists() and ge500.stat().st_size > 0:
        log(f"[SKIP] {kind}: {ge500} exists")
    else:
        filter_ge500(contigs, ge500)

def run_one_comm(
    comm_dir: Path,
    spades_bin: Path,
    gfa_script: Path,
    keep_seq: int,
    only_assembler: int,
    threads: int,
    python_bin: Path | None,
) -> None:
    log(f"[COMM] {comm_dir.name}")
    # Assemble both kinds
    for kind in ("generic", "filtered"):
        assemble_one_kind(comm_dir, kind, spades_bin, threads, only_assembler)

        # After assembly, extract graph
        asm_dir = comm_dir / f"assembly/metaSPAdes_{kind}"
        gdir = asm_dir / "graph_for_gnn"
        graph_extract(gfa_script, asm_dir, gdir, keep_seq, python_bin)

def main() -> None:
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--train-base", type=Path, default=Path("/path/to/training/train_communities"),
                    help="Root containing TrainC_* community directories")
    ap.add_argument("--target", required=True, help="Community name (e.g., TrainC_short_01) or 'ALL'")
    ap.add_argument("--spades-bin", type=Path, required=True,
                    help="Path to metaspades.py")
    ap.add_argument("--gfa-script", type=Path, required=True,
                    help="Path to gfa_extract_validate.py")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "16")))
    ap.add_argument("--keep-seq", type=int, default=0, help="Pass-through arg to gfa_extract_validate.py")
    ap.add_argument("--only-assembler", type=int, default=0, help="1 to add --only-assembler to SPAdes")
    ap.add_argument("--python-bin", type=Path, default=None, help="Python to run gfa_extract_validate.py (default: this python)")
    args = ap.parse_args()

    # Validate tools
    if not args.spades_bin.exists():
        raise SystemExit(f"[ERR] SPAdes not found: {args.spades_bin}")
    if not args.gfa_script.exists():
        raise SystemExit(f"[ERR] gfa_extract_validate.py not found: {args.gfa_script}")

    if args.target == "ALL":
        # Iterate TrainC_short_* (mirrors original)
        found = False
        for d in sorted(args.train_base.glob("TrainC_short_*")):
            if d.is_dir():
                found = True
                run_one_comm(d, args.spades_bin, args.gfa_script, args.keep_seq, args.only_assembler, args.threads, args.python_bin)
        if not found:
            log("[WARN] No TrainC_short_* directories found.")
    else:
        comm_dir = args.train_base / args.target
        if not comm_dir.is_dir():
            raise SystemExit(f"[ERR] Community not found: {comm_dir}")
        run_one_comm(comm_dir, args.spades_bin, args.gfa_script, args.keep_seq, args.only_assembler, args.threads, args.python_bin)

    log("[DONE]")

if __name__ == "__main__":
    main()
