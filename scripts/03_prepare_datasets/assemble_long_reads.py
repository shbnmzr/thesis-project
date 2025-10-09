#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Assemble ONT long reads with Flye

Usage examples:
  # 1) With a single reads file (FASTA/FASTQ, plain or .gz)
  python assemble_long_reads.py --reads /path/to/sim.fastq.gz --threads 16 --outroot /path/to/work/flye

  # 2) With a directory containing multiple read files (any mix of fasta/fastq, gz or plain)
  python assemble_long_reads.py --reads /path/to/dir --threads 32

  # 3) Auto-locate reads from a community (sim_nanosim_<split>_<scenario>/reads/sim.*)
  python assemble_long_reads.py --comm /nfsd/.../TrainC_long_01 --split train --scenario generic

Notes:
  * Temp dir is created UNDER THE CURRENT DIRECTORY:  ./.tmp_<JOBID-or-PID>
    - Cleaned on success; kept on failure for inspection.
  * Outputs are copied to --outroot (default: <comm>/work/flye_<split>_<scenario> or ./flye_out)
  * Requires Flye installed and on PATH (or pass --flye /path/to/flye)
"""

from __future__ import annotations
import argparse
import atexit
import gzip
import hashlib
import os
import shutil
import subprocess as sp
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, List, Optional

# ---------- small utils ----------

def log(msg: str) -> None:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def is_gz(p: Path) -> bool:
    return p.suffix == ".gz"

def ext_key(p: Path) -> str:
    return p.suffix.lower()

def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def which(cmd: str) -> Optional[str]:
    return shutil.which(cmd)

def read_sources_from_dir(d: Path) -> List[Path]:
    pats = [
        "*.fastq", "*.fq", "*.fasta", "*.fa", "*.fna",
        "*.fastq.gz", "*.fq.gz", "*.fasta.gz", "*.fa.gz", "*.fna.gz",
    ]
    files: List[Path] = []
    for pat in pats:
        files.extend(sorted(d.glob(pat)))
    return files

def open_maybe_gz(p: Path):
    if is_gz(p):
        return gzip.open(p, "rb")
    return p.open("rb")

# ---------- staging ----------

def stage_reads_to_single_file(reads: Path, staged_out: Path) -> None:
    """Concatenate one file or all files in a directory into staged_out (decompressing .gz on the fly)."""
    ensure_dir(staged_out.parent)
    if reads.is_dir():
        sources = read_sources_from_dir(reads)
        if not sources:
            raise FileNotFoundError(f"No read files found in directory: {reads}")
        log(f"Staging {len(sources)} files from directory: {reads}")
    else:
        if not reads.exists():
            raise FileNotFoundError(f"Reads not found: {reads}")
        sources = [reads]
        log(f"Staging single reads file: {reads}")

    # Stream-concat to minimize tmp usage
    with staged_out.open("wb") as outfh:
        for src in sources:
            log(f"  + {src}")
            with open_maybe_gz(src) as infh:
                shutil.copyfileobj(infh, outfh, length=1024 * 1024)
    if staged_out.stat().st_size == 0:
        raise RuntimeError(f"Staged reads is empty: {staged_out}")

# ---------- community auto-locate ----------

def find_sim_reads_from_comm(comm: Path, split: str, scenario: str) -> Path:
    base = comm / f"sim_nanosim_{split}_{scenario}" / "reads"
    candidates = [
        base / "sim.fastq",
        base / "sim.fastq.gz",
        base / "sim.fasta",
        base / "sim.fasta.gz",
    ]
    for c in candidates:
        if c.exists() and c.is_file():
            return c
    raise FileNotFoundError(f"Could not find simulated reads under {base} (expected sim.[fastq|fasta][.gz])")

# ---------- compression of artifacts ----------

def compress_copy(src: Path, dst_gz: Path) -> None:
    # Prefer pigz if available
    pigz = which("pigz")
    if pigz:
        with dst_gz.open("wb") as fo:
            p = sp.Popen([pigz, "-9", "-c", str(src)], stdout=fo)
            rc = p.wait()
            if rc != 0:
                raise RuntimeError(f"pigz failed (rc={rc}) on {src}")
        return
    # Fallback to Python gzip (slower)
    with src.open("rb") as fi, gzip.open(dst_gz, "wb", compresslevel=9) as fo:
        shutil.copyfileobj(fi, fo)

# ---------- manifest & checksums ----------

def write_manifest_and_checksums(outroot: Path, meta: dict, files: Iterable[Path]) -> None:
    manifest = outroot / "MANIFEST.txt"
    checks = outroot / "CHECKSUMS.sha256"
    ensure_dir(outroot)

    # Manifest
    with manifest.open("w") as f:
        f.write("# Flye minimal artifact manifest\n")
        for k, v in meta.items():
            f.write(f"{k}: {v}\n")
        f.write("\n## Files:\n")
        for p in files:
            if p.exists():
                f.write(f" - {p.name}\n")

    # Checksums
    with checks.open("w") as f:
        for p in files:
            if p.exists():
                f.write(f"{sha256_file(p)}  {p.name}\n")

# ---------- main ----------

def main() -> None:
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--reads", type=Path, help="Reads file OR directory (fasta/fastq; .gz ok)")
    ap.add_argument("--comm", type=Path, help="Community directory (auto-find simulated reads)")
    ap.add_argument("--comm-name", type=str, help="Name under BASE_COMM_ROOT to resolve --comm")
    ap.add_argument("--base-comm-root", type=Path,
                    default=Path("/path/to/train_communities"),
                    help="Root used when --comm-name is provided")
    ap.add_argument("--split", choices=["train", "val", "test"], default="train")
    ap.add_argument("--scenario", choices=["generic", "filtered"], default="generic")

    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--flye", type=str, default="flye", help="Flye executable")
    ap.add_argument("--outroot", type=Path, default=None, help="Final results directory")
    ap.add_argument("--tmpdir", type=Path, default=None,
                    help="Override temp root (default: ./.tmp_<JOBID-or-PID> under CWD)")
    ap.add_argument("--min-overlap", type=int, default=None, help="Optional Flye --min-overlap value")
    ap.add_argument("--keep-temp-on-success", action="store_true",
                    help="Do NOT delete temp dir on success (for debugging)")

    args = ap.parse_args()

    # Resolve community path if only name was provided
    comm: Optional[Path] = args.comm
    if comm is None and args.comm_name:
        comm = args.base_comm_root / args.comm_name

    # Determine reads source
    reads: Optional[Path] = args.reads
    if reads is None:
        if comm is None:
            ap.error("Provide --reads OR (--comm/--comm-name and --split/--scenario)")
        if not comm.exists():
            ap.error(f"--comm not found: {comm}")
        reads = find_sim_reads_from_comm(comm, args.split, args.scenario)

    # Determine outroot
    outroot: Path
    if args.outroot is not None:
        outroot = args.outroot
    else:
        if comm is not None:
            outroot = comm / f"work/flye_{args.split}_{args.scenario}"
        else:
            outroot = Path.cwd() / "flye_out"

    ensure_dir(outroot)
    ensure_dir(Path("logs/flye"))

    # Build temp root inside current directory
    jobid = os.environ.get("SLURM_JOB_ID")
    tag = jobid if jobid else str(os.getpid())
    tmp_root = args.tmpdir or (Path.cwd() / f".tmp_{tag}")
    stage = tmp_root / "stage"
    flye_out = tmp_root / "flye_out"
    ensure_dir(stage)
    ensure_dir(flye_out)

    # Make sure tools that respect TMPDIR use our local dir
    os.environ["TMPDIR"] = str(tmp_root)

    # Announce
    log(f"Job tag:     {tag}")
    log(f"Threads:     {args.threads}")
    log(f"Reads:       {reads}")
    if comm:
        log(f"Community:   {comm} (split={args.split}, scenario={args.scenario})")
    log(f"Outroot:     {outroot}")
    log(f"TEMP Root:   {tmp_root}")

    # Cleanup handlers
    def on_exit_success():
        if args.keep_temp_on_success:
            log(f"[keep] Temp preserved at: {tmp_root}")
        else:
            shutil.rmtree(tmp_root, ignore_errors=True)

    def on_exit_failure():
        fail_bundle = outroot / f"_failed_{tag}"
        ensure_dir(fail_bundle)
        # Try to copy useful bits
        for name in ("flye.log", "params.json"):
            src = flye_out / name
            if src.exists():
                shutil.copy2(src, fail_bundle / src.name)
        log(f"[fail] Temp kept for inspection: {tmp_root}")
        log(f"[fail] Logs copied to: {fail_bundle}")

    # We handle cleanup manually to differentiate success/failure
    atexit.register(lambda: None)

    try:
        # Stage reads to a single file
        staged_reads = stage / "reads.fq"
        stage_reads_to_single_file(reads, staged_reads)

        # Build Flye command
        cmd = [
            args.flye,
            "--nano-raw", str(staged_reads),
            "--meta",
            "--threads", str(args.threads),
            "--out-dir", str(flye_out),
        ]
        if args.min_overlap:
            cmd += ["--min-overlap", str(args.min_overlap)]

        log(f"RUN: {' '.join(cmd)}")
        rc = sp.call(cmd)
        if rc != 0:
            raise RuntimeError(f"Flye failed (rc={rc})")

        # Collect minimal artifacts
        ensure_dir(outroot)
        def copy_if_exists(name: str) -> Optional[Path]:
            src = flye_out / name
            if src.exists():
                dst = outroot / name
                shutil.copy2(src, dst)
                return dst
            return None

        copied: List[Path] = []
        for core in ("assembly.fasta", "assembly_graph.gfa", "assembly_info.txt", "params.json", "flye.log"):
            p = copy_if_exists(core)
            if p:
                copied.append(p)

        # Compress copies alongside (keep plain too)
        gz_copies: List[Path] = []
        for p in copied:
            gz = outroot / f"{p.name}.gz"
            try:
                compress_copy(p, gz)
                gz_copies.append(gz)
            except Exception as e:
                log(f"[warn] Compression failed for {p.name}: {e}")

        # Manifest + checksums
        meta = {
            "job_id": tag,
            "date_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "threads": str(args.threads),
            "reads_source": str(reads),
            "temp_root": str(tmp_root),
        }
        if comm:
            meta["community"] = str(comm)
            meta["split"] = args.split
            meta["scenario"] = args.scenario

        files_for_manifest = copied + gz_copies
        write_manifest_and_checksums(outroot, meta, files_for_manifest)

        log("[DONE] Flye (ONT+meta) finished successfully.")
        log(f"Results: {outroot}")
        on_exit_success()

    except Exception as e:
        log(f"[ERR] {e}")
        on_exit_failure()
        sys.exit(1)


if __name__ == "__main__":
    main()
