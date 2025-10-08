#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generic community maker (short/long)

- Selects per-class genomes (non-overlapping across communities by default).
- Writes genomes/*.sel per class.
- Writes abundance/*.tsv: generic, filtered, long (per class) + merged_long.tsv.
- Saves a manifest.json per community.
- Tracks used genomes in OUTROOT/.used_genomes.json to avoid reuse (unless --reuse-allowed).

Example:
  python make_communities.py \
    --outroot /path/to/training/train_communities \
    --seed 42 \
    --n-short 8 --n-long 8 \
    --root-prok /path/to/dataset/train/archaea --root-prok /path/to/dataset/train/bacteria \
    --root-euk  /path/to/dataset/train/fungi   --root-euk  /path/to/dataset/train/protozoa \
    --root-vir /path/to/dataset/train/viral \
    --root-pla  /path/to/dataset/train/plasmid \
    --glob *.fna --glob *.fa --glob *.fasta --glob *.fna.gz --glob *.fa.gz --glob *.fasta.gz \
    --generic-prok 0.56 --generic-euk 0.24 --generic-vir 0.10 --generic-pla 0.10 \
    --filtered-prok 0.14 --filtered-euk 0.06 --filtered-vir 0.40 --filtered-pla 0.40 \
    --min-prok 100 --max-prok 200 \
    --min-euk   50 --max-euk  100 \
    --min-vir  120 --max-vir  200 \
    --min-pla  120 --max-pla  200 \
    --alpha-prok 0.7 --alpha-euk 0.9 --alpha-vir 1.0 --alpha-pla 1.0
"""

from __future__ import annotations
import argparse
import hashlib
import json
import os
import random
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

CLASSES = ["prokaryote", "microeukaryote", "virus", "plasmid"]

def log(msg: str) -> None:
    print(msg, flush=True)

def md5_of_path(p: Path) -> str:
    h = hashlib.md5()
    h.update(str(p).encode("utf-8", "ignore"))
    return h.hexdigest()

def canon(p: Path) -> str | None:
    try:
        rp = p.resolve(strict=True)
    except FileNotFoundError:
        return None
    if rp.stat().st_size <= 0:
        return None
    return str(rp)

def list_fastas(roots: List[str], patterns: List[str]) -> List[str]:
    files: set[str] = set()
    for root in roots:
        base = Path(root)
        if not base.exists():
            continue
        for pat in patterns:
            for q in base.rglob(pat):
                if q.is_file():
                    rp = canon(q)
                    if rp is not None:
                        files.add(rp)
    return sorted(files)

def dirichlet_weights(n: int, alpha: float, rng: np.random.Generator) -> np.ndarray:
    if n <= 0:
        return np.zeros(0, dtype=np.float64)
    if alpha <= 0:
        return np.full(n, 1.0 / n, dtype=np.float64)
    w = rng.dirichlet([alpha] * n).astype(np.float64)
    return w / w.sum()

def write_lines(p: Path, lines: List[str]) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w") as f:
        for x in lines:
            f.write(str(x).strip() + "\n")

def write_tsv_pairs(
    p: Path, paths: List[str], weights_sum_to: float, alpha: float, rng: np.random.Generator, reuse_weights: np.ndarray | None = None
) -> np.ndarray:
    p.parent.mkdir(parents=True, exist_ok=True)
    n = len(paths)
    if n == 0:
        p.write_text("")
        return np.zeros(0, dtype=np.float64)
    if reuse_weights is None:
        w = dirichlet_weights(n, alpha, rng)
    else:
        w = np.asarray(reuse_weights, dtype=np.float64)
        if w.size != n:
            raise ValueError("reuse_weights size mismatch")
        w = np.maximum(w, 0)
        if w.sum() == 0:
            w = np.full(n, 1.0 / n, dtype=np.float64)
        w = w / w.sum()
    w_out = w * float(weights_sum_to)
    with p.open("w") as f:
        for path, wi in zip(paths, w_out):
            f.write(f"{path}\t{wi:.10g}\n")
    return w

def load_pairs(p: Path) -> List[Tuple[str, float]]:
    out: List[Tuple[str, float]] = []
    if p.exists():
        with p.open() as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                a = line.split("\t")
                if len(a) >= 2:
                    out.append((a[0], float(a[1])))
    return out

def load_used_registry(outroot: Path) -> Dict[str, List[str]]:
    used_file = outroot / ".used_genomes.json"
    if used_file.exists():
        try:
            return json.loads(used_file.read_text())
        except Exception:
            return {c: [] for c in CLASSES}
    return {c: [] for c in CLASSES}

def save_used_registry(outroot: Path, reg: Dict[str, List[str]]) -> None:
    (outroot / ".used_genomes.json").write_text(json.dumps(reg, indent=2))

def partition_nonoverlapping(files: List[str], n_bins: int, prng: random.Random) -> List[List[str]]:
    arr = list(files)
    prng.shuffle(arr)
    bins = [[] for _ in range(n_bins)]
    for i, f in enumerate(arr):
        bins[i % n_bins].append(f)
    return bins

def make_one_community(
    name: str,
    outroot: Path,
    picks: Dict[str, List[str]],
    force: bool,
    rng: np.random.Generator,
    ratios_generic: Dict[str, float],
    ratios_filtered: Dict[str, float],
    alpha: Dict[str, float],
    seed: int,
    reuse_allowed: bool,
    roots: Dict[str, List[str]],
) -> None:
    base = outroot / name
    genomes_dir = base / "genomes"
    abr = base / "abundance"
    genomes_dir.mkdir(parents=True, exist_ok=True)
    abr.mkdir(parents=True, exist_ok=True)

    must_have_sels = [genomes_dir / f"{cls}.sel" for cls in CLASSES]
    must_have_abund = [
        "prokaryote_generic",
        "microeukaryote_generic",
        "virus_generic",
        "plasmid_generic",
        "prokaryote_filtered",
        "microeukaryote_filtered",
        "virus_filtered",
        "plasmid_filtered",
        "prokaryote_long",
        "microeukaryote_long",
        "virus_long",
        "plasmid_long",
        "merged_long",
    ]

    if (
        not force
        and all(p.exists() and p.stat().st_size > 0 for p in must_have_sels)
        and all((abr / f"{k}.tsv").exists() and (abr / f"{k}.tsv").stat().st_size > 0 for k in must_have_abund)
    ):
        log(f"[SKIP] {name}: selections + abundance already present")
        return

    for cls in CLASSES:
        sel_p = genomes_dir / f"{cls}.sel"
        if force or not (sel_p.exists() and sel_p.stat().st_size > 0):
            write_lines(sel_p, picks[cls])

    # Dirichlet weights per class (reused across generic/filtered/long)
    W = {cls: dirichlet_weights(len(picks[cls]), alpha[cls], rng) for cls in CLASSES}

    # generic
    write_tsv_pairs(abr / "prokaryote_generic.tsv", picks["prokaryote"], ratios_generic["prokaryote"], alpha["prokaryote"], rng, reuse_weights=W["prokaryote"])
    write_tsv_pairs(abr / "microeukaryote_generic.tsv", picks["microeukaryote"], ratios_generic["microeukaryote"], alpha["microeukaryote"], rng, reuse_weights=W["microeukaryote"])
    write_tsv_pairs(abr / "virus_generic.tsv", picks["virus"], ratios_generic["virus"], alpha["virus"], rng, reuse_weights=W["virus"])
    write_tsv_pairs(abr / "plasmid_generic.tsv", picks["plasmid"], ratios_generic["plasmid"], alpha["plasmid"], rng, reuse_weights=W["plasmid"])

    # filtered
    write_tsv_pairs(abr / "prokaryote_filtered.tsv", picks["prokaryote"], ratios_filtered["prokaryote"], alpha["prokaryote"], rng, reuse_weights=W["prokaryote"])
    write_tsv_pairs(abr / "microeukaryote_filtered.tsv", picks["microeukaryote"], ratios_filtered["microeukaryote"], alpha["microeukaryote"], rng, reuse_weights=W["microeukaryote"])
    write_tsv_pairs(abr / "virus_filtered.tsv", picks["virus"], ratios_filtered["virus"], alpha["virus"], rng, reuse_weights=W["virus"])
    write_tsv_pairs(abr / "plasmid_filtered.tsv", picks["plasmid"], ratios_filtered["plasmid"], alpha["plasmid"], rng, reuse_weights=W["plasmid"])

    # long (sum to 1.0 within class)
    write_tsv_pairs(abr / "prokaryote_long.tsv", picks["prokaryote"], 1.0, alpha["prokaryote"], rng, reuse_weights=W["prokaryote"])
    write_tsv_pairs(abr / "microeukaryote_long.tsv", picks["microeukaryote"], 1.0, alpha["microeukaryote"], rng, reuse_weights=W["microeukaryote"])
    write_tsv_pairs(abr / "virus_long.tsv", picks["virus"], 1.0, alpha["virus"], rng, reuse_weights=W["virus"])
    write_tsv_pairs(abr / "plasmid_long.tsv", picks["plasmid"], 1.0, alpha["plasmid"], rng, reuse_weights=W["plasmid"])

    # merged_long.tsv (scaled by GENERIC class mix), normalized to 1.0
    prok = load_pairs(abr / "prokaryote_long.tsv")
    euk = load_pairs(abr / "microeukaryote_long.tsv")
    vir = load_pairs(abr / "virus_long.tsv")
    pla = load_pairs(abr / "plasmid_long.tsv")

    merged_rows: List[Tuple[str, float]] = []
    for p, w in prok:
        merged_rows.append((p, w * ratios_generic["prokaryote"]))
    for p, w in euk:
        merged_rows.append((p, w * ratios_generic["microeukaryote"]))
    for p, w in vir:
        merged_rows.append((p, w * ratios_generic["virus"]))
    for p, w in pla:
        merged_rows.append((p, w * ratios_generic["plasmid"]))

    agg: Dict[str, float] = {}
    for p, w in merged_rows:
        agg[p] = agg.get(p, 0.0) + w
    tot = sum(agg.values()) or 1.0
    with (abr / "merged_long.tsv").open("w") as f:
        for p, w in agg.items():
            f.write(f"{p}\t{w / tot:.10g}\n")

    # manifest
    manifest = {
        "name": name,
        "seed": seed,
        "reuse_allowed": bool(reuse_allowed),
        "class_roots": roots,
        "ratios_generic": ratios_generic,
        "ratios_filtered": ratios_filtered,
        "ranges": {k: [len(picks[k]), len(picks[k])] for k in CLASSES},  # actual counts (min=max=current)
        "alpha": alpha,
        "files_md5": {cls: {str(p): md5_of_path(Path(p)) for p in picks[cls]} for cls in CLASSES},
    }
    (base / "manifest.json").write_text(json.dumps(manifest, indent=2))

    log(f"[OK] {name}: prok={len(picks['prokaryote'])} euk={len(picks['microeukaryote'])} vir={len(picks['virus'])} pla={len(picks['plasmid'])}")

def main() -> None:
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--outroot", type=Path, required=True, help="Where to create communities (TrainC_* subdirs)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--force", action="store_true", help="Overwrite/refresh selections and abundance")
    ap.add_argument("--reuse-allowed", action="store_true", help="Allow reusing genomes across communities")

    ap.add_argument("--n-short", type=int, default=8, help="Number of short communities")
    ap.add_argument("--n-long", type=int, default=8, help="Number of long communities")
    ap.add_argument("--short-prefix", type=str, default="TrainC_short", help="Naming prefix for short")
    ap.add_argument("--long-prefix", type=str, default="TrainC_long", help="Naming prefix for long")

    # Roots (repeatable)
    ap.add_argument("--root-prok", action="append", default=[], help="Add a prokaryote root dir (repeatable)")
    ap.add_argument("--root-euk",  action="append", default=[], help="Add a microeukaryote root dir (repeatable)")
    ap.add_argument("--root-vir",  action="append", default=[], help="Add a virus root dir (repeatable)")
    ap.add_argument("--root-pla",  action="append", default=[], help="Add a plasmid root dir (repeatable)")

    # Glob patterns (repeatable)
    ap.add_argument("--glob", action="append", default=["*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fa.gz", "*.fasta.gz"],
                    help="File glob patterns to include (repeatable)")

    # Ranges per class
    ap.add_argument("--min-prok", type=int, default=100); ap.add_argument("--max-prok", type=int, default=200)
    ap.add_argument("--min-euk",  type=int, default=50);  ap.add_argument("--max-euk",  type=int, default=100)
    ap.add_argument("--min-vir",  type=int, default=120); ap.add_argument("--max-vir",  type=int, default=200)
    ap.add_argument("--min-pla",  type=int, default=120); ap.add_argument("--max-pla",  type=int, default=200)

    # Ratios (must sum to 1 within each group)
    ap.add_argument("--generic-prok", type=float, default=0.56)
    ap.add_argument("--generic-euk",  type=float, default=0.24)
    ap.add_argument("--generic-vir",  type=float, default=0.10)
    ap.add_argument("--generic-pla",  type=float, default=0.10)

    ap.add_argument("--filtered-prok", type=float, default=0.14)
    ap.add_argument("--filtered-euk",  type=float, default=0.06)
    ap.add_argument("--filtered-vir",  type=float, default=0.40)
    ap.add_argument("--filtered-pla",  type=float, default=0.40)

    # Dirichlet alphas
    ap.add_argument("--alpha-prok", type=float, default=0.7)
    ap.add_argument("--alpha-euk",  type=float, default=0.9)
    ap.add_argument("--alpha-vir",  type=float, default=1.0)
    ap.add_argument("--alpha-pla",  type=float, default=1.0)

    args = ap.parse_args()

    outroot: Path = args.outroot
    outroot.mkdir(parents=True, exist_ok=True)

    # Validate ratios
    generic = {
        "prokaryote": args.generic_prok,
        "microeukaryote": args.generic_euk,
        "virus": args.generic_vir,
        "plasmid": args.generic_pla,
    }
    filtered = {
        "prokaryote": args.filtered_prok,
        "microeukaryote": args.filtered_euk,
        "virus": args.filtered_vir,
        "plasmid": args.filtered_pla,
    }
    for name, vec in [("GENERIC", generic), ("FILTERED", filtered)]:
        s = sum(vec.values())
        if abs(s - 1.0) > 1e-6:
            raise SystemExit(f"[ERR] {name} ratios must sum to 1.0 (got {s})")

    ranges = {
        "prokaryote": (args.min_prok, args.max_prok),
        "microeukaryote": (args.min_euk, args.max_euk),
        "virus": (args.min_vir, args.max_vir),
        "plasmid": (args.min_pla, args.max_pla),
    }
    roots = {
        "prokaryote": args.root_prok,
        "microeukaryote": args.root_euk,
        "virus": args.root_vir,
        "plasmid": args.root_pla,
    }
    patterns = args.glob

    rng = np.random.default_rng(args.seed)
    prng = random.Random(args.seed)

    # Load registry of used paths to avoid reuse unless allowed
    USED = load_used_registry(outroot)

    # Build pools
    files_by_cls: Dict[str, List[str]] = {}
    for cls in CLASSES:
        all_files = list_fastas(roots[cls], patterns)
        if not all_files:
            raise SystemExit(f"[ERR] No FASTAs for '{cls}' under roots: {roots[cls]}")
        if not args.reuse_allowed:
            prior = set(USED.get(cls, []))
            all_files = [p for p in all_files if p not in prior]
        files_by_cls[cls] = all_files
        log(f"[POOL] {cls}: {len(all_files)} genomes available")

    # Names
    short_names = [f"{args.short_prefix}_{i:02d}" for i in range(1, args.n_short + 1)]
    long_names = [f"{args.long_prefix}_{i:02d}" for i in range(1, args.n_long + 1)]
    all_names = short_names + long_names
    if not all_names:
        log("[WARN] No communities requested")
        return

    # Partition each class pool into non-overlapping shards across communities
    shards_by_cls: Dict[str, List[List[str]]] = {
        cls: partition_nonoverlapping(files_by_cls[cls], len(all_names), prng) for cls in CLASSES
    }

    # Plan picks per community
    planned_picks: Dict[str, Dict[str, List[str]]] = {name: {cls: [] for cls in CLASSES} for name in all_names}
    for idx, name in enumerate(all_names):
        for cls in CLASSES:
            lo, hi = ranges[cls]
            shard = shards_by_cls[cls][idx] if idx < len(shards_by_cls[cls]) else []
            shard2 = list(shard)
            prng.shuffle(shard2)

            if len(shard2) < lo:
                if args.reuse_allowed:
                    pool = list_fastas(roots[cls], patterns)
                    prng.shuffle(pool)
                    extra = [p for p in pool if p not in shard2]
                    need = lo - len(shard2)
                    shard2 += extra[:need]
                else:
                    raise SystemExit(
                        f"[ERR] '{cls}' shard for {name} too small "
                        f"({len(shard2)} < {lo}) and --reuse-allowed not set."
                    )

            k = min(max(lo, min(hi, len(shard2))), len(shard2))
            picks_cls = sorted(shard2[:k])
            planned_picks[name][cls] = picks_cls

    # Write selections/abundances for each community
    for name in all_names:
        # Update USED registry
        for cls in CLASSES:
            USED.setdefault(cls, [])
            for p in planned_picks[name][cls]:
                if p not in USED[cls]:
                    USED[cls].append(p)

        make_one_community(
            name=name,
            outroot=outroot,
            picks=planned_picks[name],
            force=args.force,
            rng=rng,
            ratios_generic=generic,
            ratios_filtered=filtered,
            alpha={
                "prokaryote": args.alpha_prok,
                "microeukaryote": args.alpha_euk,
                "virus": args.alpha_vir,
                "plasmid": args.alpha_pla,
            },
            seed=args.seed,
            reuse_allowed=args.reuse_allowed,
            roots=roots,
        )

    save_used_registry(outroot, USED)
    log(f"[DONE] Communities prepared under {outroot}/")

if __name__ == "__main__":
    main()
