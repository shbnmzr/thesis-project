#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generic graph builder for GNN training data
===========================================

Replaces the bash job that built graphs for short reads. Now generic:
- Supports presets (short_metaspades, long_flye) or explicit overrides
- Builds per split (train/val/...) and per kind (generic/filtered/...)
- Writes per-run derived YAML configs from a base config
- Calls a module (default: data_prep.build_graphs_short) to materialize graphs
- Idempotent via stamp files

Examples
--------
Short-read SPAdes (like your old script):
  python build_graphs.py \
    --repo /path/to/training_kit/gnn_short_reads \
    --base-config /path/to/training_kit/config.yaml \
    --splits train val \
    --root-train /path/to/train_communities \
    --root-val  /path/to/val_communities \
    --kinds generic filtered \
    --preset short_metaspades

Long-read Flye (example layout):
  python build_graphs.py \
    --repo /path/to/training_kit/gnn_long_reads \
    --base-config /path/to/training_kit/config.yaml \
    --splits train \
    --root-train /path/to/train_communities \
    --kinds generic \
    --preset long_flye

Fully custom (explicit paths are relative to each community dir):
  python build_graphs.py \
    --repo /path/to/training_kit/gnn_any_reads \
    --base-config /path/to/training_kit/config.yaml \
    --splits train \
    --root-train /path/to/train_communities \
    --kinds filtered \
    --preset custom \
    --contigs-path "eval_filtered/contigs.fa" \
    --gfa-candidates "assembly/myasm_filtered/graph.gfa,assembly/myasm_filtered/graph.fastg" \
    --save-pt-dir "graphs/myasm_{kind}"
"""
from __future__ import annotations

import argparse
import os
import sys
import subprocess as sp
from pathlib import Path
from typing import Dict, List, Tuple


def log(*a):
    print(*a, flush=True)


def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def symlink_force(src: Path, dst: Path):
    try:
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        dst.symlink_to(src)
        return True
    except Exception:
        # Fallback to copy when symlink not allowed
        try:
            import shutil
            shutil.copyfile(src, dst)
            return True
        except Exception:
            return False


def preset_defaults(preset: str, kind: str) -> Dict[str, object]:
    if preset == "short_metaspades":
        return {
            "contigs_path": f"eval_{kind}/contigs_ge500.fasta",
            "gfa_candidates": [
                f"assembly/metaSPAdes_{kind}/assembly_graph_with_scaffolds.gfa",
                f"assembly/metaSPAdes_{kind}/assembly_graph.gfa",
                f"assembly/metaSPAdes_{kind}/assembly_graph.fastg",
            ],
            "save_pt_dir": f"graphs/short_{kind}",
        }
    if preset == "long_flye":
        return {
            "contigs_path": f"assembly/flye_{kind}/assembly.fasta",
            "gfa_candidates": [
                f"assembly/flye_{kind}/assembly_graph.gfa",
            ],
            "save_pt_dir": f"graphs/long_{kind}",
        }
    # custom: require explicit contigs/gfa via CLI (save dir still defaults)
    return {
        "contigs_path": None,
        "gfa_candidates": [],
        "save_pt_dir": f"graphs/{kind}",
    }


def write_cfg(
    base_cfg: Path,
    repo_dir: Path,
    root: Path,
    kind: str,
    out_yaml: Path,
    preset: str,
    contigs_override: str | None,
    gfa_override: str | None,
    save_dir_override: str | None,
):
    import yaml

    with base_cfg.open() as f:
        cfg = yaml.safe_load(f)

    cfg.setdefault("data", {})
    cfg["data"]["root"] = str(root)
    # pattern "" means scan all community subdirs; the builder will gate on presence
    cfg["data"]["pattern"] = cfg["data"].get("pattern", "")

    # start with preset defaults
    d = preset_defaults(preset, kind)

    # apply CLI overrides
    contigs_path = contigs_override or d["contigs_path"]
    gfa_candidates = (
        [p.strip() for p in gfa_override.split(",")] if gfa_override else d["gfa_candidates"]
    )
    save_pt_dir = (save_dir_override or d["save_pt_dir"]).replace("{kind}", kind)

    if not contigs_path:
        raise SystemExit(
            "[ERR] contigs path is required (use --contigs-path or choose a preset that defines it)"
        )

    cfg["data"]["contigs_path"] = contigs_path
    cfg["data"]["gfa_candidates"] = gfa_candidates
    cfg["data"]["save_pt_dir"] = save_pt_dir

    ensure_dir(out_yaml.parent)
    with out_yaml.open("w") as w:
        yaml.safe_dump(cfg, w, sort_keys=False)


def show_cuda_info():
    try:
        import torch  # type: ignore
        log("torch:", torch.__version__, "cuda:", getattr(torch.version, "cuda", None))
        log("CUDA available:", torch.cuda.is_available())
        if torch.cuda.is_available():
            try:
                log("Device:", torch.cuda.get_device_name(0))
            except Exception:
                pass
    except Exception as e:
        log("[INFO] PyTorch CUDA probe skipped:", e)


def build_one(
    repo_dir: Path,
    module: str,
    base_cfg: Path,
    art_dir: Path,
    tmp_cfg_dir: Path,
    split: str,
    root: Path,
    kind: str,
    preset: str,
    contigs_override: str | None,
    gfa_override: str | None,
    save_dir_override: str | None,
    force: bool,
) -> None:
    cfg_yaml = tmp_cfg_dir / f"build_{kind}_{split}.yaml"
    stamp = art_dir / f".build_{kind}_{split}_done"

    if stamp.exists() and not force:
        log(f"[SKIP] {kind} @ {split} already built ({stamp}).")
        return

    write_cfg(
        base_cfg=base_cfg,
        repo_dir=repo_dir,
        root=root,
        kind=kind,
        out_yaml=cfg_yaml,
        preset=preset,
        contigs_override=contigs_override,
        gfa_override=gfa_override,
        save_dir_override=save_dir_override,
    )

    log(f"[RUN ] Building graphs: split={split} kind={kind}")
    env = os.environ.copy()
    # ensure the repo is in PYTHONPATH so module imports work
    env["PYTHONPATH"] = str(repo_dir) + os.pathsep + env.get("PYTHONPATH", "")
    # thread hygiene if running on SLURM
    t = int(env.get("SLURM_CPUS_PER_TASK", env.get("OMP_NUM_THREADS", "8")))
    env["OMP_NUM_THREADS"] = str(t)
    env.setdefault("PYTORCH_CUDA_ALLOC_CONF", "max_split_size_mb:128,expandable_segments:True")

    cmd = [sys.executable, "-m", module, "--config", str(cfg_yaml)]
    rc = sp.call(cmd, cwd=str(repo_dir), env=env)
    if rc != 0:
        raise SystemExit(f"[ERR] Module failed with rc={rc}: {' '.join(cmd)}")

    # stamp success
    ensure_dir(art_dir)
    stamp.write_text("ok\n")
    log(f"[DONE] {kind} @ {split} built. Stamp: {stamp}")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Generic graph builder for GNN datasets (short/long/custom)."
    )
    ap.add_argument("--repo", type=Path, default=Path.cwd(), help="Path to repo root to run in")
    ap.add_argument("--base-config", type=Path, required=True, help="Base YAML config")
    ap.add_argument(
        "--module",
        type=str,
        default="data_prep.build_graphs_short",
        help="Python module to run (must accept --config YAML)",
    )

    ap.add_argument(
        "--splits",
        type=str,
        nargs="+",
        default=["train", "val"],
        help="Which splits to build (names only, e.g. train val test)",
    )
    ap.add_argument("--root-train", type=Path, help="Root of train communities")
    ap.add_argument("--root-val", type=Path, help="Root of val communities")
    ap.add_argument("--root-test", type=Path, help="Root of test communities")
    ap.add_argument(
        "--root", type=Path, help="Single root for when using just one split name (advanced)"
    )

    ap.add_argument(
        "--kinds", type=str, nargs="+", default=["generic", "filtered"], help="Kinds to build"
    )
    ap.add_argument(
        "--preset",
        choices=["short_metaspades", "long_flye", "custom"],
        default="short_metaspades",
        help="Choose defaults for contigs/GFA/save dirs; CLI overrides win",
    )
    ap.add_argument("--contigs-path", type=str, default=None,
                    help="Relative contigs path (per community); supports {kind}")
    ap.add_argument("--gfa-candidates", type=str, default=None,
                    help="Comma-separated relative paths to GFA/FASTG (per community); supports {kind}")
    ap.add_argument("--save-pt-dir", type=str, default=None,
                    help="Relative output dir for saved graphs; supports {kind}")

    ap.add_argument("--art-dir", type=Path, default=Path("runs/graphs"),
                    help="Artifacts directory under repo (stamps, tmp cfgs)")
    ap.add_argument("--force", action="store_true", help="Rebuild even if stamp exists")

    return ap.parse_args()


def main():
    args = parse_args()

    repo_dir = args.repo.resolve()
    base_cfg = args.base_config.resolve()
    if not base_cfg.exists():
        raise SystemExit(f"[ERR] Base config not found: {base_cfg}")

    # Ensure repo exists & is current working directory for module calls
    ensure_dir(repo_dir)
    # Make base config visible at repo/config.yaml (compat w/ any code expecting it there)
    cfg_link = repo_dir / "config.yaml"
    if not symlink_force(base_cfg, cfg_link):
        log("[WARN] Could not place config.yaml symlink/copy in repo; proceeding anyway.")

    # Show CUDA env info (optional)
    show_cuda_info()

    # Prepare artifact dirs
    art_dir = (repo_dir / args.art_dir).resolve()
    tmp_cfg_dir = art_dir / "_tmp_cfgs"
    ensure_dir(tmp_cfg_dir)

    # Map split -> root
    split_roots: Dict[str, Path] = {}
    for s in args.splits:
        key = f"root_{s}"
        root = getattr(args, key, None)
        if root is None and len(args.splits) == 1 and args.root:
            root = args.root
        if root is None:
            # Provide nice defaults for train/val/test if missing
            # only if the names match typical ones
            if s == "train" and args.root_train:
                root = args.root_train
            elif s == "val" and args.root_val:
                root = args.root_val
            elif s == "test" and args.root_test:
                root = args.root_test
        if root is None:
            raise SystemExit(f"[ERR] --root-{s} required for split '{s}'")
        root = Path(root).resolve()
        if not root.exists():
            raise SystemExit(f"[ERR] Split root not found for {s}: {root}")
        split_roots[s] = root

    # Build all requested (split, kind)
    for split, root in split_roots.items():
        for kind in args.kinds:
            build_one(
                repo_dir=repo_dir,
                module=args.module,
                base_cfg=base_cfg,
                art_dir=art_dir,
                tmp_cfg_dir=tmp_cfg_dir,
                split=split,
                root=root,
                kind=kind,
                preset=args.preset,
                contigs_override=(args.contigs_path.replace("{kind}", kind)
                                  if args.contigs_path else None),
                gfa_override=(args.gfa_candidates.replace("{kind}", kind)
                              if args.gfa_candidates else None),
                save_dir_override=(args.save_pt_dir.replace("{kind}", kind)
                                   if args.save_pt_dir else None),
                force=args.force,
            )

    log("[ALL DONE] Graph building complete.")


if __name__ == "__main__":
    main()
