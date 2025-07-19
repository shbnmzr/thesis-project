#!/usr/bin/env python3
"""
4CAC Pipeline Runner
--------------------
Runs the full 4CAC classification pipeline:
1. XGBoost classification
2. Graph-based refinement

Features:
- Validates inputs before execution.
- Logs to both console and file for reproducibility.
- Skips re-running if results already exist (unless forced).
"""

import argparse
import subprocess
import logging
from pathlib import Path
import sys
from scripts.helpers.logger import setup_logger

def ensure_trailing_slash(path: Path) -> str:
    """Ensure a trailing slash for paths required by 4CAC."""
    path_str = str(path)
    return path_str if path_str.endswith("/") else path_str + "/"

def run_command(cmd: list, step_name: str):
    logging.info(f"Running {step_name}: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logging.info(result.stdout)
    except subprocess.CalledProcessError as e:
        logging.error(f"{step_name} failed with return code {e.returncode}")
        logging.error(f"Standard Output:\n{e.stdout}")
        logging.error(f"Standard Error:\n{e.stderr}")
        sys.exit(1)

def validate_inputs(asmdir: Path):
    """Validate required metaSPAdes outputs."""
    required_files = ["scaffolds.fasta", "assembly_graph_with_scaffolds.gfa", "scaffolds.paths"]
    for f in required_files:
        file_path = asmdir / f
        if not file_path.exists():
            logging.error(f"Missing required file: {file_path}")
            sys.exit(1)

def skip_if_exists(asmdir: Path, force: bool) -> bool:
    """Skip execution if final 4CAC output already exists."""
    final_output = asmdir / "4CAC_classification.fasta"
    if final_output.exists() and not force:
        logging.info(f"Skipping 4CAC: {final_output} already exists. Use --force to re-run.")
        return True
    return False

def run_4cac(asmdir: Path, threads: int, python_exec: Path, tools_dir: Path, force: bool):
    """Run the full 4CAC pipeline."""
    validate_inputs(asmdir)
    if skip_if_exists(asmdir, force):
        return

    asmdir_str = ensure_trailing_slash(asmdir)
    contigs = asmdir / "scaffolds.fasta"
    xgb_out = asmdir / "scaffolds.fasta.probs_xgb_4class.out"

    # Step 1: XGBoost Classification
    classify_xgb_script = tools_dir / "classify_xgb.py"
    run_command([
        str(python_exec), str(classify_xgb_script),
        "-f", str(contigs),
        "-o", str(xgb_out),
        "-p", str(threads)
    ], "XGBoost Classification")

    # Step 2: Graph-based Refinement
    classify_4cac_script = tools_dir / "classify_4CAC.py"
    run_command([
        str(python_exec), str(classify_4cac_script),
        "--assembler", "metaSPAdes",
        "--asmdir", asmdir_str,
        "--classdir", asmdir_str,
        "--outdir", asmdir_str
    ], "Graph-based Refinement (4CAC)")

    logging.info(f"4CAC completed successfully. Final results in {asmdir}")

def main():
    parser = argparse.ArgumentParser(description="Run 4CAC (XGBoost + Graph refinement) on assembled contigs.")
    parser.add_argument("--asmdir", required=True, help="Directory containing metaSPAdes output.")
    parser.add_argument("--threads", type=int, default=8, help="Number of CPU threads.")
    parser.add_argument("--log_file", default="run_4cac.log", help="Path to the log file.")
    parser.add_argument("--python_exec", default="path/to/bin/python3.9", help="Python interpreter with 4CAC installed.")
    parser.add_argument("--tools_dir", default="path/to/4CAC", help="Path to the 4CAC code directory.")
    parser.add_argument("--force", action="store_true", help="Force re-run even if outputs exist.")
    args = parser.parse_args()

    setup_logger(Path(args.log_file))

    asmdir = Path(args.asmdir).resolve()
    python_exec = Path(args.python_exec)
    tools_dir = Path(args.tools_dir).resolve()

    logging.info("=== Starting 4CAC Pipeline ===")
    logging.info(f"Parameters:\nASMDIR={asmdir}\nTHREADS={args.threads}\nTOOLS={tools_dir}\nFORCE={args.force}")
    run_4cac(asmdir, args.threads, python_exec, tools_dir, args.force)
    logging.info("Pipeline finished successfully.")

if __name__ == "__main__":
    main()
