import subprocess
from pathlib import Path
import argparse
import time


def main(split, category):
    input_dir = Path(f"/path/to/4CAC_dataset/{split}/{category}/cleaned")
    output_dir = Path(f"/path/to/assembly_graphs/{split}/{category}")
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_dir.exists():
        print(f"[ERROR] Input directory not found: {input_dir}")
        return

    for fasta_file in input_dir.glob("*.cleaned.fna"):
        genome_id = fasta_file.stem.replace(".cleaned", "")
        genome_out_dir = output_dir / genome_id
        genome_out_dir.mkdir(parents=True, exist_ok=True)

        gfa_file = genome_out_dir / "assembly_graph_with_scaffolds.gfa"
        if gfa_file.exists():
            print(f"[SKIP] Already exists: {genome_id}")
            continue

        cmd = [
            "/path/to/venv/spades/bin/spades.py",
            "--meta",
            "--trusted-contigs", str(fasta_file),
            "-o", str(genome_out_dir),
            "--only-assembler"
        ]

        print(f"[INFO] Running metaSPAdes on: {genome_id}")
        start = time.time()
        try:
            subprocess.run(cmd, check=True)
            duration = time.time() - start
            print(f"[DONE] {genome_id} in {duration:.1f} seconds")
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] metaSPAdes failed for {genome_id}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run metaSPAdes on cleaned .fna contigs using --trusted-contigs.")
    parser.add_argument("--split", type=str, required=True, help="Data split (e.g., train, val, test)")
    parser.add_argument("--category", type=str, required=True, help="Genome category (e.g., bacteria, archaea)")
    args = parser.parse_args()

    main(args.split, args.category)
