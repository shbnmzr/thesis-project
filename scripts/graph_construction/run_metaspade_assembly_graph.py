import subprocess
from pathlib import Path
import argparse
import shutil


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

        if (genome_out_dir / "assembly_graph.gfa").exists():
            print(f"[SKIP] Already exists: {genome_id}")
            continue

        # Prepare temp .fa copy
        renamed_fasta = genome_out_dir / f"{genome_id}.fa"
        genome_out_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(fasta_file, renamed_fasta)

        cmd = [
            "/path/to/venv/spades/bin/spades.py",
            "--sc",
            "-s", str(renamed_fasta),
            "-o", str(genome_out_dir),
            "--only-assembler"
        ]

        print(f"[INFO] Running metaSPAdes on: {genome_id}")
        try:
            subprocess.run(cmd, check=True)
            print(f"[DONE] {genome_id}")
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] metaSPAdes failed for {genome_id}: {e}")
        finally:
            # Always remove the temporary .fa file
            try:
                renamed_fasta.unlink()
                print(f"[CLEANUP] Removed temporary file: {renamed_fasta}")
            except Exception as e:
                print(f"[WARNING] Could not delete temp file {renamed_fasta}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run metaSPAdes on cleaned .fna files.")
    parser.add_argument("--split", type=str, required=True, help="Data split (e.g., train, val, test)")
    parser.add_argument("--category", type=str, required=True, help="Genome category (e.g., bacteria, archaea)")
    args = parser.parse_args()

    main(args.split, args.category)
