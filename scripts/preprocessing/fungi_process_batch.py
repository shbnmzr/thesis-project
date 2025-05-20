import argparse
from Bio import SeqIO
import subprocess
from pathlib import Path
import csv
import sys
import time
import builtins
from concurrent.futures import ThreadPoolExecutor, as_completed

# --- Configurable Parameters ---
CATEGORY = "fungi"
BASE_DIR = Path("../4CAC_dataset")
GENOME_DIR = BASE_DIR / "genomes" / CATEGORY
LOG_DIR = BASE_DIR / "logs"
MIN_LEN = 1000
MAX_AMBIG = 0.05
CDHIT_ID = 0.99
MIN_CONTIGS = 10
CDHIT_THREADS = 4   # Each CD-HIT instance uses 3 threads
CDHIT_MEMORY_MB = 8000
MAX_WORKERS = 4  # Number of parallel threads per job (adjust based on --cpus-per-task)

# --- Argument Parsing ---
parser = argparse.ArgumentParser()
parser.add_argument("--start", type=int, required=True, help="Start index of genome slice")
parser.add_argument("--end", type=int, required=True, help="End index of genome slice (exclusive)")
args = parser.parse_args()

# --- Logging ---
LOG_DIR.mkdir(parents=True, exist_ok=True)
batch_log = LOG_DIR / f"{CATEGORY}_batch_{args.start}_{args.end}.log"
log_stream = open(batch_log, "w")

def dual_print(*args, **kwargs):
    builtins.print(*args, **kwargs, flush=True)
    builtins.print(*args, **kwargs, file=log_stream, flush=True)

print = dual_print

# --- Validate and Collect Files ---
if not GENOME_DIR.exists():
    print(f"[ERROR] Genome directory not found: {GENOME_DIR.resolve()}")
    sys.exit(1)

all_files = sorted(list(GENOME_DIR.glob("*.fna")))
slice_files = all_files[args.start:args.end]

print(f"[INFO] Batch range: {args.start} to {args.end} ({len(slice_files)} genomes)\n")
if not slice_files:
    print(f"[WARNING] No FASTA files found in this slice.")
    sys.exit(0)

def clean_fasta(fasta_path):
    seen = set()
    records = []
    total_before = 0
    try:
        for rec in SeqIO.parse(fasta_path, "fasta"):
            total_before += 1
            seq = str(rec.seq).upper()
            if len(seq) == 0 or seq in seen:
                continue
            seen.add(seq)
            if len(seq) < MIN_LEN or (seq.count("N") / len(seq)) > MAX_AMBIG:
                continue
            records.append(rec)
    except Exception as e:
        print(f"[ERROR] Failed to parse {fasta_path.name}: {e}")
        return None, 0, 0

    if records:
        cleaned_path = fasta_path.with_suffix(".cleaned.fna")
        SeqIO.write(records, cleaned_path, "fasta")
        return cleaned_path, total_before, len(records)
    else:
        fasta_path.unlink()
        return None, total_before, 0

def run_cd_hit(input_path, output_path):
    stdout_path = output_path.with_suffix(".cdhit_stdout.log")
    stderr_path = output_path.with_suffix(".cdhit_stderr.log")
    cmd = [
        "cd-hit",
        "-i", str(input_path),
        "-o", str(output_path),
        "-c", str(CDHIT_ID),
        "-n", "5",
        "-M", str(CDHIT_MEMORY_MB),
        "-T", str(CDHIT_THREADS)
    ]
    with open(stdout_path, "w") as out, open(stderr_path, "w") as err:
        result = subprocess.run(cmd, stdout=out, stderr=err)
    input_path.unlink()
    return result.returncode == 0

def count_contigs(path):
    return sum(1 for _ in SeqIO.parse(path, "fasta"))

def get_gc_avg_len(path):
    seqs = list(SeqIO.parse(path, "fasta"))
    if not seqs:
        return 0.0, 0.0
    total_gc = sum((str(s.seq).upper().count("G") + str(s.seq).upper().count("C")) for s in seqs)
    total_len = sum(len(s.seq) for s in seqs)
    return 100 * total_gc / total_len, total_len / len(seqs)

def process_file(fasta_file):
    print(f"[INFO] --> {fasta_file.name}")
    cleaned, before, after_filter = clean_fasta(fasta_file)

    if cleaned:
        dedup = fasta_file.with_suffix(".dedup.fna")
        if run_cd_hit(cleaned, dedup) and dedup.exists():
            after_cd = count_contigs(dedup)
            if after_cd < MIN_CONTIGS:
                dedup.unlink()
                gc, avg_len, status = 0.0, 0.0, "DELETED (too few contigs)"
            else:
                gc, avg_len = get_gc_avg_len(dedup)
                status = "KEPT"
        else:
            after_cd, gc, avg_len, status = 0, 0.0, 0.0, "FAILED (CD-HIT)"
    else:
        after_cd, gc, avg_len, status = 0, 0.0, 0.0, "DELETED (empty after filtering)"

    return (fasta_file.name, before, after_filter, after_cd, f"{gc:.2f}", f"{avg_len:.1f}", status)

# --- Multithreaded Processing ---
report = []
with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
    futures = {executor.submit(process_file, f): f for f in slice_files}
    for future in as_completed(futures):
        try:
            report.append(future.result())
        except Exception as e:
            print(f"[ERROR] Exception while processing {futures[future].name}: {e}")

# --- Save CSV Report ---
csv_path = LOG_DIR / f"{CATEGORY}_batch_{args.start}_{args.end}_report.csv"
with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Filename", "Original Contigs", "After Filtering", "After CD-HIT", "GC %", "Avg Length", "Status"])
    writer.writerows(report)

print(f"\n[INFO] Batch complete. Report saved to {csv_path}")
time.sleep(2)