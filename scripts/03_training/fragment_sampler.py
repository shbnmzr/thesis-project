#!/usr/bin/env python3
import argparse
import gzip
import os
import bisect
from pathlib import Path
import random

def fasta_record_lengths(path: Path):
    """
    Stream only lengths of FASTA records from .fna or .fna.gz without storing sequences.
    Yields (record_id, length).
    """
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        header = None
        length = 0
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, length
                header = line[1:].strip().split()[0]
                length = 0
            else:
                length += len(line.strip())
        if header is not None:
            yield header, length

def iter_fasta_records(path: Path):
    """
    Stream full FASTA records from .fna or .fna.gz.
    Yields (record_id, sequence_string_upper).
    """
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        header, seq_parts = None, []
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts).upper()
                header = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            yield header, "".join(seq_parts).upper()

def discover_records(genome_dirs, frag_len, min_len=300):
    """
    First pass: discover records that can yield fragments and gather their weights.
    Returns:
      recs: list of dicts with {file, header, gid, length, avail}
      cum_weights: cumulative sum of avail for bisect
      total_avail: total available starts across all records
    """
    recs = []
    cum_weights = []
    total = 0

    for d in genome_dirs:
        for fp in sorted(Path(d).glob("*.fna*")):
            for header, L in fasta_record_lengths(fp):
                if L < max(min_len, frag_len):
                    continue
                avail = L - frag_len + 1
                gid = f"{fp.name}:{header}"  # provenance-compatible ID
                recs.append({"file": fp, "header": header, "gid": gid, "length": L, "avail": avail})
                total += avail
                cum_weights.append(total)

    return recs, cum_weights, total

def assign_counts(recs, cum_weights, total_avail, n_target, seed=42):
    """
    Assign how many fragments each record should produce such that each possible
    start across the dataset is equally likely (exactly matches global uniform over starts).
    Returns a list 'counts' the same length as recs.
    """
    import random
    rng = random.Random(seed)
    m = len(recs)
    counts = [0] * m
    if total_avail == 0 or n_target <= 0:
        return counts
    for _ in range(n_target):
        r = rng.randrange(total_avail)  # [0, total_avail-1]
        j = bisect.bisect_right(cum_weights, r)
        if j >= m:
            j = m - 1
        counts[j] += 1
    return counts

def sample_from_sequence(seq, gid, L, k, class_name, fa, lab, rng, allow_N=False, retries=10):
    """
    For a single sequence string 'seq', emit k fragments of length L.
    Avoid Ns by retrying up to 'retries' times unless allow_N=True.
    Write to open FASTA 'fa' and labels CSV 'lab'.
    """
    n = len(seq)
    start_max = n - L
    if start_max < 0 or k <= 0:
        return 0
    written = 0
    for _ in range(k):
        ok = False
        for __ in range(retries):
            start = rng.randint(0, start_max)
            frag = seq[start:start+L]
            if allow_N or (("N" not in frag) and ("n" not in frag)):
                ok = True
                break
        if not ok:
            # Skip this fragment if we couldn't find a clean window in 'retries' attempts
            continue
        cid = f"{gid}|{L}|{start}"
        fa.write(f">{cid}\n")
        for i in range(0, L, 80):
            fa.write(frag[i:i+80] + "\n")
        lab.write(f"{cid},{class_name},{gid},{L}\n")
        written += 1
    return written

def main():
    ap = argparse.ArgumentParser(
        description="Sample artificial contigs (fragments) from genomes for a given class."
    )
    ap.add_argument("--genome_dirs", nargs="+", required=True,
                    help="One or more directories with *.fna or *.fna.gz genomes.")
    ap.add_argument("--class_name", required=True,
                    choices=["Prokaryote","Plasmid","Virus","Microeukaryote"])
    ap.add_argument("--frag_len", type=int, required=True,
                    help="Fragment length in bp (e.g., 500, 1000, 5000, 10000, 50000).")
    ap.add_argument("--n_target", type=int, required=True,
                    help="Number of fragments to sample.")
    ap.add_argument("--out_fasta", required=True,
                    help="Output FASTA file for fragments.")
    ap.add_argument("--out_labels", required=True,
                    help="Output CSV with columns: contig_id,class_name,genome_id,length")
    ap.add_argument("--seed", type=int, default=42, help="Random seed.")
    ap.add_argument("--allow_N", action="store_true", help="Allow Ns in sampled fragments.")
    ap.add_argument("--max_retries_per_fragment", type=int, default=10,
                    help="Retries per fragment to avoid Ns.")
    args = ap.parse_args()

    # Pass A: discover records and build weights
    recs, cum_weights, total_avail = discover_records(args.genome_dirs, args.frag_len)
    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_labels).parent.mkdir(parents=True, exist_ok=True)

    # If nothing is long enough, write empty CSV header and exit gracefully
    if total_avail == 0:
        with open(args.out_fasta, "w") as _fa, open(args.out_labels, "w") as lab:
            lab.write("contig_id,class_name,genome_id,length\n")
        print(f"[WARN] No sequences long enough for frag_len={args.frag_len}. Wrote 0 fragments.")
        return

    # Assign how many fragments to draw from each record (exact uniform over starts)
    counts = assign_counts(recs, cum_weights, total_avail, args.n_target, seed=args.seed)

    # Pass B: re-read sequences one-by-one and emit assigned fragments
    rng = random.Random(args.seed)
    written_total = 0
    with open(args.out_fasta, "w") as fa, open(args.out_labels, "w") as lab:
        lab.write("contig_id,class_name,genome_id,length\n")

        # Build a lookup of how many to draw per (file, header)
        need = {}
        for rec, c in zip(recs, counts):
            if c > 0:
                need.setdefault(rec["file"], {})[rec["header"]] = (rec["gid"], c)

        # Iterate files and stream records
        for fp, inner in need.items():
            for header, seq in iter_fasta_records(fp):
                entry = inner.get(header)
                if entry is None:
                    continue
                gid, k = entry
                written_total += sample_from_sequence(
                    seq, gid, args.frag_len, k, args.class_name, fa, lab,
                    rng, allow_N=args.allow_N, retries=args.max_retries_per_fragment
                )

    print(f"Wrote {written_total} fragments: class={args.class_name} len={args.frag_len} → {args.out_fasta}")

if __name__ == "__main__":
    main()
