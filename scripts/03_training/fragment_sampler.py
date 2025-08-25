#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sample artificial contigs with:
  - uniform-per-genome allocation (equal target per input FASTA file)
  - within-record without-replacement start positions
Deterministic given --seed. Avoids Ns unless --allow_N is set.

Outputs:
  - FASTA of fragments
  - CSV labels: contig_id,class_name,genome_id,length
"""
import argparse
import gzip
from pathlib import Path
import random
import bisect

def fasta_record_lengths(path: Path):
    """Yield (record_id, length) for a .fna / .fa / .fasta (optionally .gz)."""
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
    """Yield (record_id, sequence_upper) for a .fna / .fa / .fasta (optionally .gz)."""
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        header, parts = None, []
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts).upper()
                header = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if header is not None:
            yield header, "".join(parts).upper()

def discover_records(genome_dirs, frag_len, min_len=300):
    """
    Discover records that can yield fragments and gather their capacities.
    Returns:
      recs: list of dicts {file, header, gid, length, avail}
    """
    recs = []
    for d in genome_dirs:
        for fp in sorted(Path(d).glob("*")):
            if not fp.is_file():
                continue
            name = fp.name.lower()
            if not (name.endswith(".fa") or name.endswith(".fna") or name.endswith(".fasta")
                    or name.endswith(".fa.gz") or name.endswith(".fna.gz") or name.endswith(".fasta.gz")):
                continue
            for header, L in fasta_record_lengths(fp):
                if L < max(min_len, frag_len):
                    continue
                avail = L - frag_len + 1
                gid = f"{fp.name}:{header}"  # provenance
                recs.append({"file": fp, "header": header, "gid": gid, "length": L, "avail": avail})
    return recs

def assign_counts_uniform_per_genome(recs, n_target, seed=42):
    """
    Uniform-per-genome (file) allocation with water-filling and per-record proportional split.
    Returns: list 'counts' of same length as recs (fragments per record).
    """
    rng = random.Random(seed)

    # Group record indices by genome (file)
    genomes = {}
    for i, r in enumerate(recs):
        genomes.setdefault(r["file"], []).append(i)

    # Filter genomes with nonzero total capacity
    g_info = []
    for fp, idxs in genomes.items():
        cap = sum(max(0, recs[i]["avail"]) for i in idxs)
        if cap > 0:
            g_info.append((fp, idxs, cap))

    if not g_info or n_target <= 0:
        return [0] * len(recs)

    G = len(g_info)
    base = n_target // G
    rem = n_target % G

    # Initial equal allocation per genome, deterministic remainder split
    alloc = {fp: base for fp, _, _ in g_info}
    for fp, _, _ in sorted(g_info, key=lambda t: str(t[0]))[:rem]:
        alloc[fp] += 1

    # Cap by capacity and track leftover
    leftover = 0
    for fp, idxs, cap in g_info:
        if alloc[fp] > cap:
            leftover += (alloc[fp] - cap)
            alloc[fp] = cap

    # Water-fill into genomes with spare capacity
    if leftover > 0:
        pool = []
        for fp, idxs, cap in g_info:
            spare = cap - alloc[fp]
            if spare > 0:
                pool.append([fp, spare])
        pool.sort(key=lambda x: str(x[0]))
        j = 0
        while leftover > 0 and pool:
            fp, spare = pool[j]
            if spare > 0:
                alloc[fp] += 1
                spare -= 1
                leftover -= 1
                pool[j][1] = spare
            j = (j + 1) % len(pool)
            # prune exhausted
            pool = [p for p in pool if p[1] > 0]
            if pool:
                j %= len(pool)

    # Split each genome's allocation across its records by record capacity
    counts = [0] * len(recs)
    for fp, idxs, cap in g_info:
        tgt = alloc[fp]
        if tgt <= 0:
            continue
        caps = [max(0, recs[i]["avail"]) for i in idxs]
        total_cap = sum(caps)
        if total_cap == 0:
            continue
        shares = [tgt * (c / total_cap) for c in caps]
        ints = [int(s) for s in shares]
        given = sum(ints)
        rems = [(s - i, j) for j, (s, i) in enumerate(zip(shares, ints))]
        need = tgt - given
        rems.sort(reverse=True)
        for k in range(need):
            ints[rems[k][1]] += 1
        # respect per-record capacity
        for local_j, i in enumerate(idxs):
            counts[i] = min(ints[local_j], recs[i]["avail"])
    return counts

def sample_from_sequence_without_replacement(seq, gid, L, k, class_name, fa, lab, rng, allow_N=False):
    """
    Emit up to k fragments WITHOUT replacement over start positions.
    """
    n = len(seq)
    start_max = n - L
    if start_max < 0 or k <= 0:
        return 0

    starts = list(range(0, start_max + 1))
    rng.shuffle(starts)

    written = 0
    for start in starts:
        if written >= k:
            break
        frag = seq[start:start+L]
        if allow_N or (("N" not in frag) and ("n" not in frag)):
            cid = f"{gid}|{L}|{start}"
            fa.write(f">{cid}\n")
            for i in range(0, L, 80):
                fa.write(frag[i:i+80] + "\n")
            lab.write(f"{cid},{class_name},{gid},{L}\n")
            written += 1
    return written

def main():
    ap = argparse.ArgumentParser(
        description="Sample artificial contigs: uniform-per-genome + without-replacement starts."
    )
    ap.add_argument("--genome_dirs", nargs="+", required=True,
                    help="One or more directories with genomes (*.fa/.fna/.fasta[.gz]).")
    ap.add_argument("--class_name", required=True,
                    choices=["Prokaryote","Plasmid","Virus","Microeukaryote"])
    ap.add_argument("--frag_len", type=int, required=True,
                    help="Fragment length in bp (e.g., 500, 1000, 5000, 10000, 50000).")
    ap.add_argument("--n_target", type=int, required=True,
                    help="Total number of fragments to sample.")
    ap.add_argument("--out_fasta", required=True,
                    help="Output FASTA file for fragments.")
    ap.add_argument("--out_labels", required=True,
                    help="Output CSV (contig_id,class_name,genome_id,length).")
    ap.add_argument("--seed", type=int, default=42, help="Random seed.")
    ap.add_argument("--allow_N", action="store_true", help="Allow Ns in sampled fragments.")
    args = ap.parse_args()

    recs = discover_records(args.genome_dirs, args.frag_len)
    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_labels).parent.mkdir(parents=True, exist_ok=True)

    if not recs:
        with open(args.out_fasta, "w") as _fa, open(args.out_labels, "w") as lab:
            lab.write("contig_id,class_name,genome_id,length\n")
        print("[WARN] No suitable records found. Wrote 0 fragments.")
        return

    counts = assign_counts_uniform_per_genome(recs, args.n_target, seed=args.seed)
    rng = random.Random(args.seed)

    written_total = 0
    with open(args.out_fasta, "w") as fa, open(args.out_labels, "w") as lab:
        lab.write("contig_id,class_name,genome_id,length\n")

        # Build per-file intention map for deterministic iteration
        need = {}
        for rec, c in zip(recs, counts):
            if c > 0:
                need.setdefault(rec["file"], {})[rec["header"]] = (rec["gid"], c)

        for fp in sorted(need.keys(), key=lambda p: str(p)):
            inner = need[fp]
            for header, seq in iter_fasta_records(fp):
                entry = inner.get(header)
                if entry is None:
                    continue
                gid, k = entry
                written_total += sample_from_sequence_without_replacement(
                    seq, gid, args.frag_len, k, args.class_name, fa, lab, rng,
                    allow_N=args.allow_N
                )

    target_sum = sum(counts)
    if written_total < target_sum:
        print(f"[INFO] Wrote {written_total} fragments (targetable windows {target_sum}); "
              f"shortfall likely due to N filtering.")
    else:
        print(f"Wrote {written_total} fragments: class={args.class_name} len={args.frag_len} → {args.out_fasta}")

if __name__ == "__main__":
    main()
