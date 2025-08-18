#!/usr/bin/env python3
import argparse, gzip, numpy as np, pandas as pd
from pathlib import Path

DNA = "ACGT"
COMP = str.maketrans("ACGT","TGCA")

def revcomp(s): return s.translate(COMP)[::-1]
def is_ambig(s): return any(c not in "ACGT" for c in s)
def canonical(s):
    rc = revcomp(s)
    return s if s <= rc else rc

def canonical_kmer_list(k):
    import itertools
    seen = set()
    for t in itertools.product(DNA, repeat=k):
        s = "".join(t)
        seen.add(canonical(s))
    return sorted(seen)

# Precompute canonical index per k
_KIDX = {k: {km:i for i,km in enumerate(canonical_kmer_list(k))} for k in [3,4,5,6,7]}
_KCOLS = [f"k{k}:{km}" for k in [3,4,5,6,7] for km in canonical_kmer_list(k)]

def iter_fasta(p: Path):
    op = gzip.open if str(p).endswith(".gz") else open
    with op(p, "rt") as f:
        hdr, parts = None, []
        for line in f:
            if line.startswith(">"):
                if hdr: yield hdr, "".join(parts).upper()
                hdr = line[1:].strip().split()[0]; parts = []
            else:
                parts.append(line.strip())
        if hdr: yield hdr, "".join(parts).upper()

def kfreq_block(seq, k, idx):
    n=len(seq)
    if n<k:
        return np.zeros(len(idx), dtype=float)
    v=np.zeros(len(idx), dtype=float)
    valid=0
    for i in range(n-k+1):
        s=seq[i:i+k]
        if is_ambig(s): continue
        v[idx[canonical(s)]]+=1.0; valid+=1
    if valid>0: v/=float(valid)
    return v

def features_3to7(seq):
    return np.concatenate([kfreq_block(seq,k,_KIDX[k]) for k in [3,4,5,6,7]])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ids, rows = [], []
    for cid, seq in iter_fasta(Path(args.fasta)):
        ids.append(cid)
        rows.append(features_3to7(seq))
    X = np.vstack(rows) if rows else np.zeros((0,10952), dtype=float)
    if X.shape[1] != 10952:
        raise SystemExit(f"Expected 10,952 features, got {X.shape[1]}")
    df = pd.DataFrame(X, columns=_KCOLS)
    df.insert(0, "contig_id", ids)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {df.shape[0]} rows, {df.shape[1]-1} features → {args.out}")

if __name__ == "__main__":
    main()
