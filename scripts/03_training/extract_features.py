#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Streaming k-mer feature extraction:
  - Canonical k-mers for k=3..7, concatenated (dim = 10,952)
  - Per-k L1 normalization
  - Adds length and GC% (GC over ACGT only)
  - Writes TSV or TSV.GZ based on --out extension
"""
import argparse, gzip
from pathlib import Path
import numpy as np

COMP = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(s): return s.translate(COMP)[::-1]
def canonical(s):
    rc = revcomp(s)
    return s if s <= rc else rc

def canonical_kmer_list(k):
    import itertools
    seen = set()
    for t in itertools.product("ACGT", repeat=k):
        s = "".join(t)
        seen.add(canonical(s))
    return sorted(seen)

# Precompute once
_KS = [3,4,5,6,7]
_KMER_LISTS = {k: canonical_kmer_list(k) for k in _KS}
_KIDX = {k: {km:i for i,km in enumerate(_KMER_LISTS[k])} for k in _KS}
_KCOLS = [f"k{k}:{km}" for k in _KS for km in _KMER_LISTS[k]]
_EXPECTED_DIM = sum(len(_KMER_LISTS[k]) for k in _KS)  # 10952
assert len(_KCOLS) == _EXPECTED_DIM, "Column list length mismatch"

def iter_fasta(p: Path):
    op = gzip.open if str(p).endswith(".gz") else open
    with op(p, "rt") as f:
        hdr, parts = None, []
        for line in f:
            if line.startswith(">"):
                if hdr is not None:
                    seq = "".join(parts).upper()
                    yield hdr, seq
                hdr = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if hdr is not None:
            seq = "".join(parts).upper()
            yield hdr, seq

def is_ambig_window(s):
    for c in s:
        if c not in "ACGT":
            return True
    return False

def kfreq_block(seq, k, idx):
    n = len(seq)
    v = np.zeros(len(idx), dtype=np.float32)
    if n < k:  # no windows
        return v
    valid = 0
    for i in range(n - k + 1):
        s = seq[i:i+k]
        if is_ambig_window(s):
            continue
        v[idx[canonical(s)]] += 1.0
        valid += 1
    if valid > 0:
        v /= float(valid)
    return v

def features_3to7(seq):
    blocks = [kfreq_block(seq, k, _KIDX[k]) for k in _KS]
    return np.concatenate(blocks).astype(np.float32, copy=False)

def gc_and_length(seq):
    # Compute GC over A/C/G/T only; length is full sequence length
    L = len(seq)
    a = seq.count("A"); c = seq.count("C"); g = seq.count("G"); t = seq.count("T")
    acgt = a + c + g + t
    gc = (c + g) / acgt if acgt > 0 else 0.0
    return L, gc

def open_maybe_gzip(out_path: Path):
    return gzip.open(out_path, "wt") if str(out_path).endswith(".gz") else open(out_path, "wt")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch-size", type=int, default=5000)
    args = ap.parse_args()

    outp = Path(args.out)
    wrote = 0
    with open_maybe_gzip(outp) as out:
        # Add length and gc right after contig_id
        out.write("contig_id\tlength\tgc\t" + "\t".join(_KCOLS) + "\n")

        buf_ids, buf_lens, buf_gcs, buf_rows = [], [], [], []
        for cid, seq in iter_fasta(Path(args.fasta)):
            L, GC = gc_and_length(seq)
            buf_ids.append(cid)
            buf_lens.append(L)
            buf_gcs.append(GC)
            buf_rows.append(features_3to7(seq))
            if len(buf_ids) >= args.batch_size:
                for i, row in enumerate(buf_rows):
                    out.write(
                        f"{buf_ids[i]}\t{buf_lens[i]}\t{buf_gcs[i]:.6f}\t" +
                        "\t".join(f"{x:.7g}" for x in row.tolist()) + "\n"
                    )
                wrote += len(buf_ids)
                buf_ids.clear(); buf_lens.clear(); buf_gcs.clear(); buf_rows.clear()

        # flush
        for i, row in enumerate(buf_rows):
            out.write(
                f"{buf_ids[i]}\t{buf_lens[i]}\t{buf_gcs[i]:.6f}\t" +
                "\t".join(f"{x:.7g}" for x in row.tolist()) + "\n"
            )
        wrote += len(buf_ids)

    print(f"Wrote {wrote} rows, {2 + len(_KCOLS)} feature cols → {outp}")

if __name__ == "__main__":
    main()
