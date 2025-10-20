#!/usr/bin/env python3
"""
Split & Relabel to 4 Classes

Classes: prokaryote, eukaryote, plasmid, virus

Reads from:
  <base_dir>/{train,val,test}/<category>/*.(fna|fa|fasta|*.gz)

Writes per-record FASTAs to:
  <out_dir>/{train,val,test}/{prokaryote|eukaryote|plasmid|virus}/

Writes manifest to:
  <out_dir>/../metadata/<category>_refs4_manifest.csv

Usage:
  python split_and_relabel.py \
    --base_dir /path/to/my_dataset \
    --out_dir  /path/to/my_dataset/refs4 \
    --category bacteria
"""

import argparse
import csv
import gzip
import logging
import re
from pathlib import Path
from typing import List, Dict
from Bio import SeqIO

# Map coarse categories to domain classes for chromosomal records
CATEGORY_TO_DOMAIN = {
    "archaea": "prokaryote",
    "bacteria": "prokaryote",
    "fungi": "eukaryote",
    "protozoa": "eukaryote",
    "virus": "virus",
}

def replicon_type(description: str) -> str:
    """Return 'plasmid', 'virus', or 'chromosomal' from description."""
    d = (description or "").lower()
    if "plasmid" in d:
        return "plasmid"
    if any(w in d for w in ["virus", "virion", "phage", "bacteriophage"]):
        return "virus"
    return "chromosomal"

def class4_for_record(category: str, description: str) -> str:
    """Map a record to final class among {'prokaryote','eukaryote','plasmid','virus','unknown'}."""
    rtype = replicon_type(description)
    if rtype == "plasmid":
        return "plasmid"
    if rtype == "virus":
        return "virus"
    # chromosomal/scaffold → map via category
    dom = CATEGORY_TO_DOMAIN.get(category.lower())
    return dom if dom in {"prokaryote", "eukaryote", "virus"} else "unknown"

def sanitize_id(s: str) -> str:
    """Make a filesystem-safe identifier for filenames."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s or "")

def iter_fasta_records(path: Path):
    """Yield SeqRecords from FASTA or FASTA.GZ."""
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
            yield from SeqIO.parse(fh, "fasta")
    else:
        with open(path, "rt", encoding="utf-8", errors="ignore") as fh:
            yield from SeqIO.parse(fh, "fasta")

def find_fasta_files(root: Path) -> List[Path]:
    exts = ["*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fa.gz", "*.fasta.gz"]
    files: List[Path] = []
    for ext in exts:
        files.extend(root.glob(ext))
    return sorted(files)

def process_split_category(
    base_dir: Path,
    out_dir: Path,
    split: str,
    category: str,
    keep_unknown: bool = False,
) -> List[Dict]:
    """
    For one split/category, write per-record FASTAs into out_dir/<split>/<class4>/.
    Returns manifest rows.
    """
    in_dir = base_dir / split / category
    manifest_rows: List[Dict] = []

    if not in_dir.exists():
        logging.warning(f"[{split}/{category}] directory not found: {in_dir}")
        return manifest_rows

    fasta_files = find_fasta_files(in_dir)
    logging.info(f"[{split}/{category}] Found {len(fasta_files)} files")

    for fpath in fasta_files:
        try:
            for rec in iter_fasta_records(fpath):
                cls = class4_for_record(category, rec.description or rec.id)
                if cls == "unknown" and not keep_unknown:
                    continue

                out_cls_dir = out_dir / split / cls
                out_cls_dir.mkdir(parents=True, exist_ok=True)

                # Build deterministic filename per record
                src_stem = fpath.stem  # strips one suffix (works for .gz too)
                rec_id = sanitize_id(rec.id)
                out_name = f"{src_stem}__{rec_id}.fna"
                out_path = out_cls_dir / out_name

                # Write single-record FASTA
                SeqIO.write([rec], str(out_path), "fasta")

                manifest_rows.append({
                    "split": split,
                    "category": category,
                    "class4": cls,
                    "accession": rec.id,
                    "description": rec.description,
                    "path": str(out_path),
                    "source_file": str(fpath),
                })
        except Exception as e:
            logging.exception(f"Failed parsing {fpath}: {e}")

    return manifest_rows

def write_manifest(rows: List[Dict], manifest_path: Path):
    if not rows:
        logging.warning(f"No rows to write in manifest: {manifest_path}")
        return
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["split","category","class4","accession","description","path","source_file"]
    with manifest_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    logging.info(f"Wrote manifest: {manifest_path} (n={len(rows)})")

def main():
    ap = argparse.ArgumentParser(description="Split & relabel FASTA records into 4 classes.")
    ap.add_argument("--base_dir", required=True, type=Path, help="Base dir with train/val/test/<category>/*")
    ap.add_argument("--out_dir",  required=True, type=Path, help="Output dir for refs4/<split>/<class4>/")
    ap.add_argument("--category", required=True, help="archaea | bacteria | fungi | protozoa | virus")
    ap.add_argument("--keep_unknown", action="store_true", help="Keep records that cannot be classified")
    ap.add_argument("--log_level", default="INFO", choices=["DEBUG","INFO","WARNING","ERROR"])
    args = ap.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s [%(levelname)s] %(message)s")

    cat = args.category.lower()
    if cat not in {"archaea","bacteria","fungi","protozoa","virus"}:
        raise SystemExit(f"Unsupported category: {args.category}")

    all_rows: List[Dict] = []
    for split in ["train","val","test"]:
        all_rows.extend(process_split_category(args.base_dir, args.out_dir, split, cat, keep_unknown=args.keep_unknown))

    # Manifest lives under ../metadata if out_dir ends with refs4; otherwise under out_dir/metadata
    default_meta_root = args.out_dir.parent if args.out_dir.name == "refs4" else args.out_dir
    manifest_path = default_meta_root / "metadata" / f"{cat}_refs4_manifest.csv"
    write_manifest(all_rows, manifest_path)

if __name__ == "__main__":
    main()
