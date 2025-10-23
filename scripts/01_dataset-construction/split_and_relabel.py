#!/usr/bin/env python3
"""
Split & Relabel to 4 Classes

Classes: prokaryote, eukaryote, plasmid, virus

Reads:
  <base_dir>/{train,val,test}/<category>/*.(fna|fa|fasta|*.gz)

Writes per-record FASTAs:
  <out_dir>/{train,val,test}/{prokaryote|eukaryote|plasmid|virus}/

Writes manifest (no stats):
  <meta_root>/metadata/<category>_refs4_manifest.csv
  where meta_root = out_dir/.. if out_dir name == "refs4", else out_dir

Flags:
  --delete_source    delete original file after successful split (>=1 record written)
  --trash_dir DIR    move original files here instead of deleting (overrides --delete_source)
  --dry_run          do not write or delete/move anything; just log actions

Usage:
  python split_and_relabel.py \
    --base_dir /path/to/my_dataset \
    --out_dir  /path/to/my_dataset/refs4 \
    --category bacteria \
    --delete_source
"""

import argparse
import csv
import gzip
import logging
import re
import shutil
from pathlib import Path
from typing import Dict, Iterable, List
from Bio import SeqIO

CATEGORY_TO_DOMAIN = {
    "archaea": "prokaryote",
    "bacteria": "prokaryote",
    "fungi": "eukaryote",
    "protozoa": "eukaryote",
    "virus": "virus",
}

def replicon_type(description: str) -> str:
    d = (description or "").lower()
    if "plasmid" in d:
        return "plasmid"
    if any(w in d for w in ["virus", "virion", "phage", "bacteriophage"]):
        return "virus"
    return "chromosomal"

def class4_for_record(category: str, description: str) -> str:
    rtype = replicon_type(description)
    if rtype == "plasmid":
        return "plasmid"
    if rtype == "virus":
        return "virus"
    dom = CATEGORY_TO_DOMAIN.get(category.lower())
    return dom if dom in {"prokaryote", "eukaryote", "virus"} else "unknown"

def sanitize_id(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s or "")

def iter_fasta_records(path: Path):
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

def open_manifest(manifest_path: Path):
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    fh = manifest_path.open("w", newline="")
    writer = csv.DictWriter(fh, fieldnames=["split","category","class4","accession","description","path","source_file"])
    writer.writeheader()
    return fh, writer

def process_file(
    fpath: Path,
    split: str,
    category: str,
    out_dir: Path,
    writer: csv.DictWriter,
    keep_unknown: bool,
    dry_run: bool,
) -> int:
    """
    Returns number of records successfully written from this source file.
    """
    written = 0
    try:
        for rec in iter_fasta_records(fpath):
            cls = class4_for_record(category, rec.description or rec.id)
            if cls == "unknown" and not keep_unknown:
                continue

            out_cls_dir = out_dir / split / cls
            if not dry_run:
                out_cls_dir.mkdir(parents=True, exist_ok=True)

            src_stem = fpath.stem
            rec_id = sanitize_id(rec.id)
            out_name = f"{src_stem}__{rec_id}.fna"
            out_path = out_cls_dir / out_name

            if dry_run:
                logging.info(f"[DRY] Would write {out_path}")
            else:
                SeqIO.write([rec], str(out_path), "fasta")

            writer.writerow({
                "split": split,
                "category": category,
                "class4": cls,
                "accession": rec.id,
                "description": rec.description,
                "path": str(out_path),
                "source_file": str(fpath),
            })
            written += 1
    except Exception as e:
        logging.exception(f"Failed parsing {fpath}: {e}")
    return written

def delete_or_trash(fpath: Path, trash_dir: Path | None, dry_run: bool):
    if trash_dir is not None:
        trash_dir.mkdir(parents=True, exist_ok=True)
        dst = trash_dir / fpath.name
        if dry_run:
            logging.info(f"[DRY] Would move source → {dst}")
        else:
            shutil.move(str(fpath), str(dst))
            logging.info(f"Moved source → {dst}")
    else:
        if dry_run:
            logging.info(f"[DRY] Would DELETE source {fpath}")
        else:
            fpath.unlink(missing_ok=True)
            logging.info(f"Deleted source {fpath}")

def process_split_category(
    base_dir: Path,
    out_dir: Path,
    split: str,
    category: str,
    writer: csv.DictWriter,
    keep_unknown: bool,
    delete_source: bool,
    trash_dir: Path | None,
    dry_run: bool,
) -> None:
    in_dir = base_dir / split / category
    if not in_dir.exists():
        logging.warning(f"[{split}/{category}] directory not found: {in_dir}")
        return

    fasta_files = find_fasta_files(in_dir)
    logging.info(f"[{split}/{category}] Found {len(fasta_files)} files")

    for fpath in fasta_files:
        logging.info(f"Processing {fpath}")
        n = process_file(fpath, split, category, out_dir, writer, keep_unknown, dry_run)
        if n > 0 and (delete_source or trash_dir is not None):
            delete_or_trash(fpath, trash_dir, dry_run)
        elif n == 0:
            logging.warning(f"No records written from {fpath}; skipping deletion.")

def main():
    ap = argparse.ArgumentParser(description="Split & relabel FASTA records into 4 classes; optionally delete originals.")
    ap.add_argument("--base_dir", required=True, type=Path, help="Base dir with train/val/test/<category>/*")
    ap.add_argument("--out_dir",  required=True, type=Path, help="Output dir for refs4/<split>/<class4>/")
    ap.add_argument("--category", required=True, help="archaea | bacteria | fungi | protozoa | virus")
    ap.add_argument("--keep_unknown", action="store_true", help="Keep records that cannot be classified")
    ap.add_argument("--delete_source", action="store_true", help="Delete original file after successful split")
    ap.add_argument("--trash_dir", type=Path, default=None, help="Move originals here instead of deleting")
    ap.add_argument("--dry_run", action="store_true", help="Log actions without writing or deleting")
    ap.add_argument("--log_level", default="INFO", choices=["DEBUG","INFO","WARNING","ERROR"])
    args = ap.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level), format="%(asctime)s [%(levelname)s] %(message)s")

    cat = args.category.lower()
    if cat not in {"archaea","bacteria","fungi","protozoa","virus"}:
        raise SystemExit(f"Unsupported category: {args.category}")

    # Decide manifest location
    meta_root = args.out_dir.parent if args.out_dir.name == "refs4" else args.out_dir
    manifest_path = meta_root / "metadata" / f"{cat}_refs4_manifest.csv"
    if args.dry_run:
        logging.info(f"[DRY] Would write manifest to {manifest_path}")

    # Open manifest (unless dry run)
    if args.dry_run:
        fh = None
        writer = None
    else:
        fh, writer = open_manifest(manifest_path)

    try:
        for split in ["train","val","test"]:
            if args.dry_run:
                # create a temporary no-op writer that just logs
                class _Dummy:
                    def writerow(self, row): logging.debug(f"[DRY] manifest row: {row}")
                writer = _Dummy()
            process_split_category(
                base_dir=args.base_dir,
                out_dir=args.out_dir,
                split=split,
                category=cat,
                writer=writer,
                keep_unknown=args.keep_unknown,
                delete_source=args.delete_source,
                trash_dir=args.trash_dir,
                dry_run=args.dry_run,
            )
    finally:
        if not args.dry_run and fh is not None:
            fh.close()
            logging.info(f"Wrote manifest: {manifest_path}")

if __name__ == "__main__":
    main()
