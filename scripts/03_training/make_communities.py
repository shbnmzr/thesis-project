"""
Create synthetic training communities by sampling genomes from filesystem roots and
emitting (a) per-class genome selection lists (*.sel) and (b) per-genome abundance
TSVs built from Dirichlet-distributed weights, scaled by per-class abundance ratios.

Classes are:
  - prokaryote
  - microeukaryote
  - virus
  - plasmid

For each requested community (short/long), we:
  1) Randomly sample K genomes per class, where K is drawn uniformly from a per-class
     [lo, hi] range.
  2) Write the selected genome file paths to {outroot}/{community}/genomes/{class}.sel
  3) For each target abundance profile ("generic" and "filtered"), draw per-genome
     Dirichlet weights within each class and scale by the class ratio, then write a
     two-column TSV: <genome_path>\t<weight> in {outroot}/{community}/abundance/*.tsv

Notes
-----
* Files are discovered by recursively globbing under the provided root directories
  with the given filename patterns (e.g., *.fna, *.fa.gz, ...).
* Randomness is driven by --seed for reproducibility.
* If a class lacks sufficient genomes to satisfy its lower bound, the program exits
  with an error.
"""

from pathlib import Path
import argparse, os, sys, random, glob, json

# Fixed, known class labels used throughout the pipeline
CLASSES = ["prokaryote","microeukaryote","virus","plasmid"]


def list_fastas(roots, patterns):
    """Recursively list FASTA-like files under *roots* matching any of *patterns*.

    Parameters
    ----------
    roots : Iterable[str | Path]
        One or more directory roots to search.
    patterns : Iterable[str]
        Glob patterns (e.g., "*.fna", "*.fa.gz"). Each pattern is searched
        recursively from each root.

    Returns
    -------
    list[str]
        Sorted list of matching file paths (as strings) with duplicates removed.
    """
    files = set()
    for root in roots:
        p = Path(root)
        if not p.exists():
            # Silently skip non-existent roots, allowing optional inputs
            continue
        for pat in patterns:
            # Use ** to recurse; glob.glob requires a string path
            files.update(glob.glob(str(p / "**" / pat), recursive=True))
    return sorted(files)


def sample_genomes(files, k, rng):
    """Randomly sample *k* genome paths from *files* without replacement.

    Shuffles a copy of *files* in-place using the provided *rng* (random.Random)
    and returns the first *k* entries.
    """
    files = list(files)
    rng.shuffle(files)
    return files[:k]


def dirichlet_weights(n, rng, alpha=0.3, floor_frac=0.01):
    """Draw a length-*n* Dirichlet-like weight vector with optional flooring.

    Implementation detail: Instead of calling a Dirichlet sampler, we draw *n*
    independent Gamma(alpha, 1.0) variates and normalize them. Optionally, we
    floor each weight to at least *floor_frac/n* and renormalize to avoid
    extremely small values.

    Parameters
    ----------
    n : int
        Number of weights to generate.
    rng : random.Random
        RNG used for reproducible draws.
    alpha : float, default 0.3
        Concentration parameter; smaller values yield sparser weights.
    floor_frac : float, default 0.01
        Minimum total mass guaranteed across all items (distributed evenly).

    Returns
    -------
    list[float]
        A list of *n* non-negative weights summing to 1.0 (within FP error).
    """
    # Draw Gamma(alpha, 1.0) variates
    g = [rng.gammavariate(alpha, 1.0) for _ in range(n)]
    s = sum(g) or 1.0
    w = [x / s for x in g]

    if floor_frac and floor_frac > 0:
        # Apply a per-item floor, then renormalize
        floor = floor_frac / float(n)
        w = [max(x, floor) for x in w]
        s = sum(w)
        w = [x / s for x in w]
    return w


def write_sel(outdir, picks):
    """Write per-class genome selections (*.sel) under *outdir*.

    Parameters
    ----------
    outdir : Path
        Directory to create and write selection files into.
    picks : dict[str, list[str]]
        Mapping class -> list of selected genome file paths (strings).
    """
    outdir.mkdir(parents=True, exist_ok=True)
    for cls, lst in picks.items():
        with open(outdir / f"{cls}.sel", "w") as f:
            for p in lst:
                f.write(p + "\n")


def write_per_genome_tsv(out_tsv, picks, class_ratio, rng, alpha=0.3, floor_frac=0.01):
    """Write per-genome abundance weights to a TSV.

    For each class in *picks*, generate Dirichlet weights across that class's
    selected genomes, then scale each per-genome weight by the class-level
    abundance ratio in *class_ratio*.

    The final file is a two-column TSV: <path>\t<weight> (floating).
    """
    rows = []
    for cls, files in picks.items():
        n = len(files)
        if n == 0:
            # No genomes were selected for this class; skip
            continue
        # Dirichlet weights within the class
        w = dirichlet_weights(n, rng, alpha=alpha, floor_frac=floor_frac)
        # Scale by class ratio and enqueue rows
        for path, wi in zip(files, w):
            rows.append((path, class_ratio[cls] * float(wi)))

    # Ensure parent directory exists and emit TSV lines with reasonable precision
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w") as f:
        for path, w in rows:
            f.write(f"{path}\t{w:.10g}\n")


def make_one(comm_name, outroot, files_by_cls, ranges, generic, filtered, rng):
    """Assemble and write one community under {outroot}/{comm_name}.

    Steps:
      1) Draw K per-class from *ranges[cls]* and sample that many genomes
      2) Write selection lists (*.sel)
      3) Write abundance TSVs for both the "generic" and "filtered" class ratios
    """
    base = outroot / comm_name
    (base / "genomes").mkdir(parents=True, exist_ok=True)

    # Per-class genome selections
    picks = {}
    for cls in CLASSES:
        lo, hi = ranges[cls]
        if len(files_by_cls[cls]) < lo:
            sys.exit(
                f"[ERR] Not enough genomes for {cls}: need >= {lo}, found {len(files_by_cls[cls])}"
            )
        # Uniform integer draw in [lo, hi]
        k = rng.randint(lo, hi)
        picks[cls] = sample_genomes(files_by_cls[cls], k, rng)

    # Write per-class selection lists
    write_sel(base / "genomes", picks)

    # Emit per-genome abundance files for both ratio presets
    abr = base / "abundance"
    for target, ratio in [
        ("prokaryote_generic.tsv", generic),
        ("microeukaryote_generic.tsv", generic),
        ("virus_generic.tsv", generic),
        ("plasmid_generic.tsv", generic),
        ("prokaryote_filtered.tsv", filtered),
        ("microeukaryote_filtered.tsv", filtered),
        ("virus_filtered.tsv", filtered),
        ("plasmid_filtered.tsv", filtered),
    ]:
        write_per_genome_tsv(abr / target, picks, ratio, rng)


def main():
    """CLI entrypoint: parse args, discover files, and build communities."""
    ap = argparse.ArgumentParser(description="Create training communities")

    # Reproducibility and community counts
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--short", type=int, default=0, help="Number of short communities")
    ap.add_argument("--long", type=int, default=0, help="Number of long communities")

    # I/O roots and discovery patterns
    ap.add_argument("--outroot", type=Path, required=True)
    ap.add_argument("--roots-prok", nargs="+", required=True)
    ap.add_argument("--roots-euk", nargs="+", required=True)
    ap.add_argument("--roots-vir", nargs="+", required=True)
    ap.add_argument("--roots-pla", nargs="+", required=True)
    ap.add_argument(
        "--patterns",
        nargs="+",
        default=["*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fa.gz", "*.fasta.gz"],
    )

    # Sampling ranges and ratio presets
    ap.add_argument("--ranges", type=str, default="", help="JSON for per-class ranges")
    ap.add_argument(
        "--ratios-generic",
        type=str,
        required=True,
        help='JSON like {"prokaryote":0.56,"microeukaryote":0.24,"virus":0.10,"plasmid":0.10}',
    )
    ap.add_argument("--ratios-filtered", type=str, required=True)

    args = ap.parse_args()

    # RNG setup and ensure output root exists
    rng = random.Random(args.seed)
    outroot = args.outroot
    outroot.mkdir(parents=True, exist_ok=True)

    # Default per-class [lo, hi] ranges (inclusive) if none provided
    ranges = json.loads(args.ranges) if args.ranges else {
        "prokaryote": [100, 200],
        "microeukaryote": [50, 100],
        "virus": [120, 200],
        "plasmid": [120, 200],
    }
    # Normalize to tuple[int, int] for convenience/safety
    ranges = {k: (int(v[0]), int(v[1])) for k, v in ranges.items()}

    # Map classes to their respective search roots
    roots = {
        "prokaryote": args.roots_prok,
        "microeukaryote": args.roots_euk,
        "virus": args.roots_vir,
        "plasmid": args.roots_pla,
    }

    patterns = args.patterns

    # Discover candidate genome files per class
    files_by_cls = {cls: list_fastas(roots[cls], patterns) for cls in CLASSES}

    # Early validation: ensure we found at least one file per class
    for cls, lst in files_by_cls.items():
        if not lst:
            sys.exit(f"[ERR] No FASTA files for class '{cls}' under roots: {roots[cls]}")

    # Parse ratio presets (must include all classes as keys)
    generic = json.loads(args.ratios_generic)
    filtered = json.loads(args.ratios_filtered)

    # Build requested communities
    for i in range(1, args.short + 1):
        make_one(f"TrainC_short_{i:02d}", outroot, files_by_cls, ranges, generic, filtered, rng)
    for i in range(1, args.long + 1):
        make_one(f"TrainC_long_{i:02d}", outroot, files_by_cls, ranges, generic, filtered, rng)

    print(f"[OK] Training communities created under {outroot}/")


if __name__ == "__main__":
    main()
