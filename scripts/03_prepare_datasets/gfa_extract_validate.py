#!/usr/bin/env python3
import sys, os, re, gzip, json, argparse
from collections import defaultdict, Counter
from pathlib import Path

def open_any(p):
    p = str(p)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")

COV_RE = re.compile(r"(?:^|_)cov_([0-9]*\.?[0-9]+)")

def gc_of(seq):
    if not seq: return 0.0
    atgc = sum(seq.count(x) for x in "ACGTacgt")
    if atgc == 0: return 0.0
    g = seq.count("G")+seq.count("g")
    c = seq.count("C")+seq.count("c")
    return (g+c)/atgc

def choose_gfa(asm_dir: str):
    asm_dir = Path(asm_dir)
    cands = [
        asm_dir / "assembly_graph_with_scaffolds.gfa",
        asm_dir / "assembly_graph.gfa",
        asm_dir / "assembly_graph_with_scaffolds.gfa.gz",
        asm_dir / "assembly_graph.gfa.gz",
    ]
    for c in cands:
        if c.is_file() and c.stat().st_size > 0:
            return str(c)
    # fallback: pick the biggest *.gfa*
    picks = sorted(asm_dir.glob("*.gfa*"), key=lambda p: p.stat().st_size if p.exists() else -1, reverse=True)
    return str(picks[0]) if picks else None

def parse_gfa(gfa_path, keep_seq=False):
    segs = {}            # id -> dict
    edges = []           # dicts: src,dst,src_orient,dst_orient,overlap
    contig_paths = {}    # contig -> [(unitig, +/-), ...]

    def cov_from_id(sid):
        m = COV_RE.search(sid)
        return float(m.group(1)) if m else None

    def parse_cov(tags):
        cov = None
        for t in tags:
            # DP:f:xx, dp:f:xx, cv:f:xx, cv:i:xx (seen in some graphs)
            if t.startswith(("DP:f:","dp:f:","cv:f:","CV:f:")):
                try: cov = float(t.split(":")[-1]); break
                except: pass
            if t.startswith(("cv:i:","CV:i:")):
                try: cov = float(t.split(":")[-1]); break
                except: pass
        return cov

    with open_any(gfa_path) as f:
        for ln,line in enumerate(f,1):
            if not line or line[0] == '#': continue
            parts = line.rstrip("\n").split("\t")
            tag = parts[0]
            if tag == "S":
                sid = parts[1]
                seq = parts[2]
                cov = parse_cov(parts[3:]) or cov_from_id(sid)
                L = len(seq) if seq != "*" else None
                GC = gc_of(seq) if seq != "*" else None
                segs[sid] = {
                    "node_id": sid,
                    "length": L,
                    "gc": GC,
                    "cov": cov,
                    "has_seq": int(seq != "*"),
                    "seq": (seq if keep_seq and seq != "*" else "")
                }
            elif tag == "L":
                if len(parts) < 6:
                    continue
                u, ou, v, ov, ovl = parts[1:6]
                if ou not in ("+","-") or ov not in ("+","-"):
                    continue
                m = re.match(r"(\d+)", ovl)  # overlap like 10M
                ovl_bp = int(m.group(1)) if m else None
                edges.append({
                    "src": u, "dst": v,
                    "src_orient": ou, "dst_orient": ov,
                    "overlap_raw": ovl, "overlap": ovl_bp
                })
            elif tag == "P":
                contig = parts[1]
                walk = parts[2].split(",")
                parsed = []
                for w in walk:
                    m = re.match(r"(.+)([+-])$", w)
                    if m: parsed.append((m.group(1), m.group(2)))
                contig_paths[contig] = parsed
            # ignore others
    return segs, edges, contig_paths

def contig_graph_from_paths(contig_paths, kmin_shared=1):
    # connect contigs sharing >= kmin_shared unitigs
    u2c = defaultdict(set)
    for c, walk in contig_paths.items():
        for u,_ in walk:
            u2c[u].add(c)
    es = set()
    if kmin_shared <= 1:
        for u, CC in u2c.items():
            CC = sorted(CC)
            for i in range(len(CC)):
                for j in range(i+1, len(CC)):
                    es.add((CC[i], CC[j]))
    else:
        # count shared unitigs
        pairs = Counter()
        for CC in u2c.values():
            CC = sorted(CC)
            for i in range(len(CC)):
                for j in range(i+1, len(CC)):
                    pairs[(CC[i], CC[j])] += 1
        es = {ab for ab, cnt in pairs.items() if cnt >= kmin_shared}
    return sorted(es)

def write_tsv_nodes(path, segs):
    with open(path,"w") as out:
        out.write("node_id\tlength\tgc\tcov\n")
        for sid in sorted(segs.keys()):
            d = segs[sid]
            L = "" if d["length"] is None else d["length"]
            GC = "" if d["gc"] is None else f"{d['gc']:.6f}"
            COV= "" if d["cov"] is None else f"{d['cov']:.6f}"
            out.write(f"{sid}\t{L}\t{GC}\t{COV}\n")

def write_tsv_edges(path, edges):
    with open(path,"w") as out:
        out.write("src\tdst\tsrc_orient\tdst_orient\toverlap\n")
        for e in sorted(edges, key=lambda x: (x["src"], x["dst"])):
            ol = "" if e["overlap"] is None else e["overlap"]
            out.write(f"{e['src']}\t{e['dst']}\t{e['src_orient']}\t{e['dst_orient']}\t{ol}\n")

def write_tsv_contigs_nodes(path, contig_paths):
    with open(path,"w") as out:
        out.write("contig_id\n")
        for c in sorted(contig_paths.keys()):
            out.write(f"{c}\n")

def write_tsv_contigs_edges(path, contig_edges):
    with open(path,"w") as out:
        out.write("src\tdst\n")
        for a,b in sorted(contig_edges):
            out.write(f"{a}\t{b}\n")

def sanity_checks(segs, edges, contig_paths, contigs_ge500_fa=None):
    rep = {}
    rep["n_nodes"] = len(segs)
    rep["n_edges"] = len(edges)
    rep["n_paths"] = len(contig_paths)
    rep["n_nodes_no_seq_field"] = sum(1 for d in segs.values() if d["has_seq"] == 0)
    rep["n_nodes_zero_len"]     = sum(1 for d in segs.values() if (d["length"] or 0) == 0)
    bad_ori = sum(1 for e in edges if e["src_orient"] not in "+-" or e["dst_orient"] not in "+-")
    rep["n_edges_bad_orient"] = bad_ori

    missing_src = [e for e in edges if e["src"] not in segs]
    missing_dst = [e for e in edges if e["dst"] not in segs]
    rep["n_edges_missing_src"] = len(missing_src)
    rep["n_edges_missing_dst"] = len(missing_dst)

    deg = Counter()
    for e in edges:
        if e["src"] in segs: deg[e["src"]] += 1
        if e["dst"] in segs: deg[e["dst"]] += 1
    if deg:
        vals = sorted(deg.values())
        rep["degree_min"] = int(vals[0])
        rep["degree_med"] = int(vals[len(vals)//2])
        rep["degree_max"] = int(vals[-1])
    else:
        rep["degree_min"] = rep["degree_med"] = rep["degree_max"] = 0

    rep["contigs_ge500_in_paths"] = None
    if contigs_ge500_fa and os.path.isfile(contigs_ge500_fa):
        total = 0
        names = set()
        with open_any(contigs_ge500_fa) as f:
            for line in f:
                if line.startswith(">"):
                    total += 1
                    nm = line[1:].strip().split()[0]
                    names.add(nm)
        kept = sum(1 for nm in names if nm in contig_paths)
        rep["contigs_ge500_in_paths"] = {"total_headers": total, "present_in_paths": kept}

    connected = set()
    for e in edges:
        connected.add(e["src"]); connected.add(e["dst"])
    rep["n_isolated_nodes"] = sum(1 for nid in segs if nid not in connected)
    return rep

def autodiscover_contigs(asm_dir: Path) -> str | None:
    # Try to find eval*/contigs*_ge*.fa* around asm_dir
    base = asm_dir.parent
    cands = []
    for d in base.glob("eval*"):
        cands += list(d.glob("contigs*_ge*.fa*"))
    if not cands:
        # Also check asm_dir itself
        cands = list(asm_dir.glob("contigs*_ge*.fa*"))
    if not cands:
        return None
    cands.sort(key=lambda p: p.stat().st_size if p.exists() else -1, reverse=True)
    return str(cands[0])

def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--asm-dir", help="Assembly directory containing GFA(s)")
    p.add_argument("--out-dir", required=True, help="Output directory for TSV/JSON")
    p.add_argument("--gfa", default=None, help="Explicit GFA path (overrides autodetect)")
    p.add_argument("--keep-seq", action="store_true", help="Keep segment sequences in nodes.json (seq field)")
    p.add_argument("--contigs-fa", default=None, help="Path to ≥500bp contigs FASTA (for sanity check); autodetect if not provided")
    p.add_argument("--min-shared-unitigs", type=int, default=1, help="Edges between contigs require at least this many shared unitigs")
    p.add_argument("--no-contig-graph", action="store_true", help="Skip writing contig_nodes/contig_edges TSVs")
    args = p.parse_args()

    if not args.gfa and not args.asm_dir:
        print("ERROR: Provide --gfa or --asm-dir", file=sys.stderr)
        sys.exit(2)

    asm_dir = Path(args.asm_dir) if args.asm_dir else Path(Path(args.gfa).parent)
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    gfa = args.gfa or choose_gfa(str(asm_dir))
    if not gfa:
        print(f"[ERR] No GFA found under {asm_dir}", file=sys.stderr)
        sys.exit(3)

    segs, edges, contig_paths = parse_gfa(gfa, keep_seq=args.keep_seq)

    nodes_tsv = out_dir / "nodes.tsv"
    edges_tsv = out_dir / "edges.tsv"
    write_tsv_nodes(str(nodes_tsv), segs)
    write_tsv_edges(str(edges_tsv), edges)

    if not args.no_contig_graph and contig_paths:
        contig_nodes_tsv = out_dir / "contig_nodes.tsv"
        contig_edges_tsv = out_dir / "contig_edges.tsv"
        write_tsv_contigs_nodes(str(contig_nodes_tsv), contig_paths)
        c_edges = contig_graph_from_paths(contig_paths, kmin_shared=args.min_shared_unitigs)
        write_tsv_contigs_edges(str(contig_edges_tsv), c_edges)

    ge500 = args.contigs_fa or autodiscover_contigs(asm_dir)
    report = sanity_checks(segs, edges, contig_paths, ge500)
    report["asm_dir"] = str(asm_dir)
    report["gfa"] = str(gfa)
    report["contigs_ge500"] = ge500

    with open(out_dir / "graph_report.json", "w") as w:
        json.dump(report, w, indent=2)

    print("[GRAPH] Summary")
    for k,v in report.items():
        print(f"  {k}: {v}")
    print(f"[OK] nodes.tsv: {nodes_tsv}")
    print(f"[OK] edges.tsv: {edges_tsv}")
    if contig_paths and not args.no_contig_graph:
        print(f"[OK] contig_nodes.tsv: {out_dir/'contig_nodes.tsv'}")
        print(f"[OK] contig_edges.tsv: {out_dir/'contig_edges.tsv'}")

if __name__ == "__main__":
    main()
