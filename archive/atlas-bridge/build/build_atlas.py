#!/usr/bin/env python3
"""
build_atlas.py  —  Reproducible build for BioLab Atlas.

What it does (all from REAL public data — nothing invented):
  1. Downloads CDS FASTA for each genome in atlas_meta.SPECIES from NCBI RefSeq (cached).
  2. Extracts each anchor gene by annotated symbol  -> real coding sequence + locus + length.
  3. Computes GC%, protein length, and pairwise amino-acid identity (BLOSUM62 global align).
  4. Merges with curated META/expr  -> writes genes_data.js
  5. Injects genes_data.js into atlas_template.html  -> writes ../index.html

Run:   python build_atlas.py
Deps:  biopython  (pip install biopython)

INTEGRITY: this script only reports what it can derive from RefSeq + the curated META.
Presence is by annotated SYMBOL (a proxy; the pipeline's orthogroups are definitive).
Identity is a global amino-acid alignment. None of these are hand-edited.
"""
import os, re, gzip, json, sys, urllib.request
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices
import atlas_meta as M

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache")
TEMPLATE = os.path.join(HERE, "atlas_template.html")
OUT_DATA = os.path.join(HERE, "genes_data.js")
OUT_HTML = os.path.join(os.path.dirname(HERE), "index.html")
WIN = 120          # bp of the 5' coding window shown on the helix
os.makedirs(CACHE, exist_ok=True)

# ---------- fetch ----------
def ftp_cds_url(acc):
    p = acc.split("_")[1].split(".")[0]
    base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{p[0:3]}/{p[3:6]}/{p[6:9]}/"
    html = urllib.request.urlopen(base, timeout=60).read().decode()
    m = re.search(re.escape(acc) + r"_[A-Za-z0-9._]+/", html)
    if not m:
        raise RuntimeError(f"Could not resolve assembly dir for {acc} under {base}")
    d = m.group(0).strip("/")
    return base + d + "/" + d + "_cds_from_genomic.fna.gz"

def fetch_cds(acc):
    out = os.path.join(CACHE, acc + ".cds.fna")
    if os.path.exists(out):
        return out
    url = ftp_cds_url(acc)
    print(f"  downloading {acc}  ({url.split('/')[-1]})")
    raw = urllib.request.urlopen(url, timeout=180).read()
    open(out, "w").write(gzip.decompress(raw).decode())
    return out

# ---------- parse / index ----------
def parse_fasta(path):
    recs, hdr, seq = [], None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if hdr is not None: recs.append((hdr, "".join(seq)))
            hdr, seq = line, []
        else:
            seq.append(line)
    if hdr is not None: recs.append((hdr, "".join(seq)))
    return recs

def tagf(h, k):
    m = re.search(r"\[" + k + r"=([^\]]+)\]", h)
    return m.group(1) if m else None

def index_by_symbol(path):
    out = {}
    for h, s in parse_fasta(path):
        g, lt = tagf(h, "gene"), tagf(h, "locus_tag")
        if g and g not in out:
            out[g] = {"locus": lt, "seq": s}
    return out

def lookup(idx, gene):
    for sym in M.SYN.get(gene, [gene]):
        if sym in idx:
            return idx[sym]
    return None

# ---------- metrics ----------
def gc(s):
    s = s.upper(); n = len(s)
    return round(100 * (s.count("G") + s.count("C")) / n, 1) if n else None

def protein(s):
    try:
        return str(Seq(s).translate(table=11, to_stop=True))
    except Exception:
        return ""

_aligner = PairwiseAligner()
_aligner.mode = "global"
_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
_aligner.open_gap_score = -10
_aligner.extend_gap_score = -0.5

def identity(p1, p2):
    if not p1 or not p2:
        return None
    aln = _aligner.align(p1, p2)[0]
    ident = 0
    for (s1, e1), (s2, e2) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(e1 - s1):
            if p1[s1 + k] == p2[s2 + k]:
                ident += 1
    return round(100 * ident / min(len(p1), len(p2)), 1)

# ---------- build ----------
def main():
    sp_order = list(M.SPECIES.keys())
    spi = {s: i for i, s in enumerate(sp_order)}
    print("Fetching genomes…")
    idx = {s: index_by_symbol(fetch_cds(M.SPECIES[s]["acc"])) for s in sp_order}

    print("Extracting genes + computing metrics…")
    genes = []
    for g in M.ANCHORS:
        meta = M.META[g]
        sp, prots = {}, {}
        for s in sp_order:
            r = lookup(idx[s], g)
            if r:
                sp[s] = {"locus": r["locus"], "len": len(r["seq"]),
                         "seq": r["seq"][:WIN], "gc": gc(r["seq"]), "aa": len(protein(r["seq"]))}
                prots[s] = protein(r["seq"])
            else:
                sp[s] = None
        idt = {}
        for i in range(len(sp_order)):
            for j in range(i + 1, len(sp_order)):
                a, b = sp_order[i], sp_order[j]
                if a in prots and b in prots:
                    idt[a + "|" + b] = identity(prots[a], prots[b])
        rec = {"g": g, **meta, "sp": sp, "idt": idt}
        genes.append(rec)
        present = [s for s in sp_order if sp[s]]
        print(f"  {g:6s} present in {len(present)}/{len(sp_order)}: {','.join(present)}")

    open(OUT_DATA, "w").write("const GENES=" + json.dumps(genes, separators=(",", ":")) + ";")
    print(f"Wrote {OUT_DATA}")

    # Single source of truth: project SPECIES/PHENO from atlas_meta into the page
    sp_js = {s: {"n": M.SPECIES[s]["name"], "s": M.SPECIES[s]["strain"],
                 "kind": M.SPECIES[s]["kind"], "lin": M.SPECIES[s].get("lin", ""),
                 "dom": M.SPECIES[s].get("dom", ""),
                 "distant": bool(M.SPECIES[s].get("distant", False))} for s in sp_order}
    asm_js = {s: M.SPECIES[s]["acc"] for s in sp_order}
    repl = {
        "/*__GENES__*/":   open(OUT_DATA).read(),
        "/*__SPECIES__*/": json.dumps(sp_js, ensure_ascii=False),
        "/*__SP_ORDER__*/": json.dumps(sp_order),
        "/*__ASM__*/":     json.dumps(asm_js),
        "/*__PHENO__*/":   json.dumps(M.PHENO, ensure_ascii=False),
    }
    html = open(TEMPLATE).read()
    for k, v in repl.items():
        if k not in html:
            sys.exit(f"ERROR: template is missing the {k} placeholder.")
        html = html.replace(k, v)
    open(OUT_HTML, "w").write(html)
    print(f"Wrote {OUT_HTML}  ({len(html)} bytes)")
    print("Done. Deploy index.html, or `vercel --prod` from the atlas/ folder.")

if __name__ == "__main__":
    main()
