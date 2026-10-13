#!/usr/bin/env python3
"""
pipeline_adapter.py  —  the seam between the Atlas and the BioLab pipeline.

Two tiers, both INTEGRITY-bound (never fabricate; fail loud when data is absent):

  TIER A  (runs here, on the genome set the Atlas already uses):
    rbh_ortholog(gene)        reciprocal-best-hit homology across all genomes, with an
                              identity threshold and an honest confidence label. This is the
                              real fix for "archaea show 0/18 by symbol" — it finds radB/radA-
                              class homologs symbol-matching can't. It is a PROXY; it labels
                              itself as such and defers to Tier B orthogroups as definitive.
    differential_presence(..) aggregate/sift: which anchor genes are confidently present in
                              one set of genomes and absent in another -> candidate genes.

  TIER B  (the heavy BioLab pipeline — orthogroups, multi-lineage DE, convergence scoring):
    orthogroups / expression_de / convergence_candidates read the pipeline's REAL output
    files. If those outputs are not present, they raise PipelineUnavailable with the path to
    wire — they do NOT invent results. Claude Code connects these to orchestrate.py.

Nothing here returns a number it did not compute or read from real data.
"""
import os, sys, functools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "build"))   # reuse the validated build helpers
import build_atlas as B
import atlas_meta as M

# Where the heavy pipeline writes its outputs. Override with env BIOLAB_OUT.
PIPELINE_OUT = os.environ.get("BIOLAB_OUT", os.path.join(HERE, "..", "..", "pipeline", "output"))

# Honest homology thresholds (amino-acid % identity over the aligned region)
ID_HIGH = 30.0   # >= this AND reciprocal  -> confident homolog
ID_WEAK = 22.0   # 22-30% or non-reciprocal -> weak/uncertain
# below ID_WEAK -> reported as "no confident homolog"


class PipelineUnavailable(RuntimeError):
    """Raised when a Tier-B result is requested but the pipeline output isn't present.
    We raise instead of fabricating."""


# ---------------- Tier A: real homology over the Atlas genome set ----------------
def _kmers(s, k=4):
    return frozenset(s[i:i+k] for i in range(len(s) - k + 1))

@functools.lru_cache(maxsize=1)
def _proteomes():
    """Per-genome list of (locus, gene, protein, kmerset). kmers precomputed for fast prefilter."""
    P = {}
    for s in M.SPECIES:
        recs = []
        for h, seq in B.parse_fasta(B.fetch_cds(M.SPECIES[s]["acc"])):
            p = B.protein(seq)
            if len(p) >= 30:
                recs.append((B.tagf(h, "locus_tag"), B.tagf(h, "gene"), p, _kmers(p)))
        P[s] = recs
    return P

def _best_hit(query, prot_list, topN=40):
    qk = _kmers(query)
    nq = len(qk) or 1
    # cosine-like normalized k-mer overlap removes the long-protein bias of raw counts
    cand = sorted(((len(qk & km) / ((nq * len(km)) ** 0.5 or 1), i)
                   for i, (_, _, _, km) in enumerate(prot_list) if km),
                  reverse=True)[:topN]
    best = None
    for _, i in cand:
        loc, gene, p, _ = prot_list[i]
        sc = B._aligner.align(query, p)[0].score
        if best is None or sc > best[0]:
            best = (sc, B.identity(query, p), loc, gene, i)
    return best  # (score, identity, locus, gene, index) or None

def _dr1_protein(gene):
    idx = B.index_by_symbol(B.fetch_cds(M.SPECIES["DR1"]["acc"]))
    syms = M.SYN.get(gene, [gene])
    for s in syms:
        if s in idx:
            return B.protein(idx[s]["seq"])
    return None

@functools.lru_cache(maxsize=128)
def rbh_ortholog(gene):
    """Reciprocal-best-hit homology for a DR1 anchor gene across every genome.
    Returns a list of per-organism dicts with an honest confidence label."""
    q = _dr1_protein(gene)
    if not q:
        raise ValueError(f"{gene} is not an annotated DR1 anchor protein")
    P = _proteomes()
    dr1 = P["DR1"]
    out = []
    for s in M.SPECIES:
        fwd = _best_hit(q, P[s])
        if not fwd:
            out.append({"sp": s, "confidence": "none", "hit": None, "identity": None,
                        "rbh": False}); continue
        # reciprocity: does the forward hit map back to our query protein in DR1?
        back = _best_hit(P[s][fwd[4]][2], dr1)
        rbh = bool(back and back[1] is not None and
                   (P["DR1"][back[4]][2] == q))
        idt = fwd[1]
        if rbh and idt >= ID_HIGH:
            conf = "high"
        elif (rbh and idt >= ID_WEAK) or (idt >= ID_HIGH):
            conf = "weak"
        else:
            conf = "none"
        out.append({"sp": s, "hit": fwd[3] or "(unnamed)", "locus": fwd[2],
                    "identity": round(idt, 1) if idt is not None else None,
                    "rbh": rbh, "confidence": conf})
    return {"gene": gene, "method": "reciprocal-best-hit (BLOSUM62), proxy",
            "note": "Confidence is a quick-homology proxy. Definitive orthology = Tier-B "
                    "pipeline orthogroups. 'none' means no reliable homolog by this method "
                    "— not proof of true absence.",
            "results": out}

@functools.lru_cache(maxsize=128)
def _homolog_presence(gene):
    """Forward-best-hit presence per genome (lighter than full RBH; used for aggregate sifting)."""
    q = _dr1_protein(gene)
    if not q:
        raise ValueError(f"{gene} is not an annotated DR1 anchor protein")
    P = _proteomes()
    out = {}
    for s in M.SPECIES:
        bh = _best_hit(q, P[s])
        idt = bh[1] if bh else None
        out[s] = {"identity": idt, "hit": bh[3] if bh else None,
                  "present": bool(idt is not None and idt >= ID_HIGH)}
    return out

def differential_presence(present_in, absent_in, genes=None):
    """Sift: anchor genes confidently present across `present_in` and absent across `absent_in`."""
    genes = genes or M.ANCHORS
    hits = []
    for g in genes:
        try:
            pr = _homolog_presence(g)
        except ValueError:
            continue
        if (all(pr[s]["present"] for s in present_in if s in pr) and
                all(not pr[s]["present"] for s in absent_in if s in pr)):
            hits.append({"gene": g, "present": {s: pr[s]["identity"] for s in present_in if s in pr},
                         "absent": absent_in})
    return {"present_in": present_in, "absent_in": absent_in,
            "method": "forward-best-hit homology over the anchor set (proxy, id>=30%)",
            "candidates": hits}

@functools.lru_cache(maxsize=1)
def convergence_preview():
    """Lineage-aware sift over the anchor set: for each gene, how many INDEPENDENT radioresistant
    lineages carry a confident homolog, and whether the sensitive controls do too. Labels each gene
    universal / lineage-specific / convergent-candidate. A PROXY; genome-wide convergence = Tier B."""
    resist = [s for s in M.SPECIES if M.SPECIES[s]["kind"] != "sensitive control"]
    sens = [s for s in M.SPECIES if M.SPECIES[s]["kind"] == "sensitive control"]
    lin = {s: M.SPECIES[s].get("lin", "") for s in M.SPECIES}
    out = []
    for g in M.ANCHORS:
        pr = _homolog_presence(g)
        present = {s for s in M.SPECIES if pr[s]["present"]}
        rlins = sorted({lin[s] for s in present if s in resist})
        in_sens = [s for s in sens if s in present]
        n = len(rlins)
        if all(s in present for s in sens) and n >= 2:
            label = "universal"
        elif n >= 2 and not in_sens:
            label = "convergent candidate"
        elif n == 1 and not in_sens:
            label = "lineage-specific"
        elif n == 0:
            label = "sparse"
        else:
            label = "mixed"
        out.append({"gene": g, "lineages": rlins, "n_lineages": n,
                    "in_sensitive": in_sens, "label": label})
    rank = {"convergent candidate": 0, "lineage-specific": 1, "mixed": 2, "universal": 3, "sparse": 4}
    out.sort(key=lambda x: (rank[x["label"]], -x["n_lineages"]))
    return {"method": "lineage-aware forward-best-hit homology over the anchor set (proxy, id>=30%)",
            "note": "Genome-wide convergence over ALL orthogroups is a Tier-B pipeline job; this "
                    "previews the signal on the 18-gene anchor set. universal = core repair shared "
                    "by everything incl. sensitive; lineage-specific = unique to one lineage (e.g. the "
                    "Deinococcus ddr genes); convergent candidate = present across multiple independent "
                    "radioresistant lineages but absent in the sensitive controls. Convergent "
                    "candidates are LEADS to verify with Tier-B orthogroups (a >=30% best-hit can be "
                    "spurious) — not conclusions.",
            "radioresistant_lineages": sorted(set(lin[s] for s in resist)),
            "results": out}


# ---------------- Tier B: the heavy pipeline (fail loud, never fabricate) ----------------
def _require(path, what):
    if not os.path.exists(path):
        raise PipelineUnavailable(
            f"{what} not found at {path}. Wire this to the BioLab pipeline output "
            f"(set BIOLAB_OUT or generate it), or run the pipeline. Not fabricating.")
    return path

def orthogroups(genes=None):
    """Read OrthoFinder-style Orthogroups.tsv from the pipeline. Raises if absent."""
    p = _require(os.path.join(PIPELINE_OUT, "orthogroups", "Orthogroups.tsv"),
                 "Pipeline orthogroups")
    rows = [l.rstrip("\n").split("\t") for l in open(p)]
    header, data = rows[0], rows[1:]
    ogs = [dict(zip(header, r)) for r in data]
    if genes:
        gl = set(genes)
        ogs = [o for o in ogs if any(g in (o.get(c, "")) for g in gl for c in header[1:])]
    return {"source": p, "n_orthogroups": len(ogs), "orthogroups": ogs[:200]}

def expression_de(sp, geo):
    """Read a real DESeq2 result table the pipeline produced for (sp, geo). Raises if absent."""
    p = _require(os.path.join(PIPELINE_OUT, "de", f"{sp}_{geo}_deseq2.csv"),
                 f"DESeq2 result for {sp}/{geo}")
    import csv
    with open(p) as f:
        rows = list(csv.DictReader(f))
    return {"source": p, "sp": sp, "geo": geo, "n_genes": len(rows), "results": rows[:500]}

def convergence_candidates():
    """Read the pipeline's cross-lineage convergence candidate scores. Raises if absent."""
    p = _require(os.path.join(PIPELINE_OUT, "convergence", "candidates.json"),
                 "Convergence candidates")
    import json
    return json.load(open(p))


def capabilities():
    """Report what is live right now — without fabricating availability."""
    def ok(path):
        return os.path.exists(path)
    return {
        "tier_a_homology": True,  # always available: computes on the genome set
        "tier_a_convergence": True,
        "genomes": list(M.SPECIES.keys()),
        "n_lineages": len(set(M.SPECIES[s].get("lin", "") for s in M.SPECIES)),
        "tier_b_orthogroups": ok(os.path.join(PIPELINE_OUT, "orthogroups", "Orthogroups.tsv")),
        "tier_b_de_dir": ok(os.path.join(PIPELINE_OUT, "de")),
        "tier_b_convergence": ok(os.path.join(PIPELINE_OUT, "convergence", "candidates.json")),
        "pipeline_out": os.path.abspath(PIPELINE_OUT),
    }
