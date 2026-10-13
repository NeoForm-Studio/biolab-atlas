"""Measure a proteome: amino-acid composition and trait marker proteins (via local NCBI BLAST+). No network."""
import shutil
import subprocess
import tempfile
from pathlib import Path

IVYWREL = set("IVYWREL")
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def parse_fasta(text):
    """Return a list of (header, sequence) pairs."""
    records, header, chunks = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:].strip(), []
        elif header is not None:
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def composition(records):
    """Proteome-wide amino-acid signals.

    ivywrel: fraction of I, V, Y, W, R, E, L residues. Zeldovich et al. 2007 (PLoS Comput Biol 3:e5)
    found it tracks optimal growth temperature as T = 937 * F - 335 (r = 0.93 across prokaryotes).
    acidic_excess: (D + E - K - R) / residues. Salt-in halophiles have strongly acidic proteomes.
    """
    counts, total = {}, 0
    for _, seq in records:
        for aa in seq.upper():
            if aa in AMINO_ACIDS:  # skip stops, X and other ambiguity codes
                counts[aa] = counts.get(aa, 0) + 1
                total += 1
    if not total:
        return None
    f = lambda letters: sum(counts.get(a, 0) for a in letters) / total
    ivywrel = f(IVYWREL)
    return {
        "proteins": len(records),
        "ivywrel": round(ivywrel, 4),
        "predicted_ogt_c": round(937 * ivywrel - 335),
        "acidic_excess": round(f("DE") - f("KR"), 4),
    }


MIN_COVERAGE, MIN_IDENTITY, MAX_EVALUE = 70, 30, 1e-10


def best_hits(tsv):
    """blastp -outfmt '6 qseqid sseqid pident qcovs evalue bitscore' -> {query: best hit by bitscore}."""
    best = {}
    for line in tsv.splitlines():
        q, s, pident, qcov, evalue, bits = line.split("\t")
        hit = {"accession": s, "identity": float(pident), "coverage": int(qcov), "evalue": float(evalue),
               "bits": float(bits)}
        if q not in best or hit["bits"] > best[q]["bits"]:
            best[q] = hit
    return best


def passes(hit, gene):
    """A BLAST hit counts as the marker when it is significant, long enough and similar enough."""
    return bool(hit) and (hit["evalue"] <= MAX_EVALUE and hit["coverage"] >= MIN_COVERAGE
                          and hit["identity"] >= gene.get("min_identity", MIN_IDENTITY))


def find_markers(records, panels, reference_faa):
    """For each trait panel, BLAST the reference marker proteins against this proteome."""
    if not shutil.which("blastp"):
        raise RuntimeError("NCBI BLAST+ (blastp) is required; install it and make sure blastp is on PATH")
    with tempfile.TemporaryDirectory() as tmp:
        subject = Path(tmp) / "proteome.faa"
        subject.write_text("".join(f">{h.split()[0]}\n{seq}\n" for h, seq in records))
        run = subprocess.run(["blastp", "-query", str(reference_faa), "-subject", str(subject), "-evalue", "1e-5",
                              "-max_hsps", "1", "-outfmt", "6 qseqid sseqid pident qcovs evalue bitscore"],
                             capture_output=True, text=True, check=True)
    best = best_hits(run.stdout)
    out = {}
    for trait, panel in panels.items():
        genes = []
        for gene in panel["genes"]:
            hit = best.get(gene["ref"])
            found = passes(hit, gene)
            genes.append({"name": gene["name"], "found": found,
                          **({k: hit[k] for k in ("accession", "identity", "coverage")} if hit else {})})
        out[trait] = {"found": sum(g["found"] for g in genes), "of": len(genes), "genes": genes}
    return out
