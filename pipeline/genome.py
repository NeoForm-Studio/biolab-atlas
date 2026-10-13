"""Pure functions: parse a proteome and measure it. No network access here, so it is easy to test."""
import re

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


def product(header):
    """'WP_1.1 MULTISPECIES: recombinase RecA [Deinococcus]' -> 'recombinase RecA'."""
    name = re.sub(r"^\S+\s+(MULTISPECIES:\s+)?", "", header)
    return name.rsplit(" [", 1)[0].strip()


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


def find_markers(records, panels):
    """For each trait panel, list which marker proteins are annotated in this proteome."""
    products = [(h.split()[0], product(h)) for h, _ in records]
    out = {}
    for trait, panel in panels.items():
        genes = []
        for gene in panel["genes"]:
            pattern = re.compile(gene["match"], re.IGNORECASE)
            hits = [acc for acc, name in products if pattern.fullmatch(name)]
            genes.append({"name": gene["name"], "found": bool(hits), "accessions": hits[:5]})
        out[trait] = {  # labels and notes live in markers.json
            "found": sum(g["found"] for g in genes),
            "of": len(genes),
            "genes": genes,
        }
    return out
