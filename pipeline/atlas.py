"""BioLab Atlas pipeline: look up any bacterium or archaeon in public databases and profile its genome.

    python pipeline/atlas.py "Deinococcus radiodurans" "Thermus thermophilus" ...

For each name it writes library/<slug>.json, then rebuilds library/index.json.
Sources: NCBI Taxonomy + Datasets, NCBI nucleotide (16S records), UniProt, Europe PMC, BacDive.
Needs NCBI BLAST+ (blastp) on PATH; otherwise standard library only.
"""
import http.client
import io
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from datetime import date
from pathlib import Path

from genome import composition, find_markers, parse_fasta
from score import score_library

VERSION = "0.1"
HERE = Path(__file__).resolve().parent
LIBRARY = HERE.parent / "library"
PANELS = json.loads((HERE / "markers.json").read_text(encoding="utf-8"))
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BACDIVE = "https://api.bacdive.dsmz.de/v2"
DOMAINS = {2: "Bacteria", 2157: "Archaea"}
SPECIES_OR_BELOW = {"SPECIES", "SUBSPECIES", "STRAIN", "SEROTYPE", "SEROGROUP", "BIOTYPE", "BIOVAR", "ISOLATE",
                    "PATHOGROUP", "FORMA_SPECIALIS"}
LEVELS = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}


# ---------------------------------------------------------------- HTTP
def fetch(url, want="json"):
    """GET with retries. want = 'json' | 'bytes' | 'response' (body + headers)."""
    req = urllib.request.Request(url, headers={"User-Agent": f"biolab-atlas/{VERSION}"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
                time.sleep(0.35)  # stay under NCBI's 3 requests/second without an API key
                if want == "json":
                    return json.loads(body)
                return body if want == "bytes" else (body, r.headers)
        except urllib.error.HTTPError as e:
            if e.code not in (429, 500, 502, 503, 504) or attempt == 3:
                raise
        except (urllib.error.URLError, TimeoutError, ConnectionError, http.client.HTTPException):
            if attempt == 3:
                raise
        time.sleep(2 ** attempt)


def q(text):
    return urllib.parse.quote(text, safe="")


# ---------------------------------------------------------------- sources
def resolve(name):
    """NCBI Taxonomy node for a name or tax id, or None if NCBI does not recognise it."""
    nodes = fetch(f"{DATASETS}/taxonomy/taxon/{q(name)}").get("taxonomy_nodes", [])
    taxa = [n["taxonomy"] for n in nodes if "taxonomy" in n]
    # a name can exist in two kingdoms (Bacillus is also a stick insect): prefer the prokaryote
    return next((t for t in taxa if set(t.get("lineage", [])) & DOMAINS.keys()), taxa[0] if taxa else None)


def is_species_level(taxon):
    """True for a species or anything below it. Unranked taxa count only if their parent does."""
    rank = taxon.get("rank", "NO_RANK")
    if rank != "NO_RANK":
        return rank in SPECIES_OR_BELOW
    parent = resolve(str(taxon["lineage"][-1])) if taxon.get("lineage") else None
    return bool(parent) and is_species_level(parent)


def inventory(taxon):
    tax_id, name = taxon["tax_id"], taxon["organism_name"]
    counts = {c["type"]: c["count"] for c in taxon.get("counts", [])}
    _, h = fetch(f"https://rest.uniprot.org/uniprotkb/search?query=taxonomy_id:{tax_id}&size=0", "response")
    _, hr = fetch(f"https://rest.uniprot.org/uniprotkb/search?query=taxonomy_id:{tax_id}+AND+reviewed:true&size=0",
                  "response")
    term = q(f'txid{tax_id}[Organism:exp] AND "16S ribosomal RNA"[Title]')
    s16 = fetch(f"{EUTILS}/esearch.fcgi?db=nuccore&rettype=count&retmode=json&term={term}")
    pmc = fetch(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q(chr(34) + name + chr(34))}"
                f"&format=json&pageSize=1")
    return {
        "assemblies": counts.get("COUNT_TYPE_ASSEMBLY", 0),
        "uniprot_proteins": int(h.get("X-Total-Results", 0)),
        "uniprot_reviewed": int(hr.get("X-Total-Results", 0)),
        "sequences_16s": int(s16["esearchresult"]["count"]),
        "papers": pmc.get("hitCount", 0),
    }


def temperature_range(text):
    """'30-37' -> (30.0, 37.0), '26' -> (26.0, 26.0), anything else -> None."""
    m = re.fullmatch(r"\s*(-?\d+(?:\.\d+)?)\s*(?:-\s*(-?\d+(?:\.\d+)?))?\s*", str(text or ""))
    return (float(m[1]), float(m[2] or m[1])) if m else None


def phenotype(name):
    """Measured growth data from BacDive (DSMZ), preferring the type strain. None if BacDive has no record."""
    words = name.split()
    if len(words) < 2:
        return None
    page = fetch(f"{BACDIVE}/taxon/{q(words[0])}/{q(words[1])}")
    strains = []
    for _ in range(3):  # big species have thousands of strains: look at up to 300 for the type strain
        if not page.get("results"):
            break
        batch = list(fetch(f"{BACDIVE}/fetch/{';'.join(map(str, page['results']))}")["results"].values())
        strains += batch
        if any(s["Name and taxonomic classification"].get("type strain") == "yes" for s in batch) or not page.get("next"):
            break
        page = fetch(page["next"])
    if not strains:
        return None
    has_temp = lambda s: bool(s.get("Culture and growth conditions", {}).get("culture temp"))
    strain = (next((s for s in strains if s["Name and taxonomic classification"].get("type strain") == "yes"), None)
              or next((s for s in strains if has_temp(s)), strains[0]))
    listify = lambda x: x if isinstance(x, list) else [x] if x else []
    temps = [t for t in listify(strain.get("Culture and growth conditions", {}).get("culture temp"))
             if t.get("growth") == "positive"]
    span = lambda rows: (lambda r: [min(r), max(r)] if r else None)(
        [v for t in rows for v in (temperature_range(t.get("temperature")) or ())])
    salt_ranges = [temperature_range(h["concentration"].replace("%", ""))
                   for h in listify(strain.get("Physiology and metabolism", {}).get("halophily"))
                   if h.get("growth") == "positive" and h.get("salt") == "NaCl" and "%" in h.get("concentration", "")]
    salt = [v for r in salt_ranges if r and r[1] <= 30 for v in r]  # above 30% is not a NaCl percentage
    bacdive_id = strain["General"]["BacDive-ID"]
    return {
        "bacdive_id": bacdive_id,
        "url": f"https://bacdive.dsmz.de/strain/{bacdive_id}",
        "type_strain": strain["Name and taxonomic classification"].get("type strain") == "yes",
        "optimum_temp_c": span([t for t in temps if t.get("type") == "optimum"]),
        "growth_temp_c": span(temps),
        "max_nacl_percent": max(salt) if salt else None,
    }


def assemblies(tax_id, extra=""):
    url = f"{DATASETS}/genome/taxon/{tax_id}/dataset_report?page_size=100{extra}"
    return fetch(url).get("reports", [])


def pick_assembly(tax_id):
    """NCBI's reference genome if there is one, else the best annotated, most complete assembly."""
    reports = (assemblies(tax_id, "&filters.reference_only=true")
               or assemblies(tax_id, "&filters.has_annotation=true&filters.assembly_level=complete_genome"
                                     "&filters.assembly_level=chromosome")
               or assemblies(tax_id, "&filters.has_annotation=true")
               or assemblies(tax_id))
    if not reports:
        return None
    return min(reports, key=lambda r: (
        "annotation_info" not in r,
        LEVELS.get(r["assembly_info"].get("assembly_level"), 9),
        not r["accession"].startswith("GCF_"),
        -r["assembly_stats"].get("contig_n50", 0),
    ))


def replicons(accession):
    reports = fetch(f"{DATASETS}/genome/accession/{accession}/sequence_reports").get("reports", [])
    return [{
        "name": r.get("chr_name") if r.get("chr_name") not in (None, "na") else r.get("assigned_molecule_location_type"),
        "type": r.get("assigned_molecule_location_type"),
        "accession": r.get("refseq_accession") or r.get("genbank_accession"),
        "length": r.get("length"),
        "gc_percent": r.get("gc_percent"),
    } for r in reports if r.get("role") == "assembled-molecule"]


def proteome(accession):
    data = fetch(f"{DATASETS}/genome/accession/{accession}/download?include_annotation_type=PROT_FASTA", "bytes")
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        faa = [n for n in z.namelist() if n.endswith("protein.faa")]
        return parse_fasta(z.read(faa[0]).decode()) if faa else []


def relatives(tax_id):
    """Species under this taxon that have an NCBI reference genome."""
    out, seen = [], set()
    for r in assemblies(tax_id, "&filters.reference_only=true"):
        org = r["organism"]["organism_name"]
        if org not in seen:
            seen.add(org)
            out.append({"name": org, "tax_id": r["organism"]["tax_id"], "accession": r["accession"]})
        if len(out) == 5:
            break
    return out


# ---------------------------------------------------------------- profile
def profile(query):
    entry = {"query": query, "generated": date.today().isoformat(), "pipeline_version": VERSION}
    taxon = resolve(query)
    if taxon is None:
        return {**entry, "status": "unknown_name",
                "message": f"No record of '{query}' in NCBI Taxonomy. Check the spelling or try a synonym."}

    lineage = taxon.get("lineage", []) + [taxon["tax_id"]]
    domain = next((DOMAINS[t] for t in lineage if t in DOMAINS), None)
    entry["taxon"] = {"tax_id": taxon["tax_id"], "name": taxon["organism_name"],
                      "rank": taxon.get("rank", "NO_RANK"), "domain": domain}
    if domain is None:
        return {**entry, "status": "not_prokaryote",
                "message": f"{taxon['organism_name']} is not a bacterium or archaeon; the atlas covers only those."}
    if not is_species_level(taxon):
        return {**entry, "status": "not_species", "relatives": relatives(taxon["tax_id"]),
                "message": f"{taxon['organism_name']} is not a single species. Pick a species."}

    entry["inventory"] = inventory(taxon)
    try:  # optional source: a BacDive outage must not lose the genome profile
        entry["phenotype"] = phenotype(taxon["organism_name"])
    except Exception as e:
        entry["phenotype"] = None
        print(f"{'warning':15} BacDive lookup failed for {query}: {e}")
    best = pick_assembly(taxon["tax_id"])
    if best is None:
        parent = taxon["lineage"][-1] if taxon.get("lineage") else None
        return {**entry, "status": "no_genome", "relatives": relatives(parent) if parent else [],
                "message": "No public genome assembly for this taxon yet. Listed relatives are close taxa "
                           "that have one."}

    info, stats = best["assembly_info"], best["assembly_stats"]
    level = info.get("assembly_level")
    entry["genome"] = {
        "accession": best["accession"],
        "organism": best["organism"]["organism_name"],
        "level": level,
        "reference": info.get("refseq_category") == "reference genome",
        "size_bp": int(stats.get("total_sequence_length", 0)),
        "gc_percent": stats.get("gc_percent"),
        "contigs": stats.get("number_of_contigs"),
        "replicons": replicons(best["accession"]) if level in ("Complete Genome", "Chromosome") else [],
        "annotated": "annotation_info" in best,
    }
    records = proteome(best["accession"]) if entry["genome"]["annotated"] else []
    if not records:
        return {**entry, "status": "genome_only",
                "message": "A genome exists but has no protein annotation, so only genome size and GC are shown."}
    entry["proteome"] = composition(records)
    entry["markers"] = find_markers(records, PANELS, HERE / "markers.faa")
    return {**entry, "status": "ok", "message": "Genome and proteome profiled."}


# ---------------------------------------------------------------- library
def slug(text):
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


def rebuild_index():
    rows = []
    for path in sorted(LIBRARY.glob("*.json")):
        if path.name == "index.json" or path.name.startswith("_"):
            continue
        e = json.loads(path.read_text(encoding="utf-8"))
        rows.append({
            "file": path.name,
            "name": e.get("taxon", {}).get("name", e["query"]),
            "query": e["query"],
            "status": e["status"],
            "domain": e.get("taxon", {}).get("domain"),
            "gc_percent": (e.get("genome") or {}).get("gc_percent"),
            "predicted_ogt_c": (e.get("proteome") or {}).get("predicted_ogt_c"),
            "acidic_excess": (e.get("proteome") or {}).get("acidic_excess"),
            "measured_optimum_c": (e.get("phenotype") or {}).get("optimum_temp_c"),
            "measured_growth_c": (e.get("phenotype") or {}).get("growth_temp_c"),
            "max_nacl_percent": (e.get("phenotype") or {}).get("max_nacl_percent"),
            "markers": {k: f"{v['found']}/{v['of']}" for k, v in e.get("markers", {}).items()},
            "genes": {g["name"]: (round(g["identity"]) if g["found"] else None)
                      for v in e.get("markers", {}).values() for g in v["genes"]},
            "weak": {g["name"]: round(g["identity"]) for v in e.get("markers", {}).values()
                     for g in v["genes"] if not g["found"] and "identity" in g},
        })
    write_json(LIBRARY / "_validation.json", score_library(rows))
    write_json(LIBRARY / "index.json", rows)


def write_json(path, data):
    """Write via a temp file so an interrupted run never leaves half a file behind."""
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    tmp.replace(path)


def main(names):
    if not names:
        sys.exit(__doc__)
    LIBRARY.mkdir(exist_ok=True)
    try:
        for name in (n.strip() for n in names if n.strip()):
            try:
                entry = profile(name)
            except Exception as e:  # one failed lookup must not stop the rest
                print(f"{'error':15} {name}: {e}")
                continue
            if entry["status"] not in ("ok", "genome_only", "no_genome"):  # not a species: report, don't store
                print(f"{entry['status']:15} {entry['message']} {[r['name'] for r in entry.get('relatives', [])]}")
                continue
            path = LIBRARY / f"{slug(entry['taxon']['name'])}.json"
            write_json(path, entry)
            p = entry.get("proteome") or {}
            print(f"{entry['status']:15} {path.name:45} OGT~{p.get('predicted_ogt_c', '-')}C  "
                  f"acidic {p.get('acidic_excess', '-')}  " +
                  " ".join(f"{k} {v['found']}/{v['of']}" for k, v in entry.get("markers", {}).items()))
    finally:
        rebuild_index()


if __name__ == "__main__":
    main(sys.argv[1:])
