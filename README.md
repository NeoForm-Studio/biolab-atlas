# BioLab Atlas

Look up any bacterium or archaeon and see what public data exists for it and how much its genome resembles
known extremophiles. It opens with a 3D scene of *Deinococcus radiodurans*, a bacterium that survives
γ-radiation doses up to 10,000 Gy and rebuilds a genome shattered into hundreds of fragments.

![D. radiodurans scene](media/01-hero.png)

## Running locally

A static site with no build step. Serve the folder and open it:

```bash
python -m http.server 8000
# http://localhost:8000              the 3D D. radiodurans scene
# http://localhost:8000/explorer.html  the species explorer
```

## What's in here

| Path | What it is |
|------|-----------|
| `index.html` | 3D introduction: D. radiodurans genome repair |
| `explorer.html` | Species explorer: search the library, compare species |
| `library/` | One JSON profile per species, plus `index.json` |
| `pipeline/` | The lookup pipeline that builds the library |
| `assets/` | Background plates for the 3D scene, AI-generated with `tools/genplate.mjs` (needs your own `OPENAI_API_KEY`) |

## Species lookup pipeline

`pipeline/` looks up any bacterium or archaeon in public databases and profiles its genome. It needs Python 3.10+
(standard library only) and [NCBI BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK569861/) (`blastp` on PATH).
`explorer.html` shows the results: open it through a local server (`python -m http.server`), search a species, and
see its public data, genome, marker genes and where it sits among every other profiled species.

```bash
python pipeline/atlas.py "Thermus thermophilus" "Halomonas elongata"
python -m unittest discover -s pipeline
```

For each name it:

1. **Resolves** the name in NCBI Taxonomy. Unknown names, non-prokaryotes and ranks above species are reported and
   not stored. A genus lists member species that have reference genomes.
2. **Inventories** public data: NCBI genome assemblies, UniProt proteins (all and reviewed), 16S rRNA records in
   NCBI Nucleotide, Europe PMC papers mentioning the name, and measured growth data from
   [BacDive](https://bacdive.dsmz.de/) (type strain preferred: growth temperatures, optimum, NaCl tolerance).
3. **Picks a genome**: NCBI's reference genome if there is one, otherwise the most complete annotated assembly.
   Name a strain (for example `"Escherichia coli K-12"`) to get that strain.
4. **Profiles it**: replicons, size, GC, and from the annotated proteome:
   - `predicted_ogt_c`: optimal growth temperature estimated from the IVYWREL amino-acid fraction,
     T = 937F − 335 (Zeldovich et al. 2007, PLoS Comput Biol 3:e5, r = 0.93 over 204 proteomes). Treat it as
     a rough estimate: in our tests it was close for Pyrococcus furiosus (98 °C vs about 100 °C) but about 30 °C
     too high for Thermus thermophilus (about 70 °C) and about 20 °C too high for Halomonas elongata.
   - `acidic_excess`: (D + E − K − R) per residue. Salt-in halophiles stand out (Halobacterium 0.079 against
     0.00–0.03 for the others tested).
   - Trait marker panels (`pipeline/markers.json`): Deinococcus-type DNA-damage-response proteins, reverse gyrase
     (all known hyperthermophiles and some thermophiles), ectoine synthesis (salt-out osmoprotection). Each
     reference protein (`pipeline/markers.faa`, from NCBI RefSeq) is searched against the proteome with `blastp`.
     It counts as present at e-value ≤ 1e-10, ≥ 70% coverage of the reference and ≥ 30% identity. EctB needs 45%,
     because ordinary class III aminotransferases match it at 27–39% across nearly its full length, while
     distant ectoine producers (Streptomyces, Virgibacillus, Marinococcus) sit at 50–54%. PprM was left
     out: it is related to cold-shock proteins and matches mesophiles at 51–69% identity. A missing marker is not
     evidence the trait is absent.
5. **Saves** `library/<species>.json` and rebuilds `library/index.json`. A species with no genome is saved as
   `no_genome`, with relatives that do have one.

These are signals that a genome resembles known extremophiles, not a determination that an organism is one.

## Sources and credits

- Data: [NCBI](https://www.ncbi.nlm.nih.gov/) Taxonomy, Datasets and Nucleotide; [UniProt](https://www.uniprot.org/);
  [Europe PMC](https://europepmc.org/); [BacDive](https://bacdive.dsmz.de/) (DSMZ). Marker search: NCBI BLAST+.
- Growth-temperature estimate: Zeldovich, Berezovsky & Shakhnovich 2007, PLoS Comput Biol 3:e5.
- Reverse gyrase in hyperthermophiles: Forterre 2002, Trends Genet 18:236.
- D. radiodurans claims on the 3D page: genome reassembly from hundreds of fragments, Zahradka et al. 2006,
  Nature 443:569 (doi:10.1038/nature05160); γ-radiation resistance up to 10 kGy, Xie et al. 2026, Appl Environ
  Microbiol (doi:10.1128/aem.02107-25).
- Rendering: [Three.js](https://threejs.org/).
