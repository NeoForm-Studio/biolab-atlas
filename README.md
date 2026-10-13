# BioLab Atlas — the genome that rebuilds itself

An interactive, cinematic web explorer for the DNA-repair machinery of
*Deinococcus radiodurans*, the bacterium that survives radiation doses lethal to
almost all other life. Scroll through a volumetric 3D cell and explore the genes
behind its self-repairing genome — with cross-species comparison, expression
data, and plain-language explanations of what each repair factor does.

## Live site

<https://biolab-atlas-deploy.vercel.app/>

## Screenshots

The cinematic hero — a *D. radiodurans* chromosome caught mid-repair:

![Hero — the repair engine](media/01-hero.png)

The BioLab Workspace — genes, conservation, expression, and literature side by side:

![Workspace dashboard](media/02-workspace.png)

The cross-species tool — comparing a repair gene between two organisms' real sequences:

![Cross-species comparison](media/03-cross-species.png)

## Running locally

It's a fully static site — no build step, no dependencies. Serve the folder with
any static server, for example:

```bash
python -m http.server 8000
# then open http://localhost:8000
```

Opening `index.html` directly also works, though a local server is recommended so
`data.json` loads over HTTP.

## What's in here

| Path | What it is |
|------|-----------|
| `index.html` | The full experience (scene, scroll spine, gene stations) |
| `data.json` | Gene atlas — sequences, loci, cross-species identity, expression |
| `assets/` | Rendered background plates and stills |
| `tools/genplate.mjs` | Helper that generated the background plates (needs your own `OPENAI_API_KEY`) |
| `RECIPE-cinematic-volumetric-3d-scene.md` | Build notes for the 3D scene |

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
   NCBI Nucleotide, and Europe PMC papers mentioning the name.
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

## Credits

Gene and sequence references link out to [PubMed](https://pubmed.ncbi.nlm.nih.gov/)
and [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/). Rendering uses
[Three.js](https://threejs.org/) and [GSAP](https://greensock.com/gsap/).
