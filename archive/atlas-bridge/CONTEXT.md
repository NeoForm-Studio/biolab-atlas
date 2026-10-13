# BioLab Atlas — build context & handoff

A self-contained WebGL explorer of DNA-repair genes across radiation-resistant bacteria
(*D. radiodurans*, *D. deserti*, *D. gobiensis*) vs the sensitive control *E. coli*. Every
number in it is derived from real public data — RefSeq genomes, a public GEO expression
series, and published survival phenotypes. **Claude Code owns the build from here.**

---

## File map
```
atlas/
  index.html              # BUILT OUTPUT — this is the deployable app. Do not hand-edit.
  vercel.json             # static config (cleanUrls)
  CONTEXT.md              # this file
  .gitignore              # ignores build/cache/ (raw FASTAs are re-downloadable)
  build/
    atlas_meta.py         # CURATED KNOWLEDGE — edit this to add genomes/genes/expression
    build_atlas.py        # fetch RefSeq + compute GC/identity + inject -> ../index.html
    atlas_template.html   # UI/template (Three.js helix, panels, EXPL text, Compute mode)
    genes_data.js         # GENERATED data (committed for reference; rebuilt each run)
    cache/                # downloaded FASTAs — gitignored, disposable
  server/                 # LOCAL pipeline bridge (localhost only) — powers the Compute tab
    atlas_server.py       # FastAPI app on 127.0.0.1:8765
    pipeline_adapter.py   # Tier A real homology + Tier B fail-loud seams to BioLab
    requirements.txt
    README.md
```

**Source of truth split:** `atlas_meta.py` = the curated facts (which genomes, which genes,
plain-language descriptions, the real expression numbers, the survival phenotypes).
`atlas_template.html` = the UI and the lay-explanation text (the `EXPL` and `PHENO` blocks in
its `<script>`). `build_atlas.py` = the deterministic step that turns RefSeq + meta into the page.

---

## Build & backup workflow  ← do this after EVERY change
```bash
cd atlas/build
python build_atlas.py            # rebuilds ../index.html from current meta + RefSeq
# eyeball ../index.html in a browser, then back it up:
cd ..
git add atlas/ ; git commit -m "atlas: <what changed>" ; git push
```
- Rebuild before committing — never hand-edit `index.html`; it will be overwritten.
- Commit `build/*.py`, `build/genes_data.js`, `atlas_template.html`, and `index.html`.
- Do **not** commit `build/cache/` (gitignored — FASTAs re-download in minutes).
- For anything experimental, branch first: `git checkout -b atlas-<feature>`, build, commit,
  open a PR / merge when it checks out. The repo is the backup; commit small and often.

Deps: `pip install biopython`. The build needs internet (NCBI RefSeq).

---

## INTEGRITY RULES — non-negotiable
This tool's whole value is that it never fabricates. Preserve these or it's worthless:

1. **No invented data.** Sequences, GC%, protein length, identity, and presence are computed
   from RefSeq by `build_atlas.py`. Never type these in by hand.
2. **Expression is real or absent.** `expr` in `atlas_meta.py` holds DESeq2 output from a real
   GEO series (currently DR1 only, GSE176207, 10 kGy IR, 3v3). Only DR1 has it. **Never add an
   `expr` value that isn't from a real differential-expression run on real reads.** If a second
   organism gets radiation RNA-seq, add the genome, run DE, then paste the real numbers.
3. **Presence is "by annotated symbol," a proxy — not proof of absence.** A gene showing in 1/4
   genomes may simply be unannotated or renamed elsewhere. The UI says this; keep it saying it.
   Definitive orthology comes from the BioLab pipeline's orthogroups, not from this tool.
4. **Identity is a global amino-acid alignment** (BLOSUM62), not residue-by-residue curation.
   High identity ⇒ likely same function; it does not prove identical kinetics.
5. **Phenotype is cited literature, not estimates** (see provenance below).
6. **The honest ceiling, kept visible:** there is *no* evidence *D. radiodurans* makes more RecA
   or PolA than *E. coli*. Its resistance is multi-factorial — extra Deinococcus genes, multiple
   genome copies, the IrrE switch, and manganese-antioxidant proteome protection. The tool must
   never imply "more/faster repair protein" explains the difference. (This lives in `EXPL.tells`.)

---

## How to extend

**Add a genome** → in `atlas_meta.py`, add an entry to `SPECIES` with its RefSeq `acc`
(`GCF_...`), and a `PHENO` entry with a *cited* survival dose. Rebuild. It auto-joins the
species switcher, cross-species picker, presence matrix, and identity matrix.

**Add real expression** → only after a real DE run. Put `{ "log2fc":…, "padj":…, "baseMean":…,
"env":"…", "geo":"GSE…", "n":"3v3" }` into `META[gene]["expr"]`. Rebuild. (If it's a second
organism, the template's Expression view currently renders DR1; extending it to multi-species
is the next UI task — and the same real-data bar applies.)

**Add an anchor gene** → append the symbol to `ANCHORS` and add a `META[symbol]` block
(`role, scope, pathway, plain, why, tags, expr`). Add a `SYN` entry if it has alias symbols.
Rebuild; identity/GC/presence fill in automatically.

**Edit a lay explanation** → the expandable `＋` texts are the `EXPL` object in
`atlas_template.html`; per-species survival blurbs are `PHENO` there too. Edit, rebuild.

---

## Data provenance (cite when extending)
- **Genomes (RefSeq), 12 across 4 lineages + 2 sensitive controls** — all verified fetchable:
  - *Deinococcus–Thermus:* DR1 `GCF_000008565.1` · DDES `GCF_000020685.1` · DGOB `GCF_000252445.1` ·
    DGEO (D. geothermalis) `GCF_000196275.1` · DPRO (D. proteolyticus) `GCF_000190555.1` · **TTH** (T. thermophilus, *sensitive control*) `GCF_000091545.1`
  - *Actinobacteria:* RXY (R. xylanophilus) `GCF_019448335.1` · KRAD (K. radiotolerans) `GCF_000017305.1`
  - *Cyanobacteria:* CHRO (Chroococcidiopsis thermalis) `GCF_000317125.1`
  - *Bacteroidetes:* HYM (Hymenobacter yonginensis) `GCF_027625995.1`
  - *Euryarchaeota (Archaea, flagged `distant`):* TGAM (T. gammatolerans) `GCF_000022365.1` ·
    HALO (H. salinarum NRC-1) `GCF_000006805.1` · PFUR (P. furiosus) `GCF_008245085.1`
  - *Proteobacteria:* EC (E. coli K-12, *sensitive*) `GCF_000005845.2` · SHEW (S. oneidensis MR-1, *sensitive, high-Fe*) `GCF_000146165.2`
  - CDS pulled from `*_cds_from_genomic.fna.gz`. **Archaea undercount by symbol** (they use radA not recA,
    etc.) — flagged `distant:True` in SPECIES and caveated in the UI. This is expected, not a bug.

- **Expression:** *D. radiodurans* GSE176207 (ionizing radiation, 3 vs 3), DESeq2 on raw counts — DR1 only.

### Expression-data landscape (GEO scan, for adding a real second-organism layer)
Real radiation transcriptomes that EXIST and could be run through DE (Rule 2 still applies — only add
after a real DE run; most below lack biological replicates, which weakens DESeq2 + the convergence gate):
  - **DDES — GSE56058** (RNA-seq, IR) — second Deinococcus, replicate-able. Best candidate.
  - **KRAD — GSE51096** (RNA-seq, IR) — cross-*phylum* radiation RNA-seq (Actinobacteria); n=2, unreplicated.
  - **DGOB — GSE29088** (UV RNA-seq, n=2) · **PFUR — GSE5919 / GDS2469** (gamma, microarray not RNA-seq).
  - **DR1 — GSE175875** (Radiation/Desiccation Response regulon, n=41) — more DR1 depth.
  - **HALO** — large oxidative-stress transcriptome catalog (oxidative damage = radiation's proximal hit).
Gaps (no public radiation expression): **RXY, CHRO**, and notably **TGAM** (only a cadmium study) — these
are where the lab's own RNA-seq would be novel.

## Survival phenotype sources (per PHENO entry)
- DR1 5,000 Gy ~no loss / 10% at 12,000 Gy (PNAS 1009312107; mednexus radmp.2023.03.001; Wikipedia).
- DDES ~15 kGy (IJSEM 16280508) · DGOB >15 kGy, >600 J/m² UV (IJSEM 19502345) · DGEO multi-kGy (Deinococcus).
- TTH radiation-SENSITIVE, few hundred Gy (FEMS Microbiol Rev 2018) — in-phylum negative control.
- RXY >25 kGy (AEM 71/9/5225) · KRAD D10 10–12 kGy (Frontiers fmicb.2026.1777693; Phillips 2002).
- CHRO ~15 kGy (Billi 2000) · TGAM up to 30 kGy, D10 ~7.5 kGy (Wikipedia; Tapias 2009) ·
  HALO D10 ~5 kGy (PMC3505630) · PFUR multi-kGy (PMC3505630).
- EC ~150 Gy D10, 200–800 Gy lethal; 5–10 Gy lethal to humans (FEMS 2018).
- **Honest ceiling:** no elevated RecA/PolA in DR vs EC; resistance from proteome protection — PNAS 2010.

---

## Deploy
Static. From `atlas/`: `vercel --prod` (CLI), or import the repo in the Vercel dashboard
(framework: Other). `index.html` serves at root. Needs internet at view time for the Three.js
and Google-Fonts CDNs; ask if you want those vendored for offline use.

## Pipeline bridge (the Compute tab) — server/
A local backend that makes the Atlas a live front-end for the pipeline. **Localhost only**;
nothing is exposed. Run: `cd server; pip install -r requirements.txt; python atlas_server.py`,
then open the Atlas and click **Compute**.

- **Tier A also includes a convergence sift** (`/api/convergence_preview`): across all genomes,
  it labels each anchor gene universal / lineage-specific / convergent-candidate by independent
  radioresistant lineage. Honest finding it surfaces: the core repair toolkit (recA, uvrAB,
  recN, gyrA, ruvB, polA) is UNIVERSAL — shared by the sensitive controls too — so having it
  isn't the explanation; candidates (ddrC, ddrB) are leads to verify with Tier-B orthogroups.
- **Tier A (live now):** reciprocal-best-hit homology over the 12-genome set + differential
  presence. Real compute, honest confidence (`high`/`weak`/`none`), explicitly a PROXY. This is
  what recovers archaeal radA/radB (e.g. T. gammatolerans recA → radB, RBH, ~36%) that the
  by-symbol presence in the static atlas shows as absent. It defers to Tier-B orthogroups.
- **Tier B (wire to BioLab):** `orthogroups` / `expression_de` / `convergence_candidates` in
  `pipeline_adapter.py` read the pipeline's real output (default `../pipeline/output`, override
  with env `BIOLAB_OUT`). If absent they raise `PipelineUnavailable` → the API returns 404 and
  the UI says "not wired yet" — they NEVER fabricate. **To finish the connection:** point
  `BIOLAB_OUT` at the pipeline's output dir (OrthoFinder `Orthogroups.tsv`, DE CSVs,
  `convergence/candidates.json`), or replace those three adapter functions with direct calls
  into `orchestrate.py`. Keep the fail-loud contract: real data or an honest error, nothing invented.
