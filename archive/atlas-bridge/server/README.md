# Atlas ↔ pipeline bridge (localhost only)

Turns the Atlas "Compute" tab into a live front-end for the BioLab pipeline. Nothing is
exposed publicly — it binds to 127.0.0.1.

## Run
```bash
cd atlas/server
pip install -r requirements.txt
python atlas_server.py          # http://localhost:8765
```
Then open `atlas/index.html` and click the **Compute** tab. The status dot goes green when the
bridge is reachable; if the server isn't running, the rest of the atlas still works.

## What it computes
- **Tier A (works now, on the 12-genome set):** reciprocal-best-hit homology (`/api/ortholog`)
  and differential presence (`/api/differential`). Real compute; honest confidence labels;
  a PROXY that defers to Tier-B orthogroups. This is what recovers archaeal radA/radB that
  symbol-matching shows as absent.
- **Tier B (the heavy pipeline):** `/api/orthogroups`, `/api/expression`, `/api/convergence`
  read the pipeline's real output. If it isn't generated, they return 404 with the path to
  wire — they never fabricate. Point `BIOLAB_OUT` at the pipeline output dir, or wire
  `pipeline_adapter.py`'s Tier-B functions to `orchestrate.py`.

## Integrity
Every endpoint returns really-computed data or an honest error. No invented orthogroups,
expression, or phenotype — same rules as the rest of the project (see ../CONTEXT.md).
