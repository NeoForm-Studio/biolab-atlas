#!/usr/bin/env python3
"""
atlas_server.py  —  local bridge between the Atlas UI and the BioLab pipeline.

Runs ONLY on localhost. Nothing is exposed publicly. Start it, open the Atlas, and the
"Compute" tab calls these endpoints. If this server isn't running, the Atlas still works as a
static explorer — it just shows the bridge as offline.

  pip install -r requirements.txt
  python atlas_server.py            # serves http://localhost:8765

Endpoints (all return REAL computed/looked-up data, or an honest error — never fabrication):
  GET  /api/health                      what's live (Tier A always; Tier B if pipeline output exists)
  POST /api/ortholog   {gene}           reciprocal-best-hit homology across all genomes
  POST /api/differential {present_in, absent_in}   sift candidate lineage-distinguishing genes
  POST /api/orthogroups {genes?}        Tier B: pipeline orthogroups (404 if not generated)
  POST /api/expression {sp, geo}        Tier B: a real DESeq2 table (404 if not present)
  GET  /api/convergence                 Tier B: convergence candidate scores (404 if not present)
"""
import uvicorn
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional
import pipeline_adapter as PA

app = FastAPI(title="BioLab Atlas — pipeline bridge", version="1.0")
# localhost only — the Atlas may be opened via file:// or a local static server
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

@app.get("/api/health")
def health():
    return {"ok": True, "service": "biolab-atlas-bridge", **PA.capabilities()}

class GeneIn(BaseModel):
    gene: str

@app.post("/api/ortholog")
def ortholog(q: GeneIn):
    try:
        return PA.rbh_ortholog(q.gene)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

class DiffIn(BaseModel):
    present_in: List[str]
    absent_in: List[str]

@app.post("/api/differential")
def differential(q: DiffIn):
    if not q.present_in:
        raise HTTPException(status_code=400, detail="present_in must list at least one genome")
    return PA.differential_presence(q.present_in, q.absent_in)

@app.get("/api/convergence_preview")
def convergence_preview():
    return PA.convergence_preview()

class OGIn(BaseModel):
    genes: Optional[List[str]] = None

@app.post("/api/orthogroups")
def orthogroups(q: OGIn):
    try:
        return PA.orthogroups(q.genes)
    except PA.PipelineUnavailable as e:
        raise HTTPException(status_code=404, detail=str(e))

class ExprIn(BaseModel):
    sp: str
    geo: str

@app.post("/api/expression")
def expression(q: ExprIn):
    try:
        return PA.expression_de(q.sp, q.geo)
    except PA.PipelineUnavailable as e:
        raise HTTPException(status_code=404, detail=str(e))

@app.get("/api/convergence")
def convergence():
    try:
        return PA.convergence_candidates()
    except PA.PipelineUnavailable as e:
        raise HTTPException(status_code=404, detail=str(e))

if __name__ == "__main__":
    print("BioLab Atlas bridge → http://localhost:8765  (localhost only)")
    print("Tier A homology is live. Tier B reads pipeline output at:", PA.capabilities()["pipeline_out"])
    uvicorn.run(app, host="127.0.0.1", port=8765, log_level="info")
