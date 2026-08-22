"""Prokka annotation module."""

import os
import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from .config import (
    DIRS, RESULTS_DIR, CPU_THREADS, ENV_NAME,
    PROKKA_BIN, CONDA, FORCE_RERUN, run_cmd
)

log = logging.getLogger("BAAP")


def prokka_outputs_valid(outdir):
    """Check if Prokka outputs are valid."""
    required = [
        outdir / "PROKKA.faa",
        outdir / "PROKKA.tsv",
        outdir / "PROKKA.gff",
        outdir / "PROKKA.gbk",
    ]
    return all(p.exists() and p.stat().st_size > 0 for p in required)


def run_prokka(genome, outdir):
    """Run Prokka annotation."""
    outdir.mkdir(parents=True, exist_ok=True)
    
    if not FORCE_RERUN and prokka_outputs_valid(outdir):
        print("Existing Prokka outputs detected — skipping.")
        return outdir
    
    db_dir = Path(os.environ.get("PROKKA_DB_DIR", ""))
    if not db_dir.exists():
        raise RuntimeError(f"Prokka database directory not found: {db_dir}")
    
    cmd = [
        str(CONDA), "run", "-n", ENV_NAME, str(PROKKA_BIN),
        "--dbdir", str(db_dir),
        str(genome),
        "--outdir", str(outdir),
        "--prefix", "PROKKA",
        "--cpus", str(CPU_THREADS),
        "--force",
    ]
    
    env_vars = os.environ.copy()
    env_vars["JAVA_TOOL_OPTIONS"] = "-Xlog:disable"
    
    run_cmd(cmd, stderr_path=DIRS["logs"] / "prokka.stderr.log", env=env_vars)
    
    return outdir


def parse_prokka(prokka_dir):
    """Parse Prokka output files."""
    faa = prokka_dir / "PROKKA.faa"
    tsv = prokka_dir / "PROKKA.tsv"
    
    if not faa.exists():
        raise FileNotFoundError(f"Missing {faa}")
    
    rows = []
    for r in SeqIO.parse(faa, "fasta"):
        rows.append({
            "Protein_ID": r.id,
            "Prokka_Accession": r.id,
            "Protein_Length": len(r.seq),
            "Protein_Sequence": str(r.seq),
        })
    
    df = pd.DataFrame(rows)
    
    if tsv.exists():
        try:
            ann = pd.read_csv(tsv, sep="\t", dtype=str)
            ann.to_csv(RESULTS_DIR / "prokka_annotations.tsv", sep="\t", index=False)
            
            product_col = next(
                (c for c in ["product", "Product", "annotation"] if c in ann.columns),
                None
            )
            
            if product_col and "locus_tag" in ann.columns:
                df = df.merge(
                    ann[["locus_tag", product_col]].rename(
                        columns={
                            "locus_tag": "Prokka_Accession",
                            product_col: "Prokka_Annotation"
                        }
                    ),
                    on="Prokka_Accession",
                    how="left",
                )
        except Exception as e:
            log.warning("Could not parse Prokka TSV: %s", e)
    
    if "Prokka_Annotation" not in df.columns:
        df["Prokka_Annotation"] = ""
    
    return df


def extract_protein_ids(faa_path):
    """Extract protein IDs from a FASTA file."""
    return [r.id for r in SeqIO.parse(faa_path, "fasta")]


def get_protein_sequences(faa_path):
    """Get protein sequences as a dict."""
    return {r.id: str(r.seq) for r in SeqIO.parse(faa_path, "fasta")}