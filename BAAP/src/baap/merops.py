"""MEROPS database and DIAMOND screening module."""

import re
import subprocess
import logging
from pathlib import Path

import pandas as pd

from .config import (
    DIRS, RESULTS_DIR, EVALUE_CUTOFF, IDENTITY_CUTOFF,
    SUBJECT_COVERAGE_CUTOFF, QUERY_COVERAGE_CUTOFF,
    MAX_TARGET_SEQS, CPU_THREADS, stage, log
)

# MEROPS paths
MEROPS_DIR = DIRS["merops"]
MEROPS_RAW = MEROPS_DIR / "pepunit.lib"
MEROPS_FASTA = MEROPS_DIR / "merops_cleaned.fasta"
MEROPS_DB = MEROPS_DIR / "merops_db"
MEROPS_DMND = Path(str(MEROPS_DB) + ".dmnd")
DIAMOND_RESULTS = MEROPS_DIR / "diamond_results.tsv"

MEROPS_URL = "https://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib"

# Output paths
MEROPS_PASS = RESULTS_DIR / "MEROPS_PASS_hits.csv"
MEROPS_FAILED = RESULTS_DIR / "MEROPS_failed_hits.csv"
MEROPS_NOHIT = RESULTS_DIR / "MEROPS_no_hits.csv"
MEROPS_SUMMARY = RESULTS_DIR / "MEROPS_summary.csv"

DIAMOND_COLUMNS = [
    "Query_ID", "Subject_ID", "Subject_Title", "Percent_Identity",
    "Alignment_Length", "Mismatches", "Gap_Openings", "Qstart", "Qend",
    "Sstart", "Send", "Evalue", "Bitscore", "Query_Length", "Subject_Length"
]


def download_and_prepare_merops():
    """Download and prepare MEROPS database."""
    stage("PREPARING MEROPS")
    
    if not (MEROPS_RAW.exists() and MEROPS_RAW.stat().st_size > 0):
        subprocess.run([
            "wget", "-q", "--show-progress",
            "-O", str(MEROPS_RAW), MEROPS_URL
        ], check=True)
    
    if not MEROPS_FASTA.exists() or MEROPS_FASTA.stat().st_size == 0:
        count = 0
        with open(MEROPS_RAW, "r", encoding="latin-1") as infile, \
             open(MEROPS_FASTA, "w") as outfile:
            for record in SeqIO.parse(infile, "fasta"):
                record.seq = record.seq.upper()
                SeqIO.write(record, outfile, "fasta")
                count += 1
        if count == 0:
            raise RuntimeError("MEROPS database contains no FASTA sequences.")
    
    if not (MEROPS_DMND.exists() and MEROPS_DMND.stat().st_size > 0):
        subprocess.run([
            "diamond", "makedb",
            "--in", str(MEROPS_FASTA),
            "--db", str(MEROPS_DB)
        ], check=True)


def run_diamond(query_faa):
    """Run DIAMOND BLASTP against MEROPS."""
    cmd = [
        "diamond", "blastp",
        "--query", str(query_faa),
        "--db", str(MEROPS_DMND),
        "--out", str(DIAMOND_RESULTS),
        "--outfmt", "6",
        "qseqid", "sseqid", "stitle", "pident", "length",
        "mismatch", "gapopen", "qstart", "qend", "sstart",
        "send", "evalue", "bitscore", "qlen", "slen",
        "--evalue", str(EVALUE_CUTOFF),
        "--max-target-seqs", str(MAX_TARGET_SEQS),
        "--threads", str(CPU_THREADS),
    ]
    subprocess.run(cmd, check=True)


def extract_merops_family(title):
    """Extract MEROPS family from subject title."""
    if pd.isna(title):
        return "Unknown"
    t = str(title)
    m = re.search(r"\[([A-Z]\d{1,3})(?:\.\d+)?\]", t)
    if m:
        return m.group(1)
    m = re.search(r"\[([A-Z]\d{1,3})", t)
    return m.group(1) if m else "Unknown"


def extract_merops_annotation(title):
    """Extract annotation from MEROPS title."""
    if pd.isna(title):
        return "NA"
    t = str(title).strip()
    t = re.split(r"\[[A-Z]\d{1,3}(?:\.\d+)?\]", t, maxsplit=1)[0].strip()
    t = re.sub(r"^MER\d+\s*-\s*", "", t)
    return t.strip() or "NA"


def classify_merops_hit(row):
    """Classify a MEROPS hit."""
    if pd.isna(row["Evalue"]):
        return "HIGH_EVALUE"
    if row["Evalue"] > EVALUE_CUTOFF:
        return "HIGH_EVALUE"
    if row["Percent_Identity"] < IDENTITY_CUTOFF:
        return "LOW_IDENTITY"
    if row["Subject_Coverage"] < SUBJECT_COVERAGE_CUTOFF:
        return "LOW_SUBJECT_COVERAGE"
    if row["Query_Coverage"] < QUERY_COVERAGE_CUTOFF:
        return "LOW_QUERY_COVERAGE"
    return "PASS"


def screen_merops(prokka_faa, prokka_annotation_map):
    """Run MEROPS screening pipeline."""
    stage("MEROPS PROTEASE SCREENING")
    
    download_and_prepare_merops()
    run_diamond(prokka_faa)
    
    if not DIAMOND_RESULTS.exists() or DIAMOND_RESULTS.stat().st_size == 0:
        log.warning("No DIAMOND hits against MEROPS.")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    
    df = pd.read_csv(DIAMOND_RESULTS, sep="\t", names=DIAMOND_COLUMNS)
    
    # Convert numeric columns
    numeric_cols = [
        "Percent_Identity", "Alignment_Length", "Mismatches",
        "Gap_Openings", "Qstart", "Qend", "Sstart", "Send",
        "Evalue", "Bitscore", "Query_Length", "Subject_Length"
    ]
    for c in numeric_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    
    # Calculate coverage
    df["Query_Coverage"] = ((df["Qend"] - df["Qstart"] + 1) / df["Query_Length"]) * 100
    df["Subject_Coverage"] = ((df["Send"] - df["Sstart"] + 1) / df["Subject_Length"]) * 100
    
    # Add annotations
    df["MEROPS_Family"] = df["Subject_Title"].apply(extract_merops_family)
    df["MEROPS_Annotation"] = df["Subject_Title"].apply(extract_merops_annotation)
    df["Prokka_Annotation"] = df["Query_ID"].map(prokka_annotation_map).fillna("")
    df["Status"] = df.apply(classify_merops_hit, axis=1)
    
    df.to_csv(RESULTS_DIR / "MEROPS_all_hits.csv", index=False)
    
    # Best hits per query
    best_hits = (
        df.sort_values(
            by=["Query_ID", "Bitscore", "Evalue", "Percent_Identity"],
            ascending=[True, False, True, False]
        )
        .drop_duplicates(subset=["Query_ID"], keep="first")
        .copy()
    )
    
    passed = best_hits[best_hits["Status"] == "PASS"].sort_values("Bitscore", ascending=False).copy()
    failed = best_hits[best_hits["Status"] != "PASS"].sort_values("Bitscore", ascending=False).copy()
    
    # No-hit proteins
    all_query_ids = set(prokka_annotation_map.keys())
    hit_ids = set(best_hits["Query_ID"])
    no_hit_ids = all_query_ids - hit_ids
    
    nohit_rows = []
    for pid in sorted(no_hit_ids):
        nohit_rows.append({
            "Protein_ID": pid,
            "Prokka_Annotation": prokka_annotation_map.get(pid, ""),
            "Status": "NO_MEROPS_HIT"
        })
    nohit_df = pd.DataFrame(nohit_rows)
    
    # Save outputs
    passed.to_csv(MEROPS_PASS, index=False)
    failed.to_csv(MEROPS_FAILED, index=False)
    nohit_df.to_csv(MEROPS_NOHIT, index=False)
    
    # Summary
    summary = pd.DataFrame({
        "Metric": [
            "Total proteins",
            "Proteins with MEROPS hit",
            "MEROPS PASS",
            "MEROPS failed",
            "No MEROPS hit"
        ],
        "Count": [
            len(all_query_ids),
            len(hit_ids),
            len(passed),
            len(failed),
            len(no_hit_ids)
        ]
    })
    summary.to_csv(MEROPS_SUMMARY, index=False)
    
    log.info("MEROPS: %d PASS | %d failed | %d no-hit",
             len(passed), len(failed), len(no_hit_ids))
    
    return df, passed, failed