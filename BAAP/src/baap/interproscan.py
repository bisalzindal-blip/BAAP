"""InterProScan REST API integration module."""

import time
import logging
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO

from .config import (
    RESULTS_DIR, EMAIL, BATCH_SIZE, SUBMISSION_SLEEP,
    STATUS_SLEEP, MAX_WAIT_MINUTES, log
)

BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"

PROTEASE_ANALYSES = {
    "Pfam", "SMART", "CDD", "SUPERFAMILY", "Gene3D",
    "PANTHER", "ProSiteProfiles", "ProSitePatterns",
    "PRINTS", "TIGRFAM", "HAMAP"
}

IPR_COLUMNS = [
    "Protein_ID", "Sequence_MD5", "Sequence_Length",
    "Analysis", "Signature_Accession", "Signature_Description",
    "Start", "End", "Score", "Status", "Date",
    "InterPro_Accession", "InterPro_Description",
    "GO_Annotations", "Pathways"
]


def run_interproscan_rest(candidate_records):
    """Run InterProScan via REST API."""
    if EMAIL.strip() == "" or "YOUR_REAL_EMAIL" in EMAIL:
        raise RuntimeError("Set a real email address in EMAIL before running InterProScan.")
    
    if not candidate_records:
        raise RuntimeError("No candidates supplied to InterProScan.")
    
    batches = [candidate_records[i:i+BATCH_SIZE] for i in range(0, len(candidate_records), BATCH_SIZE)]
    job_ids = []
    
    for batch_num, batch in enumerate(batches, start=1):
        fasta_text = ""
        for rec in batch:
            fasta_text += f">{rec.id}\n{str(rec.seq)}\n"
        
        payload = {
            "email": EMAIL,
            "sequence": fasta_text,
            "goterms": "true",
            "pathways": "true",
        }
        
        log.info("Submitting InterProScan batch %d/%d (%d proteins)",
                 batch_num, len(batches), len(batch))
        
        try:
            r = requests.post(f"{BASE_URL}/run", data=payload, timeout=180)
        except requests.exceptions.RequestException as e:
            log.error("InterProScan submission failed: %s", e)
            continue
        
        if r.status_code != 200:
            log.error("InterProScan submission failed: %s", r.text[:1000])
            continue
        
        job_id = r.text.strip()
        job_ids.append((batch_num, job_id))
        
        if batch_num < len(batches):
            time.sleep(SUBMISSION_SLEEP)
    
    if not job_ids:
        raise RuntimeError("No InterProScan jobs were successfully submitted.")
    
    completed = {}
    start_time = time.time()
    max_wait_seconds = MAX_WAIT_MINUTES * 60
    
    while True:
        all_done = True
        for batch_num, job_id in job_ids:
            if job_id in completed:
                continue
            try:
                r = requests.get(f"{BASE_URL}/status/{job_id}", timeout=60)
                status = r.text.strip()
            except requests.exceptions.RequestException:
                all_done = False
                continue
            
            log.info("Batch %d job %s status: %s", batch_num, job_id, status)
            
            if status == "FINISHED":
                completed[job_id] = True
            elif status in {"FAILURE", "ERROR"}:
                completed[job_id] = False
            else:
                all_done = False
        
        if all_done and len(completed) == len(job_ids):
            break
        
        if time.time() - start_time > max_wait_seconds:
            raise TimeoutError(f"InterProScan jobs exceeded {MAX_WAIT_MINUTES} minutes.")
        
        time.sleep(STATUS_SLEEP)
    
    # Download results
    all_tables = []
    for batch_num, job_id in job_ids:
        if not completed.get(job_id, False):
            continue
        
        try:
            result = requests.get(f"{BASE_URL}/result/{job_id}/tsv", timeout=120)
            result.raise_for_status()
            
            batch_file = RESULTS_DIR / f"interpro_batch_{batch_num}.tsv"
            batch_file.write_text(result.text)
            
            df = pd.read_csv(
                batch_file, sep="\t", comment="#", header=None,
                names=IPR_COLUMNS, dtype=str
            )
            all_tables.append(df)
        except Exception as e:
            log.warning("Could not download batch %d: %s", batch_num, e)
    
    if not all_tables:
        raise RuntimeError("No InterProScan results were successfully downloaded.")
    
    merged = pd.concat(all_tables, ignore_index=True)
    merged.to_csv(RESULTS_DIR / "interproscan.tsv", sep="\t", index=False)
    
    return merged