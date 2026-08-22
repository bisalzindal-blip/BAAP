"""Main pipeline orchestration module."""

import sys
import time
import shutil
import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from . import __version__ as PIPELINE_VERSION
from . import __name__ as PIPELINE_NAME
from .genome import validate_genome, get_genome_stats
from .prokka import run_prokka, parse_prokka
from .merops import screen_merops
from .interproscan import run_interproscan_rest
from .validation import validate_domains
from .reporting import create_final_audit, create_complete_zip, save_manifest
from .config import (
    PROJECT_DIR, DIRS, RESULTS_DIR, LOG_FILE,
    PROKKA_BIN, install_prokka, FORCE_RERUN,
    stage, log, setup_logging
)

ROOT = Path(PROJECT_DIR)
PROJECT_DIR = ROOT


def run_pipeline(genome_path=None, upload_from_colab=True):
    """Execute the complete BAAP pipeline.
    
    Args:
        genome_path: Path to genome FASTA file (optional)
        upload_from_colab: If True, use Colab file upload (default: True)
    
    Returns:
        dict: Pipeline results and statistics
    """
    
    # Setup logging
    setup_logging()
    
    stage(f"{PIPELINE_NAME} v{PIPELINE_VERSION}")
    
    print("""
Scientific workflow:
  
Genome
  ↓
Prokka
  ↓
All predicted proteins
  ↓
MEROPS DIAMOND sensitive screening
  ↓
MEROPS PASS candidates
  ↓
InterProScan independent domain validation
  ↓
Protease-specific Pfam / InterPro evidence
  ↓
GO peptidase evidence
  ↓
Positive / negative / inactive keyword context
  ↓
Protein-level evidence integration
  ↓
TRUE / PUTATIVE / INACTIVE / NON-PROTEASE / REVIEW
""")
    
    # Install Prokka if needed
    if not PROKKA_BIN.exists():
        install_prokka()
    
    # Handle genome input
    if upload_from_colab:
        from google.colab import files
        stage("UPLOAD GENOME FASTA")
        uploaded = files.upload()
        
        fasta_candidates = [
            Path(name) for name in uploaded
            if Path(name).suffix.lower() in {".fasta", ".fa", ".fna"}
        ]
        
        if not fasta_candidates:
            raise ValueError("No .fasta/.fa/.fna genome file uploaded.")
        
        target_fasta = DIRS["input"] / fasta_candidates[0].name
        shutil.copy2(fasta_candidates[0], target_fasta)
    else:
        if genome_path is None:
            raise ValueError("genome_path required when upload_from_colab=False")
        target_fasta = Path(genome_path)
        if not target_fasta.exists():
            raise FileNotFoundError(f"Genome file not found: {target_fasta}")
        shutil.copy2(target_fasta, DIRS["input"] / target_fasta.name)
        target_fasta = DIRS["input"] / target_fasta.name
    
    # Validate genome
    stage("VALIDATING GENOME")
    records, genome_stats = validate_genome(target_fasta)
    print(pd.Series(genome_stats).to_string())
    
    pd.DataFrame([genome_stats]).to_csv(
        RESULTS_DIR / "genome_statistics.csv", index=False
    )
    
    # Run Prokka
    stage("RUNNING PROKKA")
    prokka_dir = run_prokka(target_fasta, DIRS["prokka"])
    
    # Parse Prokka
    stage("PARSING PROKKA")
    protein_df = parse_prokka(prokka_dir)
    
    protein_df.to_csv(
        RESULTS_DIR / "master_proteins.tsv",
        sep="\t", index=False
    )
    
    prokka_annotation_map = dict(zip(
        protein_df["Protein_ID"],
        protein_df["Prokka_Annotation"]
    ))
    
    prokka_faa = prokka_dir / "PROKKA.faa"
    
    # MEROPS screening
    merops_all, merops_pass, merops_failed = screen_merops(
        prokka_faa, prokka_annotation_map
    )
    
    # Get candidate records
    prokka_records = {r.id: r for r in SeqIO.parse(prokka_faa, "fasta")}
    
    candidate_ids = set(merops_pass["Query_ID"]) if not merops_pass.empty else set()
    
    candidate_records = []
    for pid in sorted(candidate_ids):
        if pid not in prokka_records:
            continue
        rec = prokka_records[pid]
        rec.id = pid
        rec.name = pid
        rec.description = ""
        candidate_records.append(rec)
    
    CANDIDATE_FASTA = RESULTS_DIR / "MEROPS_candidates.faa"
    SeqIO.write(candidate_records, CANDIDATE_FASTA, "fasta")
    
    candidate_table = pd.DataFrame({
        "Protein_ID": [r.id for r in candidate_records],
        "Protein_Length": [len(r.seq) for r in candidate_records]
    })
    candidate_table.to_csv(RESULTS_DIR / "MEROPS_candidate_proteins.csv", index=False)
    
    # Check if we have candidates
    if not candidate_records:
        print("\nNo MEROPS PASS candidates.")
        zip_path = create_complete_zip()
        
        try:
            from google.colab import files
            files.download(str(zip_path))
        except Exception:
            print(f"Download manually from: {zip_path}")
        
        return {
            "status": "no_candidates",
            "total_proteins": len(protein_df),
            "candidate_count": 0
        }
    
    print(f"\nMEROPS candidates for InterProScan: {len(candidate_records)}")
    
    # Run InterProScan
    stage("RUNNING INTERPROSCAN")
    ipr_df = run_interproscan_rest(candidate_records)
    
    # Validate domains
    validated_df = validate_domains(ipr_df, candidate_records)
    
    # Create final audit
    audit_df = create_final_audit(validated_df)
    
    # Display results
    stage("FINAL PROTEASE CLASSIFICATION")
    counts = validated_df["Final_Classification"].value_counts()
    print(counts.to_string())
    
    # Save manifest
    manifest = save_manifest(
        protein_df, candidate_records, validated_df,
        PIPELINE_NAME, PIPELINE_VERSION
    )
    
    # Create zip
    zip_path = create_complete_zip()
    
    stage("PIPELINE COMPLETE")
    
    print(f"""
All results:
  
{RESULTS_DIR}
  
Complete ZIP:
  
{zip_path}
  
Important final files:
  
1. FINAL_protease_audit_table.csv
2. BAAP_complete_annotation.xlsx
3. TRUE_PROTEASES.faa
4. PUTATIVE_PROTEASES.faa
5. INACTIVE_PROTEASES.faa
6. NON_PROTEASES.faa
7. REVIEW_PROTEINS.faa
8. validated_proteases.csv
9. interproscan_domain_evidence.csv
10. MEROPS_all_hits.csv
11. MEROPS_PASS_hits.csv
12. MEROPS_failed_hits.csv
13. MEROPS_no_hits.csv
14. FINAL_classification_summary.csv
15. pipeline_manifest.csv
""")
    
    # Download zip
    try:
        from google.colab import files
        files.download(str(zip_path))
    except Exception:
        print(f"\nZIP location: {zip_path}")
    
    return {
        "status": "complete",
        "total_proteins": len(protein_df),
        "candidate_count": len(candidate_records),
        "true_proteases": int((validated_df["Final_Classification"] == "True protease").sum()),
        "putative_proteases": int((validated_df["Final_Classification"] == "Putative protease").sum()),
        "inactive": int((validated_df["Final_Classification"] == "Inactive").sum()),
        "non_proteases": int((validated_df["Final_Classification"] == "Non-protease").sum()),
        "review": int((validated_df["Final_Classification"] == "Review").sum())
    }


if __name__ == "__main__":
    run_pipeline()