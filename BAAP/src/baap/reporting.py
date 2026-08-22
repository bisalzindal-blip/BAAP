"""Reporting and output generation module."""

import zipfile
import logging
from pathlib import Path

import pandas as pd

from .config import RESULTS_DIR, ROOT, EVALUE_CUTOFF, IDENTITY_CUTOFF, SUBJECT_COVERAGE_CUTOFF, QUERY_COVERAGE_CUTOFF, BATCH_SIZE

log = logging.getLogger("BAAP")


def create_final_audit(validated_df):
    """Create final audit table."""
    audit_columns = [
        "Protein_ID", "Protein_Length", "MEROPS_Family", "MEROPS_Annotation",
        "Percent_Identity", "Query_Coverage", "Subject_Coverage", "Evalue", "Bitscore",
        "Domain_Count", "Protease_Domain_Count", "True_Protease_Domain_Count",
        "Putative_Protease_Domain_Count", "Inactive_Domain_Count", "Non_Protease_Domain_Count",
        "Review_Domain_Count", "Signature_Accession", "InterPro_Accession",
        "Domain_Description", "Protease_Class", "Positive_Keywords", "Negative_Keywords",
        "Strong_GO_Evidence", "Maximum_Evidence_Score", "Classification_Reasons",
        "Final_Classification", "Domain_Validated"
    ]
    
    audit_columns = [c for c in audit_columns if c in validated_df.columns]
    audit_df = validated_df[audit_columns].copy()
    
    audit_df.to_csv(RESULTS_DIR / "FINAL_protease_audit_table.csv", index=False)
    
    # Excel workbook
    excel_path = RESULTS_DIR / "BAAP_complete_annotation.xlsx"
    
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        audit_df.to_excel(writer, sheet_name="Final_Audit", index=False)
        
        for sheet_name, file_name in [
            ("InterPro_Domains", "interproscan_domain_evidence.csv"),
            ("MEROPS_All_Hits", "MEROPS_all_hits.csv"),
            ("MEROPS_PASS", "MEROPS_PASS_hits.csv"),
            ("MEROPS_Failed", "MEROPS_failed_hits.csv"),
            ("Summary", "FINAL_classification_summary.csv"),
        ]:
            file_path = RESULTS_DIR / file_name
            if file_path.exists():
                pd.read_csv(file_path).to_excel(writer, sheet_name=sheet_name, index=False)
    
    return audit_df


def create_complete_zip():
    """Create a zip archive of all results."""
    stage("CREATING COMPLETE RESULTS ZIP")
    
    zip_path = ROOT / "BAAP_complete_results.zip"
    
    if zip_path.exists():
        zip_path.unlink()
    
    excluded = {zip_path.resolve()}
    
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for path in ROOT.rglob("*"):
            if not path.is_file():
                continue
            if path.resolve() in excluded:
                continue
            try:
                arcname = path.relative_to(ROOT)
                z.write(path, arcname=str(arcname))
            except Exception as e:
                log.warning("Could not add %s: %s", path, e)
    
    print(f"\nCOMPLETE ZIP:\n{zip_path}")
    return zip_path


def save_manifest(protein_df, candidate_records, validated_df, pipeline_name, pipeline_version):
    """Save pipeline manifest."""
    manifest = {
        "Pipeline": pipeline_name,
        "Version": pipeline_version,
        "Workflow": "Prokka → MEROPS DIAMOND → InterProScan → evidence integration",
        "MEROPS_Evalue": EVALUE_CUTOFF,
        "MEROPS_identity_percent": IDENTITY_CUTOFF,
        "MEROPS_subject_coverage_percent": SUBJECT_COVERAGE_CUTOFF,
        "MEROPS_query_coverage_percent": QUERY_COVERAGE_CUTOFF,
        "InterProScan_batch_size": BATCH_SIZE,
        "Total_Prokka_proteins": len(protein_df),
        "MEROPS_candidates": len(candidate_records),
        "True_proteases": int((validated_df["Final_Classification"] == "True protease").sum()),
        "Putative_proteases": int((validated_df["Final_Classification"] == "Putative protease").sum()),
        "Inactive": int((validated_df["Final_Classification"] == "Inactive").sum()),
        "Non_protease": int((validated_df["Final_Classification"] == "Non-protease").sum()),
        "Review": int((validated_df["Final_Classification"] == "Review").sum()),
    }
    
    pd.DataFrame([manifest]).T.rename(columns={0: "Value"}).to_csv(
        RESULTS_DIR / "pipeline_manifest.csv"
    )
    
    return manifest


def stage(msg):
    """Print a section header."""
    print("\n" + "=" * 80)
    print(msg)
    print("=" * 80)