"""MEROPS database and DIAMOND screening module."""

import re
import subprocess
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO

from .config import (
    DIRS,
    RESULTS_DIR,
    EVALUE_CUTOFF,
    IDENTITY_CUTOFF,
    SUBJECT_COVERAGE_CUTOFF,
    QUERY_COVERAGE_CUTOFF,
    MAX_TARGET_SEQS,
    CPU_THREADS,
    stage,
    log,
)



# MEROPS paths


MEROPS_DIR = DIRS["merops"]

MEROPS_RAW = MEROPS_DIR / "pepunit.lib"
MEROPS_FASTA = MEROPS_DIR / "merops_cleaned.fasta"

MEROPS_DB = MEROPS_DIR / "merops_db"
MEROPS_DMND = Path(str(MEROPS_DB) + ".dmnd")

DIAMOND_RESULTS = MEROPS_DIR / "diamond_results.tsv"

MEROPS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/merops/"
    "current_release/pepunit.lib"
)



# Output paths


MEROPS_PASS = RESULTS_DIR / "MEROPS_PASS_hits.csv"
MEROPS_FAILED = RESULTS_DIR / "MEROPS_failed_hits.csv"
MEROPS_NOHIT = RESULTS_DIR / "MEROPS_no_hits.csv"
MEROPS_SUMMARY = RESULTS_DIR / "MEROPS_summary.csv"

# DIAMOND output columns

DIAMOND_COLUMNS = [
    "Query_ID",
    "Subject_ID",
    "Subject_Title",
    "Percent_Identity",
    "Alignment_Length",
    "Mismatches",
    "Gap_Openings",
    "Qstart",
    "Qend",
    "Sstart",
    "Send",
    "Evalue",
    "Bitscore",
    "Query_Length",
    "Subject_Length",
]


# Download and prepare MEROPS


def download_and_prepare_merops():
    """Download MEROPS and prepare a DIAMOND-compatible database."""

    stage("PREPARING MEROPS")

    MEROPS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )


    # Download MEROPS database


    if not (
        MEROPS_RAW.exists()
        and MEROPS_RAW.stat().st_size > 0
    ):

        log.info("Downloading MEROPS database...")

        try:

            response = requests.get(
                MEROPS_URL,
                stream=True,
                timeout=120,
            )

            response.raise_for_status()

            with open(
                MEROPS_RAW,
                "wb",
            ) as outfile:

                for chunk in response.iter_content(
                    chunk_size=1024 * 1024
                ):

                    if chunk:
                        outfile.write(chunk)

        except requests.RequestException as exc:

            # Remove incomplete download
            if MEROPS_RAW.exists():
                MEROPS_RAW.unlink()

            raise RuntimeError(
                "MEROPS database download failed.\n"
                f"URL: {MEROPS_URL}\n"
                f"Reason: {exc}"
            ) from exc


    # Verify downloaded file

    if not (
        MEROPS_RAW.exists()
        and MEROPS_RAW.stat().st_size > 0
    ):

        raise RuntimeError(
            "MEROPS database file is missing or empty:\n"
            f"{MEROPS_RAW}"
        )

    log.info(
        "MEROPS raw database ready: %s",
        MEROPS_RAW,
    )


    # Convert MEROPS library to FASTA

    if not (
        MEROPS_FASTA.exists()
        and MEROPS_FASTA.stat().st_size > 0
    ):

        log.info(
            "Converting MEROPS library to FASTA..."
        )

        count = 0

        with open(
            MEROPS_RAW,
            "r",
            encoding="latin-1",
            errors="replace",
        ) as infile, open(
            MEROPS_FASTA,
            "w",
            encoding="utf-8",
        ) as outfile:

            for record in SeqIO.parse(
                infile,
                "fasta",
            ):

                record.seq = record.seq.upper()

                SeqIO.write(
                    record,
                    outfile,
                    "fasta",
                )

                count += 1

        if count == 0:
            raise RuntimeError(
                "MEROPS database contains no FASTA sequences."
            )

        log.info(
            "MEROPS FASTA sequences prepared: %d",
            count,
        )

 
    # Build DIAMOND database

    if not (
        MEROPS_DMND.exists()
        and MEROPS_DMND.stat().st_size > 0
    ):

        log.info(
            "Building DIAMOND MEROPS database..."
        )

        subprocess.run(
            [
                "diamond",
                "makedb",
                "--in",
                str(MEROPS_FASTA),
                "--db",
                str(MEROPS_DB),
            ],
            check=True,
        )

    if not (
        MEROPS_DMND.exists()
        and MEROPS_DMND.stat().st_size > 0
    ):

        raise RuntimeError(
            "DIAMOND database was not created:\n"
            f"{MEROPS_DMND}"
        )

    log.info("MEROPS database ready.")



# DIAMOND screening

def run_diamond(query_faa):
    """Run DIAMOND BLASTP against MEROPS."""

    DIAMOND_RESULTS.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    cmd = [
        "diamond",
        "blastp",

        "--query",
        str(query_faa),

        "--db",
        str(MEROPS_DMND),

        "--out",
        str(DIAMOND_RESULTS),

        "--outfmt",
        "6",

        "qseqid",
        "sseqid",
        "stitle",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qlen",
        "slen",

        "--evalue",
        str(EVALUE_CUTOFF),

        "--max-target-seqs",
        str(MAX_TARGET_SEQS),

        "--threads",
        str(CPU_THREADS),
    ]

    subprocess.run(
        cmd,
        check=True,
    )

# MEROPS family extraction

def extract_merops_family(title):
    """Extract MEROPS family from a MEROPS subject title."""

    if pd.isna(title):
        return "Unknown"

    text = str(title)

    match = re.search(
        r"\[([A-Z]\d{1,3})(?:\.\d+)?\]",
        text,
    )

    if match:
        return match.group(1)

    match = re.search(
        r"\[([A-Z]\d{1,3})",
        text,
    )

    return (
        match.group(1)
        if match
        else "Unknown"
    )


# MEROPS annotation extraction

def extract_merops_annotation(title):
    """Extract annotation from MEROPS subject title."""

    if pd.isna(title):
        return "NA"

    text = str(title).strip()

    text = re.split(
        r"\[[A-Z]\d{1,3}(?:\.\d+)?\]",
        text,
        maxsplit=1,
    )[0].strip()

    text = re.sub(
        r"^MER\d+\s*-\s*",
        "",
        text,
    )

    return text.strip() or "NA"



# MEROPS hit classification

def classify_merops_hit(row):
    """Classify an individual MEROPS hit."""

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


def screen_merops(
    prokka_faa,
    prokka_annotation_map,
):
    """Run the complete MEROPS screening pipeline."""

    stage("MEROPS PROTEASE SCREENING")

    # Prepare MEROPS
    download_and_prepare_merops()

    # Run DIAMOND
    run_diamond(prokka_faa)

    # Check DIAMOND output
 
    if (
        not DIAMOND_RESULTS.exists()
        or DIAMOND_RESULTS.stat().st_size == 0
    ):

        log.warning(
            "No DIAMOND hits against MEROPS."
        )

        return (
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
        )

    # Read DIAMOND results


    df = pd.read_csv(
        DIAMOND_RESULTS,
        sep="\t",
        names=DIAMOND_COLUMNS,
    )


    # Convert numeric columns

    numeric_cols = [
        "Percent_Identity",
        "Alignment_Length",
        "Mismatches",
        "Gap_Openings",
        "Qstart",
        "Qend",
        "Sstart",
        "Send",
        "Evalue",
        "Bitscore",
        "Query_Length",
        "Subject_Length",
    ]

    for column in numeric_cols:

        df[column] = pd.to_numeric(
            df[column],
            errors="coerce",
        )


    # Calculate coverage

    df["Query_Coverage"] = (
        (
            df["Qend"]
            - df["Qstart"]
            + 1
        )
        / df["Query_Length"]
    ) * 100

    df["Subject_Coverage"] = (
        (
            df["Send"]
            - df["Sstart"]
            + 1
        )
        / df["Subject_Length"]
    ) * 100

 
    # MEROPS annotations

    df["MEROPS_Family"] = (
        df["Subject_Title"]
        .apply(extract_merops_family)
    )

    df["MEROPS_Annotation"] = (
        df["Subject_Title"]
        .apply(extract_merops_annotation)
    )

    df["Prokka_Annotation"] = (
        df["Query_ID"]
        .map(prokka_annotation_map)
        .fillna("")
    )

  
    # Hit classification

    df["Status"] = df.apply(
        classify_merops_hit,
        axis=1,
    )

    # Save all hits
    df.to_csv(
        RESULTS_DIR / "MEROPS_all_hits.csv",
        index=False,
    )


    # Best hit per protein

    best_hits = (
        df.sort_values(
            by=[
                "Query_ID",
                "Bitscore",
                "Evalue",
                "Percent_Identity",
            ],
            ascending=[
                True,
                False,
                True,
                False,
            ],
        )
        .drop_duplicates(
            subset=["Query_ID"],
            keep="first",
        )
        .copy()
    )


    # PASS hits
    passed = (
        best_hits[
            best_hits["Status"] == "PASS"
        ]
        .sort_values(
            "Bitscore",
            ascending=False,
        )
        .copy()
    )

  
    # Failed hits
    failed = (
        best_hits[
            best_hits["Status"] != "PASS"
        ]
        .sort_values(
            "Bitscore",
            ascending=False,
        )
        .copy()
    )

   
    # No-hit proteins

    all_query_ids = set(
        prokka_annotation_map.keys()
    )

    hit_ids = set(
        best_hits["Query_ID"]
    )

    no_hit_ids = (
        all_query_ids
        - hit_ids
    )

    nohit_rows = []

    for protein_id in sorted(
        no_hit_ids
    ):

        nohit_rows.append(
            {
                "Protein_ID": protein_id,
                "Prokka_Annotation": (
                    prokka_annotation_map.get(
                        protein_id,
                        "",
                    )
                ),
                "Status": "NO_MEROPS_HIT",
            }
        )

    nohit_df = pd.DataFrame(
        nohit_rows
    )

    # Save outputs
    passed.to_csv(
        MEROPS_PASS,
        index=False,
    )

    failed.to_csv(
        MEROPS_FAILED,
        index=False,
    )

    nohit_df.to_csv(
        MEROPS_NOHIT,
        index=False,
    )

    # Summary
    summary = pd.DataFrame(
        {
            "Metric": [
                "Total proteins",
                "Proteins with MEROPS hit",
                "MEROPS PASS",
                "MEROPS failed",
                "No MEROPS hit",
            ],
            "Count": [
                len(all_query_ids),
                len(hit_ids),
                len(passed),
                len(failed),
                len(no_hit_ids),
            ],
        }
    )

    summary.to_csv(
        MEROPS_SUMMARY,
        index=False,
    )

    log.info(
        "MEROPS: %d PASS | %d failed | %d no-hit",
        len(passed),
        len(failed),
        len(no_hit_ids),
    )

    return (
        df,
        passed,
        failed,
    )
