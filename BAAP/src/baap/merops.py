"""MEROPS database and DIAMOND screening module.

BAAP — Bacterial Annotation & Analysis of Proteases

This module:

1. Downloads the MEROPS peptide-unit library.
2. Validates the downloaded FASTA.
3. Converts the library into a clean FASTA file.
4. Builds a DIAMOND protein database.
5. Screens Prokka-predicted proteins against MEROPS.
6. Calculates query and subject coverage.
7. Extracts MEROPS family and annotation information.
8. Classifies hits using configurable thresholds.
9. Produces PASS, FAILED and NO-HIT result tables.

MEROPS is used here as a sensitive protease-candidate screening
resource. Final protease classification must be supported by
independent domain/function evidence downstream.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Dict, Tuple

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


# ============================================================================
# PATHS
# ============================================================================

MEROPS_DIR = DIRS["merops"]

MEROPS_RAW = MEROPS_DIR / "pepunit.lib"
MEROPS_FASTA = MEROPS_DIR / "merops_cleaned.fasta"

MEROPS_DB = MEROPS_DIR / "merops_db"
MEROPS_DMND = Path(str(MEROPS_DB) + ".dmnd")

DIAMOND_RESULTS = MEROPS_DIR / "diamond_results.tsv"


# ============================================================================
# OUTPUT FILES
# ============================================================================

MEROPS_PASS = RESULTS_DIR / "MEROPS_PASS_hits.csv"
MEROPS_FAILED = RESULTS_DIR / "MEROPS_failed_hits.csv"
MEROPS_NOHIT = RESULTS_DIR / "MEROPS_no_hits.csv"
MEROPS_SUMMARY = RESULTS_DIR / "MEROPS_summary.csv"
MEROPS_ALL = RESULTS_DIR / "MEROPS_all_hits.csv"


# ============================================================================
# MEROPS DOWNLOAD URLS
# ============================================================================

MEROPS_URLS = [
    "https://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib",
    "https://www.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib",
]


# ============================================================================
# DIAMOND OUTPUT FORMAT
# ============================================================================

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


NUMERIC_COLUMNS = [
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


# ============================================================================
# GENERAL UTILITIES
# ============================================================================

def _check_command(command: str) -> None:
    """Ensure an external executable is available."""
    if shutil.which(command) is None:
        raise RuntimeError(
            f"Required executable '{command}' was not found in PATH."
        )


def _remove_file(path: Path) -> None:
    """Safely remove a file if it exists."""
    try:
        if path.exists():
            path.unlink()
    except OSError as exc:
        raise RuntimeError(
            f"Unable to remove file: {path}"
        ) from exc


# ============================================================================
# MEROPS DOWNLOAD
# ============================================================================

def _download_with_requests(
    url: str,
    output_path: Path,
    timeout: int = 600,
) -> bool:
    """Download MEROPS using requests."""

    try:
        log.info("Downloading MEROPS from: %s", url)

        response = requests.get(
            url,
            stream=True,
            timeout=(30, timeout),
            headers={
                "User-Agent": "BAAP/1.0 "
                              "(Bacterial Annotation & Analysis of Proteases)"
            },
        )

        response.raise_for_status()

        temporary_path = output_path.with_suffix(".download")

        with open(temporary_path, "wb") as handle:

            for chunk in response.iter_content(
                chunk_size=1024 * 1024
            ):
                if chunk:
                    handle.write(chunk)

        if (
            not temporary_path.exists()
            or temporary_path.stat().st_size == 0
        ):
            _remove_file(temporary_path)
            return False

        temporary_path.replace(output_path)

        log.info(
            "MEROPS download completed: %.2f MB",
            output_path.stat().st_size / (1024 * 1024),
        )

        return True

    except requests.RequestException as exc:

        log.warning(
            "MEROPS requests download failed: %s",
            exc,
        )

        return False

    except OSError as exc:

        log.warning(
            "MEROPS file-writing failed: %s",
            exc,
        )

        return False


def _download_with_wget(
    url: str,
    output_path: Path,
    timeout: int = 600,
) -> bool:
    """Download MEROPS using wget."""

    if shutil.which("wget") is None:
        return False

    temporary_path = output_path.with_suffix(".download")

    _remove_file(temporary_path)

    command = [
        "wget",
        "-q",
        "--show-progress",
        "--timeout",
        str(timeout),
        "--tries",
        "3",
        "-O",
        str(temporary_path),
        url,
    ]

    try:

        result = subprocess.run(
            command,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout + 60,
        )

        if (
            result.returncode == 0
            and temporary_path.exists()
            and temporary_path.stat().st_size > 0
        ):

            temporary_path.replace(output_path)

            log.info(
                "MEROPS wget download completed: %.2f MB",
                output_path.stat().st_size / (1024 * 1024),
            )

            return True

        log.warning(
            "wget failed for %s: %s",
            url,
            result.stderr.decode(
                errors="replace"
            )[-500:],
        )

    except subprocess.TimeoutExpired:

        log.warning(
            "wget timed out while downloading MEROPS."
        )

    finally:

        _remove_file(temporary_path)

    return False


def _download_with_curl(
    url: str,
    output_path: Path,
    timeout: int = 600,
) -> bool:
    """Download MEROPS using curl."""

    if shutil.which("curl") is None:
        return False

    temporary_path = output_path.with_suffix(".download")

    _remove_file(temporary_path)

    command = [
        "curl",
        "-L",
        "--fail",
        "--retry",
        "3",
        "--retry-delay",
        "3",
        "--max-time",
        str(timeout),
        "-o",
        str(temporary_path),
        url,
    ]

    try:

        result = subprocess.run(
            command,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout + 60,
        )

        if (
            result.returncode == 0
            and temporary_path.exists()
            and temporary_path.stat().st_size > 0
        ):

            temporary_path.replace(output_path)

            log.info(
                "MEROPS curl download completed: %.2f MB",
                output_path.stat().st_size / (1024 * 1024),
            )

            return True

        log.warning(
            "curl failed for %s: %s",
            url,
            result.stderr.decode(
                errors="replace"
            )[-500:],
        )

    except subprocess.TimeoutExpired:

        log.warning(
            "curl timed out while downloading MEROPS."
        )

    finally:

        _remove_file(temporary_path)

    return False


def _download_merops() -> None:
    """Download the MEROPS peptide-unit library."""

    MEROPS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    if (
        MEROPS_RAW.exists()
        and MEROPS_RAW.stat().st_size > 0
    ):
        log.info(
            "Existing MEROPS library found: %s",
            MEROPS_RAW,
        )
        return

    methods = [
        _download_with_requests,
        _download_with_wget,
        _download_with_curl,
    ]

    for method in methods:

        for url in MEROPS_URLS:

            log.info(
                "Trying %s",
                method.__name__,
            )

            if method(
                url,
                MEROPS_RAW,
            ):
                return

            time.sleep(1)

    raise RuntimeError(
        "\n"
        + "=" * 80
        + "\n"
        "MEROPS DOWNLOAD FAILED\n"
        + "=" * 80
        + "\n"
        "BAAP could not download pepunit.lib automatically.\n\n"
        f"Expected file:\n{MEROPS_RAW}\n\n"
        "Please download pepunit.lib from the MEROPS current-release "
        "directory and place it at the path above."
    )


# ============================================================================
# MEROPS VALIDATION
# ============================================================================

def _validate_merops_raw() -> None:
    """Validate that pepunit.lib contains FASTA records."""

    if (
        not MEROPS_RAW.exists()
        or MEROPS_RAW.stat().st_size == 0
    ):
        raise RuntimeError(
            f"MEROPS library is missing or empty:\n{MEROPS_RAW}"
        )

    try:

        with open(
            MEROPS_RAW,
            "r",
            encoding="latin-1",
            errors="replace",
        ) as handle:

            first_line = handle.readline().strip()

    except OSError as exc:

        raise RuntimeError(
            f"Unable to read MEROPS library:\n{MEROPS_RAW}"
        ) from exc

    if not first_line.startswith(">"):

        raise RuntimeError(
            "MEROPS pepunit.lib does not appear to be FASTA.\n"
            f"First line: {first_line[:200]}"
        )


# ============================================================================
# MEROPS FASTA PREPARATION
# ============================================================================

def _prepare_merops_fasta() -> None:
    """Convert MEROPS library into cleaned FASTA."""

    if (
        MEROPS_FASTA.exists()
        and MEROPS_FASTA.stat().st_size > 0
    ):
        log.info(
            "Existing cleaned MEROPS FASTA detected."
        )
        return

    log.info(
        "Converting MEROPS library to FASTA..."
    )

    temporary_path = (
        MEROPS_FASTA.with_suffix(".tmp")
    )

    _remove_file(temporary_path)

    sequence_count = 0

    try:

        with open(
            MEROPS_RAW,
            "r",
            encoding="latin-1",
            errors="replace",
        ) as infile, open(
            temporary_path,
            "w",
            encoding="utf-8",
        ) as outfile:

            for record in SeqIO.parse(
                infile,
                "fasta",
            ):

                sequence = str(
                    record.seq
                ).upper()

                if not sequence:
                    continue

                record.seq = record.seq.__class__(
                    sequence
                )

                SeqIO.write(
                    record,
                    outfile,
                    "fasta",
                )

                sequence_count += 1

    except Exception as exc:

        _remove_file(temporary_path)

        raise RuntimeError(
            "Failed to convert MEROPS library to FASTA."
        ) from exc

    if sequence_count == 0:

        _remove_file(temporary_path)

        raise RuntimeError(
            "MEROPS library contains no valid FASTA sequences."
        )

    temporary_path.replace(
        MEROPS_FASTA
    )

    log.info(
        "MEROPS FASTA prepared: %d sequences",
        sequence_count,
    )


# ============================================================================
# DIAMOND DATABASE
# ============================================================================

def _build_merops_database() -> None:
    """Build DIAMOND database from cleaned MEROPS FASTA."""

    _check_command("diamond")

    if (
        MEROPS_DMND.exists()
        and MEROPS_DMND.stat().st_size > 0
    ):

        log.info(
            "Existing DIAMOND MEROPS database detected."
        )

        return

    log.info(
        "Building DIAMOND MEROPS database..."
    )

    command = [
        "diamond",
        "makedb",
        "--in",
        str(MEROPS_FASTA),
        "--db",
        str(MEROPS_DB),
    ]

    try:

        result = subprocess.run(
            command,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    except OSError as exc:

        raise RuntimeError(
            "Unable to execute DIAMOND makedb."
        ) from exc

    if result.returncode != 0:

        raise RuntimeError(
            "DIAMOND database construction failed:\n"
            + result.stderr[-2000:]
        )

    if (
        not MEROPS_DMND.exists()
        or MEROPS_DMND.stat().st_size == 0
    ):

        raise RuntimeError(
            f"DIAMOND database was not created:\n{MEROPS_DMND}"
        )

    log.info(
        "DIAMOND MEROPS database ready."
    )


def download_and_prepare_merops() -> None:
    """Download, validate and prepare MEROPS for DIAMOND."""

    stage(
        "PREPARING MEROPS"
    )

    MEROPS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    _check_command("diamond")

    _download_merops()

    _validate_merops_raw()

    _prepare_merops_fasta()

    _build_merops_database()

    log.info(
        "MEROPS preparation completed successfully."
    )


# ============================================================================
# DIAMOND SCREENING
# ============================================================================

def run_diamond(
    query_faa: Path,
) -> None:
    """Run DIAMOND BLASTP against MEROPS."""

    query_faa = Path(
        query_faa
    )

    if not query_faa.exists():

        raise FileNotFoundError(
            f"Protein FASTA not found:\n{query_faa}"
        )

    if query_faa.stat().st_size == 0:

        raise ValueError(
            f"Protein FASTA is empty:\n{query_faa}"
        )

    _check_command("diamond")

    if not (
        MEROPS_DMND.exists()
        and MEROPS_DMND.stat().st_size > 0
    ):

        raise RuntimeError(
            "MEROPS DIAMOND database is unavailable."
        )

    DIAMOND_RESULTS.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    _remove_file(
        DIAMOND_RESULTS
    )

    command = [
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

    log.info(
        "Running DIAMOND BLASTP against MEROPS..."
    )

    try:

        result = subprocess.run(
            command,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    except OSError as exc:

        raise RuntimeError(
            "Unable to execute DIAMOND BLASTP."
        ) from exc

    if result.returncode != 0:

        raise RuntimeError(
            "DIAMOND BLASTP failed:\n"
            + result.stderr[-3000:]
        )

    log.info(
        "DIAMOND screening completed."
    )


# ============================================================================
# MEROPS ANNOTATION
# ============================================================================

def extract_merops_family(
    title: object,
) -> str:
    """Extract MEROPS family such as S08, M16, etc."""

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
        r"\b([A-Z]\d{1,3})(?:\.\d+)?\b",
        text,
    )

    if match:
        return match.group(1)

    return "Unknown"


def extract_merops_annotation(
    title: object,
) -> str:
    """Extract descriptive annotation from MEROPS title."""

    if pd.isna(title):
        return "NA"

    text = str(
        title
    ).strip()

    text = re.sub(
        r"\s*\[[A-Z]\d{1,3}(?:\.\d+)?\].*$",
        "",
        text,
    )

    text = re.sub(
        r"^MER\d+\s*[-:]\s*",
        "",
        text,
        flags=re.IGNORECASE,
    )

    return (
        text.strip()
        or "NA"
    )


# ============================================================================
# HIT CLASSIFICATION
# ============================================================================

def classify_merops_hit(
    row: pd.Series,
) -> str:
    """Classify one DIAMOND MEROPS hit."""

    if pd.isna(row["Evalue"]):
        return "HIGH_EVALUE"

    if row["Evalue"] > EVALUE_CUTOFF:
        return "HIGH_EVALUE"

    if pd.isna(row["Percent_Identity"]):
        return "LOW_IDENTITY"

    if row["Percent_Identity"] < IDENTITY_CUTOFF:
        return "LOW_IDENTITY"

    if pd.isna(row["Subject_Coverage"]):
        return "LOW_SUBJECT_COVERAGE"

    if row["Subject_Coverage"] < SUBJECT_COVERAGE_CUTOFF:
        return "LOW_SUBJECT_COVERAGE"

    if pd.isna(row["Query_Coverage"]):
        return "LOW_QUERY_COVERAGE"

    if row["Query_Coverage"] < QUERY_COVERAGE_CUTOFF:
        return "LOW_QUERY_COVERAGE"

    return "PASS"


# ============================================================================
# EMPTY OUTPUT HANDLING
# ============================================================================

def _write_nohit_outputs(
    protein_ids,
    annotation_map,
) -> Tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    """Write outputs when DIAMOND produces no hits."""

    rows = [
        {
            "Protein_ID": protein_id,
            "Prokka_Annotation": annotation_map.get(
                protein_id,
                "",
            ),
            "Status": "NO_MEROPS_HIT",
        }
        for protein_id in sorted(
            protein_ids
        )
    ]

    nohit_df = pd.DataFrame(
        rows
    )

    nohit_df.to_csv(
        MEROPS_NOHIT,
        index=False,
    )

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
                len(protein_ids),
                0,
                0,
                0,
                len(protein_ids),
            ],
        }
    )

    summary.to_csv(
        MEROPS_SUMMARY,
        index=False,
    )

    return (
        pd.DataFrame(),
        pd.DataFrame(),
        pd.DataFrame(),
    )


# ============================================================================
# MAIN MEROPS SCREENING
# ============================================================================

def screen_merops(
    prokka_faa,
    prokka_annotation_map: Dict[str, str],
):
    """
    Run complete MEROPS + DIAMOND screening.

    Returns
    -------
    tuple
        all_hits, passed_hits, failed_hits
    """

    stage(
        "MEROPS PROTEASE SCREENING"
    )

    RESULTS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ------------------------------------------------------------------
    # 1. Prepare MEROPS
    # ------------------------------------------------------------------

    download_and_prepare_merops()

    # ------------------------------------------------------------------
    # 2. Run DIAMOND
    # ------------------------------------------------------------------

    run_diamond(
        prokka_faa
    )

    # ------------------------------------------------------------------
    # 3. Determine query proteins
    # ------------------------------------------------------------------

    all_query_ids = set(
        prokka_annotation_map.keys()
    )

    # ------------------------------------------------------------------
    # 4. Handle no hits
    # ------------------------------------------------------------------

    if (
        not DIAMOND_RESULTS.exists()
        or DIAMOND_RESULTS.stat().st_size == 0
    ):

        log.warning(
            "DIAMOND produced no MEROPS hits."
        )

        return _write_nohit_outputs(
            all_query_ids,
            prokka_annotation_map,
        )

    # ------------------------------------------------------------------
    # 5. Read DIAMOND results
    # ------------------------------------------------------------------

    df = pd.read_csv(
        DIAMOND_RESULTS,
        sep="\t",
        names=DIAMOND_COLUMNS,
        dtype={
            "Query_ID": "string",
            "Subject_ID": "string",
            "Subject_Title": "string",
        },
    )

    if df.empty:

        return _write_nohit_outputs(
            all_query_ids,
            prokka_annotation_map,
        )

    # ------------------------------------------------------------------
    # 6. Numeric conversion
    # ------------------------------------------------------------------

    for column in NUMERIC_COLUMNS:

        df[column] = pd.to_numeric(
            df[column],
            errors="coerce",
        )

    # ------------------------------------------------------------------
    # 7. Coverage calculation
    # ------------------------------------------------------------------

    valid_query_length = (
        df["Query_Length"]
        > 0
    )

    valid_subject_length = (
        df["Subject_Length"]
        > 0
    )

    df["Query_Coverage"] = pd.NA
    df["Subject_Coverage"] = pd.NA

    df.loc[
        valid_query_length,
        "Query_Coverage",
    ] = (
        (
            df.loc[
                valid_query_length,
                "Qend",
            ]
            - df.loc[
                valid_query_length,
                "Qstart",
            ]
            + 1
        )
        / df.loc[
            valid_query_length,
            "Query_Length",
        ]
        * 100
    )

    df.loc[
        valid_subject_length,
        "Subject_Coverage",
    ] = (
        (
            df.loc[
                valid_subject_length,
                "Send",
            ]
            - df.loc[
                valid_subject_length,
                "Sstart",
            ]
            + 1
        )
        / df.loc[
            valid_subject_length,
            "Subject_Length",
        ]
        * 100
    )

    # ------------------------------------------------------------------
    # 8. MEROPS annotation
    # ------------------------------------------------------------------

    df["MEROPS_Family"] = (
        df["Subject_Title"]
        .apply(
            extract_merops_family
        )
    )

    df["MEROPS_Annotation"] = (
        df["Subject_Title"]
        .apply(
            extract_merops_annotation
        )
    )

    df["Prokka_Annotation"] = (
        df["Query_ID"]
        .map(
            prokka_annotation_map
        )
        .fillna("")
    )

    # ------------------------------------------------------------------
    # 9. Classify every hit
    # ------------------------------------------------------------------

    df["Status"] = df.apply(
        classify_merops_hit,
        axis=1,
    )

    # ------------------------------------------------------------------
    # 10. Save all hits
    # ------------------------------------------------------------------

    df.to_csv(
        MEROPS_ALL,
        index=False,
    )

    # ------------------------------------------------------------------
    # 11. Select best hit per protein
    # ------------------------------------------------------------------

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
            na_position="last",
        )
        .drop_duplicates(
            subset=["Query_ID"],
            keep="first",
        )
        .copy()
    )

    # ------------------------------------------------------------------
    # 12. PASS
    # ------------------------------------------------------------------

    passed = (
        best_hits[
            best_hits["Status"] == "PASS"
        ]
        .sort_values(
            by="Bitscore",
            ascending=False,
        )
        .copy()
    )

    # ------------------------------------------------------------------
    # 13. FAILED
    # ------------------------------------------------------------------

    failed = (
        best_hits[
            best_hits["Status"] != "PASS"
        ]
        .sort_values(
            by="Bitscore",
            ascending=False,
        )
        .copy()
    )

    # ------------------------------------------------------------------
    # 14. NO HIT
    # ------------------------------------------------------------------

    hit_ids = set(
        best_hits["Query_ID"]
    )

    no_hit_ids = (
        all_query_ids
        - hit_ids
    )

    nohit_rows = [
        {
            "Protein_ID": protein_id,
            "Prokka_Annotation": prokka_annotation_map.get(
                protein_id,
                "",
            ),
            "Status": "NO_MEROPS_HIT",
        }
        for protein_id in sorted(
            no_hit_ids
        )
    ]

    nohit_df = pd.DataFrame(
        nohit_rows
    )

    # ------------------------------------------------------------------
    # 15. Save classified results
    # ------------------------------------------------------------------

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

    # ------------------------------------------------------------------
    # 16. Summary
    # ------------------------------------------------------------------

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
