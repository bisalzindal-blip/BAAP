"""MEROPS database and DIAMOND screening module."""

import re
import shutil
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

# Alternative URLs for MEROPS download (multiple mirrors)
MEROPS_URLS = [
    "https://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib",
    "https://www.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib",
]


def _check_command(command):
    """Check whether an external command is available."""
    if shutil.which(command) is None:
        raise RuntimeError(f"Required executable '{command}' was not found in PATH.")


def _download_with_requests(url, output_path, timeout=600):
    """Download a file using requests with progress tracking."""
    try:
        log.info("Downloading from: %s", url)
        
        session = requests.Session()
        session.mount('http://', requests.adapters.HTTPAdapter(max_retries=3))
        session.mount('https://', requests.adapters.HTTPAdapter(max_retries=3))
        
        with session.get(
            url,
            stream=True,
            timeout=(30, timeout),
            headers={
                "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
            },
            verify=False
        ) as response:
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            downloaded = 0
            
            with open(output_path, "wb") as outfile:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        outfile.write(chunk)
                        downloaded += len(chunk)
                        if total_size > 0:
                            progress = (downloaded / total_size) * 100
                            if int(progress) % 10 == 0:
                                log.info("Download progress: %.0f%%", progress)
            
            if output_path.exists() and output_path.stat().st_size > 0:
                log.info("Download completed: %.2f MB", output_path.stat().st_size / (1024 * 1024))
                return True
            else:
                log.warning("Downloaded file is empty")
                return False
                
    except requests.RequestException as exc:
        log.warning("requests download failed: %s", exc)
        if output_path.exists():
            output_path.unlink()
        return False


def _download_with_wget(url, output_path, timeout=600):
    """Download a file using wget with better error handling."""
    try:
        log.info("Attempting download with wget from: %s", url)
        
        commands = [
            [
                "wget",
                "--no-check-certificate",
                "--timeout", str(timeout),
                "--tries", "5",
                "--retry-connrefused",
                "--waitretry", "5",
                "-O", str(output_path),
                url
            ],
            [
                "wget",
                "--no-check-certificate",
                "-O", str(output_path),
                url
            ]
        ]
        
        for cmd in commands:
            try:
                result = subprocess.run(cmd, check=False, capture_output=True, timeout=timeout + 30)
                
                if result.returncode == 0:
                    if output_path.exists() and output_path.stat().st_size > 0:
                        log.info("wget download completed: %.2f MB", output_path.stat().st_size / (1024 * 1024))
                        return True
                    else:
                        log.warning("wget downloaded empty file")
                        if output_path.exists():
                            output_path.unlink()
                else:
                    log.warning("wget attempt failed with code %d", result.returncode)
                    if output_path.exists():
                        output_path.unlink()
                        
            except subprocess.TimeoutExpired:
                log.warning("wget command timed out")
                if output_path.exists():
                    output_path.unlink()
                continue
                
    except Exception as exc:
        log.warning("wget download failed: %s", exc)
        if output_path.exists():
            output_path.unlink()
        return False
    
    return False


def _download_with_curl(url, output_path, timeout=600):
    """Download a file using curl as fallback."""
    try:
        log.info("Attempting download with curl from: %s", url)
        
        subprocess.run([
            "curl",
            "-L",
            "--max-time", str(timeout),
            "--retry", "5",
            "--retry-delay", "5",
            "--insecure",
            "-o", str(output_path),
            url
        ], check=True, capture_output=True)
        
        if output_path.exists() and output_path.stat().st_size > 0:
            log.info("curl download completed: %.2f MB", output_path.stat().st_size / (1024 * 1024))
            return True
        else:
            log.warning("curl downloaded empty file")
            if output_path.exists():
                output_path.unlink()
            return False
            
    except subprocess.CalledProcessError as exc:
        log.warning("curl download failed with code %d", exc.returncode)
        if output_path.exists():
            output_path.unlink()
        return False
    except FileNotFoundError:
        log.warning("curl command not found in PATH")
        return False


def _download_merops():
    """Download the MEROPS peptide-unit library safely with multiple methods."""
    
    MEROPS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Check if file already exists and is valid
    if MEROPS_RAW.exists() and MEROPS_RAW.stat().st_size > 0:
        log.info("Existing MEROPS database found at: %s", MEROPS_RAW)
        return
    
    # Try all download methods in order of preference
    download_methods = [
        _download_with_requests,
        _download_with_wget,
        _download_with_curl
    ]
    
    for method in download_methods:
        for url in MEROPS_URLS:
            log.info("Attempting download method: %s with URL: %s", method.__name__, url)
            if method(url, MEROPS_RAW):
                return
            import time
            time.sleep(1)
    
    # If all automatic downloads fail, provide manual instructions
    raise RuntimeError(
        "\n" + "=" * 80 + "\n"
        "ERROR: Unable to download the MEROPS peptide-unit library.\n"
        "=" * 80 + "\n\n"
        "All download methods failed. Please manually download the file:\n\n"
        "1. Visit: https://ftp.ebi.ac.uk/pub/databases/merops/current_release/\n"
        "2. Download the file: pepunit.lib\n"
        "3. Upload it to the following location in your Colab environment:\n"
        f"   {MEROPS_RAW}\n\n"
        "Alternatively, try downloading it directly in a Colab cell:\n\n"
        "!wget -O /content/protease_pipeline/merops/pepunit.lib "
        "https://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib\n\n"
        "After uploading the file, run the pipeline again."
    )


def _validate_merops_fasta():
    """Check that the MEROPS raw file contains FASTA records."""
    try:
        with open(MEROPS_RAW, "r", encoding="latin-1", errors="replace") as infile:
            first_line = infile.readline().strip()
    except OSError as exc:
        raise RuntimeError(f"Unable to read MEROPS database: {MEROPS_RAW}") from exc

    if not first_line.startswith(">"):
        raise RuntimeError(
            "Downloaded MEROPS file does not appear to be a FASTA file.\n"
            f"File: {MEROPS_RAW}\n"
            f"First line: {first_line[:200]}"
        )


def download_and_prepare_merops():
    """
    Download MEROPS and prepare a DIAMOND-compatible database.
    """
    stage("PREPARING MEROPS")
    MEROPS_DIR.mkdir(parents=True, exist_ok=True)
    
    # 1. Check DIAMOND
    _check_command("diamond")
    
    # 2. Download MEROPS
    if not (MEROPS_RAW.exists() and MEROPS_RAW.stat().st_size > 0):
        log.info("MEROPS database not found. Downloading MEROPS...")
        _download_merops()
    else:
        log.info("Existing MEROPS database detected. Skipping download.")
    
    # 3. Validate MEROPS raw file
    if not (MEROPS_RAW.exists() and MEROPS_RAW.stat().st_size > 0):
        raise RuntimeError(f"MEROPS database is missing or empty:\n{MEROPS_RAW}")
    
    _validate_merops_fasta()
    
    # 4. Convert MEROPS to cleaned FASTA
    if not (MEROPS_FASTA.exists() and MEROPS_FASTA.stat().st_size > 0):
        log.info("Converting MEROPS library to FASTA...")
        
        count = 0
        temporary_fasta = MEROPS_FASTA.with_suffix(".tmp")
        
        try:
            with open(MEROPS_RAW, "r", encoding="latin-1", errors="replace") as infile, \
                 open(temporary_fasta, "w", encoding="utf-8") as outfile:
                
                for record in SeqIO.parse(infile, "fasta"):
                    sequence = str(record.seq).upper()
                    if not sequence:
                        continue
                    record.seq = sequence
                    SeqIO.write(record, outfile, "fasta")
                    count += 1
            
            if count == 0:
                raise RuntimeError("MEROPS database contains no FASTA sequences.")
            
            temporary_fasta.replace(MEROPS_FASTA)
            
        except Exception as e:
            if temporary_fasta.exists():
                temporary_fasta.unlink()
            raise RuntimeError(f"Failed to convert MEROPS to FASTA: {e}")
        
        log.info("MEROPS FASTA sequences prepared: %d", count)
    else:
        log.info("Existing cleaned MEROPS FASTA detected. Skipping conversion.")
    
    # 5. Build DIAMOND database
    if not (MEROPS_DMND.exists() and MEROPS_DMND.stat().st_size > 0):
        log.info("Building DIAMOND MEROPS database...")
        
        try:
            subprocess.run(
                [
                    "diamond",
                    "makedb",
                    "--in", str(MEROPS_FASTA),
                    "--db", str(MEROPS_DB),
                ],
                check=True,
                capture_output=True,
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"DIAMOND database build failed: {e.stderr.decode() if e.stderr else str(e)}")
    else:
        log.info("Existing DIAMOND MEROPS database detected. Skipping database construction.")
    
    # 6. Final validation
    if not (MEROPS_DMND.exists() and MEROPS_DMND.stat().st_size > 0):
        raise RuntimeError(f"DIAMOND MEROPS database was not created:\n{MEROPS_DMND}")
    
    log.info("MEROPS database ready.")


def run_diamond(query_faa):
    """Run DIAMOND BLASTP against the MEROPS database."""
    query_faa = Path(query_faa)
    
    if not query_faa.exists():
        raise FileNotFoundError(f"Protein FASTA file not found: {query_faa}")
    
    if query_faa.stat().st_size == 0:
        raise ValueError(f"Protein FASTA file is empty: {query_faa}")
    
    _check_command("diamond")
    
    DIAMOND_RESULTS.parent.mkdir(parents=True, exist_ok=True)
    
    # Remove stale results
    if DIAMOND_RESULTS.exists():
        DIAMOND_RESULTS.unlink()
    
    cmd = [
        "diamond",
        "blastp",
        "--query", str(query_faa),
        "--db", str(MEROPS_DMND),
        "--out", str(DIAMOND_RESULTS),
        "--outfmt", "6",
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
        "--evalue", str(EVALUE_CUTOFF),
        "--max-target-seqs", str(MAX_TARGET_SEQS),
        "--threads", str(CPU_THREADS),
    ]
    
    log.info("Running DIAMOND BLASTP against MEROPS...")
    
    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr.decode() if e.stderr else str(e)
        raise RuntimeError(f"DIAMOND BLASTP failed: {error_msg}")


def extract_merops_family(title):
    """Extract MEROPS family identifier from a subject title."""
    if pd.isna(title):
        return "Unknown"
    
    text = str(title)
    
    match = re.search(r"\[([A-Z]\d{1,3})(?:\.\d+)?\]", text)
    if match:
        return match.group(1)
    
    match = re.search(r"\[([A-Z]\d{1,3})", text)
    if match:
        return match.group(1)
    
    return "Unknown"


def extract_merops_annotation(title):
    """Extract annotation from a MEROPS subject title."""
    if pd.isna(title):
        return "NA"
    
    text = str(title).strip()
    
    text = re.split(r"\[[A-Z]\d{1,3}(?:\.\d+)?\]", text, maxsplit=1)[0].strip()
    text = re.sub(r"^MER\d+\s*-\s*", "", text)
    
    return text.strip() or "NA"


def classify_merops_hit(row):
    """Classify an individual MEROPS hit."""
    if pd.isna(row["Evalue"]) or row["Evalue"] > EVALUE_CUTOFF:
        return "HIGH_EVALUE"
    
    if pd.isna(row["Percent_Identity"]) or row["Percent_Identity"] < IDENTITY_CUTOFF:
        return "LOW_IDENTITY"
    
    if pd.isna(row["Subject_Coverage"]) or row["Subject_Coverage"] < SUBJECT_COVERAGE_CUTOFF:
        return "LOW_SUBJECT_COVERAGE"
    
    if pd.isna(row["Query_Coverage"]) or row["Query_Coverage"] < QUERY_COVERAGE_CUTOFF:
        return "LOW_QUERY_COVERAGE"
    
    return "PASS"


def screen_merops(prokka_faa, prokka_annotation_map):
    """
    Run the complete MEROPS screening pipeline.
    
    Returns:
        tuple: (all_hits, passed_hits, failed_hits)
    """
    stage("MEROPS PROTEASE SCREENING")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # 1. Prepare MEROPS
    download_and_prepare_merops()
    
    # 2. Run DIAMOND
    run_diamond(prokka_faa)
    
    # 3. Check DIAMOND results
    if not DIAMOND_RESULTS.exists() or DIAMOND_RESULTS.stat().st_size == 0:
        log.warning("No DIAMOND hits against MEROPS.")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    
    # 4. Read DIAMOND results
    df = pd.read_csv(DIAMOND_RESULTS, sep="\t", names=DIAMOND_COLUMNS)
    
    # 5. Convert numeric columns
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
        df[column] = pd.to_numeric(df[column], errors="coerce")
    
    # 6. Calculate alignment coverage
    df["Query_Coverage"] = ((df["Qend"] - df["Qstart"] + 1) / df["Query_Length"]) * 100
    df["Subject_Coverage"] = ((df["Send"] - df["Sstart"] + 1) / df["Subject_Length"]) * 100
    
    # Remove invalid coverage values
    df["Query_Coverage"] = df["Query_Coverage"].replace([float("inf"), float("-inf")], pd.NA)
    df["Subject_Coverage"] = df["Subject_Coverage"].replace([float("inf"), float("-inf")], pd.NA)
    
    # 7. MEROPS annotations
    df["MEROPS_Family"] = df["Subject_Title"].apply(extract_merops_family)
    df["MEROPS_Annotation"] = df["Subject_Title"].apply(extract_merops_annotation)
    df["Prokka_Annotation"] = df["Query_ID"].map(prokka_annotation_map).fillna("")
    
    # 8. Classify individual hits
    df["Status"] = df.apply(classify_merops_hit, axis=1)
    
    # 9. Save all MEROPS hits
    df.to_csv(RESULTS_DIR / "MEROPS_all_hits.csv", index=False)
    
    # 10. Select best hit for each protein
    best_hits = (
        df.sort_values(
            by=["Query_ID", "Bitscore", "Evalue", "Percent_Identity"],
            ascending=[True, False, True, False]
        )
        .drop_duplicates(subset=["Query_ID"], keep="first")
        .copy()
    )
    
    # 11. PASS hits
    passed = (
        best_hits[best_hits["Status"] == "PASS"]
        .sort_values("Bitscore", ascending=False)
        .copy()
    )
    
    # 12. Failed hits
    failed = (
        best_hits[best_hits["Status"] != "PASS"]
        .sort_values("Bitscore", ascending=False)
        .copy()
    )
    
    # 13. Proteins with no MEROPS hit
    all_query_ids = set(prokka_annotation_map.keys())
    hit_ids = set(best_hits["Query_ID"])
    no_hit_ids = all_query_ids - hit_ids
    
    nohit_rows = [
        {
            "Protein_ID": protein_id,
            "Prokka_Annotation": prokka_annotation_map.get(protein_id, ""),
            "Status": "NO_MEROPS_HIT",
        }
        for protein_id in sorted(no_hit_ids)
    ]
    
    nohit_df = pd.DataFrame(nohit_rows)
    
    # 14. Save classified outputs
    passed.to_csv(MEROPS_PASS, index=False)
    failed.to_csv(MEROPS_FAILED, index=False)
    nohit_df.to_csv(MEROPS_NOHIT, index=False)
    
    # 15. Generate summary
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
    
    summary.to_csv(MEROPS_SUMMARY, index=False)
    
    # 16. Logging
    log.info("MEROPS: %d PASS | %d failed | %d no-hit", len(passed), len(failed), len(no_hit_ids))
    
    return df, passed, failed
