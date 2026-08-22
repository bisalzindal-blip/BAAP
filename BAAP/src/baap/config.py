"""Configuration and constants for BAAP pipeline."""

import os
import sys
import time
import shutil
import subprocess
import logging
from pathlib import Path

# Pipeline metadata
PIPELINE_NAME = "BAAP"
PIPELINE_VERSION = "1.0"

# Paths
PROJECT_DIR = "/content/protease_pipeline"
FORCE_RERUN = False
CPU_THREADS = 2
EMAIL = "newwapppp@example.com"

# MEROPS filtering parameters
EVALUE_CUTOFF = 1e-5
IDENTITY_CUTOFF = 30.0
SUBJECT_COVERAGE_CUTOFF = 60.0
QUERY_COVERAGE_CUTOFF = 20.0
MAX_TARGET_SEQS = 10

# InterProScan parameters
BATCH_SIZE = 100
SUBMISSION_SLEEP = 20
STATUS_SLEEP = 45
MAX_WAIT_MINUTES = 180

# Directory structure
ROOT = Path(PROJECT_DIR)

DIRS = {
    "input": ROOT / "input",
    "prokka": ROOT / "prokka",
    "results": ROOT / "results",
    "logs": ROOT / "logs",
    "merops": ROOT / "merops",
    "interproscan": ROOT / "interproscan",
}

for p in DIRS.values():
    p.mkdir(parents=True, exist_ok=True)

RESULTS_DIR = DIRS["results"]
LOG_FILE = DIRS["logs"] / "pipeline.log"

# Miniforge / Prokka paths
MINIFORGE_DIR = Path("/content/miniforge3")
ENV_NAME = "bisal-prokka"
ENV_DIR = MINIFORGE_DIR / "envs" / ENV_NAME
PROKKA_BIN = ENV_DIR / "bin" / "prokka"
CONDA = MINIFORGE_DIR / "bin" / "conda"
PROKKA_DB_DIR = ENV_DIR / "db"


def setup_logging():
    """Configure logging for the pipeline."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(LOG_FILE),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return logging.getLogger(PIPELINE_NAME)


log = setup_logging()


def stage(msg):
    """Print a section header."""
    print("\n" + "=" * 80)
    print(msg)
    print("=" * 80)


def run_cmd(cmd, cwd=None, stdout_path=None, stderr_path=None, check=True, env=None):
    """Run a command with logging."""
    cmd = [str(x) for x in cmd]
    log.info("RUN: %s", " ".join(cmd))
    
    stdout_h = open(stdout_path, "w") if stdout_path else subprocess.PIPE
    stderr_h = open(stderr_path, "w") if stderr_path else subprocess.PIPE
    
    full_env = os.environ.copy()
    if env:
        full_env.update(env)
    
    try:
        p = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=stdout_h,
            stderr=stderr_h,
            text=True,
            check=False,
            env=full_env,
        )
    finally:
        if stdout_path:
            stdout_h.close()
        if stderr_path:
            stderr_h.close()
    
    if check and p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {' '.join(cmd)}")
    
    return p


def install_prokka():
    """Install Prokka using conda."""
    stage("INSTALLING / VALIDATING PROKKA")
    
    if FORCE_RERUN and MINIFORGE_DIR.exists():
        shutil.rmtree(MINIFORGE_DIR)
    
    if not MINIFORGE_DIR.exists():
        subprocess.run([
            "wget", "-q", "-O", "/content/miniforge.sh",
            "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
        ], check=True)
        
        subprocess.run([
            "bash", "/content/miniforge.sh",
            "-b", "-p", str(MINIFORGE_DIR)
        ], check=True)
    
    subprocess.run([str(CONDA), "config", "--system", "--add", "channels", "conda-forge"], check=True)
    subprocess.run([str(CONDA), "config", "--system", "--add", "channels", "bioconda"], check=True)
    
    if not ENV_DIR.exists():
        r = subprocess.run(
            [str(CONDA), "create", "-y", "-n", ENV_NAME, "prokka"],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
        )
        if r.returncode != 0:
            print(r.stdout)
            raise RuntimeError("Failed to create Prokka environment.")
    
    hamap = PROKKA_DB_DIR / "hmm" / "HAMAP.hmm.h3i"
    if not hamap.exists():
        if PROKKA_DB_DIR.exists():
            shutil.rmtree(PROKKA_DB_DIR)
        PROKKA_DB_DIR.mkdir(parents=True, exist_ok=True)
        
        r = subprocess.run(
            [str(CONDA), "run", "-n", ENV_NAME, str(PROKKA_BIN), "--setupdb"],
            capture_output=True, text=True
        )
        if r.returncode != 0:
            print(r.stdout, r.stderr)
            raise RuntimeError("Failed to set up Prokka databases.")
    
    if not PROKKA_BIN.exists():
        raise RuntimeError(f"Prokka executable not found: {PROKKA_BIN}")
    
    os.environ["PROKKA_BIN"] = str(PROKKA_BIN)
    os.environ["PROKKA_DB_DIR"] = str(PROKKA_DB_DIR)
    os.environ["PATH"] = str(ENV_DIR / "bin") + os.pathsep + os.environ.get("PATH", "")
    
    print("Prokka installation ready.")