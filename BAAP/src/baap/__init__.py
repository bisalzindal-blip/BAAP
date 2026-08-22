"""BAAP — Bacterial Annotation & Analysis of Proteases

BAAP is a bioinformatics pipeline for identifying and classifying
proteases from bacterial genomes using Prokka, MEROPS DIAMOND
screening, and InterProScan domain validation.
"""

__version__ = "1.0"
__name__ = "BAAP"

from .pipeline import run_pipeline
from .genome import validate_genome, get_genome_stats
from .prokka import run_prokka, parse_prokka
from .merops import screen_merops
from .interproscan import run_interproscan_rest
from .classification import (
    classify_merops_hit,
    evaluate_domain_evidence,
    detect_protease_class,
    validate_domains
)
from .reporting import (
    create_final_audit,
    create_complete_zip,
    save_manifest
)

__all__ = [
    "run_pipeline",
    "validate_genome",
    "get_genome_stats",
    "run_prokka",
    "parse_prokka",
    "screen_merops",
    "run_interproscan_rest",
    "classify_merops_hit",
    "evaluate_domain_evidence",
    "detect_protease_class",
    "validate_domains",
    "create_final_audit",
    "create_complete_zip",
    "save_manifest",
]
