"""BAAP — Bacterial Annotation & Analysis of Proteases.

BAAP is a bioinformatics pipeline for identifying and classifying
proteases from bacterial genomes using genome annotation, MEROPS
screening, and InterProScan/domain evidence.
"""

__version__ = "1.0.0"
__name__ = "BAAP"

# Core pipeline
from .pipeline import run_pipeline

# Genome utilities
from .genome import (
    validate_genome,
    get_genome_stats,
)

# Prokka
from .prokka import (
    run_prokka,
    parse_prokka,
)

# MEROPS
from .merops import (
    screen_merops,
)

# InterProScan
from .interproscan import (
    run_interproscan_rest,
)

# Classification
from .classification import (
    evaluate_domain_evidence,
    detect_protease_class,
)

# Reporting
from .reporting import (
    create_final_audit,
    create_complete_zip,
    save_manifest,
)

__all__ = [
    # Package
    "__version__",
    "__name__",

    # Pipeline
    "run_pipeline",

    # Genome
    "validate_genome",
    "get_genome_stats",

    # Prokka
    "run_prokka",
    "parse_prokka",

    # MEROPS
    "screen_merops",

    # InterProScan
    "run_interproscan_rest",

    # Classification
    "evaluate_domain_evidence",
    "detect_protease_class",

    # Reporting
    "create_final_audit",
    "create_complete_zip",
    "save_manifest",
]
