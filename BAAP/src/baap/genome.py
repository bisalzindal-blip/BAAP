"""Genome validation and statistics module."""

import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO

DNA_ALPHABET = set("ACGTNRYKMSWBDHV")

log = logging.getLogger("BAAP")


def validate_genome(path):
    """Validate genome FASTA and compute statistics.
    
    Args:
        path: Path to genome FASTA file
        
    Returns:
        tuple: (records, stats_dict)
        
    Raises:
        ValueError: If invalid nucleotides are found
    """
    records = list(SeqIO.parse(path, "fasta"))
    
    if not records:
        raise ValueError("FASTA contains no sequences.")
    
    bad = []
    total = 0
    gc = 0
    
    for r in records:
        seq = str(r.seq).upper()
        invalid = sorted(set(seq) - DNA_ALPHABET)
        if invalid:
            bad.append((r.id, invalid))
        total += len(seq)
        gc += seq.count("G") + seq.count("C")
    
    if bad:
        raise ValueError(f"Invalid nucleotide symbols found: {bad[:5]}")
    
    stats = {
        "Genome_size_bp": total,
        "Number_of_contigs": len(records),
        "GC_percent": round((100 * gc / total) if total else 0, 3),
        "Input_FASTA": str(path),
    }
    
    return records, stats


def get_genome_stats(path):
    """Get genome statistics without full validation."""
    _, stats = validate_genome(path)
    return stats


def count_contigs(records):
    """Count number of contigs."""
    return len(records)


def get_total_length(records):
    """Get total genome length."""
    return sum(len(r.seq) for r in records)


def calculate_gc(records):
    """Calculate GC percentage."""
    total = 0
    gc = 0
    for r in records:
        seq = str(r.seq).upper()
        total += len(seq)
        gc += seq.count("G") + seq.count("C")
    return round((100 * gc / total) if total else 0, 3)