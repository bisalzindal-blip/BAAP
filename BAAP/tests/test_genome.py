"""Tests for genome module."""

import pytest
from pathlib import Path
from Bio import SeqIO
from baap.genome import validate_genome, get_genome_stats, DNA_ALPHABET


def test_dna_alphabet():
    """Test DNA alphabet contains expected characters."""
    assert "A" in DNA_ALPHABET
    assert "C" in DNA_ALPHABET
    assert "G" in DNA_ALPHABET
    assert "T" in DNA_ALPHABET
    assert "N" in DNA_ALPHABET


def test_validate_genome_valid(tmp_path):
    """Test validation of valid genome FASTA."""
    fasta_path = tmp_path / "test.fasta"
    with open(fasta_path, "w") as f:
        f.write(">contig1\nACGTACGT\n>contig2\nNNNNNNNN\n")
    
    records, stats = validate_genome(fasta_path)
    
    assert len(records) == 2
    assert stats["Genome_size_bp"] == 16
    assert stats["Number_of_contigs"] == 2
    assert stats["GC_percent"] == 50.0


def test_validate_genome_invalid(tmp_path):
    """Test validation of invalid genome FASTA."""
    fasta_path = tmp_path / "test.fasta"
    with open(fasta_path, "w") as f:
        f.write(">contig1\nACGTXACGT\n")
    
    with pytest.raises(ValueError, match="Invalid nucleotide symbols"):
        validate_genome(fasta_path)


def test_validate_genome_empty(tmp_path):
    """Test validation of empty FASTA."""
    fasta_path = tmp_path / "test.fasta"
    fasta_path.touch()
    
    with pytest.raises(ValueError, match="contains no sequences"):
        validate_genome(fasta_path)