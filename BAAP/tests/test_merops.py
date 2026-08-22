"""Tests for MEROPS module."""

import pytest
import pandas as pd
from baap.merops import (
    extract_merops_family,
    extract_merops_annotation,
    classify_merops_hit
)


def test_extract_merops_family():
    """Test MEROPS family extraction."""
    assert extract_merops_family("Protein [S01.001]") == "S01"
    assert extract_merops_family("Protein [S1]") == "S1"
    assert extract_merops_family("No family") == "Unknown"
    assert extract_merops_family(pd.NA) == "Unknown"


def test_extract_merops_annotation():
    """Test MEROPS annotation extraction."""
    assert extract_merops_annotation("Subtilisin [S01.001]") == "Subtilisin"
    assert extract_merops_annotation("MER12345 - Trypsin [S01]") == "Trypsin"
    assert extract_merops_annotation("No bracket") == "No bracket"
    assert extract_merops_annotation(pd.NA) == "NA"


def test_classify_merops_hit():
    """Test MEROPS hit classification."""
    # PASS case
    row = {
        "Evalue": 1e-10,
        "Percent_Identity": 50.0,
        "Subject_Coverage": 80.0,
        "Query_Coverage": 50.0,
    }
    assert classify_merops_hit(row) == "PASS"
    
    # HIGH_EVALUE
    row["Evalue"] = 1e-1
    assert classify_merops_hit(row) == "HIGH_EVALUE"
    
    # LOW_IDENTITY
    row["Evalue"] = 1e-10
    row["Percent_Identity"] = 20.0
    assert classify_merops_hit(row) == "LOW_IDENTITY"
    
    # LOW_SUBJECT_COVERAGE
    row["Percent_Identity"] = 50.0
    row["Subject_Coverage"] = 30.0
    assert classify_merops_hit(row) == "LOW_SUBJECT_COVERAGE"
    
    # LOW_QUERY_COVERAGE
    row["Subject_Coverage"] = 80.0
    row["Query_Coverage"] = 10.0
    assert classify_merops_hit(row) == "LOW_QUERY_COVERAGE"
    
    # NA Evalue
    row["Evalue"] = pd.NA
    assert classify_merops_hit(row) == "HIGH_EVALUE"