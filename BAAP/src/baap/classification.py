"""Protease classification module with domain evidence evaluation."""

import re
import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from .config import RESULTS_DIR, log
from .interproscan import PROTEASE_ANALYSES

# Protease Pfam patterns
PROTEASE_PFAM_PATTERNS = {
    # Serine proteases
    "PF00082": "Subtilase / subtilisin family S8",
    "PF00089": "Trypsin serine protease family S1",
    "PF00326": "Prolyl oligopeptidase family S9",
    "PF00574": "ClpP serine protease family S14",
    "PF10502": "Signal peptidase S26",
    "PF01252": "Lipoprotein signal peptidase II",
    "PF01694": "Rhomboid intramembrane protease S54",
    "PF00768": "D-Ala-D-Ala carboxypeptidase S11",
    "PF02113": "D-Ala-D-Ala carboxypeptidase S13",
    
    # Metalloproteases
    "PF01546": "Peptidase M20",
    "PF07687": "Peptidase M20 dimerization domain",
    "PF05193": "Peptidase M16 C-terminal domain",
    "PF00675": "Peptidase M16 N-terminal domain",
    "PF00557": "Peptidase M24",
    "PF05343": "Peptidase M42",
    "PF02163": "Peptidase M50",
    "PF04951": "D-aminopeptidase M55",
    "PF02073": "Thermophilic metalloprotease M29",
    
    # Threonine protease
    "PF00227": "Proteasome beta subunit / T1",
    
    # Other proteolytic systems
    "PF03418": "Germination protease A25",
    "PF05908": "Poly-gamma-glutamate hydrolase",
    
    # Cautionary
    "PF01965": "DJ-1/PfpI family — protease activity not assignable from Pfam alone",
}

# True protease keywords
TRUE_PROTEASE_KEYWORDS = {
    # SERINE
    "subtilisin": "Subtilisin serine protease",
    "subtilase": "Subtilase serine protease",
    "subtilisin-like": "Subtilisin-like serine protease",
    "fervidolysin": "Fervidolysin-like subtilase",
    "thermitase": "Thermitase-like subtilase",
    "trypsin": "Trypsin-family serine protease",
    "chymotrypsin": "Chymotrypsin-family serine protease",
    "elastase": "Elastase serine protease",
    "proteinase k": "Proteinase K family",
    "prolyl oligopeptidase": "Prolyl oligopeptidase",
    "prolyl endopeptidase": "Prolyl endopeptidase",
    "dipeptidyl peptidase": "Dipeptidyl peptidase",
    "acylaminoacyl peptidase": "Acylaminoacyl peptidase",
    "clpp": "ClpP serine protease",
    "clp protease": "Clp protease",
    "lon protease": "Lon ATP-dependent protease",
    "lona": "LonA protease",
    "signal peptidase": "Signal peptidase",
    "lipoprotein signal peptidase": "Lipoprotein signal peptidase II",
    "tail-specific protease": "Tail-specific protease",
    "carboxyl-terminal protease": "C-terminal processing protease",
    "sedolisin": "Sedolisin protease",
    "rhomboid": "Rhomboid intramembrane protease",
    "rhomboid-like": "Rhomboid-like protease",
    "protease iv": "Protease IV",
    "sohb": "SohB protease",
    
    # CYSTEINE
    "papain": "Papain-family cysteine protease",
    "cathepsin": "Cathepsin cysteine protease",
    "calpain": "Calpain cysteine protease",
    "legumain": "Legumain cysteine protease",
    "caspase": "Caspase cysteine protease",
    "metacaspase": "Metacaspase",
    "gingipain": "Gingipain cysteine protease",
    "deubiquitinase": "Deubiquitinase",
    "deubiquitin": "Deubiquitinase",
    "sumo protease": "SUMO-specific protease",
    "senp": "SUMO-specific protease",
    
    # ASPARTIC
    "pepsin": "Pepsin-family aspartic protease",
    "renin": "Renin-family aspartic protease",
    "bace": "BACE aspartic protease",
    "plasmepsin": "Plasmepsin aspartic protease",
    "presenilin": "Presenilin intramembrane protease",
    "prepilin peptidase": "Prepilin peptidase",
    "type 4 prepilin peptidase": "Type IV prepilin peptidase",
    "omptin": "Omptin-family protease",
    
    # METALLO
    "aminopeptidase": "Aminopeptidase",
    "aminopeptidase n": "Aminopeptidase N",
    "alanyl aminopeptidase": "Alanyl aminopeptidase",
    "lysyl aminopeptidase": "Lysyl aminopeptidase",
    "methionyl aminopeptidase": "Methionyl aminopeptidase",
    "aminopeptidase p": "Aminopeptidase P",
    "leucyl aminopeptidase": "Leucyl aminopeptidase",
    "glutamyl aminopeptidase": "Glutamyl aminopeptidase",
    "thimet oligopeptidase": "Thimet oligopeptidase",
    "neurolysin": "Neurolysin",
    "thermolysin": "Thermolysin",
    "pseudolysin": "Pseudolysin",
    "collagenase": "Collagenase",
    "astacin": "Astacin metalloprotease",
    "meprin": "Meprin metalloprotease",
    "neprilysin": "Neprilysin",
    "carboxypeptidase": "Carboxypeptidase",
    "metallocarboxypeptidase": "Metallocarboxypeptidase",
    "pitrilysin": "Pitrilysin",
    "insulysin": "Insulysin",
    "prolidase": "Prolidase",
    "peptidase m20": "Peptidase M20",
    "peptidase m24": "Peptidase M24",
    "peptidase m28": "Peptidase M28",
    "peptidase m41": "FtsH metalloprotease",
    "ftsh": "FtsH ATP-dependent metalloprotease",
    "site-2 protease": "Site-2 metalloprotease",
    "peptidase m50": "Peptidase M50",
    "d-aminopeptidase": "D-aminopeptidase",
    "jamm": "JAMM metalloprotease",
    "jamm-like": "JAMM-like metalloprotease",
    "beta-lytic metallopeptidase": "Beta-lytic metallopeptidase",
    "staphylolysin": "Staphylolysin",
    "lysostaphin": "Lysostaphin",
    
    # THREONINE
    "proteasome": "Proteasome threonine protease",
    "hslv": "HslV threonine protease",
    "ntn hydrolase": "N-terminal nucleophile hydrolase",
    "threonine protease": "Threonine protease",
    "threonine peptidase": "Threonine peptidase",
    
    # GLUTAMIC
    "glutamic protease": "Glutamic protease",
    "eqolisin": "Eqolisin glutamic protease",
    
    # OTHER PEPTIDASES
    "peptidase": "Peptidase",
    "protease": "Protease",
    "proteinase": "Proteinase",
    "endopeptidase": "Endopeptidase",
    "exopeptidase": "Exopeptidase",
    "oligopeptidase": "Oligopeptidase",
    "d-ala-d-ala carboxypeptidase": "D-Ala-D-Ala carboxypeptidase",
    "d-alanyl-d-alanine carboxypeptidase": "D-Ala-D-Ala carboxypeptidase",
}

# Protease class keywords
CLASS_KEYWORDS = {
    "Serine": [
        "serine", "subtilisin", "subtilase", "fervidolysin",
        "thermitase", "trypsin", "chymotrypsin", "elastase",
        "clpp", "clp protease", "lon protease", "signal peptidase",
        "rhomboid", "prolyl oligopeptidase", "prolyl endopeptidase"
    ],
    "Cysteine": [
        "cysteine protease", "papain", "cathepsin", "calpain",
        "legumain", "caspase", "metacaspase", "gingipain",
        "deubiquitinase", "deubiquitin", "sumo protease"
    ],
    "Aspartic": [
        "aspartic protease", "aspartyl protease", "pepsin",
        "renin", "bace", "plasmepsin", "presenilin"
    ],
    "Metalloprotease": [
        "metalloprotease", "metalloendopeptidase", "metallopeptidase",
        "thermolysin", "pseudolysin", "collagenase", "astacin",
        "meprin", "neprilysin", "aminopeptidase", "carboxypeptidase",
        "metallocarboxypeptidase", "pitrilysin", "insulysin",
        "prolidase", "ftsh", "site-2 protease", "jamm",
        "peptidase m20", "peptidase m24", "peptidase m28", "peptidase m50"
    ],
    "Threonine": [
        "threonine protease", "threonine peptidase", "proteasome", "hslv", "ntn hydrolase"
    ],
    "Glutamic": [
        "glutamic protease", "eqolisin"
    ],
    "Peptide_Lyase": [
        "peptide lyase", "asparagine peptide lyase"
    ],
}

# Non-protease keywords
NON_PROTEASE_KEYWORDS = {
    "serpin": "Serine protease inhibitor",
    "serine protease inhibitor": "Serine protease inhibitor",
    "cystatin": "Cysteine protease inhibitor",
    "cysteine protease inhibitor": "Cysteine protease inhibitor",
    "stefin": "Cysteine protease inhibitor",
    "kunitz": "Kunitz-type protease inhibitor",
    "kazal": "Kazal-type protease inhibitor",
    "bowman-birk": "Bowman-Birk protease inhibitor",
    "alpha-2-macroglobulin": "Alpha-2-macroglobulin protease inhibitor",
    "timp": "Metalloprotease inhibitor",
    "tissue inhibitor": "Tissue inhibitor of metalloproteases",
    "trypsin inhibitor": "Trypsin inhibitor",
    "chymotrypsin inhibitor": "Chymotrypsin inhibitor",
    "elastase inhibitor": "Elastase inhibitor",
    "antithrombin": "Antithrombin protease inhibitor",
    "alpha-1-antitrypsin": "Alpha-1-antitrypsin inhibitor",
    "non-peptidase homologue": "MEROPS non-peptidase homologue",
    "non peptidase homologue": "MEROPS non-peptidase homologue",
    "inactive peptidase": "Inactive peptidase",
    "inactive protease": "Inactive protease",
    "pseudopeptidase": "Pseudopeptidase",
    "pseudoenzyme": "Pseudoenzyme",
    "catalytically inactive": "Catalytically inactive",
    "degenerate peptidase": "Degenerate peptidase",
    "degenerate active site": "Degenerate active site",
    "sortase": "Sortase / transpeptidase",
    "penicillin binding protein": "Penicillin-binding protein",
    "pbp": "Penicillin-binding protein",
    "dd-transpeptidase": "DD-transpeptidase",
    "transpeptidase": "Transpeptidase",
    "transglycosylase": "Transglycosylase",
    "lytic transglycosylase": "Lytic transglycosylase",
    "muramidase": "Muramidase",
    "lysozyme": "Lysozyme",
    "peptidoglycan": "Peptidoglycan biosynthesis",
    "murein": "Murein biosynthesis",
    "ld-transpeptidase": "LD-transpeptidase",
    "amidase": "Amidase",
    "lipase": "Lipase",
    "esterase": "Esterase",
    "phospholipase": "Phospholipase",
    "thioesterase": "Thioesterase",
    "carboxylesterase": "Carboxylesterase",
    "glycosidase": "Glycosidase",
    "amylase": "Amylase",
    "cellulase": "Cellulase",
    "glucosidase": "Glucosidase",
    "galactosidase": "Galactosidase",
    "mannosidase": "Mannosidase",
    "chitinase": "Chitinase",
    "pectinase": "Pectinase",
    "pectate lyase": "Pectate lyase",
    "glycosyl hydrolase": "Glycosyl hydrolase",
    "rnase": "RNase",
    "dnase": "DNase",
    "nuclease": "Nuclease",
    "endonuclease": "Endonuclease",
    "exonuclease": "Exonuclease",
    "helicase": "Helicase",
    "topoisomerase": "Topoisomerase",
    "gyrase": "Gyrase",
    "polymerase": "Polymerase",
    "kinase": "Kinase",
    "phosphatase": "Phosphatase",
    "phosphodiesterase": "Phosphodiesterase",
    "phosphorylase": "Phosphorylase",
    "dehydrogenase": "Dehydrogenase",
    "oxidase": "Oxidase",
    "reductase": "Reductase",
    "oxygenase": "Oxygenase",
    "hydroxylase": "Hydroxylase",
    "peroxidase": "Peroxidase",
    "catalase": "Catalase",
    "isomerase": "Isomerase",
    "mutase": "Mutase",
    "racemase": "Racemase",
    "epimerase": "Epimerase",
    "synthase": "Synthase",
    "synthetase": "Synthetase",
    "ligase": "Ligase",
    "cyclase": "Cyclase",
    "amidotransferase": "Amidotransferase",
    "aminotransferase": "Aminotransferase",
    "transaminase": "Transaminase",
    "decarboxylase": "Decarboxylase",
    "methyltransferase": "Methyltransferase",
    "acetyltransferase": "Acetyltransferase",
    "glutaminase": "Glutaminase",
    "asparaginase": "Asparaginase",
    "glutamine synthetase": "Glutamine synthetase",
    "abc transporter": "ABC transporter",
    "transporter": "Transporter",
    "porin": "Porin",
    "permease": "Permease",
    "ion channel": "Ion channel",
    "ribosomal protein": "Ribosomal protein",
    "translation factor": "Translation factor",
    "elongation factor": "Elongation factor",
    "initiation factor": "Initiation factor",
    "release factor": "Release factor",
    "transcription factor": "Transcription factor",
    "transcriptional regulator": "Transcriptional regulator",
    "sigma factor": "Sigma factor",
    "chaperone": "Chaperone",
    "chaperonin": "Chaperonin",
    "heat shock protein": "Heat shock protein",
    "flagellin": "Flagellin",
    "pilin": "Pilin",
    "fimbrial": "Fimbrial protein",
    "actin": "Actin",
    "tubulin": "Tubulin",
    "collagen": "Collagen",
    "histone": "Histone",
    "immunoglobulin": "Immunoglobulin",
    "antibody": "Antibody",
    "dna-binding": "DNA-binding protein",
    "dna repair": "DNA repair protein",
    "recombination": "Recombination protein",
    "restriction enzyme": "Restriction enzyme",
    "fha domain": "FHA signaling domain",
    "ankyrin repeat": "Ankyrin repeat",
    "wd40 repeat": "WD40 repeat",
    "tpr repeat": "TPR repeat",
    "beta-lactamase": "Beta-lactamase",
}

INACTIVE_KEYWORDS = [
    "inactive domain", "inactive peptidase", "inactive protease",
    "catalytically inactive", "catalytically dead",
    "non-peptidase homologue", "non peptidase homologue",
    "would be peptidase", "pseudoenzyme", "pseudopeptidase",
    "degenerate active site", "degenerate catalytic site",
    "inactive active site"
]

UNCERTAIN_KEYWORDS = [
    "putative", "predicted", "probable", "possible",
    "hypothetical", "uncharacterized", "family protein",
    "domain-containing protein"
]

STRONG_PROTEASE_GO_TERMS = {
    "GO:0008233": "Peptidase activity",
    "GO:0008236": "Serine-type peptidase activity",
    "GO:0004252": "Serine-type endopeptidase activity",
    "GO:0004175": "Endopeptidase activity",
    "GO:0070008": "Serine-type exopeptidase activity",
    "GO:0004180": "Carboxypeptidase activity",
    "GO:0008238": "Exopeptidase activity",
    "GO:0004222": "Metalloendopeptidase activity",
    "GO:0004197": "Cysteine-type endopeptidase activity",
}

SUPPORTIVE_GO_TERMS = {
    "GO:0008270": "Zinc ion binding",
    "GO:0005509": "Calcium ion binding",
}

_KEYWORD_CACHE = {}


def match_keyword(keyword, text):
    """Match a keyword in text using regex."""
    keyword = str(keyword).lower().strip()
    text = str(text).lower()
    
    if keyword not in _KEYWORD_CACHE:
        _KEYWORD_CACHE[keyword] = re.compile(r"(?<![a-z0-9])" + re.escape(keyword) + r"(?![a-z0-9])")
    
    return bool(_KEYWORD_CACHE[keyword].search(text))


def matched_keywords(dictionary, text):
    """Find all keywords from dictionary that match text."""
    hits = []
    for keyword in dictionary:
        if match_keyword(keyword, text):
            hits.append(keyword)
    return hits


def extract_go_ids(go_text):
    """Extract GO IDs from text."""
    if pd.isna(go_text):
        return set()
    return set(re.findall(r"GO:\d{7}", str(go_text)))


def strong_go_hits(go_text):
    """Get strong GO hits."""
    ids = extract_go_ids(go_text)
    return sorted(ids.intersection(STRONG_PROTEASE_GO_TERMS))


def supportive_go_hits(go_text):
    """Get supportive GO hits."""
    ids = extract_go_ids(go_text)
    return sorted(ids.intersection(SUPPORTIVE_GO_TERMS))


def detect_protease_class(text):
    """Detect protease class from text."""
    text = str(text).lower()
    for cls, keywords in CLASS_KEYWORDS.items():
        for keyword in keywords:
            if match_keyword(keyword, text):
                return cls
    return "Unknown"


def evaluate_domain_evidence(annotation_text, signature_accession, analysis,
                             interpro_accession, go_annotations):
    """Evaluate domain evidence for protease classification."""
    text = str(annotation_text).lower().strip()
    pfam_id = str(signature_accession).upper().strip()
    analysis = str(analysis).strip()
    interpro = str(interpro_accession).strip()
    
    strong_go = strong_go_hits(go_annotations)
    supportive_go = supportive_go_hits(go_annotations)
    
    positive_keywords = matched_keywords(TRUE_PROTEASE_KEYWORDS, text)
    negative_keywords = matched_keywords(NON_PROTEASE_KEYWORDS, text)
    
    inactive_hits = [x for x in INACTIVE_KEYWORDS if match_keyword(x, text)]
    uncertain_hits = [x for x in UNCERTAIN_KEYWORDS if match_keyword(x, text)]
    
    pfam_specific = pfam_id in PROTEASE_PFAM_PATTERNS
    pfam_description = PROTEASE_PFAM_PATTERNS.get(pfam_id, "")
    
    # DJ-1/PfpI special case
    if pfam_id == "PF01965":
        return {
            "is_protease_evidence": False,
            "strength": "REVIEW",
            "classification": "Review",
            "reason": "DJ-1/PfpI family detected; protease activity cannot be assigned from Pfam family alone.",
            "protease_class": "Review",
            "positive_keywords": "; ".join(positive_keywords),
            "negative_keywords": "; ".join(negative_keywords),
            "inactive_keywords": "; ".join(inactive_hits),
            "strong_GO": "; ".join(strong_go),
            "supportive_GO": "; ".join(supportive_go),
            "evidence_score": 0,
        }
    
    # Inactive evidence
    if inactive_hits:
        return {
            "is_protease_evidence": False,
            "strength": "INACTIVE",
            "classification": "Inactive",
            "reason": "Inactive/pseudoenzyme evidence: " + "; ".join(inactive_hits),
            "protease_class": "Inactive",
            "positive_keywords": "; ".join(positive_keywords),
            "negative_keywords": "; ".join(negative_keywords),
            "inactive_keywords": "; ".join(inactive_hits),
            "strong_GO": "; ".join(strong_go),
            "supportive_GO": "; ".join(supportive_go),
            "evidence_score": 0,
        }
    
    # Pfam-specific protease signature
    if pfam_specific:
        evidence_score = 5
        if strong_go:
            evidence_score += 2
        if positive_keywords:
            evidence_score += 1
        
        return {
            "is_protease_evidence": True,
            "strength": "STRONG",
            "classification": "True protease",
            "reason": f"Protease-associated Pfam signature {pfam_id}: {pfam_description}" +
                     ("; GO peptidase activity also present." if strong_go else ""),
            "protease_class": detect_protease_class(text + " " + pfam_description),
            "positive_keywords": "; ".join(positive_keywords),
            "negative_keywords": "; ".join(negative_keywords),
            "inactive_keywords": "; ".join(inactive_hits),
            "strong_GO": "; ".join(strong_go),
            "supportive_GO": "; ".join(supportive_go),
            "evidence_score": evidence_score,
        }
    
    # Negative keywords with GO conflict
    if negative_keywords:
        if strong_go:
            return {
                "is_protease_evidence": True,
                "strength": "CONFLICT",
                "classification": "Review",
                "reason": f"Conflicting evidence: non-protease annotation ({'; '.join(negative_keywords)}) but strong GO peptidase evidence ({'; '.join(strong_go)}).",
                "protease_class": detect_protease_class(text),
                "positive_keywords": "; ".join(positive_keywords),
                "negative_keywords": "; ".join(negative_keywords),
                "inactive_keywords": "; ".join(inactive_hits),
                "strong_GO": "; ".join(strong_go),
                "supportive_GO": "; ".join(supportive_go),
                "evidence_score": 1,
            }
        
        return {
            "is_protease_evidence": False,
            "strength": "NEGATIVE",
            "classification": "Non-protease",
            "reason": "Non-protease annotation: " + "; ".join(negative_keywords),
            "protease_class": "Non-protease",
            "positive_keywords": "; ".join(positive_keywords),
            "negative_keywords": "; ".join(negative_keywords),
            "inactive_keywords": "; ".join(inactive_hits),
            "strong_GO": "; ".join(strong_go),
            "supportive_GO": "; ".join(supportive_go),
            "evidence_score": 0,
        }
    
    # Strong GO + positive keywords
    if strong_go and positive_keywords:
        return {
            "is_protease_evidence": True,
            "strength": "STRONG",
            "classification": "True protease",
            "reason": "Protease-associated annotation supported by GO peptidase activity: " + "; ".join(strong_go),
            "protease_class": detect_protease_class(text),
            "positive_keywords": "; ".join(positive_keywords),
            "negative_keywords": "",
            "inactive_keywords": "",
            "strong_GO": "; ".join(strong_go),
            "supportive_GO": "; ".join(supportive_go),
            "evidence_score": 4,
        }
    
    # Specific keywords (not generic)
    specific_keywords = [k for k in positive_keywords if k not in {
        "protease", "proteinase", "peptidase", "endopeptidase",
        "exopeptidase", "oligopeptidase"
    }]
    
    if specific_keywords:
        return {
            "is_protease_evidence": True,
            "strength": "MODERATE",
            "classification": "Putative protease",
            "reason": f"Protease-associated annotation keyword(s): {'; '.join(specific_keywords)}; no protease-specific InterPro/Pfam signature sufficient for a True call.",
            "protease_class": detect_protease_class(text),
            "positive_keywords": "; ".join(positive_keywords),
            "negative_keywords": "",
            "inactive_keywords": "",
            "strong_GO": "",
            "supportive_GO": "; ".join(supportive_go),
            "evidence_score": 2,
        }
    
    # Generic protease wording + GO
    if strong_go and positive_keywords:
        return {
            "is_protease_evidence": True,
            "strength": "STRONG",
            "classification": "True protease",
            "reason": "Generic protease/peptidase annotation supported by GO peptidase activity.",
            "protease_class": detect_protease_class(text),
            "positive_keywords": "; ".join(positive_keywords),
            "negative_keywords": "",
            "inactive_keywords": "",
            "strong_GO": "; ".join(strong_go),
            "supportive_GO": "; ".join(supportive_go),
            "evidence_score": 4,
        }
    
    # Generic protease wording only
    if positive_keywords:
        return {
            "is_protease_evidence": True,
            "strength": "WEAK",
            "classification": "Putative protease",
            "reason": "Generic protease/peptidase wording detected without sufficiently specific domain/GO evidence.",
            "protease_class": detect_protease_class(text),
            "positive_keywords": "; ".join(positive_keywords),
            "negative_keywords": "",
            "inactive_keywords": "",
            "strong_GO": "",
            "supportive_GO": "; ".join(supportive_go),
            "evidence_score": 1,
        }
    
    # GO only
    if strong_go:
        return {
            "is_protease_evidence": True,
            "strength": "MODERATE",
            "classification": "Putative protease",
            "reason": "GO peptidase activity detected without protease-specific text/domain evidence.",
            "protease_class": detect_protease_class(text),
            "positive_keywords": "",
            "negative_keywords": "",
            "inactive_keywords": "",
            "strong_GO": "; ".join(strong_go),
            "supportive_GO": "; ".join(supportive_go),
            "evidence_score": 3,
        }
    
    # Uncertain
    if uncertain_hits:
        return {
            "is_protease_evidence": False,
            "strength": "REVIEW",
            "classification": "Review",
            "reason": "Uncertain annotation: " + "; ".join(uncertain_hits),
            "protease_class": "Review",
            "positive_keywords": "; ".join(positive_keywords),
            "negative_keywords": "",
            "inactive_keywords": "",
            "strong_GO": "",
            "supportive_GO": "; ".join(supportive_go),
            "evidence_score": 0,
        }
    
    return {
        "is_protease_evidence": False,
        "strength": "NONE",
        "classification": "Review",
        "reason": "No sufficiently informative protease or non-protease evidence.",
        "protease_class": "Unknown",
        "positive_keywords": "",
        "negative_keywords": "",
        "inactive_keywords": "",
        "strong_GO": "",
        "supportive_GO": "; ".join(supportive_go),
        "evidence_score": 0,
    }