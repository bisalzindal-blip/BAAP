"""Domain validation and protein classification module."""

import pandas as pd
from Bio import SeqIO

from .config import RESULTS_DIR
from .classification import evaluate_domain_evidence
from .interproscan import PROTEASE_ANALYSES


VALIDATED_CSV = RESULTS_DIR / "validated_proteases.csv"
VALIDATED_FASTA = RESULTS_DIR / "validated_proteases.faa"

FINAL_TRUE_FASTA = RESULTS_DIR / "TRUE_PROTEASES.faa"
FINAL_PUTATIVE_FASTA = RESULTS_DIR / "PUTATIVE_PROTEASES.faa"
FINAL_REVIEW_FASTA = RESULTS_DIR / "REVIEW_PROTEINS.faa"
FINAL_INACTIVE_FASTA = RESULTS_DIR / "INACTIVE_PROTEASES.faa"
FINAL_NONPROTEASE_FASTA = RESULTS_DIR / "NON_PROTEASES.faa"


def validate_domains(ipr_df, candidate_records):
    """Validate InterProScan domains and classify proteins."""

    stage("INTERPROSCAN DOMAIN-LEVEL VALIDATION")

    ipr_df = ipr_df.copy()

   
    required_columns = [
        "Protein_ID",
        "Analysis",
        "Signature_Accession",
    ]

    missing = [c for c in required_columns if c not in ipr_df.columns]

    if missing:
        raise ValueError(
            "InterProScan dataframe is missing required columns: "
            + ", ".join(missing)
        )

   
    for col in ["Sequence_Length", "Start", "End", "Score"]:
        if col in ipr_df.columns:
            ipr_df[col] = pd.to_numeric(
                ipr_df[col],
                errors="coerce"
            )

  
    domain_df = ipr_df[
        ipr_df["Analysis"].isin(PROTEASE_ANALYSES)
    ].copy()

    domain_df = domain_df[
        domain_df["Signature_Accession"].notna()
    ]

    domain_df = domain_df[
        domain_df["Signature_Accession"]
        .astype(str)
        .str.strip()
        .ne("")
    ]

    optional_columns = [
        "Signature_Description",
        "InterPro_Description",
        "GO_Annotations",
        "Pathways",
        "InterPro_Accession",
    ]

    for col in optional_columns:
        if col not in domain_df.columns:
            domain_df[col] = ""

    domain_df["Domain_Description"] = (
        domain_df["InterPro_Description"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    missing_description = (
        domain_df["Domain_Description"] == ""
    )

    domain_df.loc[
        missing_description,
        "Domain_Description"
    ] = (
        domain_df.loc[
            missing_description,
            "Signature_Description"
        ]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    domain_df["Pfam_ID"] = ""

    pfam_mask = (
        domain_df["Analysis"]
        .astype(str)
        .str.lower()
        .eq("pfam")
    )

    domain_df.loc[
        pfam_mask,
        "Pfam_ID"
    ] = (
        domain_df.loc[
            pfam_mask,
            "Signature_Accession"
        ]
        .fillna("")
        .astype(str)
        .str.upper()
        .str.strip()
    )

  
    duplicate_columns = [
        "Protein_ID",
        "Analysis",
        "Signature_Accession",
        "Start",
        "End",
    ]

    duplicate_columns = [
        c for c in duplicate_columns
        if c in domain_df.columns
    ]

    domain_df = (
        domain_df
        .drop_duplicates(subset=duplicate_columns)
        .reset_index(drop=True)
    )

 
    evaluated = []

    for _, row in domain_df.iterrows():

        combined_text = " ".join([
            str(row.get("Signature_Description", "")),
            str(row.get("InterPro_Description", "")),
            str(row.get("Domain_Description", "")),
        ])

        evaluation = evaluate_domain_evidence(
            annotation_text=combined_text,
            signature_accession=row.get("Pfam_ID", ""),
            analysis=row.get("Analysis", ""),
            interpro_accession=row.get(
                "InterPro_Accession",
                ""
            ),
            go_annotations=row.get(
                "GO_Annotations",
                ""
            ),
        )

        evaluated.append(evaluation)

    evaluated_df = pd.DataFrame(evaluated)

    domain_df = pd.concat(
        [
            domain_df.reset_index(drop=True),
            evaluated_df.reset_index(drop=True),
        ],
        axis=1,
    )

    domain_df.to_csv(
        RESULTS_DIR / "interproscan_domain_evidence.csv",
        index=False,
    )


    protein_rows = []

    for protein_id, group in domain_df.groupby("Protein_ID"):

        true_domains = group[
            group["classification"] == "True protease"
        ]

        putative_domains = group[
            group["classification"] == "Putative protease"
        ]

        inactive_domains = group[
            group["classification"] == "Inactive"
        ]

        nonprotease_domains = group[
            group["classification"] == "Non-protease"
        ]

        review_domains = group[
            group["classification"] == "Review"
        ]

    
        if len(true_domains) > 0:
            final_classification = "True protease"

        elif len(putative_domains) > 0:
            final_classification = "Putative protease"

        elif len(inactive_domains) > 0:
            final_classification = "Inactive"

        elif len(nonprotease_domains) > 0:
            final_classification = "Non-protease"

        else:
            final_classification = "Review"

      
        protease_domain_rows = group[
            group["is_protease_evidence"] == True
        ]

        classes = sorted(
            set(
                protease_domain_rows["protease_class"]
                .astype(str)
            )
            - {
                "Unknown",
                "Review",
                "Non-protease",
            }
        )

      
        signature_ids = sorted(
            set(
                group["Signature_Accession"]
                .dropna()
                .astype(str)
            )
        )

        interpro_ids = sorted(
            set(
                group["InterPro_Accession"]
                .dropna()
                .astype(str)
            )
            - {""}
        )

        descriptions = sorted(
            set(
                group["Domain_Description"]
                .dropna()
                .astype(str)
            )
            - {""}
        )

        reasons = sorted(
            set(
                group["reason"]
                .dropna()
                .astype(str)
            )
        )

        positive_keywords = sorted(
            set(
                kw.strip()
                for x in group["positive_keywords"].fillna("")
                for kw in str(x).split(";")
                if kw.strip()
            )
        )

        negative_keywords = sorted(
            set(
                kw.strip()
                for x in group["negative_keywords"].fillna("")
                for kw in str(x).split(";")
                if kw.strip()
            )
        )

        strong_go = sorted(
            set(
                go.strip()
                for x in group["strong_GO"].fillna("")
                for go in str(x).split(";")
                if go.strip()
            )
        )

        max_score = pd.to_numeric(
            group["evidence_score"],
            errors="coerce",
        ).max()

       
        protein_rows.append({
            "Protein_ID": protein_id,
            "Domain_Count": len(group),
            "Protease_Domain_Count": int(
                group["is_protease_evidence"].sum()
            ),
            "True_Protease_Domain_Count": len(
                true_domains
            ),
            "Putative_Protease_Domain_Count": len(
                putative_domains
            ),
            "Inactive_Domain_Count": len(
                inactive_domains
            ),
            "Non_Protease_Domain_Count": len(
                nonprotease_domains
            ),
            "Review_Domain_Count": len(
                review_domains
            ),
            "Signature_Accession": "; ".join(
                signature_ids
            ),
            "InterPro_Accession": "; ".join(
                interpro_ids
            ),
            "Domain_Description": "; ".join(
                descriptions
            ),
            "Protease_Class": "; ".join(
                classes
            ),
            "Positive_Keywords": "; ".join(
                positive_keywords
            ),
            "Negative_Keywords": "; ".join(
                negative_keywords
            ),
            "Strong_GO_Evidence": "; ".join(
                strong_go
            ),
            "Maximum_Evidence_Score": (
                int(max_score)
                if pd.notna(max_score)
                else 0
            ),
            "Classification_Reasons": (
                " | ".join(reasons)[:5000]
            ),
            "Final_Classification": final_classification,
            "Domain_Validated": (
                final_classification
                in {
                    "True protease",
                    "Putative protease",
                }
            ),
        })

    protein_classification = pd.DataFrame(
        protein_rows
    )

 
    candidate_df = pd.DataFrame({
        "Protein_ID": [
            r.id for r in candidate_records
        ],
        "Protein_Length": [
            len(r.seq) for r in candidate_records
        ],
    })

    validated_df = candidate_df.merge(
        protein_classification,
        on="Protein_ID",
        how="left",
    )

    numeric_cols = [
        "Domain_Count",
        "Protease_Domain_Count",
        "True_Protease_Domain_Count",
        "Putative_Protease_Domain_Count",
        "Inactive_Domain_Count",
        "Non_Protease_Domain_Count",
        "Review_Domain_Count",
        "Maximum_Evidence_Score",
    ]

    for col in numeric_cols:

        if col not in validated_df.columns:
            validated_df[col] = 0

        validated_df[col] = (
            pd.to_numeric(
                validated_df[col],
                errors="coerce",
            )
            .fillna(0)
            .astype(int)
        )

    if "Final_Classification" not in validated_df.columns:
        validated_df["Final_Classification"] = "Review"

    validated_df["Final_Classification"] = (
        validated_df["Final_Classification"]
        .fillna("Review")
    )

    if "Domain_Validated" not in validated_df.columns:
        validated_df["Domain_Validated"] = False

    validated_df["Domain_Validated"] = (
        validated_df["Domain_Validated"]
        .fillna(False)
        .astype(bool)
    )


    merops_pass_file = (
        RESULTS_DIR / "MEROPS_PASS_hits.csv"
    )

    if merops_pass_file.exists():

        merops_pass = pd.read_csv(
            merops_pass_file
        )

        if not merops_pass.empty:

            merops_columns = [
                "Query_ID",
                "Subject_ID",
                "Subject_Title",
                "Percent_Identity",
                "Alignment_Length",
                "Qstart",
                "Qend",
                "Sstart",
                "Send",
                "Evalue",
                "Bitscore",
                "Query_Length",
                "Subject_Length",
                "Query_Coverage",
                "Subject_Coverage",
                "MEROPS_Family",
                "MEROPS_Annotation",
            ]

            merops_columns = [
                c
                for c in merops_columns
                if c in merops_pass.columns
            ]

            merops_subset = (
                merops_pass[merops_columns]
                .rename(
                    columns={
                        "Query_ID": "Protein_ID"
                    }
                )
            )

            validated_df = validated_df.merge(
                merops_subset,
                on="Protein_ID",
                how="left",
            )

 
    validated_df.to_csv(
        VALIDATED_CSV,
        index=False,
    )


    record_map = {
        r.id: r
        for r in candidate_records
    }

    classification_files = {
        "True protease": FINAL_TRUE_FASTA,
        "Putative protease": FINAL_PUTATIVE_FASTA,
        "Review": FINAL_REVIEW_FASTA,
        "Inactive": FINAL_INACTIVE_FASTA,
        "Non-protease": FINAL_NONPROTEASE_FASTA,
    }

    for classification, filepath in classification_files.items():

        # CORRECT pandas indexing
        ids = set(
            validated_df.loc[
                validated_df["Final_Classification"]
                == classification,
                "Protein_ID",
            ]
        )

        records = [
            record_map[protein_id]
            for protein_id in ids
            if protein_id in record_map
        ]

        SeqIO.write(
            records,
            filepath,
            "fasta",
        )



    # CORRECT pandas indexing
    validated_ids = set(
        validated_df.loc[
            validated_df["Final_Classification"].isin(
                [
                    "True protease",
                    "Putative protease",
                ]
            ),
            "Protein_ID",
        ]
    )

    validated_records = [
        record_map[protein_id]
        for protein_id in validated_ids
        if protein_id in record_map
    ]

    SeqIO.write(
        validated_records,
        VALIDATED_FASTA,
        "fasta",
    )

 

    summary = (
        validated_df["Final_Classification"]
        .value_counts()
        .rename_axis("Classification")
        .reset_index(name="Count")
    )

    summary.to_csv(
        RESULTS_DIR / "FINAL_classification_summary.csv",
        index=False,
    )

    return validated_df


def stage(msg):
    """Print a section header."""

    print("\n" + "=" * 80)
    print(msg)
    print("=" * 80)
