# BAAP 

## Bacterial Annotation & Analysis of Proteases 
BAAP is a bioinformatics pipeline for genome-wide identification, annotation, and evidence-based classification of bacterial proteases from whole-genome nucleotide FASTA sequences.
 
### Workflow 
Genome FASTA → Genome validation → Prokka annotation → Protein prediction → MEROPS DIAMOND screening → InterProScan domain validation → Pfam/InterPro evidence → GO peptidase evidence → Evidence integration → Protease classification → Comprehensive results 

### Main classifications 
- True protease 
- Putative protease 
- Inactive 
- Non-protease 
- Review 

### Input 
BAAP accepts a bacterial whole-genome nucleotide FASTA file: 
- `.fasta` 
- `.fa` 
- `.fna` 

### Major tools and databases 
- Prokka 
- DIAMOND 
- MEROPS 
- InterProScan 
- Pfam 
- Gene Ontology 

### Output 
BAAP generates detailed CSV/TSV files, FASTA files, an Excel annotation workbook, classification summaries, an audit table, and a complete ZIP archive of the analysis results. For detailed installation and usage instructions, see the complete documentation in this repository.