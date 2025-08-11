# Gene Table Documentation

## Overview
This comprehensive gene table contains information for 19,911 human protein-coding genes with 96 columns of diverse genomic, transcriptomic, and proteomic data.

## Table Dimensions
- **Rows**: 19,911 (human protein-coding genes)
- **Columns**: 96 (diverse annotations and expression data)
- **File Size**: ~27 MB
- **Format**: CSV (comma-separated values)

## Column Descriptions

### Basic Gene Information (Columns 1-5)

#### 1. `gene_id`
- **Data Type**: Integer/String
- **Source**: NCBI Entrez Gene database
- **Description**: Unique NCBI Gene ID (formerly Entrez Gene ID)
- **Example**: "7157" (for TP53)
- **Generation**: Retrieved from Ensembl BioMart query

#### 2. `gene_name`
- **Data Type**: String
- **Source**: HUGO Gene Nomenclature Committee (HGNC)
- **Description**: Official gene symbol
- **Example**: "TP53", "BRCA1"
- **Generation**: Retrieved from Ensembl BioMart with HGNC symbols

#### 3. `chromosome`
- **Data Type**: String
- **Source**: Human genome assembly (GRCh38/hg38)
- **Description**: Chromosome location
- **Format**: "chr1", "chr2", ..., "chrX", "chrY"
- **Generation**: Retrieved from Ensembl BioMart genomic coordinates

#### 4. `start_position`
- **Data Type**: Integer
- **Source**: Human genome assembly (GRCh38/hg38)
- **Description**: Genomic start position (1-based coordinate)
- **Example**: 7661779
- **Generation**: Retrieved from Ensembl BioMart genomic coordinates

#### 5. `end_position`
- **Data Type**: Integer
- **Source**: Human genome assembly (GRCh38/hg38)
- **Description**: Genomic end position (1-based coordinate)
- **Example**: 7687538
- **Generation**: Retrieved from Ensembl BioMart genomic coordinates

### Additional Identifiers (Columns 6-10)

#### 6. `ensembl_id`
- **Data Type**: String
- **Source**: Ensembl database
- **Description**: Ensembl gene identifier
- **Format**: "ENSG" followed by 8 digits
- **Example**: "ENSG00000141510"
- **Generation**: Simulated based on Ensembl ID format

#### 7. `refseq_id`
- **Data Type**: String
- **Source**: NCBI RefSeq database
- **Description**: Reference sequence identifier for mRNA
- **Format**: "NM_" followed by 6 digits
- **Coverage**: ~90% of genes have RefSeq IDs
- **Generation**: Simulated; empty for ~10% of genes

#### 8. `uniprot_id`
- **Data Type**: String
- **Source**: UniProt database
- **Description**: UniProt protein identifier
- **Format**: Letter (P/Q/O) followed by 5 digits
- **Coverage**: ~80% of genes have UniProt IDs
- **Generation**: Simulated based on UniProt format

#### 9. `hgnc_id`
- **Data Type**: String
- **Source**: HGNC database
- **Description**: HGNC database identifier
- **Format**: "HGNC:" followed by number
- **Coverage**: ~85% of genes
- **Generation**: Simulated based on HGNC format

#### 10. `omim_id`
- **Data Type**: String
- **Source**: OMIM (Online Mendelian Inheritance in Man)
- **Description**: OMIM gene identifier
- **Format**: 6-digit number
- **Coverage**: ~30% of genes (disease-associated genes)
- **Generation**: Simulated for disease-relevant genes

### Gene Characteristics (Columns 11-20)

#### 11. `gene_type`
- **Data Type**: String
- **Values**: "protein_coding" (all entries)
- **Description**: Gene biotype classification
- **Generation**: Set to "protein_coding" for all entries (filtered dataset)

#### 12. `gene_length`
- **Data Type**: Integer
- **Description**: Total gene length in base pairs
- **Calculation**: `end_position - start_position`
- **Range**: Typically 1,000 to 2,000,000 bp

#### 13. `exon_count`
- **Data Type**: Integer
- **Description**: Number of exons in the gene
- **Range**: 1-50
- **Generation**: Simulated based on typical gene structures

#### 14. `transcript_count`
- **Data Type**: Integer
- **Description**: Number of known transcript isoforms
- **Range**: 1-20
- **Generation**: Simulated based on typical alternative splicing patterns

#### 15. `gc_content`
- **Data Type**: Float
- **Description**: GC content (proportion of G and C nucleotides)
- **Range**: 0.3-0.7 (30%-70%)
- **Generation**: Simulated based on human genome GC distribution

#### 16. `strand`
- **Data Type**: String
- **Values**: "+" (forward strand) or "-" (reverse strand)
- **Description**: DNA strand orientation
- **Generation**: Randomly assigned (approximately 50/50 distribution)

#### 17. `is_essential`
- **Data Type**: Boolean/NULL
- **Values**: True, False, or NULL (unknown)
- **Description**: Whether gene is essential for cell survival
- **Generation**: Simulated based on essential gene databases

#### 18. `is_housekeeping`
- **Data Type**: Boolean
- **Description**: Whether gene is a housekeeping gene (constitutively expressed)
- **Coverage**: ~15% marked as housekeeping
- **Generation**: Simulated (15% probability)

#### 19. `paralogs_count`
- **Data Type**: Integer
- **Description**: Number of paralogous genes (gene duplications)
- **Range**: 0-10
- **Generation**: Simulated based on gene family sizes

#### 20. `orthologs_count`
- **Data Type**: Integer
- **Description**: Number of orthologous genes across species
- **Range**: 0-500
- **Generation**: Simulated based on conservation patterns

### GTEx Expression Data (Columns 21-74)

These 54 columns contain gene expression levels from the GTEx (Genotype-Tissue Expression) project version 8, covering all major human tissues. Expression values are in TPM (Transcripts Per Million).

#### 21-74. `gtex_[tissue_name]`
- **Data Type**: Float
- **Unit**: TPM (Transcripts Per Million)
- **Source**: GTEx v8 database (simulated)
- **Description**: Median gene expression in specific tissue
- **Range**: 0 to 50,000 TPM

**Complete tissue list:**
- `gtex_adipose_subcutaneous`: Subcutaneous adipose tissue
- `gtex_adipose_visceral_omentum`: Visceral adipose tissue
- `gtex_adrenal_gland`: Adrenal gland
- `gtex_artery_aorta`: Aorta
- `gtex_artery_coronary`: Coronary artery
- `gtex_artery_tibial`: Tibial artery
- `gtex_bladder`: Bladder
- `gtex_brain_amygdala`: Brain amygdala
- `gtex_brain_anterior_cingulate_cortex_ba24`: Anterior cingulate cortex
- `gtex_brain_caudate_basal_ganglia`: Caudate (basal ganglia)
- `gtex_brain_cerebellar_hemisphere`: Cerebellar hemisphere
- `gtex_brain_cerebellum`: Cerebellum
- `gtex_brain_cortex`: Brain cortex
- `gtex_brain_frontal_cortex_ba9`: Frontal cortex (BA9)
- `gtex_brain_hippocampus`: Hippocampus
- `gtex_brain_hypothalamus`: Hypothalamus
- `gtex_brain_nucleus_accumbens_basal_ganglia`: Nucleus accumbens
- `gtex_brain_putamen_basal_ganglia`: Putamen (basal ganglia)
- `gtex_brain_spinal_cord_cervical_c1`: Spinal cord
- `gtex_brain_substantia_nigra`: Substantia nigra
- `gtex_breast_mammary_tissue`: Breast tissue
- `gtex_cells_cultured_fibroblasts`: Cultured fibroblasts
- `gtex_cells_ebv_transformed_lymphocytes`: EBV-transformed lymphocytes
- `gtex_cervix_ectocervix`: Ectocervix
- `gtex_cervix_endocervix`: Endocervix
- `gtex_colon_sigmoid`: Sigmoid colon
- `gtex_colon_transverse`: Transverse colon
- `gtex_esophagus_gastroesophageal_junction`: Gastroesophageal junction
- `gtex_esophagus_mucosa`: Esophagus mucosa
- `gtex_esophagus_muscularis`: Esophagus muscularis
- `gtex_fallopian_tube`: Fallopian tube
- `gtex_heart_atrial_appendage`: Heart atrial appendage
- `gtex_heart_left_ventricle`: Heart left ventricle
- `gtex_kidney_cortex`: Kidney cortex
- `gtex_kidney_medulla`: Kidney medulla
- `gtex_liver`: Liver
- `gtex_lung`: Lung
- `gtex_minor_salivary_gland`: Salivary gland
- `gtex_muscle_skeletal`: Skeletal muscle
- `gtex_nerve_tibial`: Tibial nerve
- `gtex_ovary`: Ovary
- `gtex_pancreas`: Pancreas
- `gtex_pituitary`: Pituitary gland
- `gtex_prostate`: Prostate
- `gtex_skin_not_sun_exposed_suprapubic`: Non-sun-exposed skin
- `gtex_skin_sun_exposed_lower_leg`: Sun-exposed skin
- `gtex_small_intestine_terminal_ileum`: Small intestine (ileum)
- `gtex_spleen`: Spleen
- `gtex_stomach`: Stomach
- `gtex_testis`: Testis
- `gtex_thyroid`: Thyroid
- `gtex_uterus`: Uterus
- `gtex_vagina`: Vagina
- `gtex_whole_blood`: Whole blood

### Protein Information (Columns 75-84)

#### 75. `protein_length`
- **Data Type**: Integer/Float
- **Unit**: Amino acids
- **Description**: Length of the protein product
- **Range**: 50-5,000 amino acids
- **Coverage**: ~90% of genes
- **Generation**: Simulated based on typical protein sizes

#### 76. `molecular_weight`
- **Data Type**: Float
- **Unit**: Daltons (Da)
- **Description**: Molecular weight of the protein
- **Range**: 5,000-500,000 Da
- **Coverage**: ~90% of genes
- **Generation**: Simulated based on protein length

#### 77. `isoelectric_point`
- **Data Type**: Float
- **Description**: pH at which protein has no net charge
- **Range**: 4.0-11.0
- **Coverage**: ~85% of genes
- **Generation**: Simulated based on typical protein pI distribution

#### 78. `instability_index`
- **Data Type**: Float
- **Description**: Protein stability score (>40 suggests instability)
- **Range**: 10-80
- **Coverage**: ~80% of genes
- **Generation**: Simulated

#### 79. `aliphatic_index`
- **Data Type**: Float
- **Description**: Relative volume of aliphatic side chains
- **Range**: 40-120
- **Coverage**: ~80% of genes
- **Generation**: Simulated based on typical values

#### 80. `gravy_score`
- **Data Type**: Float
- **Description**: Grand Average of Hydropathy
- **Range**: -2.0 to 2.0 (negative = hydrophilic, positive = hydrophobic)
- **Coverage**: ~80% of genes
- **Generation**: Simulated

#### 81. `domains_count`
- **Data Type**: Integer
- **Description**: Number of protein domains
- **Range**: 0-15
- **Generation**: Simulated based on domain architecture patterns

#### 82. `transmembrane_domains`
- **Data Type**: Integer
- **Description**: Number of transmembrane helices
- **Range**: 0-12
- **Coverage**: ~30% have at least one TM domain
- **Generation**: Simulated based on membrane protein frequency

#### 83. `signal_peptide`
- **Data Type**: Boolean/NULL
- **Description**: Presence of signal peptide for secretion
- **Values**: True, False, or NULL
- **Generation**: Simulated

#### 84. `subcellular_location`
- **Data Type**: String
- **Values**: "nucleus", "cytoplasm", "membrane", "mitochondria", "ER", "golgi", "extracellular", "unknown"
- **Description**: Primary subcellular localization
- **Generation**: Randomly assigned based on typical distributions

### Expression Statistics (Columns 85-93)

#### 85. `mean_expression_all_tissues`
- **Data Type**: Float
- **Unit**: TPM
- **Description**: Mean expression across all 54 GTEx tissues
- **Calculation**: Average of all gtex_* columns
- **Typical Range**: 0-100 TPM

#### 86. `max_expression_all_tissues`
- **Data Type**: Float
- **Unit**: TPM
- **Description**: Maximum expression in any tissue
- **Calculation**: Maximum value across all gtex_* columns

#### 87. `min_expression_all_tissues`
- **Data Type**: Float
- **Unit**: TPM
- **Description**: Minimum expression across tissues
- **Calculation**: Minimum value across all gtex_* columns

#### 88. `std_expression_all_tissues`
- **Data Type**: Float
- **Unit**: TPM
- **Description**: Standard deviation of expression across tissues
- **Calculation**: Standard deviation of all gtex_* columns

#### 89. `cv_expression`
- **Data Type**: Float
- **Description**: Coefficient of variation (std/mean)
- **Calculation**: `std_expression_all_tissues / mean_expression_all_tissues`
- **Interpretation**: Higher values indicate more variable expression

#### 90. `tissues_expressed_count`
- **Data Type**: Integer
- **Description**: Number of tissues with expression >1 TPM
- **Range**: 0-54
- **Calculation**: Count of tissues where expression > 1 TPM

#### 91. `highly_expressed_tissues`
- **Data Type**: Integer
- **Description**: Number of tissues with expression >100 TPM
- **Range**: 0-54
- **Calculation**: Count of tissues where expression > 100 TPM

#### 92. `tissue_specificity_tau`
- **Data Type**: Float
- **Description**: Tau score for tissue specificity (Yanai et al., 2005)
- **Range**: 0-1 (0 = ubiquitous, 1 = tissue-specific)
- **Calculation**: Based on expression distribution across tissues

#### 93. `max_expression_tissue`
- **Data Type**: String
- **Description**: Tissue with highest expression
- **Values**: Name of GTEx tissue (without "gtex_" prefix)
- **Calculation**: Tissue corresponding to max_expression_all_tissues

### Metadata (Columns 94-96)

#### 94. `data_source`
- **Data Type**: String
- **Value**: "Simulated_GTEx_v8"
- **Description**: Source of the data
- **Note**: Indicates data generation method

#### 95. `last_updated`
- **Data Type**: DateTime
- **Format**: YYYY-MM-DD HH:MM:SS
- **Description**: Timestamp of data generation
- **Example**: "2024-01-15 14:30:00"

#### 96. `version`
- **Data Type**: String
- **Value**: "1.0"
- **Description**: Version of the data table
- **Purpose**: Track data updates and schema changes

## Data Generation Methods

### Real Data Sources
1. **Gene coordinates**: Retrieved from Ensembl BioMart (GRCh38/hg38)
2. **Gene IDs and symbols**: From NCBI Entrez and HGNC databases

### Simulated Data
1. **Expression values**: Generated using log-normal distributions with tissue-specific patterns
2. **Protein properties**: Based on typical distributions in UniProt
3. **Functional annotations**: Simulated identifiers following standard formats

## Usage Notes

### Expression Data Interpretation
- **TPM (Transcripts Per Million)**: Normalized expression measure
  - TPM < 1: Not expressed or very low expression
  - TPM 1-10: Low expression
  - TPM 10-100: Moderate expression
  - TPM > 100: High expression
  - TPM > 1000: Very high expression

### Tissue Specificity
- **Tau score interpretation**:
  - τ < 0.3: Ubiquitously expressed
  - τ 0.3-0.7: Moderately specific
  - τ > 0.7: Tissue-specific
  - τ > 0.9: Highly tissue-specific

### Missing Data
- NaN or empty strings indicate missing/unavailable data
- Coverage varies by data type (noted in each column description)

## File Format
- **Format**: CSV (Comma-Separated Values)
- **Encoding**: UTF-8
- **Header**: First row contains column names
- **Delimiter**: Comma (,)
- **Missing values**: Empty string or NaN

## Citation
When using this data, please note that:
- Gene coordinates are from Ensembl (GRCh38/hg38)
- Expression patterns are simulated based on GTEx v8 distributions
- This is a demonstration dataset for analysis pipelines

## Updates and Versions
- **Version 1.0**: Initial release with 96 columns including 54 GTEx tissues
- **Last Updated**: See `last_updated` column for generation timestamp

## Contact
For questions about this gene table, please refer to the generation scripts:
- `generate_table_with_expression.py`: Main generation script
- `human_genes.csv`: Source gene list
- `gtex_expression_fetcher.py`: GTEx data retrieval utilities