import pandas as pd
import numpy as np
import random
from datetime import datetime
import os
import requests
import sys

def fetch_human_genes():
    """
    Fetch comprehensive list of human protein-coding genes from Ensembl BioMart
    """
    
    # BioMart REST API endpoint
    biomart_url = "http://www.ensembl.org/biomart/martservice"
    
    # XML query for human genes with coordinates
    xml_query = """<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6">
        <Dataset name="hsapiens_gene_ensembl" interface="default">
            <Attribute name="entrezgene_id"/>
            <Attribute name="hgnc_symbol"/>
            <Attribute name="chromosome_name"/>
            <Attribute name="start_position"/>
            <Attribute name="end_position"/>
            <Attribute name="gene_biotype"/>
        </Dataset>
    </Query>"""
    
    print("Fetching human genes from Ensembl BioMart...")
    
    # Make request
    response = requests.post(biomart_url, data={'query': xml_query})
    
    if response.status_code != 200:
        print(f"Error: Failed to fetch data. Status code: {response.status_code}")
        return None
    
    # Parse TSV response
    lines = response.text.strip().split('\n')
    headers = lines[0].split('\t')
    
    genes = []
    for line in lines[1:]:
        values = line.split('\t')
        if len(values) == len(headers):
            gene_data = dict(zip(headers, values))
            
            # Filter for protein-coding genes on standard chromosomes
            if (gene_data.get('Gene type') == 'protein_coding' and 
                gene_data.get('Chromosome/scaffold name', '').replace('chr', '') in 
                ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']):
                
                # Only include genes with valid NCBI gene ID and HGNC symbol
                if gene_data.get('NCBI gene (formerly Entrezgene) ID') and gene_data.get('HGNC symbol'):
                    genes.append({
                        'gene_id': gene_data['NCBI gene (formerly Entrezgene) ID'],
                        'gene_name': gene_data['HGNC symbol'],
                        'chrom': f"chr{gene_data['Chromosome/scaffold name']}",
                        'start': int(gene_data['Gene start (bp)']),
                        'end': int(gene_data['Gene end (bp)'])
                    })
    
    print(f"Retrieved {len(genes)} protein-coding genes")
    
    # Convert to DataFrame and save
    df = pd.DataFrame(genes)
    df.to_csv('human_genes.csv', index=False)
    print(f"Saved gene list to human_genes.csv")
    
    return df

def fetch_from_hgnc():
    """
    Alternative: Fetch from HGNC (HUGO Gene Nomenclature Committee)
    """
    print("Fetching human genes from HGNC...")
    
    # HGNC provides a comprehensive gene list
    url = "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_entrez_id&col=gd_locus_type&col=gd_pub_chrom_map&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
    
    response = requests.get(url)
    
    if response.status_code != 200:
        print(f"Error fetching from HGNC: {response.status_code}")
        return None
    
    # Parse TSV
    lines = response.text.strip().split('\n')
    headers = lines[0].split('\t')
    
    genes = []
    for line in lines[1:]:
        values = line.split('\t')
        if len(values) >= 5:
            # Filter for protein-coding genes with Entrez IDs
            if 'protein-coding gene' in values[3] and values[2]:  # Has Entrez ID
                genes.append({
                    'gene_id': values[2],
                    'gene_name': values[1],
                    'hgnc_id': values[0],
                    'locus_type': values[3],
                    'cyto_location': values[4] if len(values) > 4 else ''
                })
    
    print(f"Retrieved {len(genes)} protein-coding genes from HGNC")
    
    df = pd.DataFrame(genes)
    df.to_csv('hgnc_genes.csv', index=False)
    return df


def generate_gtex_like_expression(gene_name, tissue_type):
    """
    Generate realistic expression values based on gene and tissue type
    """
    # Set seed based on gene name for reproducibility
    random.seed(hash(gene_name + tissue_type) % 2**32)
    np.random.seed(hash(gene_name + tissue_type) % 2**32)
    
    # Housekeeping genes have higher baseline expression
    housekeeping_genes = ['GAPDH', 'ACTB', 'B2M', 'HPRT1', 'RPL13A', 'YWHAZ']
    is_housekeeping = any(hk in gene_name for hk in housekeeping_genes)
    
    # Tissue-specific patterns
    tissue_patterns = {
        'brain': ['GFAP', 'MAP2', 'SYN', 'RBFOX', 'NRGN'],
        'liver': ['ALB', 'CYP', 'HNF', 'SERPINA', 'FGB'],
        'heart': ['MYH6', 'MYH7', 'TNNT2', 'NPPA', 'TTN'],
        'muscle': ['MYH1', 'MYH2', 'ACTA1', 'CKM', 'DMD'],
        'blood': ['HBB', 'HBA', 'CD4', 'CD8', 'CD19'],
        'pancreas': ['INS', 'GCG', 'AMY', 'PRSS1', 'CPA1']
    }
    
    # Base expression level
    if is_housekeeping:
        base_expr = np.random.lognormal(4, 0.5)  # Higher for housekeeping
    else:
        base_expr = np.random.lognormal(1, 1.5)  # Lower for others
    
    # Adjust for tissue specificity
    for tissue_key, genes in tissue_patterns.items():
        if tissue_key in tissue_type.lower() and any(g in gene_name for g in genes):
            base_expr *= np.random.uniform(5, 20)  # Boost for tissue-specific genes
            break
    
    # Add noise
    expression = base_expr * np.random.uniform(0.5, 1.5)
    
    # Cap at reasonable TPM values
    return min(expression, 50000)

def generate_complete_table_with_simulated_gtex():
    """
    Generate comprehensive gene table with simulated GTEx-like expression data
    """
    
    print("=" * 80)
    print("COMPREHENSIVE GENE TABLE WITH SIMULATED GTEx EXPRESSION DATA")
    print("=" * 80)
    
    # Load gene list
    print("\n1. Loading gene list...")
    genes_df = fetch_human_genes()

    total_genes = len(genes_df)
    print(f"   Loaded {total_genes:,} genes")
    
    # GTEx tissues (54 tissues from GTEx v8)
    gtex_tissues = [
        'Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Adrenal_Gland',
        'Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial',
        'Bladder', 'Brain_Amygdala', 'Brain_Anterior_cingulate_cortex_BA24',
        'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere', 'Brain_Cerebellum',
        'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus',
        'Brain_Hypothalamus', 'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia',
        'Brain_Spinal_cord_cervical_c1', 'Brain_Substantia_nigra',
        'Breast_Mammary_Tissue', 'Cells_Cultured_fibroblasts', 'Cells_EBV_transformed_lymphocytes',
        'Cervix_Ectocervix', 'Cervix_Endocervix', 'Colon_Sigmoid', 'Colon_Transverse',
        'Esophagus_Gastroesophageal_Junction', 'Esophagus_Mucosa', 'Esophagus_Muscularis',
        'Fallopian_Tube', 'Heart_Atrial_Appendage', 'Heart_Left_Ventricle',
        'Kidney_Cortex', 'Kidney_Medulla', 'Liver', 'Lung',
        'Minor_Salivary_Gland', 'Muscle_Skeletal', 'Nerve_Tibial',
        'Ovary', 'Pancreas', 'Pituitary', 'Prostate',
        'Skin_Not_Sun_Exposed_Suprapubic', 'Skin_Sun_Exposed_Lower_leg',
        'Small_Intestine_Terminal_Ileum', 'Spleen', 'Stomach',
        'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood'
    ]
    
    print(f"\n2. Generating comprehensive table with {len(gtex_tissues)} GTEx tissues...")
    
    # Initialize data dictionary
    data = {
        # Basic gene information (5 columns)
        'gene_id': genes_df['gene_id'],
        'gene_name': genes_df['gene_name'],
        'chromosome': genes_df['chrom'],
        'start_position': genes_df['start'],
        'end_position': genes_df['end'],
        
        # Additional identifiers (5 columns)
        'ensembl_id': [f"ENSG{random.randint(10000000, 99999999)}" for _ in range(total_genes)],
        'refseq_id': [f"NM_{random.randint(100000, 999999)}" if random.random() > 0.1 else '' for _ in range(total_genes)],
        'uniprot_id': [f"{random.choice(['P','Q','O'])}{random.randint(10000, 99999)}" if random.random() > 0.2 else '' for _ in range(total_genes)],
        'hgnc_id': [f"HGNC:{random.randint(1000, 50000)}" if random.random() > 0.15 else '' for _ in range(total_genes)],
        'omim_id': [str(random.randint(100000, 620000)) if random.random() > 0.7 else '' for _ in range(total_genes)],
        
        # Gene characteristics (10 columns)
        'gene_type': ['protein_coding' for _ in range(total_genes)],
        'gene_length': genes_df['end'] - genes_df['start'],
        'exon_count': [random.randint(1, 50) for _ in range(total_genes)],
        'transcript_count': [random.randint(1, 20) for _ in range(total_genes)],
        'gc_content': [random.uniform(0.3, 0.7) for _ in range(total_genes)],
        'strand': [random.choice(['+', '-']) for _ in range(total_genes)],
        'is_essential': [random.choice([True, False, None]) for _ in range(total_genes)],
        'is_housekeeping': [random.random() > 0.85 for _ in range(total_genes)],
        'paralogs_count': [random.randint(0, 10) if random.random() > 0.5 else 0 for _ in range(total_genes)],
        'orthologs_count': [random.randint(0, 500) if random.random() > 0.3 else 0 for _ in range(total_genes)]
    }
    
    # Add GTEx expression columns (54 columns)
    print("\n3. Generating GTEx expression data for all tissues...")
    for tissue in gtex_tissues:
        tissue_col = f'gtex_{tissue.lower()}'
        expression_values = []
        
        for gene_name in genes_df['gene_name']:
            if pd.notna(gene_name):
                expr_value = generate_gtex_like_expression(str(gene_name), tissue)
                expression_values.append(expr_value)
            else:
                expression_values.append(np.nan)
        
        data[tissue_col] = expression_values
        
        # Progress indicator
        if len([k for k in data.keys() if k.startswith('gtex_')]) % 10 == 0:
            print(f"   Generated expression for {len([k for k in data.keys() if k.startswith('gtex_')])} tissues...")
    
    # Additional columns
    print("\n4. Adding protein and functional annotations...")
    
    # Protein information (10 columns)
    data.update({
        'protein_length': [random.randint(50, 5000) if random.random() > 0.1 else np.nan for _ in range(total_genes)],
        'molecular_weight': [random.uniform(5000, 500000) if random.random() > 0.1 else np.nan for _ in range(total_genes)],
        'isoelectric_point': [random.uniform(4, 11) if random.random() > 0.15 else np.nan for _ in range(total_genes)],
        'instability_index': [random.uniform(10, 80) if random.random() > 0.2 else np.nan for _ in range(total_genes)],
        'aliphatic_index': [random.uniform(40, 120) if random.random() > 0.2 else np.nan for _ in range(total_genes)],
        'gravy_score': [random.uniform(-2, 2) if random.random() > 0.2 else np.nan for _ in range(total_genes)],
        'domains_count': [random.randint(0, 15) if random.random() > 0.3 else 0 for _ in range(total_genes)],
        'transmembrane_domains': [random.randint(0, 12) if random.random() > 0.7 else 0 for _ in range(total_genes)],
        'signal_peptide': [random.choice([True, False]) if random.random() > 0.2 else None for _ in range(total_genes)],
        'subcellular_location': [random.choice(['nucleus', 'cytoplasm', 'membrane', 'mitochondria', 'ER', 'golgi', 'extracellular', 'unknown']) for _ in range(total_genes)]
    })
    
    # Create DataFrame
    result_df = pd.DataFrame(data)
    
    # Calculate expression statistics
    print("\n5. Calculating expression statistics...")
    gtex_cols = [col for col in result_df.columns if col.startswith('gtex_')]
    
    result_df['mean_expression_all_tissues'] = result_df[gtex_cols].mean(axis=1)
    result_df['max_expression_all_tissues'] = result_df[gtex_cols].max(axis=1)
    result_df['min_expression_all_tissues'] = result_df[gtex_cols].min(axis=1)
    result_df['std_expression_all_tissues'] = result_df[gtex_cols].std(axis=1)
    result_df['cv_expression'] = result_df['std_expression_all_tissues'] / result_df['mean_expression_all_tissues']
    result_df['tissues_expressed_count'] = (result_df[gtex_cols] > 1).sum(axis=1)
    result_df['highly_expressed_tissues'] = (result_df[gtex_cols] > 100).sum(axis=1)
    result_df['tissue_specificity_tau'] = result_df[gtex_cols].apply(calculate_tau_score, axis=1)
    
    # Identify tissue with highest expression
    result_df['max_expression_tissue'] = result_df[gtex_cols].idxmax(axis=1).str.replace('gtex_', '')
    
    # Add metadata
    result_df['data_source'] = 'Simulated_GTEx_v8'
    result_df['last_updated'] = datetime.now()
    result_df['version'] = '1.0'
    
    # Save results
    output_file = 'comprehensive_gene_table_with_gtex_expression.csv'
    result_df.to_csv(output_file, index=False)
    
    print("\n" + "=" * 80)
    print("GENERATION COMPLETE!")
    print("=" * 80)
    
    print(f"\n📊 Table Statistics:")
    print(f"   Total rows: {len(result_df):,}")
    print(f"   Total columns: {len(result_df.columns)}")
    print(f"   GTEx tissue columns: {len(gtex_cols)}")
    print(f"   File saved: {output_file}")
    print(f"   File size: {os.path.getsize(output_file) / (1024**2):.2f} MB")
    
    return result_df

def calculate_tau_score(expression_row):
    """
    Calculate tau score for tissue specificity (Yanai et al., 2005)
    Returns value between 0 (ubiquitous) and 1 (tissue-specific)
    """
    values = expression_row.dropna()
    if len(values) == 0 or values.max() == 0:
        return np.nan
    
    # Normalize to max expression
    normalized = values / values.max()
    n = len(normalized)
    
    # Calculate tau
    tau = (n - normalized.sum()) / (n - 1) if n > 1 else 0
    
    return tau

def big_gene_table():
    _, output_file_name = sys.argv
    df = generate_complete_table_with_simulated_gtex()
    df.to_csv(output_file_name, index=False)


if __name__ == "__main__":
    df = generate_complete_table_with_simulated_gtex()