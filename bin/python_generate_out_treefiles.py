#!/usr/bin/env python

import os
import sys
import json
import glob
import logging
import itertools
import numpy as np
import pandas as pd
from pybedtools import BedTool


## for nextflow
script_path = os.path.realpath(__file__)
bin_dir = os.path.dirname(script_path)
project_root = os.path.dirname(bin_dir)
utils_dir = os.path.join(project_root, 'utils')
sys.path.append(utils_dir)

GENE_INFO_FILE = os.path.join(utils_dir, 'gencode.v45.chr_patch_hapl_scaff.basic.annotation.bed')


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def process_tree_summary(treesum_file, out_path):
    """Parses the tree summary JSON and exports a TSV."""
    with open(treesum_file) as f:
        data = json.load(f)['trees']
        
    rows = []
    for model_idx, model_data in data.items():
        pops = pd.DataFrame.from_dict(model_data['populations'], orient='index')
        pops['model_index'] = int(model_idx) + 1
        
        # Expand cellular prevalence
        cp_df = pd.DataFrame(pops['cellular_prevalence'].tolist(), index=pops.index)
        cp_df.columns = [f'cp_{i+1}' for i in range(cp_df.shape[1])]
        
        pops = pd.concat([pops.drop(columns=['cellular_prevalence']), cp_df], axis=1)
        
        # Add tree metrics
        for metric in ['clustering_index', 'branching_index', 'linearity_index', 'llh']:
            pops[metric] = model_data[metric]
            
        pops = pops.reset_index().rename(columns={'index': 'population'})
        
        # Handle tree structure
        structure_row = [(np.nan, '0')]
        for parent, children in model_data['structure'].items():
            for child in children:
                structure_row.append((parent, str(child)))
        
        struct_df = pd.DataFrame(structure_row, columns=['parent', 'child'])
        pops = pd.merge(pops, struct_df, left_on='population', right_on='child', how='left')
        rows.append(pops)
        
    df_tree_sum = pd.concat(rows, ignore_index=True)
    df_tree_sum = df_tree_sum.sort_values(by=['llh', 'population'], ascending=[False, True]).reset_index(drop=True)
    df_tree_sum.to_csv(out_path, sep='\t', index=False)
    return df_tree_sum

def parse_ssm_cnv_files(ssm_file, cnvs_file):
    """Cleans and prepares SSM and CNV dataframes."""
    # Process SSMs
    ssm_df = pd.read_csv(ssm_file, sep='\t')[['id', 'gene']].rename(columns={'id':'ssms'})
    coords = ssm_df['gene'].str.split('_', expand=True)
    ssm_df['chrom'] = 'chr' + coords[0]
    ssm_df['end'] = coords[1].astype(int)
    ssm_df['start'] = ssm_df['end'] - 1
    ssm_df = ssm_df.drop('gene', axis=1)

    # Process CNVs
    cnvs_raw = pd.read_csv(cnvs_file, sep='\t')[['cnv', 'physical_cnvs']]
    cnv_list = []
    for _, row in cnvs_raw.iterrows():
        for item in row['physical_cnvs'].split(';'):
            parts = {p.split('=')[0]: p.split('=')[1] for p in item.split(',')}
            cnv_list.append({
                'cnvs': row['cnv'],
                'chrom': 'chr' + parts['chrom'],
                'start': int(parts['start']),
                'end': int(parts['end'])
            })
    return ssm_df, pd.DataFrame(cnv_list)

def get_mutation_assignments(mutass_dir, model_index):
    """Loads mutation assignments from the top-scoring model."""
    path = f"{mutass_dir}/{model_index}.json"
    with open(path) as f:
        data = json.load(f)['mut_assignments']
    
    dfs = []
    for pop_idx, content in data.items():
        # Ensure equal length for DataFrame creation
        max_l = max(len(v) for v in content.values())
        content_padded = {k: v + [np.nan]*(max_l - len(v)) for k, v in content.items()}
        
        df = pd.DataFrame(content_padded)
        df['population'] = pop_idx
        dfs.append(df[['population', 'ssms', 'cnvs']])
    
    return pd.concat(dfs)

def run_genomic_intersection(mutass_df, ssm_df, cnvs_df, df_anno, df_gene_of_interest, out_path):
    """
    Handles BedTool logic: intersects CNVs and SSMs, applies gene annotations,
    and filters for genes of interest.
    """
    # Prepare CP DataFrame
    cp_df = mutass_df.copy()
    cp_df = cp_df[['population'] + [col for col in cp_df.columns if 'cp_' in col]]
    cp_df['population'] = cp_df['population'].astype(int)
    cp_df = cp_df.drop_duplicates()
    
    # Prepare Bed-ready DataFrames
    m_ssms = pd.merge(mutass_df[['population', 'ssms']], ssm_df, on='ssms').drop_duplicates()
    m_cnvs = pd.merge(mutass_df[['population', 'cnvs']], cnvs_df, on='cnvs').drop_duplicates()
    
    # Define columns for BedTool conversion
    ssm_cols = ['chrom', 'start', 'end', 'population', 'ssms']
    cnv_cols = ['chrom', 'start', 'end', 'population', 'cnvs']
    
    # convert into bed format
    bt_ssms = BedTool.from_dataframe(m_ssms[ssm_cols])
    bt_cnvs = BedTool.from_dataframe(m_cnvs[cnv_cols])
    
    # Find Overlaps (CNVs and SSMs in same regions)
    bt_intersect = bt_cnvs.intersect(bt_ssms, wa=True, wb=True)
    
    # Reconstruct headers for the intersection result
    intersect_cols = [f"{c}_cnvs" for c in cnv_cols] + [f"{c}_ssms" for c in ssm_cols]
    df_intersect = bt_intersect.to_dataframe(names=intersect_cols)
    
    # Process Overlapping Sets
    bt_overlap_cnvs = df_intersect[['chrom_cnvs','start_cnvs','end_cnvs','population_cnvs','cnvs_cnvs']] \
        .drop_duplicates(subset=['population_cnvs', 'start_cnvs']) \
        .rename(columns={'chrom_cnvs':'chrom','start_cnvs':'start', 'end_cnvs':'end', 
                         'population_cnvs':'population', 'cnvs_cnvs':'ssms_cnvs'})

    bt_overlap_ssms = df_intersect[['chrom_ssms','start_ssms','end_ssms','population_ssms','ssms_ssms']] \
        .drop_duplicates(subset=['population_ssms', 'start_ssms']) \
        .rename(columns={'chrom_ssms':'chrom','start_ssms':'start', 'end_ssms':'end', 
                         'population_ssms':'population', 'ssms_ssms':'ssms_cnvs'})
    
    # Prepare SSM and CNV bed-ready data
    m_ssms = pd.merge(mutass_df[['population', 'ssms']], ssm_df, on='ssms').drop_duplicates()
    m_cnvs = pd.merge(mutass_df[['population', 'cnvs']], cnvs_df, on='cnvs').drop_duplicates()
    
    bt_ssms = BedTool.from_dataframe(m_ssms[['chrom', 'start', 'end', 'population', 'ssms']])
    bt_cnvs = BedTool.from_dataframe(m_cnvs[['chrom', 'start', 'end', 'population', 'cnvs']])
    
    # Handle Non-overlapping Mutations
    bt_cnvs_only = bt_cnvs.intersect(bt_ssms, v=True).to_dataframe(names=['chrom', 'start', 'end', 'population', 'ssms_cnvs'])
    bt_ssms_only = bt_ssms.intersect(bt_cnvs, v=True).to_dataframe(names=['chrom', 'start', 'end', 'population', 'ssms_cnvs'])

    # Combine all mutation types
    combined_muts = pd.concat([bt_overlap_cnvs, bt_overlap_ssms, bt_cnvs_only, bt_ssms_only])
    combined_muts['gene_coordinate'] = (combined_muts['chrom'].astype(str) + '_' + 
                                        combined_muts['start'].astype(str) + '_' + 
                                        combined_muts['end'].astype(str))
    
    # Annotation Logic
    # Merge with priority annotation first
    annotated = pd.merge(combined_muts, df_anno, on='gene_coordinate', how='left')
    annotated = annotated.loc[~((annotated['vaf'].isna()) & (annotated['ssms_cnvs'].str.startswith('s')))]
    
    # For rows missing priority annotation, use the BedTool gene file
    missing_anno = annotated[annotated['gene_names'].isnull()].copy()
    
    if not missing_anno.empty:
        # Load global gene info
        df_gene_info = pd.read_csv(GENE_INFO_FILE, sep='\t')
        gene_info_cols = df_gene_info.columns.tolist()
        bt_gene_info = BedTool.from_dataframe(df_gene_info)
        
        bt_missing = BedTool.from_dataframe(missing_anno)
        # wa=True, wb=True here adds the gene info to the mutation record
        bt_extra_anno = bt_gene_info.intersect(bt_missing, wa=True, wb=True)
        
        df_extra_anno = bt_extra_anno.to_dataframe(names=gene_info_cols + list(missing_anno.columns))
        
        # Combine the custom annotations with the Bed-derived ones
        annotated = pd.concat([df_extra_anno, annotated[annotated['gene_names'].notnull()]])
        annotated['GENE_NAME'] = annotated['GENE_NAME'].fillna(annotated['gene_names'])

    # Final Filter & Export
    # Keep only genes of interest
    final_df = pd.merge(annotated, df_gene_of_interest, on='GENE_NAME', how='inner')
    final_df = final_df.drop_duplicates(subset=['population', 'gene_coordinate'])
    
    # Split back into ssms/cnvs columns based on ID prefix
    final_df['ssms'] = final_df["ssms_cnvs"].apply(lambda x: x if str(x).startswith("s") else "")
    final_df['cnvs'] = final_df["ssms_cnvs"].apply(lambda x: x if str(x).startswith("c") else "")
    
    final_df = final_df[['chrom', 'start', 'end', 'population', 'ssms', 'cnvs', 'gene_coordinate', 'GENE_NAME', 'vaf', 't_alt_count', 't_ref_count']]
    final_df = final_df.sort_values(by=['population', 'chrom', 'start']).reset_index(drop=True)
    final_df = pd.merge(final_df, cp_df, on='population', how='inner') # Merge CP info to final df
    final_df.to_csv(out_path, sep='\t', index=False)
    
    return final_df

def main(sample_id, annotation_file, treesum_file, mutass_dir, ssm_file, cnvs_file, goi_file, out_dir):

    logger.info("Processing PhyloWGS outfiles %s", sample_id)
    
    os.makedirs(out_dir, exist_ok=True)
    
    # Process Tree Summary
    df_tree_sum = process_tree_summary(
        treesum_file, f"{out_dir}/{sample_id}__treesum.tsv")
    logger.info("Wrote phyloWGS tree summary: %s/%s__treesum.tsv", out_dir, sample_id)

    # Get Mutations for Top Model
    top_model = df_tree_sum['model_index'].iloc[0]
    mutass_df = get_mutation_assignments(mutass_dir, top_model)
    # Merge CP info to mutass_df
    df_tree_sum = df_tree_sum.loc[df_tree_sum['model_index']==top_model][['population'] + [col for col in df_tree_sum.columns if 'cp_' in col]]
    mutass_df = pd.merge(mutass_df, df_tree_sum, on='population', how='inner')

    # Parse Metadata
    ssm_df, cnvs_df = parse_ssm_cnv_files(ssm_file, cnvs_file)
    
    # Load Priority Annotations and Gene List
    df_anno = pd.read_csv(annotation_file, sep='\t')
    df_anno['gene_coordinate'] = 'chr' + df_anno['chrom'].astype(str) + '_' + \
                                 df_anno['start'].astype(str) + '_' + df_anno['end'].astype(str)
    df_anno = df_anno.rename(columns={'GENE_NAME':'gene_names'})[['gene_coordinate', 'gene_names', 'vaf', 't_alt_count', 't_ref_count']]
    
    # df_goi = pd.DataFrame(GENES_OF_INTEREST, columns=['GENE_NAME'])
    df_goi = pd.read_csv(goi_file, sep='\t', names=['GENE_NAME'], header=None)

    # Run Intersection and Save
    out_path = f"{out_dir}/{sample_id}__mutass_with_geneinfo.tsv"
    run_genomic_intersection(mutass_df, ssm_df, cnvs_df, df_anno, df_goi,out_path)

    logger.info("Wrote phyloWGS mutation annotation file: %s", out_path)
    logger.info("Processing phyloWGS outfiles completed for sample %s", sample_id)
    
    
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description="Process PhyloWGS output with Gene Annotations")
    
    parser.add_argument("--sample_id", type=str, help="sample ID", required=True)
    parser.add_argument("--anno_file", type=str, help="Annotation file", required=True)
    parser.add_argument("--treesum_file", type=str, help="treesum_file from PhyloWGS", required=True)
    parser.add_argument("--mutass_dir", type=str, help="Mutass dirctory from PhyloWGS", required=True)
    parser.add_argument("--ssm_file", type=str, help="SSM file", required=True)
    parser.add_argument("--cnvs_file", type=str, help="CNVs file", required=True)
    parser.add_argument("--goi_file", type=str, help="Gene of interest", required=True)
    parser.add_argument("--outdir", type=str, help="outdir", required=True)

    args = parser.parse_args()
    
    main(sample_id=args.sample_id,
                    annotation_file=args.anno_file,
                    treesum_file=args.treesum_file,
                    mutass_dir=args.mutass_dir,
                    ssm_file=args.ssm_file,
                    cnvs_file=args.cnvs_file,
                    goi_file=args.goi_file,
                    out_dir=args.outdir)






