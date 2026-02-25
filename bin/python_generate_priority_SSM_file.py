#!/usr/bin/env python

from pybedtools import BedTool
import pandas as pd
import glob
import sys
import os


## for nextflow
script_path = os.path.realpath(__file__)
bin_dir = os.path.dirname(script_path)
project_root = os.path.dirname(bin_dir)
utils_dir = os.path.join(project_root, 'utils')
sys.path.append(utils_dir)


def priority_ssm(sample_vcffile_str, sample_id, goi_file, outdir):

    gene_info_file = os.path.join(utils_dir, 'gencode.v45.chr_patch_hapl_scaff.basic.annotation.bed')
    df_gene = pd.read_csv(gene_info_file, sep='\t')
    df_gene['CHROM'] = df_gene['CHROM'].str.extract(r'chr(\w+)')
    df_gene_col= df_gene.columns.tolist()
    df_gene = BedTool.from_dataframe(df_gene) # Convert df to BedTool
    
    sample_file_str = sample_vcffile_str
    
    sample_vcf_file_list = []
    # sample_anno_info_list = []
    for i in sample_file_str.split(' '):
        sample_vcf_file_list.append(i.split('=')[-1])
        # sample_anno_info_list.append('%s_annotation_info.tsv'%(i.split('=')[-1].split('_reformated.vcf')[0]))
        
        
    def vcf_to_bed(vcf_path):
        bed_records = []
        with open(vcf_path, "r") as f:
            for line in f:
                if line.startswith("#"):  # skip header
                    continue
                parts = line.strip().split("\t")
                chrom = parts[0]
                pos = int(parts[1])   # VCF POS (1-based)
                ref = parts[3]
                alt = parts[4]
                # BED format is 0-based start, half-open intervals
                start = pos - 1
                end = pos - 1 + len(ref)  # handle indels properly
                ssm = chrom + '_' + str(pos)
                info = parts[7]
                vaf = info.split(';')[0].split('=')[-1]
                t_alt_count = info.split(';')[1].split('=')[-1]
                t_ref_count = info.split(';')[2].split('=')[-1]
                # func_refGene = parts[8]
                # caller = parts[9]
                bed_records.append([chrom, start, end, ref, alt, ssm, vaf, t_alt_count, t_ref_count])
        # Convert to DataFrame
        bed_df = pd.DataFrame(bed_records, columns=["chrom", "start", "end", "ref", "alt", "ssm", 'vaf', 't_alt_count', 't_ref_count'])
        return bed_df
    
    
    dfs_vcf = []
    for vcf_file in sample_vcf_file_list:
        df_vcf = vcf_to_bed(vcf_file)
        # print (df_vcf)
        dfs_vcf.append(df_vcf)
    
    df_vcf = pd.concat(dfs_vcf).reset_index(drop=True)
    # print (df_vcf)
    
    df_vcf_cols = df_vcf.columns.tolist()
    df_vcf_bedtools = BedTool.from_dataframe(df_vcf) # # Convert df to BedTool
    bt_annotated = df_gene.intersect(df_vcf_bedtools, wa=True, wb=True) # Intersect (add gene annotation to each mutation)
    
    df_annotated_cols = df_gene_col + df_vcf_cols
    df_vcf_annotated = bt_annotated.to_dataframe(names=df_annotated_cols)
    # print (df_vcf_annotated)

    
    ### # Keep only genes of interest
    df_gene_of_intreset = pd.read_csv(goi_file, sep='\t', names=['GENE_NAME'], header=None)
    df_filtered = pd.merge(df_vcf_annotated, df_gene_of_intreset, on='GENE_NAME', how='inner')
    df_filtered = df_filtered.drop_duplicates(subset=['ssm'])
    # print (df_filtered)
    ###########################################
          
    
    df_filtered.to_csv('%s/%s_priority_ssm_with_geneinfo.tsv'%(outdir, sample_id), sep='\t', index=False)
    
    df_filtered = df_filtered[['ssm']]
    df_filtered.to_csv('%s/%s_priority_ssm_file.tsv'%(outdir, sample_id), sep='\t', index=False, header=None)
    
    print ('Done')
    return None

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--sample_id", type=str, help="sample ID")
    parser.add_argument("--vcf_files", type=str, help="string of vcf file")
    parser.add_argument("--goi_file", type=str, help="gene of interest file")
    parser.add_argument("--outdir", type=str, help="outdir")
    
    args = parser.parse_args()
    
    priority_ssm(sample_id=args.sample_id,
                            sample_vcffile_str=args.vcf_files,
                            goi_file=args.goi_file,
                            outdir=args.outdir)



