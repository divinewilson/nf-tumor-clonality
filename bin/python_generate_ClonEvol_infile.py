#!/usr/bin/env python

import os
import json
import logging
import pandas as pd


MIN_CP_THRESHOLD = 0.01
CI_MARGIN = 5
MIN_BRANCH_SIZE = 2
MAX_BRANCH_SIZE = 12

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s")

def load_sample_names(treesum_json):
    with open(treesum_json) as f:
        data = json.load(f)
    sample_names = [ f"{k} : {v.split('__cnvs.txt')[0]}" for k, v in data['params']['cnv_files'].items()]
    return ", ".join(sorted(sample_names))


def format_ci_string(row, cp_cols):
    ci_parts = []

    for cp_col in cp_cols:
        cp_val = row[cp_col]
        if cp_val > MIN_CP_THRESHOLD:
            cp_pct = round(cp_val * 100, 1)
            lower = max(0, round(cp_pct - CI_MARGIN, 1))
            upper = min(100, round(cp_pct + CI_MARGIN, 1))
            sample = f"S{cp_col.split('_')[-1]}"
            ci_parts.append(f"{sample} : {lower}-{upper}%(0.95)/p=0.000")
    return ",".join(ci_parts) if ci_parts else "NA"


def scale_branch_size(series):
    if series.nunique() == 1:
        return pd.Series([MAX_BRANCH_SIZE] * len(series), index=series.index)
    return ((series - series.min()) / (series.max() - series.min())) * (MAX_BRANCH_SIZE - MIN_BRANCH_SIZE) + MIN_BRANCH_SIZE



def generate_clonevol_infiles(sample_id, mutass_file, treesum_file, treesum_json, outdir):
    os.makedirs(outdir, exist_ok=True)
    logging.info("Processing ClonEvol outfiles sample %s", sample_id)

    sample_name_str = load_sample_names(treesum_json)

    df_mutass = pd.read_csv(mutass_file, sep="\t")
    df_tree_sum = pd.read_csv(treesum_file, sep="\t")

    if df_tree_sum.empty:
        raise ValueError("treesum_file is empty")

    top_model = df_tree_sum['model_index'].iloc[0]
    df_tree_sum = df_tree_sum[df_tree_sum['model_index'] == top_model]

    # ---------- Variant table ----------
    df_variant = (
        df_mutass[['population', 'gene_coordinate', 'GENE_NAME']]
        .drop(columns='gene_coordinate')
        .assign(population=lambda x: x['population'].astype(int))
    )

    cp_cols = [c for c in df_tree_sum.columns if c.startswith('cp_')]
    df_cp_info = df_tree_sum[['population'] + cp_cols].copy()
    df_cp_info['population'] = df_cp_info['population'].astype(int)

    df_variant = df_variant.merge(df_cp_info, on='population', how='left')
    df_variant['is.driver'] = df_variant['GENE_NAME'].notna().map({True: 'TRUE', False: 'FALSE'})
    df_variant['GENE_NAME'] = df_variant['GENE_NAME'].fillna('-')

    df_variant = df_variant.rename(
        columns={'population': 'cluster', 'GENE_NAME': 'gene'}
    )

    variant_out = os.path.join(outdir, f"{sample_id}_variant.tsv")
    df_variant.to_csv(variant_out, sep="\t", index=False)
    logging.info("Wrote ClonEvol variant file: %s", variant_out)

    # ---------- Gene events ----------
    df_gene_event = (
        df_variant[df_variant['gene'] != '-'][['cluster', 'gene']]
        .drop_duplicates()
        .groupby('cluster', as_index=False)
        .agg(event=('gene', lambda x: ",".join(sorted(x))))
        .rename(columns={'cluster': 'clone'})
    )

    # ---------- Tree table ----------
    df_tree = df_tree_sum[['child', 'parent'] + cp_cols].dropna()
    df_tree['parent'] = df_tree['parent'].replace(0, -1)

    df_tree['sample.with.nonzero.cell.frac.ci'] = df_tree.apply(
        format_ci_string, axis=1, cp_cols=cp_cols
    )

    df_tree['sample_info'] = sample_name_str
    df_tree['branch_size'] = scale_branch_size(df_tree[cp_cols].max(axis=1))

    df_tree = df_tree.drop(columns=cp_cols).rename(columns={'child': 'clone'})
    df_tree = df_tree.merge(df_gene_event, on='clone', how='left')

    tree_out = os.path.join(outdir, f"{sample_id}_tree.tsv")
    df_tree.to_csv(tree_out, sep="\t", index=False)
    logging.info("Wrote ClonEvol tree file: %s", tree_out)
    logging.info("Processing ClonEvol outfiles completed for sample %s", sample_id)


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate Clonevol input files from PhyloWGS output")
    
    parser.add_argument("--sample_id", type=str, help="sample ID", required=True)
    parser.add_argument("--mutass_file", type=str, help="Mutation Annotation file with gene info", required=True)
    parser.add_argument("--treesum_file", type=str, help="treesum_file in tsv format", required=True)
    parser.add_argument("--treesum_json", type=str, help="treesum_file in json format", required=True)
    parser.add_argument("--outdir", type=str, help="outdir", required=True)

    args = parser.parse_args()
    
    generate_clonevol_infiles(sample_id=args.sample_id,
                              mutass_file=args.mutass_file,
                              treesum_file=args.treesum_file,
                              treesum_json=args.treesum_json,
                              outdir=args.outdir)