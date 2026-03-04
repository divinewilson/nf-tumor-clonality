#!/usr/bin/env python

import numpy as np
import pandas as pd
import plotly.graph_objects as go


def pops_bubble_plot(sample_id, annotation_file, out_dir):
    ##Load and Prepare Data
    df = pd.read_csv(annotation_file, sep='\t')
    df['vaf'] = pd.to_numeric(df['vaf'], errors='coerce')
    df['t_alt_count'] = pd.to_numeric(df['t_alt_count'], errors='coerce')
    df['t_ref_count'] = pd.to_numeric(df['t_ref_count'], errors='coerce')
    df['total_depth'] = df['t_alt_count'] + df['t_ref_count']
    df = df.dropna(subset=['vaf', 'population', 'GENE_NAME', 'total_depth'])
    
    ##Weighted Aggregation Logic
    def weighted_mean(vafs, weights):
        return np.average(vafs, weights=weights)
    
    agg_data = []
    for (pop, gene), group in df.groupby(['population', 'GENE_NAME']):
        count = len(group)
        w_vaf = weighted_mean(group['vaf'], group['total_depth'])
        vaf_std = group['vaf'].std() if count > 1 else 0  # Standard Deviation = Heterogeneity
        med_depth = group['total_depth'].median()
        
        agg_data.append({'population': pop,
                         'GENE_NAME': gene,
                         'w_vaf': w_vaf,
                         'vaf_std': vaf_std,
                         'mut_count': count,
                         'med_depth': med_depth})
        
    agg_df = pd.DataFrame(agg_data)
    agg_df['pop_jittered'] = agg_df['population'] + np.random.uniform(-0.15, 0.15, size=len(agg_df))
    
    ##Create Heterogeneity-Aware Bubble Plot
    fig = go.Figure()
    genes = sorted(agg_df['GENE_NAME'].unique())
    
    for gene in genes:
        gdf = agg_df[agg_df['GENE_NAME'] == gene]
        # Size = Mutation Count
        bubble_sizes = np.interp(gdf['mut_count'], (agg_df['mut_count'].min(), agg_df['mut_count'].max()), (12, 55))
        
        fig.add_trace(go.Scatter(x=gdf['pop_jittered'],
                                 y=gdf['w_vaf'],
                                 name=gene,
                                 mode='markers',
                                 marker=dict(size=bubble_sizes, opacity=0.8,
                                             # color=gdf['med_depth'], colorscale='Viridis',
                                             showscale=False, line=dict(width=2, color='#4C4B4B')),
                                 error_y=dict(type='data', array=gdf['vaf_std'], visible=True,
                                              thickness=1, width=5, color='rgba(50, 50, 50, 0.4)'),
                                 customdata=np.stack((gdf['GENE_NAME'],
                                                      gdf['mut_count'],
                                                      gdf['med_depth'],
                                                      gdf['vaf_std']), axis=-1),
                                 hovertemplate=("<b>%{customdata[0]}</b><br>Pop: %{x:.0f}<br>Weighted VAF: %{y:.3f}<br>" +
                                                "Heterogeneity (±Std): %{customdata[3]:.3f}<br>" +
                                                "Mutations: %{customdata[1]:.0f}<br>Median Depth: %{customdata[2]:.0f}<extra></extra>")))
    
    
    ##Layout
    fig.update_layout(xaxis=dict(title="Population", tickvals=sorted(agg_df['population'].unique()),
                                 ticktext=[f"Pop {int(p)}" for p in sorted(agg_df['population'].unique())], gridcolor='whitesmoke'),
                      yaxis_title="Weighted Mean VAF", template="plotly_white",
                      legend=dict(title="Genes", orientation="v", y=1,yanchor="top"),
                      margin=dict(l=50, r=50, t=50, b=50), height=650)
    
    fig.write_html('%s/%s__pops_bubbleplot_interactive.html'%(out_dir, sample_id))
    
    return None

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description="Genome wide VAF scatter plot")
    
    parser.add_argument("--sample_id", type=str, help="sample ID", required=True)
    parser.add_argument("--anno_file", type=str, help="Annotation file", required=True)
    parser.add_argument("--outdir", type=str, help="outdir", required=True)

    args = parser.parse_args()
    
    pops_bubble_plot(sample_id=args.sample_id,
                    annotation_file=args.anno_file,
                    out_dir=args.outdir)
