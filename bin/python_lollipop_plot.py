#!/usr/bin/env python

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go


def genome_wide_VAF_scatterplot(sample_id, annotation_file, out_dir):
    # 1. Load your data
    df = pd.read_csv(annotation_file, sep='\t')
    
    ## prepare the data, clean VAF column
    df['vaf'] = pd.to_numeric(df['vaf'], errors='coerce')
    df = df.dropna(subset=['vaf'])
    
    ##Setup Genome-Wide X-Axis
    # Sort chromosomes naturally (chr1, chr2... chrX)
    def chrom_key(c):
        val = str(c).replace('chr', '')
        if val == 'X': return 23
        if val == 'Y': return 24
        try: return int(val)
        except: return 100
        
    chrom_list = sorted(df['chrom'].unique(), key=chrom_key)
    
    # Actual hg38 chromosome lengths (approximate in bp)
    # These define the relative width of each "bin" on the X-axis
    HG38_LENS = {
        'chr1': 248956422,  'chr2': 242193529,  'chr3': 198295559,  'chr4': 190214555,
        'chr5': 181538259,  'chr6': 170805979,  'chr7': 159345973,  'chr8': 145138636,
        'chr9': 138394717,  'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
        'chr17': 83257020,  'chr18': 80373285,  'chr19': 58617616,  'chr20': 64444167,
        'chr21': 46709983,  'chr22': 50818468,  'chrX': 156040895,  'chrY': 57227415}
    
    
    # Calculate offsets based on actual lengths
    offsets = {}
    total_offset = 0
    chrom_midpoints = {}
    
    for chrom in chrom_list:
        length = HG38_LENS.get(chrom, 100_000_000) # Default if not found
        offsets[chrom] = total_offset
        chrom_midpoints[chrom] = total_offset + (length / 2)
        total_offset += length + 5_000_000 # Scale width: Relative space = Actual bp + small padding for visual separation
        
    df['genomic_pos'] = df.apply(lambda r: r['start'] + offsets[r['chrom']], axis=1)
    
    ##Create the Figure
    genes = sorted(df['GENE_NAME'].unique())
    pops = sorted(df['population'].unique())
    fig = go.Figure()
    
    # Add a trace for every gene
    for gene in genes:
        gene_df = df[df['GENE_NAME'] == gene]
        fig.add_trace(go.Scatter(
            x=gene_df['genomic_pos'],
            y=gene_df['vaf'],
            name=gene,
            mode='markers',
            marker=dict(size=9, opacity=0.7, line=dict(width=1, color='white')),
            # Attach metadata for hover
            customdata=np.stack((gene_df['chrom'], gene_df['start'], gene_df['population'], 
                                 gene_df['t_alt_count'], gene_df['t_ref_count']), axis=-1),
            hovertemplate=("<b>Gene: " + gene + "</b><br>" +
                            "Chromosome: %{customdata[0]}<br>" +
                            "Position: %{customdata[1]:,}<br>" +
                            "VAF: %{y}<br>" +
                            "Population: %{customdata[2]}<br>" +
                            "Counts (Alt/Ref): %{customdata[3]}/%{customdata[4]}<extra></extra>")))
    
    ##Add Chromosome Background Shading
    for i, chrom in enumerate(chrom_list):
        if i % 2 == 1: # Shade every other chromosome
            start_shading = offsets[chrom]
            end_shading = offsets[chrom_list[i+1]] if i+1 < len(chrom_list) else total_offset
            fig.add_vrect(x0=start_shading, x1=end_shading, fillcolor="gray", 
                          opacity=0.1, layer="below", line_width=0)
            
    ##Build the Filter Dropdown
    buttons = []
    
    # "All Populations" Option
    buttons.append(dict(
        method="update", label="All Populations",
        args=[{"x": [df[df['GENE_NAME'] == g]['genomic_pos'] for g in genes],
                "y": [df[df['GENE_NAME'] == g]['vaf'] for g in genes],
                "customdata": [np.stack((df[df['GENE_NAME'] == g]['chrom'], 
                                        df[df['GENE_NAME'] == g]['start'], 
                                        df[df['GENE_NAME'] == g]['population'],
                                        df[df['GENE_NAME'] == g]['t_alt_count'], 
                                        df[df['GENE_NAME'] == g]['t_ref_count']), axis=-1) for g in genes],
                "visible": [True] * len(genes)}]))
    
    # Specific Population Options
    for p in pops:
        x_data, y_data, custom_data, visibility = [], [], [], []
        for g in genes:
            sub = df[(df['GENE_NAME'] == g) & (df['population'] == p)]
            x_data.append(sub['genomic_pos'])
            y_data.append(sub['vaf'])
            custom_data.append(np.stack((sub['chrom'], sub['start'], sub['population'], 
                                         sub['t_alt_count'], sub['t_ref_count']), axis=-1) if not sub.empty else [])
            visibility.append(len(sub) > 0)
        
        buttons.append(dict(
            method="update", label=f"Population {p}",
            args=[{"x": x_data, "y": y_data, "customdata": custom_data, "visible": visibility}]))
        
    ##Final Layout Customization
    fig.update_layout(
        xaxis=dict(
            tickvals=list(chrom_midpoints.values()),
            ticktext=list(chrom_midpoints.keys()),
            title="Chromosome",
            showgrid=False),
        yaxis_title="VAF",
        template="plotly_white",
        updatemenus=[{"buttons": buttons, "direction": "down", "showactive": True, "x": 0.01, "y": 1.15}],
        legend=dict(title="Genes", orientation="v", y=1,yanchor="top"))
    
    
    fig.write_html('%s/%s__pops_VAFs_interactive.html'%(out_dir, sample_id))
    
    return None



if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description="Genome wide VAF scatter plot")
    
    parser.add_argument("--sample_id", type=str, help="sample ID", required=True)
    parser.add_argument("--anno_file", type=str, help="Annotation file", required=True)
    parser.add_argument("--outdir", type=str, help="outdir", required=True)

    args = parser.parse_args()
    
    genome_wide_VAF_scatterplot(sample_id=args.sample_id,
                                annotation_file=args.anno_file,
                                out_dir=args.outdir)