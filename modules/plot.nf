process GENOME_VAF_plot {
    tag "$sample_id"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(mutass_geneinfo), path(treesum_tsv), path(treesum_json)
    
    output:
    path ("*")
    
    script:
    """
    python_lollipop_plot.py \
        --sample_id "${sample_id}" \
        --anno_file "${mutass_geneinfo}" \
        --outdir ./
    """
    
}

process POPS_VAF_plot {
    tag "$sample_id"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(mutass_geneinfo), path(treesum_tsv), path(treesum_json)
    
    output:
    path ("*")
    
    script:
    """
    python_bubble_plot.py \
        --sample_id "${sample_id}" \
        --anno_file "${mutass_geneinfo}" \
        --outdir ./
    """
    
}