process MAKE_CLONEVOL_INFILES {
    tag "$sample_id"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(mutass_geneinfo), path(treesum_tsv), path(treesum_json)

    output:
    tuple val(sample_id), path("${sample_id}_variant.tsv"), path("${sample_id}_tree.tsv"), emit: clonevol_inputs

    script:
    """
    python_generate_ClonEvol_infile.py \
        --sample_id "${sample_id}" \
        --mutass_file "${mutass_geneinfo}" \
        --treesum_file "${treesum_tsv}" \
        --treesum_json "${treesum_json}"  \
        --outdir ./
    """
}

process CLONEVOL {
    tag "$sample_id"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(variant), path(tree)
    
    output:
    path ("*")
    
    script:
    """
    Rscript_clonevol.R "${tree}" "${variant}" "${sample_id}"_clonality_tree.pdf "${params.tree_width}" "${params.tree_height}"
    """
    
}