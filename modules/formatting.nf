process WRITE_RESULTS {
    tag "$sample_id"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(trees_zip)

    output:
    tuple val(sample_id), path("${sample_id}.treesumm.json.gz"), path("${sample_id}.muts.json.gz"), path("${sample_id}.mutass.zip"), emit: phylowgs_outputs

    script:
    """
    write_results.py \
        "${sample_id}" \
        "${trees_zip}" \
        "${sample_id}.treesumm.json.gz" \
        "${sample_id}.muts.json.gz" \
        "${sample_id}.mutass.zip"
    """
}

process FORMAT_RESULTS {
    tag "$sample_id"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(ssm_geneinfo), path(ssm), path(cnv), \
          path(treesum_json), path(muts_json), path(mutass_zip), path(goi_file)
    
    output:
    tuple val(sample_id), path("${sample_id}__mutass_with_geneinfo.tsv"), \
            path("${sample_id}__treesum.tsv"), path("${sample_id}.treesumm.json"), emit: phylowgs_format_outputs
    
    script:
    """
    # Decompress & Unzip
    gunzip -c "${treesum_json}" > "${treesum_json.baseName}"
    gunzip -c "${muts_json}" > "${muts_json.baseName}"
    mkdir -p mutass
    unzip -o "${mutass_zip}" -d mutass
    
    # Final Formatting
    python_generate_out_treefiles.py \
        --sample_id "${sample_id}" \
        --anno_file "${ssm_geneinfo}" \
        --treesum_file "${treesum_json.baseName}" \
        --mutass_dir mutass \
        --ssm_file "${ssm}" \
        --cnvs_file "${cnv}" \
        --goi_file "${goi_file}" \
        --outdir ./
    """
}


    
    
