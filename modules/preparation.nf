process GENERATE_PRIORITY_SSM {
    tag "$sample_id"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf_files), path(goi_file)

    output:
    tuple val(sample_id), path("${sample_id}_priority_ssm_file.tsv"), path("${sample_id}_priority_ssm_with_geneinfo.tsv"), emit: priority_ssm

    script:
    """
    python_generate_priority_SSM_file.py \
        --sample_id "${sample_id}" \
        --vcf_files "${vcf_files}" \
        --goi_file "${goi_file}" \
        --outdir ./
    """
}

process PREPARE_PHYLOWGS_INPUTS {
    tag "$sample_id"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(priority_ssm), path(priority_ssm_geneinfo)
    path cnvs
    path vcfs
    
    output:
    tuple val(sample_id), path("${sample_id}_ssm_data.txt"), path("${sample_id}_cnv_data.txt"), path("${sample_id}_params.json"), emit: phylowgs_inputs

    script:
    // 1. Build the dynamic strings for CNVs and VCFs
    // This creates a list like ["S1=file1", "S2=file2"]
    def cnv_args = []
    def vcf_type_args = []
    def vcf_args = []
    
    cnvs.eachWithIndex { file, i -> 
        cnv_args << "--cnvs S${i+1}=${file}" }
        
    cnvs.eachWithIndex { file, i ->
        vcf_type_args << "--vcf-type S${i+1}=${params.vcf_types}"}
        
    vcfs.eachWithIndex { file, i -> 
        // Handles the positional S1=file1 S2=file2 at the end of your command
        vcf_args << "S${i+1}=${file}"}
    
    // 2. Join the lists into a single string space-separated
    def cnv_string = cnv_args.join(' ')
    def vcf_type_string = vcf_type_args.join(' ')
    def vcf_string = vcf_args.join(' ')
    
    """
    create_phylowgs_inputs.py \
        --sex male \
        --sample-size "${params.subsample_num}" \
        --priority-ssms "${priority_ssm}" \
        --output-cnvs "${sample_id}_cnv_data.txt" \
        --output-variants "${sample_id}_ssm_data.txt" \
        --output-params "${sample_id}_params.json" \
        ${cnv_string} \
        ${vcf_type_string} \
        ${vcf_string}
    """
}