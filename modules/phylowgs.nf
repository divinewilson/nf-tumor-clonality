process RUN_MULTIEVOLVE {
    tag "$sample_id"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(ssm), path(cnv), path(params_json)

    output:
    tuple val(sample_id), path("trees.zip"), emit: trees

    script:
    """
    mkdir temp_dir
    multievolve.py \
        --num-chains 4 \
        --ssms "${ssm}" \
        --cnvs "${cnv}" \
        --output-dir ./ \
        --params "${params_json}" \
        --tmp-dir temp_dir \
        --burnin-samples 500 \
        --mcmc-samples 1450
    """
}