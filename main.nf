nextflow.enable.dsl = 2

// 1. Import Modules
//include { GENERATE_PRIORITY_SSM}                           from './modules/preparation.nf'
include { GENERATE_PRIORITY_SSM; PREPARE_PHYLOWGS_INPUTS } from './modules/preparation.nf'
include { RUN_MULTIEVOLVE }                                from './modules/phylowgs.nf'
include { WRITE_RESULTS; FORMAT_RESULTS }                  from './modules/formatting.nf'
include {MAKE_CLONEVOL_INFILES; CLONEVOL }                 from './modules/clonevol.nf'

// 2. Define Workflow
workflow {

    // Create a channel for input data
    cnv_files_ch = Channel.fromPath(params.cnv_files).collect().map { it.sort() }
    vcf_files_ch = Channel.fromPath(params.vcf_files).collect().map { it.sort() }
    generate_priority_ssm_ch = vcf_files_ch.map { filenames -> [ params.sample_id, filenames, file(params.goi_file) ] }
    //generate_priority_ssm_ch.view()

    // --- EXECUTION ---

    GENERATE_PRIORITY_SSM(generate_priority_ssm_ch)
    
    PREPARE_PHYLOWGS_INPUTS(
                            GENERATE_PRIORITY_SSM.out.priority_ssm,
                            cnv_files_ch,
                            vcf_files_ch )

    RUN_MULTIEVOLVE(PREPARE_PHYLOWGS_INPUTS.out.phylowgs_inputs)

    WRITE_RESULTS(RUN_MULTIEVOLVE.out.trees)
    
    //make channel for data formating
    priority_ssm_ch = GENERATE_PRIORITY_SSM.out.priority_ssm
        .map { sample_id, tsv, ssm_geneinfo -> [ sample_id, ssm_geneinfo ] }
    ssm_cnv_ch = PREPARE_PHYLOWGS_INPUTS.out.phylowgs_inputs
        .map { sample_id, ssm, cnv, params_json -> [ sample_id, ssm, cnv ] }
    combined_ch = priority_ssm_ch.join(ssm_cnv_ch)
    format_result_ch = combined_ch.join(WRITE_RESULTS.out.phylowgs_outputs)
    format_result_ch_ = format_result_ch.combine( Channel.fromPath(params.goi_file) )
    
    
    FORMAT_RESULTS(format_result_ch_)
    
    MAKE_CLONEVOL_INFILES(FORMAT_RESULTS.out.phylowgs_format_outputs)
    
    CLONEVOL(MAKE_CLONEVOL_INFILES.out.clonevol_inputs)
}