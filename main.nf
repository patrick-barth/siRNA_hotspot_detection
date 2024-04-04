#!/usr/bin/env nextflow

import groovy.json.JsonOutput // used for parameter output

nextflow.enable.dsl=2

include{
    collect_metadata
    get_md5sum
    multiqc
    collect_versions
} from './modules/default_processes.nf'

include{
    quality_control
    quality_control_2
    adapter_removal
    length_filter
    filter_bacterial_contamination
} from './modules/read_processing.nf'

/*
 * Prints help and exits workflow afterwards when parameter --help is set to true
 */

if ( params.help ) {
    help = """main.nf: A description of your script and maybe some examples of how
                |                to run the script
                |Required arguments:
                |   --reads         Location of the input file file (FASTQ).
                |
                |Optional arguments:
                |   --min_length    Minimum length for reads after adapter trimming.
                |                   [default: ${params.min_length}]
                |  -w            The NextFlow work directory. Delete the directory once the process
                |                is finished [default: ${workDir}]""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

//preparation for workflow

/*
 * Welcome log to be displayed before workflow
 */
log.info """\
        ${params.manifest.name} v${params.manifest.version}
        ==========================
        reads        : ${params.reads}
        reference    : ${params.reference}
        output to    : ${params.output_dir}
        --
        run as       : ${workflow.commandLine}
        started at   : ${workflow.start}
        config files : ${workflow.configFiles}
        """
        .stripIndent()

//essential input files
input_reads     = Channel.fromPath( params.reads )
input_reference = Channel.fromPath( params.reference )

// Collect all input files
input_files = input_reads.concat(input_reference)
                .flatten().toList()

/*
 * Starting subworkflow descriptions
 */
kraken_db = file(params.kraken_db_dir).toAbsolutePath()
workflow preprocessing {
    take: 
        input_reads
        kraken_db
    main:
        quality_control(input_reads)
        adapter_removal(input_reads)
        length_filter(adapter_removal.out.fastq_length_filtered)
        if(params.filter_bac_cont){
            filter_bacterial_contamination(length_filter.out.fastq_length_filtered,kraken_db)
        }
        processed_reads = params.filter_bac_cont ? filter_bacterial_contamination.out.fastq_bac_cont_filtered : length_filter.out.fastq_length_filtered
        quality_control_2(processed_reads)

        versions = quality_control.out.version.first()
                    .concat(quality_control_2.out.version.first())
                    .concat(adapter_removal.out.version.first())
                    .concat(length_filter.out.version.first())

        versions = params.filter_bac_cont ? filter_bacterial_contamination.out.version : versions

    emit:
        //data for multiqc
        multiqc_quality_control                     = quality_control.out
        multiqc_quality_control_post_preprocessing  = quality_control_2.out
        multiqc_adapter_removal                     = adapter_removal.out.report_trimming
        multiqc_bac_filter  = params.filter_bac_cont ? filter_bacterial_contamination.out.report : Channel.empty()

        versions = versions

        // data for downstream processes
        fastq_reads = processed_reads
}



/*
 * Actual workflow connecting subworkflows
 */
workflow {
    preprocessing(input_reads)

    // Collect metadata
    collect_metadata()
    get_md5sum(input_files)
    collect_versions(collect_metadata.out.version
                        .concat(get_md5sum.out.version)
                        .unique()
                        .flatten().toList()
    )
}



/*
 * Prints complection status to command line
 */
workflow.onComplete{
	println "Pipeline completed at: $workflow.complete"
	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
	println "Something went wrong :("
	println "Pipeline execution stopped with following error message: ${workflow.errorMessage}"
}