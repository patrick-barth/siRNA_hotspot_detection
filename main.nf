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

include{
    build_index
    mapping
    extract_perfect_hits
} from './modules/alignment.nf'

include{
    sort_alignments
    index_alignments
    transform_to_bed
    find_potential_hotspots
    generate_R_plots
} from './modules/coverage_visualization.nf'

include{
    get_read_names
    collect_reads
    extract_reads
} from './modules/read_extraction.nf'

include{
    count_RNAs
} from './modules/get_best_hits.nf'

include{
    get_length_distribution
    calc_percent
    visualize_length_dist
} from './modules/length_distribution.nf'

include{
    get_length
    get_seq_only
    get_nucleotide_distribution
    visualize_nuc_distri
} from './modules/nucleotide_distribution.nf'

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

//TODO: make correct channel for all lengths
tmp_length = Channel.of(params.length_of_interest)

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
        length_filter(adapter_removal.out.fastq_trimmed)
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
        multiqc_quality_control                     = quality_control.out.output
        multiqc_quality_control_post_preprocessing  = quality_control_2.out.output
        multiqc_adapter_removal                     = adapter_removal.out.report
        multiqc_bac_filter  = params.filter_bac_cont ? filter_bacterial_contamination.out.report : Channel.empty()

        versions = versions

        // data for downstream processes
        fastq_reads = processed_reads
}

workflow alignment {
    take:
        reads
        reference
    main:
        build_index(reference)
        mapping(build_index.out.index.first(),
        reads)
        extract_perfect_hits(mapping.out.bam_alignments)

        versions = build_index.out.version.first()
            .concat(mapping.out.version.first())
            .concat(extract_perfect_hits.out.version.first())

    emit:
        all_alignments = mapping.out.bam_alignments
        perfect_alignments = extract_perfect_hits.out.bam_alignments

        versions = versions
        report = mapping.out.report
}

workflow coverage_visualization {
    take:
        alignments
        reference
    main:
        sort_alignments(alignments)
        index_alignments(sort_alignments.out.bam_sorted)
        transform_to_bed(index_alignments.out.bam_bai_index)
        find_potential_hotspots(transform_to_bed.out.bed_coverage)
        generate_R_plots(transform_to_bed.out.bed_coverage)

        versions = sort_alignments.out.version.first()
            .concat(index_alignments.out.version.first())
            .concat(transform_to_bed.out.version.first())
            .concat(find_potential_hotspots.out.version.first())
            .concat(generate_R_plots.out.version.first())
    
    emit:
        versions = versions
}

workflow read_extraction {
    take:
        alignments
        reads
    main:
        get_read_names(alignments)
        collect_reads(reads.collect()) //If several samples should be analysed in parallel this needs to be adapted
        extract_reads(get_read_names.out.names
            .combine(collect_reads.out.reads))

        versions = get_read_names.out.version.first()
            .concat(extract_reads.out.version.first())

    emit:
        reads = extract_reads.out.reads
        versions = versions
}

workflow get_best_hits {
    take:
        reads
    main:
        get_seq_only(reads)
        count_RNAs(get_seq_only.out.txt_reads_only)
}

workflow length_distribution {
    take:
        reads
    main:
        get_length_distribution(reads)
        calc_percent(get_length_distribution.out.distribution)
        visualize_length_dist(calc_percent.out.percent)

        versions = get_length_distribution.out.version.first()
            .concat(calc_percent.out.version.first())
            .concat(visualize_length_dist.out.version.first())

    emit:
        percentages = calc_percent.out.percent
        versions = versions
}

workflow nucleotide_distribution {
    take:
        reads
        length
    main:
        get_length(reads.combine(length))
        get_seq_only(get_length.out.fastq_reads_len_filtered)
        get_nucleotide_distribution(get_seq_only.out.txt_reads_only)
        visualize_nuc_distri(get_nucleotide_distribution.out.nuc_percent)

        versions = get_length.out.version.first()
            .concat(get_nucleotide_distribution.out.version.first())
            .concat(visualize_nuc_distri.out.version.first())

    emit:
        versions = versions
}

/*
 * Actual workflow connecting subworkflows
 */
workflow {
    preprocessing(input_reads,kraken_db)
    alignment(preprocessing.out.fastq_reads,input_reference)
    read_extraction(alignment.out.all_alignments,
        preprocessing.out.fastq_reads)
    length_distribution(read_extraction.out.reads)
    coverage_visualization(alignment.out.perfect_alignments,input_reference)
    nucleotide_distribution(read_extraction.out.reads,tmp_length)

    // Collect metadata
    collect_metadata()
    get_md5sum(input_files)
    collect_versions(collect_metadata.out.version
        .concat(preprocessing.out.versions)
        .concat(alignment.out.versions)
        .concat(coverage_visualization.out.versions)
        .concat(get_md5sum.out.version)
        .concat(read_extraction.out.versions)
        .concat(length_distribution.out.versions)
        .concat(nucleotide_distribution.out.versions)
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