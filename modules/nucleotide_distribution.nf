process get_length {
    tag{query.simpleName}

    input:
    tuple path(query), val(length)

    output:
    tuple path("${query.simpleName}.${length}.fastq"), val("${length}"), emit: fastq_reads_len_filtered
    path("${task.process}.version.txt"), 	emit: version

    script:
    """
    cutadapt -m ${length} -M ${length} -o ${query.simpleName}.${length}.fastq ${query}

    echo -e "${task.process}\tcutadapt\t\$(cutadapt --version)" > ${task.process}.version.txt
    """
}

process get_seq_only {
    tag{query.baseName}

    input:
    tuple path(query), val(length)

    output:
    tuple path("${query.baseName}.txt"), val("${length}"), emit: txt_reads_only

    """
    awk '{if(NR%4==2) print \$1}' ${query} > ${query.baseName}.txt
    """
    
}

process get_nucleotide_distribution {
    publishDir "${params.output_dir}/nucleotide_distribution", mode: 'copy', pattern: "${query.baseName}.nuc_*.tsv"
    tag{query.baseName}

    input:
    tuple path(query), val(length)

    output:
    path("${query.baseName}.nuc_dist.tsv"), emit: nuc_dist
    path("${query.baseName}.nuc_percent.tsv"), emit: nuc_percent
    path("${task.process}.version.txt"), 	emit: version

    script:
    """
    calc-nucleotide-distribution.py \
        --sequences ${query} \
        --output ${query.baseName}.nuc_dist.tsv \
        --output_percent ${query.baseName}.nuc_percent.tsv \
        --length ${length}

    echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" > ${task.process}.version.txt
    """
}

process visualize_nuc_distri {
    tag{query.simpleName}
    publishDir "${params.output_dir}/nucleotide_distribution/visualization", mode: 'copy', pattern: "${query.baseName}.pdf"

    input:
    path(query)

    output:
    path("${query.baseName}.pdf"), emit: nuc_dist_visualization
    path("${task.process}.version.txt"), 	emit: version

    script:
    """
    visualize_nuc_cov.R --input ${query} \
        --output ${query.baseName} \
        --type pdf

    echo -e "${task.process}\tR\t\$(Rscript --version 2>&1 | cut -d' ' -f5)" > ${task.process}.version.txt
    """
}