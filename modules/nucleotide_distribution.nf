//TODO: Add version
process get_length {
    tag{query.simpleName}

    input:
    tuple path(query), val(length)

    output:
    tuple path("${query.simpleName}.${length}.fastq"), val("${length}"), emit: fastq_reads_len_filtered

    script:
    """
    cutadapt -m ${length} -M ${length} -o ${query.simpleName}.${length}.fastq ${query}
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

//TODO: Add version
process get_nucleotide_distribution {
    publishDir "${params.output_dir}/nucleotide_distribution", mode: 'copy', pattern: "${query.baseName}.nuc_*.tsv"
    tag{query.baseName}

    input:
    tuple path(query), val(length)

    output:
    path("${query.baseName}.nuc_dist.tsv"), emit: nuc_dist
    path("${query.baseName}.nuc_percent.tsv"), emit: nuc_percent

    script:
    """
    calc-nucleotide-distribution.py \
        --sequences ${query} \
        --output ${query.baseName}.nuc_dist.tsv \
        --output_percent ${query.baseName}.nuc_percent.tsv \
        --length ${length}
    """
}

//TODO: Add version of R
process visualize_nuc_distri {
    tag{query.simpleName}
    publishDir "${params.output_dir}/nucleotide_distribution/visualization", mode: 'copy', pattern: "${query.baseName}.pdf"

    input:
    path(query)

    output:
    path("${query.baseName}.pdf"), emit: nuc_dist_visualization

    script:
    """
    visualize_nuc_cov.R --input ${query} \
        --output ${query.baseName} \
        --type pdf
    """
}