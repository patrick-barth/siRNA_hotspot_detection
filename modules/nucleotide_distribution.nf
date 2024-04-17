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

//TODO: Add version
process get_nucleotide_distribution {
    publishDir "${params.output_dir}/nucleotide_distribution", mode: 'copy', pattern: "${query.baseName}.nuc_dist.tsv"
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
        --output_percentage ${query.baseName}.nuc_percent.tsv \
        --length ${length}
    """
}

process visualize_nuc_distri {
    tag{query.simpleName}

    input:
    path(query)

    output:

    script:
    """
    """
}