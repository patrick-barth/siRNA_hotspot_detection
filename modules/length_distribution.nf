process get_length_distribution {
    tag{query.simpleName}

    input:
    path(query)

    output:
    path("${query.simpleName}.length-distribution.txt"), emit: distribution
    path("${task.process}.version.txt"), 	emit: version

    """
	awk '{if(NR%4==2) print length(\$1)}' ${query} | sort -n | uniq -c > ${query.simpleName}.length-distribution.txt

    echo -e "${task.process}\tawk\t\$(awk -W version | head -1)" > ${task.process}.version.txt
	"""
}

//TODO: Add container
process calc_percent {
    tag{query.simpleName}
    publishDir "${params.output_dir}/length-distribution", mode: 'copy', pattern: "${query.simpleName}.length_dist.txt"

    input:
    path(query)

    output:
    path("${query.simpleName}.length_dist.txt"), emit: percent
    path("${task.process}.version.txt"), 	emit: version

    script:
	"""
	percentage-for-length-distribution.py --input ${query} --output ${query.simpleName}.length_dist.txt

    echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" > ${task.process}.version.txt
	"""
}

process visualize_length_dist {
    tag{query.simpleName}
    publishDir "${params.output_dir}/length-distribution/visualization", mode: 'copy', pattern: "${query.simpleName}.pdf"

    input:
    path(query)

    output:
    path("${query.simpleName}.pdf"), emit: len_dist_visualization
    path("${task.process}.version.txt"), 	emit: version

    script:
    """
    visualize_len_dis.R --input ${query} \
        --output ${query.simpleName} \
        --type pdf

        echo -e "${task.process}\tR\t\$(Rscript --version 2>&1 | cut -d' ' -f5)" > ${task.process}.version.txt
    """
}