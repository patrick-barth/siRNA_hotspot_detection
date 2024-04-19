process count_RNAs {
    tag{query.simpleName}

    input:
    path(query)

    output:
    path("${query.simpleName}.length-distribution.txt"), emit: counted_hits

    """
    uniq -c ${query} | sort -nr > ${query.simpleName}.counts.txt
	"""
}