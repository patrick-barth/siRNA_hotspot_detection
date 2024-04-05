process build_index {

	input:
	path(ref)

	output:
	tuple path("${ref}"), path("${ref}.*"), emit: index
	path("${task.process}.version.txt"), 	emit: version

	"""
	bowtie2-build ${ref} ${ref}
	
	echo -e "${task.process}\tbowtie2\t\$(bowtie2-build --version | head -1 | rev | cut -f1 -d' ' | rev)" > ${task.process}.version.txt
	"""
}

process mapping {
	tag {query.simpleName}
	publishDir "${params.output_dir}/alignments", mode: 'copy', pattern: "${query.simpleName}.bam"
	publishDir "${params.output_dir}/statistics", mode: 'copy', pattern: "${query.simpleName}.statistics.txt"

	input:
	tuple path(ref), path(index)
	path(query)

	output:
	path("${query.simpleName}.bam"), 			emit: bam_alignments
	path("${query.simpleName}.statistics.txt"), emit: report
	path("${task.process}.version.txt"), 		emit: version

	script:
	def all_alignments = params.report_all_alignments ? '-a' : ''
	def some_alignments = params.max_alignments && !params.report_all_alignments ? "-k " + params.max_alignments : ''

	"""
	bowtie2 --no-unal \
		--very-sensitive \
		-L 10 \
		-q \
		${all_alignments} \
		${some_alignments} \
		-p ${task.cpus} \
		--seed 0 \
		-U ${query} \
		-x ${ref} \
		2> ${query.simpleName}.statistics.txt | samtools view -bS - > ${query.simpleName}.bam

	echo -e "${task.process}\tbowtie2\t\$(bowtie2 --version | head -1 | rev | cut -f1 -d' ' | rev)" > ${task.process}.version.txt
	echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" >> ${task.process}.version.txt
	"""
}

process filter_perfect_hits {
	tag {query.simpleName}
	publishDir "${params.output_dir}/perfect-alignments", mode: 'copy', pattern: "${bam.simpleName}.perfect_hits.bam"

	input:
	path(query)

	output:
	path("${bam.simpleName}.perfect_hits.bam"), emit: bam_perfect_alignments
    path("${task.process}.version.txt"), 		emit: version

	script:
	"""
	samtools view -F 4 -h ${query} | awk '/^@/ || !/XS:i:/' | samtools view -bS - > ${query.simpleName}.perfect_hits.bam
    
    echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" >> ${task.process}.version.txt
	"""
}