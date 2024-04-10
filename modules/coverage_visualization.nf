process sort_alignments{
	tag {query.simpleName}

	input:
	path(query)

	output:
	path("${query.baseName}.sorted.bam"), 	emit: bam_sorted
	path("${task.process}.version.txt"), 	emit: version

	"""
	samtools sort ${query} -o ${query.baseName}.sorted.bam

	echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" >> ${task.process}.version.txt
	"""
}

process index_alignments {
	tag {query.simpleName}

	input:
	path(query)

	output:
	tuple path("${query}"), path("${query}.bai"),   emit: bam_bai_index
    path("${task.process}.version.txt"), 		    emit: version

	script:
	"""
	samtools index ${query}

    echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" > ${task.process}.version.txt
	"""
}

//TODO: get version
process transform_to_bed {
	tag {query}
	publishDir "${params.output_dir}/bed-files", mode: 'copy'

	input:
	tuple path(query), path(index) 

	output:
	tuple path("${query.simpleName}.for.bed"), path("${query.simpleName}.rev.bed"), emit: bed_coverage

	script:
	"""
    genomeCoverageBed -ibam ${query} -bga -strand + > ${query.simpleName}.for.bed
	genomeCoverageBed -ibam ${query} -bga -strand - > ${query.simpleName}.rev.bed
	"""
}

//TODO: Check how the script needs to be adapted for forward reverse split
//TODO: get version
process find_potential_hotspots {
	publishDir "${params.output_dir}/potential-hotspots", mode: 'copy', pattern: "*.hotspots.*.txt"

	input:
	tuple path(query_for), path(query_rev)

	output:
	tuple path("${query_for.simpleName}.hotspots.for.txt"), path("${query_rev.simpleName}.hotspots.rev.txt")

	script:
	def percentile = params.percentile_hotspot ? '--percentile ' + params.percentile_hotspot : ''
	def cutoff = params.cutoff_hotspot ? '--cutoff ' + params.cutoff_hotspot : ''
	"""
	hotspot-detection.py --input ${query_for} ${percentile} ${cutoff} > ${query_for.simpleName}.hotspots.for.txt
	hotspot-detection.py --input ${query_rev} --cutoff ${cutoff} > ${query_rev.simpleName}.hotspots.rev.txt
	"""
}

process generate_R_plots {
	tag {query_for.simpleName}
	publishDir "${params.output_dir}", mode: 'copy', pattern: 'igv-session.xml'

	input:
	tuple path(query_for), path(query_rev)

	"""
	visualize_coverage.R \
		--input_for ${query_for} \
		--input_rev ${query_rev} \
		--output ${query_for.simpleName}
	"""
}