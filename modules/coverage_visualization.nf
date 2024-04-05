process index_alignments {
	tag {query}

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
	tuple path("${bam.simpleName}.for.bed"), path("${bam.simpleName}.rev.bed"), emit: bed_coverage

	script:
	"""
    genomeCoverageBed -ibam ${bam} -bg -strand + > ${bam.simpleName}.for.bed
	genomeCoverageBed -ibam ${bam} -bg -strand - > ${bam.simpleName}.rev.bed
	"""
}

//TODO: Check how the script needs to be adapted for forward reverse split
//TODO: get version
process find_potential_hotspots {
	publishDir "${params.output}/potential-hotspots", mode: 'copy'

	input:
	tuple path(query_for), path(query_rev)
    file(reference)

	output:
	path("${query_for.simpleName}.hotspots.txt")

	script:
	"""
	hotspot-detection.py --input ${bed} --cutoff ${params.cutoff} > ${bed.simpleName}.hotspots.txt
	"""
}

process generate_R_plots {
	publishDir "${params.output}", mode: 'copy', pattern: 'igv-session.xml'

	input:
	tuple path(query_for), path(query_rev)
	path(ref)

	"""
	echo "hi"
	"""
}