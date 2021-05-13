#!/usr/bin/env nextflow

nextflow.enable.dsl=1

version = '1.0'

params.extra_fasta = '../data-genome/dmel_repbase_lib.fasta'
params.genome_fasta = 'ftp://ftp.flybase.net/genomes/dmel/current/fasta/dmel-all-chromosome-r6.39.fasta.gz'
params.genome_gtf = 'ftp://ftp.flybase.net/releases/current/dmel_r6.39/gtf/dmel-all-r6.39.gtf.gz'
params.prefix = './plus-repeats'


extra_fasta_ch = Channel.fromPath(params.extra_fasta)
genome_fasta_ch = Channel.fromPath(params.genome_fasta)
genome_gtf_ch = Channel.fromPath(params.genome_gtf)

process sed_strip_rpm_weirdness {
	input:
	file extra_fasta from extra_fasta_ch

	output:
	file "clean_extra.fasta" into clean_extra_fasta_ch, clean_extra_fasta_ch_2

	"""
	sed -r '/^>/s/#.+//' $extra_fasta > clean_extra.fasta
	"""
}


process samtools_faidx_get_seqlens {
		input:
		file fasta from clean_extra_fasta_ch

		output:
		file "extra.fasta.fai" into fai_ch

    """
    samtools faidx $fasta && \
			mv '$fasta'.fai extra.fasta.fai
    """
}

process bedtools_makewindows_get_bed {
		input:
		file fai from fai_ch

		output:
		file "extra.bed" into bed_ch

		"""
		bedtools makewindows -g $fai -n 1 -i src > extra.bed
		"""
}

process bedtogenepred_get_genepred {
	conda 'ucsc-bedtogenepred'

	input:
	file bed from bed_ch

	output:
	file "extra.gp" into genepred_ch

	"""
	bedToGenePred $bed stdout | awk 'BEGIN{FS=OFS="\t"} {\$3="+"} 1' > extra.gp
	"""
}

process genepredtogtf_get_gtf {
	conda 'ucsc-genepredtogtf gffread'

	input:
	file gp from genepred_ch

	output:
	file "extra.gtf" into extra_gtf_ch

	"""
	genePredToGtf file $gp stdout | gffread -T > "extra.gtf"
	"""
}


process combine_gtfs {
	publishDir 'results', mode: 'copy'
	conda 'gffread'

	input:
	file gtf from extra_gtf_ch
	file g_gtf from genome_gtf_ch

	output:
	file params.prefix + ".gtf"

	script:
	is_gz = {assert '.gz' =~ params.genome_gtf}
	of = params.prefix + ".gtf"
	if( is_gz )
		"""
		cat $gtf <(zcat $g_gtf) | gffread -T > $of
		"""
	else
		"""
		cat $gtf <(cat $g_gtf) | gffread -T > $of
		"""
}

process repeatmasker_mask_extra {
	publishDir "results", pattern: '*.fasta.gz', mode: 'copy'
	conda 'repeatmasker'

	cpus 12

	input:
	file genome from genome_fasta_ch
	file extra from clean_extra_fasta_ch_2

	output:
	file params.prefix + ".fasta.gz"

	script:
	is_gz = {assert '.gz' =~ params.genome_fasta}
	of = params.prefix + ".fasta.gz"
	if( is_gz )
		"""
		gunzip $genome -c > genome.fasta

		RepeatMasker -e ncbi -pa ${task.cpus} -s -lib $extra -no_is -nolow -dir . genome.fasta

		cat $extra genome.fasta.masked  | gzip -c > $of
		"""
	else
		"""
		cp $genome genome.fasta

		RepeatMasker -e ncbi -pa ${task.cpus} -s -lib $extra -no_is -nolow -dir . genome.fasta

		cat $extra genome.fasta.masked  | gzip -c > $of
		"""
}
