#!/usr/bin/env nextflow

nextflow.enable.dsl=1

version = '0.2.0'

params.extra_fasta = 'repeats.fasta'
params.genome_fasta = 'http://ftp.flybase.net/releases/FB2021_04/dmel_r6.41/fasta/dmel-all-chromosome-r6.41.fasta.gz'
params.genome_gtf = 'http://ftp.flybase.net/releases/FB2021_04/dmel_r6.41/gtf/dmel-all-r6.41.gtf.gz'
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
		conda 'samtools=1.16.1'
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
		conda 'bedtools=2.30.0'
		input:
		file fai from fai_ch

		output:
		file "extra.bed" into bed_ch

		"""
		bedtools makewindows -g $fai -n 1 -i src > extra.bed
		"""
}

process bedtogenepred_get_genepred {
	conda 'ucsc-bedtogenepred=377'

	input:
	file bed from bed_ch

	output:
	file "extra.gp" into genepred_ch

	"""
	bedToGenePred $bed stdout | awk 'BEGIN{FS=OFS="\t"} {\$3="+"} 1' > extra.gp
	"""
}

process genepredtogtf_get_gtf {
	conda 'ucsc-genepredtogtf=377 gffread=0.12.7'

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
	conda 'gffread=0.12.7'

	input:
	file gtf from extra_gtf_ch
	file g_gtf from genome_gtf_ch

	output:
	file params.prefix + ".gtf" into combined_gtf_ch

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
	publishDir "results", mode: 'copy'
	conda 'repeatmasker=4.1.2 samtools=1.16.1'

	cpus 12
	time '1h'

	input:
	file genome from genome_fasta_ch
	file extra from clean_extra_fasta_ch_2

	output:
	file params.prefix + ".repeatmasked.fasta.gz" into combined_fa_ch
	file params.prefix + ".repeatmasked.gff" into combined_gff_ch

	script:
	is_gz = {assert '.gz' =~ params.genome_fasta}
	of = params.prefix + ".repeatmasked.fasta.gz"
	ofgff = params.prefix + ".repeatmasked.gff"
	if( is_gz )
		"""
		gunzip $genome -c > genome.fasta

		RepeatMasker -e ncbi -pa ${task.cpus} -s -lib $extra -no_is -gff -dir . genome.fasta

		cat $extra genome.fasta.masked  | bgzip -c > $of
		cat genome.fasta.out.gff > $ofgff
		"""
	else
		"""
		cp $genome genome.fasta

		RepeatMasker -e ncbi -pa ${task.cpus} -s -lib $extra -no_is -gff -dir . genome.fasta

		cat $extra genome.fasta.masked  | bgzip -c > $of
		cat genome.fasta.out.gff > $ofgff
		"""
}