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
		conda 'samtools'
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
		conda 'bedtools'
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
	publishDir "results", pattern: '*.fasta.gz', mode: 'copy'
	conda 'repeatmasker'

	cpus 12

	input:
	file genome from genome_fasta_ch
	file extra from clean_extra_fasta_ch_2

	output:
	file params.prefix + ".fasta.gz" into combined_fa_ch

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

process end2end_ltrs {
	publishDir "results", mode: 'copy'
	conda 'bioconductor-rtracklayer r-tidyverse'

	input:
	file fa from combined_fa_ch
	file gtf from combined_gtf_ch

	output:
	file "end2end-ltrs.*" into end2end_ch

	"""
	#!/usr/bin/env Rscript

	library(rtracklayer)
	library(tidyverse)

	fa <- import("$fa")

	gtf <- import("$gtf")

	internals <- fa[str_detect(names(fa),"[-_]I")]
	ltrs <- fa[str_detect(names(fa),"[-_]LTR")]
	others <- fa[!names(fa) %in% c(names(internals),names(ltrs))]

	# find tes that have an internal and ltr
	internals_names <- str_extract(string = names(internals), pattern = "^.+(?=[-_]I)")
	ltrs_names <- str_extract(string =  names(ltrs), pattern = "^.+(?=[-_]LTR)")
	mergeable <- intersect(internals_names, ltrs_names)
	names(internals_names) <-  names(internals)
	names(ltrs_names) <-  names(ltrs)

	# find Tes that only have 1 or the other
	unmergeable <- c(names(ltrs)[!ltrs_names %in% internals_names],
	  names(internals)[!internals_names %in% ltrs_names])

	# name the seqs so we can access them by the merged names
	names(internals) <- internals_names
	names(ltrs) <- ltrs_names

	res <- list()

	for (m in mergeable) {
  	in_seq <- internals[[m]]
  	ltr_seq <- ltrs[[m]]

  	mrg <- c(ltr_seq,in_seq,ltr_seq)

  	res[[m]] <- mrg
	}

	res <- c(Biostrings::DNAStringSet(res), fa[!names(fa) %in% mergeable])

	mergeable_gtf <- gtf[(str_extract(string = gtf\$gene_id, pattern = "^.+(?=[-_]I)") %in% mergeable) | (str_extract(string = gtf\$gene_id, pattern = "^.+(?=[-_]LTR)") %in% mergeable)]
	unmergeable_gtf <- gtf[!gtf\$gene_id %in% mergeable_gtf\$gene_id]

	gtf_df <- as_tibble(mergeable_gtf)

	gtf_df <- gtf_df %>%
	  mutate(seqnames=if_else(str_detect(seqnames,"[-_]I"),str_extract(string = seqnames, pattern = "^.+(?=[-_]I)"), as.character(seqnames))) %>%
	  mutate(seqnames=if_else(str_detect(seqnames,"[-_]LTR"),str_extract(string = seqnames, pattern = "^.+(?=[-_]LTR)"), as.character(seqnames)))


	res_gtf <- gtf_df %>%
	  mutate(start=if_else(str_detect(gene_id,"[-_]I") & seqnames %in% mergeable, start + width(ltrs[seqnames]), start)) %>%
	  mutate(end=if_else(str_detect(gene_id,"[-_]I") & seqnames %in% mergeable, end + width(ltrs[seqnames]), end)) %>%
	  dplyr::select(-width) %>%
	  GRanges() %>%
	  c(.,unmergeable_gtf)

	res_gtf <- sortSeqlevels(res_gtf)

	res_gtf <- sort(res_gtf)

	export(res_gtf,"end2end-ltrs.gtf")

	export(res,"end2end-ltrs.fasta.gz")

	"""

}
