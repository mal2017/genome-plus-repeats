library(rtracklayer)
library(tidyverse)
#https://www.biostars.org/p/464924/

fix_names <- function(nms) {
  w0 <- gsub("[\t\n\r\v\f\a\b ]","",x=gsub("(?<=[\t\n\r\v\f\a\b ]).+","",x = nms, perl = T))
  w <- gsub(".+\\|","", x=w0)
  w1 <- gsub("#.+","",x = w)
  w1
}

te_fa <- import("repeats.fasta")

names(te_fa) <- names(te_fa) %>% fix_names()

# --------------------------

genome_and_tes_fa <- import("results/plus-repeats.repeatmasked.fasta.gz")

names(genome_and_tes_fa) <- names(genome_and_tes_fa) %>% fix_names()


stopifnot(all(names(te_fa) %in% names(genome_and_tes_fa)))


# -----

tes_gtf <- import("results/plus-repeats.gtf")

seqlevels(tes_gtf) <- seqlevels(tes_gtf) %>% fix_names()

tes_gtf$transcript_id <- tes_gtf$transcript_id %>% fix_names()
tes_gtf$gene_id <- tes_gtf$gene_id %>% fix_names()

stopifnot(all(names(te_fa) %in% seqnames(tes_gtf)))
stopifnot(all(seqnames(tes_gtf) %in% names(genome_and_tes_fa)))

# -----

masked_gff <- import("results/plus-repeats.repeatmasked.gff")

a0 <- masked_gff$Target %>% fix_names()
a2 <- str_remove(a0,"\"")
a3 <- str_remove(a2,"Motif:")

masked_gff$Target <- a3
masked_gff$name <- a3

# ----
# have to match the ucsc track format: https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema
rpm_cn <- c("chrom","start","end","te","smith.waterman","strand","subs.pct","del.pct","ins.pct","bases.past.match","class","bases.in.cons.complement","match.start","match.end","ins.id","has.higher.match")
masked_outf <- read_table(pipe("rmsk2bed < results/plus-repeats.repeatmasked.out"), col_names = rpm_cn)

masked_outf<-masked_outf %>%
  mutate(bin = 1) %>%
  relocate(bin) %>%
  select(bin,smith.waterman, subs.pct, del.pct, ins.pct, chrom, start, end, bases.in.cons.complement, strand, te, class, match.start, match.end, bases.past.match) %>%
  mutate(id = row_number()) %>%
  mutate(class2 = class) %>%
  dplyr::relocate(class2,.after=class) 

# sanity check for parsing issues.
#masked_outf %>% filter(!str_detect(te,regex("\\([ACTG]+\\)n"))) %>% pull(te) %>% unique()

masked_outf$te <- masked_outf$te %>% fix_names()

# ----

export(te_fa,"results/Tidalbase_transposon_sequence.filtered.fasta")
export(genome_and_tes_fa,"results/plus-repeats.repeatmasked.fixednames.fasta")
export(tes_gtf, "results/plus-repeats.fixednames.gtf")
export(masked_gff, "results/plus-repeats.repeatmasked.fixednames.bed")
write_tsv(masked_outf, "results/plus-repeats.repeatmasked.txt", col_names = F) # this is the ucsc track format

# these must be bgzipped afterwards and the Tidalbase seqs must be samtools faidx'd as well
