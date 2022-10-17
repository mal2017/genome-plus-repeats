library(rtracklayer)
library(tidyverse)
#https://www.biostars.org/p/464924/

te_fa <- import("repeats.fasta")

names(te_fa) <- names(te_fa) %>% str_replace_all("\\\\","_")

x0 <- gsub("[\t\n\r\v\f\a\b ]","",x=gsub("(?<=[\t\n\r\v\f\a\b ]).+","",x = names(te_fa), perl = T))
x <- gsub(".+\\|","", x=x0)
x1 <- gsub("#.+","",x = x)

names(te_fa) <- gsub("gb\\|.+\\|","",x = x1)

# --------------------------

genome_and_tes_fa <- import("results/plus-repeats.repeatmasked.fasta.gz")

names(genome_and_tes_fa) <- names(genome_and_tes_fa) %>% str_replace_all("\\\\","_")

z0 <- gsub("[\t\n\r\v\f\a\b ]","",x=gsub("(?<=[\t\n\r\v\f\a\b ]).+","",x = names(genome_and_tes_fa), perl = T))
z <- gsub(".+\\|","", x=z0)
z1 <- gsub("#.+","",x = z)

names(genome_and_tes_fa) <- gsub("gb\\|.+\\|","",x = z1)


stopifnot(all(names(te_fa) %in% names(genome_and_tes_fa)))


# -----

tes_gtf <- import("results/plus-repeats.gtf")

seqlevels(tes_gtf) <- seqlevels(tes_gtf) %>% str_replace_all("\\\\","_")

a0 <- gsub("[\t\n\r\v\f\a\b ]","",x=gsub("(?<=[\t\n\r\v\f\a\b ]).+","",x = seqlevels(tes_gtf), perl = T))
a <- gsub(".+\\|","", x=a0)
a1 <- gsub("#.+","",x = a)

seqlevels(tes_gtf) <- gsub("gb\\|.+\\|","",x = a1)

tes_gtf$transcript_id <- seqnames(tes_gtf)
tes_gtf$gene_id <- seqnames(tes_gtf)

stopifnot(all(names(te_fa) %in% seqnames(tes_gtf)))
stopifnot(all(seqnames(tes_gtf) %in% names(genome_and_tes_fa)))

# -----

masked_gff <- import("results/plus-repeats.repeatmasked.gff")

a0 <- gsub("[\t\n\r\v\f\a\b ]","",x=gsub("(?<=[\t\n\r\v\f\a\b ]).+","",x = masked_gff$Target, perl = T))
a <- gsub(".+\\|","", x=a0)
a1 <- gsub("#.+","",x = a)
a2 <- str_remove(a1,"\"")
a3 <- str_remove(a2,"Motif:")

masked_gff$Target <- a3
masked_gff$name <- a3

# ----

export(te_fa,"results/Tidalbase_transposon_sequence.filtered.fasta")
export(genome_and_tes_fa,"results/plus-repeats.repeatmasked.fixednames.fasta")
export(tes_gtf, "results/plus-repeats.fixednames.gtf")
export(masked_gff, "results/plus-repeats.repeatmasked.fixednames.bed")


# these must be bgzipped afterwards and the Tidalbase seqs must be samtools faidx'd as well
