rm(list = ls())

library(dada2)
library(ShortRead)
library(Biostrings)
library(phyloseq)
library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)
set.seed(123456)


setwd("set/to/work/dir")
unite.ref <- "path/to/unite/database/sh_general_release_dynamic_s_all_25.07.2023.fasta"  # CHANGE ME to location on your machine
path <- "path/to/all/raw/reads" # To raw reads
fnFs <- sort(list.files(path, pattern = "L001_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "L001_R2.fastq.gz", full.names = TRUE))

# Primers for removal
FWD <- "GTGARTCATCGAATCTTTG"
REV <- "TCCTCCGCTTATTGATATGC"


allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# Put N-filtered files in filtN/ subdirectory
fnFs.filtN <- file.path("path/to/filtering/dir/filtN", basename(fnFs)) 
fnRs.filtN <- file.path("path/to/filtering/dir/filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# Counts number of reads in which the primer is found
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
# SANITY CHECK: primer counts in reads
primer_count <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
                      FWD.ReverseReads = sapply(FWD.orients,primerHits, fn = fnRs.filtN[[1]]),
                      REV.ForwardReads = sapply(REV.orients, primerHits,fn = fnFs.filtN[[1]]),
                      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
print("SANITY CHECK: primer_count")
print(primer_count)

# Run version on cut adapt on your machine
cutadapt <- "path/to/cutadapt/bin/cutadapt" 
system2(cutadapt, args = "--version")


path.cut <- file.path("path/to/cutadapt/output/cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# SANITY CHECK: primer_removal
primer_removal <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
                        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
                        REV.ForwardReads = sapply(REV.orients, primerHits,fn = fnFs.cut[[1]]),
                        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
print("SANITY CHECK: primer_removal")
print(primer_removal)

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "L001_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "L001_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
#get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
get.sample.name <- function(fname) {
  parts <- strsplit(basename(fname), "_")[[1]]
  sample_id <- sub("^[0-9]+-", "", parts[1])  # Remove the leading numbers and hyphen
  paste(sample_id, parts[2], sep = "_")
}

sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Assign outputs of trimmed and filtered files
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

# Filter and trim
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(4, 4), truncQ = 2,
                     minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
# SANITY CHECK: filtered_reads
print("SANITY CHECK: filtered reads (out)")
print(head(out))

saveRDS(out, file = "save/an/RDS_backup/filtered_reads.rds")
out <- readRDS("load/RDS_backup/filtered_reads.rds")

########################
## Learn error rates  ##
########################

# Run from here if need to restart
path.cut <- file.path("path/to/cutadapt/output/cutadapt")
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

# If some have 0 reads after filtering put in direct path
filtFs <- sort(list.files("path/to/cutadapt/output/cutadapt/filtered", 
                          pattern = "L001_R1.fastq.gz", full.names = TRUE)) 
filtRs <- sort(list.files("path/to/cutadapt/output/cutadapt/filtered", 
                          pattern = "L001_R2.fastq.gz", full.names = TRUE))

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, minOverlap = 12)

# Construct a sequence table
seqtab <- makeSequenceTable(mergers)

# SANITY  CHECK: seqtab
print("SANITY CHECK: inspecting sequence table (seqtab)")
print(dim(seqtab))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
print("SANITY CHECK: inspecting sequence lengths after chimera removal (seqtab.nochim)")
print(table(nchar(getSequences(seqtab.nochim))))

# Track reads through pipeline (ultimate SANITY CHECK)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN),
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
print(head(track))
saveRDS(track, file = "save/to/RDS_backup/tracked_filtering.rds")

################
## Save RDS  ##
###############
saveRDS(seqtab.nochim, file = "save/to/RDS_backup/seqtab_nochim_wellhomes_all.rds")

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)
saveRDS(taxa, file = "save/to/RDS_backup/taxa_wellhomes_all.rds")

taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)

samples.out <- rownames(seqtab.nochim)
location <- samples.out %>% str_replace("^[0-9]+-", "") %>% str_extract("^[A-Z0-9]+(?:_[0-9])?")
sample_df <- data.frame(Location = location)
row.names(sample_df) <- samples.out

ps <- phyloseq(tax_table(taxa), sample_data(sample_df), otu_table(seqtab.nochim, taxa_are_rows = FALSE))
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
standard_plot <- plot_bar(ps.top20, x="Location", fill = "Genus")


######################################################
## Pull out ASVs so dataframes can be analysis in R ##
######################################################

## Code adapted from:
## https://astrobiomike.github.io/amplicon/dada2_workflow_ex#analysis-in-r

## giving sequences headers more manageable names (ASV_1, ASV_2...) ##
asv_seqs <- colnames(seqtab.nochim) 
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character") 
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
} 

## Make and write out a fasta file containing the final asv sequences
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "save/output/to/ASV_files/asv.fa")

## make a count table ##
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "save/output/to/ASV_files/ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# Split the data into ASV identifiers and sequences
ASV <- asv_fasta[seq(1, length(asv_fasta), by = 2)]
sequences <- asv_fasta[seq(2, length(asv_fasta), by = 2)]

# Create a data frame with ASV and sequences
ASV_to_amplicon.tib <- data.frame(ASV = ASV, Sequence = sequences) %>% as_tibble()

taxa.df <- as.data.frame(taxa)
taxa.df$Sequence <- row.names(taxa.df)
taxa.tib <- taxa.df %>% as_tibble()

ASVs_taxonomy.tib <- inner_join(taxa.df, ASV_to_amplicon.tib,  join_by("Sequence")) %>%
  as_tibble() %>% select(-Sequence)

ASVs_taxonomy.tib
write.table(ASVs_taxonomy.tib, "save/output/to/ASV_files/ASVs_taxonomy.tsv", sep="\t", quote = F, col.names = NA)
