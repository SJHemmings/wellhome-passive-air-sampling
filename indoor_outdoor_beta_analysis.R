#Clear R
rm(list = ls())

##Load libraries##
library(tidyverse)
library(readxl)
library(phyloseq)
library(ggtext)
library(RColorBrewer)
library(ggsci)
library(vegan)
library(glue)
library(devtools)
library(pairwiseAdonis)

set.seed(19980106) 
setwd("set/as/work/dir") 

#####################################
##1. Clean Data into tibble format ##
##################################### 

## Read in metadata set as tibble
METADATA <- 
  read.csv("path/to/Wellhomes_metadata.csv") %>% as.tibble() %>% 
  rename_all(tolower) %>% 
  filter(sequenced %in% c("Panel 1", "Panel 2", "Panel 3")) %>% 
  rename(sample_id = id)

## Read in the ASV count table as a tibble 
ASV_COUNTS <- 
  read.csv("path/to/dada2/output/ASVs_counts.csv") %>% 
  as_tibble()

## Read in Taxonomy table
TAXONOMY <- 
  read_tsv("path/to/dada2/output/ASVs_taxonomy.tsv") %>% as.tibble() %>% select(-"...1") %>% 
  mutate_at(vars(-ASV), ~ ifelse(!is.na(.), str_sub(., start = 4), .)) %>% rename_all(tolower) %>% 
  mutate(asv = str_replace_all(asv, ">", ""))

## Pivot ASV counts into long format
ASV_COUNTS_PIVOT <- 
  ASV_COUNTS %>% 
  pivot_longer(-asv, names_to = "sample_id", values_to = "count") %>%
  rename_all(tolower)  

# Find out how many reads are in dataset, to decide what threshold to set #
ASV_COUNTS %>%
  pivot_longer(-"asv") %>% 
  rename(sample_id = name,
         reads = value) %>% 
  inner_join(., TAXONOMY) %>% 
  filter(kingdom == "Fungi") %>% 
  select(asv, sample_id, reads) %>% 
  group_by(sample_id) %>% 
  mutate(total_reads = sum(reads)) %>% 
  group_by(sample_id, total_reads) %>% 
  summarise() %>% 
  arrange(total_reads) 

## Remove samples below 6,000 reads, therefore losing samples:
## NEG_1, NEG_2, NEG_3, WH010_1, WH046_2, A014 and drop standards
FILT_ASV_COUNTS <- 
  ASV_COUNTS %>% 
  pivot_longer(-"asv") %>%  
  rename(
    sample_id = name,
    reads = value) %>% 
  inner_join(., TAXONOMY) %>% 
  filter(kingdom == "Fungi") %>% 
  select(asv, sample_id, reads) %>%
  inner_join(., METADATA) %>% filter(type != "Standard") %>%
  select(asv, sample_id, reads) %>%  
  group_by(sample_id) %>% 
  mutate(total_reads = sum(reads)) %>% 
  filter(total_reads > 6000) %>% 
  select(-total_reads) %>% ungroup(sample_id) %>%  
  pivot_wider(names_from = asv, values_from = reads) %>% 
  as.data.frame()

##########
## PCoA ##
##########

## Put data into a matrix 
rownames(FILT_ASV_COUNTS) <- FILT_ASV_COUNTS$sample_id
FILT_ASV_COUNTS <- FILT_ASV_COUNTS[, -1] 
ENV.matrix <- as.matrix(FILT_ASV_COUNTS) 

#########################
## vegan PCoA analysis ##
#########################

## Use vegdist to calculate a distance matrix using bray-curtis, this will produce a lower triangle distance matrix
## Make PCoA from distance matrix, and eig dist
ENV.dist <- avgdist(ENV.matrix, dmethod = "bray", sample = 15569) #Sample = the lowest number of counts in a sample
pcoa <- cmdscale(ENV.dist, k=2, eig = TRUE, add = TRUE) 

## Make PCoA tibble from data frame
POSITIONS <- pcoa$points 
colnames(POSITIONS) <- c("PCoA1", "PCoA2") 

## Load and clean metadata
ENV_META_DATA <- 
  POSITIONS %>% as_tibble(rownames = "sample_id") %>% 
  inner_join(METADATA, by = "sample_id") %>%  
  select(sample_id, type, season, start_date, collection_date, long, lat, campaign) #Make metadata table only containing sample location and Environment

## adonis analysis
test <- adonis2(ENV.dist~ENV_META_DATA$type, permutations = 1e4) 
p_value <- test$`Pr(>F)`[1]
R2 <- test$R2[1]
F_value <- test$F[1]

## Check Distribution 
BETADIS <- betadisper(ENV.dist, ENV_META_DATA$type)
anova(BETADIS)
permutest(BETADIS)

ENV_PCOA.tibble <- 
  POSITIONS %>% as_tibble(rownames = "sample_id") %>% 
  inner_join(METADATA, by = "sample_id") %>%  
  select(sample_id, type, PCoA1, PCoA2, season, start_date, collection_date, long, lat, campaign)

##########
##ggplot##
##########

percent_explained <- 100 * pcoa$eig / sum(pcoa$eig) 
round_percent <- format(round(percent_explained [1:2], digits = 1), nsmall = 1, trim = TRUE) #percent for x & y axis

labs <- c(glue("PCo 1 ({round_percent[1]}%)"),
          glue("PCo 2 ({round_percent[2]}%)"))

CENTROID <- 
  ENV_PCOA.tibble %>% 
  group_by(type) %>% summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2)) %>% 
  rename(Environment = type)

IN_OUT_PCOA <-
  ENV_PCOA.tibble %>%
  rename(Environment = type) %>% 
  ggplot(aes(x = PCoA1, y = PCoA2, color = Environment)) +
  geom_point(size = 2) +
  geom_point(data = CENTROID, size = 4, shape = 21, color = "black", aes(fill = Environment), show.legend = FALSE) +
  scale_color_manual(values = c("#9b72cb", "#7BB662")) +
  scale_fill_manual(values = c("#9b72cb", "#7BB662")) +
  labs(x = labs[1], y = labs[2]) +
  stat_ellipse(aes(fill = Environment), geom = "polygon", alpha = 0.25, show.legend = FALSE) +  # Using stat_ellipse with geom="polygon"
  theme_classic() +
  theme(panel.background = element_rect(colour = "black", size = 0.1)) + 
  scale_fill_manual(values = c("Indoor" = "#9b72cb", "Outdoor" = "#7BB662"))

##################
##STATS READ OUT##
##################
p_value 
R2 
F_value 
IN_OUT_PCOA
