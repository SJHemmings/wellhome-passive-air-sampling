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

##Set Working Dir##
set.seed(123456)
setwd("set/to/working/dir")

#####################################
## Clean Data into tibble format ##
##################################### 

METADATA <- 
  read.csv("path/to/metadata/Wellhomes_metadata.csv") %>% as.tibble() %>% 
  rename_all(tolower) %>% 
  filter(sequenced %in% c("Panel 1", "Panel 2", "Panel 3")) %>% 
  rename(sample_id = id)

## Read in ASV count table
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

# Find out how many reads are in dataset, decide what to filter out below #
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

############################################
## PCoA of indoor vs outdoor environments ##
############################################
## Put data into a matrix ##
rownames(FILT_ASV_COUNTS) <- FILT_ASV_COUNTS$sample_id 
FILT_ASV_COUNTS <- FILT_ASV_COUNTS[, -1] 
ENV.matrix <- as.matrix(FILT_ASV_COUNTS) 

#######################
##vegan PCoA analysis##
#######################
#Use vegdist to calculate a distance matrix using bray-curtis, this will produce a lower triangle distance matrix
ENV.dist <- avgdist(ENV.matrix, dmethod = "bray", sample = 15569) #Sample = the lowest number of counts in a sample
pcoa <- cmdscale(ENV.dist, k=2, eig = TRUE, add = TRUE) 

#Make PCoA tibble from data frame
POSITIONS <- pcoa$points 
colnames(POSITIONS) <- c("PCoA1", "PCoA2")

ENV_META_DATA <- 
  POSITIONS %>% as_tibble(rownames = "sample_id") %>% 
  inner_join(METADATA, by = "sample_id") %>%  
  select(sample_id, type, season, start_date, collection_date, long, lat, campaign) 

##adonis analysis##
test <- adonis2(ENV.dist~ENV_META_DATA$type, permutations = 1e4) 
p_value <- test$`Pr(>F)`[1]
R2 <- test$R2[1]
F_value <- test$F[1]

## Check Distribution ##
BETADIS <- betadisper(ENV.dist, ENV_META_DATA$type)
anova(BETADIS)
permutest(BETADIS)

ENV_PCOA.tibble <- 
  POSITIONS %>% as_tibble(rownames = "sample_id") %>% 
  inner_join(METADATA, by = "sample_id") %>%  
  select(sample_id, type, PCoA1, PCoA2, season, start_date, collection_date, long, lat, campaign)

#######################################
## ggplot for indoor vs outdoor PCoA ##
#######################################

percent_explained <- 100 * pcoa$eig / sum(pcoa$eig) 
round_percent <- format(round(percent_explained [1:2], digits = 1), nsmall = 1, trim = TRUE) 

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

##STATS READ OUT##
p_value 
R2 
F_value 
IN_OUT_PCOA

############################################
## PCoA for Indoor & Outdoor Seasonality  ##
############################################

## Filter for indoor ##
INDOOR_FILT_ASV_COUNTS <- 
  ASV_COUNTS %>% 
  pivot_longer(-"asv") %>%  
  rename(
    sample_id = name,
    reads = value) %>% 
  inner_join(., TAXONOMY) %>% 
  filter(kingdom == "Fungi") %>% 
  select(asv, sample_id, reads) %>% 
  inner_join(., METADATA) %>% filter(type != "Standard") %>%
  filter(type == "Indoor") %>% 
  select(asv, sample_id, reads) %>%  
  group_by(sample_id) %>% 
  mutate(total_reads = sum(reads)) %>%
  filter(total_reads > 6000) %>% 
  select(-total_reads) %>% ungroup(sample_id) %>%  
  pivot_wider(names_from = asv, values_from = reads) %>% 
  as.data.frame()

## Filter for outdoor ##
OUTDOOR_FILT_ASV_COUNTS <-
  ASV_COUNTS %>% 
  pivot_longer(-"asv") %>%  
  rename(
    sample_id = name,
    reads = value) %>% 
  inner_join(., TAXONOMY) %>% 
  filter(kingdom == "Fungi") %>% 
  select(asv, sample_id, reads) %>%
  inner_join(., METADATA) %>% filter(type != "Standard") %>% 
  filter(type == "Outdoor") %>% 
  select(asv, sample_id, reads) %>%  
  group_by(sample_id) %>% 
  mutate(total_reads = sum(reads)) %>%
  filter(total_reads > 6000) %>% 
  select(-total_reads) %>% ungroup(sample_id) %>%  
  pivot_wider(names_from = asv, values_from = reads) %>% 
  as.data.frame()

#################
## Indoor PCoA ##
#################
## Put data into a matrix ##
rownames(INDOOR_FILT_ASV_COUNTS) <- INDOOR_FILT_ASV_COUNTS$sample_id 
INDOOR_FILT_ASV_COUNTS <- INDOOR_FILT_ASV_COUNTS[, -1] 
INDOOR_ENV.matrix <- as.matrix(INDOOR_FILT_ASV_COUNTS) 

##vegan PCoA analysis##
#Use vegdist to calculate a distance matrix using bray-curtis, this will produce a lower triangle distance matrix
INDOOR_ENV.dist <- avgdist(INDOOR_ENV.matrix, dmethod = "bray", sample = 15569) #Sample = the lowest number of counts in a sample
pcoa <- cmdscale(INDOOR_ENV.dist, k=2, eig = TRUE, add = TRUE) 

#Make PCoA tibble from data frame
POSITIONS <- pcoa$points 
colnames(POSITIONS) <- c("PCoA1", "PCoA2") 

ENV_META_DATA <- 
  POSITIONS %>% as_tibble(rownames = "sample_id") %>% 
  inner_join(METADATA, by = "sample_id") %>%  
  select(sample_id, type, season, start_date, collection_date, long, lat, campaign) 

## Calculate indoor Statistics ##
## PERMANOVA ##
test <- adonis2(INDOOR_ENV.dist~ENV_META_DATA$season, permutations = 1e4) 
p_value <- test$`Pr(>F)`[1]
R2 <- test$R2[1]
F_value <- test$F[1]
test

##Pairwise PERMOANOVA##
INDOOR_META_DATA <- ENV_META_DATA %>% filter(type == "Indoor")
Indoor.pairwise.Perm.result <- pairwise.adonis(INDOOR_ENV.dist, INDOOR_META_DATA$season, sim.method = "bray", p.adjust.m = "bonferroni")
Indoor.pairwise.Perm.result

## PERMDISP ##
BETADIS <- betadisper(INDOOR_ENV.dist, ENV_META_DATA$season)
anova(BETADIS)
permutest(BETADIS)

ENV_PCOA.tibble <- 
  POSITIONS %>% as_tibble(rownames = "sample_id") %>% 
  inner_join(METADATA, by = "sample_id") %>%  
  select(sample_id, type, PCoA1, PCoA2, season, start_date, collection_date, long, lat, campaign)

###########################
## ggplot of Indoor PCoA ##
###########################
percent_explained <- 100 * pcoa$eig / sum(pcoa$eig) 
round_percent <- format(round(percent_explained [1:2], digits = 1), nsmall = 1, trim = TRUE) 

labs <- c(glue("PCoA 1 ({round_percent[1]}%)"),
          glue("PCoA 2 ({round_percent[2]}%)"))

# Ensure that the seasons are in the desired order (Winter, Spring, Summer, Autumn)
ENV_PCOA.tibble <- ENV_PCOA.tibble %>%
  mutate(season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")))

CENTROID <- 
  ENV_PCOA.tibble %>% 
  group_by(season) %>% 
  summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2))

INDOOR_SEASONALITY <-
  ENV_PCOA.tibble %>%
  rename(Environment = type) %>% 
  ggplot(aes(x = PCoA1, y = PCoA2, color = season)) +
  stat_ellipse(aes(fill = season), geom = "polygon", alpha = 0.2, show.legend = FALSE) + 
  geom_point(size = 2) +
  geom_point(data = CENTROID, size = 4, shape = 21, color = "black", aes(fill = season), show.legend = FALSE) +
  scale_color_manual(values = c("Winter" = "#5BBCD6", "Spring" = "#2BB187", "Summer" = "#F2AD00", "Autumn" = "#FF0000")) +
  scale_fill_manual(values = c("Winter" = "#5BBCD6", "Spring" = "#2BB187", "Summer" = "#F2AD00", "Autumn" = "#FF0000")) +
  labs(x = labs[1], y = labs[2]) +
  theme_classic() +
  theme(panel.background = element_rect(colour = "black", size = 0.1), legend.position = "none")

##############################
## Outdoor Seasonality PCoA ##
##############################
## Put data into a matrix ##
rownames(OUTDOOR_FILT_ASV_COUNTS) <- OUTDOOR_FILT_ASV_COUNTS$sample_id 
OUTDOOR_FILT_ASV_COUNTS <- OUTDOOR_FILT_ASV_COUNTS[, -1] 
OUTDOOR_ENV.matrix <- as.matrix(OUTDOOR_FILT_ASV_COUNTS) 

## vegan PCoA analysis ##
OUTDOOR_ENV.dist <- avgdist(OUTDOOR_ENV.matrix, dmethod = "bray", sample = 15569) #Sample = the lowest number of counts in a sample
pcoa <- cmdscale(OUTDOOR_ENV.dist, k=2, eig = TRUE, add = TRUE) 

# Make PCoA tibble from data frame
POSITIONS <- pcoa$points 
colnames(POSITIONS) <- c("PCoA1", "PCoA2")

ENV_META_DATA <- 
  POSITIONS %>% as_tibble(rownames = "sample_id") %>% 
  inner_join(METADATA, by = "sample_id") %>%  
  select(sample_id, type, season, start_date, collection_date, long, lat, campaign) 

## PERMANOVA analysis##
test <- adonis2(OUTDOOR_ENV.dist~ENV_META_DATA$season, permutations = 1e4) 
p_value <- test$`Pr(>F)`[1]
R2 <- test$R2[1]
F_value <- test$F[1]

## PERMDISP analysis ##
BETADIS <- betadisper(OUTDOOR_ENV.dist, ENV_META_DATA$season)
anova(BETADIS)
permutest(BETADIS)

ENV_PCOA.tibble <- 
  POSITIONS %>% as_tibble(rownames = "sample_id") %>% 
  inner_join(METADATA, by = "sample_id") %>%  
  select(sample_id, type, PCoA1, PCoA2, season, start_date, collection_date, long, lat, campaign)

###########################
##ggplot for Outdoor PCoA ##
############################

percent_explained <- 100 * pcoa$eig / sum(pcoa$eig) 
round_percent <- format(round(percent_explained [1:2], digits = 1), nsmall = 1, trim = TRUE) 

labs <- c(glue("PCoA 1 ({round_percent[1]}%)"),
          glue("PCoA 2 ({round_percent[2]}%)"))

# Ensure that the seasons are in the desired order (Winter, Spring, Summer, Autumn)
ENV_PCOA.tibble <- ENV_PCOA.tibble %>%
  mutate(season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")))

CENTROID <- 
  ENV_PCOA.tibble %>% 
  group_by(season) %>% 
  summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2))

OUTDOOR_SEASONALITY <-
  ENV_PCOA.tibble %>%
  rename(Environment = type) %>% 
  ggplot(aes(x = PCoA1, y = PCoA2, color = season)) +
  stat_ellipse(aes(fill = season), geom = "polygon", alpha = 0.2, show.legend = FALSE) + 
  geom_point(size = 2) +
  geom_point(data = CENTROID, size = 4, shape = 21, color = "black", aes(fill = season), show.legend = FALSE) +
  scale_color_manual(values = c("Winter" = "#5BBCD6", "Spring" = "#2BB187", "Summer" = "#F2AD00", "Autumn" = "#FF0000")) +
  scale_fill_manual(values = c("Winter" = "#5BBCD6", "Spring" = "#2BB187", "Summer" = "#F2AD00", "Autumn" = "#FF0000")) +
  labs(x = labs[1], y = labs[2]) +
  theme_classic() +
  theme(panel.background = element_rect(colour = "black", size = 0.1), legend.position = "none")

## Outdoor pairwise PERMANOVA ##
OUTDOOR_META_DATA <- ENV_META_DATA %>% filter(type == "Outdoor")
Outdoor.pairwise.Perm.result <- pairwise.adonis(OUTDOOR_ENV.dist, OUTDOOR_META_DATA$season, sim.method = "bray", p.adjust.m = "bonferroni")
Outdoor.pairwise.Perm.result
Indoor.pairwise.Perm.result