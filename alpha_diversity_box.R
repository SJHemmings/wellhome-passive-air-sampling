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
library(patchwork)

set.seed(20220927)
setwd("path/to/working/dir")

###################################
## Clean Data into tibble format ##
###################################

METADATA <- 
  read.csv("path/to/metadata/Wellhomes_metadata.csv") %>% 
  as.tibble() %>% 
  rename_all(tolower) %>%  
  filter(sequenced  %in% c("Panel 1", "Panel 2", "Panel 3")) %>% 
  rename(sample_id = id)

ASV_COUNTS <- 
  read.csv("path/to/dada2/output/ASVs_counts.csv") %>% 
  as_tibble() 

TAXONOMY <- 
  read_tsv("path/to/dada2/output/ASVs_taxonomy.tsv") %>% 
  as.tibble() %>%
  select(-"...1") %>%
  mutate_at(vars(-ASV), ~ ifelse(!is.na(.), str_sub(., start = 4), .)) %>% 
  rename_all(tolower) %>% 
  mutate(asv = str_replace_all(asv, ">", "")) 

#How many reads are in each sample?
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
  arrange(total_reads) %>% print(n = 94)

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
  select(-total_reads) %>% ungroup(sample_id) 

#########################################################################
##Stats: Richness, Shannon, Simpson, inv. Simpson, total indv in sample##
#########################################################################

RICHNESS <- 
  function(x){
    sum(x>0)
  }

SHANNON <- 
  function(x){
    rabund <- x[x>0]/sum(x)
    -sum(rabund*log(rabund))
  }


SIMPSON <-
  function(x){
    n <- sum(x)
    sum(x*(x - 1)/(n*(n - 1)))
  }

###################
##Alpha div table##
####################
METADATA_ALPHA_DIVERSITY <- 
  FILT_ASV_COUNTS %>%
  group_by(sample_id) %>% 
  summarise(sobs = RICHNESS(reads),
            shannon = SHANNON(reads),
            simpson = SIMPSON(reads), 
            invsimpson = 1/simpson,
            n = sum(reads)) %>% 
  inner_join(METADATA, by = "sample_id") %>% 
  select(sample_id, sobs, shannon, simpson, invsimpson, n, type) %>% arrange(sobs)


#########
##Plots##
#########
## Number of Reads After DADA2 ##
READS_PLOT <-
  METADATA_ALPHA_DIVERSITY %>% 
  select(n, type) %>% 
  ggplot(aes(x=type, y=n, fill=type)) +
  scale_fill_d3(alpha = 0.7) +
  geom_boxplot(show.legend = FALSE, width=0.5, outlier.shape = NA, alpha=0.65,
               coef=1.5)+
  geom_jitter(show.legend = FALSE, width=0.25, shape=21, color="black", size=2) +
  coord_cartesian(ylim = c(50000, 200000)) +
  
  labs(x="Environment", y="Number of Reads") +
  theme(
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16), 
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18)  
  ) +
  theme(axis.text.x = element_markdown()) +
  scale_fill_manual(values = c("Indoor" = "#9b72cb", "Outdoor" = "#7BB662"))

## sobs ##
SOBS <-
  METADATA_ALPHA_DIVERSITY %>%
  select(sobs, type) %>% 
  ggplot(aes(x=type, y=sobs, fill=type)) +
  scale_fill_d3(alpha = 0.7) +
  geom_boxplot(show.legend = FALSE, width=0.5, outlier.shape = NA, alpha=0.65,
               coef=1.5)+
  geom_jitter(show.legend = FALSE, width=0.25, shape=21, color="black", size=2) +
  coord_cartesian(ylim = c(0, 1500)) +
  
  labs(x="Environment", y="Number of Observed ASVs") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 16),  
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18)  
  ) +
  scale_fill_manual(values = c("Indoor" = "#9b72cb", "Outdoor" = "#7BB662"))

## Shannon ##
SHANNON_PLOT <-
  METADATA_ALPHA_DIVERSITY %>%
  select(shannon, type) %>% 
  ggplot(aes(x=type, y=shannon, fill=type)) +
  scale_fill_d3(alpha = 0.7) +
  geom_boxplot(show.legend = FALSE, width=0.5, outlier.shape = NA, alpha=0.65,
               coef=1.5)+
  geom_jitter(show.legend = FALSE, width=0.25, shape=21, color="black", size=2) +
  coord_cartesian(ylim = c(0, 8)) +
  
  labs(x="Environment", y="Shannon Diversity") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16), 
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18)  
  ) +
  scale_fill_manual(values = c("Indoor" = "#9b72cb", "Outdoor" = "#7BB662"))

## inv-simpson ##
INV_SIMPSON_PLOT <- 
  METADATA_ALPHA_DIVERSITY %>% 
  select(invsimpson, type) %>% 
  ggplot(aes(x=type, y=invsimpson, fill=type)) +
  scale_fill_d3(alpha = 0.7) +
  geom_boxplot(show.legend = FALSE, width=0.5, outlier.shape = NA, alpha=0.65,
               coef=1.5)+
  geom_jitter(show.legend = FALSE, width=0.25, shape=21, color="black", size=2) +
  coord_cartesian(ylim = c(0, 125)) +
  
  labs(x="Environment", y="Inverse Simpson Diversity") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 16),  
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18)  
  ) +
  scale_fill_manual(values = c("Indoor" = "#9b72cb", "Outdoor" = "#7BB662"))

## simpson ##
SIMPSON_PLOT <- 
  METADATA_ALPHA_DIVERSITY %>% 
  select(simpson, type) %>% 
  ggplot(aes(x=type, y=simpson, fill=type)) +
  scale_fill_d3(alpha = 0.7) +
  geom_boxplot(show.legend = FALSE, width=0.5, outlier.shape = NA, alpha=0.65,
               coef=1.5)+
  geom_jitter(show.legend = FALSE, width=0.25, shape=21, color="black", size=2) +
  coord_cartesian(ylim = c(0, 0.4)) +
  
  labs(x="Environment", y="Simpson Diversity") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 16),  
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18)  
  ) +
  scale_fill_manual(values = c("Indoor" = "#9b72cb", "Outdoor" = "#7BB662"))


READS_SIG <- 
  METADATA_ALPHA_DIVERSITY %>% 
  t.test(n ~ type, data = .)
READS_SIG$p.value

SOBS_SIG <- 
  METADATA_ALPHA_DIVERSITY %>% 
  t.test(sobs ~ type, data = .)
SOBS_SIG$p.value

SHANNON_SIG <- 
  METADATA_ALPHA_DIVERSITY %>% 
  t.test(shannon ~ type, data = .)
SHANNON_SIG$p.value

INV_SIMPSON_SIG <- 
  METADATA_ALPHA_DIVERSITY %>% 
  t.test(invsimpson ~ type, data = .)
INV_SIMPSON_SIG$p.value

SIMPSON_SIG <- 
  METADATA_ALPHA_DIVERSITY %>% 
  t.test(simpson ~ type, data = .)
SIMPSON_SIG$p.value


(SOBS | SHANNON_PLOT | INV_SIMPSON_PLOT)
