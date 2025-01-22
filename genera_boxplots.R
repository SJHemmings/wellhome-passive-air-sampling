# Clear R
rm(list = ls())

## Load libraries
library(tidyverse)
library(readxl)
library(phyloseq)
library(ggtext)
library(RColorBrewer)
library(ggsci)

## Set Working Dir
setwd("/set/to/working/dir") # Set working directory as project folder

#####################################
## 1. Clean Data into tibble format ##
#####################################

# Read in metadata
METADATA <- 
  read.csv("path/to/Wellhomes_metadata.csv") %>% as_tibble() %>% 
  rename_all(tolower) %>% 
  filter(sequenced %in% c("Panel 1", "Panel 2", "Panel 3")) %>% 
  rename(sample_id = id)

# Read in ASV count table
ASV_COUNTS <- 
  read.csv("path/to/dada2/output/ASVs_counts.csv") %>% 
  as_tibble()

# Read in Taxonomy table
TAXONOMY <- 
  read_tsv("path/to/dada2/output/ASVs_taxonomy.tsv") %>% as.tibble() %>% select(-"...1") %>% 
  mutate_at(vars(-ASV), ~ ifelse(!is.na(.), str_sub(., start = 4), .)) %>% rename_all(tolower) %>% 
  mutate(asv = str_replace_all(asv, ">", ""))

# Pivot ASV counts into long format
ASV_COUNTS_PIVOT <- 
  ASV_COUNTS %>% 
  pivot_longer(-asv, names_to = "sample_id", values_to = "count") %>% 
  rename_all(tolower)

# Filter samples
FILT_ASV_COUNTS <- 
  ASV_COUNTS %>% 
  pivot_longer(-"asv") %>%  
  rename(sample_id = name, reads = value) %>% 
  inner_join(., TAXONOMY) %>% 
  filter(kingdom == "Fungi") %>% 
  select(asv, sample_id, reads) %>% 
  inner_join(., METADATA) %>% filter(type != "Standard") %>% 
  select(asv, sample_id, reads) %>%  
  group_by(sample_id) %>% 
  mutate(total_reads = sum(reads)) %>% 
  filter(total_reads > 6000) %>% 
  select(-total_reads) %>% ungroup(sample_id)

# Calculate relative abundance
ASV_REL_ABUND <- 
  inner_join(FILT_ASV_COUNTS, METADATA, by = "sample_id") %>% 
  select(asv, sample_id, reads, home, campaign, group, type, sequenced, season) %>% 
  inner_join(TAXONOMY, by = "asv") %>% 
  filter(kingdom == "Fungi") %>% 
  group_by(sample_id) %>%  
  mutate(rel_abund = reads / sum(reads)) %>% 
  ungroup(sample_id) %>% 
  pivot_longer(cols = c("kingdom", "phylum", "class", "order", "family", "genus", "species", "asv"), 
               names_to = "level", 
               values_to = "taxon") %>% 
  mutate(taxon = if_else(is.na(taxon), "Unknown Fungal Classification", taxon))

######################################
## Further cleaning for plotting ##
######################################

# Select Genus Level
GENUS_REL_ABUND <- 
  ASV_REL_ABUND %>% 
  filter(level == "genus") %>% 
  group_by(sample_id, taxon, type, campaign, season) %>% 
  summarise(rel_abund = 100 * sum(rel_abund), .groups = "drop") %>% 
  mutate(taxon = str_replace(taxon, "^(\\S*)$", "*\\1*"))

# Pool genera with <1% mean abundance
GENUS_POOL <- 
  GENUS_REL_ABUND %>% 
  group_by(taxon, type) %>% 
  summarise(mean = mean(rel_abund), .groups = "drop") %>% 
  group_by(taxon) %>% 
  summarise(pool = max(mean) < 1, 
            mean = mean(mean), 
            .groups = "drop")

# Identify significant genera
TEST_SIG <- 
  GENUS_POOL %>% filter(pool == FALSE) %>% 
  select(taxon) %>% unique()

#####################
## ggplot Box Plot ##
#####################

ALL_PLOT <- 
  inner_join(GENUS_POOL, GENUS_REL_ABUND) %>% 
  mutate(taxon = if_else(pool, "Other", taxon)) %>% 
  group_by(type, taxon, sample_id) %>% 
  summarise(rel_abund = sum(rel_abund), 
            mean = sum(mean), 
            .groups = "drop") %>%  
  mutate(taxon = factor(taxon), 
         taxon = fct_reorder(taxon, mean, .desc = TRUE)) %>% 
  rename(Genus = taxon) %>% 
  rename(Environment = type) %>% 
  ggplot(aes(fill = Environment, x = rel_abund, y = Genus)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.size = 0, outlier.shape = 1) +
  scale_fill_d3(alpha = 0.7) +
  scale_x_continuous(expand = c(0.02, 0), limits = c(0, 60)) +
  labs(y = "Fungal Genus", x = "Relative Abundance (%)") +
  theme_classic() +   
  theme(axis.text.y = element_markdown(),
        legend.position = c(0.97, 0.57)) +
  scale_fill_manual(values = c("Indoor" = "#9b72cb", "Outdoor" = "#7BB662"))

#####################
## Significance Testing ##
#####################

p_values <- 
  TEST_SIG %>% 
  pull(taxon) %>% 
  map_df(~ {
    taxon <- .
    test_result <- GENUS_REL_ABUND %>% 
      filter(taxon == !!taxon) %>% 
      t.test(rel_abund ~ type, data = .)
    tibble(
      Genus = taxon,
      p_value = test_result$p.value
    )
  })

print(p_values) %>% print(n=28)
