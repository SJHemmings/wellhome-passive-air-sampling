#Clear R
rm(list = ls())

##Load libraries##
library(tidyverse)
library(readxl)
library(phyloseq)
library(ggtext)
library(RColorBrewer)
library(ggsci)

##Set Working Dir##
setwd("path/to/work/dir")

#####################################
##1. Clean Data into tibble format ##
##################################### 
METADATA <- 
  read.csv("path/to/metadata/Wellhomes_metadata.csv") %>% as.tibble() %>% 
  rename_all(tolower) %>% 
  filter(sequenced %in% c("Panel 1", "Panel 2", "Panel 3")) %>% 
  rename(sample_id = id)

ASV_COUNTS <- 
  read.csv("path/to/dada2/output/ASVs_counts.csv") %>% 
  as_tibble()

TAXONOMY <- 
  read_tsv("path/to/dada2/output/ASVs_taxonomy.tsv") %>% as.tibble() %>% select(-"...1") %>% 
  mutate_at(vars(-ASV), ~ ifelse(!is.na(.), str_sub(., start = 4), .)) %>% rename_all(tolower) %>% 
  mutate(asv = str_replace_all(asv, ">", ""))

ASV_COUNTS_PIVOT <- 
  ASV_COUNTS %>% 
  pivot_longer(-asv, names_to = "sample_id", values_to = "count") %>%
  rename_all(tolower)  

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
  select(-total_reads) %>% ungroup(sample_id)

## Find out the relative abundance of the filtered asv counts
ASV_REL_ABUND <- 
  inner_join(FILT_ASV_COUNTS, METADATA, by="sample_id") %>%
  select(asv, sample_id, reads, home, campaign, group, type, sequenced, season) %>% 
  inner_join(TAXONOMY, by="asv") %>% 
  filter(kingdom == "Fungi") %>% 
  group_by(sample_id) %>%  
  mutate(rel_abund = reads / sum(reads)) %>% 
  ungroup(sample_id) %>%
  pivot_longer(cols = c("kingdom", "phylum", "class", "order", "family", "genus", "species", "asv"), 
               names_to = "level", 
               values_to = "taxon") %>% 
  mutate(taxon = if_else(is.na(taxon), "Unknown Fungal Classification", taxon))

########################################
## Plot relative abundance as Stacked ##
########################################
## Filter to just look at Genus but can change this to any classification
GENUS_REL_ABUND <- 
  ASV_REL_ABUND %>% 
  filter(level=="genus") %>%
  group_by(sample_id, taxon, type, campaign, season) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>% 
  mutate(taxon = str_replace(taxon, 
                             "^(\\S*)$", "*\\1*"))

# Plotting all genus is abit much just plot over 2% rel_abund
GENUS_POOL <- 
  GENUS_REL_ABUND %>%
  group_by(taxon, type) %>% 
  summarise(mean = mean(rel_abund), .groups = "drop") %>%
  group_by(taxon) %>% 
  summarise(pool = max(mean) < 2, 
            mean = mean(mean), 
            .groups = "drop") %>% 
  arrange(desc(mean)) %>% print(n=30)

## Make stacked bar charts for camp1, camp2 and outdoors
camp1 <-
  inner_join(GENUS_POOL, GENUS_REL_ABUND) %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%  
  group_by(type, taxon, sample_id, campaign, season) %>%
  summarise(rel_abund = sum(rel_abund), 
            mean = sum(mean),
            .groups = "drop") %>% 
  rename(Genus = taxon) %>%  
  mutate(campaign = factor(campaign, levels = c("Campaign 1", "Priority 1", "Campaign 2", "Priority 2", "Outdoor"))) %>%  
  filter(campaign %in% c("Campaign 1", "Priority 1")) %>%
  ggplot(aes(x = sample_id, y = rel_abund, fill = Genus)) + 
  geom_col(width = 0.92) + 
  facet_grid(~campaign, scales = "free_x", space = "free", switch = "x") + 
  scale_fill_d3(palette = "category20") +
  labs(x = NULL, y = "Relative Abundance (%)") +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 1, size = 14),
    legend.text = element_markdown(),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing.x = unit(0.5, "cm"),
    legend.position = "bottom", 
    legend.box = "horizontal",  
    legend.box.just = "center" 
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

camp2 <-
  inner_join(GENUS_POOL, GENUS_REL_ABUND) %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%  
  group_by(type, taxon, sample_id, campaign, season) %>%
  summarise(rel_abund = sum(rel_abund), 
            mean = sum(mean),
            .groups = "drop") %>% 
  rename(Genus = taxon) %>% 
  mutate(campaign = factor(campaign, levels = c("Campaign 1", "Priority 1", "Campaign 2", "Priority 2", "Outdoor"))) %>% 
  filter(campaign %in% c("Campaign 2", "Priority 2")) %>%
  ggplot(aes(x = sample_id, y = rel_abund, fill = Genus)) +
  geom_col(width = 0.92) + 
  facet_grid(~campaign, scales = "free_x", space = "free", switch = "x") + 
  scale_fill_d3(palette = "category20") +
  labs(x = NULL, y = "Relative Abundance (%)") +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 1, size = 14),
    legend.text = element_markdown(),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing.x = unit(0.5, "cm"),
    legend.position = "bottom",  
    legend.box = "horizontal", 
    legend.box.just = "center"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

outdoor <-
  inner_join(GENUS_POOL, GENUS_REL_ABUND) %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>% 
  group_by(type, taxon, sample_id, campaign, season) %>%
  summarise(rel_abund = sum(rel_abund), 
            mean = sum(mean),
            .groups = "drop") %>%
  rename(Genus = taxon) %>%
  mutate(campaign = factor(campaign, levels = c("Campaign 1", "Priority 1", "Campaign 2", "Priority 2", "Outdoor"))) %>%
  filter(campaign %in% c("Outdoor")) %>%
  ggplot(aes(x = sample_id, y = rel_abund, fill = Genus)) +
  geom_col(width = 0.92) + 
  facet_grid(~campaign, scales = "free_x", space = "free", switch = "x") +
  scale_fill_d3(palette = "category20") +
  labs(x = NULL, y = "Relative Abundance (%)") +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 1, size = 14),
    legend.text = element_markdown(),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing.x = unit(0.5, "cm"),
    legend.position = "bottom", 
    legend.box = "horizontal",  
    legend.box.just = "center"  
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
