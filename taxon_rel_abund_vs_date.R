#Clear R
rm(list = ls())

##Load libraries##
library(tidyverse)
library(readxl)
library(phyloseq)
library(ggtext)
library(RColorBrewer)
library(ggsci)
library(wesanderson)

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
  select(-total_reads) %>% ungroup(sample_id) 

## Find out the relative abundance of the filtered asv counts
ASV_REL_ABUND <- 
  inner_join(FILT_ASV_COUNTS, METADATA, by="sample_id") %>% 
  select(asv, sample_id, reads, home, campaign, group, type, sequenced, season, collection_date) %>% 
  inner_join(TAXONOMY, by="asv") %>% 
  filter(kingdom == "Fungi") %>% 
  group_by(sample_id) %>%  
  mutate(rel_abund = reads / sum(reads)) %>% 
  ungroup(sample_id) %>%
  pivot_longer(cols = c("kingdom", "phylum", "class", "order", "family", "genus", "species", "asv"), 
               names_to = "level", 
               values_to = "taxon") %>%
  mutate(taxon = if_else(is.na(taxon), "Unknown Fungal Classification", taxon))

#########################
## Plot Phylum by date ##
#########################
## Filter to just look at Genus but can change this to any classification
PHY_REL_ABUND <- 
  ASV_REL_ABUND %>% 
  filter(level=="phylum") %>% 
  group_by(sample_id, taxon, type, campaign, season) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>% 
  mutate(taxon = str_replace(taxon, 
                             "^(\\S*)$", "*\\1*")) 

# Plotting all genus is abit much just plot over 2% rel_abund
PHY_POOL <- 
  PHY_REL_ABUND %>%
  group_by(taxon, type) %>% 
  summarise(mean = mean(rel_abund), .groups = "drop") %>%
  group_by(taxon) %>% 
  summarise(pool = max(mean) < 25, 
            mean = mean(mean), 
            .groups = "drop") %>% 
  arrange(desc(mean)) %>% print(n=30)

###################################
## ggplot phylum by data: Indoor ##
###################################

PHY_DATE <- 
  inner_join(PHY_REL_ABUND, PHY_POOL) %>%
  inner_join(METADATA) %>% 
  select(taxon, sample_id, type, campaign, rel_abund, collection_date, pool, mean) %>% 
  filter(type == "Indoor") %>% 
  mutate(collection_date = as.Date(collection_date, format = "%d/%m/%Y")) %>% 
  mutate(taxon = if_else(pool, "Other", taxon)) %>%  
  group_by(taxon, sample_id, collection_date) %>%
  summarise(rel_abund = sum(rel_abund), 
            mean = sum(mean),
            .groups = "drop") %>%  
  rename(Phylum = taxon) %>% 
  filter(Phylum %in% c("*Ascomycota*", "*Basidiomycota*"))

Indoor_phylum_vs_date <-
ggplot(PHY_DATE, aes(x = collection_date, y = rel_abund, color = Phylum)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_point(size = 2, alpha = 0.7, shape = 21, stroke = 0.5, color = "black") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +  
  scale_x_date(date_labels = "%b-%Y", date_breaks = "1 month") +  
  labs(x = "Collection Date", y = "Relative Abundance", color = "Phylum") +
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.position = "none") +  
  scale_color_manual(values = wes_palette("Darjeeling1"))

###################################
## ggplot phylum by data: Outdoor ##
###################################

PHY_DATE <- 
  inner_join(PHY_REL_ABUND, PHY_POOL) %>%
  inner_join(METADATA) %>% 
  select(taxon, sample_id, type, campaign, rel_abund, collection_date, pool, mean) %>% 
  filter(type == "Outdoor") %>% 
  mutate(collection_date = as.Date(collection_date, format = "%d/%m/%Y")) %>% 
  mutate(taxon = if_else(pool, "Other", taxon)) %>%  
  group_by(taxon, sample_id, collection_date) %>%
  summarise(rel_abund = sum(rel_abund), 
            mean = sum(mean),
            .groups = "drop") %>% 
  rename(Phylum = taxon) %>% 
  filter(Phylum %in% c("*Ascomycota*", "*Basidiomycota*"))

Outdoor_phylum_vs_date <-
ggplot(PHY_DATE, aes(x = collection_date, y = rel_abund, color = Phylum)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_point(size = 2, alpha = 0.7, shape = 21, stroke = 0.5, color = "black") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
  scale_x_date(date_labels = "%b-%Y", date_breaks = "1 month") + 
  labs(x = "Collection Date", y = "Relative Abundance", color = "Phylum") +
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.position = "none") +  
  scale_color_manual(values = wes_palette("Darjeeling1"))

GENUS_REL_ABUND <- 
  ASV_REL_ABUND %>% 
  filter(level=="genus") %>% 
  group_by(sample_id, taxon, type, campaign, season) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>% 
  mutate(taxon = str_replace(taxon, 
                             "^(\\S*)$", "*\\1*"))

###########################
## Indoor genera vs date ##
###########################

GENUS_POOL <- 
  GENUS_REL_ABUND %>%
  filter(type == "Indoor") %>% 
  group_by(taxon, type) %>%
  summarise(mean = mean(rel_abund), .groups = "drop") %>%
  group_by(taxon) %>% 
  summarise(pool = max(mean) < 1, 
            mean = mean(mean), 
            .groups = "drop") %>% 
  arrange(desc(mean))

GENUS_DATE <-
  inner_join(GENUS_REL_ABUND, GENUS_POOL) %>%
  inner_join(METADATA) %>%
  select(taxon, sample_id, type, campaign, rel_abund, collection_date, pool, mean) %>%
  filter(type == "Indoor") %>%
  mutate(collection_date = as.Date(collection_date, format = "%d/%m/%Y")) %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>% 
  group_by(taxon, sample_id, collection_date) %>%
  summarise(rel_abund = sum(rel_abund),
            mean = sum(mean),
            .groups = "drop") %>% 
  rename(Genus = taxon) %>%
  filter(!(Genus %in% c("Other", "Unknown Fungal Classification")))

Indoor_genera_vs_date <-
ggplot(GENUS_DATE, aes(x = collection_date, y = rel_abund, color = Genus)) +
  geom_point(size = 2, alpha = 0.3) + 
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, color = "black") +  
  geom_vline(xintercept = as.Date(c("2023-06-21", "2023-12-21", "2022-12-21", "2024-06-21")), linetype = "dashed", color = "black") +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "2 months") + 
  labs(x = "Collection Date", y = "Genus Relative Abundance (%)", color = "Genus") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.text = element_markdown(),
    strip.text = element_markdown(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    strip.background = element_rect(color = "black", fill = "lightgray")
  ) + 
  facet_wrap(~ Genus, scales = "free_y") +  
  scale_y_continuous(limits = c(0, NA))

###########################
## Outdoor genera vs date ##
###########################

GENUS_POOL <- 
  GENUS_REL_ABUND %>%
  filter(type == "Outdoor") %>% 
  group_by(taxon, type) %>% 
  summarise(mean = mean(rel_abund), .groups = "drop") %>%
  group_by(taxon) %>% 
  summarise(pool = max(mean) < 1, 
            mean = mean(mean), 
            .groups = "drop") %>% 
  arrange(desc(mean))

GENUS_DATE <-
  inner_join(GENUS_REL_ABUND, GENUS_POOL) %>%
  inner_join(METADATA) %>%
  select(taxon, sample_id, type, campaign, rel_abund, collection_date, pool, mean) %>%
  filter(type == "Outdoor") %>%
  mutate(collection_date = as.Date(collection_date, format = "%d/%m/%Y")) %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>% 
  group_by(taxon, sample_id, collection_date) %>%
  summarise(rel_abund = sum(rel_abund),
            mean = sum(mean),
            .groups = "drop") %>%  
  rename(Genus = taxon) %>%
  filter(!(Genus %in% c("Other", "Unknown Fungal Classification")))

Outdoor_genera_vs_date <-
  ggplot(GENUS_DATE, aes(x = collection_date, y = rel_abund, color = Genus)) +
  geom_point(size = 2, alpha = 0.3) +  
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, color = "black") + 
  geom_vline(xintercept = as.Date(c("2023-06-21", "2023-12-21", "2022-12-21", "2024-06-21")), linetype = "dashed", color = "black") +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "2 months") +
  labs(x = "Collection Date", y = "Genus Relative Abundance (%)", color = "Genus") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.text = element_markdown(),
    strip.text = element_markdown(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(color = "black", fill = "lightgray")
  ) + 
  facet_wrap(~ Genus, scales = "free_y") + 
  scale_y_continuous(limits = c(0, NA))

