#Clear R
rm(list = ls())


##Load libraries##
library(tidyverse)
library(readxl)
library(phyloseq)
library(ggtext)
library(RColorBrewer)
library(ggsci)
library(patchwork)

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
  read_tsv("path/to/dada2/output/ASVs_taxonomy.tsv") %>% 
  as.tibble() %>% select(-"...1") %>% 
  mutate_at(vars(-ASV), ~ ifelse(!is.na(.), str_sub(., start = 4), .)) %>% 
  rename_all(tolower) %>% 
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
  select(asv, sample_id, reads, home, campaign, group, type, sequenced, season) %>% 
  inner_join(TAXONOMY, by="asv") %>%
  filter(kingdom == "Fungi") %>%
  group_by(sample_id) %>%  
  mutate(rel_abund = reads / sum(reads)) %>% 
  ungroup(sample_id) %>%
  pivot_longer(cols = c("kingdom", "phylum", "class", "order", "family", "genus", "species", "asv"), 
               names_to = "level", 
               values_to = "taxon") %>% 
  mutate(taxon = if_else(is.na(taxon), "Unknown_Fungal_Classification", taxon))

## Organise dataset so tsv will run correctly in FUNGuilds package
wellhome_input_guilds <- 
  inner_join(FILT_ASV_COUNTS, METADATA, by="sample_id") %>% 
  inner_join(TAXONOMY, by="asv") %>% 
  filter(kingdom == "Fungi") %>% 
  group_by(sample_id) %>% 
  mutate(species = if_else(is.na(species), "Unknown_Fungal_Classification", species),
         species = if_else(str_detect(species, "Incertae_sedis$"), "Unknown_Fungal_Classification", species)) %>% 
  mutate(species = paste0("s__", genus, "_", species)) %>% 
  pivot_longer(cols = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), 
               names_to = "level", 
               values_to = "taxon") %>% 
  filter(level == "species") %>% 
  select(asv, sample_id, reads, taxon) %>%
  rename(OTU = asv,
         taxonomy = taxon) %>% 
  pivot_wider(., names_from = sample_id, values_from = reads)
write_tsv(wellhome_input_guilds, "path/to/work/dir/wellhome_input_guilds.tsv")

## Stop at this part of the pipeline and save wellhome_input_guilds as a .tsv file ##
## This can be run through the package FUNGuilds: https://github.com/UMNFuN/FUNGuild
## save output as wellhome_input_guilds.guilds.txt

guilds <- read_tsv("Metadata/wellhome_input_guilds.guilds.txt")
guilds %>% colnames()
guilds.edit <- 
  guilds %>% select(-"Citation/Source", -"Notes", -"Confidence Ranking", -"Trait", 
                    -"Growth Morphology", -"Trophic Mode", -"Taxon Level", -"Taxon")
write_tsv(guilds.edit, "path/to/work/dir/wellhome_guild_final_edit.tsv") 
## edit file in bash so that the guild names are seperated by tabs and columns are titled Guild1,2,3,4...7 
## save as wellhome_guild_final_edit.tsv

guilds <- read_tsv("path/to/work/dir/wellhome_guild_final_edit.tsv")
guild_counts <- 
  guilds %>% select("OTU", "taxonomy", "Guild1", "Guild2", "Guild3", "Guild4", "Guild5", "Guild6", "Guild7") %>% 
  pivot_longer(cols = c("Guild1", "Guild2", "Guild3", "Guild4", "Guild5", "Guild6", "Guild7"),
               names_to = "guild_level",
               values_to = "Guild") %>% 
  drop_na(Guild) %>% 
  distinct() %>%
  rename(asv = OTU) %>% 
  inner_join(., FILT_ASV_COUNTS) %>% 
  rename(count = reads)

Indoor_guild_percentages <- 
  guild_counts %>% 
  inner_join(., METADATA) %>% 
  filter(type == "Indoor") %>% 
  select("asv", "taxonomy", "guild_level", "Guild", "sample_id", "count") %>% 
  group_by(sample_id, Guild) %>% 
  summarise(guild_count = sum(count)) %>% 
  ungroup() %>% 
  group_by(sample_id) %>% 
  mutate(sample_avg = guild_count / sum(guild_count)) %>% 
  ungroup() %>% 
  select(Guild, sample_avg) %>%  
  group_by(Guild) %>% 
  summarise(guild_avg = sum(sample_avg)) %>% 
  mutate(percentage = 100 * guild_avg / sum(guild_avg)) %>% 
  mutate(Guild = if_else(percentage < 2, "Other Guild >2%", Guild)) %>% 
  group_by(Guild) %>%                                                 
  summarise(percentage = sum(percentage)) %>%                         
  mutate(Guild = fct_reorder(Guild, percentage, .desc = FALSE)) %>%   
  arrange(desc(percentage))

Outdoor_guild_percentages <- 
  guild_counts %>% 
  inner_join(., METADATA) %>% 
  filter(type == "Outdoor") %>% 
  select("asv", "taxonomy", "guild_level", "Guild", "sample_id", "count") %>% 
  group_by(sample_id, Guild) %>% 
  summarise(guild_count = sum(count)) %>% 
  ungroup() %>% 
  group_by(sample_id) %>% 
  mutate(sample_avg = guild_count / sum(guild_count)) %>% 
  ungroup() %>% 
  select(Guild, sample_avg) %>%  
  group_by(Guild) %>% 
  summarise(guild_avg = sum(sample_avg)) %>% 
  mutate(percentage = 100 * guild_avg / sum(guild_avg)) %>% 
  mutate(Guild = if_else(percentage < 2, "Other Guild >2%", Guild)) %>% 
  group_by(Guild) %>%                                                 
  summarise(percentage = sum(percentage)) %>%                         
  mutate(Guild = fct_reorder(Guild, percentage, .desc = FALSE)) %>%   
  arrange(desc(percentage))

# Combine both datasets to ensure consistent colors
all_guilds <- bind_rows(
  Indoor_guild_percentages %>% mutate(source = "Indoor"),
  Outdoor_guild_percentages %>% mutate(source = "Outdoor")
)

# Define a consistent color palette for all Guilds
unique_guilds <- unique(all_guilds$Guild)
guild_colors <- scale_fill_manual(
  values = setNames(ggsci::pal_d3("category20")(length(unique_guilds)), unique_guilds)
)

# Indoor plot
Indoor_pie_plot <- ggplot(Indoor_guild_percentages, aes(x = "", y = percentage, fill = Guild)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  theme_void() +
  guides(fill = guide_legend(reverse = TRUE, title = "Guild")) +
  guild_colors

# Outdoor plot
Outdoor_pie_plot <- ggplot(Outdoor_guild_percentages, aes(x = "", y = percentage, fill = Guild)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  theme_void() +
  guides(fill = guide_legend(reverse = TRUE, title = "Guild")) +
  guild_colors

# Combine the two plots
Indoor_pie_plot | Outdoor_pie_plot

Percentage_of_Guilds <- # Edit as required to view % of in/outdoor guilds 
  guild_counts %>% 
  inner_join(., METADATA) %>% 
  filter(type == "Outdoor") %>% 
  select("asv", "taxonomy", "guild_level", "Guild", "sample_id", "count") %>% 
  group_by(sample_id, Guild) %>% 
  summarise(guild_count = sum(count)) %>% 
  ungroup() %>% 
  group_by(sample_id) %>% 
  mutate(sample_avg = guild_count / sum(guild_count)) %>% 
  ungroup() %>% 
  select(Guild, sample_avg) %>%  
  group_by(Guild) %>% 
  summarise(guild_avg = sum(sample_avg)) %>% 
  mutate(percentage = 100 * guild_avg / sum(guild_avg)) %>% 
  mutate(Guild = if_else(percentage < 2, "Other Guild >2%", Guild)) %>% 
  group_by(Guild) %>%                                                  
  summarise(percentage = sum(percentage)) %>%                          
  mutate(Guild = fct_reorder(Guild, percentage, .desc = FALSE)) %>%    
  arrange(desc(percentage))
