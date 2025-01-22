#Clear R and load libraries
rm(list = ls())
library(tidyverse)
library(readxl)
library(phyloseq)
library(ComplexHeatmap)
library(vegan)
library(dendsort)
library(colorRamp2)

# Read and preprocess metadata
metadata <- read.csv("Wellhomes_metadata.csv") %>%
  as_tibble() %>%
  rename_all(tolower) %>%
  filter(sequenced %in% c("Panel 1", "Panel 2", "Panel 3")) %>%
  rename(sample_id = id)

# Read ASV count table
asv_counts <- read.csv("path/to/dada2/output/ASVs_counts.csv") %>%
  as_tibble()

# Read and preprocess taxonomy table
taxonomy <- read_tsv("path/to/dada2/output/ASVs_taxonomy.tsv") %>%
  as_tibble() %>%
  select(-"...1") %>%
  mutate_at(vars(-ASV), ~ ifelse(!is.na(.), str_sub(., start = 4), .)) %>%
  rename_all(tolower) %>%
  mutate(asv = str_replace_all(asv, ">", ""))

# Filter and process ASV counts
filtered_asv <- asv_counts %>%
  pivot_longer(-asv, names_to = "sample_id", values_to = "reads") %>%
  rename_all(tolower) %>%
  inner_join(taxonomy, by = "asv") %>%
  filter(kingdom == "Fungi") %>%
  inner_join(metadata, by = "sample_id") %>%
  filter(type != "Standard") %>%
  group_by(sample_id) %>%
  mutate(total_reads = sum(reads)) %>%
  filter(total_reads > 6000) %>%
  ungroup() %>%
  select(asv, sample_id, reads)

# Calculate relative abundance
asv_rel_abund <- filtered_asv %>%
  inner_join(metadata, by = "sample_id") %>%
  inner_join(taxonomy, by = "asv") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = reads / sum(reads)) %>%
  ungroup() %>%
  pivot_longer(
    cols = c("kingdom", "phylum", "class", "order", "family", "genus", "species", "asv"),
    names_to = "level",
    values_to = "taxon"
  ) %>%
  mutate(taxon = if_else(is.na(taxon), "Unknown Fungal Classification", taxon))

# Aggregate data at the genus level
genus_rel_abund <- asv_rel_abund %>%
  filter(level == "genus") %>%
  group_by(sample_id, taxon, type, campaign, season) %>%
  summarise(rel_abund = 100 * sum(rel_abund), .groups = "drop")

# Identify top genera and pool others
genus_pool <- genus_rel_abund %>%
  group_by(taxon, type) %>%
  summarise(mean = mean(rel_abund), .groups = "drop") %>%
  group_by(taxon) %>%
  summarise(pool = max(mean) < 1, .groups = "drop")

heatmap_data <- genus_rel_abund %>%
  inner_join(genus_pool, by = "taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  filter(!(taxon %in% c("Unknown Fungal Classification", "Other"))) %>%
  group_by(taxon, sample_id) %>%
  summarise(rel_abund = sum(rel_abund), .groups = "drop") %>%
  pivot_wider(names_from = sample_id, values_from = rel_abund, values_fill = 0)

# Prepare heatmap matrix
row_names <- heatmap_data$taxon
heatmap_matrix <- as.matrix(heatmap_data[-1])
rownames(heatmap_matrix) <- row_names

# Generate dendrograms
row_dend <- as.dendrogram(hclust(vegdist(heatmap_matrix, method = "bray")))
col_dend <- dendsort(hclust(vegdist(t(heatmap_matrix), method = "bray")))

# Define color function
max_value <- max(heatmap_matrix, na.rm = TRUE)
col_fun <- colorRamp2(c(0, 2.5, max_value), c("white", "lavender", "#8E44AD"))

# Define sample colors
sample_colors <- ifelse(grepl("^A", colnames(heatmap_matrix)), "#7BB662", "#9b72cb")
names(sample_colors) <- colnames(heatmap_matrix)
column_annotation <- HeatmapAnnotation(Environment = anno_simple(colnames(heatmap_matrix), col = sample_colors))

# Create heatmap
Heatmap(
  heatmap_matrix,
  name = "Rel_abund",
  col = col_fun,
  cluster_columns = col_dend,
  cluster_rows = row_dend,
  column_dend_height = unit(2, "cm"),
  row_dend_width = unit(2, "cm"),
  row_split = 3,
  column_split = 3,
  row_names_gp = gpar(fontsize = 16),
  column_names_gp = gpar(fontsize = 6),
  column_labels = as.expression(lapply(colnames(heatmap_matrix), function(x) bquote(italic(.(x))))),
  top_annotation = column_annotation,
  heatmap_legend_param = list(title = "Rel_abund", position = "left")
)
