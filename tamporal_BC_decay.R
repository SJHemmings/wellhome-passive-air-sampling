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
library(geosphere)
library(mgcv)
library(broom) 

## Set Working Dir ##
set.seed(123456) 
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
  select(-total_reads) %>% ungroup(sample_id) %>%  
  pivot_wider(names_from = asv, values_from = reads) %>% 
  as.data.frame()

############################
## Put data into a matrix ##
############################

rownames(FILT_ASV_COUNTS) <- FILT_ASV_COUNTS$sample_id 
FILT_ASV_COUNTS <- FILT_ASV_COUNTS[, -1]
ENV.matrix <- as.matrix(FILT_ASV_COUNTS) 

################################
## Temporal Bray-Curtis Decay ##
################################

ENV.dist <- avgdist(ENV.matrix, dmethod = "bray", sample = 15569)
pcoa <- cmdscale(ENV.dist, k=2, eig = TRUE, add = TRUE)

ENV.dist
METADATA <- 
  METADATA %>% 
  select(sample_id, home, collection_date, long, lat, type) %>% glimpse()

## Convert collection_date to Date format ##
METADATA <- 
  METADATA %>%
  filter(type %in% c("Indoor", "Outdoor")) %>% 
  mutate(collection_date = as.Date(collection_date, format = "%d/%m/%Y"))

# Calculate the number of days since the earliest date
earliest_date <- min(METADATA$collection_date, na.rm = TRUE)
METADATA <- METADATA %>%
  mutate(days = as.numeric(difftime(collection_date, earliest_date, units = "days")))
glimpse(METADATA)

## Create a matrix of latitude and longitude ##
coords <- METADATA %>% select(long, lat)

# Calculate the pairwise distance matrix (in meters) between all locations
geo_dist <- distm(coords, fun = distHaversine) / 1000  

# Convert to a dataframe and add sample IDs as rownames and colnames
geo_dist <- as.data.frame(geo_dist)
rownames(geo_dist) <- METADATA$sample_id
colnames(geo_dist) <- METADATA$sample_id

#############################################################
## Combine Temporal and Geographical data with bray curtis ##
#############################################################

# Convert the Bray-Curtis distance matrix (ENV.dist) into a long format
bray_curtis_long <-
  as.data.frame(as.table(as.matrix(ENV.dist))) %>%
  rename(sample_id1 = Var1, sample_id2 = Var2, bray_curtis = Freq)

# Convert the geographic distance matrix into a long format
geo_dist_long <-
  as.data.frame(as.table(as.matrix(geo_dist))) %>%
  rename(sample_id1 = Var1, sample_id2 = Var2, geo_dist_km = Freq)

# Join both dataframes
combined_df <-
  bray_curtis_long %>%
  left_join(geo_dist_long, by = c("sample_id1", "sample_id2"))

# Now add the temporal differences (days apart)
temporal_diff <- METADATA %>%
  select(sample_id, days)

combined_df <- combined_df %>%
  left_join(temporal_diff, by = c("sample_id1" = "sample_id")) %>%
  rename(days_sample1 = days) %>%
  left_join(temporal_diff, by = c("sample_id2" = "sample_id")) %>%
  rename(days_sample2 = days) %>%
  mutate(days_diff = abs(days_sample1 - days_sample2))

# Add sample types
sample_types <- METADATA %>%
  select(sample_id, type)

combined_df <- combined_df %>%
  left_join(sample_types, by = c("sample_id1" = "sample_id")) %>%
  rename(type_sample1 = type) %>%
  left_join(sample_types, by = c("sample_id2" = "sample_id")) %>%
  rename(type_sample2 = type)

# Filter out duplicate comparisons (same sample pairs)
combined_df <- combined_df %>%
  filter(sample_id1 != sample_id2)

# Glimpse the final combined dataframe
glimpse(combined_df)

####################
## r and r2 stats ##
####################
# Subset for indoor samples
indoor_data <- combined_df %>%
  filter(type_sample1 == "Indoor" & type_sample2 == "Indoor")

# Subset for outdoor samples
outdoor_data <- combined_df %>%
  filter(type_sample1 == "Outdoor" & type_sample2 == "Outdoor")

#####################################################################################
# Create Bray-Curtis distance matrices for indoor and outdoor samples Geographical ##
#####################################################################################

indoor_bray_matrix <- reshape2::acast(indoor_data, sample_id1 ~ sample_id2, value.var = "bray_curtis")
outdoor_bray_matrix <- reshape2::acast(outdoor_data, sample_id1 ~ sample_id2, value.var = "bray_curtis")

#########################################################
## Subset Indoor and Outdoor Samples for Temporal Data ##
#########################################################

# Subset for indoor samples based on temporal data
indoor_data_temporal <- combined_df %>%
  filter(type_sample1 == "Indoor" & type_sample2 == "Indoor")

# Subset for outdoor samples based on temporal data
outdoor_data_temporal <- combined_df %>%
  filter(type_sample1 == "Outdoor" & type_sample2 == "Outdoor")

##################################################
## Calculate Mantel Statistic for Temporal Data ##
##################################################

# Create Bray-Curtis distance matrices for indoor and outdoor samples (temporal)
indoor_bray_matrix_temporal <- reshape2::acast(indoor_data_temporal, sample_id1 ~ sample_id2, value.var = "bray_curtis")
outdoor_bray_matrix_temporal <- reshape2::acast(outdoor_data_temporal, sample_id1 ~ sample_id2, value.var = "bray_curtis")

# Create Temporal distance matrices for indoor and outdoor samples
indoor_temporal_matrix <- reshape2::acast(indoor_data_temporal, sample_id1 ~ sample_id2, value.var = "days_diff")
outdoor_temporal_matrix <- reshape2::acast(outdoor_data_temporal, sample_id1 ~ sample_id2, value.var = "days_diff")

# Mantel test for indoor samples (temporal)
mantel_indoor_temporal <- vegan::mantel(as.dist(indoor_bray_matrix_temporal), as.dist(indoor_temporal_matrix))
mantel_indoor_temporal  

# Mantel test for outdoor samples (temporal)
mantel_outdoor_temporal <- vegan::mantel(as.dist(outdoor_bray_matrix_temporal), as.dist(outdoor_temporal_matrix))
mantel_outdoor_temporal  

################################################################
## Calculate Pearson Correlation (R) and R² for Temporal Data ##
################################################################

# Pearson correlation (R) and R² for indoor samples (temporal)
cor_indoor_temporal <- cor(indoor_data_temporal$days_diff, indoor_data_temporal$bray_curtis, method = "pearson")
r_squared_indoor_temporal <- cor_indoor_temporal^2  

# Pearson correlation (R) and R² for outdoor samples (temporal)
cor_outdoor_temporal <- cor(outdoor_data_temporal$days_diff, outdoor_data_temporal$bray_curtis, method = "pearson")
r_squared_outdoor_temporal <- cor_outdoor_temporal^2  

# Print R and R² for both indoor and outdoor samples (temporal)
cat("Indoor Samples (Temporal): R =", cor_indoor_temporal, ", R² =", r_squared_indoor_temporal, "\n")
cat("Outdoor Samples (Temporal): R =", cor_outdoor_temporal, ", R² =", r_squared_outdoor_temporal, "\n")

#############################
## ggplot BC decay vs time ##
#############################

# Indoor and Outdoor colors
color_indoor <- "#9b72cb"
color_outdoor <- "#7BB662"

# Bray-Curtis vs Temporal Differences for Indoor and Outdoor (with fitted gam)
Bray_Curtis_Temporal_Indoor <- 
  combined_df %>%
  filter(type_sample1 == "Indoor" & type_sample2 == "Indoor") %>%
  ggplot(aes(x = days_diff, y = bray_curtis, color = type_sample1)) +
  geom_point(alpha = 0.6, size = 0.5) +
  scale_color_manual(values = c("Indoor" = color_indoor, "Outdoor" = color_outdoor)) +
  geom_smooth(method = "gam", 
              formula = y ~ s(x, k = 8), 
              method.args = list(family = quasibinomial()),
              color = "black") + 
  geom_vline(xintercept = 365, linetype = "dashed", color = "black") +
  labs(x = "Days Difference", y = "Bray-Curtis Distance", color = "Sample Type") +
  theme_classic() +
  theme(panel.background = element_rect(colour = "black", size = 0.1), legend.position = "none")

Bray_Curtis_Temporal_Outdoor <- 
  combined_df %>%
  filter(type_sample1 == "Outdoor" & type_sample2 == "Outdoor") %>%  
  ggplot(aes(x = days_diff, y = bray_curtis, color = type_sample1)) +
  geom_point(alpha = 0.6, size = 0.5) +
  scale_color_manual(values = c("Indoor" = color_indoor, "Outdoor" = color_outdoor)) +
  geom_smooth(method = "gam", 
              formula = y ~ s(x, k = 8), 
              method.args = list(family = quasibinomial()),
              color = "black") +  
  geom_vline(xintercept = 365, linetype = "dashed", color = "black") +
  labs(x = "Days Difference", y = "Bray-Curtis Distance", color = "Sample Type") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
   theme(panel.background = element_rect(colour = "black", size = 0.1), legend.position = "none")

###############################################################################
###############################################################################

#################################
## Outdoor temporal Statistics ##
#################################

# Filter the data
outdoor_data <- combined_df %>%
  filter(type_sample1 == "Outdoor" & type_sample2 == "Outdoor")  

# Fit the GAM model
gam_outdoor_model <- gam(bray_curtis ~ s(days_diff, k=8), data = outdoor_data, family = quasibinomial())

# Summary of the model to get R-squared
outdoor_model_summary <- summary(gam_outdoor_model)

# Calculate R-squared
outdoor_r_squared <- outdoor_model_summary$r.sq

# Calculate Pearson correlation coefficient (R)
predicted_values <- predict(gam_outdoor_model)
outdoor_r <- cor.test(outdoor_data$bray_curtis, predicted_values)
outdoor_r <- cor(outdoor_data$bray_curtis, predicted_values)

# # Calculate Pearson correlation coefficient (R)
# mantel(indoor_data$bray_curtis, predicted_values)

# Print the results
cat("R-squared:", outdoor_r_squared, "\n") #r2 outdoor
mantel_outdoor_temporal
outdoor_model_summary

#################################################################################
#################################################################################

################################
## Indoor temporal Statistics ##
################################

# Filter the data
outdoor_data <- combined_df %>%
  filter(type_sample1 == "Indoor" & type_sample2 == "Indoor")  

# Fit the GAM model
gam_indoor_model <- gam(bray_curtis ~ s(days_diff, k=8), data = indoor_data, family = quasibinomial())

# Summary of the model to get R-squared
indoor_model_summary <- summary(gam_indoor_model)

# Calculate R-squared
indoor_r_squared <- indoor_model_summary$r.sq

# Calculate Pearson correlation coefficient (R)
predicted_values <- predict(gam_indoor_model)
outdoor_r <- cor.test(indoor_data$bray_curtis, predicted_values)
outdoor_r <- cor(indoor_data$bray_curtis, predicted_values)

# Print the results
cat("R-squared:", indoor_r_squared, "\n") #r2 outdoor
mantel_indoor_temporal
indoor_model_summary
