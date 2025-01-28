##Clear R##
rm(list = ls())

##Load libraries##
library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)
library(ggsci)
library(pals)
library(colorspace)
library(pals)
library(scales)

setwd("path/to/work/dir")

# Load qPCR results in
plate1 <- read_csv("path/to/metadata/2024-05-16_Panel1.csv") %>% as.tibble() 
plate2 <- read_csv("path/to/metadata/2024-05-17_plate2_123602.csv") %>% as.tibble() 
plate3 <- read.csv("path/to/metadata/2024-05-20_plate3_100322.csv") %>% as.tibble()
plate4 <- read.csv("path/to/metadata/2024-05-20_plate4_150825.csv") %>% as.tibble()
plate5 <- read.csv("path/to/metadata/2024-05-21_plate5_141219.csv") %>% as.tibble()
plate6 <- read.csv("path/to/metadata/2024-05-24_plate6_181023.csv") %>% as.tibble() 
plate7 <- read.csv("path/to/metadata/2024-05-29_plate7_124201.csv") %>% as.tibble()
plate8 <- read.csv("path/to/metadata/2024-07-11_plate8_170330.csv") %>% as.tibble()
plate9 <- read.csv("path/to/metadata/2024-07-15_plate9_131558.csv") %>% as.tibble()
plate10 <- read.csv("path/to/metadata/2024-07-16_plate10_163652.csv") %>% as.tibble()

plate1 <-
  plate1 %>% filter(Task == "UNKNOWN") %>%
  rename(sample_id = `Sample Name`,
         quantity_mean = `Quantity Mean`) %>%
  group_by(sample_id, quantity_mean) %>% summarise() %>% ungroup() %>%
  mutate(`Total GE` = 18 * quantity_mean) %>% select(-quantity_mean) %>% 
  mutate(across(everything(), ~ replace_na(., 0)))

plate2 <-
  plate2 %>% filter(Task == "UNKNOWN") %>%
  rename(sample_id = `Sample Name`,
         quantity_mean = `Quantity Mean`) %>%
  group_by(sample_id, quantity_mean) %>% summarise() %>% ungroup() %>%
  mutate(`Total GE` = 18 * quantity_mean) %>% select(-quantity_mean) %>%
  mutate(across(everything(), ~ replace_na(., 0)))

plate3 <-
  plate3 %>% filter(Task == "UNKNOWN") %>%
  rename(sample_id = `Sample.Name`,
         quantity_mean = `Quantity.Mean`) %>%
  group_by(sample_id, quantity_mean) %>% summarise() %>% ungroup() %>%
  mutate(`Total GE` = 18 * quantity_mean) %>% select(-quantity_mean) %>% 
  mutate(across(everything(), ~ replace_na(., 0)))

plate4 <- 
  plate4 %>% filter(Task == "UNKNOWN") %>%
  rename(sample_id = `Sample.Name`,
         quantity_mean = `Quantity.Mean`) %>%
  group_by(sample_id, quantity_mean) %>% summarise() %>% ungroup() %>% 
  mutate(`Total GE` = 18 * quantity_mean) %>% select(-quantity_mean) %>%
  mutate(across(everything(), ~ replace_na(., 0)))

plate5 <- 
  plate5 %>% filter(Task == "UNKNOWN") %>%
  rename(sample_id = `Sample.Name`,
         quantity_mean = `Quantity.Mean`) %>%
  group_by(sample_id, quantity_mean) %>% summarise() %>% ungroup() %>% 
  mutate(`Total GE` = 18 * quantity_mean) %>% select(-quantity_mean) %>% 
  mutate(across(everything(), ~ replace_na(., 0)))

plate6 <- 
  plate6 %>% filter(Task == "UNKNOWN") %>%
  rename(sample_id = `Sample.Name`,
         quantity_mean = `Quantity.Mean`) %>%
  group_by(sample_id, quantity_mean) %>% summarise() %>% ungroup() %>% 
  mutate(`Total GE` = 18 * quantity_mean) %>% select(-quantity_mean) %>%
  mutate(across(everything(), ~ replace_na(., 0)))

plate7 <- 
  plate7 %>% filter(Task == "UNKNOWN") %>%
  rename(sample_id = `Sample.Name`,
         quantity_mean = `Quantity.Mean`) %>%
  group_by(sample_id, quantity_mean) %>% summarise() %>% ungroup() %>% 
  mutate(`Total GE` = 18 * quantity_mean) %>% select(-quantity_mean) %>% 
  mutate(across(everything(), ~ replace_na(., 0)))

plate8 <- 
  plate8 %>% filter(Task == "UNKNOWN") %>%
  rename(sample_id = `Sample.Name`,
         quantity_mean = `Quantity.Mean`) %>%
  group_by(sample_id, quantity_mean) %>% summarise() %>% ungroup() %>% 
  mutate(`Total GE` = 18 * quantity_mean) %>% select(-quantity_mean) %>% 
  mutate(across(everything(), ~ replace_na(., 0)))

plate9 <- 
  plate9 %>% filter(Task == "UNKNOWN") %>%
  rename(sample_id = `Sample.Name`,
         quantity_mean = `Quantity.Mean`) %>%
  group_by(sample_id, quantity_mean) %>% summarise() %>% ungroup() %>% 
  mutate(`Total GE` = 18 * quantity_mean) %>% select(-quantity_mean) %>% 
  mutate(across(everything(), ~ replace_na(., 0)))

plate10 <- 
  plate10 %>% filter(Task == "UNKNOWN") %>%
  rename(sample_id = `Sample.Name`,
         quantity_mean = `Quantity.Mean`) %>%
  group_by(sample_id, quantity_mean) %>% summarise() %>% ungroup() %>% 
  mutate(`Total GE` = 18 * quantity_mean) %>% select(-quantity_mean) %>% #Multiply for full volume (36ul elution)
  mutate(across(everything(), ~ replace_na(., 0)))


metadata <- 
  read.csv("Path/to/metadata/Wellhomes_metadata.csv") %>% as.tibble() %>% 
  select(ID, Type, qPCR_plate, Sequenced, Campaign, Home, Season) %>% 
  filter(qPCR_plate %in% c("PLATE1", "PLATE2", "PLATE3", "PLATE4", "PLATE5", 
                           "PLATE6", "PLATE7", "PLATE8", "PLATE9", "PLATE10")) %>% 
  rename(sample_id = ID)

all_plates <- rbind(plate1, plate2, plate3, plate4, plate5, plate6, plate7, plate8, plate9, plate10) 
all_plates <- inner_join(all_plates, metadata) 
all_plates$sample_id

###############
## Formating ##
###############

all_plates <- 
  all_plates %>% 
  filter(Sequenced %in% c("Panel 1", "Panel 2", "Panel 3")) %>% filter(Type == "Indoor")

## Box Plot ##
qPCR_BOX <-
  all_plates %>%
  filter(Type == "Indoor") %>%
  ggplot(aes(x = Type, y = `Total GE`, fill = Type)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, outlier.shape = NA, alpha = 0.65, coef = 0) +
  geom_jitter(show.legend = FALSE, width = 0.25, shape = 21, color = "black", size = 2) +
  theme_classic() +
  labs(title = NULL, x = NULL, y = "Total Fungal GE") +
  theme(
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),  
    axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text = element_text(size = 14),  
    strip.background = element_blank(),
  ) +
  scale_y_continuous(
    trans = log2_trans(),
    labels = scales::comma,
    breaks = trans_breaks("log2", function(x) 2^x, n = 9),
    expand = expansion(mult = c(0, 0.05)) 
  ) +
  scale_fill_manual(values = c("Indoor" = "#9b72cb"))

summary_stats <- 
  all_plates %>% filter(Type == "Indoor") %>% 
  summarise(
    Min = min(`Total GE`),
    Q1 = quantile(`Total GE`, 0.25),
    Median = median(`Total GE`),
    Q3 = quantile(`Total GE`, 0.75),
    Max = max(`Total GE`),
    Mean = mean(`Total GE`),
    SD = sd(`Total GE`),
    IQR = IQR(`Total GE`),
    Lower_Bound = Q1 - 1.5 * IQR(`Total GE`),
    Upper_Bound = Q3 + 1.5 * IQR(`Total GE`)
  )
print(summary_stats)

# Save genome equivalence so absolute abundance can be calculated 
Metadata_GE <- all_plates
saveRDS(Metadata_GE, "path/to/work/dir/Metadata_GE.rds")


