# BINF6210 Assignment 1: Exploration of of Drosophila using BOLD ----
# Author: Iroayo Toki
# Published: November 4, 2025 
# Last updated: November 12, 2025
# Model Organism: Drosophilia


##Packages----
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(vegan)
library(viridis)
library(ggplot2)

#Assignment 1----
#1. Read in and transform----
#Reading in the .tsv file for drosophilia genus
df_bold <- read_tsv(file = "../data/Drosophilia_BOLD_data.tsv")
df_bold

names(df_bold)

#Seperates the coordinates(coord) into lat and long and converts to numeric values, regex searches for all characters that are digits, minus sign or decimal points and removes the others 
df_bold <- df_bold %>%
  separate(coord, into = c("Lat","Long"), sep = ",", convert = FALSE) %>%
  mutate(
    Lat  = as.numeric(str_replace_all(Lat, "[^0-9.-]", "")),
    Long = as.numeric(str_replace_all(Long, "[^0-9.-]", ""))
  )
class(df_bold$Lat)

#reduce the number of variables into a new dataframe
df_bold_sub <- df_bold[, c("processid", "bin_uri", "species", "country/ocean", 
                           "Lat", "Long", "marker_code", "nuc")]
df_bold_sub


#Create table with lat ranges and their corresponding regions, indicate whether or not to include the upper and lower bounds.
df_lat_regions <- data.frame(
  region = c("Polar", "Temperate", "Tropical", "Temperate", "Polar"),
  lat_max = c(-66.5, -23.5, 23.5, 66.5, 90),
  lat_min = c(-90, -66.5, -23.5, 23.5, 66.5),
  include_max = c(FALSE, FALSE, TRUE, TRUE, TRUE),
  include_min = c(TRUE, TRUE, TRUE, FALSE, TRUE)
)

#Create a new variable grouping latitude into Tropical, Temperate and Polar regions
df_bold_sub <- df_bold_sub %>% 
  rowwise() %>% 
  mutate(
    latgrouped = df_lat_regions$region[
      which(
        (ifelse(df_lat_regions$include_min, Lat >= df_lat_regions$lat_min, Lat > df_lat_regions$lat_min)) &
          (ifelse(df_lat_regions$include_max, Lat <= df_lat_regions$lat_max, Lat < df_lat_regions$lat_max))
      )
    ][1]) %>% 
  ungroup()
#Check the count for each region
#Tropical = 20004, Temperate = 8846, NA = 8002, Polar = 23
df_bold_sub %>% 
  count(latgrouped, sort = TRUE)

#Dataset showing grouping by BINS and their average latitude (median used to protect against outliers). 
df_bin_summary <- df_bold_sub %>%
  group_by(bin_uri) %>%
  summarise(avg_lat = median(Lat, na.rm = TRUE)) %>%
  rowwise() %>% 
  mutate(
    latgrouped = df_lat_regions$region[
      which(
        (ifelse(df_lat_regions$include_min, avg_lat >= df_lat_regions$lat_min, avg_lat > df_lat_regions$lat_min))
        &
          (ifelse(df_lat_regions$include_max, avg_lat <= df_lat_regions$lat_max, avg_lat < df_lat_regions$lat_max
                  )
          )
      )
    ][1]) %>% 
  ungroup()

sum(is.na(df_bin_summary$avg_lat)) # 173 BINs have no latitiude values


#2. Exploration to understand parameters----

length(unique(df_bold_sub$species))
length(unique(df_bold_sub$bin_uri))
length(unique(df_bold_sub$bin_uri))/length(unique(df_bold_sub$species)) #More species than BINs 

sum(!is.na(df_bold_sub$bin_uri))/length((df_bold_sub$bin_uri)) # 91.1% of samples have BINs

#checking that latitude values are within the right range (-90 to +90). Found that data is oversampled between 10-20° (Sampling bias may be at play, this is key for regression analysis and we try to reduce its effects by grouping data into latitude bands ).
hist(df_bold_sub$Lat)

#3 Data exploration to study relationship between latitudinal regions and biodiversity---- 
#Figures show increased biodiversity in tropical regions compared to temperate regions based on available data, this could be due to stabler climates and warmer temperatures and increased vegetation, we will further check for this using a statistical test

#Bar chart showing samples per region
ggplot(df_bold_sub) +
  geom_bar(aes(x = latgrouped), fill = "blue") +
  labs(title = "Samples per Region", x = "Latitude Group(Region)", y = "Count")

#Bar chart showing  Number of BINs per region
ggplot(df_bin_summary) +
  geom_bar(aes(x = latgrouped), fill = "brown") +
  labs(title = "BINs per Region", x = "Latitude Group(Region)", y = "Numbers of BINs found")

#Scatterplot showing Unique BINs latitude distribution
ggplot(df_bin_summary) +
  geom_point(mapping = aes(x = avg_lat, y = bin_uri, colour = latgrouped, shape = latgrouped), size = 2, na.rm = TRUE) +
  labs(title = "Distribution of Drosophilia BINs by latitude ", x = " Average Latitude", y = "Unique BINs") + scale_y_discrete(labels= NULL) + scale_color_viridis_d(option= "inferno") 

#4. Statistical test and visualization for BIN richness by latitude/geographic region----

#grouped into latiude band groups of every 10 degrees away from the equator in either direction with unique BINs and avg latitude for regression analysis
#Mean used for average since data points are grouped within a 10 degree range reducing effect of outliers
df_diversity_bins <- df_bold_sub %>%
  mutate(lat_band = cut(Lat, breaks = seq(-90, 90, by = 10))) %>%
  group_by(lat_band) %>%
  summarise(BIN_richness = length(unique(bin_uri)), avg_lat = mean(Lat, na.rm = TRUE)) %>%  mutate(latgrouped = case_when(
    avg_lat >= -23.5 & avg_lat <= 23.5 ~ "Tropical",
    (avg_lat > 23.5 & avg_lat <= 66.5) | avg_lat < (-23.5 & avg_lat >= -66.5) ~ "Temperate",avg_lat > 66.5 | avg_lat < -66.5 ~ "Polar")) %>% drop_na()

#moderately negative correlation the farther from the equator the less unique BINs we encounter
cor(df_diversity_bins$BIN_richness, abs(df_diversity_bins$avg_lat))

#p-value 0.09453 abs to account for positive and negative values
model <- lm(BIN_richness ~ abs(avg_lat), data = df_diversity_bins)
summary(model)

#Regression curve showing we encounter less unique BINs as move away from the equator with outliers close to the outer boundaries of tropical zone  around 20°(recall outliers from histogram)
ggplot(df_diversity_bins, aes(x = abs(avg_lat), y = BIN_richness, shape= latgrouped, color = latgrouped )) +
  geom_point(size = 2 ) +
  geom_smooth(method = "lm", se = TRUE, color = "black", aes(group = 1)) +
  labs(title = "Relationship between Latitude and BIN Richness",x = "Absolute Latitude (°)", y = "Number of Unique BINs") +  theme_minimal()

