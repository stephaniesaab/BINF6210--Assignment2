# Add a header with some details about the file and project
# Title: 
# Author:
# Published:
# Last updated:
#
#




##Packages----
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(vegan)
library(viridis)
library(ggplot2)

#Assignment 1----
#1. Read in and transform----
df_bold <- read_tsv(file = "../data/Drosophilia_BOLD_data.tsv")
df_bold
#Reading in the .tsv file for drosophilia genus

names(df_bold)

df_bold <- df_bold %>%
  separate(coord, into = c("Lat","Long"), sep = ",", convert = FALSE) %>%
  mutate(
    Lat  = as.numeric(str_replace_all(Lat, "[^0-9.-]", "")),
    Long = as.numeric(str_replace_all(Long, "[^0-9.-]", ""))
  )
class(df_bold$Lat)
#Seperates the coordinates(coord) into lat and long and converts to numeric values, regex searches for all characters that are digits, minus sign or decimal points and removes the others 

df_bold.sub <- df_bold[, c("processid", "bin_uri", "species", "country/ocean", "Lat", "Long", "marker_code", "nuc")]
df_bold.sub
#reduce the number of variables into a new dataframe

df_bold.sub <- df_bold.sub %>%
  mutate(
    latgrouped=case_when(
      Lat <= 23.5 & Lat >= -23.5 ~ "Tropical",
      (Lat> 23.5 & Lat <= 66.5) | (Lat < -23.5 & Lat >= -66.5) ~ "Temperate",
      Lat > 66.5 | Lat < -66.5~ "Polar"))
df_bold.sub %>% count(latgrouped, sort = TRUE)
#Create a new variable grouping latitude into Tropical, Temperate and Polar regions

df_bin_summary <- df_bold.sub %>%
  group_by(bin_uri) %>%
  summarise(avg_lat = median(Lat, na.rm = TRUE)) %>%
  mutate(latgrouped = case_when(
      avg_lat >= -23.5 & avg_lat <= 23.5 ~ "Tropical",
      (avg_lat > 23.5 & avg_lat <= 66.5) | avg_lat < (-23.5 & avg_lat >= -66.5) ~ "Temperate",
      avg_lat > 66.5 | avg_lat < -66.5 ~ "Polar"))
#Dataset showing grouping by BINS and their average latitude (median used to protect against outliers). 

sum(is.na(df_bin_summary$avg_lat))
# 173 BINs have no latitiude values


#2. Exploration to understand parameters----

length(unique(df_bold.sub$species))
length(unique(df_bold.sub$bin_uri))
length(unique(df_bold.sub$bin_uri))/length(unique(df_bold.sub$species))
#More species than BINs 

sum(!is.na(df_bold.sub$bin_uri))/length((df_bold.sub$bin_uri))
# 91.1% of samples have BINs


hist(df_bold.sub$Lat)
#checking that latitude values are within the right range (-90 to +90). Data is oversampled between 10-20° (Sampling bias may be at play, this is key for regression analysis and we try to reduce its effects by grouping data into latitude bands ).
#3 Data exploration to study relationship between latitudinal regions and biodiversity---- 
ggplot(df_bold.sub) +
  geom_bar(aes(x = latgrouped), fill = "blue") +
  labs(title = "Samples per Region", x = "Latitude Group(Region)", y = "Count")
#Bar chart showing samples per region

ggplot(df_bin_summary) +
  geom_bar(aes(x = latgrouped), fill = "brown") +
  labs(title = "BINs per Region", x = "Latitude Group(Region)", y = "Numbers of BINs found")
#Bar chart showing  Number of BINs per region

ggplot(df_bin_summary) +
  geom_point(mapping = aes(x = avg_lat, y = bin_uri, colour = latgrouped, shape = latgrouped), size = 2, na.rm = TRUE) +
  labs(title = "Distribution of Drosophilia BINs by latitude ", x = " Average Latitude", y = "Unique BINs") + scale_y_discrete(labels= NULL) + scale_color_viridis_d(option= "inferno") 
#Scatterplot showing Unique BINs latitude distribution

#Figures show increased biodiversity in tropical regions compared to temperate regions based on available data, this could be due to stabler climates and warmer temperatures and increased vegetation, we will further check for this using a statistical test

#4. Statistical test and visualization for BIN richness by latitude/geographic region----
df_diversity_bins <- df_bold.sub %>%
  mutate(lat_band = cut(Lat, breaks = seq(-90, 90, by = 10))) %>%
  group_by(lat_band) %>%
  summarise(BIN_richness = length(unique(bin_uri)), avg_lat = mean(Lat, na.rm = TRUE)) %>%  mutate(latgrouped = case_when(
    avg_lat >= -23.5 & avg_lat <= 23.5 ~ "Tropical",
    (avg_lat > 23.5 & avg_lat <= 66.5) | avg_lat < (-23.5 & avg_lat >= -66.5) ~ "Temperate",avg_lat > 66.5 | avg_lat < -66.5 ~ "Polar")) %>% drop_na()
#grouped into latiude band groups of every 10 degrees away from the equator in either direction with unique BINs and avg latitude for regression analysis

#Mean used for average since data points are grouped within a 10 degree range reducing effect of outliers
cor(df_diversity_bins$BIN_richness, abs(df_diversity_bins$avg_lat))
#moderately negative correlation the farther from the equator the less unique BINs we encounter

model <- lm(BIN_richness ~ abs(avg_lat), data = df_diversity_bins)
summary(model)
#p-value 0.09453 abs to account for positive and negative values

ggplot(df_diversity_bins, aes(x = abs(avg_lat), y = BIN_richness, shape= latgrouped, color = latgrouped )) +
  geom_point(size = 2 ) +
  geom_smooth(method = "lm", se = TRUE, color = "black", aes(group = 1)) +
  labs(title = "Relationship between Latitude and BIN Richness",x = "Absolute Latitude (°)", y = "Number of Unique BINs") +  theme_minimal()
#Regression curve showing we encounter less unique BINs as move away from the equator with outliers close to the outer boundaries of tropical zone  around 20°(recall outliers from histogram)
