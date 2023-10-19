## ---------------------------
##
## Script: filter_feature_table_138_r1-r5
##
## Purpose: filter feature table
##
## Author: Ikaia Leleiwi
##
## Date Created: December 1st, 2022
##
## Copyright (c) Ikaia Leleiwi, 2022
## Email: ileleiwi@gmail.com
##
## ---------------------------
##
## Notes:
##  #find matches in metadata sample ids and feature table columns
##  #remove all samples that have 0 in every row
##  #remove all rows that have 0 in every column
##  #remove any samples with <1000 total counts
##  #remove any rows that don't have >0 in at least 3 samples
##  #remove mitochondria and chloroplast 
## ---------------------------

## set working directory

setwd("~/Salmonella-Multiomics-Paper")

## ---------------------------

library(tidyverse)

#Data

feature_table <- read_tsv("Data/Amplicon/feature_table_138_r1-r5_git.tsv") 


taxonomy <- read_tsv("Data/Amplicon/taxonomy_138_r1-r5_git.tsv")


#find matches in metadata sample ids and feature table columns
feature_table_clean <- feature_table %>%
#remove all samples that have 0 in every ASV
  column_to_rownames(var = "asv_id") %>%
  t() %>%
  as.data.frame() %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  filter(total > 0) %>%
  select(-total) %>%
#remove all ASVs that have 0 in every sample
  t() %>%
  as.data.frame() %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  filter(total > 0) %>%
  select(-total) %>%
#remove any samples with <1000 total counts
  t() %>%
  as.data.frame() %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  filter(total >= 1000) %>%
  select(-total) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "asv_id")
 
#remove mitochondria and chloroplast and Unassigned
asv_ids_mito_chloro <- taxonomy %>%
  filter(str_detect(Taxon, "Mitochondria|Chloroplast|Unassigned|d__Eukaryota")) %>%
  pull(`Feature ID`)

feature_table_clean <- feature_table_clean %>%
  filter(!(asv_id %in% asv_ids_mito_chloro))

taxonomy_clean <- taxonomy %>%
  filter(!(`Feature ID` %in% asv_ids_mito_chloro))

write_tsv(feature_table_clean, "Data/Amplicon/feature_table_clean_138_r1-r5_git.tsv")
write_tsv(taxonomy_clean, "Data/Amplicon/taxonomy_clean_138_r1-r5_git.tsv")
