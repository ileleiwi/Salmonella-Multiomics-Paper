## ---------------------------
##
## Script: 16S_lefse_lda_genus_data
##
## Purpose: create df of genus level 16S for lefse analysis
##
## Author: Ikaia Leleiwi
##
## Date Created: October 19th, 2022
##
## Copyright (c) Ikaia Leleiwi, 2022
## Email: ileleiwi@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory

setwd("~/Salmonella-Multiomics-Paper")

## ---------------------------

##Libraries

library(tidyverse)

##Data
mice <- read_tsv("Data/Amplicon/mice_in_paper.tsv")

keep_mice <- mice %>%
  filter(analysis_type %in% c("metat", "metab"),
         true_false == 1) %>%
  pull(mouse_id) %>%
  unique()

keep_samples <- mice %>%
  filter(analysis_type == "16s",
         day %in% c(10, 11, 12),
         mouse_id %in% keep_mice)

metadata <- read_tsv("Data/Amplicon/metadata_r1-r5.tsv") %>%
  rename("sample" = "#SampleID") %>%
  as.data.frame()

taxonomy <- read_tsv("Data/Amplicon/taxonomy_clean_138_r1-r5_git.tsv") %>%
  select("asv_id" = "Feature ID", Taxon)

ft <- read_tsv("Data/Amplicon/feature_table_clean_138_r1-r5_git.tsv") 


##Functions##
relabund <- function(df, columns = c(NA)) 
  #takes feature table and calculates relative abundance of columns, omitting NA's
  #considers only columns listed in columns argument (character vector of column names). 
  #columns (default) = all columns
{
  if (NA %in% columns){
    df <- sweep(df, 2, colSums(df, na.rm = TRUE), '/')
  }
  else {
    df_relabund <- df %>% select(all_of(columns)) %>%
      sweep(., 2, colSums(., na.rm = TRUE), '/')
    df[,columns] <- df_relabund
  }
  
  return(df)
}


# salmonella ASV 4cbfff144d4e7a4e0f4619ed505be070

##Clean Data##
#filter feature table and calculate relative abundance 

ft_c <- ft %>%
  left_join(taxonomy,
            by = "asv_id") %>%
  separate(col = "Taxon",
           into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
           remove = FALSE,
           sep = ";") %>%
  mutate(genus = ifelse(asv_id == "4cbfff144d4e7a4e0f4619ed505be070",
                        "Salmonella",
                        genus)) %>%
  drop_na(genus) %>%
  filter(!(genus %in% c("g__unculturedbacterium", "g__unclutured", 
                        "g__unidentified", "g__Undibacterium", 
                        "g__unidentified", "g__uncultured",
                        "g__unculturedmicroorganism"))) %>%
  select(-c(Taxon, domain, phylum, class, order, family, genus, species)) %>%
  column_to_rownames(var = "asv_id") %>%
  relabund() 
  

ft_cc <- ft_c %>%
  rownames_to_column(var = "asv_id") %>%
  left_join(taxonomy,
            by = "asv_id") %>%
  separate(col = "Taxon",
           into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
           remove = FALSE,
           sep = ";") %>%
  mutate(genus = ifelse(asv_id == "4cbfff144d4e7a4e0f4619ed505be070",
                        "Salmonella",
                        genus)) %>%
  mutate(genus = str_remove(genus, "g__"),
         genus = str_remove_all(genus, "\\[|\\]"),
         asv_genus = paste(asv_id, genus, sep = "_")) %>%
  pivot_longer(cols = 2:42,
               names_to = "Sample",
               values_to = "relative_abundance") %>%
  group_by(genus, Sample) %>%
  summarise(relabund_sum = sum(relative_abundance)) %>%
  pivot_wider(names_from = "Sample",
              values_from = relabund_sum) %>%
  rename("X.SampleID" = "genus")

vect <- colnames(ft_cc)
vect2 <- c()
for(i in vect){
  if(i == "X.SampleID"){
    vect2 <- c(vect2, "Treatment")
  }else{
    temp <- metadata[metadata$sample == i,"treatment"]
    vect2 <- c(vect2, temp)
  }
}

vdf <- as.data.frame(vect2) %>% t() %>% as.data.frame()
vdf1 <- as.data.frame(vect) %>% t() %>% as.data.frame()
colnames(vdf) <- colnames(ft_cc)
colnames(vdf1) <- colnames(ft_cc)

df_out <- rbind(vdf, vdf1, ft_cc)

write_tsv(df_out, "Data/Lefse_input_genus_16S_138_metab_metat.tsv", col_names = FALSE)

# Run in terminal
# conda install -c biobakery lefse
# conda activate lefse
# format_input.py Lefse_input_genus_16S_138_metab_metat.tsv Lefse_input_genus_16S_138_metab_metat.tsv.in -c 1 -u 2 -o 1000000
# run_lefse.py Lefse_input_genus_16S_138_metab_metat.tsv.in Lefse_input_genus_16S_138_metab_metat.tsv.in.res
# plot_res.py Lefse_input_genus_16S_138_metab_metat.tsv.in.res Lefse_input_genus_16S_138_metab_metat.tsv.in.res.png --dpi 300
