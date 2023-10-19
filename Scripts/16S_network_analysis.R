## ---------------------------
##
## Script: 16S_network_analysis
##
## Purpose: run a correlation network analysis on high responders
##
## Author: Ikaia Leleiwi
##
## Date Created: September 08 2023
##
## Copyright (c) Ikaia Leleiwi, 2023
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
library(igraph)
library(Hmisc)
library(Matrix)
library(viridis)
library(ggraph)

##Functions
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

##Data

#metadata
meta <- read_tsv("Data/Amplicon/metadata_r1-r5.tsv")

#taxonomy
tax <- read_tsv("Data/Amplicon/taxonomy_clean_138_r1-r5_git.tsv") %>%
  select(-Confidence)
colnames(tax) <- c("feature_id", "taxa")
tax <- tax %>%
  mutate(taxa = str_remove_all(taxa, ".__"))

#ASV feature table
ft <- read_tsv("Data/Amplicon/feature_table_clean_138_r1-r5_git.tsv")

#relative abundance table
rel <- ft %>%
  column_to_rownames(var = "asv_id") %>%
  relabund() %>%
  rownames_to_column(var = "asv_id")

salm_abund <- rel %>%
  filter(asv_id == "4cbfff144d4e7a4e0f4619ed505be070") %>%
  select(-asv_id) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample_id") %>%
  rename("salm_abund" = "V1") %>%
  mutate(salm_abund = salm_abund*100)

#clean data
#build unique ids
meta2 <- meta %>%
  mutate(Cage = as.character(Cage),
         Cage = ifelse(is.na(Cage), "X", Cage),
         new_samp_id = paste(Round, Cage, Mouse, Day, sep = "_")) %>%
  as.data.frame()

#build new column names with new unique ids
new_col_names_ft <- c("asv_id")
for(i in colnames(ft)[-1]){
  temp <- meta2 %>%
    filter(`#SampleID` == i) %>%
    pull(new_samp_id)
  new_col_names_ft <- c(new_col_names_ft, temp)
}

#assign new sample names
colnames(ft) <- new_col_names_ft

#join with taxonomy and combine by genus
ft_tax <- ft %>%
  left_join(tax, 
            by = c("asv_id" = "feature_id")) %>%
  separate(col = "taxa",
           into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
           remove = FALSE,
           sep = ";") %>%
  mutate(genus = ifelse(asv_id == "4cbfff144d4e7a4e0f4619ed505be070",
                        "Salmonella",
                        genus)) %>%
  pivot_longer(cols = any_of(meta2$new_samp_id),
               names_to = "sample",
               values_to = "count") 

ft_tax_genus <- ft_tax %>%
  group_by(genus, sample) %>%
  summarise(count = sum(count)) %>%
  pivot_wider(names_from = "sample",
              values_from = "count",
              values_fill = 0) %>%
  mutate(genus = str_remove(genus, " ")) #remove leading space

#filter to genera with at least 5000 total counts across all samples
ft_tax_genus_5000 <- ft_tax_genus %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  filter(!is.na(genus),
         total > 5000,
         !(genus %in% c("uncultured", "uncultured bacterium", "unidentified"))) %>%
  select(-total)


#make sample keep list with hr mice days 10-12 where salm abundance is >=25
#filter out gener a with 0 in all samples
#and filter out genera with 1000 or fewer counts

salm_25 <- salm_abund %>%
  filter(salm_abund >= 25) %>%#all these samples are from days 10-12
  pull(sample_id)

day10_12 <- meta2 %>%
  filter(Day %in% c(10, 11, 12),
         treatment == "salmonella") %>%
  pull(new_samp_id)

ft_10_12_genus <- ft_tax_genus_5000 %>%
  select(genus, any_of(day10_12)) %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  filter(total > 1000) %>%
  select(-total)


#make highresponder relabund at genus level
rel_gen <- rel %>%
  left_join(tax, by = c("asv_id" = "feature_id")) %>%
  separate(col = "taxa",
           into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
           remove = FALSE,
           sep = ";") %>%
  mutate(genus = ifelse(asv_id == "4cbfff144d4e7a4e0f4619ed505be070",
                        "Salmonella",
                        genus),
         genus = str_remove(genus, "g__"),
         genus = str_replace(genus, " ", "")) %>%
  pivot_longer(cols = any_of(meta2$new_samp_id),
               names_to = "sample",
               values_to = "relabund") %>%
  left_join(meta2, by = c("sample"="new_samp_id")) %>%
  select(sample, genus, relabund) %>%
  filter(sample %in% day10_12,
         genus != "uncultured") %>%
  group_by(genus, sample) %>%
  summarise(relabund = sum(relabund)) %>%
  pivot_wider(names_from = "sample",
              values_from = "relabund",
              values_fill = 0) %>%
  filter(genus %in% ft_10_12_genus$genus)

#correlation analysis
ft_10_12_genus_mat <- ft_10_12_genus %>%
  column_to_rownames(var = "genus") %>%
  as.matrix() %>%
  t()

gen_cor <- rcorr(ft_10_12_genus_mat, type = "spearman")

#significance
gen_pval <- forceSymmetric(gen_cor$P) %>% as.matrix()
p_sig_TF <- gen_pval < 0.05 #significant correlation = TRUE, NS correlation = F

#correlation
gen_r_val <- gen_cor$r
gen_r_val_TF <- gen_r_val > 0 #positive correlation = TRUE, negative correlation = F
gen_r_val_pos <- gen_r_val_TF * gen_r_val #remove negative correlations, all negative == 0
p_sig_r <- gen_r_val_pos  *p_sig_TF #remove non significant correlations, all NS == 0

#adjacency matrix
adjm <- as.matrix(p_sig_r)

#find bugs that interact with salmonella and bugs that interact with those bugs with at least as correlation coefficient of 0.5
salm_interactions <- c(which(adjm[,"Salmonella"] != 0), which(is.na(adjm[,"Salmonella"])))
salm_interactions

entero_interactions <- c(which(adjm[,"Enterococcus"] != 0), which(is.na(adjm[,"Enterococcus"])))
entero_interactions

lacto_interactions <- c(which(adjm[,"Lactobacillus"] != 0), which(is.na(adjm[,"Lactobacillus"])))
lacto_interactions

keep_interactions <- unique(c(names(salm_interactions), 
                              names(entero_interactions), 
                              names(lacto_interactions)))

other_interactions <- colnames(adjm)[which(!(colnames(adjm) %in% keep_interactions))]


adjm_sig <- adjm[keep_interactions, keep_interactions]
adjm_other <- adjm[other_interactions, other_interactions]


#graph objects 
net.grph <- graph.adjacency(adjm, mode = "min", 
                            weighted = TRUE, diag = FALSE)
edgew <- E(net.grph)$weight
no_interaction_nodes <- V(net.grph)[degree(net.grph) == 0]

net.grph <- delete.vertices(net.grph, no_interaction_nodes)

#salmonella objects
net.grph.salm <- delete.vertices(net.grph, !(V(net.grph)$name %in% keep_interactions))
#size vector
size_order_salm <- V(net.grph.salm)$name

rel_abund_vect_salm <- rel_gen %>%
  filter(genus %in% size_order_salm) %>% 
  mutate(mean_relabund = rowMeans(across(where(is.numeric))),
         mean_relabund = mean_relabund*100) %>%
  pull(mean_relabund) 

rel_abund_df_salm <- rel_gen %>%
  filter(genus %in% size_order_salm) %>% 
  mutate(mean_relabund = rowMeans(across(where(is.numeric))),
         mean_relabund = mean_relabund*100) %>%
  select(genus, mean_relabund)

relabund_factor_salm <- rel_abund_df_salm %>%
  arrange((mean_relabund)) %>%
  pull(genus)

#graph objects other
net.grph.other <- delete.vertices(net.grph, !(V(net.grph)$name %in% other_interactions))
#size vector
size_order_other <- V(net.grph.other)$name

rel_abund_vect_other <- rel_gen %>%
  filter(genus %in% size_order_other) %>% 
  mutate(mean_relabund = rowMeans(across(where(is.numeric))),
         mean_relabund = mean_relabund*100) %>%
  pull(mean_relabund) 

rel_abund_df_other<- rel_gen %>%
  filter(genus %in% size_order_other) %>% 
  mutate(mean_relabund = rowMeans(across(where(is.numeric))),
         mean_relabund = mean_relabund*100) %>%
  select(genus, mean_relabund)


relabund_factor_other <- rel_abund_df_other %>%
  arrange((mean_relabund)) %>%
  pull(genus)

##plots salm

#min max scaling function
min_max <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

rel_abund_vect_min_max_salm <- (min_max(rel_abund_vect_salm) + .5)*20
#plot

#network plot salmonella correleated
plot(net.grph.salm, 
     vertex.size=rel_abund_vect_min_max_salm,
     vertex.frame.color="black",
     edge.curved=F,
     edge.width=1.5,
     layout=layout_on_grid,
     edge.color=ifelse(edgew < 0,"red","blue"),
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.font=2)
legend("topright",legend = "vertex")

#Heatmap of other correlated bacteria in community
rel_abund_df_other %>%
  mutate(xcol = 1,
         genus = factor(genus, levels = relabund_factor_other)) %>%
  ggplot(aes(x = xcol,y = genus)) +
  geom_tile(aes(fill = mean_relabund)) +
  geom_text(aes(label = round(mean_relabund,2),
                color = ifelse(mean_relabund > 5,"dark gray","light gray")),
            size = 2,
            show.legend = FALSE) +
  coord_equal() +
  scale_fill_viridis(option = "rocket") +
  scale_color_manual(values = c("black", "light gray")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom")


