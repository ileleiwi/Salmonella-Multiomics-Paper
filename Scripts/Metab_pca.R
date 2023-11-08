## ---------------------------
##
## Script: Metab_pca
##
## Purpose: make pca plot from metabolite data
##
## Author: Ikaia Leleiwi
##
## Date Created: September 19, 2022
##
## Copyright (c) Ikaia Leleiwi, 2022
## Email: ileleiwi@gmail.com
##
## ---------------------------
##
## Notes: https://www.intechopen.com/chapters/52527
##   
##
## ---------------------------

## set working directory

setwd("~/Salmonella-Multiomics-Paper")

## ---------------------------

##Libraries

library(tidyverse)
library(vegan)
library(RColorBrewer)
library(ggnewscale)
library(gridExtra)

##Data
compound_metadata <- read_tsv("data/Metabolomics/Metabolite_metadata_group.tsv")

mouse_metadata <- read_tsv("data/Metabolomics/metab_mouse_metadata.tsv") 

ft <- read_tsv("data/Metabolomics/metab_feature_table_clean.tsv")

OUT_PATH <- "figures/PCA"

#PCA data cleaning functions
replacezero <- function(x) "[<-" (x, !x | is.na(x), min(x[x > 0], na.rm = TRUE) / 2)

paretoscale <- function(z){
  rowmean <- apply(z, 1, mean) # row means
  rowsd <- apply(z, 1, sd) #row standard deviation
  rowsqrtsd <- sqrt(rowsd) #sqrt of sd
  rv <- sweep(z, 1, rowmean, "-") #mean center
  rv <- sweep(rv, 1, rowsqrtsd, "/") #divide by sqrtsd
  return(rv)
}

#PCA
df_pca <- ft %>%
  column_to_rownames(var = "Compound") %>%
  mutate_if(is.character, as.numeric) %>%
  select(-starts_with("MV"), -starts_with("QC"), -starts_with("PRE"))

veg_ft <- ft %>%
  column_to_rownames(var = "Compound") %>%
  mutate_if(is.character, as.numeric) 

veg_pca <- rda(decostand(veg_ft, method = "hellinger"), scale = TRUE, center = TRUE)

#replace zeros with lowest value in matrix divided by 2
df_pca_nozero <- apply(df_pca, 1, replacezero)

log_df_pca <- log(df_pca_nozero, 2)
pareto_log_df_pca <- paretoscale(log_df_pca)

pca <- prcomp(pareto_log_df_pca, center = TRUE, scale = FALSE)
pcaresults <- summary(pca)

#permanova on bray distances
df_pca_meta <- df_pca %>%
  rownames_to_column(var = "compound") %>%
  pivot_longer(cols = -compound,
               names_to = "sample",
               values_to = "abund") %>%
  mutate(treatment = ifelse(str_detect(sample, "Salmonella"), "infected", "uninfected")) %>%
  select(-compound, -abund) %>%
  unique()


bray_dist <- vegdist(t(df_pca), method = "bray")
adonis2(bray_dist~treatment,data = df_pca_meta, method = "bray", permutations = 999)

#scores plot
scores <- as.data.frame(pcaresults$x) %>%
  rownames_to_column(var = "Compound") %>%
  separate(col = "Compound",
           into = c("sample", "mode", "rep", "treatment"),
           sep = "_")

scores %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = treatment), size = 3, shape = 21) +
  geom_text(aes(label = sample), nudge_x = 1, size = 3) +
  stat_ellipse(aes(fill = treatment, color = treatment), alpha = 0.2, geom = "polygon") +
  labs(x = "PC1 (63%)",
       y = "PC2 (7.5%)") +
  theme_bw()

scores %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = treatment), size = 3, shape = 21) +
  #geom_text(aes(label = sample), nudge_x = 1, size = 3) +
  stat_ellipse(aes(fill = treatment, color = treatment), alpha = 0.2, geom = "polygon")

#loadings plot
loadings <- as.data.frame(pcaresults$rotation) %>%
  mutate(pc1.change = ifelse(PC1 > 0.09, "up",
                             ifelse(PC1 < -0.09, "down",
                                    "nochange")))

significant_loadings <- loadings %>%
  subset(PC1 > 0.09 | PC1 < -0.09 | PC2 > 0.09 | PC2 < -0.09 )

sig_metadata <- compound_metadata %>%
  filter(Compound %in% rownames(significant_loadings))

loadings %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = pc1.change)) +
  labs(x = "PC1 (63%)",
       y = "PC2 (7.5%)") +
  theme_bw()

#combined
sig_compounds <- compound_metadata %>%
  filter(`Anova (p)` <= 0.05) %>%
  pull(Compound)


compound_names <- compound_metadata %>%
  select(Compound, c_name = `Accepted Annotation`, Group) 

#function to replace -Inf with lowest value in colum/2
repl_zero <- function(x){
  m <- min(x[x > 0])
  m_2 <- m/2
  x[x==-Inf] <- m_2
  return(x)
}

#euclidian distance funciton to get magnitude of vector arrows from 0,0
distance <- function(x, y, home = c(0,0)) {
  sqrt((x-home[1])^2 + (y-home[2])^2)
}


metab_abund <- df_pca %>%
  t() %>%
  log(., 2)

new <- matrix(nrow = nrow(metab_abund), ncol = ncol(metab_abund))
for(i in 1:ncol(metab_abund)){
  new[,i] <- repl_zero(metab_abund[,i])
}
colnames(new) <- colnames(metab_abund)
rownames(new) <- rownames(metab_abund)
metab_abund <- new

#loadings vectors
vl <- envfit(pca, metab_abund, perm = 1000, p.max = 0.005)

#create loadings df for plot arrows
v_scrs <- as.data.frame(scores(vl, display = "vectors"))
v_scrs <- cbind(v_scrs*ordiArrowMul(vl), Compound = rownames(v_scrs))

v_scrs <- v_scrs %>%
  filter(Compound %in% rownames(significant_loadings),
         Compound%in% sig_compounds) %>%
  left_join(compound_names, by = "Compound") %>%
  mutate(dist_pc1pc2 = distance(PC1, PC2), #distance in space from center
         scaled_dist = dist_pc1pc2 * 10)



#get top 10 compounds based on distance in space from center
top_dist_compounds <- v_scrs %>%
  arrange(desc(scaled_dist)) %>%
  slice_head(n = 10) %>%
  mutate(num_lab = 1:10) %>%
  select(Compound, c_name, num_lab)


#make plotting columns for loadings df
v_scrs <- v_scrs %>%
  mutate(top_c = ifelse(Compound %in% top_dist_compounds$Compound,"top","other"),
         alpha_factor = ifelse(top_c == "top", 1, 0.5),
         plot_pc1 = ifelse(top_c == "top", PC1, 0),
         plot_pc2 = ifelse(top_c == "top", PC2, 0),
         top_label = ifelse(top_c == "top", c_name, "")) %>%
  left_join(top_dist_compounds, by = c("Compound","c_name"))


#Group pallett == Set1

#Main PCA plot
pca_plot <- scores %>%
  ggplot(aes(x = PC1, y = PC2)) +
  coord_fixed() +
  geom_segment(data = v_scrs,
               aes(x= 0, xend = PC1*20, 
                   y = 0, yend = PC2*20, 
                   alpha = alpha_factor,
                   color = Group),
               arrow = arrow(length = unit(0.5, "cm")),
               inherit.aes = FALSE) +
  scale_color_brewer(palette = "Paired") +
  new_scale_color() +
  geom_text(data = v_scrs,
            aes(label = num_lab,
                x = PC1*20,
                y = PC2*20),
            size = 3) +
  geom_point(aes(fill = treatment), size = 3, shape = 21) +
  scale_fill_manual(values = c("#247169", "#8c2c26")) +
  scale_color_manual(values = c("#247169", "#8c2c26")) +
  stat_ellipse(aes(fill = treatment, color = treatment), alpha = 0.2, geom = "polygon") +
  labs(x = "PC1 (63%)",
       y = "PC2 (7.5%)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank())

#pdf(paste0(OUT_PATH,"/metab_pca_blank.pdf"), width = 7.5, height = 7.5)
pca_plot
#dev.off()

# pdf(paste0(OUT_PATH,"/metab_pca_legend.pdf"),
#     width = 7.5, height = 7.5)
# pca_plot
# dev.off()


limits_pc1 <- v_scrs %>%
  arrange(PC1) %>%
  pull(c_name)

#pc1 loadings barchart
pc1_loadings_bar <- v_scrs %>%
  ggplot(aes(x = reorder(c_name,PC1), y = abs(PC1))) +
  geom_bar(stat = "identity",
           aes(fill = Group)) +
  geom_text(aes(label = top_label)) +
  labs(y = "Scaled PC1 Loadings of Significant Compounds",
       x = NULL) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.position = "none")
# theme(legend.position = "none",
#       axis.text.x = element_blank(),
#       axis.text.y = element_blank(),
#       axis.title = element_blank())

#pdf(paste0(OUT_PATH, "/pc1_loadings_bar.pdf"), width = 7.5, height = 1.5)
pc1_loadings_bar
#dev.off()

#pc2 loadings barchart
pc2_loadings_bar <- v_scrs %>%
  ggplot(aes(x = abs(PC2), y = reorder(c_name,PC2))) +
  geom_bar(stat = "identity",
           aes(fill = Group)) +
  geom_text(aes(label = num_lab)) +
  labs(y = "Scaled PC2 Loadings of Significant Compounds",
       x = NULL) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        legend.position = "none")
# theme(legend.position = "none",
#       axis.text.x = element_blank(),
#       axis.text.y = element_blank(),
#       axis.title = element_blank())

#pdf(paste0(OUT_PATH, "/pc2_loadings_bar.pdf"), width = 1.5, height = 5.725)
pc2_loadings_bar
#dev.off()


empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())

#grid.arrange(top_bar, empty, pca_plot, nrow = 2, ncol = 2, widths=c(4, 1), heights=c(1, 4))
