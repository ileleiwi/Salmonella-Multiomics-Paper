## ---------------------------
##
## Script: 16S_stacked_bar_charts_metat_metab
##
## Purpose: make stacked barcharts for mice with metabolomics and metaT
##
## Author: Ikaia Leleiwi
##
## Date Created: October 18th, 2022
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

setwd(paste0("~/Salmonella-Multiomics-Paper"))

## ---------------------------

##Libraries

library(tidyverse)

##Data

ra <- read_tsv("Data/Amplicon/relabund_metab_metat_138.tsv")

metabo_mice <- c("3_X_13","3_X_15","3_X_16","3_X_20","3_X_27","3_X_28",
                 "3_X_44","3_X_45","3_X_46","3_X_47","3_X_48","3_X_49","3_X_50")
metat_mice <- c("3_X_13","3_X_15","3_X_16","3_X_28","3_X_32",
                "3_X_45","3_X_46","3_X_47","3_X_49","3_X_50")

union(metabo_mice, metat_mice) %>% length()

meta <- read_tsv("Data/Amplicon/metadata_r1-r5_paperid.tsv")
paper_id_join <- meta %>% select(sample, paper_id)

##Clean Data

ra_c <- ra %>%
  mutate(metat_metab_mice = case_when(mouse_id %in% union(metabo_mice, metat_mice) ~ "metat_metab",
                                      mouse_id %in% metat_mice ~ "metat",
                                      mouse_id %in% metabo_mice ~ "metab",
                                      T ~ "other"),
         Day = factor(Day, levels= c("-2","-1","0","10","11","12"))) %>%
  left_join(paper_id_join, by = "sample")


top_classes_late <- ra_c  %>%
  filter(Day %in% c("10","11","12")) %>%
  group_by(class) %>%
  summarise(mean_relabund = mean(relative_abundance)) %>%
  arrange(desc(mean_relabund)) %>%
  slice_head(n=8) %>%
  mutate(class = ifelse(class == "c__Lentisphaeria", "Other", class))

top_classes_early <- ra_c  %>%
  filter(Day %in% c("-2","-1","0")) %>%
  mutate(class = case_when(!(class %in% top_classes_late$class) ~ "Other",
                           T ~ class)) %>%
  group_by(class) %>%
  summarise(mean_relabund = mean(relative_abundance)) %>%
  arrange(desc(mean_relabund))


late_factor_levels <- top_classes_early %>%
  mutate(class = str_remove(class, "c__")) %>%
  arrange(desc(mean_relabund)) %>%
  pull(class)

metat_metab_join <- unique(ra_c[,c("paper_id", "metat_metab_mice", "treatment")])

plot_df <- ra_c %>%
  mutate(class = case_when(!(class %in% top_classes_late$class) ~ "Other",
                           T ~ class),
         class = str_remove(class, "c__"),
         class = factor(class, levels = (late_factor_levels))) %>%
  group_by(paper_id, Day, class) %>%
  summarise(mean_relabund = mean(relative_abundance), .groups = "drop") %>%
  left_join(metat_metab_join, by = "paper_id") %>%
  arrange(metat_metab_mice) %>%
  mutate(paper_id = factor(paper_id),
         paper_id = fct_reorder(paper_id, metat_metab_mice)) 



Verrucomicrobiae = "#f6b5d1"
#Lentisphaeria = "#ff6db6"
Bacteroidia = "#b6dbff"
Clostridia = "#490092"
Bacilli = "#006ddb"
Gammaproteobacteria = "#8b2f24"
Mollicutes = "#b66dff"
Dehalobacteriia = "#009292"
Other = "#000000"
Negativicutes = "#004949"
Syntrophia = "#009292"

Coriobacteriia = "#924900"

plot_df$class %>% levels()

pal <- c(Bacteroidia, Clostridia, Verrucomicrobiae, Bacilli, 
         Gammaproteobacteria, Negativicutes, Other, Syntrophia)
##Plot
plot_df %>%
  filter(treatment == "control") %>%
  ggplot(aes(x = paper_id, y = mean_relabund*100, fill = class)) +
  geom_bar(position = "fill", stat = "identity", color = "black") +
  scale_fill_manual(values = pal) +
  theme_classic() +
  theme(legend.position = "bottom") +
  facet_wrap(~Day, nrow = 1)
#dev.off()

plot_df %>%
  filter(treatment == "salmonella") %>%
  ggplot(aes(x = paper_id, y = mean_relabund*100, fill = class)) +
  geom_bar(position = "fill", stat = "identity", color = "black") +
  scale_fill_manual(values = pal) +
  theme_classic() +
  theme(legend.position = "bottom") +
  facet_wrap(~Day, nrow = 1)

## sample `48-(11)-C-F-LN` has only 75 total counts in asv feature table so it was lost in filtering


###Stats early late and treatment by class

stat_df <- plot_df %>%
  mutate(el = ifelse(Day %in% c(-2, -1, 0), "Early", "Late"))

stat_df_early <- stat_df %>%
  filter(el == "Early")

#independent class

willcox_stats <- function(df = stat_df, c, trt_el = T){
  
  if(trt_el){
    out <- df %>%
      filter(class == c) %>%
      group_by(treatment) %>%
      nest() %>%
      mutate(wilx_tst = map(data, 
                            ~wilcox.test(mean_relabund ~ el, 
                                         data = .x))) %>%
      mutate(summ = map(wilx_tst, broom::glance)) %>%
      unnest(summ) %>%
      ungroup()
  }else{
    out <- df %>%
      filter(class == c) %>%
      group_by(el) %>%
      nest() %>%
      mutate(wilx_tst = map(data, 
                            ~wilcox.test(mean_relabund ~ treatment, 
                                         data = .x))) %>%
      mutate(summ = map(wilx_tst, broom::glance)) %>%
      unnest(summ) %>%
      ungroup()
  }
  
  return(out)
  
}

#vector of classes
class_vec <- stat_df %>% pull(class) %>% unique() %>% as.character()


#early only, group by early late, relabund ~ treatment
wilx_early <- map(class_vec, ~willcox_stats(df = stat_df_early,
                                            c = .x,
                                            trt_el = F))
names(wilx_early) <- class_vec

#groupby treatment, relabund ~ early late
stat_wilx_trt_el <- map(class_vec, ~willcox_stats(c = .x))
names(stat_wilx_trt_el) <- class_vec

#groupby early late, relabund ~ treatment
class_vec <- stat_df %>% pull(class) %>% unique() %>% as.character()
stat_wilx_el_trt <- map(class_vec, ~willcox_stats(c = .x, trt_el = F))
names(stat_wilx_el_trt) <- class_vec

combine_dfs <- function(ls, grp){
  
  out_df <- data.frame(el = character(),
                       pval = numeric(),
                       metric = character())
  for(i in names(ls)){
    add_df <- data.frame(group = ls[[i]][[grp]],
                         pval = ls[[i]][["p.value"]],
                         metric = c(i,i))
    
    out_df <- rbind(out_df, add_df)
  }
  return(out_df)
}

diff_salm_vs_ctrl_df <- combine_dfs(stat_wilx_el_trt, grp = "el")
diff_early_vs_late_df <- combine_dfs(stat_wilx_trt_el, grp = "treatment")

#groupby early late, relabund ~ treatment
diff_salm_vs_ctrl_df %>%
  mutate(pval = format(pval, scientific = FALSE, big.mark = ","))

#groupby treatment, relabund ~ early late
diff_early_vs_late_df %>%
  mutate(pval = format(pval, scientific = FALSE, big.mark = ","))
