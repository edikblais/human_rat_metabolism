# This script reproduces parts of Figures 2 and 3
# of the manuscript. 

library(devtools)



# All helper functions defined in this project
# have the prefix: ef (edik function)
source("ef_helper_functions.R")


#######################################################
################### Load packages  ####################
#######################################################

library(dplyr)
library(reshape2)
library(ggplot2)

# project.folder = "C:/Dropbox/ratcon1_appendix/rdy/"

#######################################################
############# Load Expanded View Data  ##########
#######################################################

# curated reaction and gpr info annotated in Table EV1
rxn.info = c("MSB-15-6464_TableEV1.txt") %>% 
  read.table(sep = "\t", quote = "", header = T) %>% ef_df
gpr.size.curated = rxn.info %>% 
  select(rxn_id,n_hsa = n_gene_hsa, n_rno = n_gene_rno)

# orthology data available in Table EV3
gpr.orthology.data.load = c("MSB-15-6464_TableEV3.txt") %>% 
  read.table(sep = "\t", quote = "", header = T) %>% ef_df

# Calculate gpr sizes obtained by replacing all human genes with rat orthologs
gpr.size.all = gpr.orthology.data.load %>% 
  group_by(rxn_id) %>% 
  summarize(n_hsa = length(unique(hsa[!hsa %in% c("0","",NA)])),
            n_rno = length(unique(rno[!rno %in% c("0","",NA)]))) %>% 
  ungroup

# average GPR size ratio = 1.36
# between number of rat genes / human genes
# for all shared reactions in curated models
gpr.size.all %>% 
  filter(n_rno > 0, n_hsa > 0) %>% 
  summarize(mean_gpr_ratio = mean(n_rno / n_hsa),
            mean_gpr_difference = mean(n_rno - n_hsa))

# average GPR size ratio = 1.04 
# between number of rat genes / human genes
# for all shared reactions in curated models
gpr.size.curated %>% 
  filter(n_rno > 0, n_hsa > 0) %>% 
  summarize(mean_gpr_ratio = mean(n_rno / n_hsa),
            mean_gpr_difference = mean(n_rno - n_hsa))

# average GPR size ratio = 1.04 
# between number of rat genes / human genes
# for shared reactions that had GPR info in draft models
gpr.size.curated %>% 
  filter(rxn_id %in% gpr.size.all[["rxn_id"]]) %>%
  filter(n_rno > 0, n_hsa > 0) %>% 
  summarize(mean_gpr_ratio = mean(n_rno / n_hsa),
            mean_gpr_difference = mean(n_rno - n_hsa))

gpr.orthology.data.load %>% select(rxn_id) %>% distinct


# GPR conversion algorithm filters orthology annotations
# based on two parameters:
########################## database.minimum
# database.minimum = the mininum number of database 
# in which a unique ortholog pair must occur.
# orthologs were obtained from 5 databases,
# so we assigned a higher confidence to 
# human genes mapped to rat genes in multiple databases
# this value ranged from 0 - 5

########################## rank.maximum
# rank.maximum = maximum number of allowed 
# rat orthologs to be replaced by a human gene
# the score that was used to rank rat genes 
# was determined by the following equation:
# gene score = 5*(status+annotation+evidence+type)
# this value ranged from 0 - 100

# The prioritization of orthologs depended on bulk downloads 
# from multiple databases. The final scores are included, 
# but we did not receive authorization to distribute 
# these databases in whole so the methods that we used are
# described below:
###### type = NCBI gene type
# protein-coding evidence = 5
# no protein-coding evidence = 0
###### status = Ensembl gene summary status
# KNOWN = 5
# KNOWN_BY_PROJECTION = 4
# NOVEL = 3
# PUTATIVE = 1
# no data = 0
###### annotation = UniProt annotation score 
# value assigned was same as UniProt annotation score, otherwise 0
###### evidence = UniProt gene evidence
# Evidence at protein level = 5
# Evidence at transcript level = 4
# Inferred from homology = 3
# Predicted = 2
# Uncertain = 1
# no data = 0

hmr2.rxns = rxn.info %>% select(rxn_id,hmr2_id) %>% filter(nchar(hmr2_id) > 0)
# Apply a reasonable range of cutoffs for our two parameters:
gpr.conversion.result = lapply(c(1:9),function(rank.maximum) {
  lapply(c(1:5),function(database.minimum) {
    hmr2.rxns %>% select(rxn_id) %>% distinct %>% 
      mutate(cutoff_rank = rank.maximum, 
             cutoff_database = database.minimum) %>% 
      left_join(gpr.orthology.data.load %>% 
                  mutate(rno = ifelse(hsa_rank <= rank.maximum, rno, "0"),
                         rno = ifelse(orthology_count >= database.minimum, rno, "0")) %>% 
                  select(rxn_id, hsa, rno, hsa_score,rno_score), by = "rxn_id")}) %>% 
    rbind_all}) %>% rbind_all %>%
  group_by(rxn_id) %>% mutate(n_hsa = length(unique(setdiff(hsa,c(NA, "0", ""))))) %>% ungroup

gpr.all.orthology = hmr2.rxns %>% select(rxn_id) %>% distinct %>% 
  mutate(cutoff_rank = Inf, cutoff_database = 0) %>% 
  left_join(gpr.orthology.data.load %>% select(rxn_id, hsa, rno, hsa_score,rno_score), by = "rxn_id") %>% 
  group_by(rxn_id) %>% 
  mutate(n_hsa = length(unique(setdiff(hsa,c(NA, "0", "")))),
         n_rno = length(unique(setdiff(rno,c(NA, "0", ""))))) %>% ungroup

gpr.size.limit = 5
gpr.conversion.count = gpr.conversion.result %>% 
  select(rxn_id,n_hsa,cutoff_database,cutoff_rank) %>% distinct %>% 
  left_join(gpr.conversion.result %>% 
              mutate(rno = ifelse(!is.na(rno) & nchar(rno)>0,rno,"0")) %>%
              filter(rno != "0") %>% 
              group_by(rxn_id,n_hsa,cutoff_database,cutoff_rank) %>% 
              summarize(n_rno = n_distinct(rno)) %>% ungroup) %>%
  mutate(n_rno = ifelse(is.na(n_rno),0,n_rno),n_hsa = ifelse(is.na(n_hsa),0,n_hsa),
         n_diff = n_hsa - n_rno,
         facet_database = paste0(cutoff_database, "+ databases"), 
         facet_rank = paste0("Top ", cutoff_rank, " orthologs"),
         rno = paste0(pmin(n_rno, gpr.size.limit), ifelse(n_rno >= gpr.size.limit, "+", "")),
         hsa = paste0(pmin(n_hsa, gpr.size.limit), ifelse(n_hsa >= gpr.size.limit, "+", "")),
         rxn_organism = ifelse(n_hsa > 0 & n_rno > 0, "shared",ifelse(
           n_rno > 0, "rat-specific",ifelse(n_hsa > 0, "human-specific", "non-enzymatic"))),
         n_bias = ifelse(n_rno > 0 & n_hsa > 0, sign(n_diff), 0))

gpr.conversion.summary = gpr.conversion.count %>% 
  filter(n_rno > 0,n_hsa > 0) %>% 
  select(facet_database,facet_rank,rxn_id,n_hsa,n_rno) %>% distinct %>%
  group_by(facet_database,facet_rank) %>% 
  summarize(mean_gpr_ratio = mean(n_rno / n_hsa),
            mean_gpr_difference = mean(n_rno - n_hsa)) %>%
  ungroup

# The first row of this solution was considered the optimal result:
gpr.conversion.summary %>% 
  arrange(abs(mean_gpr_difference) + 
            abs(log2(mean_gpr_ratio)))

# The optimal cutoff was identified as 2+ databases and up to 2 orthologs
gpr.conversion.threshold = gpr.conversion.count %>% 
  mutate(shared = (n_rno > 0) & (n_hsa > 0)) %>% 
  group_by(facet_database,facet_rank) %>% 
  mutate(gpr_shared = sum(shared),
         gpr_ratio = mean(n_rno[shared] / n_hsa[shared]),
         gpr_diff = mean(n_rno[shared] - n_hsa[shared]),
         gpr_human_specific = sum(n_rno == 0 & n_hsa > 0),
         threshold_label = paste0("",round(gpr_ratio,2)),
         threshold_label_x = ifelse(gpr_ratio < 1,"1","4"),
         threshold_label_y = ifelse(gpr_ratio < 1,"5+","0")) %>% ungroup %>% 
  arrange(rxn_id,cutoff_database, cutoff_rank) %>% 
  mutate(cutoff_optimal = cutoff_database == 2 & cutoff_rank == 2)

gpr.color = c(`non-enzymatic` = "#616062", shared = "#674EA7", `rat-specific` = "#CC0000", `human-specific` = "#3C78D8",
              rat = "#CC0000", human = "#3C78D8", rno = "#CC0000", hsa = "#3C78D8", none = "#616062", both = "#674EA7",
              similar = "#674EA7", opposite = "#f79646", ctl = "pink", trt = "green")

# Create plot labels for different thresholds
gpr.threshold.label = gpr.conversion.threshold %>% 
  filter(cutoff_rank <= 5) %>% 
  select(facet_database,
         facet_rank,
         hsa = threshold_label_x,
         rno = threshold_label_y,
         label = threshold_label) %>% 
  distinct
gpr.threshold.label
gpr.size.optimal = gpr.conversion.threshold %>% 
  filter(cutoff_database == 2, cutoff_rank == 2)


gpr.conversion.plot = gpr.conversion.threshold %>% 
  filter(cutoff_rank <= 5) %>%
  ggplot(aes(x = factor(hsa),y = factor(rno))) + 
  geom_point(alpha = 0.3,aes(color = rxn_organism),size = 1,position = "jitter") +
  geom_abline(linetype = "dashed",alpha= 0.5) +
  geom_text(data = gpr.threshold.label,aes(label = label),size = 2.5) +
  scale_color_manual(values = gpr.color) +
  theme_bw(base_size = 12) + 
  theme(legend.position = "none",strip.background = element_rect(fill = 'white')) + 
  xlab("Human GPR size") + ylab("Rat GPR size") +
  facet_grid(facet_rank~facet_database,as.table = F)
gpr.conversion.plot
ggsave(filename = "Reproduce_Appendix_Figure_S1.pdf",
       plot = gpr.conversion.plot,width = 6.7, height = 6.7)

# orthologous gene pairs present in at least 2 of the 5 orthology databases were considered for conversion
# up to 2 rat orthologs were allowed to be converted from each human gene

##################################

# See ef_helper_functions for plotting parameters used to recreate figure 3
gpr.all.orthology.plot = gpr.all.orthology %>% 
  select(rxn_id,n_hsa,n_rno) %>% distinct %>%
  ef_reproduce_figure3(9)
gpr.optimal.plot = gpr.size.optimal %>%
  ef_reproduce_figure3(9)
gpr.curated.plot = gpr.size.curated %>% distinct %>% 
  ef_reproduce_figure3(9)

gpr.all.orthology.plot
gpr.optimal.plot
gpr.curated.plot

ggsave(filename = "Reproduce_Figure3A.pdf",plot = gpr.optimal.plot,width = 2.2, height = 2.2)
ggsave(filename = "Reproduce_Figure3B.pdf",plot = gpr.all.orthology.plot,width = 2.2, height = 2.2)
ggsave(filename = "Reproduce_Figure3C.pdf",plot = gpr.curated.plot,width = 2.2, height = 2.2)

figure3.data = rbind_all(list(
  gpr.size.all %>% select(rxn_id,n_rno,n_hsa) %>% mutate(figure = "3A"),
  gpr.size.optimal %>% select(rxn_id,n_rno,n_hsa) %>% mutate(figure = "3B"),
  gpr.size.curated %>% select(rxn_id,n_rno,n_hsa) %>% mutate(figure = "3C"))) %>% 
  mutate(rno = pmin(n_rno,9),hsa = pmin(n_hsa,9),
         shared = rno > 0 & hsa > 0,
         organism = ifelse(shared,"shared",ifelse(
           rno > 0,"rat-specific",ifelse(hsa > 0,"human-specific","non-enzymatic"))),
         color = gpr.color[organism]) %>%
  group_by(figure) %>% mutate(gpr_ratio = mean(n_rno[shared] / n_hsa[shared])) %>% ungroup

figure3.data %>% filter(shared) %>% count(figure,gpr_ratio)
figure3.data %>% write.table("2015_08_03_sourceDataForFigure3.txt",sep = "\t",quote = F,row.names = F)

figure2a.plot = figure3.data %>% filter(figure == "3C") %>% ef_reproduce_figure2
figure2a.plot
ggsave(filename = "Reproduce_Figure2A.pdf",plot = figure2a.plot,width = 2.2, height = 2.2)
figure2a.data = figure2a.plot$data %>% select(rno,hsa,reaction_density = n_gpr, organism,color) 
figure2a.data %>% write.table("2015_08_03_sourceDataForFigure2A.txt",sep = "\t",quote = F,row.names = F)

# sessionInfo()
# R version 3.2.0 (2015-04-16)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 8 x64 (build 9200)
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_1.0.1  reshape2_1.4.1 dplyr_0.4.1   
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_0.11.6      digest_0.6.8     assertthat_0.1   MASS_7.3-40      grid_3.2.0       plyr_1.8.2       gtable_0.1.2    
# [8] DBI_0.3.1        magrittr_1.5     scales_0.2.4     stringi_0.4-1    lazyeval_0.1.10  proto_0.3-10     tools_3.2.0     
# [15] stringr_1.0.0    munsell_0.4.2    parallel_3.2.0   colorspace_1.2-6