
########################################################
#################  Load packages  ######################
########################################################


library(glpkAPI)
library(sybil)
library(dplyr)
library(ggplot2)
library(reshape2)

# The default objective function for performing 
# flux balance analysis is the biomass reaction
maximize_objective = function(model,objective = "RCR99999") {
  stopifnot(require(sybil))
  opt_result = model %>% changeObjFunc(objective = objective) %>% 
    optimizeProb(algorithm = "fba", retOptSol = T)
  opt_result@lp_obj
}
# Import sybil models of iRno and iHsa
irno_model = readTSVmod("2015_07_rno","tsv",
                        balanceReact=F,def_bnd = 10000)
ihsa_model = readTSVmod("2015_07_hsa","tsv",
                        balanceReact=F,def_bnd = 10000)

# Two reactions allow for vitamin C input/output 
# at the extracellular boundary (exchange) and
# output within the cyotsol (demand)
vitaminc.exchange = data_frame(
  rxn_id = "RCR30122", rxn_name = "ascorbate exchange",
  lb = -1, ub = 1000)
vitaminc.demand = data_frame(
  rxn_id = "RCR90119", rxn_name = "vitamin C demand", 
  lb = 0, ub = 1000)

# Both iRno and iHsa can utilize vitamin C 
# under (default) physiological conditions:
irno_model %>% 
  maximize_objective(vitaminc.demand$rxn_id) # 473.3072
ihsa_model %>% 
  maximize_objective(vitaminc.demand$rxn_id) # 1

# Removing vitamin C from physiological conditions 
# blocks vitamin C utilization in iHsa but not iRno
irno_model %>% changeBounds(
  react = vitaminc.exchange$rxn_id,lb = 0,ub = 1000) %>% 
  maximize_objective(vitaminc.demand$rxn_id) # 472.3072
ihsa_model %>% changeBounds(
  react = vitaminc.exchange$rxn_id,lb = 0,ub = 1000) %>% 
  maximize_objective(vitaminc.demand$rxn_id) # 0

##################################################################

# growth was tested while vitamin C uptake was varied between zero 
# and the default physiological value of 1 fmol / cell / hr
vitaminc.range = c(1,5*c(0:6)/100)

vitaminc.biomass = vitaminc.range %>% lapply(function(x,ex) {
  ex$vitaminc = x
  ex$rno = irno_model %>% changeBounds(react = ex$rxn_id,lb = ex$lb * x,ub = ex$ub) %>% 
    maximize_objective # default objective is biomass
  ex$hsa = ihsa_model %>% changeBounds(react = ex$rxn_id,lb = ex$lb * x,ub = ex$ub) %>% 
    maximize_objective
  print(ex) # print progress to console
  ex
},vitaminc.exchange) %>% rbind_all

# the uptake rate of hepatocytes was estimated to be:
# 1 fmol / cell / hr 
# this is consistent with being about 1 order of magnitude above the 
# rate at which vitaminc could be the rate-limiting factor.
# this reveals insights into how robust cells might be to decreased 
# availability of this essential nutrient in normal physiological conditions.

figure4c = vitaminc.biomass %>% filter(vitaminc <= .5) %>% 
  rbind(vitaminc.biomass %>% filter(vitaminc == 1) %>% mutate(vitaminc = 0.35)) %>%
  select(rxn_name,vitaminc,rno,hsa) %>% 
  melt(c("rxn_name","vitaminc")) %>% as.tbl %>%
  group_by(variable) %>% mutate(nvalue = value / max(value)) %>% ungroup %>% 
  ggplot(aes(x = vitaminc, y = value, color = variable,group = variable)) + 
  geom_hline(yintercept = 0,size = 2,alpha = 0.3) + #alpha 0.3,size = 2,
  geom_line(alpha = 0.8,size = 2) + geom_point(size = 4) + 
  geom_vline(xintercept = 0.325,size = 1,alpha = 0.3,linetype = "dashed") + #alpha 0.3,size = 2,
  theme_bw(base_size = 12) + scale_color_manual(values = c(rno = "#CC0000",hsa = "#3C78D8")) + 
  xlab("vitamin C uptake limit\n(fmol / cell / hour)") + 
  ylab("predicted growth rate\n(1 / hour)") +
  scale_x_reverse(breaks = .05*c(0:7),
                  labels = setNames(gsub("0.35","1",.05*c(0:7),fixed = T),.05*c(0:7))) + 
  theme(title = element_text(size = 12),legend.position="none")

figure4c
ggsave(filename = "Reproduce_Figure4c.pdf",plot = figure4c,width = 3.22, height = 2.5)
figure4c.data = figure4c$data %>% select(vitaminc_max_uptake = vitaminc, organism = variable,max_growth = value, normalized_growth = nvalue)
figure4c.data %>% write.table("2015_08_03_sourceDataForFigure4C.txt",sep = "\t",quote = F,row.names = F)

# based on biomass reaction, 
# the estimated concentration of vitamin C 
# in one hepatocyte was about 6 fmol / cell.
