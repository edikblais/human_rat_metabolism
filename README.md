# human_rat_metabolism
Reconciled reconstructions of rat and human metabolism. Scripts used during the reconstruction of iRno and iHsa are included.

The first script (in R), compares the complexity of gene-protein-reaction (GPR) directly between rat and human metabolic models as described in Figures 2 and 3. This script inputs gene-protein-reaction (GPR) relationship information for various versions of iRno and iHsa to compare the global distribution of GPR sizes between the rat and human models. The goal of this script was to demonstrate that  GPR sizes were consistent across species after prioritizing orthologs during the optimized conversion of iHsa to iRno as well as after manual curation.

The second script (in R), validates knowns differences between rats and humans related to vitamin C synthesis as described in Figure 4. First, this script reproduces the known result that rats can produce vitamin C de novo and that humans cannot. This script also reproduces the prediction that humans rely on dietary sources of vitamin C under physiological conditions, consistent with the known effects of scurvy. Using the sybil toolbox in R, Flux balance analysis was applied to manually curated versions of iRno and iHsa to test growth under various uptake rates for vitamin C.

The third script (in R), validates known differences between rats and human related to bile acid synthesis as described in Figure 5.

The fourth script (in MATLAB), integrates gene expression changes into rat and human metabolic networks to produce toxicology subnetworks that are used to make predictions described in Figures 6 and 7. COBRA toolbox models iRno and iHsa were loaded and converted into TIGER models. Gene expression fold changes and FDR-adjusted p-values were integrated into each model independently to generate subnetworks consistent with toxicant-treated and control-treated experimental conditions. 119 toxicants were initially examined and 70 that induced significant metabolic expression changes in both rat and human hepatocytes were selected for further analysis.

The fifth script (in R), uses models generated in the fourth script to generate reaction-state predictions that were used to make predictions described in Figures 6 and 7.

The sixth script (in MATLAB) validates rat and human metabolic models with the ability to perform biologically relevant tasks. New tasks to validate species-specific differences and non-specific functions were added to a compendium of tasks that were adapted from HMR2 and Recon2.


