
% This script integrates gene expression data 
% from the Japanese Toxicogenomics database 
% into human and rat metabolic network models
% This script utilizes the COBRA and TIGER toolboxes
% and an academic version of GUROBI 6.04 (license required).

%% Setup GUROBI solver
addpath(genpath('C:/gurobi604/win64/matlab'));
% Add a path to the GUROBI license file (gurobi.lic):
addpath('C:/Users/Edik/');

%% Initialize COBRA toolbox
% COBRA toolbox is necessary to load models
addpath(genpath('C:/Users/Edik/Documents/MATLAB/cobra/'));
initCobraToolbox % this only needs to be run once
changeCobraSolver('gurobi5');
%% Initialize TIGER toolbox
% TIGER toolbox is necessary to integrate gene expression data
addpath(genpath('C:/Users/Edik/Documents/MATLAB/tiger/'));
start_tiger
set_solver('gurobi');
set_solver_option('MaxTime',60*60);
set_solver_option('IntFeasTol',1e-8);

%% Load human-specific and rat-specific metabolic models
% rno = rattus norvegicus
rno_cobra_load = xls2model('iRnoCOBRA.xlsx'); 
% hsa = homo sapiens
hsa_cobra_load = xls2model('iHsaCOBRA.xlsx'); 

% Each model contains the same sets reactions initially.
% Reactions conisdered "off" have lb and ub equal to 0
off_in_rno = rno_cobra_load.rxns(...
    rno_cobra_load.lb == 0 & rno_cobra_load.ub == 0);
off_in_hsa = hsa_cobra_load.rxns(...
    hsa_cobra_load.lb == 0 & hsa_cobra_load.ub == 0);

% Note that some reactions present in HMR2 were 
% disabled in both models by default:
off_in_both = intersect(off_in_rno,off_in_hsa);

% Remove reactions disabled in each model and set the objective to biomass
rno_cobra = changeObjective(removeRxns(...
    rno_cobra_load,off_in_rno),'RCR99999');
hsa_cobra = changeObjective(removeRxns(...
    hsa_cobra_load,off_in_hsa),'RCR99999');

% Save these COBRA models since xls2model takes a while to import a model
save cobra_models.mat rno_cobra hsa_cobra;

% Convert each COBRA model to a TIGER model
rno_tiger = cobra_to_tiger(rno_cobra,'fast_gpr',true,'status',true);
hsa_tiger = cobra_to_tiger(hsa_cobra,'fast_gpr',true,'status',true);

% Save these TIGER models since conversion takes a while
save tiger_models.mat rno_tiger hsa_tiger;

%% Load saved models
load cobra_models.mat;
load tiger_models.mat;
%% Test each model's ability to produce biomass
% Because exchange fluxes were formulated as fmol / cell / hour and
% biomass was formulated in units of fmol / cell, 
% the growth rate is measured as doublings per hour (1 / hour):

% Values should be the same for cobra and tiger models of the same organism
optimizeCbModel(rno_cobra) % 0.0480 doublings / hour
fba(rno_tiger) % 0.0480 doublings / hour
optimizeCbModel(hsa_cobra) % 0.0403 doublings / hour
fba(hsa_tiger) % 0.0403 doublings / hour

%% Specify a minimum growth rate (biomass flux) for each model
rno_obj_flux_result = fba(rno_tiger);
rno_obj_flux = rno_obj_flux_result.val;
hsa_obj_flux_result = fba(hsa_tiger);
hsa_obj_flux = hsa_obj_flux_result.val;
rno_obj_flux * 24 % about 1.1514 doubling per day
hsa_obj_flux * 24 % about 0.9679 doubling per day

% We want both models to be able to double at least once per week:
obj_value_desired = 1 / 24 / 7;
rno_obj_frac = obj_value_desired / rno_obj_flux;
hsa_obj_frac = obj_value_desired / hsa_obj_flux;

%% Load gene expression changes from the Japanese toxicogenomics database
% Column 1 contains gene identifiers
% Columns 2:120 contain gene expression fold changes or p-values
rno_logfc = readtable('tiger_rno_logfc.txt','Delimiter','\t');
rno_fdr = readtable('tiger_rno_fdr.txt','Delimiter','\t');
hsa_logfc = readtable('tiger_hsa_logfc.txt','Delimiter','\t');
hsa_fdr = readtable('tiger_hsa_fdr.txt','Delimiter','\t');

%% Identify usable datasets from the Japanese toxicogenomics database
% Only integrate expression data for compounds that induced 
% significant (false-discovery rate-adjusted p-value less than 0.1) 
% expression changes in both rat and human hepatocytes
fdr_threshold = 0.1;
rat_index_ok = find(sum(table2array(...
    rno_fdr(:,2:120)) < fdr_threshold, 1) > 0);
human_index_ok = find(sum(table2array(...
    hsa_fdr(:,2:120)) < fdr_threshold, 1) > 0);
index_ok = intersect(human_index_ok,rat_index_ok); 
% 88 compounds identified

%% Integrate gene expression changes across 119 potential compounds
for iii = 1:119
    if (ismember(iii,index_ok))
        rno_genes_input = table2cell(rno_fdr(:,1));
        rno_fc_input = 2.^table2array(rno_logfc(:,iii+1));
        rno_fdr_input = table2array(rno_fdr(:,iii+1));
        rno_genes_ok = (rno_fdr_input < fdr_threshold) & ...
            ismember(rno_genes_input, rno_tiger.varnames);
        % Only pass significantly altered genes into MADE
        % otherwise MADE will sometimes disable genes in both conditions
        sum(rno_genes_ok) 
        rno_made_result = made(rno_tiger,...
            rno_fc_input(rno_genes_ok,1),...
            rno_fdr_input(rno_genes_ok,1),...
            'gene_names',rno_genes_input(rno_genes_ok,1),...
            'obj_frac',rno_obj_frac,...
            'p_thresh',fdr_threshold,...
            'log_fold_change',false,...
            'set_IntFeasTol',1e-8,...
            'round_states',false);
        
        % Export gene state changes as a table
        % NOTE: round_states was set to FALSE so 
        % numbers can be small but non-zero
        rno_made_tbl = table(rno_made_result.genes,...
            rno_made_result.gene_states,...
            rno_made_result.opt_states,...
            repmat(rno_made_result.adj_vals(1),length(rno_made_result.genes),1),...
            repmat(rno_made_result.adj_vals(2),length(rno_made_result.genes),1));
        writetable(rno_made_tbl,['tiger_made_rno_gene_cpd',num2str(iii),'.txt'],'Delimiter','\t');
    end
    if (ismember(iii,index_ok))
        hsa_genes_input = table2cell(hsa_fdr(:,1));
        hsa_fc_input = 2.^table2array(hsa_logfc(:,iii+1));
        hsa_fdr_input = table2array(hsa_fdr(:,iii+1));
        hsa_genes_ok = (hsa_fdr_input < fdr_threshold) & ...
            ismember(hsa_genes_input, hsa_tiger.varnames);
        % Only pass significantly altered genes into MADE
        % otherwise MADE will sometimes disable genes in both conditions
       sum(hsa_genes_ok)
        hsa_made_result = made(hsa_tiger,...
            hsa_fc_input(hsa_genes_ok,1),...
            hsa_fdr_input(hsa_genes_ok,1),...
            'gene_names',hsa_genes_input(hsa_genes_ok,1),...
            'obj_frac',hsa_obj_frac,...
            'p_thresh',fdr_threshold,...
            'log_fold_change',false,...
            'set_IntFeasTol',1e-8,...
            'round_states',false);
        
        % Export gene state changes as a table
        % NOTE: round_states was set to FALSE so 
        % numbers can be small but non-zero
        hsa_made_tbl = table(hsa_made_result.genes,...
            hsa_made_result.gene_states,...
            hsa_made_result.opt_states,...
            repmat(hsa_made_result.adj_vals(1),length(hsa_made_result.genes),1),...
            repmat(hsa_made_result.adj_vals(2),length(hsa_made_result.genes),1));
        writetable(hsa_made_tbl,['tiger_made_hsa_gene_cpd',num2str(iii),'.txt'],'Delimiter','\t');
    end
end

%% Gene states were applied to SYBIL models of iHsa and iRno within R.
% end
