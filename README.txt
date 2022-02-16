IOMA - Integrative Omics-Metabolic Analysis (DOI: 10.1093/bioinformatics/btq183)

The input to the method include the following fields:

% model - metabolic model in standard format
% exp_mat - expression matrix in which row j and column i has [e^j_i]/[e^ref_i]
% met_frw_mat - saturation values matrix in which row j and column i has a+^j_i
% met_bck_mat - saturation values matrix in which row j and column i has a-^j_i
% rxns_idx - A vector with  reactions indices corresponding to the core set
% used in the problem
% ko_rxns - A matrix where for each row (knockout j), the indices of the
% knocked-out reactions are specify
% ex_flux - exchange fluxes matrix where each row represent a different exchange reaction
% and each column represent a differen knockout
% ex_rxns - A vector with the corresponding exchange reaction indices 
% GR - A vector of growth rates values in each one of the knockouts

The output to the method include the following fields:
% V_v = predicted flux rates