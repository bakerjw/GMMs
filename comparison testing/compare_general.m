function [max_err,error_test,periods] = compare_general(OQ_file,gmm)
%COMPARE_GENERAL test file to compare OpenQuake outputs to those from the
%corresponding Matlab GMM function, for median values
%   Extracts input parameters from the output OpenQuake file and compares
%   values to those from a Matlab function output
%% Comprehensive Testfile Recreation and Comparison
% global OQ_table 
OQ_table = readtable(OQ_file,"ReadVariableNames",true);%,"NumHeaderLines",1);
get_periods  = readtable(OQ_file,"ReadVariableNames",false); % to keep the period values

% Set input periods to be the same as OpenQuake
[pga_flag,inda] = ismember('pga',OQ_table.Properties.VariableNames);
if ismember('pgv',OQ_table.Properties.VariableNames)
    [pgv_flag,indv] = ismember('pgv',OQ_table.Properties.VariableNames);
else
    pgv_flag = 0; 
    indv = 0; 
end
if inda~= 0 && indv ~= 0
    ind = min(inda,indv);
elseif indv == 0
    ind = inda; 
end
OQ_results = OQ_table{:,ind:end}; % Array of results
periods = get_periods{1,ind:end}; % Array of periods, with NaN for pga or pgv columns, if applicable
if pgv_flag == 1 && pga_flag == 1
    if inda < indv
        periods(1:2) = [0 -1]; % These replace NaN with the appropriate input for pga or pgv
    else
        periods(1:2) = [-1 0];
    end
elseif pga_flag == 1
    periods(1) = 0;
end

% Inputs to run GMM
addpath ..\GitHub\gmms
median_GMM = zeros(height(OQ_table),width(OQ_results));
for i = 1:height(OQ_table)
    rup.h_eff = []; % Use default value (overwritten)
    rup.HW = 0; % Default not hanging wall
    rup.AS = 0; % Default not aftershock
    if ismember('rup_mag',OQ_table.Properties.VariableNames)
        rup.M = OQ_table.rup_mag(i);
        rup.h_eff = max(1, 10.^(-1.72+0.43.*rup.M)); % Old equation from a_2015 which OpenQuake still uses
    end
    if ismember('rup_dip',OQ_table.Properties.VariableNames)
        rup.delta = OQ_table.rup_dip(i);
    else
        rup.delta = [];
    end
    if ismember('rup_ztor',OQ_table.Properties.VariableNames)
        rup.Ztor = OQ_table.rup_ztor(i);
    else 
        rup.Ztor = [];
    end
    if ismember('rup_width',OQ_table.Properties.VariableNames)
        rup.W = OQ_table.rup_width(i);
    else 
        rup.W = [];
    end
    if ismember('rup_rake',OQ_table.Properties.VariableNames)
        rup.lambda = OQ_table.rup_rake(i); 
    else
        rup.lambda = 0; % Strike-slip as default
    end
    if ismember('dist_rrup',OQ_table.Properties.VariableNames)
        rup.Rrup = OQ_table.dist_rrup(i);
    else
        rup.Rrup = []; 
    end
    if ismember('dist_rjb',OQ_table.Properties.VariableNames)
        rup.Rjb = OQ_table.dist_rjb(i);
    else 
        rup.Rjb = []; 
    end
    if ismember('dist_rx',OQ_table.Properties.VariableNames)
        rup.Rx = OQ_table.dist_rx(i);
        if rup.delta ~= 90 %&& rup.Rx >0
            rup.HW = 1; % Hanging wall
        end
    else 
        rup.Rx = []; 
    end
    if ismember('dist_ry0',OQ_table.Properties.VariableNames)
        rup.Ry0 = OQ_table.dist_ry0(i);
    else 
        rup.Ry0 = []; 
    end
    if ismember('dist_rhypo',OQ_table.Properties.VariableNames)
        rup.Rhyp = OQ_table.dist_rhypo(i);
    else 
        rup.Rhyp = []; 
    end
    if ismember('rup_hypo_depth',OQ_table.Properties.VariableNames)
        rup.Zhyp = OQ_table.rup_hypo_depth(i);
    else 
        rup.Zhyp = []; 
    end
    if ismember('site_vs30',OQ_table.Properties.VariableNames)
        site.Vs30 = OQ_table.site_vs30(i);
    else 
        site.Vs30 = []; 
    end
    if ismember('site_vs30measured',OQ_table.Properties.VariableNames)
        site.fvs30 = OQ_table.site_vs30measured(i);
    else 
        site.fvs30 = []; 
    end
    if ismember('site_z1pt0',OQ_table.Properties.VariableNames)
        site.Z10 = OQ_table.site_z1pt0(i); 
    else 
        site.Z10 = [];
    end
    if ismember('site_z2pt5',OQ_table.Properties.VariableNames)
        site.Z25 = OQ_table.site_z2pt5(i);
    else 
        site.Z25 = []; 
    end
    site.Zbot = []; % Not provided by OpenQuake
    site.region = 0; % Global as default (can change depending on input- separate .csv files for different regions from OpenQuake)
    site.is_soil =[]; % Not provided by OpenQuake
    

    for j = 1:length(periods)
        T = periods(j);
        [median_GMM(i,j),~,~] = active_gmms(T,rup,site,gmm);
    end
end

error_test = median_GMM - OQ_results;
max_err = max(max(abs(error_test)));

%% Running with values from spreadsheet- correct for ask2014 (uncomment and pause at end to view error matrices)
% inds_ask = [1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	23	24];
% baseline = readmatrix('ask_spreadsheet.csv');
% oq_err = baseline - OQ_results(53:89,inds_ask);
% mat_err = baseline - median_GMM(53:89,inds_ask);
%% Running with values from spreadsheet--6 correct for cb2014 (uncomment and pause at end to view error matrices)
% baseline = readmatrix('cb_spreadsheet.csv');
% oq_err = baseline - OQ_results(53:89,:);
% mat_err = baseline - median_GMM(53:89,:);
end

