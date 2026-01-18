% run_pr3_f5_f10_before.m
% 
% This script should reside in:
% opensha_checks/gmms/as_2008_active/matlab_runners/
%
% Run version of MATLAB function before PR3 fix; use Git to point
% to the specific version before running MATLAB.
% Load IMTs and PR3 scenario, run MATLAB function, and write results.
clear; clc;

%% Parameters
PR_str = 'pr3'; % Which PR-specific input test files?

%% Paths
repo_root = '../../../../';
dir_implementation = fullfile(repo_root, 'gmms');
dir_tests = fullfile(repo_root, 'opensha_checks', 'gmms', 'as_2008_active');
dir_output = fullfile(dir_tests, 'outputs_matlab', 'pr3_f5_f10_before');
if ~exist(dir_output, 'dir'); mkdir(dir_output); end

file_IMTs = fullfile(dir_tests, 'inputs', 'IMTs', 'as2008_IMTs.csv');
file_scenarios = fullfile(dir_tests, 'inputs', 'scenarios', 'test_scenarios.csv');
file_output = fullfile(dir_output, [PR_str '_results.csv']);

%% Import specific version of AS2008 (after using Git)
addpath(dir_implementation);

%% Load IMTs
% Read CSV file of IMTs
opts = detectImportOptions(file_IMTs);
opts = setvartype(opts, 'IMT', 'string'); % Force IMT column to be text
df_IMTs = readtable(file_IMTs, opts);   

% Get IMTs for joining later
IMT = df_IMTs.IMT;

% Get row vector of IMTs (0 = PGA, -1 = PGV, others SA(T))
T = df_IMTs.T_equiv';

%% Load scenarios
df_scenarios = readtable(file_scenarios);

%% Identify input test scenario of interest
s = df_scenarios(strcmp(df_scenarios.scenario_id, PR_str), :);

%% Run MATLAB function for given scenario
% Build rup structure
rup.M    = s.M;
rup.Rrup = s.Rrup;
rup.Rjb  = s.Rjb;
rup.Rx   = s.Rx;
rup.delta  = s.dip;
rup.lambda = s.rake;
rup.Ztor = s.Ztor;
rup.W    = s.W;
rup.AS   = s.AS;
rup.HW   = s.HW;

% Build site structure
site.Vs30  = s.Vs30;
site.fvs30 = s.fvs30;
site.Z10  = s.Z1P0;

% Call MATLAB implementation (after using Git)
[median_vals, sigma_vals] = as_2008_active(T, rup, site);

%% Convert MATLAB output to column vectors and summarize
T_out = T';
med_out = median_vals';
sig_out = sigma_vals';
df_out = table(IMT, T_out, med_out, sig_out, ...
'VariableNames', {'IMT', 'T_equiv', 'median', 'sigma'});

%% Save results
% Before saving, ensure correct version of MATLAB function was used
% writetable(df_out, file_output);

%% Echo end
fprintf('Finished running specific version of MATLAB function.\n');