% compare_pr2.m
% 
% This script should reside in:
% opensha_checks/gmms/as_2008_active/matlab_runners/
%
% Compare medians from:
%   (1) OpenSHA expected
%   (2) MATLAB after PR2 fix
%
% Print comparison table and max abs diff. Plot SA medians.
clear; clc;

%% Paths
repo_root = '../../../../';
dir_tests = fullfile(repo_root, 'opensha_checks', 'gmms', 'as_2008_active');

file_opensha = fullfile(dir_tests, 'outputs_opensha', 'pr2_opensha_expected.csv');
file_matlab = fullfile(dir_tests, 'outputs_matlab', 'pr2_f10_after', 'pr2_results.csv');

%% Load table: OpenSHA
opts_opensha = detectImportOptions(file_opensha);
opts_opensha = setvartype(opts_opensha, 'IMT', 'string'); % Force IMT column to be text
df_opensha = readtable(file_opensha, opts_opensha);   

%% Load table: MATLAB
opts_matlab = detectImportOptions(file_matlab);
opts_matlab = setvartype(opts_matlab, 'IMT', 'string'); % Force IMT column to be text
df_matlab = readtable(file_matlab, opts_matlab);   

%% Subset and rename columns
% Subset
df_opensha_sub = df_opensha(:, {'IMT', 'T_equiv', 'median'});
df_matlab_sub = df_matlab(:, {'IMT', 'T_equiv', 'median'});

% Rename
df_opensha_sub.Properties.VariableNames{'median'} = 'median_opensha';
df_matlab_sub.Properties.VariableNames{'median'} = 'median_matlab';


%% Comparison table
df_out = join(df_opensha_sub, df_matlab_sub, ...
    'Keys',{'IMT', 'T_equiv'});

%% Compute absolute percent differences
df_out.perc_diff = (df_out.median_matlab - df_out.median_opensha) ./ df_out.median_opensha * 100;
df_out.abs_perc_diff = abs(df_out.perc_diff);

% Summarize
fprintf('Maximum absolute percent difference in medians: %.2g percent.\n', max(df_out.abs_perc_diff));

% Inspect df_out for non-SA IMTs
open df_out

%% Plot only SAs
% Get subset
isSA = df_out.T_equiv > 0;
df_SA = df_out(isSA, :);


% Plot
figure; hold on; grid on; box on;

plot(df_SA.T_equiv, df_SA.median_opensha, 'bo-', 'LineWidth', 1.5, ...
    'DisplayName', 'OpenSHA');
plot(df_SA.T_equiv, df_SA.median_matlab, 'rx--', 'LineWidth', 1.5, ...
    'DisplayName', 'MATLAB (after PR2 fix)');

xlabel('Period T (s)');
ylabel('Median SA (g)');
title('PR2 Comparison: SA Medians (OpenSHA vs MATLAB)');
legend('Location','best');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

hold off;
