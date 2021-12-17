% compare_script.m to run the function compare_general for various models,
% comparing the outputs of OpenQuake to the results from the corresponding
% GMM Matlab function.
% Created by Emily Mongold
% 19 August 2021

clear;clc;
addpath ../gmms
%%

OQ_file = "IDRISS_2014_MEAN.csv";
gmm = 'i_2014';
[max_err,error_test,periods] = compare_general(OQ_file,gmm);
OQ_table = readtable(OQ_file,"ReadVariableNames",true);

%% 

OQ_file = "ATKINSON2015_MEAN.csv";
gmm = 'a_2015';
[max_err,error_test,periods] = compare_general(OQ_file,gmm);
OQ_table = readtable(OQ_file,"ReadVariableNames",true);

%% 
% Note that at the time of upload, the OpenQuake output file is
% inconsistent with the GMPE Spreadsheet and GMM Matlab function.
OQ_file = "ASK14_ResMEAN_RegCAL.csv";
gmm = 'ask_2014';
[max_err,error_test,periods] = compare_general(OQ_file,gmm);
OQ_table = readtable(OQ_file,"ReadVariableNames",true);

%% 

OQ_file = "CB2014_MEAN.csv";
gmm = 'cb_2014';
[max_err,error_test,periods] = compare_general(OQ_file,gmm);
OQ_table = readtable(OQ_file,"ReadVariableNames",true);

%% 

OQ_file = "AS97_MEAN_SS.csv";
gmm = 'as_1997';
[max_err,error_test,periods] = compare_general(OQ_file,gmm);
OQ_table = readtable(OQ_file,"ReadVariableNames",true);

%%

OQ_file = "BSSA_2014_CALIFORNIA_MEAN.csv";
gmm = 'bssa_2014';
[max_err,error_test,periods] = compare_general(OQ_file,gmm);
OQ_table = readtable(OQ_file,"ReadVariableNames",true);

%%
% To make an output csv file with absolute errors
[~,inda] = ismember('pga',OQ_table.Properties.VariableNames);
if ismember('pgv',OQ_table.Properties.VariableNames)
    [~,indv] = ismember('pgv',OQ_table.Properties.VariableNames);
else
    indv = 0; 
end
if inda~= 0 && indv ~= 0
    ind = min(inda,indv);
elseif indv == 0
    ind = inda; 
end
full_out = OQ_table;
y = ones(1,length(error_test)); 
x = ones(1,width(error_test));
full_out(:,ind:end) = mat2cell(abs(error_test),y,x);
filename = append("..\",gmm,"_output.csv");
writetable(full_out,filename);
clear x y ind inda indv
