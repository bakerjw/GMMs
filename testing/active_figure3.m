% Script to create a plot of all active models
% Created by Emily Mongold, 1/25/21
%
clear rup
clear site
clc
addpath('../gmms/')
addpath('../')

% call a helper script to build arrays of GMM names and plotting parameters
specify_gmms

%% Setting rupture and site object values
% Rupture inputs
M = 8;
R = [];
Rrup = 25;
Rjb = 25;
Rx = 25;
Ry0 = [];
HW = 0;             % No hanging wall
AS = 0;             % Not an aftershock
Ztor = 0;
Zhyp = 8;
Rhyp = sqrt(Rrup^2+Zhyp^2);
h_eff = [];
W = 10;
delta = 90;
lambda = 0;     % Strike-slip

% Site inputs
is_soil = 1;    % Soft rock
Vs30 = 760;
fvs30 = 0;      % Inferred
Z25 = 0.8;
Z10 = 0.05;
Zbot = 15;
region = 0;     % Global, no corrections applied

%% Create rupture and site objects, other inputs
rup =       rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
site =      site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);
D = 25;
% A1100 (might be optional for CB2008)
MS = 0;         % Does not use magnitude-squared term
sub_ind = [];   % Using for an active event, not subduction

T = 1000;
%% Function Calls
medvec = cell(1,length(gmm_name)-1);
periodvec = cell(1,length(gmm_name)-1);

% a_2015 has maximum magnitude of 6 and was an outlier in the plot
for i = 2:length(gmm_name)
    [medvec{i},~,periodvec{i}] = active_gmms(T,rupt,sitevar,gmm_name{i});
end


%% Create Figure
limits = [0.01 10 .99e-5 3e-2];

figure('Name','Active Figure 3')
for i = 2:length(gmm_name)
    loglog(periodvec{i},medvec{i},line_style{i},'color',line_color{i},'Linewidth',1)
    hold on
end
grid on
axis(limits)
xlabel('Period [sec]')
ylabel('SA [g]')
title("SA versus T")
legend(gmm_name{2:end},'Location','Southwest','Interpreter','none', 'fontsize', 7)

% set figure size
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [4.5 4]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 4.5 4]);
    
%% Save Figure
saveas(gcf,'../figures/Active Figure 3.pdf')
