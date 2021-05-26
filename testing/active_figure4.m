% Script to create plot of standard deviation versus period
% Based on Figure 10 from Gregor et. al. (2014)
%
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
M = 5;           
R = [];     
Rrup = 100;      
Rjb = 100;        
Rhyp = 150;      
Rx = 100;        
Ry0 = [];       
HW = 0;             % No hanging wall
AS = 0;             % Not an aftershock
Ztor = 6;       
Zhyp = 8;
h_eff = [];
W = 20;
delta = 90;
lambda = 0;         % Strike-slip

% Site inputs
is_soil = 0;
Vs30 = 760;
fvs30 = 0;          % Inferred
Z25 = 2.0;       
Z10 = 0.5;          
Zbot = 40;      
region = 0;         % Global

% Create rupture and site objects
rup = rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
site = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

T = 1000; % For all available periods
%% Function Calls
medianvec = cell(1,length(gmm_name));
periodvec = cell(1,length(gmm_name));
for i = 1:length(gmm_name)
    [medianvec{i},~,periodvec{i}] = active_gmms(T,rup,site,gmm_name{i});
end
%% Create Figure
limits = [0.01 10 0.0001 0.1];
figure('Name','Active Figure 4','NumberTitle','off')
for i = 1:length(gmm_name)
    loglog(periodvec{i},medianvec{i},line_style{i},'color',line_color{i},'LineWidth',1)
    hold on
end
grid on 
axis(limits)
xlabel('Period [sec]')
ylabel('SA [g]')
title("SA versus T")
legend(gmm_name,'Location','Northeast','Interpreter','none', 'fontsize', 7)



% set figure size
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [5 5]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 5 5]);
    
%% Save Figure
saveas(gcf,'../figures/Active Figure 4.pdf')