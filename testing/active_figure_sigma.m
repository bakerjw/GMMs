% Script to create plot of standard deviation versus period
% Based on Figure 10 from Gregor et. al. (2014)
%
% Created by Emily Mongold, 1/24/21
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
M = 6;           
R = [];     
Rrup = 30;      
Rjb = 30;        
Rhyp = 50;      
Rx = 30;        
Ry0 = [];       
HW = 0;             % No hanging wall
AS = 0;             % Not an aftershock
Ztor = 3;       
Zhyp = 8;
h_eff = [];
W = 10;
delta = 15;
lambda = 0;         % Strike-slip

% Site inputs
is_soil = 1;
Vs30 = 760;
fvs30 = 0;          % Inferred
Z25 = 2.0;       
Z10 = 0.5;          
Zbot = 10;      
region = 0;         % Global

% Create rupture and site objects
rup = rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
site = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

T = 1000;
%% Function Calls
sigmavec = cell(1,length(gmm_name));
periodvec = cell(1,length(gmm_name));

for i = 1:length(gmm_name)
    [~,sigmavec{i},periodvec{i}] = active_gmms(T,rup,site,gmm_name{i});
end
%% Create Figure
limits = [0.01 10 0 1];
figure('Name',"St.Dev.",'NumberTitle','off')
for i = 1:length(gmm_name)
    semilogx(periodvec{i},sigmavec{i},line_style{i},'color',line_color{i},'LineWidth',1);
    hold on
end
grid on 
axis(limits)
xlabel('Period [sec]')
ylabel('Sigma [Ln Units]')
title("Sigma versus T")
legend(gmm_name,'Location','SouthEast','Interpreter','none', 'fontsize', 7)
hold off


% set figure size
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [5 5]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 5 5]);
    
%% Save Figure
saveas(gcf,'../figures/Active Figure sigma.pdf')