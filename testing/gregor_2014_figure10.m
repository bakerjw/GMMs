% Script to recreate Figure 10 from Gregor et. al. (2014)
% Using the following five GMMs
%   ask_2014_active
%   bssa_2014_active
%   cb_2014_active
%   cy_2014_active
%   i_2014_active
% Created by Emily Mongold, 12/18/20
%
clear
clc; addpath('../gmms/')

% call a helper script to build arrays of GMM names and plotting parameters
specify_gmms
%% Setting rupture and site object values
% Rupture inputs
M = [5 7];           
R = [];     
Rrup = 30;      
Rjb = 30;        
Rhyp = [];      
Rx = 30;        
Ry0 = [];       
HW = 0;             % Hanging wall indicator = 1 for hanging wall, 0 otherwise
AS = 0;             % Aftershock indicator = 1 for aftershock, 0 otherwise
Ztor = 1;           % Different for each magnitude earthquake
Zhyp = 8;
h_eff = [];
W = [];
delta = 15;
lambda = 0;         % Strike-slip

% Site inputs
is_soil = [];       % Not used in these models
Vs30 = 760;
fvs30 = 0;          % 0 for inferred and 1 for measured
Z25 = 1.9826;       % Used in CB model 
Zbot = 10;      
region = 0;         % = 0 for global

% Create rupture and site objects
rupt = rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
sitevar = site(is_soil,Vs30,fvs30,Z25,[],Zbot,region);

T_vec = 0.01:0.01:10;  % Independent Variable
%% Function Calls
sigmaASK2 = zeros(2,length(T_vec));
sigmaBSSA2 = zeros(2,length(T_vec));
sigmaCB2 = zeros(2,length(T_vec));
sigmaCY2 = zeros(2,length(T_vec));
sigmaI2 = zeros(2,length(T_vec));
for j = 1:2
    rupt.M = M(j); 
    for n = 1:length(T_vec)
        [~,sigmaASK2(j,n),~] = active_gmms(T_vec(n),rupt,sitevar,'ask_2014');
        [~,sigmaBSSA2(j,n),~] = active_gmms(T_vec(n),rupt,sitevar,'bssa_2014');
        [~,sigmaCB2(j,n),~] = active_gmms(T_vec(n),rupt,sitevar,'cb_2014');
        [~,sigmaCY2(j,n),~] = active_gmms(T_vec(n),rupt,sitevar,'cy_2014');
        [~,sigmaI2(j,n),~] = active_gmms(T_vec(n),rupt,sitevar,'i_2014');
    end
end
%% Figure 10
limits = [0.01 10 0 1];
titles = ["Mag = 5, Rrup = 30 km";"Mag = 7, Rrup = 30 km"];
figure('Name','Gregor Figure 10','NumberTitle','off','Position',[10 10 600 400])
for n = 1:2
    subplot(1,2,n)
    semilogx(T_vec, sigmaASK2(n,:),'-r',T_vec, sigmaBSSA2(n,:),'-g',T_vec, sigmaCB2(n,:),'-b',T_vec, sigmaCY2(n,:),'-m',T_vec, sigmaI2(n,:),'-c','LineWidth',1)
    grid on 
    axis(limits)
    xlabel('Period [sec]')
    ylabel('Sigma [Ln Units]')
    title(titles(n))
    if n == 2
        legend(gmm_name{nga2west},'Location','Southeast','Interpreter','none')
    end
end

% set figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [5 3]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 5 3]);

%% Save Figure
saveas(gcf,'../figures/gregor10.pdf')