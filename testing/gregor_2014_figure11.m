% Script to recreate Figure 11 from Gregor et. al. (2014)
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
M = 3:0.5:8.5;           
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
is_soil = [];
Vs30 = 760;
fvs30 = 0;          % 0 for inferred and 1 for measured
Z25 = 1.9826;       % Used in CB model 
% Z10 = 0;            % Used in BSSA model
% Z10_cy = 0.4794;     
% Z10_ask = 0.4704;
Zbot = 10;      
region = 0;         % Global


% Create rupture and site objects
rupt = rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
sitevar = site(is_soil,Vs30,fvs30,Z25,[],Zbot,region);

T_vec = [0 1];      % PGA or 1.0 sec
%% Function Calls
sigmaASK2 = zeros(2,length(M));
sigmaBSSA2 = zeros(2,length(M));
sigmaCB2 = zeros(2,length(M));
sigmaCY2 = zeros(2,length(M));
for j = 1:2
    for n = 1:length(M)
        rupt.M = M(n); rup_cy.M = M(n);
        [~,sigmaASK2(j,n),~] = active_gmms(T_vec(j),rupt,sitevar,'ask_2014');
        [~,sigmaBSSA2(j,n),~] = active_gmms(T_vec(j),rupt,sitevar,'bssa_2014');
        [~,sigmaCB2(j,n),~] = active_gmms(T_vec(j),rupt,sitevar,'cb_2014');
        [~,sigmaCY2(j,n),~] = active_gmms(T_vec(j),rupt,sitevar,'cy_2014');
        if M(n)>= 5
            [~,sigmaI2(j,n),~] = active_gmms(T_vec(j),rupt,sitevar,'i_2014');
        end
    end
end
%% Figure 11
limits = [3 8.5 0 1];
titles = ["PGA";"T = 1.0 sec"];
figure('Name','Gregor Figure 11','NumberTitle','off','Position',[10 10 600 400])
for n = 1:2
    subplot(1,2,n)
    plot(M,sigmaASK2(n,:),'-r',M, sigmaBSSA2(n,:),'-g',M, sigmaCB2(n,:),'-b',M, sigmaCY2(n,:),'-m',M(5:12),sigmaI2(n,5:12),'-c','LineWidth',1)
    grid on 
    grid minor
    axis(limits)
    xlabel('Magnitude [g]')
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
saveas(gcf,'../figures/gregor11.pdf')