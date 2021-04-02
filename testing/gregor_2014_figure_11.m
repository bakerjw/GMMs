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
Ztor = 1;   % Different for each magnitude earthquake
Zhyp = 8;
h_eff = [];
W = [];
delta = 15;
lambda = 0;     % Strike-slip

% Site inputs
is_soil = [];
% soil type     = 0 for soil
%               = 1 for soft rock
%               = 2 for hard rock
Vs30 = 760;
fvs30 = 0;      % 0 for inferred and 1 for measured
Z25 = 1.9826;       % Used in CB model 
Z10 = 0;            % Used in BSSA model
Z10_cy = 0.4794;     
Z10_ask = 0.4704;
Zbot = 10;      
region = 0;    % = 0 for global
%                = 1 for California
%                = 2 for Japan
%                = 3 for China
%                = 4 for Italy
%                = 5 for Turkey
%                = 6 for Taiwan

% Create rupture and site objects
rup_cy = rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
rup = rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
site_cy = site(is_soil,Vs30,fvs30,Z25,Z10_cy,Zbot,region);
site_ask = site(is_soil,Vs30,fvs30,Z25,Z10_ask,Zbot,region);
site = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

T_vec = [0 1];      % PGA or 1.0 sec
%% Function Calls
sigmaASK = zeros(2,length(M));
sigmaBSSA = zeros(2,length(M));
sigmaCB = zeros(2,length(M));
sigmaCY = zeros(2,length(M));
for j = 1:2
    for n = 1:length(M)
        rup.M = M(n); rup_cy.M = M(n);
        [~,sigmaASK(j,n),~] = ask_2014_active(T_vec(j),rup,site_ask);
        [~,sigmaBSSA(j,n),~] = bssa_2014_active(T_vec(j),rup,site);
        [~,sigmaCB(j,n),~] = cb_2014_active(T_vec(j),rup,site);
        [~,sigmaCY(j,n),~] = cy_2014_active(T_vec(j),rup_cy,site_cy);
        if M(n)>= 5
            [~,sigmaI(j,n-4),~] = i_2014_active(T_vec(j),rup,site);
        end
    end
end
%% Figure 11
limits = [3 8.5 0 1];
titles = ["PGA";"T = 1.0 sec"];
figure('Name','Gregor Figure 11','NumberTitle','off','Position',[10 10 600 400])
for n = 1:2
    subplot(1,2,n)
    plot(M,sigmaASK(n,:),'-r',M, sigmaBSSA(n,:),'-g',M, sigmaCB(n,:),'-b',M, sigmaCY(n,:),'-m',M(5:12),sigmaI(n,:),'-c','LineWidth',1)
    grid on 
    grid minor
    axis(limits)
    xlabel('Magnitude [g]')
    ylabel('Sigma [Ln Units]')
    title(titles(n))
    if n == 2
        legend('ASK','BSSA','CB','CY','I','Location','Southeast')
    end
end
%% Save Figure
saveas(gcf,'../figures/gregor11.jpg')