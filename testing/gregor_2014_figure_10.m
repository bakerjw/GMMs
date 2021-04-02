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
Z10 = 0;            % Used in BSSA model
Z10_cy = 0.4794;     
Z10_ask = 0.4704;
Zbot = 10;      
region = 0;         % = 0 for global

% Create rupture and site objects
rup_cy = rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
rup = rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
site_cy = site(is_soil,Vs30,fvs30,Z25,Z10_cy,Zbot,region);
site_ask = site(is_soil,Vs30,fvs30,Z25,Z10_ask,Zbot,region);
site = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

T_vec = 0.01:0.01:10;  % Independent Variable
%% Function Calls
sigmaASK = zeros(2,length(T_vec));
sigmaBSSA = zeros(2,length(T_vec));
sigmaCB = zeros(2,length(T_vec));
sigmaCY = zeros(2,length(T_vec));
sigmaI = zeros(2,length(T_vec));
for j = 1:2
    rup.M = M(j); rup_cy.M = M(j);
    for n = 1:length(T_vec)
        [~,sigmaASK(j,n),~] = ask_2014_active(T_vec(n),rup,site_ask);
        [~,sigmaBSSA(j,n),~] = bssa_2014_active(T_vec(n),rup,site);
        [~,sigmaCB(j,n),~] = cb_2014_active(T_vec(n),rup,site);
        [~,sigmaCY(j,n),~] = cy_2014_active(T_vec(n),rup_cy,site_cy);
        [~,sigmaI(j,n),~] = i_2014_active(T_vec(n),rup,site);
    end
end
%% Figure 10
limits = [0.01 10 0 1];
titles = ["Mag = 5, Rrup = 30 km";"Mag = 7, Rrup = 30 km"];
figure('Name','Gregor Figure 10','NumberTitle','off','Position',[10 10 600 400])
for n = 1:2
    subplot(1,2,n)
    semilogx(T_vec, sigmaASK(n,:),'-r',T_vec, sigmaBSSA(n,:),'-g',T_vec, sigmaCB(n,:),'-b',T_vec, sigmaCY(n,:),'-m',T_vec, sigmaI(n,:),'-c','LineWidth',1)
    grid on 
    axis(limits)
    xlabel('Period [sec]')
    ylabel('Sigma [Ln Units]')
    title(titles(n))
    if n == 2
        legend('ASK','BSSA','CB','CY','I','Location','Southeast')
    end
end
%% Save Figure
saveas(gcf,'../figures/gregor10.jpg')