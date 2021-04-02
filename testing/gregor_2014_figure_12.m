% Script to recreate Figure 12 from Gregor et. al. (2014)
% Using the following five GMMs
%   ask_2014_active
%   bssa_2014_active
%   cb_2014_active
%   cy_2014_active
% Created by Emily Mongold, 12/18/20
%
clear
clc; addpath('../gmms/')
%% Setting rupture and site object values
% Rupture inputs
M = 7;           
R = [];     
Rrup = [];      
Rjb = [];        
Rhyp = [];      
Rx = 30;        
Ry0 = [];       
HW = 0;             % No Hanging Wall
AS = 0;             % Not an aftershock
Ztor = 1;       
Zhyp = 8;
h_eff = [];
W = [];
delta = 15;
lambda = 0;         % Strike-slip

% Site inputs
is_soil = [];       % Not used in these models
Vs30 = 270;
fvs30 = 0;          % 0 for inferred and 1 for measured
Z25 = 1.9826;       % Used in CB model 
Z10 = 0;            % Used in BSSA model
Z10_cy = 0.4794;     
Z10_ask = 0.4704;
Zbot = 10;      
region = 0;         % = 0 for global

% Create rupture and site objects
rup_cy = rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
rup = rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
site_cy = site(is_soil,Vs30,fvs30,Z25,Z10_cy,Zbot,region);
site_ask = site(is_soil,Vs30,fvs30,Z25,Z10_ask,Zbot,region);
site = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

T_vec = [0 1];      % PGA or 1.0 sec
rjb_vec = 1:1:250;
%% Function Calls
sigmaASK = zeros(2,length(rjb_vec));
sigmaBSSA = zeros(2,length(rjb_vec));
sigmaCB = zeros(2,length(rjb_vec));
sigmaCY = zeros(2,length(rjb_vec));
for j = 1:2
    for n = 1:length(rjb_vec)
        rup.Rjb = rjb_vec(n);
        rup_cy.Rjb = rjb_vec(n);
        rup_cy.Rx = rup_cy.Rjb;
        rup_cy.Rhyp = rup_cy.Rjb;
        rup_cy.Rrup = rup_cy.Rjb;
        rup.Rhyp = rup.Rjb;
        rup.Rrup = rup.Rjb;
        [~,sigmaASK(j,n),~] = ask_2014_active(T_vec(j),rup,site_ask);
        [~,sigmaBSSA(j,n),~] = bssa_2014_active(T_vec(j),rup,site);
        [~,sigmaCB(j,n),~] = cb_2014_active(T_vec(j),rup,site);
        [~,sigmaCY(j,n),~] = cy_2014_active(T_vec(j),rup_cy,site_cy);
    end
end
%% Figure 12
limits = [1 250 0 1];
titles = ["PGA";"T = 1.0 sec"];
figure('Name','Gregor Figure 12','NumberTitle','off','Position',[10 10 600 400])
for n = 1:2
    subplot(1,2,n)
    semilogx(rjb_vec,sigmaASK(n,:),'-r',rjb_vec, sigmaBSSA(n,:),'-g',rjb_vec, sigmaCB(n,:),'-b',rjb_vec, sigmaCY(n,:),'-m','LineWidth',1)
    grid on 
    axis(limits)
    xlabel('RJB Distance [km]')
    ylabel('Sigma [Ln Units]')
    title(titles(n))
    if n == 2
        legend('ASK','BSSA','CB','CY','Location','Southeast')
    end
end
%% Save Figure
saveas(gcf,'../figures/gregor12.jpg')