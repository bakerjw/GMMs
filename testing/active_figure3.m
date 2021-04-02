% Script to create a plot of all active models
% Created by Emily Mongold, 1/25/21
%
clear
clc; addpath('../gmms/')
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
    [medveca,~,perioda] = a_2015_active(T,rup);
    [medvecas97,~,periodas97] = as_1997_active(T,rup,site);
    [medvecas8,~,periodas8] = as_2008_active(T,rup,site);
    [medvecask14,~,periodask] = ask_2014_active(T,rup,site);
    [medvecba8,~,periodba] = ba_2008_active(T,rup,site);
    [medvecbjf97,~,periodbjf] = bjf_1997_active(T,rup,site);
    [medvecbssa14,~,periodbssa] = bssa_2014_active(T,rup,site);
    [medvecc97,~,periodc] = c_1997_active(T,rup,site,D);
    [medveccb8,~,periodcb8] = cb_2008_active(T,rup,site);
    [medveccb14,~,periodcb14] = cb_2014_active(T,rup,site);
    [medveccy8,~,periodcy8] = cy_2008_active(T,rup,site);
    [medveccy14,~,periodcy14] = cy_2014_active(T,rup,site);
    [medveci8,~,periodi8] = i_2008_active(T,rup,site);
    [medveci14,~,periodi14] = i_2014_active(T,rup,site);
    [medvecscemy97,~,periodscemy] = scemy_1997_active(T,rup,site);
    [medvecz,~,periodz] = z_2006_active(T,rup,site,sub_ind,MS);

%% Create Figure
limits = [0.01 10 0.01 1];
ascolor = [0,0.4470,0.7410]; bacolor = [0.85,0.3250,0.0980];
cbcolor = [0.929,0.694,0.1250]; cycolor = [0.494,0.184,0.556];
icolor = [0.466,0.674,0.188]; excolor = [0.301,0.745,0.933]; othercolor = [0.635,0.078,0.184];
figure('Name','Active Figure 3')
loglog(perioda, medveca,':','color',bacolor,'Linewidth',1)
hold on
loglog(periodas97, medvecas97,'-.','color',ascolor,'Linewidth',1)
loglog(periodas8, medvecas8,'--','color',ascolor,'Linewidth',1)
loglog(periodask, medvecask14,'-','color',ascolor,'Linewidth',1)
loglog(periodba, medvecba8,'--','color',bacolor,'Linewidth',1)
loglog(periodbjf, medvecbjf97,'-.','color',bacolor,'Linewidth',1)
loglog(periodbssa, medvecbssa14,'-','color',bacolor,'Linewidth',1)
loglog(periodc, medvecc97,'-.','color',excolor,'Linewidth',1)
loglog(periodcb8, medveccb8,'--','color',cbcolor,'Linewidth',1)
loglog(periodcb14, medveccb14,'color',cbcolor,'Linewidth',1)
loglog(periodcy8, medveccy8,'--','color',cycolor,'Linewidth',1)
loglog(periodcy14, medveccy14,'color',cycolor,'Linewidth',1)
loglog(periodi8,medveci8,'--','color',icolor,'Linewidth',1)
loglog(periodi14, medveci14,'color',icolor,'Linewidth',1)
loglog(periodscemy, medvecscemy97,'-.','color',othercolor,'Linewidth',1)
loglog(periodz, medvecz,':','color',excolor,'LineWidth',1)
grid on
axis(limits)
xlabel('Period [sec]')
ylabel('PSA [g]')
title("PSA versus T")
legend('a2015','as1997','as2008','ask2014','ba2008','bjf1997','bssa2014','c1997','cb2008','cb2014','cy2008','cy2014','i2008','i2014','scemy1997','z2006','Location','Southwest')

%% Save Figure
saveas(gcf,'../figures/Active Figure 3.jpg')
