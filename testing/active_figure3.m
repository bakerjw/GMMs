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

% [medvecas97,~,periodas97] = as_1997_active(T,rup,site);
%     [medvecas8,~,periodas8] = as_2008_active(T,rup,site);
%     [medvecask14,~,periodask] = ask_2014_active(T,rup,site);
%     [medvecba8,~,periodba] = ba_2008_active(T,rup,site);
%     [medvecbjf97,~,periodbjf] = bjf_1997_active(T,rup,site);
%     [medvecbssa14,~,periodbssa] = bssa_2014_active(T,rup,site);
%     [medvecc97,~,periodc] = c_1997_active(T,rup,site,D);
%     [medveccb8,~,periodcb8] = cb_2008_active(T,rup,site);
%     [medveccb14,~,periodcb14] = cb_2014_active(T,rup,site);
%     [medveccy8,~,periodcy8] = cy_2008_active(T,rup,site);
%     [medveccy14,~,periodcy14] = cy_2014_active(T,rup,site);
%     [medveci8,~,periodi8] = i_2008_active(T,rup,site);
%     [medveci14,~,periodi14] = i_2014_active(T,rup,site);
%     [medvecscemy97,~,periodscemy] = scemy_1997_active(T,rup,site);
%     [medvecz,~,periodz] = z_2006_active(T,rup,site,sub_ind,MS);

%% Create Figure
% limits = [0.01 10 0.01 1];

figure('Name','Active Figure 3')
for i = 2:length(gmm_name)
    loglog(periodvec{i},medvec{i},line_style{i},'color',line_color{i},'Linewidth',1)
    hold on
end
% loglog(periodas8, medvecas8,'--','color',ascolor,'Linewidth',1)
% loglog(periodask, medvecask14,'-','color',ascolor,'Linewidth',1)
% loglog(periodba, medvecba8,'--','color',bacolor,'Linewidth',1)
% loglog(periodbjf, medvecbjf97,'-.','color',bacolor,'Linewidth',1)
% loglog(periodbssa, medvecbssa14,'-','color',bacolor,'Linewidth',1)
% loglog(periodc, medvecc97,'-.','color',excolor,'Linewidth',1)
% loglog(periodcb8, medveccb8,'--','color',cbcolor,'Linewidth',1)
% loglog(periodcb14, medveccb14,'color',cbcolor,'Linewidth',1)
% loglog(periodcy8, medveccy8,'--','color',cycolor,'Linewidth',1)
% loglog(periodcy14, medveccy14,'color',cycolor,'Linewidth',1)
% loglog(periodi8,medveci8,'--','color',icolor,'Linewidth',1)
% loglog(periodi14, medveci14,'color',icolor,'Linewidth',1)
% loglog(periodscemy, medvecscemy97,'-.','color',othercolor,'Linewidth',1)
% loglog(periodz, medvecz,':','color',excolor,'LineWidth',1)
grid on
% axis(limits)
xlabel('Period [sec]')
ylabel('SA [g]')
title("SA versus T")
legend(gmm_name{2:end},'Location','Southwest','Interpreter','none')
%% Save Figure
saveas(gcf,'../figures/Active Figure 3new.jpg')
