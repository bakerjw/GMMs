% Script to create plot of standard deviation versus period
% Based on Figure 10 from Gregor et. al. (2014)
%
% Created by Emily Mongold, 1/25/21
%
clear
clc; addpath('../gmms/')
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
D = 40;
% A1100 (might be optional for CB2008)
MS = 0;         % Does not use magnitude-squared term
sub_ind = [];   % Using for an active event, not subduction

T_vec = 0.01:0.01:10;  % Independent Variable
T_short = 0.05:0.01:4; % For scripts which only go to 4s
T_bjf = 0.01:0.01:2; % BJF model only goes from T = 0.001 to 2s
%% Function Calls
medianvec = zeros(16,length(T_vec));
medianvec2 = zeros(16,length(T_short));
medianvecbjf = zeros(16,length(T_bjf));
for n = 1:length(T_vec)
    [medianvec(3,n),~,~] = as_2008_active(T_vec(n),rup,site);
    [medianvec(4,n),~,~] = ask_2014_active(T_vec(n),rup,site);
    [medianvec(5,n),~,~] = ba_2008_active(T_vec(n),rup,site);
    [medianvec(7,n),~,~] = bssa_2014_active(T_vec(n),rup,site);
    [medianvec(9,n),~,~] = cb_2008_active(T_vec(n),rup,site);
    [medianvec(10,n),~,~] = cb_2014_active(T_vec(n),rup,site);
    [medianvec(11,n),~,~] = cy_2008_active(T_vec(n),rup,site);
    [medianvec(12,n),~,~] = cy_2014_active(T_vec(n),rup,site);
    [medianvec(13,n),~,~] = i_2008_active(T_vec(n),rup,site);
    [medianvec(14,n),~,~] = i_2014_active(T_vec(n),rup,site);
end
for n = 1:length(T_short)
   [medianvec2(1,n),~,~] = a_2015_active(T_short(n),rup); % T only goes
%        0.03 to 5
   [medianvec2(2,n),~,~] = as_1997_active(T_short(n),rup,site); % T only
%        goes from 0.01 to 5
   [medianvec2(8,n),~,~] = c_1997_active(T_short(n),rup,site,D); % T only
%        goes from 0.05 to 4
   [medianvec2(15,n),~,~] = scemy_1997_active(T_short(n),rup,site); % T
%        only goes from 0.001 to 4
   [medianvec2(16,n),~,~] = z_2006_active(T_short(n),rup,site,sub_ind,MS); % T
%        only goes from 0 to 5
end
for n = 1:length(T_bjf)
    [medianvecbjf(6,n),~,~] = bjf_1997_active(T_bjf(n),rup,site); % T only
%        goes from 0.001 to 2
end
%% Create Figure
limits = [0.01 10 0.0001 0.1];
figure('Name','Active Figure 4','NumberTitle','off')
ascolor = [0,0.4470,0.7410]; bacolor = [0.85,0.3250,0.0980]; 
cbcolor = [0.929,0.694,0.1250]; cycolor = [0.494,0.184,0.556]; 
icolor = [0.466,0.674,0.188]; excolor = [0.301,0.745,0.933]; othercolor = [0.635,0.078,0.184];
loglog(T_short, medianvec2(1,:),':','color',bacolor,'LineWidth',1)
hold on
loglog(T_short, medianvec2(2,:),'-.','color',ascolor,'LineWidth',1)
loglog(T_vec, medianvec(3,:),'--','color',ascolor,'LineWidth',1)
loglog(T_vec, medianvec(4,:),'color',ascolor,'LineWidth',1)
loglog(T_vec, medianvec(5,:),'--','color',bacolor,'LineWidth',1)
loglog(T_bjf, medianvecbjf(6,:),'-.','color',bacolor,'LineWidth',1)
loglog(T_vec, medianvec(7,:),'color',bacolor,'LineWidth',1)
loglog(T_short, medianvec2(8,:),'-.','color',excolor,'LineWidth',1)
loglog(T_vec, medianvec(9,:),'--','color',cbcolor,'LineWidth',1)
loglog(T_vec, medianvec(10,:),'color',cbcolor,'LineWidth',1)
loglog(T_vec, medianvec(11,:),'--','color',cycolor,'LineWidth',1)
loglog(T_vec, medianvec(12,:),'color',cycolor,'LineWidth',1)
loglog(T_vec, medianvec(13,:),'--','color',icolor,'LineWidth',1)
loglog(T_vec, medianvec(14,:),'color',icolor,'LineWidth',1)
loglog(T_short, medianvec2(15,:),'-.','color',othercolor,'LineWidth',1)
loglog(T_short, medianvec2(16,:),':','color',excolor,'LineWidth',1)
grid on 
axis(limits)
xlabel('Period [sec]')
ylabel('PSA [g]')
title("PSA versus T")
legend('a2015','as1997','as2008','ask2014','ba2008','bjf1997','bssa2014','c1997','cb2008','cb2014','cy2008','cy2014','i2008','i2014','scemy1997','z2006','Location','Northeast')

%% Save Figure
saveas(gcf,'../figures/Active Figure 4.jpg')