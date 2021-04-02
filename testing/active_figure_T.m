% Script to create plot of standard deviation versus period
% Based on Figure 10 from Gregor et. al. (2014)
%
% Created by Emily Mongold, 1/24/21
%
clear
clc; addpath('../gmms/')
%% Setting rupture and site object values
% Rupture inputs
M = 7;           
R = [];     
Rrup = 30;      
Rjb = 30;        
Rhyp = 50;      
Rx = 30;        
Ry0 = [];       
HW = 0;             % No hanging wall
AS = 0;             % Not an aftershock
Ztor = 1;       
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
D = 25;
% A1100 (might be optional for CB2008)
MS = 0;         % Does not use magnitude-squared term
sub_ind = [];   % Using for an active event, not subduction

T = 1000;
%% Function Calls
    [~,sigmavec1,period1] = a_2015_active(T,rup);
    [~,sigmavec2,period2] = as_1997_active(T,rup,site); 
    [~,sigmavec3,period3] = as_2008_active(T,rup,site);
    [~,sigmavec4,period4] = ask_2014_active(T,rup,site);
    [~,sigmavec5,period5] = ba_2008_active(T,rup,site);
    [~,sigmavec6,period6] = bjf_1997_active(T,rup,site); 
    [~,sigmavec7,period7] = bssa_2014_active(T,rup,site);
    [~,sigmavec8,period8] = c_1997_active(T,rup,site,D); 
    [~,sigmavec9,period9] = cb_2008_active(T,rup,site);
    [~,sigmavec10,period10] = cb_2014_active(T,rup,site);
    [~,sigmavec11,period11] = cy_2008_active(T,rup,site);
    [~,sigmavec12,period12] = cy_2014_active(T,rup,site);
    [~,sigmavec13,period13] = i_2008_active(T,rup,site);
    [~,sigmavec14,period14] = i_2014_active(T,rup,site);
    [~,sigmavec15,period15] = scemy_1997_active(T,rup,site);
    [~,sigmavec16,period16] = z_2006_active(T,rup,site,sub_ind,MS); 

%% Create Figure
limits = [0.01 10 0 1];
as = [0,0.4470,0.7410]; ba = [0.85,0.3250,0.0980]; 
cb = [0.929,0.694,0.1250]; cy = [0.494,0.184,0.556]; i = [0.466,0.674,0.188]; 
ex = [0.301,0.745,0.933]; other = [0.635,0.078,0.184];
figure('Name',"St.Dev.",'NumberTitle','off')
semilogx(period1, sigmavec1,':','color',ba,'LineWidth',1)
hold on
semilogx(period2, sigmavec2,'-.','color',as,'LineWidth',1)
semilogx(period3, sigmavec3,'--','color',as,'LineWidth',1)
semilogx(period4, sigmavec4,'-','color',as,'LineWidth',1)
semilogx(period5, sigmavec5,'--','color',ba,'LineWidth',1)
semilogx(period6, sigmavec6,'-.','color',ba,'LineWidth',1)
semilogx(period7, sigmavec7,'-','color',ba,'LineWidth',1)
semilogx(period8, sigmavec8,'-.','color',ex,'LineWidth',1)
semilogx(period9, sigmavec9,'--','color',cb,'LineWidth',1)
semilogx(period10, sigmavec10,'-','color',cb,'LineWidth',1)
semilogx(period11, sigmavec11,'--','color',cy,'LineWidth',1)
semilogx(period12, sigmavec12,'-','color',cy,'LineWidth',1)
semilogx(period13, sigmavec13,'--','color',i,'LineWidth',1)
semilogx(period14, sigmavec14,'-','color',i,'LineWidth',1)
semilogx(period15, sigmavec15,'-.','color',other,'LineWidth',1)
semilogx(period16, sigmavec16,':','color',ex,'LineWidth',1)
grid on 
axis(limits)
xlabel('Period [sec]')
ylabel('Sigma [Ln Units]')
title("Sigma versus T")
legend('a2015','as1997','as2008','ask2014','ba2008','bjf1997','bssa2014','c1997','cb2008','cb2014','cy2008','cy2014','i2008','i2014','scemy1997','z2006','Location','SouthEast')
hold off
%% Save Figure
saveas(gcf,'../figures/Active Figure sigma.jpg')