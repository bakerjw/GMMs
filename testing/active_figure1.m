% Script to create a plot of all active models
% Based off of Figure 1 from Gregor et. al. (2014)
% 
% Created by Emily Mongold, 1/22/21
%
% clear
clc
addpath('../gmms/')
addpath('../')

% call a helper script to build arrays of GMM names and plotting parameters
specify_gmms

%% Setting rupture and site object values
% Rupture inputs
M = 6;           
R = [];     
Rrup = [];      
Rjb = [];        
Rhyp = [];      
Rx = 10;        
Ry0 = [];       
HW = 0;             % No hanging wall
AS = 0;             % Not an aftershock
Ztor = 3; 
Zhyp = 8;
h_eff = [];
W = 10;
delta = 90;
lambda = 0;         % Strike-slip

% Site inputs
is_soil = 1;        % Soft rock
Vs30 = 760;
fvs30 = 0;          % Inferred
Z25 = 0.8;      
Z10 = 0.05;     
Zbot = 15;      
region = 0;         % Global, no corrections applied

%% Create rupture and site objects, other inputs
rupt =       rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
sitevar =      site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);
D = [];
% A1100 (might be optional for CB2008)
% MS = 0;         % Does not use magnitude-squared term
% sub_ind = [];   % Using for an active event, not subduction

T = 0;     % PGA calculation
rjb_vec = 1:1:250;
%% Function Calls
pga_gmms = [1 3 4 5 7 8 9 10 11 13 15];
medianvec = zeros(16,length(rjb_vec));
for n = 1:length(rjb_vec)
        rupt.Rjb = rjb_vec(n);
        rupt.Rrup = sqrt(rupt.Rjb^2 + rupt.Ztor^2);
        rupt.Rhyp = sqrt(rupt.Rjb^2 + (rupt.Ztor+15)^2);
        for i = pga_gmms
            [medianvec(i,n),~,~] = active_gmms(T,rupt,sitevar,gmm_name{i});
        end
end
%% Create Figure
limits = [1 250 0.001 1];
figure('Name','Active Figure 1')
for i = pga_gmms
    loglog(rjb_vec, medianvec(i,:),line_style{i},'color',line_color{i},'Linewidth',1)
    hold on
end
grid on 
axis(limits)
xlabel('R_{JB} Distance [km]')
ylabel('PGA [g]')
title("PGA versus R_{JB} Distance of a M6.0 Earthquake")
legend(gmm_name{pga_gmms},'Location','Southwest','Interpreter','none')

%% Save Figure
saveas(gcf,'../figures/Active Figure 1.jpg')
