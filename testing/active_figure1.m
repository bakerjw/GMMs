% Script to create a plot of all active models
% Based off of Figure 1 from Gregor et. al. (2014)
% 
% Created by Emily Mongold, 1/22/21
%
clear
clc; addpath('../gmms/')
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
rup =       rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
site =      site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);
D = [];
% A1100 (might be optional for CB2008)
MS = 0;         % Does not use magnitude-squared term
sub_ind = [];   % Using for an active event, not subduction

T =0;     % PGA calculation
rjb_vec = 1:1:250;
%% Function Calls
medianvec = zeros(16,length(rjb_vec));
for n = 1:length(rjb_vec)
        rup.Rjb = rjb_vec(n);
        rup.Rrup = sqrt(rup.Rjb^2 + rup.Ztor^2);
        rup.Rhyp = sqrt(rup.Rjb^2 + (rup.Ztor+15)^2);
        [medianvec(1,n),~,~] = a_2015_active(T,rup);
%         [medianvec(2,n),~,~] = as_1997_active(T,rup,site);
        [medianvec(3,n),~,~] = as_2008_active(T,rup,site);
        [medianvec(4,n),~,~] = ask_2014_active(T,rup,site);
        [medianvec(5,n),~,~] = ba_2008_active(T,rup,site);
%         [medianvec(6,n),~,~] = bjf_1997_active(T,rup,site);
        [medianvec(7,n),~,~] = bssa_2014_active(T,rup,site);
%         [medianvec(8,n),~,~] = c_1997_active(T,rup,site,D);
        [medianvec(9,n),~,~] = cb_2008_active(T,rup,site);
        [medianvec(10,n),~,~] = cb_2014_active(T,rup,site);
        [medianvec(11,n),~,~] = cy_2008_active(T,rup,site);
        [medianvec(12,n),~,~] = cy_2014_active(T,rup,site);
%         [medianvec(13,n),~,~] = i_2008_active(T,rup,site);
        [medianvec(14,n),~,~] = i_2014_active(T,rup,site);
%         [medianvec(15,n),~,~] = scemy_1997_active(T,rup,site);
        [medianvec(16,n),~,~] = z_2006_active(T,rup,site,sub_ind,MS);
end
%% Create Figure
limits = [1 250 0.001 1];
ascolor = [0,0.4470,0.7410]; bacolor = [0.85,0.3250,0.0980]; 
cbcolor = [0.929,0.694,0.1250]; cycolor = [0.494,0.184,0.556]; 
icolor = [0.466,0.674,0.188]; excolor = [0.301,0.745,0.933]; othercolor = [0.635,0.078,0.184];
figure('Name','Active Figure 1')
loglog(rjb_vec, medianvec(1,:),':','color',bacolor,'Linewidth',1)
hold on
loglog(rjb_vec, medianvec(3,:),'--','color',ascolor,'Linewidth',1)
loglog(rjb_vec, medianvec(4,:),'-','color',ascolor,'Linewidth',1)
loglog(rjb_vec, medianvec(5,:),'--','color',bacolor,'Linewidth',1)
loglog(rjb_vec, medianvec(7,:),'-','color',bacolor,'Linewidth',1)
loglog(rjb_vec, medianvec(9,:),'--','color',cbcolor,'Linewidth',1)
loglog(rjb_vec, medianvec(10,:),'-','color',cbcolor,'Linewidth',1)
loglog(rjb_vec, medianvec(11,:),'--','color',cycolor,'Linewidth',1)
loglog(rjb_vec, medianvec(12,:),'-','color',cycolor,'Linewidth',1)
loglog(rjb_vec, medianvec(14,:),'color',icolor,'Linewidth',1)
loglog(rjb_vec, medianvec(16,:),':','color',excolor,'LineWidth',1);
grid on 
axis(limits)
xlabel('RJB Distance [km]')
ylabel('PGA [g]')
title("PGA versus Rjb Distance of a M6.0 Earthquake")
legend('a2015','as2008','ask2014','ba2008','bssa2014','cb2008','cb2014','cy2008','cy2014','i2014','z2006','Location','Southwest')

%% Save Figure
saveas(gcf,'../figures/Active Figure 1.jpg')
