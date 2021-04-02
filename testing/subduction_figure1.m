% Script to create a plot of all subduction models
% Based off of Figure 1 from Gregor et. al. (2014)
% 
% Created by Emily Mongold, 2/1/21
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
W = 7;
delta = [];
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
MS = 0;         % does not include magnitude-squared term
Zl = 0;         % for general
sub_ind = [0 1];% for interface and intraslab
F_faba = 1;     % for forearc or unknown sites
T =0;           % PGA calculation
rjb_vec = 1:1:250;
%% Function Calls
medianvec = zeros(5,length(rjb_vec));
for n = 1:length(rjb_vec)
        rup.Rjb = rjb_vec(n); 
        rup.Rrup = sqrt(rup.Rjb^2 + rup.Ztor^2);
        rup.R = rup.Rrup;
        rup.Rhyp = sqrt(rup.Rjb^2 + (rup.Ztor+15)^2);
    for m = 1:length(sub_ind)
        [medianvec(1,n),~,~] = ab_2003_subduction(T,rup,site,sub_ind(m),Zl);
        [medianvec(2,n),~,~] = aga_2016_subduction(T,rup,site,sub_ind(m),F_faba);
        [medianvec(3,n),~,~] = ycsh_1997_subduction(T,rup,site,sub_ind(m));
        [medianvec(4,n),~,~] = gswy_2002_subduction(0.001,rup,site);
        [medianvec(5,n),~,~] = z_2006_active(T,rup,site,sub_ind(m),MS);
    end
end
%% Create Figure
limits = [1 250 0.001 1];
titles = ["PGA versus Rjb, M6.0 Interface" "PGA versus Rjb, M6.0 Intraslab"];
figure('Name','Subduction Figure 1','NumberTitle','off','Position',[10 10 600 400])
% zcolor = [0.4660, 0.6740, 0.1880];
for m = 1:length(sub_ind)
    subplot(1,2,m)
    loglog(rjb_vec, medianvec(1,:),rjb_vec, medianvec(2,:),rjb_vec, medianvec(3,:),rjb_vec, medianvec(4,:),rjb_vec, medianvec(5,:),'LineWidth',1);
%     hold on
%     loglog(rjb_vec, medianvec(5,:),'color',zcolor,'LineWidth',1);
    grid on 
    axis(limits)
    xlabel('RJB Distance [km]')
    ylabel('PGA [g]')
    title(titles(m))
    if m == 2
        legend('ab2003','aga2016','ycsh1997','gswy2002','z2006','Location','Southwest')
    end
end
%% Save Figure
saveas(gcf,'../figures/Subduction Figure 1.jpg')
