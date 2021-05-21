% Script to create a plot of all active models
% Created by Emily Mongold, 1/25/21
%
clear
clc; addpath('../gmms/')
%% Setting rupture and site object values
% Rupture inputs
M = 8;           
R = 25;     
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
MS = 0;         % does not include magnitude-squared term
Zl = 0;         % for general
sub_ind = [0 1];% for interface and intraslab
F_faba = 1;     % for forearc or unknown sites

T = 1000;
%% Function Calls
for m = 1:length(sub_ind)
    [medvecab(:,m),~,periodab(:,m)] = ab_2003_subduction(T,rup,site,sub_ind(m),Zl);
    [medvecaga(:,m),~,periodaga(:,m)] = aga_2016_subduction(T,rup,site,sub_ind(m),F_faba);
    [medvecycsh(:,m),~,periodycsh(:,m)] = ycsh_1997_subduction(T,rup,site,sub_ind(m));
    [medvecgswy(:,m),~,periodgswy(:,m)] = gswy_2002_subduction(T,rup,site);
    [medvecz(:,m),~,periodz(:,m)] = z_2006_active(T,rup,site,sub_ind(m),MS);
end
%% Create Figure
limits = [0.01 10 0.01 2];
titles = ["SA versus T, Interface" "SA versus T, Intraslab"];
figure('Name','Subduction Figure 4','NumberTitle','off','Position',[10 10 600 400])
for m = 1:length(sub_ind)
    subplot(1,2,m)
    loglog(periodab(:,m), medvecab(:,m),periodaga(:,m), medvecaga(:,m),periodycsh(:,m), medvecycsh(:,m),periodgswy(:,m), medvecgswy(:,m),periodz(:,m), medvecz(:,m),'LineWidth',1)
    grid on 
    axis(limits)
    xlabel('Period [sec]')
    ylabel('SA [g]')
    title(titles(m))
    if m == 2
        legend('ab2003','aga2016','ycsh1997','gswy2002','z2006','Location','Southwest')
    end
end
%% Save Figure
saveas(gcf,'../figures/Subduction Figure 4.jpg')