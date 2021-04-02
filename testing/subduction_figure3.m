% Script to create a plot of all subduction models
% Based off of Figure 1 from Gregor et. al. (2014)
% 
% Created by Emily Mongold, 2/1/21
%
clear
clc; addpath('../gmms/')
%% Setting rupture and site object values
% Rupture inputs
M = 4:0.5:8;           
R = 80;     
Rrup = 80;      
Rjb = 100;        
Rhyp = 90;      
Rx = 100;        
Ry0 = [];       
HW = 0;             % No hanging wall
AS = 0;             % Not an aftershock
Ztor = 0; 
Zhyp = 8;
h_eff = [];
W = 7;
delta = [];
lambda = 0;         % Strike-slip

% Site inputs
is_soil = 0;        % Soil
Vs30 = 270;
fvs30 = 0;          % Inferred
Z25 = 0.8;      
Z10 = 0.05;     
Zbot = 15;      
region = 0;         % Global, no corrections applied

%% Create rupture and site objects, other inputs
rup =       rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
site =      site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);
MS = 0;         % does not include magnitude-squared term
Zl = 0;         % for general
sub_ind = [0 1];% for interface and intraslab
F_faba = 0;     % for forearc or unknown sites
T =1;           % PSA calculation
%% Function Calls
medianvec = zeros(5,length(M));
for n = 1:length(M)
        rup.M = M(n); 
    for m = 1:length(sub_ind)
        [medianvec(1,n,m),~,~] = ab_2003_subduction(T,rup,site,sub_ind(m),Zl);
        [medianvec(2,n,m),~,~] = aga_2016_subduction(T,rup,site,sub_ind(m),F_faba);
        [medianvec(3,n,m),~,~] = ycsh_1997_subduction(T,rup,site,sub_ind(m));
        [medianvec(4,n,m),~,~] = gswy_2002_subduction(T,rup,site);
        [medianvec(5,n,m),~,~] = z_2006_active(T,rup,site,sub_ind(m),MS);
    end
end
%% Create Figure
limits = [4 8 0.001 1];
titles = ["PGA versus M, Intraslab" "PGA versus M, Intraslab"];
figure('Name','Subduction Figure 3','NumberTitle','off','Position',[10 10 600 400])
for m = 1:2
    subplot(1,2,m)
    semilogy(M, medianvec(1,:,m),M, medianvec(2,:,m),M, medianvec(3,:,m),M, medianvec(4,:,m),M, medianvec(5,:,m),'LineWidth',1);
    grid on 
    axis(limits)
    xlabel('Magnitude')
    ylabel('PGA [g]')
    title(titles(m))
    if m ==2
        legend('ab2003','aga2016','ycsh1997','gswy2002','z2006','Location','Northwest')
    end
end

%% Save Figure
saveas(gcf,'../figures/Subduction Figure 3.jpg')