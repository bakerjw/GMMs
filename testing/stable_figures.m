% Script to create figures for the stable scripts
% Created by Emily Mongold 3/2/2021
%
clear
clc
addpath('../gmms/')
%% Setting rupture and site object values
% Rupture inputs
M = [5 6 7 8];           
R = [];        
Rrup = 1:1:250;     % Used in ab_2006_stable
Rjb = 1:1:250;      % Used in sp_2016_stable
Rhyp = [];      
Rx = [];        
Ry0 = [];       
HW = [];             
AS = [];
Ztor = [];       
Zhyp = [];
h_eff = [];
W = [];
delta = [];
lambda = [];         

% Site inputs
is_soil = [];        
Vs30 = 760;         % Used by ab_2006_stable
fvs30 = [];          
Z25 = [];       
Z10 = [];          
Zbot = [];      
region = [];         

% Create rupture and site objects
rup = rup([],R,[],[],Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
site = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

stress = 140;                   % Used by ab_2006_stable, between 35 bars and 560 bars
T = 0.01:0.01:5; 
%% Creation of Stable Figures
median = zeros(length(Rrup),length(M),2);
for n = 1:length(M)
    rup.M = M(n);
    for j = 1:length(Rrup)
        rup.Rrup = Rrup(j);
        rup.Rjb = Rjb(j);
        [median(j,n,1), sigma(j,n,1),period(1,:)] = ab_2006_stable(1,rup,site,stress); % T = 1.0s
        [median(j,n,2), sigma(j,n,2),period(2,:)] = sp_2016_stable(1,rup); % T = 1.0s
    end
end
% Create Stable Figure 1
figure('Name','Stable Figure 1','NumberTitle','off','Position',[10 10 600 400])
titles = ["M5.0","M6.0","M7.0","M8.0"];
limits = [1 250 0.001 1];
for n = 1:length(M)
    subplot(2,2,n)
    loglog(Rrup,median(:,n,1),Rrup,median(:,n,2),'LineWidth',1); 
    title(titles(n));
    grid on
    axis(limits)
    xlabel('Distance [km]')
    ylabel('PSA [g]')
    if n == 4
        legend("ab2006","sp2016","Location","NorthEast");
    end
end
hold off
% Save Figure
saveas(gcf,'../figures/Stable Figure 1.jpg')

% Inputs and function call for Figure 2
median = zeros(length(M),length(T),2);
rup.Rrup = 100;
rup.Rjb = 100;
for n = 1:length(M)
    rup.M = M(n);
    for i = 1:length(T)
        [median(n,i,1), sigma(n,i,1),period(1,:)] = ab_2006_stable(T(i),rup,site,stress);
        [median(n,i,2), sigma(n,i,2),period(2,:)] = sp_2016_stable(T(i),rup);
    end
end
% Create Stable Figure 2
figure('Name','Stable Figure 2','NumberTitle','off','Position',[10 10 600 400])
titles = ["M5.0","M6.0","M7.0","M8.0"];
limits = [0.01 5 0.0001 0.5];
for n = 1:length(M)
    subplot(2,2,n)
    loglog(T,median(n,:,1),T,median(n,:,2),'LineWidth',1);
    title(titles(n));
    grid on
    axis(limits)
    xlabel('Period [s]')
    ylabel('PSA [g]')
    if n == 4
        legend("ab2006","sp2016","Location","NorthEast");
    end
end
hold off
% Save Figure
saveas(gcf,'../figures/Stable Figure 2.jpg')

% Create Stable Figure 3
figure('Name','Stable Figure 3','NumberTitle','off','Position',[10 10 600 400])
limits = [0.5 4 0 1];
semilogx(T,sigma(3,:,1),T,sigma(3,:,2),'LineWidth',1);
title("Standard Deviation vs Period for M7.0");
grid on 
axis(limits)
xlabel('Period [s]')
ylabel('Sigma [ln units]')
legend("as1997","c1997","Location","NorthEast");
% Save Figure
saveas(gcf,'../figures/Stable Figure sigma.jpg')
