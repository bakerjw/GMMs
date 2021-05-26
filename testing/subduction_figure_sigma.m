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
R = 30;     
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
delta = [];
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
MS = 0; % does not include magnitude-squared term
Zl = 0; % for general
sub_ind = [0 1]; % for interface and intraslab
F_faba = 1; % for forearc or unknown sites
T = 1000;
%% Function Calls
for i = 1:length(sub_ind)
    [~,sigmavec1(i,:),period1(i,:)] = ab_2003_subduction(T,rup,site,sub_ind(i),Zl);
    [~,sigmavec2(i,:),period2(i,:)] = aga_2016_subduction(T,rup,site,sub_ind(i),F_faba);
    [~,sigmavec3(i,:),period3(i,:)] = ycsh_1997_subduction(T,rup,site,sub_ind(i));
    [~,sigmavec4(i,:),period4(i,:)] = gswy_2002_subduction(T,rup,site);
    [~,sigmavec5(i,:),period5(i,:)] = z_2006_active(T,rup,site,sub_ind(i),MS);
end
%% Create Figure
limits = [0.01 10 0 1];
titles = ["Sigma versus T, Interface M7.0" "Sigma versus T, Intraslab M7.0"];
figure('Name',"St.Dev.",'NumberTitle','off','Position',[10 10 600 400])
for n = 1:2
    subplot(1,2,n)
    semilogx(period1(n,:),sigmavec1(n,:),period2(n,:),sigmavec2(n,:),period3(n,:),sigmavec3(n,:),period4(n,:),sigmavec4(n,:),period5(n,:),sigmavec5(n,:),"Linewidth",1)
    grid on 
    axis(limits)
    xlabel('Period [sec]')
    ylabel('Sigma [Ln Units]')
    title(titles(n))
    if n == 2
        legend('ab2003','aga2016','ycsh1997','gswy2002','z2006','Location','Southwest')
    end
end

% set figure size
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 4]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 6 4]);

%% Save Figure
saveas(gcf,'../figures/Subduction Figure sigma.pdf')