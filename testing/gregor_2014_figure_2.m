% Script to recreate Figure 2 from Gregor et. al. (2014)
% Using the following GMMs
%   ask_2014_active
%   bssa_2014_active
%   cb_2014_active
%   cy_2014_active
% Created by Emily Mongold, 12/18/20
%
clear
clc; addpath('../gmms/')
%% Setting rupture and site object values
% Rupture inputs
M = [5 6 7 8];     
R = [];     
Rrup = [];      
Rjb = [];        
Rhyp = [];      
Rx = 10;        
Ry0 = [];       
HW = 0;             % Hanging wall indicator = 1 for hanging wall, 0 otherwise
AS = 0;             % Aftershock indicator = 1 for aftershock, 0 otherwise
Ztor = [6 3 1 0];   % Different for each magnitude earthquake
Zhyp = 8;
h_eff = [];
W = [];
delta = 90;
lambda = 0;     % Strike-slip

% Site inputs
is_soil = [];
% soil type         = 0 for soil
%                   = 1 for soft rock
%                   = 2 for hard rock
Vs30 = 270;
fvs30 = 0;          % 0 for inferred and 1 for measured
Z25 = 1.9826;       % Used in CB model 
Z10 = 0;            % Used in BSSA model
Z10_cy = 0.4794;     
Z10_ask = 0.4704;
Zbot = 15;      
region = 0;     
% Region            = 0 for global
%                   = 1 for California
%                   = 2 for Japan
%                   = 3 for China
%                   = 4 for Italy
%                   = 5 for Turkey
%                   = 6 for Taiwan

rup = rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
site_cy = site(is_soil,Vs30,fvs30,Z25,Z10_cy,Zbot,region);
site_ask = site(is_soil,Vs30,fvs30,Z25,Z10_ask,Zbot,region);
site = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

T = 1.0;  % PSA calculation
rjb_vec = 1:1:250;
%% Function Calls
medianASK = zeros(4,length(rjb_vec));
medianBSSA = zeros(4,length(rjb_vec));
medianCB = zeros(4,length(rjb_vec));
medianCY = zeros(4,length(rjb_vec));
for j = 1:4
    rup.M = M(j);
    rup.Ztor = Ztor(j);
    for n = 1:length(rjb_vec)
        rup.Rjb = rjb_vec(n);
        rup.Rx = rup.Rjb;
        rup.Rrup = sqrt(rup.Rjb^2 + rup.Ztor^2);
        [medianASK(j,n),~,~] = ask_2014_active(T,rup,site_ask);
        [medianBSSA(j,n),~,~] = bssa_2014_active(T,rup,site);
        [medianCB(j,n),~,~] = cb_2014_active(T,rup,site);
        [medianCY(j,n),~,~] = cy_2014_active(T,rup,site_cy);
    end
end
%% Figure 2
limits = [1 250 0.001 3];
titles = ["Mag = 5";"Mag = 6";"Mag = 7";"Mag = 8"];
figure('Name','Gregor Figure 2','Position',[0 0 600 900],'NumberTitle','off')
for n = 1:4
    subplot(2,2,n)
    loglog(rjb_vec, medianASK(n,:),'-r',rjb_vec, medianBSSA(n,:),'-g',rjb_vec, medianCB(n,:),'-b',rjb_vec, medianCY(n,:),'-m')
    grid on 
    axis(limits)
    xlabel('RJB Distance [km]')
    ylabel('PSA [g]')
    title(titles(n))
    if n == 4
        legend('ASK','BSSA','CB','CY','Location','Southwest')
    end
end
%% Save Figure
saveas(gcf,'../figures/gregor2.jpg')