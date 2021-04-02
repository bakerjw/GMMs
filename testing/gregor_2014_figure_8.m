% Script to recreate Figure 8 from Gregor et. al. (2014)
% Using the following five GMMs
%   ask_2014_active
%   bssa_2014_active
%   cb_2014_active
%   cy_2014_active
%   i_2014_active
% Created by Emily Mongold, 12/18/20
%
clear
clc; addpath('../gmms/')
%% Setting rupture and site object values
% Rupture inputs
M = [5 6 7 8];           
R = [];     
Rrup = 10;      
Rjb = 10;        
Rhyp = [];      
Rx = 10;        
Ry0 = [];       
HW = 0;             % Hanging wall indicator = 1 for hanging wall, 0 otherwise
AS = 0;             % Aftershock indicator = 1 for aftershock, 0 otherwise
Ztor = [6 3 1 0];   % Different for each magnitude earthquake
Zhyp = 8;
h_eff = [];
W = [];
delta = 15;
lambda = 0;         % Strike-slip

% Site inputs
is_soil = [];       % Not used in these models
Vs30 = 760;
fvs30 = 0;          % 0 for inferred and 1 for measured
Z25 = 1.9826;       % Used in CB model 
Z10 = 0;            % Used in BSSA model
Z10_cy = 0.4794;     
Z10_ask = 0.4704;
Zbot = 10;      
region = 0;         % = 0 for global


% Create rupture and site objects
rup_cy = rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
rup = rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
site_cy = site(is_soil,Vs30,fvs30,Z25,Z10_cy,Zbot,region);
site_ask = site(is_soil,Vs30,fvs30,Z25,Z10_ask,Zbot,region);
site = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

T_vec = 0.01:0.01:10;  % Independent variable
%% Function Calls
medianASK = zeros(4,length(T_vec));
medianBSSA = zeros(4,length(T_vec));
medianCB = zeros(4,length(T_vec));
medianCY = zeros(4,length(T_vec));
medianI = zeros(4,length(T_vec));
for j = 1:4
    rup.M = M(j); rup_cy.M = M(j);
    rup.Ztor = Ztor(j); rup_cy.Ztor = Ztor(j);
    for n = 1:length(T_vec)
        [medianASK(j,n),~,~] = ask_2014_active(T_vec(n),rup,site_ask);
        [medianBSSA(j,n),~,~] = bssa_2014_active(T_vec(n),rup,site);
        [medianCB(j,n),~,~] = cb_2014_active(T_vec(n),rup,site);
        [medianCY(j,n),~,~] = cy_2014_active(T_vec(n),rup_cy,site_cy);
        [medianI(j,n),~,~] = i_2014_active(T_vec(n),rup,site);
    end
end
%% Figure 8
limits = [0.01 10 0.001 1];
titles = ["Mag = 5.0";"Mag = 6.0";"Mag = 7.0";"Mag = 8.0"];
figure('Name','Gregor Figure 8','Position',[0 0 600 900],'NumberTitle','off')
for n = 1:4
    subplot(2,2,n)
    loglog(T_vec, medianASK(n,:),'-r',T_vec, medianBSSA(n,:),'-g',T_vec, medianCB(n,:),'-b',T_vec, medianCY(n,:),'-m',T_vec, medianI(n,:),'-c')
    grid on 
    axis(limits)
    xlabel('Period [sec]')
    ylabel('PSA [g]')
    title(titles(n))
    if n == 3
        legend('ASK','BSSA','CB','CY','I','Location','Southwest')
    end
end
%% Save Figure
saveas(gcf,'../figures/gregor8.jpg')