% Script to recreate Figure 9 from Gregor et. al. (2014)
% Using the following five GMMs
%   ask_2014_active
%   bssa_2014_active
%   cb_2014_active
%   cy_2014_active
%   i_2014_active
% Created by Emily Mongold, 12/18/20
%
% NOTE: This file does not recreate the figure as it appears in the text
% for BSSA2014, Vs30 = 270, default Z10,Z25. The outputs of that script do 
% not return the same values as the spreadsheet, but the error was not
% identified as of 21 Jan 2021.
%
clear
clc; addpath('../gmms/')
%% Setting rupture and site object values
% Rupture inputs
M = 7;           
R = [];          
Rjb = 10;        
Rhyp = [];      
Rx = 10;        
Ry0 = [];       
HW = 0;             % Hanging wall indicator = 1 for hanging wall, 0 otherwise
AS = 0;             % Aftershock indicator = 1 for aftershock, 0 otherwise
Ztor = 1;    
Rrup = sqrt(Rjb^2+Ztor^2); 
Zhyp = 8;
h_eff = [];
W = [];
delta = 90;
lambda = 0;     % Strike-slip

% Site inputs
is_soil = [];
Vs30 = 270;         % for first three plots
fvs30 = 0;          % 0 for inferred
% Default values for Vs30 = 760
Z25_760 = 0.6068;       % Used in CB model 
Z10_cy760 = 0.0413;     
Z10_ask760 = 0.0481;
% Default values for Vs30 = 270
Z25_270 = 1.9826;       % Used in CB model 
Z10_cy270 = 0.4794;     
Z10_ask270 = 0.4704;

Z10 = 0;                % Used in BSSA model
Zbot = 10;      
region = 0;             % = 0 for global

% Create rupture and site objects
rup_cy = rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
rup = rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
site_cy = site(is_soil,Vs30,fvs30,[],[],Zbot,region);
site_ask = site(is_soil,Vs30,fvs30,[],[],Zbot,region);
site = site(is_soil,Vs30,fvs30,[],Z10,Zbot,region);

T_vec = 0.01:0.01:10;  % Independent Variable
%% Function Calls
medianASK = zeros(4,length(T_vec));
medianBSSA = zeros(4,length(T_vec));
medianCB = zeros(4,length(T_vec));
medianCY = zeros(4,length(T_vec));
medianI = zeros(4,length(T_vec));
for j = 1:4
    if j == 4
        site.Vs30 = 760; site_ask.Vs30 = 760; site_cy.Vs30 = 760;
        site.Z10 = Z10; site_ask.Z10 = Z10_ask760; site_cy.Z10 = Z10_cy760;
        site.Z25 = Z25_760; site_ask.Z25 = Z25_760; site_cy.Z25 = Z25_760;
    end
    if j == 1
        site.Z10 = Z10; site_ask.Z10 = Z10_ask270; site_cy.Z10 = Z10_cy270;
        site.Z25 = Z25_270; site_ask.Z25 = Z25_270; site_cy.Z25 = Z25_270;
    elseif j == 2
        site.Z10 = 0.1; site_ask.Z10 = 0.1; site_cy.Z10 = 0.1;
        site.Z25 = 0.9; site_ask.Z25 = 0.9; site_cy.Z25 = 0.9;
    elseif j == 3
        site.Z10 = 1.2; site_ask.Z10 = 1.2; site_cy.Z10 = 1.2;
        site.Z25 = 4.8; site_ask.Z25 = 4.8; site_cy.Z25 = 4.8;
    end
    for n = 1:length(T_vec)
        [medianASK(j,n),~,~] = ask_2014_active(T_vec(n),rup,site_ask);
        [medianBSSA(j,n),~,~] = bssa_2014_active(T_vec(n),rup,site);
        [medianCB(j,n),~,~] = cb_2014_active(T_vec(n),rup,site);
        [medianCY(j,n),~,~] = cy_2014_active(T_vec(n),rup_cy,site_cy);
        if j == 4
            [medianI(j,n),~,~] = i_2014_active(T_vec(n),rup,site);
        end
    end
end
%% Figure 9
limits = [0.01 10 0.01 1];
titles = ["Vs30 = 270m/sec, Default Z1, Z25";"Vs30 = 270m/sec, Z1=0.1km, Z25=0.9km";"Vs30 = 270m/sec, Z1=1.2km, Z25=4.8km";"Vs30 = 760m/sec, Default Z1, Z25"];
figure('Name','Gregor Figure 9','Position',[0 0 600 900],'NumberTitle','off')
for n = 1:4
    subplot(2,2,n)
    loglog(T_vec, medianASK(n,:),'-r',T_vec, medianBSSA(n,:),'-g',T_vec, medianCB(n,:),'-b',T_vec, medianCY(n,:),'-m',T_vec, medianI(n,:),'-c')
    grid on 
    axis(limits)
    xlabel('Period [sec]')
    ylabel('PSA [g]')
    title(titles(n))
    if n == 4
        legend('ASK','BSSA','CB','CY','I','Location','Southwest')
    end
end
%% Save Figure
saveas(gcf,'../figures/gregor9.jpg')