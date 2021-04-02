% Script to recreate Figure 7 from Gregor et. al. (2014)
% Using the following GMMs
%   ask_2014_active
%   bssa_2014_active
%   cb_2014_active
%   cy_2014_active
% Created by Emily Mongold, 12/18/20
%
% NOTE: This file does not recreate the figure as it appears in the text
% for CY2014, Japan/Italy. The outputs of that script return the same
% values as the spreadsheet, and no changes in unspecified input variables
% changed the curve to match as it appears in Gregor et. al. (2014). 
%
clear
clc; addpath('../gmms/')
%% Setting rupture and site object values
% Rupture inputs
M = 7;     
R = [];     
Rrup = [];      
Rjb = [];        
Rhyp = [];      
Rx = [];        
Ry0 = [];       
HW = 1;             % Hanging wall indicator = 1 for hanging wall, 0 otherwise
AS = 0;             % Aftershock indicator = 1 for aftershock, 0 otherwise
Ztor = 1;           % Different for each magnitude earthquake
Zhyp = 8;
h_eff = [];
W = [];
delta = 90;
lambda = 0;         % Strike-slip

% Site inputs
is_soil = [];       % Not used for these models
Vs30 = 760;
fvs30 = 0;          % 0 for inferred and 1 for measured
Z25 = 0.6068;       % Used in CB model 
Z10 = 0;            % Used in BSSA model
Z10_cy = [];        % To use default values
Z10_ask = 0.0481;   % For Vs30 = 760
Zbot = 15;      
region = [0 3 2 6];     
% Region            = 0 for global
%                   = 2 for Japan
%                   = 3 for China
%                   = 6 for Taiwan

rup = rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
site_cy = site(is_soil,Vs30,fvs30,Z25,Z10_cy,Zbot,region);
site_ask = site(is_soil,Vs30,fvs30,Z25,Z10_ask,Zbot,region);
site = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

T = 0;  % PGA calculation
rjb_vec = 1:1:250;
%% Function Calls
medianASK = zeros(4,length(rjb_vec));
medianBSSA = zeros(4,length(rjb_vec));
medianCB = zeros(4,length(rjb_vec));
medianCY = zeros(4,length(rjb_vec));
for j = 1:4
    site.region = region(j); site_cy.region = region(j);
    site_ask.region = region(j);
    for n = 1:length(rjb_vec)
        rup.Rjb = rjb_vec(n);
        rup.Rx = rup.Rjb;
        rup.Rrup = sqrt(rup.Rjb^2+rup.Ztor^2);
        [medianASK(j,n),~,~] = ask_2014_active(T,rup,site_ask);
        [medianBSSA(j,n),~,~] = bssa_2014_active(T,rup,site);
        [medianCB(j,n),~,~] = cb_2014_active(T,rup,site);
        [medianCY(j,n),~,~] = cy_2014_active(T,rup,site_cy);
    end
end
%% Figure 7
limits = [1 250 0.001 1];
figure('Name','Gregor Figure 7','Position',[0 0 600 900],'NumberTitle','off')
subplot(2,2,1)
loglog(rjb_vec, medianASK(1,:),'-r',rjb_vec, medianASK(4,:),'-g',rjb_vec, medianASK(2,:),'-b',rjb_vec, medianASK(3,:),'-m')
grid on 
axis(limits)
xlabel('RJB Distance [km]')
ylabel('PGA [g]')
title("ASK")
legend('ASK','Taiwan','China','Japan','Location','Southwest')

subplot(2,2,2)
loglog(rjb_vec, medianBSSA(1,:),'-r',rjb_vec, medianBSSA(2,:),'-g',rjb_vec, medianBSSA(3,:),'-b')
grid on 
axis(limits)
xlabel('RJB Distance [km]')
ylabel('PGA [g]')
title("BSSA")
legend('BSSA','China/Turkey','Italy/Japan','Location','Southwest')

subplot(2,2,3)
loglog(rjb_vec, medianCB(1,:),'-r',rjb_vec, medianCB(3,:),'-g',rjb_vec, medianCB(2,:),'-b')
grid on 
axis(limits)
xlabel('RJB Distance [km]')
ylabel('PGA [g]')
title("CB")
legend('CB','Japan/Italy','China','Location','Southwest')

subplot(2,2,4)
loglog(rjb_vec, medianCY(1,:),'-r',rjb_vec, medianCY(3,:),'-g',rjb_vec, medianCY(2,:),'-b')
grid on 
axis(limits)
xlabel('RJB Distance [km]')
ylabel('PGA [g]')
title("CY")
legend('CY','Japan/Italy','Wenchuan','Location','Southwest')

%% Save Figure
saveas(gcf,'../figures/gregor7.jpg')
