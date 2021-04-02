% Script to recreate Figure 1 from Gregor et. al. (2014)
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
lambda = 0;         % Strike-slip

% Site inputs
is_soil = [];       % Not used in these models
Vs30 = 760;
fvs30 = 0;          % 0 for inferred and 1 for measured
Z25 = 0.6068;       % Used in CB model 
Z10 = 0;            % Used in BSSA model
Z10_cy = 0.0413;     
Z10_ask = 0.0481;
Zbot = 10;      
region = 0;         % = 0 for global

% Create rupture and site objects
rup_cy =    rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
rup =       rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
site_cy =   site(is_soil,Vs30,fvs30,Z25,Z10_cy,Zbot,region);
site_ask =  site(is_soil,Vs30,fvs30,Z25,Z10_ask,Zbot,region);
site =      site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

T =0;  % PGA calculation
rjb_vec = 1:1:250;
%% Function Calls
medianASK = zeros(4,length(rjb_vec));
medianBSSA = zeros(4,length(rjb_vec));
medianCB = zeros(4,length(rjb_vec));
medianCY = zeros(4,length(rjb_vec));
medianI = zeros(4,length(rjb_vec));
for j = 1:4
    rup.M = M(j); rup_cy.M = M(j);
    rup.Ztor = Ztor(j); rup_cy.Ztor = Ztor(j);
    for n = 1:length(rjb_vec)
        rup.Rjb = rjb_vec(n);
        rup_cy.Rjb = rjb_vec(n);
        rup_cy.Rx = rup_cy.Rjb;
        rup_cy.Rrup = sqrt(rup_cy.Rjb^2 + rup_cy.Ztor^2);
        rup.Rrup = sqrt(rup.Rjb^2 + rup.Ztor^2);
        [medianASK(j,n),~,~] = ask_2014_active(T,rup,site_ask);
        [medianBSSA(j,n),~,~] = bssa_2014_active(T,rup,site);
        [medianCB(j,n),~,~] = cb_2014_active(T,rup,site);
        [medianCY(j,n),~,~] = cy_2014_active(T,rup_cy,site_cy);
        [medianI(j,n),~,~] = i_2014_active(T,rup,site);
    end
end
%% Figure 1
limits = [1 250 0.001 1];
titles = ["Mag = 5";"Mag = 6";"Mag = 7";"Mag = 8"];
figure('Name','Gregor Figure 1','Position',[0 0 600 900],'NumberTitle','off')
for n = 1:4
    subplot(2,2,n)
    loglog(rjb_vec, medianASK(n,:),'-r',rjb_vec, medianBSSA(n,:),'-g',rjb_vec, medianCB(n,:),'-b',rjb_vec, medianCY(n,:),'-m',rjb_vec, medianI(n,:),'-c')
    grid on 
    axis(limits)
    xlabel('RJB Distance [km]')
    ylabel('PGA [g]')
    title(titles(n))
    if n == 4
        legend('ASK','BSSA','CB','CY','I','Location','Southwest')
    end
end

%% Save Figure
saveas(gcf,'../figures/gregor1.jpg')
