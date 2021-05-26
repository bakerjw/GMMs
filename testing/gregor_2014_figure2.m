% Script to recreate Figure 2 from Gregor et. al. (2014)
% Using the following GMMs
%   ask_2014_active
%   bssa_2014_active
%   cb_2014_active
%   cy_2014_active
% Created by Emily Mongold, 12/18/20
%
clear rup
clear site
clc
addpath('../gmms/')
addpath('../')

% call a helper script to build arrays of GMM names and plotting parameters
specify_gmms

%% Setting rupture and site object values
% Rupture inputs
M = [5 6 7 8];     
R = [];  Rrup = [];   Rjb = [];   Rhyp = [];  Ry0 = [];       
Rx = 10;        
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
Vs30 = 270;
fvs30 = 0;          % 0 for inferred and 1 for measured
Z25 = 1.9826;       % Used in CB model 
Z10 = [];           % Changed for each model
Zbot = 15;      
region = 0;         %  = 0 for global

rupt = rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
sitevar = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

T = 1.0;  % PSA calculation
rjb_vec = 1:1:250;
%% Function Calls
medianASK2 = zeros(4,length(rjb_vec));
medianBSSA2 = zeros(4,length(rjb_vec));
medianCB2 = zeros(4,length(rjb_vec));
medianCY2 = zeros(4,length(rjb_vec));
for j = 1:4
    rupt.M = M(j);
    rupt.Ztor = Ztor(j);
    for n = 1:length(rjb_vec)
        rupt.Rjb = rjb_vec(n);
        rupt.Rx = rupt.Rjb;
        rupt.Rrup = sqrt(rupt.Rjb^2 + rupt.Ztor^2);
        sitevar.Z10 = 0.4704;
        [medianASK2(j,n),~,~] = active_gmms(T,rupt,sitevar,'ask_2014');
        sitevar.Z10 = 0; 
        [medianBSSA2(j,n),~,~] = active_gmms(T,rupt,sitevar,'bssa_2014');
        [medianCB2(j,n),~,~] = active_gmms(T,rupt,sitevar,'cb_2014');
        sitevar.Z10 = 0.4794;
        [medianCY2(j,n),~,~] = active_gmms(T,rupt,sitevar,'cy_2014');
    end
end
%% Figure 2
limits = [1 250 0.001 3];
titles = ["Mag = 5";"Mag = 6";"Mag = 7";"Mag = 8"];
figure('Name','Gregor Figure 2','Position',[0 0 600 900],'NumberTitle','off')
for n = 1:4
    subplot(2,2,n)
    loglog(rjb_vec, medianASK2(n,:),'-r',rjb_vec, medianBSSA2(n,:),'-g',rjb_vec, medianCB2(n,:),'-b',rjb_vec, medianCY2(n,:),'-m')
    grid on 
    axis(limits)
    xlabel('RJB Distance [km]')
    ylabel('PSA [g]')
    title(titles(n))
    if n == 4
        legend(gmm_name{nga2west((1:4))},'Location','Southwest','Interpreter','none')
    end
end

% set figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [5 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 5 6]);

%% Save Figure
saveas(gcf,'../figures/gregor2.pdf')