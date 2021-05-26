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
clear rup
clear site
clc; 
addpath('../gmms/')
addpath('../')

% call a helper script to build arrays of GMM names and plotting parameters
specify_gmms
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
% Z10 = 0;            % Used in BSSA model
% Z10_cy = [];        % To use default values
% Z10_ask = 0.0481;   % For Vs30 = 760
Z10 = [0.0481 0 0 0.0413 0];    % Values corresponding with the 5 models, 0 if not used (CB and I)
Zbot = 15;      
region = [0 3 2 6];     
% Region            = 0 for global
%                   = 2 for Japan
%                   = 3 for China
%                   = 6 for Taiwan

rupt = rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
sitevar = site(is_soil,Vs30,fvs30,Z25,[],Zbot,region);

T = 0;  % PGA calculation
rjb_vec = 1:1:250;
%% Function Calls
medianvec = cell(1,length(nga2west));
% medianASK2 = zeros(4,length(rjb_vec));
% medianBSSA2 = zeros(4,length(rjb_vec));
% medianCB2 = zeros(4,length(rjb_vec));
% medianCY2 = zeros(4,length(rjb_vec));
for j = 1:4
    sitevar.region = region(j);
    for n = 1:length(rjb_vec)
        rupt.Rjb = rjb_vec(n);
        rupt.Rx = rupt.Rjb;
        rupt.Rrup = sqrt(rupt.Rjb^2+rupt.Ztor^2);
        for i = 1:length(nga2west)
            [medianvec{i}(j,n),~,~] = active_gmms(T,rupt,sitevar,gmm_name{nga2west(i)});
        end
%         [medianASK2(j,n),~,~] = active_gmms(T,rupt,sitevar,'ask_2014');
%         [medianBSSA2(j,n),~,~] = active_gmms(T,rupt,sitevar,'bssa_2014');
%         [medianCB2(j,n),~,~] = active_gmms(T,rupt,sitevar,'cb_2014');
%         [medianCY2(j,n),~,~] = active_gmms(T,rupt,sitevar,'cy_2014');
    end
end
%% Figure 7
limits = [1 250 0.001 1];
gregor_color = ['r','g','b','m','c'];
ask_ind = [1 4 2 3]; bssa_ind = [1 2 3]; cy_ind = [1 3 2];
figure('Name','Gregor Figure 7','Position',[0 0 600 900],'NumberTitle','off')
subplot(2,2,1)
for i = 1:4 % ASK
    loglog(rjb_vec,medianvec{1}(ask_ind(i),:),'color',gregor_color(i));
    hold on
end
grid on
axis(limits)
xlabel('RJB Distance [km]')
ylabel('PGA [g]')
title("ASK")
legend('ASK','Taiwan','China','Japan','Location','Southwest')
hold off

subplot(2,2,2)
for i = 1:3 % BSSA
    loglog(rjb_vec,medianvec{2}(bssa_ind(i),:),'color',gregor_color(i));
    hold on
end
grid on 
axis(limits)
xlabel('RJB Distance [km]')
ylabel('PGA [g]')
title("BSSA")
legend('BSSA','China/Turkey','Italy/Japan','Location','Southwest')
hold off

subplot(2,2,3)
for i = 1:3 % CB
    loglog(rjb_vec,medianvec{3}(bssa_ind(i),:),'color',gregor_color(i));
    hold on
end
grid on 
axis(limits)
xlabel('RJB Distance [km]')
ylabel('PGA [g]')
title("CB")
legend('CB','Japan/Italy','China','Location','Southwest')
hold off

subplot(2,2,4)
for i = 1:3 % CY
    loglog(rjb_vec,medianvec{4}(cy_ind(i),:),'color',gregor_color(i));
    hold on
end
grid on 
axis(limits)
xlabel('RJB Distance [km]')
ylabel('PGA [g]')
title("CY")
legend('CY','Japan/Italy','Wenchuan','Location','Southwest')

% set figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [5 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 5 6]);

%% Save Figure
saveas(gcf,'../figures/gregor7.pdf')
