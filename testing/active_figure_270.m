% Script to create a plot of all active models
% Created by Emily Mongold, 1/22/21
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
M = 6;           
R = [];     
Rrup = [];      
Rjb = [];        
Rhyp = [];      
Rx = 10;        
Ry0 = [];       
HW = 0;             % No hanging wall
AS = 0;             % Not an aftershock
Ztor = 3; 
Zhyp = 8;
h_eff = [];
W = 10;
delta = 90;
lambda = 0;     % Strike-slip

% Site inputs
is_soil = 1;    % Soft rock
Vs30 = 270;
fvs30 = 0;      % inferred
Z25 = 0.8;      
Z10 = 0.05;     
Zbot = 15;      
region = 0;     % global, no corrections applied

%% Create rupture and site objects, other inputs
rup =       rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
site =      site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

T =1;     % PSA calculation
rjb_vec = 1:1:250;
%% Function Calls
gmm_vec = [1:11 13:15]; % Excluding I2008 because Vs30 is out of its range
medianvec = cell(1,length(gmm_name));
for n = 1:length(rjb_vec)
        rup.Rjb = rjb_vec(n);
        rup.Rrup = sqrt(rup.Rjb^2 + rup.Ztor^2);
        rup.Rhyp = sqrt(rup.Rjb^2 + (rup.Ztor+15)^2);
        for i = gmm_vec
            medianvec{i}(n) = active_gmms(T,rup,site,gmm_name{i});
        end
end
%% Create Figure
limits = [1 250 0.001 1];
figure('Name','Active Figure 270')
for i = gmm_vec
    loglog(rjb_vec,medianvec{i},line_style{i},'color',line_color{i},'Linewidth',1);
    hold on
end
grid on 
axis(limits)
xlabel('RJB Distance [km]')
ylabel('SA [g]')
title("SA versus Rjb Distance")
legend(gmm_name{gmm_vec},'Location','Southwest','Interpreter','none')


% set figure size
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [5 5]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 5 5]);
    
%% Save Figure
saveas(gcf,'../figures/Active Figure 270.pdf')