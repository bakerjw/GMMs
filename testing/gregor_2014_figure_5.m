% Script to recreate Figure 5 from Gregor et. al. (2014)
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
M = 7;  
R = [];     
Rrup = 10;      
Rjb = 10;        
Rhyp = [];      
Rx = 10;        
Ry0 = [];       
HW = 0;         % Hanging wall indicator = 1 for hanging wall, 0 otherwise
AS = 0;         % Aftershock indicator = 1 for aftershock, 0 otherwise
Ztor = 1;
Zhyp = 8;
h_eff = [];
W = [];
delta = 0;
lambda = 0;     % Strike-slip

% Site inputs
is_soil = [];
% soil type     = 0 for soil
%               = 1 for soft rock
%               = 2 for hard rock
Vs30 = 190:10:1500;
fvs30 = 0;      % 0 for inferred and 1 for measured
Z25 = [];   % Used in CB model 
Z10 = [];        % Used in BSSA model
Z10_cy = 0.4794;     
Z10_ask = [];
Zbot = 15;      
region = 0;     
% Region        = 0 for global
%               = 1 for California
%               = 2 for Japan
%               = 3 for China
%               = 4 for Italy
%               = 5 for Turkey
%               = 6 for Taiwan
        
rup = rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);    
site = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

period = [0 0.2 1 3];  % Periods for the four plots
%% Function Calls
medianASK = zeros(4,length(Vs30));
medianBSSA = zeros(4,length(Vs30));
medianCB = zeros(4,length(Vs30));
medianCY = zeros(4,length(Vs30));
medianI = zeros(4,length(Vs30));
for j = 1:length(period)
    T = period(j);
    for n = 1:length(Vs30)
        site.Vs30 = Vs30(n);
        [medianASK(j,n),~,~] = ask_2014_active(T,rup,site);
        [medianBSSA(j,n),~,~] = bssa_2014_active(T,rup,site);
        [medianCB(j,n),~,~] = cb_2014_active(T,rup,site);
        [medianCY(j,n),~,~] = cy_2014_active(T,rup,site);
        if Vs30(n) >= 450
           [medianI(j,n),~,~] = i_2014_active(T,rup,site);
        end
    end
end

%% Figure 5
limits = [100 2000 0.01 1];
y_labels = ["PGA [g]";"PSA [g]";"PSA [g]";"PSA [g]"];
titles = ["PGA";"T = 0.2s";"T = 1.0s";"T = 3.0s"];

figure('Name','Gregor Figure 5','Position',[0 0 600 900],'NumberTitle','off')
format shortG
for n = 1:4
    subplot(2,2,n)
    loglog(Vs30,medianASK(n,:),'-r',Vs30, medianBSSA(n,:),'-g',Vs30, medianCB(n,:),'-b',Vs30, medianCY(n,:),'-m',Vs30, medianI(n,:),'-c')
    grid on 
    axis(limits)
    xlabel('Vs30')
    ylabel(y_labels(n))
    title(titles(n))
    if n == 4
        legend('ASK','BSSA','CB','CY','I','Location','Northeast')
    end
end
%% Save Figure
saveas(gcf,'../figures/gregor5.jpg')