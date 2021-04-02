% Script to recreate Figure 3 from Gregor et. al. (2014)
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
M = 3:1:9;  
R = [];     
Rrup = 30;      
Rjb = 30;        
Rhyp = [];      
Rx = 10;        
Ry0 = [];       
HW = 0;         % Hanging wall indicator = 1 for hanging wall, 0 otherwise
AS = 0;         % Aftershock indicator = 1 for aftershock, 0 otherwise
Ztor = [];
Zhyp = 8;
h_eff = [];
W = [];
delta = 15;
lambda = 0;     % Strike-slip

% Site inputs
is_soil = [];
% soil type     = 0 for soil
%               = 1 for soft rock
%               = 2 for hard rock
Vs30 = 760;
fvs30 = 0;      % 0 for inferred and 1 for measured
Z25 = 0.6068;   % Used in CB model 
Z10 = 0;        % Used in BSSA model
Z10_cy = 0.0413;     
Z10_ask = 0.0481;
Zbot = 15;      
region = 2;     
% Region        = 0 for global
%               = 1 for California
%               = 2 for Japan
%               = 3 for China
%               = 4 for Italy
%               = 5 for Turkey
%               = 6 for Taiwan

rup_i = rup(3,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
rup = rup(3,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
site_cy = site(is_soil,Vs30,fvs30,Z25,Z10_cy,Zbot,region);
site_ask = site(is_soil,Vs30,fvs30,Z25,Z10_ask,Zbot,region);
site = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

period = [0 0.2 1 3];  % Periods for the four plots
%% Function Calls
% PGA
medianASK = zeros(4,length(M));
medianBSSA = zeros(4,length(M));
medianCB = zeros(4,length(M));
medianCY = zeros(4,length(M));
medianI = zeros(4,length(M));
for j = 1:4
    T = period(j);
    for n = 1:length(M)
        if M(n)>=7
            rup.Zhyp = 12;rup_i.Zhyp = 12;
        end
        rup.M = M(n); rup_i.M = M(n);

        rup.Ztor = (max(2.673-1.136*max(rup.M-4.970,0),0))^2;
        rup_i.Ztor = rup.Ztor;

        [medianASK(j,n),~,~] = ask_2014_active(T,rup,site_ask);
        [medianBSSA(j,n),~,~] = bssa_2014_active(T,rup,site);
        [medianCB(j,n),~,~] = cb_2014_active(T,rup,site);
        [medianCY(j,n),~,~] = cy_2014_active(T,rup,site_cy);
        if M(n) >= 5
            [medianI(j,n),~,~] = i_2014_active(T,rup_i,site);
        end
    end
end
%% Figure 3
limits = [3 8.5 0.0001 1];
titles = ["PGA" "T = 0.2s" "T = 1.0s" "T = 3.0s"];
figure('Name','Gregor Figure 3','Position',[0 0 600 900],'NumberTitle','off')
for n = 1:4
    subplot(2,2,n)
    semilogy(M, medianASK(n,:),'-r',M, medianBSSA(n,:),'-g',M, medianCB(n,:),'-b',M, medianCY(n,:),'-m',M, medianI(n,:),'-c')
    grid on
    grid minor
    axis(limits)
    xlabel('Magnitude')
    if n == 1
        ylabel('PGA [g]')
    else
        ylabel('PSA [g]')
    end
    title(titles(n))
    if n == 4
        legend('ASK','BSSA','CB','CY','I','Location','Northwest')
    end
end

%% Save Figure
saveas(gcf,'../figures/gregor3.jpg')