% Script to recreate Figure 6 from Gregor et. al. (2014)
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
M = 6.7;  
R = [];     
Rrup = [];      
Rjb = [];        
Rhyp = [];      
Rx = -30:1:40;        
Ry0 = [];       
HW = 1;         % Hanging wall
AS = 0;         % Not an aftershock
Ztor = [];      % 0 or 6
Zhyp = 8;
h_eff = [];
W = 21.97;      
delta = 45;     % Dip angle of 45 degrees
lambda = [];    % Normal (top) and reverse (bottom)

% Site inputs
is_soil = [];   % Not used for these models
Vs30 = 760;
fvs30 = 0;      % 0 for inferred and 1 for measured
Z25 = 1.9826;   % Used in CB model 
Z10 = 0;        % Used in BSSA model
Z10_cy = 0.4794;     
Z10_ask = 0.4704;
Zbot = 15;      
region = 0;     % Global
        
rup = rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
site_cy = site(is_soil,Vs30,fvs30,Z25,Z10_cy,Zbot,region);
site_ask = site(is_soil,Vs30,fvs30,Z25,Z10_ask,Zbot,region);
site = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

T = 0;      % PGA
%% Function Calls
medianASK = zeros(4,length(Rx));
medianBSSA = zeros(4,length(Rx));
medianCB = zeros(4,length(Rx));
medianCY = zeros(4,length(Rx));
medianI = zeros(4,length(Rx));
for j = 1:4
    % Surface or Buried
    if j==2 || j==4
        rup.Ztor = 6;
    else 
        rup.Ztor = 0;
    end
    % Normal or Reverse Faulting
    if j <= 2
        rup.lambda = -90;
    else
        rup.lambda = 90;
    end
    for n = 1:length(Rx)
        rup.Rx = Rx(n); 
        if Rx(n)<0
            rup.Rrup = sqrt(Rx(n)^2 + rup.Ztor^2); rup.Rjb = abs(Rx(n)); 
        elseif Rx(n)>=0   
            rup.Rrup = (Rx(n) + rup.Ztor)*cos(rup.delta);
            if Rx(n)<=16
                rup.Rjb = 0; 
            else
                rup.Rjb = Rx(n)-16;
            end
        end
        [medianASK(j,n),~,~] = ask_2014_active(T,rup,site_ask);
        [medianBSSA(j,n),~,~] = bssa_2014_active(T,rup,site);
        [medianCB(j,n),~,~] = cb_2014_active(T,rup,site);
        [medianCY(j,n),~,~] = cy_2014_active(T,rup,site_cy);
        [medianI(j,n),~,~] = i_2014_active(T,rup,site);
    end
end

%% Figure 6
limits = [-30 40 0.05 1];
titles = ["Surface, Normal";"Buried, Normal";"Surface, Reverse";"Buried, Reverse"];

figure('Name','Gregor Figure 6','Position',[0 0 600 900],'NumberTitle','off')
for n = 1:4
    subplot(2,2,n)
    semilogy(Rx,medianASK(n,:),'-r',Rx, medianBSSA(n,:),'-g',Rx, medianCB(n,:),'-b',Rx, medianCY(n,:),'-m',Rx, medianI(n,:),'-c')
    grid on 
    axis(limits)
    xlabel('Rx')
    ylabel("PGA [g]")
    title(titles(n))
    if n == 4
        legend('ASK','BSSA','CB','CY','I','Location','South')
    end
end
%% Save Figure
saveas(gcf,'../figures/gregor6.jpg')