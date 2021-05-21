% Script to recreate Figure 6 from Gregor et. al. (2014)
% Using the following five GMMs
%   ask_2014_active
%   bssa_2014_active
%   cb_2014_active
%   cy_2014_active
%   i_2014_active
% Created by Emily Mongold, 12/18/20
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
Z10 = [0.0481 0 0 0.0413 0];    % Values corresponding with the 5 models, 0 if not used (CB and I)
Zbot = 15;      
region = 0;     % Global
        
rupt = rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
sitevar = site(is_soil,Vs30,fvs30,Z25,[],Zbot,region);

T = 0;      % PGA
%% Function Calls
medianvec = cell(1,length(nga2west));

for j = 1:4
    % Surface or Buried
    if j==2 || j==4
        rupt.Ztor = 6;
    else 
        rupt.Ztor = 0;
    end
    % Normal or Reverse Faulting
    if j <= 2
        rupt.lambda = -90;
    else
        rupt.lambda = 90;
    end
    for n = 1:length(Rx)
        rupt.Rx = Rx(n); 
        if Rx(n)<0
            rupt.Rrup = sqrt(Rx(n)^2 + rupt.Ztor^2); rupt.Rjb = abs(Rx(n)); 
        elseif Rx(n)>=0   
            rupt.Rrup = (Rx(n) + rupt.Ztor)*cos(rupt.delta);
            if Rx(n)<=16
                rupt.Rjb = 0; 
            else
                rupt.Rjb = Rx(n)-16;
            end
        end
        for i = 1:length(nga2west)
            [medianvec{i}(j,n),~,~] = active_gmms(T,rupt,sitevar,gmm_name{nga2west(i)});
        end
    end
end

%% Figure 6
limits = [-30 40 0.05 1];
gregor_color = ['r','g','b','m','c'];
titles = ["Surface, Normal";"Buried, Normal";"Surface, Reverse";"Buried, Reverse"];
figure('Name','Gregor Figure 6','Position',[0 0 600 900],'NumberTitle','off')
for n = 1:4
    subplot(2,2,n)
    for i = 1:length(nga2west)
        semilogy(Rx,medianvec{i}(n,:),'color',gregor_color(i));
        hold on
    end
    grid on 
    axis(limits)
    xlabel('Rx')
    ylabel("PGA [g]")
    title(titles(n))
    if n == 4
        legend(gmm_name{nga2west},'Location','South','Interpreter','none')
    end
end
%% Save Figure
saveas(gcf,'../figures/gregor6.jpg')