% Script to recreate Figure 3 from Gregor et. al. (2014)
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
Vs30 = 760;
fvs30 = 0;      % 0 for inferred and 1 for measured
Z25 = 0.6068;   % Used in CB model 
Z10 = [0.0481 0 0 0.0413 0];    % Values corresponding with the 5 models, 0 if not used (CB and I)
Zbot = 15;      
region = 2;     %   = 2 for Japan

rupt = rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
sitevar = site(is_soil,Vs30,fvs30,Z25,[],Zbot,region);

period = [0 0.2 1 3];  % Periods for the four plots
%% Function Calls
medianvec = cell(1,length(nga2west));
for j = 1:4
    T = period(j);
    for n = 1:length(M)
        if M(n)>=7
            rupt.Zhyp = 12;  
        end
        rupt.M = M(n); 
        rupt.Ztor = (max(2.673-1.136*max(rupt.M-4.970,0),0))^2;
        for i = 1:length(nga2west)
            if i ~=5
            medianvec{i}(j,n) = active_gmms(T,rupt,sitevar,gmm_name{nga2west(i)});
            else
                if M(n) >= 5
                    medianvec{i}(j,n) = active_gmms(T,rupt,sitevar,gmm_name{nga2west(i)});
                end
            end
        end
    end
end
%% Figure 3
limits = [3 8.5 0.0001 1];
titles = ["PGA" "T = 0.2s" "T = 1.0s" "T = 3.0s"];
gregor_color = ['r','g','b','m','c'];
figure('Name','Gregor Figure 3','Position',[0 0 600 900],'NumberTitle','off')
for n = 1:4
    subplot(2,2,n)
    for i = 1:length(nga2west)
        semilogy(M,medianvec{i}(n,:),'color',gregor_color(i));
        hold on
    end
%     semilogy(M, medianASK2(n,:),'-r',M, medianBSSA2(n,:),'-g',M, medianCB2(n,:),'-b',M, medianCY2(n,:),'-m',M, medianI2(n,:),'-c')
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
        legend(gmm_name{nga2west},'Location','Northwest','Interpreter','none')
    end
end

%% Save Figure
saveas(gcf,'../figures/gregor3.jpg')