% Script to recreate Figure 4 from Gregor et. al. (2014)
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
M = 7;  
R = [];     
Rrup = 100;      
Rjb = 100;        
Rhyp = [];      
Rx = 100;        
Ry0 = [];       
HW = 0;         % Hanging wall indicator = 1 for hanging wall, 0 otherwise
AS = 0;         % Aftershock indicator = 1 for aftershock, 0 otherwise
Ztor = 1;
Zhyp = 8;
h_eff = [];
W = [];
delta = 90;
lambda = 0;     % Strike-slip

% Site inputs
is_soil = [];   % Not used in these models
Vs30 = 190:10:1500;
fvs30 = 0;      % 0 for inferred and 1 for measured
Z25 = [];       % Used in CB model 
Z10 = [];       % Values change with Vs30: left empty so each GMM fills in accordingly
Zbot = 15;      
region = 0;     % 0 for global
     
rupt = rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
sitevar = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

period = [0 0.2 1 3];  % Periods for the four plots
%% Function Calls
medianvec = cell(1,length(nga2west));
for j = 1:length(period)
    T = period(j);
    for n = 1:length(Vs30)
        sitevar.Vs30 = Vs30(n);
        for i = 1:length(nga2west)
            if i ~=5
            [medianvec{i}(j,n),~,~] = active_gmms(T,rupt,sitevar,gmm_name{nga2west(i)});
            else
                if Vs30(n) >= 450
                    [medianvec{i}(j,n),~,~] = active_gmms(T,rupt,sitevar,gmm_name{nga2west(i)});
                end
            end
        end
    end
end

%% Figure 4
limits = [100 2000 0.001 0.1];
gregor_color = ['r','g','b','m','c'];
y_labels = ["PGA [g]";"PSA [g]";"PSA [g]";"PSA [g]"];
titles = ["PGA";"T = 0.2s";"T = 1.0s";"T = 3.0s"];

figure('Name','Gregor Figure 4','Position',[0 0 600 900],'NumberTitle','off')
format shortG
for n = 1:4
    if n >2
        limits = [100 2000 0.001 0.1];
    else
        limits = [100 2000 0.01 1];
    end
    subplot(2,2,n)
        for i = 1:length(nga2west)
            loglog(Vs30,medianvec{i}(n,:),'color',gregor_color(i));
            hold on
        end
    grid on 
    axis(limits)
    xlabel('Vs30')
    ylabel(y_labels(n))
    title(titles(n))
    if n == 4
        legend(gmm_name{nga2west},'Location','NorthEast','Interpreter','none')
    end
end
%% Save Figure
saveas(gcf,'../figures/gregor4new.jpg')
% save("ASK_SAs.mat","medianASK2")
% save("BSSA_SAs.mat","medianBSSA2")
% save("CB_SAs.mat","medianCB2")
% save("CY_SAs.mat","medianCY2")
% save("I_SAs.mat","medianI2")