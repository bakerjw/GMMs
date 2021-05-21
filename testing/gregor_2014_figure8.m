% Script to recreate Figure 8 from Gregor et. al. (2014)
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
M = [5 6 7 8];           
R = [];     
Rrup = 10;      
Rjb = 10;        
Rhyp = [];      
Rx = 10;        
Ry0 = [];       
HW = 0;             % Hanging wall indicator = 1 for hanging wall, 0 otherwise
AS = 0;             % Aftershock indicator = 1 for aftershock, 0 otherwise
Ztor = [6 3 1 0];   % Different for each magnitude earthquake
Zhyp = 8;
h_eff = [];
W = [];
delta = 15;
lambda = 0;         % Strike-slip

% Site inputs
is_soil = [];       % Not used in these models
Vs30 = 760;
fvs30 = 0;          % 0 for inferred and 1 for measured
Z25 = 1.9826;       % Used in CB model 
Z10 = [0.4704 0 0 0.4794 0];    % Values corresponding with the 5 models, 0 if not used (CB and I)
Zbot = 10;      
region = 0;         % = 0 for global


% Create rupture and site objects
rupt = rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
sitevar = site(is_soil,Vs30,fvs30,Z25,[],Zbot,region);

T_vec = 0.01:0.01:10;  % Independent variable
%% Function Calls
medianvec = cell(1,length(nga2west));
for j = 1:4
    rupt.M = M(j);
    rupt.Ztor = Ztor(j); 
    for n = 1:length(T_vec)
        for i = 1:length(nga2west)
            sitevar.Z10 = Z10(i);
            [medianvec{i}(j,n),~,~] = active_gmms(T_vec(n),rupt,sitevar,gmm_name{nga2west(i)});
        end
    end
end
%% Figure 8
limits = [0.01 10 0.001 1];
gregor_color = ['r','g','b','m','c'];
titles = ["Mag = 5.0";"Mag = 6.0";"Mag = 7.0";"Mag = 8.0"];
figure('Name','Gregor Figure 8','Position',[0 0 600 900],'NumberTitle','off')
for n = 1:4
    subplot(2,2,n)
    for i = 1:length(nga2west)
        loglog(T_vec,medianvec{i}(n,:),'color',gregor_color(i));
        hold on
    end
    grid on 
    axis(limits)
    xlabel('Period [sec]')
    ylabel('PSA [g]')
    title(titles(n))
    if n == 3
        legend(gmm_name{nga2west},'Location','Southwest','Interpreter','none')
    end
end
%% Save Figure
saveas(gcf,'../figures/gregor8.jpg')