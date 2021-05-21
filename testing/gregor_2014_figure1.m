% Script to recreate Figure 1 from Gregor et. al. (2014)
% Using the following five GMMs
%   ask_2014_active
%   bssa_2014_active
%   cb_2014_active
%   cy_2014_active
%   i_2014_active
% Created by Emily Mongold, 12/18/20
%
% clear
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
Rrup = [];      
Rjb = [];        
Rhyp = [];      
Rx = 10;        
Ry0 = [];       
HW = 0;             % Hanging wall indicator = 1 for hanging wall, 0 otherwise
AS = 0;             % Aftershock indicator = 1 for aftershock, 0 otherwise
Ztor = [6 3 1 0];   % Different for each magnitude earthquake
Zhyp = 8;
h_eff = [];
W = [];
delta = 90;
lambda = 0;         % Strike-slip

% Site inputs
is_soil = [];                   % Not used in these models
Vs30 = 760;
fvs30 = 0;                      % 0 for inferred and 1 for measured
Z25 = 0.6068;                   % Used in CB model 
Z10 = [0.0481 0 0 0.0413 0];    % Values corresponding with the 5 models, 0 if not used (CB and I)
Zbot = 10;      
region = 0;                     % = 0 for global

% Create rupture and site objects
rupt =       rup([],R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
sitevar =    site(is_soil,Vs30,fvs30,Z25,[],Zbot,region);

T =0;  % PGA calculation
rjb_vec = 1:1:250;
%% Function Calls
medianvec = cell(1,length(nga2west));
for j = 1:4
    rupt.M = M(j);
    rupt.Ztor = Ztor(j);
    for n = 1:length(rjb_vec)
        rupt.Rjb = rjb_vec(n);
        rupt.Rrup = sqrt(rupt.Rjb^2 + rupt.Ztor^2);
        for i = 1:length(nga2west)
            [medianvec{i}(j,n),~,~] = active_gmms(T,rupt,sitevar,gmm_name{nga2west(i)});
        end
    end
end
%% Figure 1
limits = [1 250 0.001 1];
gregor_color = ['r','g','b','m','c'];
titles = ["Mag = 5";"Mag = 6";"Mag = 7";"Mag = 8"];
figure('Name','Gregor Figure 1','Position',[0 0 600 900],'NumberTitle','off')
for n = 1:4
    subplot(2,2,n)
    for i = 1:length(nga2west)
        loglog(rjb_vec,medianvec{i}(n,:),'color',gregor_color(i));
        hold on
    end
    grid on 
    axis(limits)
    xlabel('R_{JB} Distance [km]')
    ylabel('PGA [g]')
    title(titles(n))
    if n == 4
        legend(gmm_name{nga2west},'Location','Southwest','Interpreter','none')
    end
end

%% Save Figure
saveas(gcf,'../figures/gregor1.jpg')
