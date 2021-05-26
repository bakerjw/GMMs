% Script to recreate Figure 9 from Gregor et. al. (2014)
% Using the following five GMMs
%   ask_2014_active
%   bssa_2014_active
%   cb_2014_active
%   cy_2014_active
%   i_2014_active
% Created by Emily Mongold, 12/18/20
%
% NOTE: This file does not recreate the figure as it appears in the text
% for BSSA2014, Vs30 = 270, default Z10,Z25. The outputs of that script do 
% not return the same values as the spreadsheet, but the error was not
% identified as of 21 Jan 2021.
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
M = 7;           
R = [];          
Rjb = 10;        
Rhyp = [];      
Rx = 10;        
Ry0 = [];       
HW = 0;             % Hanging wall indicator = 1 for hanging wall, 0 otherwise
AS = 0;             % Aftershock indicator = 1 for aftershock, 0 otherwise
Ztor = 1;    
Rrup = sqrt(Rjb^2+Ztor^2); 
Zhyp = 8;
h_eff = [];
W = [];
delta = 90;
lambda = 0;     % Strike-slip

% Site inputs
is_soil = [];
Vs30 = [270 270 270 760];       % Changes for each subplot, same for all GMMs
fvs30 = 0;          % 0 for inferred
Z25 = [1.9826 0.9 4.8 0.6068];  % Changes for each subplot, only used in CB
Z10 = [0.4704 0 0 0.4794 0;
    0.1 0.1 0.1 0.1 0.1;
    1.2 1.2 1.2 1.2 1.2;
    0.0481 0 0 0.0413 0];       % Changes for each subplot and each GMM
Zbot = 10;      
region = 0;             % = 0 for global

% Create rupture and site objects
rupt = rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,[],Zhyp,h_eff,W,delta,lambda);
sitevar = site(is_soil,Vs30,fvs30,[],[],Zbot,region);

T_vec = 0.01:0.01:10;  % Independent Variable
%% Function Calls
medianvec = cell(1,length(nga2west));
for j = 1:4
    sitevar.Vs30 = Vs30(j);
    sitevar.Z25 = Z25(j); 
    for n = 1:length(T_vec)
        for i = 1:length(nga2west)
            if i ~= 5
            sitevar.Z10 = Z10(j,i); 
            [medianvec{i}(j,n),~,~] = active_gmms(T_vec(n),rupt,sitevar,gmm_name{nga2west(i)});
            elseif i == 5 && j == 4
                [medianvec{i}(j,n),~,~] = active_gmms(T_vec(n),rupt,sitevar,gmm_name{nga2west(i)});
            end
        end
    end
end
%% Figure 9
limits = [0.01 10 0.01 1];
gregor_color = ['r','g','b','m','c'];
titles = ["Vs30 = 270m/sec, Default Z1, Z25";"Vs30 = 270m/sec, Z1=0.1km, Z25=0.9km";"Vs30 = 270m/sec, Z1=1.2km, Z25=4.8km";"Vs30 = 760m/sec, Default Z1, Z25"];
figure('Name','Gregor Figure 9','Position',[0 0 600 900],'NumberTitle','off')
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

% set figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [5 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 5 6]);


%% Save Figure
saveas(gcf,'../figures/gregor9.pdf')