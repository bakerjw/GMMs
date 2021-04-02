% Script to create figures for the vertical and duration scripts
% Created by Emily Mongold 2/23/2021
%
clear
clc
addpath('../gmms/')
%% Setting rupture and site object values
% Rupture inputs
M = [5 6 7 8];           
R = 1:1:250;        % Used in as_1996_duration
Rrup = 1:1:250;     % Used in ab_2006_stable and c_ and as_1997_vert
Rjb = 1:1:250;      % Used in sp_2016_stable
Rhyp = [];      
Rx = [];        
Ry0 = [];       
HW = 0;             % No hanging wall, used by as_1997_vert
AS = [];
Ztor = [];       
Zhyp = [];
h_eff = [];
W = [];
delta = [];
lambda = 0;         % Strike-slip, used by vertical models

% Site inputs
is_soil = 0;        % Used by as_1996_duration, c_ and as_1997_vert
Vs30 = 760;         % Used by ab_2006_stable
fvs30 = [];          
Z25 = [];       
Z10 = [];          
Zbot = [];      
region = [];         

% Create rupture and site objects
rup = rup([],[],[],[],Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
site = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

input_durations = [1 2 3 4];    % Used by as_1996_duration
D = 80;                         % Used by c_1997_vert
T = 1.0; 

%% Creation of Duration Figure
median = zeros(length(R),length(M),length(input_durations));
for n = 1:length(M)
    rup.M = M(n);
    for m = 1:length(input_durations)
        dur_type = input_durations(m);
        for j = 1:length(R)
            rup.R = R(j);
            [median(j,n,m), sigma] = as_1996_duration(rup, site, dur_type);
        end
    end
end
% Create Duration Figure
% limits = [1 250 0.001 1];
figure('Name','Significant Durations','NumberTitle','off','Position',[10 10 600 600])
titles = ["M5.0","M6.0","M7.0","M8.0"];
for n = 1:length(M)
    subplot(3,2,n)
    plot(R,median(:,n,1),R,median(:,n,2),R,median(:,n,3),R,median(:,n,4),'LineWidth',1);
    title(titles(n));
    grid on 
    % axis(limits)
    xlabel('Distance [km]')
    ylabel('Median Duration [s]')
end
lgd = legend("5-75% horizontal significant duration","5-75% vertical significant duration","5-95% horizontal significant duration","5-95% vertical significant duration","Location","southoutside");
hL = subplot(3,2,5:6);
posL = get(hL,'position');
set(lgd,'position',posL);
axis(hL,'off');
hold off
% Save Figure
saveas(gcf,'../figures/Duration Figure.jpg')

%% Creation of Vertical Figures
median = zeros(length(Rrup),length(M),2);
for n = 1:length(M)
    rup.M = M(n);
    for j = 1:length(Rrup)
        rup.Rrup = Rrup(j);
        [median(j,n,1), sigma(j,n,1),period(1,:)] = as_1997_vert(T,rup,site); % T = 1.0s
        [median(j,n,2), sigma(j,n,2),period(2,:)] = c_1997_vert(T,rup,site,D); % T = 1.0s
    end
end
% Create Vertical Figure 1
figure('Name','Vertical Figure 1','NumberTitle','off','Position',[10 10 600 400])
titles = ["M5.0","M6.0","M7.0","M8.0"];
limits = [1 250 0.001 1];
for n = 1:length(M)
    subplot(2,2,n)
    loglog(R,median(:,n,1),R,median(:,n,2),'LineWidth',1);
    title(titles(n));
    grid on
    axis(limits)
    xlabel('Distance [km]')
    ylabel('PSA [g]')
    if n == 4
        legend("as1997","c1997","Location","NorthEast");
    end
end
hold off
% Save Figure
saveas(gcf,'../figures/Vertical Figure 1.jpg')

% Function calls for Figure 2
median = zeros(length(M),length(T),2);
T = 0.05:0.01:4;
rup.Rrup = 100; rup.Rjb = 100; 
for n = 1:length(M)
    rup.M = M(n);
        for i = 1:length(T)
            [median(n,i,1), sigma(n,i,1),period(1,:)] = as_1997_vert(T(i),rup,site);
            [median(n,i,2), sigma(n,i,2),period(2,:)] = c_1997_vert(T(i),rup,site,D);
        end
end
% Create Vertical Figure 2
figure('Name','Vertical Figure 2','NumberTitle','off','Position',[10 10 600 400])
titles = ["M5.0","M6.0","M7.0","M8.0"];
limits = [0.05 4 0.001 0.2];
for n = 1:length(M)
    subplot(2,2,n)
    loglog(T,median(n,:,1),T,median(n,:,2),'LineWidth',1);
    title(titles(n));
    grid on 
    axis(limits)
    xlabel('Period [s]')
    ylabel('PSA [g]')
    if n == 4
        legend("as1997","c1997","Location","NorthEast");
    end
end
hold off
% Save Figure
saveas(gcf,'../figures/Vertical Figure 2.jpg')

% Create Vertical Figure 3
figure('Name','Vertical Figure 3','NumberTitle','off','Position',[10 10 600 400])
limits = [0.5 4 0 1];
semilogx(T,sigma(3,:,1),T,sigma(3,:,2),'LineWidth',1);
title("Standard Deviation vs Period for M7.0");
grid on 
axis(limits)
xlabel('Period [s]')
ylabel('Sigma [ln units]')
legend("as1997","c1997","Location","NorthEast");
% Save Figure
saveas(gcf,'../figures/Vertical Figure sigma.jpg')
