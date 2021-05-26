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
Rrup = 1:1:250;     % Used in ab_2006_stable and c_ and as_1997_vert, as_2016_duration and bsa_2009_duration
Rjb = 1:1:250;      % Used in sp_2016_stable
Rhyp = [];      
Rx = [];        
Ry0 = [];       
HW = 0;             % No hanging wall, used by as_1997_vert
AS = [];
Ztor = 0;           % Used in bsa_2009_duration
Zhyp = [];
h_eff = [];
W = [];
delta = [];
lambda = 0;         % Strike-slip, used by vertical models

% Site inputs
is_soil = 0;        % Used by as_1996_duration, c_ and as_1997_vert
Vs30 = [760 270];         % Used by ab_2006_stable, as_2016_duration, and bsa_2009_duration
fvs30 = [];          
Z25 = [];       
Z10 = 0.05;         % Used in as_2016_duration
Zbot = [];      
region = 1;         % Used in as_2016_duration

% Create rupture and site objects
rup = rup([],[],[],[],Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda);
site = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region);

input_durations = [1 3 2 4 5];      % Used by as_1996_duration
D = 80;                             % Used by c_1997_vert
T = 1.0; 

%% Creation of Duration Figures
% median = zeros(3,length(R),length(M),length(input_durations));
% sigma = zeros(3,length(R),length(M),length(input_durations));
% med_vert = zeros(3,length(R),length(M),length(input_durations));
% sigma_vert = zeros(3,length(R),length(M),length(input_durations));
% med_2080 = zeros(3,length(R),length(M),length(input_durations));
% sigma_2080 = zeros(3,length(R),length(M),length(input_durations));
for i = 1:length(Vs30)
    site.Vs30 = Vs30(i); 
    for n = 1:length(M)
        rup.M = M(n);
        for m = 1:length(input_durations)
            dur_type = input_durations(m);
            if dur_type == 1 || dur_type == 3
                for j = 1:length(R)
                    rup.R = R(j); rup.Rrup = R(j);
                    [median(1,j,n,m,i), sigma(1,j,n,m,i)] = as_1996_duration(rup, site, dur_type);
                    [median(2,j,n,m,i), sigma(2,j,n,m,i),~,~] = as_2016_duration(rup, site, dur_type);
                    [median(3,j,n,m,i), sigma(3,j,n,m,i),~,~] = bsa_2009_duration(rup, site, dur_type);
                end
            elseif dur_type == 2 || dur_type == 4
                for j = 1:length(R)
                    rup.R = R(j); rup.Rrup = R(j);
                    [med_vert(j,n,m-2,i), sigma_vert(j,n,m-2,i)] = as_2016_duration(rup, site, dur_type);
                end
            elseif dur_type == 5
                for j = 1:length(R)
                    rup.R = R(j);rup.Rrup = R(j);
                    [med_2080(j,n,m-4,i), sigma_2080(j,n,m-4,i)] = as_2016_duration(rup, site, dur_type);
                end
            end
        end
    end
end
% Create Duration Figure 1
color1 = [0,0.4470,0.7410]; color2 = [0.85,0.3250,0.0980]; 
color3 = [0.929,0.694,0.1250]; color4 = [0.494,0.184,0.556]; 
% limits = [1 250 0.001 1];
figure('Name','Significant Durations','NumberTitle','off','Position',[10 10 600 600])
titles = ["M5.0","M6.0","M7.0","M8.0"];
for n = 1:length(M)
    subplot(3,2,n)
    plot(R,med_vert(:,n,1,1),'-','color',color1,'Linewidth',1);
    hold on
    plot(R,med_vert(:,n,2,1),'--','color',color1,'Linewidth',1);
    title(titles(n));
    grid on 
    % axis(limits)
    xlabel('Distance [km]')
    ylabel('Median Duration [s]')
    %sgtitle('Vertical Significant Durations from as1996')
end
lgd = legend("5-75% vertical duration","5-95% vertical duration","Location","southoutside");
hL = subplot(3,2,5:6);
posL = get(hL,'position');
set(lgd,'position',posL);
axis(hL,'off');
hold off
% set figure size
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 6]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 6 6]);
% Save Figure
saveas(gcf,'../figures/Duration Figure 1.pdf')

% Create Duration Figure 2
figure('Name','Significant Durations','NumberTitle','off','Position',[10 10 600 600])
titles = ["M5.0, Vs30 = 760m/s","M5.0, Vs30 = 270m/s","M7.0, Vs30 = 760m/s","M7.0, Vs30 = 270m/s"];
for n = [1 3]
    for i = 1:length(Vs30)
        subplot(3,2,(n+i-1))
        plot(R,median(1,:,n,1,i),'-','color',color1,'Linewidth',1);
        hold on
        plot(R,median(1,:,n,2,i),'--','color',color1,'Linewidth',1);
        plot(R,median(2,:,n,1,i),'-','color',color2,'Linewidth',1);
        plot(R,median(2,:,n,2,i),'--','color',color2,'Linewidth',1);
        plot(R,median(3,:,n,1,i),'-','color',color3,'Linewidth',1);
        plot(R,median(3,:,n,2,i),'--','color',color3,'Linewidth',1);

        title(titles(n+i-1));
        grid on 
        xlabel('Distance [km]')
        ylabel('Median Duration [s]')
    end
end
lgd = legend("5-75% as1996","5-95% as1996","5-75% as2016","5-95% as2016","5-75% bsa2009","5-95% bsa2009","Location","southoutside");
hL = subplot(3,2,5:6);
posL = get(hL,'position');
set(lgd,'position',posL);
axis(hL,'off');
hold off
% set figure size
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 6]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 6 6]);
% Save Figure
saveas(gcf,'../figures/Duration Figure 2.pdf')

%% Creation of Vertical Figures
median = zeros(length(Rrup),length(M),2);
site.Vs30 = 760; 
M = [5 6 7 8]; 
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
    ylabel('SA [g]')
    if n == 4
        legend("as1997","c1997","Location","Southwest");
    end
end
hold off
% Save Figure
% set figure size
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 6]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 6 6]);
saveas(gcf,'../figures/Vertical Figure 1.pdf')

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
    ylabel('SA [g]')
    if n == 4
        legend("as1997","c1997","Location","Southwest");
    end
end
hold off
% set figure size
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 6]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 6 6]);
% Save Figure
saveas(gcf,'../figures/Vertical Figure 2.pdf')

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
% set figure size
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [4.5 4]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 4.5 4]);
% Save Figure
saveas(gcf,'../figures/Vertical Figure sigma.pdf')
