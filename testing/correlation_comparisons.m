% Plots of correlation models for testing
% Created by Jack Baker
% 5/7/2021

clear; close all; clc;
addpath('../correlations/')

%% setup for calculations 

% specify periods of interest
Periods = [0.010,0.020,0.022,0.025,0.029,0.030,0.032,0.035,0.036,0.040,0.042,0.044,0.045,0.046,0.048,0.050,0.055,0.060,0.065,0.067,0.070,0.075,0.080,0.085,0.090,0.095,0.10,0.11,0.12,0.13,0.133,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.22,0.24,0.25,0.26,0.28,0.29,0.30,0.32,0.34,0.35,0.36,0.38,0.40,0.42,0.44,0.45,0.46,0.48,0.50,0.55,0.60,0.65,0.667000000000000,0.70,0.75,0.80,0.85,0.90,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.5,2.6,2.8,3,3.2,3.4,3.5,3.6,3.8,4,4.2,4.4,4.6,4.8,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10];

% specify models of interest
corr_name{1}  = 'a_2011_corr';
corr_name{2}  = 'asa_2014_corr';
corr_name{3}  = 'bb_2017_corr';
corr_name{4}  = 'bc_2006_corr';
corr_name{5}  = 'bj_2008_corr';
corr_name{6}  = 'ga_2009_corr';

% specify line styles for plotting
linespec{1}  = '-';
linespec{2}  = '-.';
linespec{3}  = '--';
linespec{4}  = '-';
linespec{5}  = '-.';
linespec{6}  = '-';

% define colors for plotting
colorspec{1} = [56 95 150]/255;
colorspec{2} = [207 89 33]/255;
colorspec{3} = [0 0 0]/255;
colorspec{4} = [231 184 0]/255;
colorspec{5} = [128 0 0]/255;
colorspec{6} = [158 184 219]/255;

% make strings for labeling the subfigures
subLabel = [{'(a)'}, {'(b)'}, {'(c)'}, {'(d)'}];


%% compute correlations

for i = 1:length(corr_name)
    rho{i} = sa_correlations(Periods, Periods, corr_name{i});
end



%% plots for a given T2

tStar = [ 0.1 0.3 1 3]; % periods to consider

figure
for j=1:length(tStar)
    tIdx = find(Periods == tStar(j));
    
    subplot(2, 2, j, 'FontSize', 9)

    for i = 1:length(corr_name)
        semilogx(Periods, rho{i}(tIdx, :), linespec{i}, 'color', colorspec{i}, 'linewidth', 2);
        hold on
    end
    
    set(gca, 'ylim', [0 1])
    axis([0.01 10 0 1])
    set(gca, 'xtick', 10.^[-2:1])
    set(gca,'xticklabel', [0.01 0.1 1 10])
    xlabel('T_1 (s)');
    ylabel('\rho', 'FontSize', 9);
    text(-0.20,-0.18, subLabel{j} ,'Units', 'Normalized', 'VerticalAlignment', 'Top')
        
    if j==2
        legend(corr_name, 'location', 'southwest', 'interpreter', 'none', 'FontSize', 8);
    end
end

% fix figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [7 7]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 7 7]);

% save figure
print('-dpdf', ['../figures/compareCorrelation.pdf']); % save the figure to a file

%% contour plots


labelSize = 10;

for i = 1:length(rho)
    v = [0.1:0.1:1];
    
    Tmin = Periods(find(~isnan(rho{i}(30,:)), 1 ));
    Tmax = Periods(find(~isnan(rho{i}(30,:)), 1, 'last' ));
    

    
    figure
    [C,h] = contour(Periods, Periods, rho{i}, v, 'k');
    hold on
    plot([Tmin Tmax Tmax Tmin Tmin], [Tmin Tmin Tmax Tmax Tmin], ':k')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca,'XTick',[0.01;0.1;1;10])
    set(gca,'YTick',[0.01;0.1;1;10])
    set(gca,'XTickLabel',[0.01;0.1;1;10])
    set(gca,'YTickLabel',[0.01;0.1;1;10])
    xlabel('T_1', 'fontsize', labelSize);
    ylabel('T_2', 'fontsize', labelSize);
    title(corr_name{i}, 'fontsize', labelSize, 'interpreter', 'none');
    set(gca, 'fontsize', labelSize)
    axis square
    clabel(C,h, 'fontsize', 7); %add labels to the contours
    
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [4.5 4]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 4.5 4]);

    print('-dpdf', ['../figures/correlationContours' num2str(i) '.pdf']); % save the figure to a file
end
