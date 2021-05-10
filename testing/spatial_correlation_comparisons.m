% Plot IM correlation versus separation distance
% Jack Baker
% 1/18/2016
% edited 2/6/2018
% color figure option added 2/18/2021

clear; close all; clc
addpath('../correlations/')

%% setup for calculations 

h = [0 0.001 0.5:0.5:40]; % distances of interest

% specify models of interest
corr_name{1}  = 'gh_2008_spatial_corr';
corr_name{2}  = 'hm_2019_spatial_corr';
corr_name{3}  = 'jb_2009_spatial_corr';
corr_name{4}  = 'lb_2013_spatial_corr';
corr_name{5}  = 'mcb_2018_spatial_corr';

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
subLabel = [{'(a)'}, {'(b)'}, {'(c)'}, {'(d)'}, {'(e)'}, {'(f)'}];


%% compute correlations

T = [0.1 0.3 1 5]; % period of interest

for i = 1:length(corr_name)
    for j = 1:length(T)
        rho{i,j} = spatial_correlations(T(j), h, corr_name{i});
    end
end

%% Plot - each period in one subfigure

figure
for j=1:length(T)
    
    subplot(2, 2, j, 'FontSize', 9)

    for i = 1:length(corr_name)
        plot(h, rho{i,j}, linespec{i}, 'color', colorspec{i}, 'linewidth', 2);
        hold on
    end
    
    set(gca, 'ylim', [0 1])
    xlabel('h [km]', 'FontSize', 9);
    ylabel('\rho', 'FontSize', 9);
    title(['T = ' num2str(T(j)) ' s'])
    text(-0.20,-0.18, subLabel{j} ,'Units', 'Normalized', 'VerticalAlignment', 'Top')        
    if j==4
        legend(corr_name, 'location', 'northeast', 'interpreter', 'none', 'FontSize', 8);
    end

end

% fix figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [7 7]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 7 7]);

% save figure
print('-dpdf', ['../figures/compareSpatialCorrelation.pdf']); % save the figure to a file

%% each model in one subfigure

% define a sequential color scheme for plotting
colorspec{1} = [179,205,227]/255;
colorspec{2} = [140,150,198]/255;
colorspec{3} = [136,86,167]/255;
colorspec{4} = [129,15,124]/255;


% make legend text - flip order and add text 
for j = 1:length(T)
    legendtext{j} = ['T = ' num2str(T(length(T) - j + 1)) ' s'];
end


figure
for i=1:length(corr_name)
    
    subplot(3, 2, i, 'FontSize', 9)

    for j = length(T):-1:1
        plot(h, rho{i,j}, linespec{j}, 'color', colorspec{j}, 'linewidth', 2);
        hold on
    end
    
    set(gca, 'ylim', [0 1])
    xlabel('h [km]', 'FontSize', 9);
    ylabel('\rho', 'FontSize', 9);
    title(corr_name{i}, 'interpreter', 'none')
    if i==1
        legend(legendtext, 'location', 'northeast', 'interpreter', 'none', 'FontSize', 8);
    end

end

% fix figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [7 7]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 7 7]);

% save figure
print('-dpdf', ['../figures/compareSpatialCorrelation2.pdf']); % save the figure to a file


