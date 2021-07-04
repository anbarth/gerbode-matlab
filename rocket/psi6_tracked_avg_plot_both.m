% To plot grain and granule data together, you first must run
% psi6_tracked_avg_plot_grain and psi6_tracked_avg_plot_granule,
% which will set up the correct variables in your workspace

% Script for plotting the average psi6 phases of the grain and granules
% Written by Cora Payne '23, Inq Soncharoen '23, Maya Martinez '20

figure('Name','Average psi6 phase by grain and granule');

numG=2;
for g=1:numG
    errorbar(psi6_tracked_avg_grain_plotting(1:22,2*numG+1), ...
        psi6_tracked_avg_grain_plotting(1:22,g), psi6_tracked_avg_grain_plotting(1:22,g+numG), ...
       '-o', 'MarkerSize',6, 'MarkerEdgeColor',[grain_rgb(g,1) grain_rgb(g,2) grain_rgb(g,3)],...
    'MarkerFaceColor',[grain_rgb(g,1) grain_rgb(g,2) grain_rgb(g,3)],...
    'Color', [grain_rgb(g,1) grain_rgb(g,2) grain_rgb(g,3)], 'LineWidth', 3);
    hold on;
end


numG=10;
for g=1:numG
    errorbar(psi6_tracked_avg_granule_plotting(1:22,2*numG+1), ...
        psi6_tracked_avg_granule_plotting(1:22,g), psi6_tracked_avg_granule_plotting(1:22,g+numG), ...
        ':o','MarkerSize',4, 'MarkerEdgeColor',[granule_rgb(g,1) granule_rgb(g,2) granule_rgb(g,3)],...
    'MarkerFaceColor',[granule_rgb(g,1) granule_rgb(g,2) granule_rgb(g,3)],...
    'Color', [granule_rgb(g,1) granule_rgb(g,2) granule_rgb(g,3)]);
    hold on;
end

%legend('Grain 1', 'Grain 2', 'Granule 1', 'Granule 2', 'Granule 3', 'Granule 4', 'Granule 5', 'Granule 6', 'Granule 7', 'Granule 8', 'Granule 9', 'Granule 10');
xlabel('Time Elapsed (sec)');
ylabel('Average \Psi_6 phase');