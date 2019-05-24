

% load stuff
root_gvs = '/Volumes/GoogleDrive/Team Drives/CoBaL/Pilots/cognitiveLoadGvs';
root_vis = '/Volumes/GoogleDrive/Team Drives/CoBaL/Pilots/cognitiveLoadVision';

results_gvs = load([root_gvs filesep 'COGGP01/countingResults.mat']);
results_vis = load([root_vis filesep 'COGVP01/countingResults.mat']);
plot_settings_gvs = SettingsCustodian([root_gvs filesep 'plotSettings.txt']); 
plot_settings_vis = SettingsCustodian([root_vis filesep 'plotSettings.txt']); 
colors_gvs = plot_settings_gvs.get('colors_comparison');
colors_vis = plot_settings_vis.get('colors_comparison');


% define stuff
linewidth = 3;
shift_gvs = -0.1;
shift_vis = 0.1;
box_width = 1.5;

% count number figure
count_figure = figure; count_axes = axes; hold on
ylim([47, 61]);
plot((1:10)+shift_vis, results_vis.numbers_counted(3:end), 'd-', ...
    'linewidth', linewidth, 'color', colors_vis(2, :), ...
    'markersize', 15, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors_vis(2, :) ...
    )
plot((1:10)+shift_gvs, results_gvs.numbers_counted(3:end), 'd-', ...
    'linewidth', linewidth, 'color', colors_gvs(2, :), ...
    'markersize', 15, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors_gvs(2, :) ...
    )
singleBoxPlot(results_vis.numbers_counted(3:end), ...
    'abscissa', 11.25, 'width', box_width, 'FaceColor', lightenColor(colors_vis(2, :), 0.5), 'EdgeColor', colors_vis(2, :), ...
    'MarkerColor', colors_vis(2, :), 'WiskColor', colors_vis(2, :), 'EdgeLinewidth', 3, 'WiskLinewidth', 3, ...
    'FaceAlpha', 0.3, 'MeanLinewidth', 5, 'MeanColor', colors_vis(2, :), 'plotMedian', 0, 'ShowData', true, 'ShowOutliers', false ...
  )
singleBoxPlot(results_gvs.numbers_counted(3:end), ...
    'abscissa', 13, 'width', box_width, 'FaceColor', lightenColor(colors_gvs(2, :), 0.5), 'EdgeColor', colors_gvs(2, :), ...
    'MarkerColor', colors_gvs(2, :), 'WiskColor', colors_gvs(2, :), 'EdgeLinewidth', 3, 'WiskLinewidth', 3, ...
    'FaceAlpha', 0.3, 'MeanLinewidth', 5, 'MeanColor', colors_gvs(2, :), 'plotMedian', 0, 'ShowData', true, 'ShowOutliers', false ...
  )


error_figure = figure; error_axes = axes; hold on
ylim([-0.5, 3.5]);
plot((1:10)+shift_vis, results_vis.errors(3:end), 's-', ...
    'linewidth', linewidth, 'color', colors_vis(2, :), ...
    'markersize', 15, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors_vis(2, :) ...
    )
plot((1:10)+shift_gvs, results_gvs.errors(3:end), 's-', ...
    'linewidth', linewidth, 'color', colors_gvs(2, :), ...
    'markersize', 15, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors_gvs(2, :) ...
    )
singleBoxPlot(results_vis.errors(3:end), ...
    'abscissa', 11.25, 'width', box_width, 'FaceColor', lightenColor(colors_vis(2, :), 0.5), 'EdgeColor', colors_vis(2, :), ...
    'MarkerColor', colors_vis(2, :), 'WiskColor', colors_vis(2, :), 'EdgeLinewidth', 3, 'WiskLinewidth', 3, ...
    'FaceAlpha', 0.3, 'MeanLinewidth', 5, 'MeanColor', colors_vis(2, :), 'plotMedian', 0, 'ShowData', true, 'ShowOutliers', false ...
  )
singleBoxPlot(results_gvs.errors(3:end), ...
    'abscissa', 13, 'width', box_width, 'FaceColor', lightenColor(colors_gvs(2, :), 0.5), 'EdgeColor', colors_gvs(2, :), ...
    'MarkerColor', colors_gvs(2, :), 'WiskColor', colors_gvs(2, :), 'EdgeLinewidth', 3, 'WiskLinewidth', 3, ...
    'FaceAlpha', 0.3, 'MeanLinewidth', 5, 'MeanColor', colors_gvs(2, :), 'plotMedian', 0, 'ShowData', true, 'ShowOutliers', false ...
  )

% save figures
filename = '~/Dropbox/Writing/Proposals/2019_R01/figures/raw/withLabels/count.tiff';
saveas(count_figure, filename, 'tiff')
filename = '~/Dropbox/Writing/Proposals/2019_R01/figures/raw/withLabels/errors.tiff';
saveas(error_figure, filename, 'tiff')
            
% remove text and marks to save graphs only
set(get(count_axes, 'xaxis'), 'visible', 'off');
set(get(count_axes, 'yaxis'), 'visible', 'off');
set(get(count_axes, 'xlabel'), 'visible', 'off');
set(get(count_axes, 'ylabel'), 'visible', 'off');
set(get(count_axes, 'title'), 'visible', 'off');
set(count_axes, 'xticklabel', '');
set(count_axes, 'yticklabel', '');
set(count_axes, 'position', [0 0 1 1]);
filename = '~/Dropbox/Writing/Proposals/2019_R01/figures/raw/noLabels/count.tiff';
saveas(count_figure, filename, 'tiff')

set(get(error_axes, 'xaxis'), 'visible', 'off');
set(get(error_axes, 'yaxis'), 'visible', 'off');
set(get(error_axes, 'xlabel'), 'visible', 'off');
set(get(error_axes, 'ylabel'), 'visible', 'off');
set(get(error_axes, 'title'), 'visible', 'off');
set(error_axes, 'xticklabel', '');
set(error_axes, 'yticklabel', '');
set(error_axes, 'position', [0 0 1 1]);
filename = '~/Dropbox/Writing/Proposals/2019_R01/figures/raw/noLabels/error.tiff';
saveas(error_figure, filename, 'tiff')











