
check_out_origin                    = 0;

do_cop_plots_single                 = 0;
do_cop_plots_absolute               = 1;
do_step_parameter_plots             = 1;

save_figures                        = 0;

load subjectInfo.mat;
load(makeFileName(date, subject_id, 'gaitParametersConditions'));
load(makeFileName(date, subject_id, 'gaitParametersForceplate'));
load(makeFileName(date, subject_id, 'gaitParametersMarker'));

color_control = [0.3 0.1 1];
color_right = [1 0.3 .1] * 0.7;
color_left = [0.3 1 0.1] * 0.7;



%% check out origin
if check_out_origin
    condition = conditions_stanceL;
    
    data = lcop_x_normalized_total(:, condition);
    time_data = repmat(time_normalized', 1, size(data, 2));
    index_data = repmat((1 : size(data, 2)), length(time_normalized), 1);
    
    origin_trial_list_condition = origin_trial_list_total(condition);
    origin_start_time_list_condition = origin_start_time_list_total(condition);
    origin_end_time_list_condition = origin_end_time_list_total(condition);
    
    figure; axes; hold on; title('origin'); set(gca, 'Fontsize', 12)
    plot3(time_data, data, index_data)
return    
 
end

%% cop plots
if do_cop_plots_single
    figure; axes; hold on; title(['left foot medial-lateral CoP, N = ' num2str(sum(conditions_stanceL))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['right foot medial-lateral CoP, N = ' num2str(sum(conditions_stanceR))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR), 2), 'color', color_control, 'linewidth', 5);
end

if do_cop_plots_absolute
    figure; axes; hold on; title('right foot medial-lateral CoP, 0ms'); set(gca, 'Fontsize', 12)
    left_plot = shadedErrorBar(time_normalized, lcop_x_mean_stanceL, lcop_x_std_stanceL, {'color', color_left, 'linewidth', 5}, 1);
    right_plot = shadedErrorBar(time_normalized, rcop_x_mean_stanceR, rcop_x_std_stanceR, {'color', color_right, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
    this_legend = legend([left_plot.mainLine right_plot.mainLine], 'left CoP, mean \pm STD', 'right CoP, mean \pm STD');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'cop_absolute.eps', 'epsc2')
    end
    
end

%% step parameter plots
if do_step_parameter_plots
    step_width_figure = figure; axes('fontsize', 12); hold on; title('step width', 'fontsize', 16);
    step_width_left_histogram = histogram(step_width_total(conditions_stanceR), 'Normalization', 'pdf', 'displaystyle', 'bar', 'edgecolor', color_left, 'facecolor', lightenColor(color_left, 0.5));
    step_width_right_histogram = histogram(step_width_total(conditions_stanceL), 'Normalization', 'pdf', 'displaystyle', 'bar', 'edgecolor', color_right, 'facecolor', lightenColor(color_right, 0.5));
    
    step_length_figure = figure; axes('fontsize', 12); hold on; title('step length', 'fontsize', 16);
    step_length_left_histogram = histogram(step_length_total(conditions_stanceR), 'Normalization', 'pdf', 'displaystyle', 'bar', 'edgecolor', color_left, 'facecolor', lightenColor(color_left, 0.5));
    step_length_right_histogram = histogram(step_length_total(conditions_stanceL), 'Normalization', 'pdf', 'displaystyle', 'bar', 'edgecolor', color_right, 'facecolor', lightenColor(color_right, 0.5));
    
end


