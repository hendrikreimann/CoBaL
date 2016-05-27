
check_out_origin                    = 0;

do_cop_plots_single                 = 0;
do_cop_plots_absolute_right         = 1;
do_cop_plots_absolute_left          = 0;
do_cop_plots_response_right         = 1;
do_cop_plots_response_left          = 0;

do_heel_plots_single                = 0;
do_heel_plots_absolute_right        = 1;
do_heel_plots_absolute_left         = 0;
do_heel_plots_response_right        = 1;
do_heel_plots_response_left         = 0;

do_phase_delayed_plots              = 1;

save_figures                        = 0;

load subjectInfo.mat;
load(makeFileName(date, subject_id, 'resultsConditions'));
load(makeFileName(date, subject_id, 'resultsForceplate'));
load(makeFileName(date, subject_id, 'resultsMarker'));

color_control = [0.3 0.1 1];
color_positive = [1 0.3 .1] * 0.7;
color_negative = [0.3 1 0.1] * 0.7;



%% check out origin
if check_out_origin
%     condition = conditions_stanceL_stimPos_0ms;
%     condition = conditions_stanceL_stimNo;
    condition = conditions_stanceL_stimNo;
    
%     data = rheel_x_pos_normalized_total(:, condition);
    data = rheel_x_pos_normalized_total(:, condition);
    
    time_data = repmat(time_normalized', 1, size(data, 2));
    index_data = repmat((1 : size(data, 2)), length(time_normalized), 1);
    
    origin_trial_list_condition = origin_trial_list_total(condition);
    origin_start_time_list_condition = origin_start_time_list_total(condition);
    origin_end_time_list_condition = origin_end_time_list_total(condition);
    
    figure; axes; hold on; title('origin finder'); set(gca, 'Fontsize', 12)
    plot3(time_data, data, index_data)
return    
 
end

%% cop plots
if do_cop_plots_single
    figure; axes; hold on; title(['left foot medial-lateral CoP, N = ' num2str(sum(conditions_stanceL_stimNo))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimNo), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['right foot medial-lateral CoP, N = ' num2str(sum(conditions_stanceR_stimNo))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimNo), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['left foot medial-lateral CoP, 0ms, N = ' num2str(sum(conditions_stanceL_0ms))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
    plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
    plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    figure; axes; hold on; title(['right foot medial-lateral CoP, 0ms, N = ' num2str(sum(conditions_stanceR_0ms))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
    plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
    plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    if do_phase_delayed_plots
        figure; axes; hold on; title(['right foot medial-lateral CoP, 150ms, N = ' num2str(sum(conditions_stanceR_150ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['left foot medial-lateral CoP, 450ms, N = ' num2str(sum(conditions_stanceL_450ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);
        
        figure; axes; hold on; title(['left foot medial-lateral CoP, 150ms, N = ' num2str(sum(conditions_stanceL_150ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['right foot medial-lateral CoP, 450ms, N = ' num2str(sum(conditions_stanceR_450ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);
    end
end

if do_cop_plots_absolute_right
    % 0ms
    figure; axes; hold on; title('right foot medial-lateral CoP, 0ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, rcop_x_mean_stanceR_control, rcop_x_civ_stanceR_control, {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, rcop_x_mean_stanceR_stimPos_0ms, rcop_x_civ_stanceR_stimPos_0ms, {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, rcop_x_mean_stanceR_stimNeg_0ms, rcop_x_civ_stanceR_stimNeg_0ms, {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'cop_absolute_stanceR_0ms.eps', 'epsc2')
    end
    
    if do_phase_delayed_plots
        % 150ms
        figure; axes; hold on; title('right foot medial-lateral CoP, 150ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, rcop_x_mean_stanceR_control, rcop_x_civ_stanceR_control, {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, rcop_x_mean_stanceR_stimPos_150ms, rcop_x_civ_stanceR_stimPos_150ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, rcop_x_mean_stanceR_stimNeg_150ms, rcop_x_civ_stanceR_stimNeg_150ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_absolute_stanceR_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title('left foot medial-lateral CoP, 450ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, lcop_x_mean_stanceL_control, lcop_x_civ_stanceL_control, {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, lcop_x_mean_stanceL_stimPos_450ms, lcop_x_civ_stanceL_stimPos_450ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, lcop_x_mean_stanceL_stimNeg_450ms, lcop_x_civ_stanceL_stimNeg_450ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_absolute_stanceL_450ms.eps', 'epsc2')
        end
    end
end

if do_cop_plots_absolute_left
    figure; axes; hold on; title('left foot medial-lateral CoP, 0ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, lcop_x_mean_stanceL_control, lcop_x_civ_stanceL_control, {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, lcop_x_mean_stanceL_stimPos_0ms, lcop_x_civ_stanceL_stimPos_0ms, {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, lcop_x_mean_stanceL_stimNeg_0ms, lcop_x_civ_stanceL_stimNeg_0ms, {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'cop_absolute_stanceL_0ms.eps', 'epsc2')
    end
    
    if do_phase_delayed_plots
        % 150ms
        figure; axes; hold on; title('left foot medial-lateral CoP, 150ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, lcop_x_mean_stanceL_control, lcop_x_civ_stanceL_control, {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, lcop_x_mean_stanceL_stimPos_150ms, lcop_x_civ_stanceL_stimPos_150ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, lcop_x_mean_stanceL_stimNeg_150ms, lcop_x_civ_stanceL_stimNeg_150ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_absolute_stanceL_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title('right foot medial-lateral CoP, 450ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, rcop_x_mean_stanceR_control, rcop_x_civ_stanceR_control, {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, rcop_x_mean_stanceR_stimPos_450ms, rcop_x_civ_stanceR_stimPos_450ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, rcop_x_mean_stanceR_stimNeg_450ms, rcop_x_civ_stanceR_stimNeg_450ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_absolute_stanceR_450ms.eps', 'epsc2')
        end
    end
end

if do_cop_plots_response_right
    % 0ms
    figure; axes; hold on; title('right foot medial-lateral CoP response, 0ms'); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, rcop_x_response_mean_stanceR_stimPos_0ms, rcop_x_response_civ_stanceR_stimPos_0ms, {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, rcop_x_response_mean_stanceR_stimNeg_0ms, rcop_x_response_civ_stanceR_stimPos_0ms, {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-0.015 0.015]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'cop_response_stanceR_0ms.eps', 'epsc2')
    end
    
    if do_phase_delayed_plots
        % 150ms
        figure; axes; hold on; title('right foot medial-lateral CoP response, 150ms'); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, rcop_x_response_mean_stanceR_stimPos_150ms, rcop_x_response_civ_stanceR_stimPos_150ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, rcop_x_response_mean_stanceR_stimNeg_150ms, rcop_x_response_civ_stanceR_stimPos_150ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-0.015 0.015]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_response_stanceR_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title('left foot medial-lateral CoP response, 450ms'); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, lcop_x_response_mean_stanceL_stimPos_450ms, lcop_x_response_civ_stanceL_stimPos_450ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, lcop_x_response_mean_stanceL_stimNeg_450ms, lcop_x_response_civ_stanceL_stimPos_450ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-0.015 0.015]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_response_stanceL_450ms.eps', 'epsc2')
        end
    end
end

if do_cop_plots_response_left
    figure; axes; hold on; title('left foot medial-lateral CoP response, 0ms'); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, lcop_x_response_mean_stanceL_stimPos_0ms, lcop_x_response_civ_stanceL_stimPos_0ms, {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, lcop_x_response_mean_stanceL_stimNeg_0ms, lcop_x_response_civ_stanceL_stimPos_0ms, {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-0.015 0.015]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'cop_response_stanceL_0ms.eps', 'epsc2')
    end
    
    if do_phase_delayed_plots
        % 150ms
        figure; axes; hold on; title('left foot medial-lateral CoP response, 150ms'); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, lcop_x_response_mean_stanceL_stimPos_150ms, lcop_x_response_civ_stanceL_stimPos_150ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, lcop_x_response_mean_stanceL_stimNeg_150ms, lcop_x_response_civ_stanceL_stimPos_150ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-0.015 0.015]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_response_stanceL_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title('right foot medial-lateral CoP response, 450ms'); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, rcop_x_response_mean_stanceR_stimPos_450ms, rcop_x_response_civ_stanceR_stimPos_450ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, rcop_x_response_mean_stanceR_stimNeg_450ms, rcop_x_response_civ_stanceR_stimPos_450ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-0.015 0.015]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_response_stanceR_450ms.eps', 'epsc2')
        end
    end
end

%% heel plots
if do_heel_plots_single
    figure; axes; hold on; title(['left foot medial-lateral heel pos, N = ' num2str(sum(conditions_stanceL_stimNo))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, rheel_x_pos_normalized_total(:, conditions_stanceL_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(rheel_x_pos_normalized_total(:, conditions_stanceL_stimNo), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['right foot medial-lateral heel pos, N = ' num2str(sum(conditions_stanceR_stimNo))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, lheel_x_pos_normalized_total(:, conditions_stanceR_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(lheel_x_pos_normalized_total(:, conditions_stanceR_stimNo), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['left foot medial-lateral heel pos, 0ms, N = ' num2str(sum(conditions_stanceR_0ms))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, lheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
    plot(time_normalized, lheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
    plot(time_normalized, mean(lheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(lheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    figure; axes; hold on; title(['right foot medial-lateral heel pos, 0ms, N = ' num2str(sum(conditions_stanceL_0ms))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, rheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
    plot(time_normalized, rheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
    plot(time_normalized, mean(rheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(rheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    if do_phase_delayed_plots
        figure; axes; hold on; title(['right foot medial-lateral heel pos, 150ms, N = ' num2str(sum(conditions_stanceR_150ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, lheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, lheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(lheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(lheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['left foot medial-lateral heel pos, 450ms, N = ' num2str(sum(conditions_stanceL_450ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, rheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(rheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(rheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);    
    end    
    
    if do_phase_delayed_plots
        figure; axes; hold on; title(['left foot medial-lateral CoP, 150ms, N = ' num2str(sum(conditions_stanceR_150ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['right foot medial-lateral CoP, 450ms, N = ' num2str(sum(conditions_stanceL_450ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);
    end
end

if do_heel_plots_absolute_right
    % 0ms
    figure; axes; hold on; title('left foot medial-lateral heel pos, 0ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, lheel_x_pos_mean_stanceR_control, lheel_x_pos_civ_stanceR_control, {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, lheel_x_pos_mean_stanceR_stimPos_0ms, lheel_x_pos_civ_stanceR_stimPos_0ms, {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, lheel_x_pos_mean_stanceR_stimNeg_0ms, lheel_x_pos_civ_stanceR_stimNeg_0ms, {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'heel_absolute_stanceR_0ms.eps', 'epsc2')
    end
    
    if do_phase_delayed_plots
        % 150ms
        figure; axes; hold on; title('left foot medial-lateral heel pos, 150ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, lheel_x_pos_mean_stanceR_control, lheel_x_pos_civ_stanceR_control, {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, lheel_x_pos_mean_stanceR_stimPos_150ms, lheel_x_pos_civ_stanceR_stimPos_150ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, lheel_x_pos_mean_stanceR_stimNeg_150ms, lheel_x_pos_civ_stanceR_stimNeg_150ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'heel_absolute_stanceR_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title('right foot medial-lateral heel pos, 450ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, rheel_x_pos_mean_stanceL_control, rheel_x_pos_civ_stanceL_control, {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, rheel_x_pos_mean_stanceL_stimPos_450ms, rheel_x_pos_civ_stanceL_stimPos_450ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, rheel_x_pos_mean_stanceL_stimNeg_450ms, rheel_x_pos_civ_stanceL_stimNeg_450ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'heel_absolute_stanceL_450ms.eps', 'epsc2')
        end
    end
end

if do_heel_plots_absolute_left
    % 0ms
    figure; axes; hold on; title('right foot medial-lateral heel pos, 0ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, rheel_x_pos_mean_stanceL_control, rheel_x_pos_civ_stanceL_control, {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, rheel_x_pos_mean_stanceL_stimPos_0ms, rheel_x_pos_civ_stanceL_stimPos_0ms, {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, rheel_x_pos_mean_stanceL_stimNeg_0ms, rheel_x_pos_civ_stanceL_stimNeg_0ms, {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'heel_absolute_stanceL_0ms.eps', 'epsc2')
    end
    
    if do_phase_delayed_plots
        % 150ms
        figure; axes; hold on; title('right foot medial-lateral heel pos, 150ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, rheel_x_pos_mean_stanceL_control, rheel_x_pos_civ_stanceL_control, {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, rheel_x_pos_mean_stanceL_stimPos_150ms, rheel_x_pos_civ_stanceL_stimPos_150ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, rheel_x_pos_mean_stanceL_stimNeg_150ms, rheel_x_pos_civ_stanceL_stimNeg_150ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'heel_absolute_stanceL_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title('left foot medial-lateral heel pos, 450ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, lheel_x_pos_mean_stanceR_control, lheel_x_pos_civ_stanceR_control, {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, lheel_x_pos_mean_stanceR_stimPos_450ms, lheel_x_pos_civ_stanceR_stimPos_450ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, lheel_x_pos_mean_stanceR_stimNeg_450ms, lheel_x_pos_civ_stanceR_stimNeg_450ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'heel_absolute_stanceR_450ms.eps', 'epsc2')
        end
    end
end

if do_heel_plots_response_right
    % 0ms
    figure; axes; hold on; title('left foot medial-lateral heel pos response, 0ms'); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, lheel_x_pos_response_mean_stanceR_stimPos_0ms, lheel_x_pos_response_civ_stanceR_stimPos_0ms, {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, lheel_x_pos_response_mean_stanceR_stimNeg_0ms, lheel_x_pos_response_civ_stanceR_stimPos_0ms, {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-0.015 0.015]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'heel_response_stanceR_0ms.eps', 'epsc2')
    end
    
    if do_phase_delayed_plots
        % 150ms
        figure; axes; hold on; title('left foot medial-lateral heel pos response, 150ms'); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, lheel_x_pos_response_mean_stanceR_stimPos_150ms, lheel_x_pos_response_civ_stanceR_stimPos_150ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, lheel_x_pos_response_mean_stanceR_stimNeg_150ms, lheel_x_pos_response_civ_stanceR_stimPos_150ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-0.015 0.015]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'heel_response_stanceR_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title('right foot medial-lateral heel pos response, 450ms'); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, rheel_x_pos_response_mean_stanceL_stimPos_450ms, rheel_x_pos_response_civ_stanceL_stimPos_450ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, rheel_x_pos_response_mean_stanceL_stimNeg_450ms, rheel_x_pos_response_civ_stanceL_stimPos_450ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-0.015 0.015]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'heel_response_stanceL_450ms.eps', 'epsc2')
        end
    end
end

if do_heel_plots_response_left
    % 0ms
    figure; axes; hold on; title('right foot medial-lateral heel pos response, 0ms'); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, rheel_x_pos_response_mean_stanceL_stimPos_0ms, rheel_x_pos_response_civ_stanceL_stimPos_0ms, {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, rheel_x_pos_response_mean_stanceL_stimNeg_0ms, rheel_x_pos_response_civ_stanceL_stimPos_0ms, {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-0.015 0.015]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'heel_response_stanceL_0ms.eps', 'epsc2')
    end
    
    if do_phase_delayed_plots
        % 150ms
        figure; axes; hold on; title('right foot medial-lateral heel pos response, 150ms'); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, rheel_x_pos_response_mean_stanceL_stimPos_150ms, rheel_x_pos_response_civ_stanceL_stimPos_150ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, rheel_x_pos_response_mean_stanceL_stimNeg_150ms, rheel_x_pos_response_civ_stanceL_stimPos_150ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-0.015 0.015]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'heel_response_stanceL_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title('left foot medial-lateral heel pos response, 450ms'); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, lheel_x_pos_response_mean_stanceR_stimPos_450ms, lheel_x_pos_response_civ_stanceR_stimPos_450ms, {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, lheel_x_pos_response_mean_stanceR_stimNeg_450ms, lheel_x_pos_response_civ_stanceR_stimPos_450ms, {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-0.015 0.015]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'heel_response_stanceR_450ms.eps', 'epsc2')
        end
    end
end














return
