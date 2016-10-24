
check_out_origin                    = 0;

% CoP
do_cop_plots_single                 = 0;
do_cop_plots_absolute_right         = 1;
do_cop_plots_absolute_left          = 0;
do_cop_plots_response_right         = 0;
do_cop_plots_response_left          = 0;

% fx
do_fx_plots_single                  = 0;
do_fx_plots_absolute_right          = 0;
do_fx_plots_absolute_left           = 0;
do_fx_plots_response_right          = 0;
do_fx_plots_response_left           = 0;

% fz
do_fz_plots_single                  = 0;
do_fz_plots_absolute_right          = 0;
do_fz_plots_absolute_left           = 0;
do_fz_plots_response_right          = 0;
do_fz_plots_response_left           = 0;

% my
do_my_plots_single                  = 0;
do_my_plots_absolute_right          = 0;
do_my_plots_absolute_left           = 0;
do_my_plots_response_right          = 0;
do_my_plots_response_left           = 0;

% heel
do_heel_plots_single                = 0;
do_heel_plots_absolute_right        = 0;
do_heel_plots_absolute_left         = 0;
do_heel_plots_response_right        = 0;
do_heel_plots_response_left         = 0;

% pelvis_pos
do_pelvis_pos_plots_single          = 0;
do_pelvis_pos_plots_absolute_right  = 0;
do_pelvis_pos_plots_response_right  = 0;

% pelvis_vel
do_pelvis_vel_plots_single          = 0;
do_pelvis_vel_plots_absolute_right  = 0;
do_pelvis_vel_plots_response_right  = 0;

% pelvis_acc
do_pelvis_acc_plots_single          = 0;
do_pelvis_acc_plots_absolute_right  = 0;
do_pelvis_acc_plots_response_right  = 0;

% emg
do_emg_plots_single                 = 0;
do_emg_plots_absolute_right         = 1;
do_emg_plots_absolute_left          = 1;
do_emg_plots_response_right         = 0;
do_emg_plots_response_left          = 0;

% step
do_step_response_plots              = 0;

do_phase_delayed_plots              = 1;

save_figures                        = 1;

load subjectInfo.mat;
load(makeFileName(date, subject_id, 'resultsConditions'));
load(makeFileName(date, subject_id, 'resultsForceplate'));
load(makeFileName(date, subject_id, 'resultsMarker'));
load(makeFileName(date, subject_id, 'resultsEmg'));

color_control = [0.3 0.1 1];
color_positive = [1 0.3 .1] * 0.7;
color_negative = [0.3 1 0.1] * 0.7;

color_control = [1 1 1] * 0.7;
color_positive = [82 79 161] * 1/255;
color_negative = [253 185 19] * 1/255;


%% check out origin
if check_out_origin
    condition = conditions_stanceR_stimPos_0ms;
%     condition = conditions_stanceL_stimNo;
%     condition = conditions_stanceL_stimNo;
    
%     data = rheel_x_pos_normalized_total(:, condition);
    data = rheel_x_pos_normalized_total(:, condition);
    data = pelvis_x_pos_normalized_total(:, condition);
    
    time_data = repmat(time_normalized', 1, size(data, 2));
    index_data = repmat((1 : size(data, 2)), length(time_normalized), 1);
    
    origin_trial_list_condition = origin_trial_list_total(condition);
    origin_start_time_list_condition = origin_start_time_list_total(condition);
    origin_end_time_list_condition = origin_end_time_list_total(condition);
    
    figure; axes; hold on; title(['origin finder - ' subject_id]); set(gca, 'Fontsize', 12)
    plot3(time_data, data, index_data)
return    
 
end

%% cop plots
if do_cop_plots_single
    figure; axes; hold on; title(['left foot medial-lateral CoP, N = ' num2str(sum(conditions_stanceL_stimNo)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimNo), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['right foot medial-lateral CoP, N = ' num2str(sum(conditions_stanceR_stimNo)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimNo), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['left foot medial-lateral CoP, 0ms, N = ' num2str(sum(conditions_stanceL_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    if any(conditions_stanceL_stimNeg_0ms) plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5)); end
    if any(conditions_stanceL_stimPos_0ms) plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimPos_0ms), 'color', lightenColor(color_positive, 0.5)); end
    plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    figure; axes; hold on; title(['right foot medial-lateral CoP, 0ms, N = ' num2str(sum(conditions_stanceR_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
    plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
    plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    if do_phase_delayed_plots
        figure; axes; hold on; title(['right foot medial-lateral CoP, 150ms, N = ' num2str(sum(conditions_stanceR_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['left foot medial-lateral CoP, 450ms, N = ' num2str(sum(conditions_stanceL_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);
        
        figure; axes; hold on; title(['left foot medial-lateral CoP, 150ms, N = ' num2str(sum(conditions_stanceL_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        if any(conditions_stanceL_stimNeg_150ms) plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5)); end
        if any(conditions_stanceL_stimPos_150ms) plot(time_normalized, lcop_x_normalized_total(:, conditions_stanceL_stimPos_150ms), 'color', lightenColor(color_positive, 0.5)); end
        plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['right foot medial-lateral CoP, 450ms, N = ' num2str(sum(conditions_stanceR_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        if any(conditions_stanceR_stimNeg_450ms) plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5)); end
        if any(conditions_stanceR_stimPos_450ms) plot(time_normalized, rcop_x_normalized_total(:, conditions_stanceR_stimPos_450ms), 'color', lightenColor(color_positive, 0.5)); end
        plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);
    end
end

if do_cop_plots_absolute_right
    % 0ms
    figure; axes; hold on; title(['right foot medial-lateral CoP, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(rcop_x_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), cinv(rcop_x_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), cinv(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['right foot medial-lateral CoP, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(rcop_x_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), cinv(rcop_x_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), cinv(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['left foot medial-lateral CoP, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(lcop_x_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), cinv(lcop_x_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), cinv(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
    figure; axes; hold on; title(['left foot medial-lateral CoP, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(lcop_x_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), cinv(lcop_x_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), cinv(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['left foot medial-lateral CoP, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(lcop_x_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), cinv(lcop_x_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), cinv(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['right foot medial-lateral CoP, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(rcop_x_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), cinv(rcop_x_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), cinv(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
    figure; axes; hold on; %title(['right foot medial-lateral CoP response, 0ms - ' subject_id]); 
    set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(rcop_x_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), cinv(rcop_x_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(rcop_x_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), cinv(rcop_x_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-0.015 0.015]);
%     this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
%     set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'cop_response_stanceR_0ms.eps', 'epsc2')
    end
    
    if do_phase_delayed_plots
        % 150ms
        figure; axes; hold on; %title(['right foot medial-lateral CoP response, 150ms - ' subject_id]); 
        set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(rcop_x_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), cinv(rcop_x_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rcop_x_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), cinv(rcop_x_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-0.015 0.015]);
%         this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
%         set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_response_stanceR_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; %title(['left foot medial-lateral CoP response, 450ms - ' subject_id]); 
        set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), cinv(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), cinv(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-0.015 0.015]);
%         this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
%         set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_response_stanceL_450ms.eps', 'epsc2')
        end
    end
end

if do_cop_plots_response_left
    figure; axes; hold on; title(['left foot medial-lateral CoP response, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2), cinv(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2), cinv(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['left foot medial-lateral CoP response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2), cinv(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2), cinv(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['right foot medial-lateral CoP response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(rcop_x_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2), cinv(rcop_x_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rcop_x_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2), cinv(rcop_x_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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

%% fx plots
if do_fx_plots_single
    figure; axes; hold on; title(['left foot f_x, N = ' num2str(sum(conditions_stanceL_stimNo))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, fxl_normalized_total(:, conditions_stanceL_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimNo), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['right foot f_x, N = ' num2str(sum(conditions_stanceR_stimNo))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, fxr_normalized_total(:, conditions_stanceR_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimNo), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['left foot f_x, 0ms, N = ' num2str(sum(conditions_stanceL_0ms))]); set(gca, 'Fontsize', 12)
    if any(conditions_stanceL_stimNeg_0ms) plot(time_normalized, fxl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5)); end
    if any(conditions_stanceL_stimPos_0ms) plot(time_normalized, fxl_normalized_total(:, conditions_stanceL_stimPos_0ms), 'color', lightenColor(color_positive, 0.5)); end
    plot(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    figure; axes; hold on; title(['right foot f_x, 0ms, N = ' num2str(sum(conditions_stanceR_0ms))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, fxr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
    plot(time_normalized, fxr_normalized_total(:, conditions_stanceR_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
    plot(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    if do_phase_delayed_plots
        figure; axes; hold on; title(['right foot f_x, 150ms, N = ' num2str(sum(conditions_stanceR_150ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, fxr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, fxr_normalized_total(:, conditions_stanceR_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['left foot f_x, 450ms, N = ' num2str(sum(conditions_stanceL_450ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, fxl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, fxl_normalized_total(:, conditions_stanceL_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);
        
        figure; axes; hold on; title(['left foot f_x, 150ms, N = ' num2str(sum(conditions_stanceL_150ms))]); set(gca, 'Fontsize', 12)
        if any(conditions_stanceL_stimNeg_150ms) plot(time_normalized, fxl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5)); end
        if any(conditions_stanceL_stimPos_150ms) plot(time_normalized, fxl_normalized_total(:, conditions_stanceL_stimPos_150ms), 'color', lightenColor(color_positive, 0.5)); end
        plot(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['right foot f_x, 450ms, N = ' num2str(sum(conditions_stanceR_450ms))]); set(gca, 'Fontsize', 12)
        if any(conditions_stanceR_stimNeg_450ms) plot(time_normalized, fxr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5)); end
        if any(conditions_stanceR_stimPos_450ms) plot(time_normalized, fxr_normalized_total(:, conditions_stanceR_stimPos_450ms), 'color', lightenColor(color_positive, 0.5)); end
        plot(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);
    end
end

if do_fx_plots_absolute_right
    % 0ms
    figure; axes; hold on; title(['right foot f_x, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(fxr_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), cinv(fxr_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), cinv(fxr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['right foot f_x, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(fxr_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), cinv(fxr_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), cinv(fxr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['left foot f_x, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(fxl_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), cinv(fxl_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), cinv(fxl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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

if do_fx_plots_absolute_left
    figure; axes; hold on; title(['left foot f_x, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(fxl_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), cinv(fxl_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), cinv(fxl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['left foot f_x, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(fxl_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), cinv(fxl_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fxl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), cinv(fxl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['right foot f_x, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(fxr_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), cinv(fxr_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fxr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), cinv(fxr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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

if do_fx_plots_response_right
    % 0ms
    figure; axes; hold on; title(['right foot f_x response, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(fxr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), cinv(fxr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(fxr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), cinv(fxr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-15 15]);
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
        figure; axes; hold on; title(['right foot f_x response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(fxr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), cinv(fxr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fxr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), cinv(fxr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-15 15]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_response_stanceR_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title(['left foot f_x response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(fxl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), cinv(fxl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fxl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), cinv(fxl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-15 15]);
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

if do_fx_plots_response_left
    figure; axes; hold on; title(['left foot f_x response, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(fxl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2), cinv(fxl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(fxl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2), cinv(fxl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-15 15]);
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
        figure; axes; hold on; title(['left foot f_x response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(fxl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2), cinv(fxl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fxl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2), cinv(fxl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-15 15]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_response_stanceL_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title(['right foot f_x response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(fxr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2), cinv(fxr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fxr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2), cinv(fxr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-15 15]);
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

%% fz plots
if do_fz_plots_single
    figure; axes; hold on; title(['left foot f_z, N = ' num2str(sum(conditions_stanceL_stimNo))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, fzl_normalized_total(:, conditions_stanceL_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimNo), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['right foot f_z, N = ' num2str(sum(conditions_stanceR_stimNo))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, fzr_normalized_total(:, conditions_stanceR_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimNo), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['left foot f_z, 0ms, N = ' num2str(sum(conditions_stanceL_0ms))]); set(gca, 'Fontsize', 12)
    if any(conditions_stanceL_stimNeg_0ms) plot(time_normalized, fzl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5)); end
    if any(conditions_stanceL_stimPos_0ms) plot(time_normalized, fzl_normalized_total(:, conditions_stanceL_stimPos_0ms), 'color', lightenColor(color_positive, 0.5)); end
    plot(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    figure; axes; hold on; title(['right foot f_z, 0ms, N = ' num2str(sum(conditions_stanceR_0ms))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, fzr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
    plot(time_normalized, fzr_normalized_total(:, conditions_stanceR_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
    plot(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    if do_phase_delayed_plots
        figure; axes; hold on; title(['right foot f_z, 150ms, N = ' num2str(sum(conditions_stanceR_150ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, fzr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, fzr_normalized_total(:, conditions_stanceR_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['left foot f_z, 450ms, N = ' num2str(sum(conditions_stanceL_450ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, fzl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, fzl_normalized_total(:, conditions_stanceL_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);
        
        figure; axes; hold on; title(['left foot f_z, 150ms, N = ' num2str(sum(conditions_stanceL_150ms))]); set(gca, 'Fontsize', 12)
        if any(conditions_stanceL_stimNeg_150ms) plot(time_normalized, fzl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5)); end
        if any(conditions_stanceL_stimPos_150ms) plot(time_normalized, fzl_normalized_total(:, conditions_stanceL_stimPos_150ms), 'color', lightenColor(color_positive, 0.5)); end
        plot(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['right foot f_z, 450ms, N = ' num2str(sum(conditions_stanceR_450ms))]); set(gca, 'Fontsize', 12)
        if any(conditions_stanceR_stimNeg_450ms) plot(time_normalized, fzr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5)); end
        if any(conditions_stanceR_stimPos_450ms) plot(time_normalized, fzr_normalized_total(:, conditions_stanceR_stimPos_450ms), 'color', lightenColor(color_positive, 0.5)); end
        plot(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);
    end
end

if do_fz_plots_absolute_right
    % 0ms
    figure; axes; hold on; title(['right foot f_z, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(fzr_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), cinv(fzr_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), cinv(fzr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['right foot f_z, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(fzr_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), cinv(fzr_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), cinv(fzr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['left foot f_z, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(fzl_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), cinv(fzl_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), cinv(fzl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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

if do_fz_plots_absolute_left
    figure; axes; hold on; title(['left foot f_z, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(fzl_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), cinv(fzl_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), cinv(fzl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['left foot f_z, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(fzl_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), cinv(fzl_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fzl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), cinv(fzl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['right foot f_z, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(fzr_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), cinv(fzr_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fzr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), cinv(fzr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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

if do_fz_plots_response_right
    % 0ms
    figure; axes; hold on; title(['right foot f_z response, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(fzr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), cinv(fzr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(fzr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), cinv(fzr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-30 30]);
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
        figure; axes; hold on; title(['right foot f_z response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(fzr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), cinv(fzr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fzr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), cinv(fzr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-30 30]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_response_stanceR_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title(['left foot f_z response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(fzl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), cinv(fzl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fzl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), cinv(fzl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-30 30]);
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

if do_fz_plots_response_left
    figure; axes; hold on; title(['left foot f_z response, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(fzl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2), cinv(fzl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(fzl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2), cinv(fzl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-30 30]);
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
        figure; axes; hold on; title(['left foot f_z response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(fzl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2), cinv(fzl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fzl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2), cinv(fzl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-30 30]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_response_stanceL_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title(['right foot f_z response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(fzr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2), cinv(fzr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(fzr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2), cinv(fzr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-30 30]);
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

%% my plots
if do_my_plots_single
    figure; axes; hold on; title(['left foot m_y, N = ' num2str(sum(conditions_stanceL_stimNo))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, myl_normalized_total(:, conditions_stanceL_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimNo), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['right foot m_y, N = ' num2str(sum(conditions_stanceR_stimNo))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, myr_normalized_total(:, conditions_stanceR_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimNo), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['left foot m_y, 0ms, N = ' num2str(sum(conditions_stanceL_0ms))]); set(gca, 'Fontsize', 12)
    if any(conditions_stanceL_stimNeg_0ms) plot(time_normalized, myl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5)); end
    if any(conditions_stanceL_stimPos_0ms) plot(time_normalized, myl_normalized_total(:, conditions_stanceL_stimPos_0ms), 'color', lightenColor(color_positive, 0.5)); end
    plot(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    figure; axes; hold on; title(['right foot m_y, 0ms, N = ' num2str(sum(conditions_stanceR_0ms))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, myr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
    plot(time_normalized, myr_normalized_total(:, conditions_stanceR_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
    plot(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    if do_phase_delayed_plots
        figure; axes; hold on; title(['right foot m_y, 150ms, N = ' num2str(sum(conditions_stanceR_150ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, myr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, myr_normalized_total(:, conditions_stanceR_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['left foot m_y, 450ms, N = ' num2str(sum(conditions_stanceL_450ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, myl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, myl_normalized_total(:, conditions_stanceL_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);
        
        figure; axes; hold on; title(['left foot m_y, 150ms, N = ' num2str(sum(conditions_stanceL_150ms))]); set(gca, 'Fontsize', 12)
        if any(conditions_stanceL_stimNeg_150ms) plot(time_normalized, myl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5)); end
        if any(conditions_stanceL_stimPos_150ms) plot(time_normalized, myl_normalized_total(:, conditions_stanceL_stimPos_150ms), 'color', lightenColor(color_positive, 0.5)); end
        plot(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['right foot m_y, 450ms, N = ' num2str(sum(conditions_stanceR_450ms))]); set(gca, 'Fontsize', 12)
        if any(conditions_stanceR_stimNeg_450ms) plot(time_normalized, myr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5)); end
        if any(conditions_stanceR_stimPos_450ms) plot(time_normalized, myr_normalized_total(:, conditions_stanceR_stimPos_450ms), 'color', lightenColor(color_positive, 0.5)); end
        plot(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);
    end
end

if do_my_plots_absolute_right
    % 0ms
    figure; axes; hold on; title(['right foot m_y, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(myr_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), cinv(myr_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), cinv(myr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['right foot m_y, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(myr_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), cinv(myr_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), cinv(myr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['left foot m_y, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(myl_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), cinv(myl_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), cinv(myl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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

if do_my_plots_absolute_left
    figure; axes; hold on; title(['left foot m_y, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(myl_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), cinv(myl_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), cinv(myl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['left foot m_y, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(myl_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), cinv(myl_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(myl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), cinv(myl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['right foot m_y, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(myr_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), cinv(myr_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(myr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), cinv(myr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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

if do_my_plots_response_right
    % 0ms
    figure; axes; hold on; title(['right foot m_y response, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(myr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), cinv(myr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(myr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), cinv(myr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-20 20]);
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
        figure; axes; hold on; title(['right foot m_y response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(myr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), cinv(myr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(myr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), cinv(myr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-20 20]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_response_stanceR_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title(['left foot m_y response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(myl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), cinv(myl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(myl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), cinv(myl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-20 20]);
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

if do_my_plots_response_left
    figure; axes; hold on; title(['left foot m_y response, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(myl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2), cinv(myl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(myl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2), cinv(myl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-20 20]);
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
        figure; axes; hold on; title(['left foot m_y response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(myl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2), cinv(myl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(myl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2), cinv(myl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-20 20]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'cop_response_stanceL_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title(['right foot m_y response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(myr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2), cinv(myr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(myr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2), cinv(myr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)], 'ylim', [-20 20]);
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
    if any(conditions_stanceL_stimNeg_0ms) plot(time_normalized, rheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5)); end
    if any(conditions_stanceL_stimPos_0ms) plot(time_normalized, rheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_0ms), 'color', lightenColor(color_positive, 0.5)); end
    plot(time_normalized, mean(rheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(rheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    if do_phase_delayed_plots
        figure; axes; hold on; title(['right foot medial-lateral heel pos, 150ms, N = ' num2str(sum(conditions_stanceL_150ms))]); set(gca, 'Fontsize', 12)
        if any(conditions_stanceL_stimNeg_150ms) plot(time_normalized, lheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5)); end
        if any(conditions_stanceL_stimPos_150ms) plot(time_normalized, lheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_150ms), 'color', lightenColor(color_positive, 0.5)); end
        plot(time_normalized, mean(lheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(lheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['left foot medial-lateral heel pos, 450ms, N = ' num2str(sum(conditions_stanceR_450ms))]); set(gca, 'Fontsize', 12)
        if any(conditions_stanceR_stimNeg_450ms) plot(time_normalized, rheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5)); end
        if any(conditions_stanceR_stimPos_450ms) plot(time_normalized, rheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_450ms), 'color', lightenColor(color_positive, 0.5)); end
        plot(time_normalized, mean(rheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(rheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);    
    end    
end

if do_heel_plots_absolute_right
    % 0ms
    figure; axes; hold on; title(['left foot medial-lateral heel pos, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(lheel_x_pos_normalized_total(:, conditions_stanceR_0ms), 2), cinv(lheel_x_pos_normalized_total(:, conditions_stanceR_0ms), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(lheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), cinv(lheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(lheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), cinv(lheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['left foot medial-lateral heel pos, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(lheel_x_pos_normalized_total(:, conditions_stanceR_150ms), 2), cinv(lheel_x_pos_normalized_total(:, conditions_stanceR_150ms), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(lheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), cinv(lheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(lheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), cinv(lheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['right foot medial-lateral heel pos, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(rheel_x_pos_normalized_total(:, conditions_stanceL_450ms), 2), cinv(rheel_x_pos_normalized_total(:, conditions_stanceL_450ms), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(rheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), cinv(rheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), cinv(rheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
    figure; axes; hold on; title(['right foot medial-lateral heel pos, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
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
        figure; axes; hold on; title(['right foot medial-lateral heel pos, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
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
        figure; axes; hold on; title(['left foot medial-lateral heel pos, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
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
    figure; axes; hold on; title(['left foot medial-lateral heel pos response, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), cinv(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), cinv(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['left foot medial-lateral heel pos response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), cinv(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), cinv(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['right foot medial-lateral heel pos response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), cinv(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), cinv(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
    figure; axes; hold on; title(['right foot medial-lateral heel pos response, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
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
        figure; axes; hold on; title(['right foot medial-lateral heel pos response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
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
        figure; axes; hold on; title(['left foot medial-lateral heel pos response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
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

%% pelvis pos plots
if do_pelvis_pos_plots_single
    figure; axes; hold on; title(['medial-lateral pelvis pos, N = ' num2str(sum(conditions_stanceL_stimNo))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, pelvis_x_pos_normalized_total(:, conditions_stanceL_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceL_stimNo), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['medial-lateral pelvis pos, 0ms, N = ' num2str(sum(conditions_stanceR_0ms))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, pelvis_x_pos_normalized_total(:, conditions_stanceR_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
    plot(time_normalized, pelvis_x_pos_normalized_total(:, conditions_stanceR_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
    plot(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    if do_phase_delayed_plots
        figure; axes; hold on; title(['medial-lateral pelvis pos, 150ms, N = ' num2str(sum(conditions_stanceR_150ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, pelvis_x_pos_normalized_total(:, conditions_stanceR_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, pelvis_x_pos_normalized_total(:, conditions_stanceR_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['left foot medial-lateral heel pos, 450ms, N = ' num2str(sum(conditions_stanceL_450ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, pelvis_x_pos_normalized_total(:, conditions_stanceL_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, pelvis_x_pos_normalized_total(:, conditions_stanceL_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);    
    end    
end

if do_pelvis_pos_plots_absolute_right
    % 0ms
    figure; axes; hold on; title(['medial-lateral pelvis pos, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceR_0ms), 2), cinv(pelvis_x_pos_normalized_total(:, conditions_stanceR_0ms), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), cinv(pelvis_x_pos_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), cinv(pelvis_x_pos_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['medial-lateral pelvis pos, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceR_150ms), 2), cinv(pelvis_x_pos_normalized_total(:, conditions_stanceR_150ms), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), cinv(pelvis_x_pos_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), cinv(pelvis_x_pos_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['medial-lateral pelvis pos, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceL_450ms), 2), cinv(pelvis_x_pos_normalized_total(:, conditions_stanceL_450ms), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), cinv(pelvis_x_pos_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_pos_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), cinv(pelvis_x_pos_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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

if do_pelvis_pos_plots_response_right
    % 0ms
    figure; axes; hold on; title(['left foot medial-lateral heel pos response, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_pos_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), cinv(pelvis_x_pos_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), cinv(pelvis_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['left foot medial-lateral heel pos response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_pos_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), cinv(pelvis_x_pos_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), cinv(pelvis_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['right foot medial-lateral heel pos response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_pos_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), cinv(pelvis_x_pos_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), cinv(pelvis_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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

%% pelvis vel plots
if do_pelvis_vel_plots_single
    figure; axes; hold on; title(['medial-lateral pelvis vel, N = ' num2str(sum(conditions_stanceL_stimNo))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, pelvis_x_vel_normalized_total(:, conditions_stanceL_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceL_stimNo), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['medial-lateral pelvis vel, 0ms, N = ' num2str(sum(conditions_stanceR_0ms))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, pelvis_x_vel_normalized_total(:, conditions_stanceR_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
    plot(time_normalized, pelvis_x_vel_normalized_total(:, conditions_stanceR_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
    plot(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    if do_phase_delayed_plots
        figure; axes; hold on; title(['medial-lateral pelvis vel, 150ms, N = ' num2str(sum(conditions_stanceR_150ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, pelvis_x_vel_normalized_total(:, conditions_stanceR_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, pelvis_x_vel_normalized_total(:, conditions_stanceR_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['left foot medial-lateral heel vel, 450ms, N = ' num2str(sum(conditions_stanceL_450ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, pelvis_x_vel_normalized_total(:, conditions_stanceL_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, pelvis_x_vel_normalized_total(:, conditions_stanceL_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);    
    end    
end

if do_pelvis_vel_plots_absolute_right
    % 0ms
    figure; axes; hold on; title(['medial-lateral pelvis vel, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceR_0ms), 2), cinv(pelvis_x_vel_normalized_total(:, conditions_stanceR_0ms), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), cinv(pelvis_x_vel_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), cinv(pelvis_x_vel_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['medial-lateral pelvis vel, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceR_150ms), 2), cinv(pelvis_x_vel_normalized_total(:, conditions_stanceR_150ms), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), cinv(pelvis_x_vel_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), cinv(pelvis_x_vel_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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
        figure; axes; hold on; title(['medial-lateral pelvis vel, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceL_450ms), 2), cinv(pelvis_x_vel_normalized_total(:, conditions_stanceL_450ms), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), cinv(pelvis_x_vel_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_vel_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), cinv(pelvis_x_vel_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
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

if do_pelvis_vel_plots_response_right
    % 0ms
    figure; axes; hold on; title(['medial-lateral pelvis vel response, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_vel_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), cinv(pelvis_x_vel_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_vel_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), cinv(pelvis_x_vel_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
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
        figure; axes; hold on; title(['medial-lateral pelvis vel response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_vel_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), cinv(pelvis_x_vel_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_vel_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), cinv(pelvis_x_vel_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'heel_response_stanceR_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title(['medial-lateral pelvis vel response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_vel_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), cinv(pelvis_x_vel_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_vel_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), cinv(pelvis_x_vel_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
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

%% pelvis acc plots
if do_pelvis_acc_plots_single
    figure; axes; hold on; title(['medial-lateral pelvis acc, N = ' num2str(sum(conditions_stanceL_stimNo))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, pelvis_x_acc_normalized_total(:, conditions_stanceL_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceL_stimNo), 2), 'color', color_control, 'linewidth', 5);
    
    figure; axes; hold on; title(['medial-lateral pelvis acc, 0ms, N = ' num2str(sum(conditions_stanceR_0ms))]); set(gca, 'Fontsize', 12)
    plot(time_normalized, pelvis_x_acc_normalized_total(:, conditions_stanceR_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
    plot(time_normalized, pelvis_x_acc_normalized_total(:, conditions_stanceR_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
    plot(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    plot(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    
    if do_phase_delayed_plots
        figure; axes; hold on; title(['medial-lateral pelvis acc, 150ms, N = ' num2str(sum(conditions_stanceR_150ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, pelvis_x_acc_normalized_total(:, conditions_stanceR_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, pelvis_x_acc_normalized_total(:, conditions_stanceR_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['medial-lateral pelvis acc, 450ms, N = ' num2str(sum(conditions_stanceL_450ms))]); set(gca, 'Fontsize', 12)
        plot(time_normalized, pelvis_x_acc_normalized_total(:, conditions_stanceL_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, pelvis_x_acc_normalized_total(:, conditions_stanceL_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        plot(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);    
    end    
end

if do_pelvis_acc_plots_absolute_right
    % 0ms
    figure; axes; hold on; title(['medial-lateral pelvis acc, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceR_0ms), 2), cinv(pelvis_x_acc_normalized_total(:, conditions_stanceR_0ms), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), cinv(pelvis_x_acc_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), cinv(pelvis_x_acc_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'pelvis_absolute_stanceR_0ms.eps', 'epsc2')
    end
    
    if do_phase_delayed_plots
        % 150ms
        figure; axes; hold on; title(['medial-lateral pelvis acc, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceR_150ms), 2), cinv(pelvis_x_acc_normalized_total(:, conditions_stanceR_150ms), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), cinv(pelvis_x_acc_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), cinv(pelvis_x_acc_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'pelvis_absolute_stanceR_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title(['medial-lateral pelvis acc, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceL_450ms), 2), cinv(pelvis_x_acc_normalized_total(:, conditions_stanceL_450ms), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), cinv(pelvis_x_acc_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_acc_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), cinv(pelvis_x_acc_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'pelvis_absolute_stanceL_450ms.eps', 'epsc2')
        end
    end
end

if do_pelvis_acc_plots_response_right
    % 0ms
    figure; axes; hold on; title(['medial-lateral pelvis acc response, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_acc_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), cinv(pelvis_x_acc_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_acc_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), cinv(pelvis_x_acc_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'pelvis_response_stanceR_0ms.eps', 'epsc2')
    end
    
    if do_phase_delayed_plots
        % 150ms
        figure; axes; hold on; title(['medial-lateral pelvis acc response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_acc_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), cinv(pelvis_x_acc_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_acc_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), cinv(pelvis_x_acc_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'pelvis_response_stanceR_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title(['medial-lateral pelvis acc response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(pelvis_x_acc_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), cinv(pelvis_x_acc_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(pelvis_x_acc_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), cinv(pelvis_x_acc_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'pelvis_response_stanceL_450ms.eps', 'epsc2')
        end
    end
end



%% emg plots
if do_emg_plots_single
    % no stim
    figure; axes; hold on; title(['left gluteus medius, left stance, control, N = ' num2str(sum(conditions_stanceL_stimNo)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    plot(time_normalized, lglutmed_normalized_total(:, conditions_stanceL_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceL_stimNo), 2), 'color', color_control, 'linewidth', 5);

    figure; axes; hold on; title(['left tibialis anterior, left stance, control, N = ' num2str(sum(conditions_stanceL_stimNo)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    plot(time_normalized, ltibiant_normalized_total(:, conditions_stanceL_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceL_stimNo), 2), 'color', color_control, 'linewidth', 5);

    figure; axes; hold on; title(['left peroneus longus, left stance, control, N = ' num2str(sum(conditions_stanceL_stimNo)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    plot(time_normalized, lperolng_normalized_total(:, conditions_stanceL_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceL_stimNo), 2), 'color', color_control, 'linewidth', 5);

    figure; axes; hold on; title(['right gluteus medius, left stance, control, N = ' num2str(sum(conditions_stanceL_stimNo)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    plot(time_normalized, rglutmed_normalized_total(:, conditions_stanceL_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceL_stimNo), 2), 'color', color_control, 'linewidth', 5);

    figure; axes; hold on; title(['right tibialis anterior, left stance, control, N = ' num2str(sum(conditions_stanceL_stimNo)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    plot(time_normalized, rtibiant_normalized_total(:, conditions_stanceL_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceL_stimNo), 2), 'color', color_control, 'linewidth', 5);

    figure; axes; hold on; title(['right peroneus longus, left stance, control, N = ' num2str(sum(conditions_stanceL_stimNo)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    plot(time_normalized, rperolng_normalized_total(:, conditions_stanceL_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceL_stimNo), 2), 'color', color_control, 'linewidth', 5);

    figure; axes; hold on; title(['left gluteus medius, right stance, control, N = ' num2str(sum(conditions_stanceR_stimNo)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    plot(time_normalized, lglutmed_normalized_total(:, conditions_stanceR_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceR_stimNo), 2), 'color', color_control, 'linewidth', 5);

    figure; axes; hold on; title(['left tibialis anterior, right stance, control, N = ' num2str(sum(conditions_stanceR_stimNo)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    plot(time_normalized, ltibiant_normalized_total(:, conditions_stanceR_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceR_stimNo), 2), 'color', color_control, 'linewidth', 5);

    figure; axes; hold on; title(['left peroneus longus, right stance, control, N = ' num2str(sum(conditions_stanceR_stimNo)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    plot(time_normalized, lperolng_normalized_total(:, conditions_stanceR_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceR_stimNo), 2), 'color', color_control, 'linewidth', 5);

    figure; axes; hold on; title(['right gluteus medius, right stance, control, N = ' num2str(sum(conditions_stanceR_stimNo)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    plot(time_normalized, rglutmed_normalized_total(:, conditions_stanceR_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceR_stimNo), 2), 'color', color_control, 'linewidth', 5);

    figure; axes; hold on; title(['right tibialis anterior, right stance, control, N = ' num2str(sum(conditions_stanceR_stimNo)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    plot(time_normalized, rtibiant_normalized_total(:, conditions_stanceR_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceR_stimNo), 2), 'color', color_control, 'linewidth', 5);

    figure; axes; hold on; title(['right peroneus longus, right stance, control, N = ' num2str(sum(conditions_stanceR_stimNo)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
    plot(time_normalized, rperolng_normalized_total(:, conditions_stanceR_stimNo), 'color', lightenColor(color_control, 0.5));
    plot(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceR_stimNo), 2), 'color', color_control, 'linewidth', 5);

    % pos
    if any(conditions_stanceL_stimPos_0ms)
        figure; axes; hold on; title(['left gluteus medius, left stance, 0ms, N = ' num2str(sum(conditions_stanceL_stimPos_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, lglutmed_normalized_total(:, conditions_stanceL_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
        
        figure; axes; hold on; title(['left tibialis anterior, left stance, 0ms, N = ' num2str(sum(conditions_stanceL_stimPos_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, ltibiant_normalized_total(:, conditions_stanceL_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['left peroneus longus, left stance, 0ms, N = ' num2str(sum(conditions_stanceL_stimPos_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, lperolng_normalized_total(:, conditions_stanceL_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['right gluteus medius, left stance, 0ms, N = ' num2str(sum(conditions_stanceL_stimPos_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rglutmed_normalized_total(:, conditions_stanceL_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['right tibialis anterior, left stance, 0ms, N = ' num2str(sum(conditions_stanceL_stimPos_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rtibiant_normalized_total(:, conditions_stanceL_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['right peroneus longus, left stance, 0ms, N = ' num2str(sum(conditions_stanceL_stimPos_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rperolng_normalized_total(:, conditions_stanceL_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceL_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    end
    
    if any(conditions_stanceR_stimPos_0ms)
        figure; axes; hold on; title(['left gluteus medius, right stance, 0ms, N = ' num2str(sum(conditions_stanceR_stimPos_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, lglutmed_normalized_total(:, conditions_stanceR_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['left tibialis anterior, right stance, 0ms, N = ' num2str(sum(conditions_stanceR_stimPos_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, ltibiant_normalized_total(:, conditions_stanceR_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['left peroneus longus, right stance, 0ms, N = ' num2str(sum(conditions_stanceR_stimPos_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, lperolng_normalized_total(:, conditions_stanceR_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['right gluteus medius, right stance, 0ms, N = ' num2str(sum(conditions_stanceR_stimPos_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rglutmed_normalized_total(:, conditions_stanceR_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['right tibialis anterior, right stance, 0ms, N = ' num2str(sum(conditions_stanceR_stimPos_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rtibiant_normalized_total(:, conditions_stanceR_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);

        figure; axes; hold on; title(['right peroneus longus, right stance, 0ms, N = ' num2str(sum(conditions_stanceR_stimPos_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rperolng_normalized_total(:, conditions_stanceR_stimPos_0ms), 'color', lightenColor(color_positive, 0.5));
        plot(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_positive, 'linewidth', 5);
    end
    
    % neg
    if any(conditions_stanceL_stimPos_0ms)
        figure; axes; hold on; title(['left gluteus medius, left stance, 0ms, N = ' num2str(sum(conditions_stanceL_stimNeg_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, lglutmed_normalized_total(:, conditions_stanceL_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
        
        figure; axes; hold on; title(['left tibialis anterior, left stance, 0ms, N = ' num2str(sum(conditions_stanceL_stimNeg_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, ltibiant_normalized_total(:, conditions_stanceL_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);

        figure; axes; hold on; title(['left peroneus longus, left stance, 0ms, N = ' num2str(sum(conditions_stanceL_stimNeg_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, lperolng_normalized_total(:, conditions_stanceL_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);

        figure; axes; hold on; title(['right gluteus medius, left stance, 0ms, N = ' num2str(sum(conditions_stanceL_stimNeg_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rglutmed_normalized_total(:, conditions_stanceL_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);

        figure; axes; hold on; title(['right tibialis anterior, left stance, 0ms, N = ' num2str(sum(conditions_stanceL_stimNeg_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rtibiant_normalized_total(:, conditions_stanceL_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);

        figure; axes; hold on; title(['right peroneus longus, left stance, 0ms, N = ' num2str(sum(conditions_stanceL_stimNeg_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rperolng_normalized_total(:, conditions_stanceL_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    end
    
    if any(conditions_stanceR_stimNeg_0ms)
        figure; axes; hold on; title(['left gluteus medius, right stance, 0ms, N = ' num2str(sum(conditions_stanceR_stimNeg_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, lglutmed_normalized_total(:, conditions_stanceR_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);

        figure; axes; hold on; title(['left tibialis anterior, right stance, 0ms, N = ' num2str(sum(conditions_stanceR_stimNeg_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, ltibiant_normalized_total(:, conditions_stanceR_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);

        figure; axes; hold on; title(['left peroneus longus, right stance, 0ms, N = ' num2str(sum(conditions_stanceR_stimNeg_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, lperolng_normalized_total(:, conditions_stanceR_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);

        figure; axes; hold on; title(['right gluteus medius, right stance, 0ms, N = ' num2str(sum(conditions_stanceR_stimNeg_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rglutmed_normalized_total(:, conditions_stanceR_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);

        figure; axes; hold on; title(['right tibialis anterior, right stance, 0ms, N = ' num2str(sum(conditions_stanceR_stimNeg_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rtibiant_normalized_total(:, conditions_stanceR_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);

        figure; axes; hold on; title(['right peroneus longus, right stance, 0ms, N = ' num2str(sum(conditions_stanceR_stimNeg_0ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
        plot(time_normalized, rperolng_normalized_total(:, conditions_stanceR_stimNeg_0ms), 'color', lightenColor(color_negative, 0.5));
        plot(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_negative, 'linewidth', 5);
    end    
    
    if do_phase_delayed_plots
        % 150ms
        % pos
        if any(conditions_stanceL_stimPos_150ms)
            figure; axes; hold on; title(['left gluteus medius, left stance, 150ms, N = ' num2str(sum(conditions_stanceL_stimPos_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lglutmed_normalized_total(:, conditions_stanceL_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['left tibialis anterior, left stance, 150ms, N = ' num2str(sum(conditions_stanceL_stimPos_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, ltibiant_normalized_total(:, conditions_stanceL_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['left peroneus longus, left stance, 150ms, N = ' num2str(sum(conditions_stanceL_stimPos_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lperolng_normalized_total(:, conditions_stanceL_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['right gluteus medius, left stance, 150ms, N = ' num2str(sum(conditions_stanceL_stimPos_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rglutmed_normalized_total(:, conditions_stanceL_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['right tibialis anterior, left stance, 150ms, N = ' num2str(sum(conditions_stanceL_stimPos_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rtibiant_normalized_total(:, conditions_stanceL_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['right peroneus longus, left stance, 150ms, N = ' num2str(sum(conditions_stanceL_stimPos_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rperolng_normalized_total(:, conditions_stanceL_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceL_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);
        end

        if any(conditions_stanceR_stimPos_150ms)
            figure; axes; hold on; title(['left gluteus medius, right stance, 150ms, N = ' num2str(sum(conditions_stanceR_stimPos_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lglutmed_normalized_total(:, conditions_stanceR_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['left tibialis anterior, right stance, 150ms, N = ' num2str(sum(conditions_stanceR_stimPos_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, ltibiant_normalized_total(:, conditions_stanceR_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['left peroneus longus, right stance, 150ms, N = ' num2str(sum(conditions_stanceR_stimPos_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lperolng_normalized_total(:, conditions_stanceR_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['right gluteus medius, right stance, 150ms, N = ' num2str(sum(conditions_stanceR_stimPos_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rglutmed_normalized_total(:, conditions_stanceR_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['right tibialis anterior, right stance, 150ms, N = ' num2str(sum(conditions_stanceR_stimPos_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rtibiant_normalized_total(:, conditions_stanceR_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['right peroneus longus, right stance, 150ms, N = ' num2str(sum(conditions_stanceR_stimPos_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rperolng_normalized_total(:, conditions_stanceR_stimPos_150ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), 'color', color_positive, 'linewidth', 5);
        end

        % neg
        if any(conditions_stanceL_stimPos_150ms)
            figure; axes; hold on; title(['left gluteus medius, left stance, 150ms, N = ' num2str(sum(conditions_stanceL_stimNeg_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lglutmed_normalized_total(:, conditions_stanceL_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['left tibialis anterior, left stance, 150ms, N = ' num2str(sum(conditions_stanceL_stimNeg_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, ltibiant_normalized_total(:, conditions_stanceL_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['left peroneus longus, left stance, 150ms, N = ' num2str(sum(conditions_stanceL_stimNeg_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lperolng_normalized_total(:, conditions_stanceL_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['right gluteus medius, left stance, 150ms, N = ' num2str(sum(conditions_stanceL_stimNeg_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rglutmed_normalized_total(:, conditions_stanceL_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['right tibialis anterior, left stance, 150ms, N = ' num2str(sum(conditions_stanceL_stimNeg_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rtibiant_normalized_total(:, conditions_stanceL_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['right peroneus longus, left stance, 150ms, N = ' num2str(sum(conditions_stanceL_stimNeg_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rperolng_normalized_total(:, conditions_stanceL_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        end

        if any(conditions_stanceR_stimNeg_150ms)
            figure; axes; hold on; title(['left gluteus medius, right stance, 150ms, N = ' num2str(sum(conditions_stanceR_stimNeg_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lglutmed_normalized_total(:, conditions_stanceR_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['left tibialis anterior, right stance, 150ms, N = ' num2str(sum(conditions_stanceR_stimNeg_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, ltibiant_normalized_total(:, conditions_stanceR_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['left peroneus longus, right stance, 150ms, N = ' num2str(sum(conditions_stanceR_stimNeg_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lperolng_normalized_total(:, conditions_stanceR_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['right gluteus medius, right stance, 150ms, N = ' num2str(sum(conditions_stanceR_stimNeg_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rglutmed_normalized_total(:, conditions_stanceR_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['right tibialis anterior, right stance, 150ms, N = ' num2str(sum(conditions_stanceR_stimNeg_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rtibiant_normalized_total(:, conditions_stanceR_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['right peroneus longus, right stance, 150ms, N = ' num2str(sum(conditions_stanceR_stimNeg_150ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rperolng_normalized_total(:, conditions_stanceR_stimNeg_150ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), 'color', color_negative, 'linewidth', 5);
        end  
    
        % 450ms
        % pos
        if any(conditions_stanceL_stimPos_450ms)
            figure; axes; hold on; title(['left gluteus medius, left stance, 450ms, N = ' num2str(sum(conditions_stanceL_stimPos_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lglutmed_normalized_total(:, conditions_stanceL_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['left tibialis anterior, left stance, 450ms, N = ' num2str(sum(conditions_stanceL_stimPos_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, ltibiant_normalized_total(:, conditions_stanceL_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['left peroneus longus, left stance, 450ms, N = ' num2str(sum(conditions_stanceL_stimPos_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lperolng_normalized_total(:, conditions_stanceL_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['right gluteus medius, left stance, 450ms, N = ' num2str(sum(conditions_stanceL_stimPos_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rglutmed_normalized_total(:, conditions_stanceL_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['right tibialis anterior, left stance, 450ms, N = ' num2str(sum(conditions_stanceL_stimPos_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rtibiant_normalized_total(:, conditions_stanceL_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['right peroneus longus, left stance, 450ms, N = ' num2str(sum(conditions_stanceL_stimPos_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rperolng_normalized_total(:, conditions_stanceL_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);
        end

        if any(conditions_stanceR_stimPos_450ms)
            figure; axes; hold on; title(['left gluteus medius, right stance, 450ms, N = ' num2str(sum(conditions_stanceR_stimPos_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lglutmed_normalized_total(:, conditions_stanceR_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['left tibialis anterior, right stance, 450ms, N = ' num2str(sum(conditions_stanceR_stimPos_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, ltibiant_normalized_total(:, conditions_stanceR_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['left peroneus longus, right stance, 450ms, N = ' num2str(sum(conditions_stanceR_stimPos_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lperolng_normalized_total(:, conditions_stanceR_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['right gluteus medius, right stance, 450ms, N = ' num2str(sum(conditions_stanceR_stimPos_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rglutmed_normalized_total(:, conditions_stanceR_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['right tibialis anterior, right stance, 450ms, N = ' num2str(sum(conditions_stanceR_stimPos_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rtibiant_normalized_total(:, conditions_stanceR_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);

            figure; axes; hold on; title(['right peroneus longus, right stance, 450ms, N = ' num2str(sum(conditions_stanceR_stimPos_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rperolng_normalized_total(:, conditions_stanceR_stimPos_450ms), 'color', lightenColor(color_positive, 0.5));
            plot(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceR_stimPos_450ms), 2), 'color', color_positive, 'linewidth', 5);
        end

        % neg
        if any(conditions_stanceL_stimPos_450ms)
            figure; axes; hold on; title(['left gluteus medius, left stance, 450ms, N = ' num2str(sum(conditions_stanceL_stimNeg_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lglutmed_normalized_total(:, conditions_stanceL_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['left tibialis anterior, left stance, 450ms, N = ' num2str(sum(conditions_stanceL_stimNeg_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, ltibiant_normalized_total(:, conditions_stanceL_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['left peroneus longus, left stance, 450ms, N = ' num2str(sum(conditions_stanceL_stimNeg_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lperolng_normalized_total(:, conditions_stanceL_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['right gluteus medius, left stance, 450ms, N = ' num2str(sum(conditions_stanceL_stimNeg_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rglutmed_normalized_total(:, conditions_stanceL_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['right tibialis anterior, left stance, 450ms, N = ' num2str(sum(conditions_stanceL_stimNeg_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rtibiant_normalized_total(:, conditions_stanceL_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['right peroneus longus, left stance, 450ms, N = ' num2str(sum(conditions_stanceL_stimNeg_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rperolng_normalized_total(:, conditions_stanceL_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        end

        if any(conditions_stanceR_stimNeg_450ms)
            figure; axes; hold on; title(['left gluteus medius, right stance, 450ms, N = ' num2str(sum(conditions_stanceR_stimNeg_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lglutmed_normalized_total(:, conditions_stanceR_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['left tibialis anterior, right stance, 450ms, N = ' num2str(sum(conditions_stanceR_stimNeg_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, ltibiant_normalized_total(:, conditions_stanceR_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['left peroneus longus, right stance, 450ms, N = ' num2str(sum(conditions_stanceR_stimNeg_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, lperolng_normalized_total(:, conditions_stanceR_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['right gluteus medius, right stance, 450ms, N = ' num2str(sum(conditions_stanceR_stimNeg_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rglutmed_normalized_total(:, conditions_stanceR_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['right tibialis anterior, right stance, 450ms, N = ' num2str(sum(conditions_stanceR_stimNeg_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rtibiant_normalized_total(:, conditions_stanceR_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);

            figure; axes; hold on; title(['right peroneus longus, right stance, 450ms, N = ' num2str(sum(conditions_stanceR_stimNeg_450ms)) ' - ' subject_id]); set(gca, 'Fontsize', 12)
            plot(time_normalized, rperolng_normalized_total(:, conditions_stanceR_stimNeg_450ms), 'color', lightenColor(color_negative, 0.5));
            plot(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2), 'color', color_negative, 'linewidth', 5);
        end  
    end
end

if do_emg_plots_absolute_right
    % 0ms
    figure; axes; hold on; %title(['left gluteus medius, right stance, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(lglutmed_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), cinv(lglutmed_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), cinv(lglutmed_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%     this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'lglutmed_absolute_stanceR_0ms.eps', 'epsc2')
    end
    
    figure; axes; hold on; %title(['left tibialis anterior, right stance, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(ltibiant_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), cinv(ltibiant_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), cinv(ltibiant_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%     this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'ltibiant_absolute_stanceR_0ms.eps', 'epsc2')
    end
    
    figure; axes; hold on; %title(['left peroneus longus, right stance, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(lperolng_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), cinv(lperolng_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), cinv(lperolng_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%     this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'lperolng_absolute_stanceR_0ms.eps', 'epsc2')
    end

    figure; axes; hold on; %title(['right gluteus medius, right stance, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(rglutmed_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), cinv(rglutmed_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), cinv(rglutmed_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%     this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'rglutmed_absolute_stanceR_0ms.eps', 'epsc2')
    end
    
    figure; axes; hold on; %title(['right tibialis anterior, right stance, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(rtibiant_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), cinv(rtibiant_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), cinv(rtibiant_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%     this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'rtibiant_absolute_stanceR_0ms.eps', 'epsc2')
    end
    
    figure; axes; hold on; %title(['right peroneus longus, right stance, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(rperolng_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), cinv(rperolng_normalized_total(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), cinv(rperolng_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%     this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'rperolng_absolute_stanceR_0ms.eps', 'epsc2')
    end

    
    
    
    
    
    
    if do_phase_delayed_plots
        % 150ms
        figure; axes; hold on; %title(['left gluteus medius, right stance, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(lglutmed_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), cinv(lglutmed_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), cinv(lglutmed_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%         this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'lglutmed_absolute_stanceR_150ms.eps', 'epsc2')
        end

        figure; axes; hold on; %title(['left tibialis anterior, right stance, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(ltibiant_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), cinv(ltibiant_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), cinv(ltibiant_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%         this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'ltibiant_absolute_stanceR_150ms.eps', 'epsc2')
        end

        figure; axes; hold on; %title(['left peroneus longus, right stance, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(lperolng_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), cinv(lperolng_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), cinv(lperolng_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%         this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'lperolng_absolute_stanceR_150ms.eps', 'epsc2')
        end

        figure; axes; hold on; %title(['right gluteus medius, right stance, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(rglutmed_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), cinv(rglutmed_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), cinv(rglutmed_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%         this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'rglutmed_absolute_stanceR_150ms.eps', 'epsc2')
        end

        figure; axes; hold on; %title(['right tibialis anterior, right stance, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(rtibiant_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), cinv(rtibiant_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), cinv(rtibiant_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%         this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'rtibiant_absolute_stanceR_150ms.eps', 'epsc2')
        end

        figure; axes; hold on; %title(['right peroneus longus, right stance, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceR_stimNo), 2), cinv(rperolng_normalized_total(:, conditions_stanceR_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), cinv(rperolng_normalized_total(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), cinv(rperolng_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%         this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'rperolng_absolute_stanceR_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; %title(['left gluteus medius, right stance, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(lglutmed_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), cinv(lglutmed_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(lglutmed_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), cinv(lglutmed_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%         this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'lglutmed_absolute_stanceL_450ms.eps', 'epsc2')
        end

        figure; axes; hold on; %title(['left tibialis anterior, right stance, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(ltibiant_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), cinv(ltibiant_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(ltibiant_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), cinv(ltibiant_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%         this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'ltibiant_absolute_stanceL_450ms.eps', 'epsc2')
        end

        figure; axes; hold on; %title(['left peroneus longus, right stance, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(lperolng_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), cinv(lperolng_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(lperolng_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), cinv(lperolng_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%         this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'lperolng_absolute_stanceL_450ms.eps', 'epsc2')
        end

        figure; axes; hold on; %title(['right gluteus medius, right stance, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(rglutmed_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), cinv(rglutmed_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rglutmed_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), cinv(rglutmed_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%         this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'rglutmed_absolute_stanceL_450ms.eps', 'epsc2')
        end

        figure; axes; hold on; %title(['right tibialis anterior, right stance, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(rtibiant_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), cinv(rtibiant_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rtibiant_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), cinv(rtibiant_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%         this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'rtibiant_absolute_stanceL_450ms.eps', 'epsc2')
        end

        figure; axes; hold on; %title(['right peroneus longus, right stance, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceL_stimNo), 2), cinv(rperolng_normalized_total(:, conditions_stanceL_stimNo), 2), {'color', color_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), cinv(rperolng_normalized_total(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rperolng_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), cinv(rperolng_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]); 
%         this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'rperolng_absolute_stanceL_450ms.eps', 'epsc2')
        end
        
    end
end

if do_emg_plots_response_right
    % 0ms
    figure; axes; hold on; title(['left gluteus medius response, right stance, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(lglutmed_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), cinv(lglutmed_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(lglutmed_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), cinv(lglutmed_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'lglutmed_response_stanceR_0ms.eps', 'epsc2')
    end
    
    figure; axes; hold on; title(['left tibialis anterior response, right stance, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(ltibiant_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), cinv(ltibiant_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(ltibiant_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), cinv(ltibiant_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'ltibiant_response_stanceR_0ms.eps', 'epsc2')
    end
    
    figure; axes; hold on; title(['left peroneus longus response, right stance, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(lperolng_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), cinv(lperolng_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(lperolng_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), cinv(lperolng_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'lperolng_response_stanceR_0ms.eps', 'epsc2')
    end
    
    
    figure; axes; hold on; title(['right gluteus medius response, right stance, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(rglutmed_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), cinv(rglutmed_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(rglutmed_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), cinv(rglutmed_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'lglutmed_response_stanceR_0ms.eps', 'epsc2')
    end
    
    figure; axes; hold on; title(['right tibialis anterior response, right stance, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(rtibiant_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), cinv(rtibiant_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(rtibiant_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), cinv(rtibiant_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'ltibiant_response_stanceR_0ms.eps', 'epsc2')
    end
    
    figure; axes; hold on; title(['right peroneus longus response, right stance, 0ms - ' subject_id]); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, mean(rperolng_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), cinv(rperolng_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, mean(rperolng_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), cinv(rperolng_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right'); set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%     text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    if save_figures
        saveas(gcf, 'lperolng_response_stanceR_0ms.eps', 'epsc2')
    end
    
    
    
    
    if do_phase_delayed_plots
        % 150ms
        figure; axes; hold on; title(['left gluteus medius response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(lglutmed_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), cinv(lglutmed_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(lglutmed_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), cinv(lglutmed_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'lglutmed_response_stanceR_150ms.eps', 'epsc2')
        end

        figure; axes; hold on; title(['left tibialis anterior response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(ltibiant_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), cinv(ltibiant_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(ltibiant_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), cinv(ltibiant_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'ltibiant_response_stanceR_150ms.eps', 'epsc2')
        end

        figure; axes; hold on; title(['left peroneus longus response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(lperolng_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), cinv(lperolng_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(lperolng_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), cinv(lperolng_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'lperolng_response_stanceR_150ms.eps', 'epsc2')
        end

        figure; axes; hold on; title(['right gluteus medius response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(rglutmed_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), cinv(rglutmed_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rglutmed_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), cinv(rglutmed_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'rglutmed_response_stanceR_150ms.eps', 'epsc2')
        end

        figure; axes; hold on; title(['right tibialis anterior response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(rtibiant_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), cinv(rtibiant_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rtibiant_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), cinv(rtibiant_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'rtibiant_response_stanceR_150ms.eps', 'epsc2')
        end

        figure; axes; hold on; title(['right peroneus longus response, 150ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(rperolng_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), cinv(rperolng_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rperolng_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), cinv(rperolng_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'rperolng_response_stanceR_150ms.eps', 'epsc2')
        end

        % 450ms
        figure; axes; hold on; title(['left gluteus medius response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(lglutmed_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), cinv(lglutmed_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(lglutmed_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), cinv(lglutmed_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'lglutmed_response_stanceL_450ms.eps', 'epsc2')
        end

        figure; axes; hold on; title(['left tibialis anterior response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(ltibiant_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), cinv(ltibiant_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(ltibiant_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), cinv(ltibiant_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'ltibiant_response_stanceL_450ms.eps', 'epsc2')
        end

        figure; axes; hold on; title(['left peroneus longus response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(lperolng_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), cinv(lperolng_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(lperolng_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), cinv(lperolng_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'lperolng_response_stanceL_450ms.eps', 'epsc2')
        end

        figure; axes; hold on; title(['right gluteus medius response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(rglutmed_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), cinv(rglutmed_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rglutmed_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), cinv(rglutmed_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'rglutmed_response_stanceL_450ms.eps', 'epsc2')
        end

        figure; axes; hold on; title(['right tibialis anterior response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(rtibiant_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), cinv(rtibiant_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rtibiant_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), cinv(rtibiant_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'rtibiant_response_stanceL_450ms.eps', 'epsc2')
        end

        figure; axes; hold on; title(['right peroneus longus response, 450ms - ' subject_id]); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, mean(rperolng_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), cinv(rperolng_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2), {'color', color_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, mean(rperolng_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), cinv(rperolng_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2), {'color', color_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [time_normalized(1) time_normalized(end)]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
%         text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        if save_figures
            saveas(gcf, 'rperolng_response_stanceL_450ms.eps', 'epsc2')
        end
        
        
        
        
        
        
    end
end

%% do_step_response_plots
if do_step_response_plots
    
    
    figure; difference_axes = axes('Fontsize', 16); hold on
    linewidth = 1;
    civ_marker_height = .05;
    
    y_loc = 5.1;
    step_response_bar_pos_0ms = barh(y_loc, mean(step_response_total(conditions_stanceR_stimPos_0ms)), 'facecolor', color_positive);
    plot(mean(step_response_total(conditions_stanceR_stimPos_0ms)) + cinv(step_response_total(conditions_stanceR_stimPos_0ms))*[-1 1], [y_loc y_loc], 'k', 'linewidth', linewidth);
    plot(mean(step_response_total(conditions_stanceR_stimPos_0ms)) + cinv(step_response_total(conditions_stanceR_stimPos_0ms))*[-1 -1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    plot(mean(step_response_total(conditions_stanceR_stimPos_0ms)) + cinv(step_response_total(conditions_stanceR_stimPos_0ms))*[1 1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);

    y_loc = 5.9;
    step_response_bar_neg_0ms = barh(y_loc, mean(step_response_total(conditions_stanceR_stimNeg_0ms)), 'facecolor', color_negative);
    plot(mean(step_response_total(conditions_stanceR_stimNeg_0ms)) + cinv(step_response_total(conditions_stanceR_stimNeg_0ms))*[-1 1], [y_loc y_loc], 'k', 'linewidth', linewidth);
    plot(mean(step_response_total(conditions_stanceR_stimNeg_0ms)) + cinv(step_response_total(conditions_stanceR_stimNeg_0ms))*[-1 -1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    plot(mean(step_response_total(conditions_stanceR_stimNeg_0ms)) + cinv(step_response_total(conditions_stanceR_stimNeg_0ms))*[1 1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    
    y_loc = 3.1;
    step_response_bar_pos_150ms = barh(y_loc, mean(step_response_total(conditions_stanceR_stimPos_150ms)), 'facecolor', color_positive);
    plot(mean(step_response_total(conditions_stanceR_stimPos_150ms)) + cinv(step_response_total(conditions_stanceR_stimPos_150ms))*[-1 1], [y_loc y_loc], 'k', 'linewidth', linewidth);
    plot(mean(step_response_total(conditions_stanceR_stimPos_150ms)) + cinv(step_response_total(conditions_stanceR_stimPos_150ms))*[-1 -1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    plot(mean(step_response_total(conditions_stanceR_stimPos_150ms)) + cinv(step_response_total(conditions_stanceR_stimPos_150ms))*[1 1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);

    y_loc = 3.9;
    step_response_bar_neg_150ms = barh(y_loc, mean(step_response_total(conditions_stanceR_stimNeg_150ms)), 'facecolor', color_negative);
    plot(mean(step_response_total(conditions_stanceR_stimNeg_150ms)) + cinv(step_response_total(conditions_stanceR_stimNeg_150ms))*[-1 1], [y_loc y_loc], 'k', 'linewidth', linewidth);
    plot(mean(step_response_total(conditions_stanceR_stimNeg_150ms)) + cinv(step_response_total(conditions_stanceR_stimNeg_150ms))*[-1 -1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    plot(mean(step_response_total(conditions_stanceR_stimNeg_150ms)) + cinv(step_response_total(conditions_stanceR_stimNeg_150ms))*[1 1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    
    y_loc = 1.1;
    step_response_bar_pos_450ms = barh(y_loc, mean(step_response_total(conditions_stanceL_stimPos_450ms)), 'facecolor', color_positive);
    plot(mean(step_response_total(conditions_stanceL_stimPos_450ms)) + cinv(step_response_total(conditions_stanceL_stimPos_450ms))*[-1 1], [y_loc y_loc], 'k', 'linewidth', linewidth);
    plot(mean(step_response_total(conditions_stanceL_stimPos_450ms)) + cinv(step_response_total(conditions_stanceL_stimPos_450ms))*[-1 -1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    plot(mean(step_response_total(conditions_stanceL_stimPos_450ms)) + cinv(step_response_total(conditions_stanceL_stimPos_450ms))*[1 1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);

    y_loc = 1.9;
    step_response_bar_neg_450ms = barh(y_loc, mean(step_response_total(conditions_stanceL_stimNeg_450ms)), 'facecolor', color_negative);
    plot(mean(step_response_total(conditions_stanceL_stimNeg_450ms)) + cinv(step_response_total(conditions_stanceL_stimNeg_450ms))*[-1 1], [y_loc y_loc], 'k', 'linewidth', linewidth);
    plot(mean(step_response_total(conditions_stanceL_stimNeg_450ms)) + cinv(step_response_total(conditions_stanceL_stimNeg_450ms))*[-1 -1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    plot(mean(step_response_total(conditions_stanceL_stimNeg_450ms)) + cinv(step_response_total(conditions_stanceL_stimNeg_450ms))*[1 1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    
    set(difference_axes, 'ytick', [1.5, 3.5 5.5]);
    set(difference_axes, 'yticklabel', {'450ms', '150ms', '0ms'});
    set(gca, 'xlim', [-0.022 0.022])
%     set(difference_axes, 'xticklabel', [], 'yticklabel', []);
    
%     difference_axes.YAxis.Visible = 'off';
%     set(difference_axes, 'Fontsize', 16);
%     difference_axes.YTickLabelRotation = 0;    
    if save_figures
        saveas(gcf, 'step_response.eps', 'epsc2')
    end
    
    
    % stim response
    figure; difference_axes = axes('Fontsize', 16); hold on
    linewidth = 1;
    civ_marker_height = .05;
    
    y_loc = 1.1;
    stim_response_bar_pos_0ms = barh(y_loc, mean(stim_response_total(conditions_stanceR_stimPos_0ms)), 'facecolor', color_positive);
    plot(mean(stim_response_total(conditions_stanceR_stimPos_0ms)) + cinv(stim_response_total(conditions_stanceR_stimPos_0ms))*[-1 1], [y_loc y_loc], 'k', 'linewidth', linewidth);
    plot(mean(stim_response_total(conditions_stanceR_stimPos_0ms)) + cinv(stim_response_total(conditions_stanceR_stimPos_0ms))*[-1 -1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    plot(mean(stim_response_total(conditions_stanceR_stimPos_0ms)) + cinv(stim_response_total(conditions_stanceR_stimPos_0ms))*[1 1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);

    y_loc = 1.9;
    stim_response_bar_neg_0ms = barh(y_loc, mean(stim_response_total(conditions_stanceR_stimNeg_0ms)), 'facecolor', color_negative);
    plot(mean(stim_response_total(conditions_stanceR_stimNeg_0ms)) + cinv(stim_response_total(conditions_stanceR_stimNeg_0ms))*[-1 1], [y_loc y_loc], 'k', 'linewidth', linewidth);
    plot(mean(stim_response_total(conditions_stanceR_stimNeg_0ms)) + cinv(stim_response_total(conditions_stanceR_stimNeg_0ms))*[-1 -1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    plot(mean(stim_response_total(conditions_stanceR_stimNeg_0ms)) + cinv(stim_response_total(conditions_stanceR_stimNeg_0ms))*[1 1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    
    y_loc = 3.1;
    stim_response_bar_pos_150ms = barh(y_loc, mean(stim_response_total(conditions_stanceR_stimPos_150ms)), 'facecolor', color_positive);
    plot(mean(stim_response_total(conditions_stanceR_stimPos_150ms)) + cinv(stim_response_total(conditions_stanceR_stimPos_150ms))*[-1 1], [y_loc y_loc], 'k', 'linewidth', linewidth);
    plot(mean(stim_response_total(conditions_stanceR_stimPos_150ms)) + cinv(stim_response_total(conditions_stanceR_stimPos_150ms))*[-1 -1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    plot(mean(stim_response_total(conditions_stanceR_stimPos_150ms)) + cinv(stim_response_total(conditions_stanceR_stimPos_150ms))*[1 1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);

    y_loc = 3.9;
    stim_response_bar_neg_150ms = barh(y_loc, mean(stim_response_total(conditions_stanceR_stimNeg_150ms)), 'facecolor', color_negative);
    plot(mean(stim_response_total(conditions_stanceR_stimNeg_150ms)) + cinv(stim_response_total(conditions_stanceR_stimNeg_150ms))*[-1 1], [y_loc y_loc], 'k', 'linewidth', linewidth);
    plot(mean(stim_response_total(conditions_stanceR_stimNeg_150ms)) + cinv(stim_response_total(conditions_stanceR_stimNeg_150ms))*[-1 -1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    plot(mean(stim_response_total(conditions_stanceR_stimNeg_150ms)) + cinv(stim_response_total(conditions_stanceR_stimNeg_150ms))*[1 1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    
    y_loc = 5.1;
    stim_response_bar_pos_450ms = barh(y_loc, mean(stim_response_total(conditions_stanceL_stimPos_450ms)), 'facecolor', color_positive);
    plot(mean(stim_response_total(conditions_stanceL_stimPos_450ms)) + cinv(stim_response_total(conditions_stanceL_stimPos_450ms))*[-1 1], [y_loc y_loc], 'k', 'linewidth', linewidth);
    plot(mean(stim_response_total(conditions_stanceL_stimPos_450ms)) + cinv(stim_response_total(conditions_stanceL_stimPos_450ms))*[-1 -1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    plot(mean(stim_response_total(conditions_stanceL_stimPos_450ms)) + cinv(stim_response_total(conditions_stanceL_stimPos_450ms))*[1 1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);

    y_loc = 5.9;
    stim_response_bar_neg_450ms = barh(y_loc, mean(stim_response_total(conditions_stanceL_stimNeg_450ms)), 'facecolor', color_negative);
    plot(mean(stim_response_total(conditions_stanceL_stimNeg_450ms)) + cinv(stim_response_total(conditions_stanceL_stimNeg_450ms))*[-1 1], [y_loc y_loc], 'k', 'linewidth', linewidth);
    plot(mean(stim_response_total(conditions_stanceL_stimNeg_450ms)) + cinv(stim_response_total(conditions_stanceL_stimNeg_450ms))*[-1 -1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    plot(mean(stim_response_total(conditions_stanceL_stimNeg_450ms)) + cinv(stim_response_total(conditions_stanceL_stimNeg_450ms))*[1 1], y_loc + [-1 1]*civ_marker_height, 'k', 'linewidth', linewidth);
    
    set(difference_axes, 'ytick', [1.5, 3.5 5.5]);
    set(difference_axes, 'yticklabel', {'0ms', '150ms', '450ms'});
    
    
    
    
%     group_order = ...
%       { ...
%         'RIGHT, NEGATIVE, 0ms', ...
%         'RIGHT, POSITIVE, 0ms', ...
%         'RIGHT, NEGATIVE, 150ms', ...
%         'RIGHT, POSITIVE, 150ms', ...
%         'LEFT, NEGATIVE, 450ms', ...
%         'LEFT, POSITIVE, 450ms' ...
%       }
%     group_order = ...
%       { ...
%         'RIGHT', 'LEFT' ...
%       };
%     
%     step_response_figure = figure; axes('fontsize', 12); hold on; 
%     my_boxplot = ...
%       boxplot ...
%         ( ...
%           step_response_total(conditions_stim), ...
%           {condition_stance_foot_list_total(conditions_stim)}, ...
%           'GroupOrder', group_order ...
%         );
%     
    
%           {condition_stance_foot_list_total(conditions_stim) condition_polarity_list_total(conditions_stim) condition_delay_list_total(conditions_stim)}, ...
    
    
    

    
end