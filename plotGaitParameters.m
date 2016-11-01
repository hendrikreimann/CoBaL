
plot_detailed       = 0;
plot_shaded         = 1;

show_legend         = 1;
dictate_axes        = 0;

save_figures        = 1;

error_shades = 'cinv';
error_shades = 'std';

continuous_variable_info = {};
discrete_variable_info = {};

% markers
% continuous_variable_info = [continuous_variable_info; {'linclination_normalized_all', 'left arm inclination angle', 'angle (deg)', 'linclination'}];
% continuous_variable_info = [continuous_variable_info; {'rinclination_normalized_all', 'right arm inclination angle', 'angle (deg)', 'rinclination'}];

discrete_variable_info = [discrete_variable_info; {'step_width_all', 'step width'}];
discrete_variable_info = [discrete_variable_info; {'step_length_all', 'step length'}];

load subjectInfo.mat;
load(['analysis' filesep makeFileName(date, subject_id, 'resultsConditions')]);
load(['analysis' filesep makeFileName(date, subject_id, 'resultsBalance')]);
load(['analysis' filesep makeFileName(date, subject_id, 'resultsArmswing')]);
load(['analysis' filesep makeFileName(date, subject_id, 'resultsForceplate')]);



%% define comparisons
%
% specify which comparisons should be made by giving a condition label.
% All different conditions of that label will then be plotted in the same figure

% define conditions
% order of condition is STANCE_FOOT, PERTURBATION, DELAY, INDEX
% STANCE_FOOT: RIGHT, LEFT
% PERTURBATION: POSITIVE, NEGATIVE
% DELAY: 0, 150, 450
% INDEX: ONE, TWO
%
% this should be loaded from a file

condition_labels = {'stance foot', 'experimental'};
conditions_to_analyze = ...
  {
    'STANCE_LEFT', 'baseline'; ...
    'STANCE_RIGHT', 'baseline'; ...
    'STANCE_LEFT', 'feedback'; ...
    'STANCE_RIGHT', 'feedback'; ...
    'STANCE_LEFT', 'post'; ...
    'STANCE_RIGHT', 'post'; ...
  };

number_of_conditions_to_analyze = size(conditions_to_analyze, 1);

% define comparisons
comparison_to_make = 2; % experimental only
comparison_indices = {};
% control_conditions_for_comparisons = [];
conditions_already_compared = [];
while length(conditions_already_compared) < number_of_conditions_to_analyze
    % start with the first available condition
    i_condition = 1;
    while ismember(i_condition, conditions_already_compared)
        i_condition = i_condition + 1;
    end
    
    this_comparison = i_condition;
    % search for conditions that differ from this one only in the 
    for j_condition = 1 : number_of_conditions_to_analyze
        if i_condition ~= j_condition
            % check which conditions labels agree between these two conditions
            comparison_table = zeros(1, length(condition_labels));
            for i_label = 1 : length(condition_labels)
                comparison_table(i_label) = strcmp(conditions_to_analyze{i_condition, i_label}, conditions_to_analyze{j_condition, i_label});
            end

            comparison_table_relevant = comparison_table;
            comparison_table_relevant(comparison_to_make) = [];
            if all(comparison_table_relevant)
                this_comparison = [this_comparison, j_condition];
            end
        end
    end
    comparison_indices = [comparison_indices; this_comparison];
    conditions_already_compared = [conditions_already_compared this_comparison];
end


%% do plots
color_control = [0.3 0.1 1];
colors_comparison = ...
  [ ...
    [1 0.3 0.1] * 0.7; ...
    [0.3 1 0.1] * 0.7; ...
    [0.1 0.3 1] * 0.7; ...
  ]; % should have one row per condition in the comparison


%% plot detailed
if plot_detailed
    % single variables
    for i_variable = 1 : size(discrete_variable_info, 1)
        evalstring = ['points_to_plot = ' discrete_variable_info{i_variable, 1} ';'];
        eval(evalstring);
        for i_comparison = 1 : length(comparison_indices);
            figure; axes; hold on;
            this_comparison = comparison_indices{i_comparison};
            for i_condition = 1 : length(this_comparison)
                condition_indicator = conditions_to_analyze_indicators(:, this_comparison(i_condition));
                label_string = conditions_to_analyze{comparison_indices{i_comparison}(i_condition), comparison_to_make};
                histogram ...
                  ( ...
                    points_to_plot(:, condition_indicator), ...
                    'edgecolor', colors_comparison(i_condition, :), ...
                    'facecolor', lightenColor(colors_comparison(i_condition, :), 0.5), ...
                    'DisplayName', label_string ...
                  )
            end

            % annotate
            if show_legend
                legend('toggle')
            end
            title_string = continuous_variable_info{i_variable, 2};
            for i_label = 1 : length(condition_labels);
                if i_label ~= comparison_to_make
                    title_string = [title_string ' - ' conditions_to_analyze{comparison_indices{i_comparison}(1), i_label}];
                end
            end
            title(title_string); set(gca, 'Fontsize', 12)
        end
    end    
    
    
    % continuous variables
    for i_variable = 1 : size(continuous_variable_info, 1)
        evalstring = ['trajectories_to_plot = ' continuous_variable_info{i_variable, 1} ';'];
        eval(evalstring);
        for i_comparison = 1 : length(comparison_indices);
            figure; axes; hold on;
            this_comparison = comparison_indices{i_comparison};
            condition_mean_plots = zeros(1, length(this_comparison));
            for i_condition = 1 : length(this_comparison)
                condition_indicator = conditions_to_analyze_indicators(:, this_comparison(i_condition));
                plot ...
                  ( ...
                    time_normalized, trajectories_to_plot(:, condition_indicator), ...
                    'HandleVisibility', 'off', ...
                    'color', lightenColor(colors_comparison(i_condition, :), 0.5) ...
                  );
                label_string = conditions_to_analyze{comparison_indices{i_comparison}(i_condition), comparison_to_make};
                condition_mean_plots(i_condition) = plot ...
                  ( ...
                    time_normalized, mean(trajectories_to_plot(:, condition_indicator), 2), ...
                    'DisplayName', label_string, ...
                    'linewidth', 5, ...
                    'color', colors_comparison(i_condition, :) ...
                  );

            end

            % reorder
            for i_condition = 1 : length(this_comparison)
                uistack(condition_mean_plots(i_condition), 'top');
            end

            % annotate
            if show_legend
                legend('toggle')
            end
            title_string = continuous_variable_info{i_variable, 2};
            for i_label = 1 : length(condition_labels);
                if i_label ~= comparison_to_make
                    title_string = [title_string ' - ' conditions_to_analyze{comparison_indices{i_comparison}(1), i_label}];
                end
            end
            title(title_string); set(gca, 'Fontsize', 12)
        end
    end
end

%% plot shaded
if plot_shaded
    % single variables
    for i_variable = 1 : size(discrete_variable_info, 1)
        evalstring = ['points_to_plot = ' discrete_variable_info{i_variable, 1} ';'];
        eval(evalstring);
        
        
        for i_comparison = 1 : length(comparison_indices);
            this_comparison = comparison_indices{i_comparison};
            % make condition labels for box plot
            condition_labels_for_boxplot = cell(length(points_to_plot), 1);
            for i_condition = 1 : length(this_comparison)
                % figure out condition label
                condition_indicator = conditions_to_analyze_indicators(:, this_comparison(i_condition));
                label_string = conditions_to_analyze{comparison_indices{i_comparison}(i_condition), comparison_to_make};
                % place condition label into label cell array
                for i_point = 1 : length(condition_labels_for_boxplot)
                    if condition_indicator(i_point)
                        condition_labels_for_boxplot{i_point} = label_string;
                    end
                end
            end            
            % prune data points that we don't want to look at in this plot
            conditions_for_this_comparison = conditions_to_analyze_indicators(:, this_comparison);
            indices_for_this_comparison = any(conditions_for_this_comparison, 2);
            conditions_pruned = conditions_for_this_comparison(indices_for_this_comparison, :);
            points_to_plot_pruned = points_to_plot(indices_for_this_comparison);
            condition_labels_for_boxplot_pruned = condition_labels_for_boxplot(indices_for_this_comparison);
            
            figure; axes; hold on;
            boxplot(points_to_plot_pruned, condition_labels_for_boxplot_pruned);

            % annotate
            title_string = discrete_variable_info{i_variable, 2};
            for i_label = 1 : length(condition_labels);
                if i_label ~= comparison_to_make
                    title_string = [title_string ' - ' strrep(conditions_to_analyze{comparison_indices{i_comparison}(1), i_label}, '_', ' ')];
                    filename_string = [filename_string '_' conditionStringToFilename(conditions_to_analyze{comparison_indices{i_comparison}(1), i_label})];
                end
            end
            title(title_string, 'interpreter', 'LaTeX'); set(gca, 'Fontsize', 12)
        end
    end      
    
    for i_variable = 1 : size(continuous_variable_info, 1)
        evalstring = ['trajectories_to_plot = ' continuous_variable_info{i_variable, 1} ';'];
        eval(evalstring);
        for i_comparison = 1 : length(comparison_indices);
            legend_handles = [];
            legend_data = {};
            figure; axes; hold on;
            % plot stimulus
            this_comparison = comparison_indices{i_comparison};
            condition_mean_plots = zeros(1, length(this_comparison));
            for i_condition = 1 : length(this_comparison)
                condition_indicator = conditions_to_analyze_indicators(:, this_comparison(i_condition));
                mean_curves = mean(trajectories_to_plot(:, condition_indicator), 2);
                if strcmp(error_shades, 'cinv')
                    error_curves = cinv(trajectories_to_plot(:, condition_indicator), 2);
                elseif strcmp(error_shades, 'std')
                    error_curves = std(trajectories_to_plot(:, condition_indicator), 1, 2);
                end
                current_plots = shadedErrorBar ...
                  ( ...
                    time_normalized, ...
                    mean_curves, ...
                    error_curves, ...
                    { ...
                      'color', colors_comparison(i_condition, :), ...
                      'linewidth', 3 ...
                    }, ...
                    1 ...
                  );
                legend_handles = [legend_handles, current_plots.mainLine];
                label_string = strrep(conditions_to_analyze{comparison_indices{i_comparison}(i_condition), comparison_to_make}, '_', ' ');
                legend_data = [legend_data, label_string];
                
            end

            % annotate
            if show_legend
                this_legend = legend(legend_handles, legend_data);
            end
            title_string = continuous_variable_info{i_variable, 2};
            filename_string = continuous_variable_info{i_variable, 4};
            for i_label = 1 : length(condition_labels);
                if i_label ~= comparison_to_make
%                     title_string = [title_string ' - ' conditionStringToTitle(conditions_to_analyze{comparison_indices{i_comparison}(1), i_label})];
                    title_string = [title_string ' - ' strrep(conditions_to_analyze{comparison_indices{i_comparison}(1), i_label}, '_', ' ')];
%                     title_string = strrep(title_string, '_', ' ');
                    filename_string = [filename_string '_' conditionStringToFilename(conditions_to_analyze{comparison_indices{i_comparison}(1), i_label})];
                end
            end
            title(title_string, 'interpreter', 'LaTeX'); set(gca, 'Fontsize', 12)
            
            set(gca, 'xlim', [time_normalized(1), time_normalized(end)]);
            xlabel('normalized time (s)');
            ylabel(continuous_variable_info{i_variable, 3});
            
            if dictate_axes
                set(gca, 'xlim', [time_normalized(1), time_normalized(end)]);
                set(gca, 'ylim', [-continuous_variable_info{i_variable, 5}, continuous_variable_info{i_variable, 5}]);
            end
            
            xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
%             text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right $\rightarrow$', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right', 'interpreter', 'LaTeX')
%             text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '$\leftarrow$ left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left', 'interpreter', 'LaTeX')
            
            % save
            if save_figures
                
                % figure out folders
                if ~exist('figures', 'dir')
                    mkdir('figures')
                end
                
                filename = ['figures' filesep filename_string '.eps'];
                saveas(gcf, filename, 'epsc2')
                
            end
        end
    end    
    
    
    
    
    
end











if false
    check_out_origin                        = 0;

    do_cop_plots_single                     = 0;
    do_cop_plots_absolute                   = 0;
    do_step_parameter_plots                 = 0;

    do_cop_plots_experimental               = 0;
    do_cop_std_plots_experimental           = 1;
    do_step_parameter_plots_experimental    = 0;


    save_figures                            = 1;

    load subjectInfo.mat;
    load(makeFileName(date, subject_id, 'gaitParametersConditions'));
    load(makeFileName(date, subject_id, 'gaitParametersForceplate'));
    load(makeFileName(date, subject_id, 'gaitParametersMarker'));

    color_control = [0.3 0.1 1];
    color_right = [1 0.3 .1] * 0.7;
    color_left = [0.3 1 0.1] * 0.7;

    color_list = ...
      [ ...
             0    0.4470    0.7410
        0.8500    0.3250    0.0980
        0.9290    0.6940    0.1250
        0.4940    0.1840    0.5560
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840
      ];


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
        left_plot = shadedErrorBar(time_normalized, mean(lcop_x_normalized_total(:, conditions_stanceL), 2), std(lcop_x_normalized_total(:, conditions_stanceL), 0, 2), {'color', color_left, 'linewidth', 5}, 1);
        right_plot = shadedErrorBar(time_normalized, mean(rcop_x_normalized_total(:, conditions_stanceR), 2), std(rcop_x_normalized_total(:, conditions_stanceR), 0, 2), {'color', color_right, 'linewidth', 5}, 1);
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

    %% experimental condition

    if do_cop_plots_experimental
        number_of_experimental_conditions = length(experimental_conditions);

        % left
        figure; axes; hold on; title('LEFT foot medial-lateral CoP - relative to stance foot'); set(gca, 'Fontsize', 12)
        left_plots = zeros(1, number_of_experimental_conditions);
        for i_condition = 1 : number_of_experimental_conditions
            new_plot = shadedErrorBar(time_normalized, mean(lcop_x_normalized_total(:, strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, experimental_conditions{i_condition})), 2), std(lcop_x_normalized_total(:, strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, experimental_conditions{i_condition})), 0, 2), {'color', color_list(i_condition, :), 'linewidth', 5}, 1);
            left_plots(i_condition) = new_plot.mainLine;
        end
        legend(left_plots, experimental_conditions);
        if save_figures
            saveas(gcf, 'cop_left.eps', 'epsc2')
        end

        % right
        figure; axes; hold on; title('RIGHT foot medial-lateral CoP - relative to stance foot'); set(gca, 'Fontsize', 12)
        right_plots = zeros(1, number_of_experimental_conditions);
        for i_condition = 1 : number_of_experimental_conditions
            new_plot = shadedErrorBar(time_normalized, mean(rcop_x_normalized_total(:, strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, experimental_conditions{i_condition})), 2), std(rcop_x_normalized_total(:, strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, experimental_conditions{i_condition})), 0, 2), {'color', color_list(i_condition, :), 'linewidth', 5}, 1);
            right_plots(i_condition) = new_plot.mainLine;
        end
        legend(right_plots, experimental_conditions);
        if save_figures
            saveas(gcf, 'cop_right.eps', 'epsc2')
        end
    end

    if do_cop_std_plots_experimental

        % run quick and dirty f-tests
        cop_left_on = lcop_x_normalized_total(:, strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, 'ON'));
        cop_left_off = lcop_x_normalized_total(:, strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, 'OFF'));
        cop_left_post = lcop_x_normalized_total(:, strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, 'POST'));
        cop_right_on = rcop_x_normalized_total(:, strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, 'ON'));
        cop_right_off = rcop_x_normalized_total(:, strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, 'OFF'));
        cop_right_post = rcop_x_normalized_total(:, strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, 'POST'));
        p_left_cop_onOff = zeros(number_of_time_steps_normalized, 1);
        for i_time = 1 : number_of_time_steps_normalized
            [h, p_left_cop_onOff(i_time)] = vartest2(cop_left_on(i_time, :), cop_left_off(i_time, :));


        end

        % t-tests with FDR limitation after Benyamini and Hochberg (1995)
        false_discovery_rate = 0.05;
        p_left_cop_onOff = zeros(number_of_time_steps_normalized, 1);
        for i_time = 1 : number_of_time_steps_normalized
            [~, p_left_cop_onOff(i_time)] = vartest2(cop_left_on(i_time, :), cop_left_off(i_time, :));
        end
        first_time_step_to_test = 1;
        time_steps_to_test = first_time_step_to_test : number_of_time_steps_normalized;
        number_of_hypotheses_cop_std = length(p_left_cop_onOff);
        h_indices_cop_ml = 1 : number_of_hypotheses_cop_std;
        hypothesis_matrix_cop_ml = [p_left_cop_onOff, h_indices_cop_ml'];
        hypothesis_matrix_ascending_cop_ml = sortrows(hypothesis_matrix_cop_ml, 1);
        decision_threshold_cop_ml = (1:number_of_hypotheses_cop_std)' * 1/number_of_hypotheses_cop_std * false_discovery_rate;
        hypothesis_rejection_index_cop_ml = find(hypothesis_matrix_ascending_cop_ml(:, 1) < decision_threshold_cop_ml, 1, 'last');
        h_results_cop_ml = [ones(1, hypothesis_rejection_index_cop_ml) zeros(1, number_of_hypotheses_cop_std-hypothesis_rejection_index_cop_ml)];
        result_matrix_ascending_cop_ml = [hypothesis_matrix_ascending_cop_ml h_results_cop_ml'];
        result_matrix_left_cop_onOff = sortrows(result_matrix_ascending_cop_ml, 2);
        cop_left_ftest_result_onOff = [-1*ones(1, first_time_step_to_test-1) result_matrix_left_cop_onOff(1:number_of_hypotheses_cop_std, 3)']';

        % t-tests with FDR limitation after Benyamini and Hochberg (1995)
        false_discovery_rate = 0.05;
        p_right_cop_onOff = zeros(number_of_time_steps_normalized, 1);
        for i_time = 1 : number_of_time_steps_normalized
            [~, p_right_cop_onOff(i_time)] = vartest2(cop_right_on(i_time, :), cop_right_off(i_time, :));
        end
        first_time_step_to_test = 1;
        time_steps_to_test = first_time_step_to_test : number_of_time_steps_normalized;
        number_of_hypotheses_cop_std = length(p_right_cop_onOff);
        h_indices_cop_ml = 1 : number_of_hypotheses_cop_std;
        hypothesis_matrix_cop_ml = [p_right_cop_onOff, h_indices_cop_ml'];
        hypothesis_matrix_ascending_cop_ml = sortrows(hypothesis_matrix_cop_ml, 1);
        decision_threshold_cop_ml = (1:number_of_hypotheses_cop_std)' * 1/number_of_hypotheses_cop_std * false_discovery_rate;
        hypothesis_rejection_index_cop_ml = find(hypothesis_matrix_ascending_cop_ml(:, 1) < decision_threshold_cop_ml, 1, 'last');
        h_results_cop_ml = [ones(1, hypothesis_rejection_index_cop_ml) zeros(1, number_of_hypotheses_cop_std-hypothesis_rejection_index_cop_ml)];
        result_matrix_ascending_cop_ml = [hypothesis_matrix_ascending_cop_ml h_results_cop_ml'];
        result_matrix_right_cop_onOff = sortrows(result_matrix_ascending_cop_ml, 2);
        cop_right_ftest_result_onOff = [-1*ones(1, first_time_step_to_test-1) result_matrix_right_cop_onOff(1:number_of_hypotheses_cop_std, 3)']';


        significance_marker_height_onOff_fdr = 0.002;


        number_of_experimental_conditions = length(experimental_conditions);

        % left
        figure; axes; hold on; title('Standard deviation of the LEFT foot medial-lateral CoP'); set(gca, 'Fontsize', 12)
        left_plots = zeros(1, number_of_experimental_conditions);
        for i_condition = 1 : number_of_experimental_conditions
            left_plots(i_condition) = plot(time_normalized, std(lcop_x_normalized_total(:, strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, experimental_conditions{i_condition})), 0, 2), 'color', color_list(i_condition, :), 'linewidth', 5);
        end
        legend(left_plots, experimental_conditions);
        plot(time_normalized(cop_left_ftest_result_onOff==1), significance_marker_height_onOff_fdr, '*', 'color', color_list(4, :), 'markersize', 8, 'linewidth', 8);
        if save_figures
            saveas(gcf, 'cop_std_left.eps', 'epsc2')
        end

        % right
        figure; axes; hold on; title('Standard deviation of the RIGHT foot medial-lateral CoP'); set(gca, 'Fontsize', 12)
        right_plots = zeros(1, number_of_experimental_conditions);
        for i_condition = 1 : number_of_experimental_conditions
            right_plots(i_condition) = plot(time_normalized, std(rcop_x_normalized_total(:, strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, experimental_conditions{i_condition})), 0, 2), 'color', color_list(i_condition, :), 'linewidth', 5);
        end
        legend(right_plots, experimental_conditions);
        plot(time_normalized(cop_right_ftest_result_onOff==1), significance_marker_height_onOff_fdr, '*', 'color', color_list(4, :), 'markersize', 8, 'linewidth', 8);
        if save_figures
            saveas(gcf, 'cop_std_right.eps', 'epsc2')
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

    %% step parameter plots
    if do_step_parameter_plots_experimental
    %     % width histogram
    %     bin_edges = linspace(-0.15, 0.15, 39);
    %     step_width_figure = figure; axes('fontsize', 12); hold on; title('step width', 'fontsize', 16);
    %     step_width_left_histograms = zeros(1, number_of_experimental_conditions);
    %     step_width_right_histograms = zeros(1, number_of_experimental_conditions);
    %     for i_condition = 1 : number_of_experimental_conditions
    %        step_width_left_histograms(i_condition) = histogram(step_width_total(strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, experimental_conditions{i_condition})), 'Normalization', 'pdf', 'displaystyle', 'bar', 'edgecolor', color_list(i_condition, :), 'facecolor', lightenColor(color_list(i_condition, :), 0.5));        
    %        step_width_right_histograms(i_condition) = histogram(step_width_total(strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, experimental_conditions{i_condition})), 'Normalization', 'pdf', 'displaystyle', 'bar', 'edgecolor', color_list(i_condition, :), 'facecolor', lightenColor(color_list(i_condition, :), 0.5));        
    %     end
    %     legend(step_width_left_histograms, experimental_conditions);
    %     
    %     % length histogram
    %     bin_edges = linspace(-0.75, -0.45, 39);
    %     step_length_figure = figure; axes('fontsize', 12); hold on; title('step length', 'fontsize', 16);
    %     step_length_left_histograms = zeros(1, number_of_experimental_conditions);
    %     step_length_right_histograms = zeros(1, number_of_experimental_conditions);
    %     for i_condition = 1 : number_of_experimental_conditions
    %        step_length_left_histograms(i_condition) = histogram(step_length_total(strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, experimental_conditions{i_condition})), 'Normalization', 'pdf', 'displaystyle', 'bar', 'edgecolor', color_list(i_condition, :), 'facecolor', lightenColor(color_list(i_condition, :), 0.5));        
    %        step_length_right_histograms(i_condition) = histogram(step_length_total(strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, experimental_conditions{i_condition})), 'Normalization', 'pdf', 'displaystyle', 'bar', 'edgecolor', color_list(i_condition, :), 'facecolor', lightenColor(color_list(i_condition, :), 0.5));        
    %     end
    %     legend(step_length_left_histograms, experimental_conditions);

        % width box plots
        group_order = cell(2*number_of_experimental_conditions, 1);
        for i_condition = 1 : number_of_experimental_conditions
            left_string = ['LEFT,' experimental_conditions{i_condition}];
            group_order{i_condition} = left_string;
            right_string = ['RIGHT,' experimental_conditions{i_condition}];
            group_order{number_of_experimental_conditions + i_condition} = right_string;
        end

        step_width_figure = figure; axes('fontsize', 12); hold on; title('step width', 'fontsize', 16);
        box_plot_data = boxplot(step_width_total, {condition_stance_foot_list_total condition_experimental_list_total}, 'grouporder', group_order);
        setBoxPlotColors(gca, box_plot_data, group_order, experimental_conditions, color_list);
        if save_figures
            saveas(gcf, 'step_width.eps', 'epsc2')
        end


        step_length_figure = figure; axes('fontsize', 12); hold on; title('step length', 'fontsize', 16);
        box_plot_data = boxplot(step_length_total, {condition_stance_foot_list_total condition_experimental_list_total});
        setBoxPlotColors(gca, box_plot_data, group_order, experimental_conditions, color_list);
        if save_figures
            saveas(gcf, 'step_length.eps', 'epsc2')
        end

        % run quick and dirty ttests
        [h_left_width_onOff, p_left_width_onOff] = ttest2(step_width_total(strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, 'ON')), step_width_total(strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, 'OFF')));
        [h_left_width_onPost, p_left_width_onPost] = ttest2(step_width_total(strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, 'ON')), step_width_total(strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, 'POST')));
        [h_left_width_postOff, p_left_width_postOff] = ttest2(step_width_total(strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, 'POST')), step_width_total(strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, 'OFF')));

        [h_right_width_onOff, p_right_width_onOff] = ttest2(step_width_total(strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, 'ON')), step_width_total(strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, 'OFF')));
        [h_right_width_onPost, p_right_width_onPost] = ttest2(step_width_total(strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, 'ON')), step_width_total(strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, 'POST')));
        [h_right_width_postOff, p_right_width_postOff] = ttest2(step_width_total(strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, 'POST')), step_width_total(strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, 'OFF')));

        [h_left_length_onOff, p_left_length_onOff] = ttest2(step_length_total(strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, 'ON')), step_length_total(strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, 'OFF')));
        [h_left_length_onPost, p_left_length_onPost] = ttest2(step_length_total(strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, 'ON')), step_length_total(strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, 'POST')));
        [h_left_length_postOff, p_left_length_postOff] = ttest2(step_length_total(strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, 'POST')), step_length_total(strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_experimental_list_total, 'OFF')));

        [h_right_length_onOff, p_right_length_onOff] = ttest2(step_length_total(strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, 'ON')), step_length_total(strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, 'OFF')));
        [h_right_length_onPost, p_right_length_onPost] = ttest2(step_length_total(strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, 'ON')), step_length_total(strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, 'POST')));
        [h_right_length_postOff, p_right_length_postOff] = ttest2(step_length_total(strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, 'POST')), step_length_total(strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_experimental_list_total, 'OFF')));



    end
end















