

plot_single         = 0;
plot_absolute       = 1;

plot_rheel_x_pos    = 1;


% variable infor contains the following columns
% variable name | display name | y-label with unit
variable_info = ...
  { ...
    'rcop_x_normalized_all', 'right foot CoP, ml', 'CoP (m)'; ...
  };
%     'rcop_x_response', 'right foot CoP response, ml', 'm'; ...
%     'lcop_x_response'; ...
%     'lcop_x_normalized_all'; ...
%     'lheel_x_pos_normalized_all'; ...
%     'rheel_x_pos_normalized_all'; ...

%% load data
load subjectInfo.mat;
load(makeFileName(date, subject_id, 'resultsConditions'));
load(makeFileName(date, subject_id, 'resultsMarker'));
load(makeFileName(date, subject_id, 'resultsForceplate'));
% load(makeFileName(date, subject_id, 'resultsEmg'));



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

condition_labels = {'stance foot', 'perturbation', 'delay', 'index'};

conditions_control = ...
  {
    'LEFT', 'CONTROL', 'CONTROL', 'CONTROL'; ...
    'RIGHT', 'CONTROL', 'CONTROL', 'CONTROL'; ...
  };

conditions_to_analyze = ...
  {
    'LEFT', 'POSITIVE', '0ms', 'TWO'; ...
    'LEFT', 'POSITIVE', '150ms', 'TWO'; ...
    'LEFT', 'POSITIVE', '450ms', 'TWO'; ...
    'LEFT', 'NEGATIVE', '0ms', 'TWO'; ...
    'LEFT', 'NEGATIVE', '150ms', 'TWO'; ...
    'LEFT', 'NEGATIVE', '450ms', 'TWO'; ...
    'RIGHT', 'POSITIVE', '0ms', 'ONE'; ...
    'RIGHT', 'POSITIVE', '150ms', 'ONE'; ...
    'RIGHT', 'POSITIVE', '450ms', 'ONE'; ...
    'RIGHT', 'NEGATIVE', '0ms', 'ONE'; ...
    'RIGHT', 'NEGATIVE', '150ms', 'ONE'; ...
    'RIGHT', 'NEGATIVE', '450ms', 'ONE'; ...
  };



number_of_conditions_control = size(conditions_control, 1);
number_of_conditions_to_analyze = size(conditions_to_analyze, 1);
    
applicable_control_condition_indices = zeros(number_of_conditions_to_analyze, 1);
for i_condition = 1 : number_of_conditions_to_analyze
    if strcmp(conditions_to_analyze(i_condition, 1), 'LEFT')
        applicable_control_condition_indices(i_condition) = 1;
    elseif strcmp(conditions_to_analyze(i_condition, 1), 'RIGHT')
        applicable_control_condition_indices(i_condition) = 2;
    end
end


% define comparisons
comparison_to_make = 2; % perturbation only
% comparison_to_make = 3; % delay only
comparison_indices = {};
control_conditions_for_comparisons = [];
conditions_already_compared = [];
while length(conditions_already_compared) < number_of_conditions_to_analyze
    % start with the first available condition
    i_condition = 1;
    while ismember(i_condition, conditions_already_compared)
        i_condition = i_condition + 1;
    end
    
    this_comparison = i_condition;
    control_conditions_for_comparisons = [control_conditions_for_comparisons; applicable_control_condition_indices(i_condition)];
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


color_control = [0.3 0.1 1];
colors_comparison = ...
  [ ...
    [1 0.3 0.1] * 0.7; ...
    [0.3 1 0.1] * 0.7; ...
    [1 0.3 0.0] * 1.0; ...
  ]; % should have one row per condition in the comparison

% color_control = [1 1 1] * 0.7;
% color_positive = [82 79 161] * 1/255;
% color_negative = [253 185 19] * 1/255;



%% plot single
if plot_single
    for i_variable = 1 : size(variable_info, 1)
        evalstring = ['trajectories_to_plot = ' variable_info{i_variable, 1} ';'];
        eval(evalstring);
        for i_comparison = 1 : length(comparison_indices);
            figure; axes; hold on;
            % plot control
            condition_indicator = conditions_control_indicators(:, control_conditions_for_comparisons(i_comparison));
            plot ...
              ( ...
                time_normalized, trajectories_to_plot(:, condition_indicator), ...
                'HandleVisibility', 'off', ...
                'color', lightenColor(color_control, 0.5) ...
              );
            control_mean_plot = plot ...
              ( ...
                time_normalized, mean(trajectories_to_plot(:, condition_indicator), 2), ...
                'DisplayName', 'control', ...
                'linewidth', 5, ...
                'color', color_control ...
              );
            % plot stimulus
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
            uistack(control_mean_plot, 'top');
            for i_condition = 1 : length(this_comparison)
                uistack(condition_mean_plots(i_condition), 'top');
            end

            % annotate
            legend('toggle')
            title_string = variable_info{i_variable, 2};
            for i_label = 1 : length(condition_labels);
                if i_label ~= comparison_to_make
                    title_string = [title_string ' - ' conditions_to_analyze{comparison_indices{i_comparison}(1), i_label}];
                end
            end
            title(title_string); set(gca, 'Fontsize', 12)
        end
    end
end






%% plot absolute
if plot_absolute
    for i_variable = 1 : size(variable_info, 1)
        evalstring = ['trajectories_to_plot = ' variable_info{i_variable, 1} ';'];
        eval(evalstring);
        for i_comparison = 1 : length(comparison_indices);
            figure; axes; hold on;
            % plot control
            condition_indicator = conditions_control_indicators(:, control_conditions_for_comparisons(i_comparison));
            current_plots = shadedErrorBar ...
              ( ...
                time_normalized, ...
                mean(trajectories_to_plot(:, condition_indicator), 2), ...
                cinv(trajectories_to_plot(:, condition_indicator), 2), ...
                { ...
                  'color', color_control, ...
                  'linewidth', 3 ...
                }, ...
                1 ...
              );
            legend_handles = current_plots.mainLine;
            legend_data = {'control'};
            % plot stimulus
            this_comparison = comparison_indices{i_comparison};
            condition_mean_plots = zeros(1, length(this_comparison));
            for i_condition = 1 : length(this_comparison)
                condition_indicator = conditions_to_analyze_indicators(:, this_comparison(i_condition));
                current_plots = shadedErrorBar ...
                  ( ...
                    time_normalized, ...
                    mean(trajectories_to_plot(:, condition_indicator), 2), ...
                    cinv(trajectories_to_plot(:, condition_indicator), 2), ...
                    { ...
                      'color', colors_comparison(i_condition, :), ...
                      'linewidth', 3 ...
                    }, ...
                    1 ...
                  );
                legend_handles = [legend_handles, current_plots.mainLine];
                label_string = clarifyConditionString(conditions_to_analyze{comparison_indices{i_comparison}(i_condition), comparison_to_make});
                legend_data = [legend_data, label_string];
            end

            % annotate
            this_legend = legend(legend_handles, legend_data);
            title_string = variable_info{i_variable, 2};
%             title_string = strrep(variables_to_plot{i_variable}, '_', ' ');
%             title_string = strrep(title_string, ' normalized all', '');
            for i_label = 1 : length(condition_labels);
                if i_label ~= comparison_to_make
                    title_string = [title_string ' - ' clarifyConditionString(conditions_to_analyze{comparison_indices{i_comparison}(1), i_label})];
                end
            end
            title(title_string); set(gca, 'Fontsize', 12)
            
            set(gca, 'xlim', [time_normalized(1), time_normalized(end)]);
            xlabel('normalized time (s)');
            ylabel(variable_info{i_variable, 3});
            
            xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
            text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
            text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
        end
    end    
    
    
    
    
    
end











