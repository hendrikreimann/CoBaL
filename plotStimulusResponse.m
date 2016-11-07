

plot_detailed       = 1;
plot_overview       = 0;

show_legend         = 0;
dictate_axes        = 0;

save_figures        = 0;

% variable info contains the following columns
% variable name | display name | y-label with unit | save label
continuous_variable_info = {};
discrete_variable_info = {};

% markers
% the cell array should have the following entries in each line: variable name, variable label, unit, file label for saving, forced axis scale, positive direction label, negative direction label
% continuous_variable_info = [continuous_variable_info; {'lheel_x_pos_normalized_all', 'left heel pos, ml', 'heel pos (m)', 'lheelpos', 0, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'rheel_x_pos_normalized_all', 'right heel pos, ml', 'heel pos (m)', 'rheelpos', 0, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'trunk_angle_ml_normalized_all', 'trunk angle, ml', 'angle (deg)', 'trunkangleml', 0, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'lleg_angle_ml_normalized_all', 'left leg angle, ml', 'angle (deg)', 'llegangleml', 0, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'rleg_angle_ml_normalized_all', 'right leg angle, ml', 'angle (deg)', 'rlegangleml', 0, 'right', 'left'}];

% continuous_variable_info = [continuous_variable_info; {'lheel_x_pos_response', 'left heel pos response, ml', 'heel pos (m)', 'lheelposRsp', 0.05, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'rheel_x_pos_response', 'right heel pos response, ml', 'heel pos (m)', 'rheelposRsp', 0.05, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'trunk_angle_ml_response', 'trunk angle response, ml', 'angle (deg)', 'trunkanglemlRsp', 2, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'lleg_angle_ml_response', 'left leg angle response, ml', 'angle (deg)', 'lleganglemlRsp', 2, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'rleg_angle_ml_response', 'right leg angle response, ml', 'angle (deg)', 'rleganglemlRsp', 2, 'right', 'left'}];

% forceplate
% continuous_variable_info = [continuous_variable_info; {'cop_x_normalized_all', 'total CoP, ml', 'CoP (m)', 'copx', 0, '?', '?'}];
% continuous_variable_info = [continuous_variable_info; {'cop_x_response', 'total CoP response, ml', 'CoP (m)', 'copxRsp', 0.025, '?', '?'}];
% continuous_variable_info = [continuous_variable_info; {'f_x_normalized_all', 'total force, ml', 'f (N)', 'fx', 0, '?', '?'}];
% continuous_variable_info = [continuous_variable_info; {'f_x_response', 'total force response, ml', 'f (N)', 'fxRsp', 0.025, '?', '?'}];
% continuous_variable_info = [continuous_variable_info; {'f_z_normalized_all', 'total vertical force', 'f (N)', 'fz', 0, '?', '?'}];
% continuous_variable_info = [continuous_variable_info; {'f_z_response', 'total force vertical response', 'f (N)', 'fzRsp', 0.025, '?', '?'}];
% continuous_variable_info = [continuous_variable_info; {'m_y_normalized_all', 'total moment, ml', 'm (Nm)', 'my', 0, '?', '?'}];
% continuous_variable_info = [continuous_variable_info; {'m_y_response', 'total moment response, ml', 'm (Nm)', 'myRsp', 0.025, '?', '?'}];

% continuous_variable_info = [continuous_variable_info; {'lcop_x_normalized_all', 'left foot CoP, ml', 'CoP (m)', 'lcopx', 0, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'rcop_x_normalized_all', 'right foot CoP, ml', 'CoP (m)', 'rcopx', 0, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'lcop_x_response', 'left foot CoP response, ml', 'CoP (m)', 'lcopxRsp', 0, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'rcop_x_response', 'right foot CoP response, ml', 'CoP (m)', 'rcopxRsp', 0, 'right', 'left'}];

% EMG
% continuous_variable_info = [continuous_variable_info; {'lglutmed_normalized_all', 'left Gluteus Medius', 'EMG', 'lglutmed'}];
% continuous_variable_info = [continuous_variable_info; {'ltibiant_normalized_all', 'left Tibialis Anterior', 'EMG', 'ltibiant'}];
% continuous_variable_info = [continuous_variable_info; {'lgastroc_normalized_all', 'left Gastrocnemius Medialis', 'EMG', 'lgastroc'}];
% continuous_variable_info = [continuous_variable_info; {'lperolng_normalized_all', 'left Peroneus Longus', 'EMG', 'lperolng'}];
% continuous_variable_info = [continuous_variable_info; {'rglutmed_normalized_all', 'right Gluteus Medius', 'EMG', 'rglutmed'}];
% continuous_variable_info = [continuous_variable_info; {'rtibiant_normalized_all', 'right Tibialis Anterior', 'EMG', 'rtibiant'}];
% continuous_variable_info = [continuous_variable_info; {'rgastroc_normalized_all', 'right Gastrocnemius Medialis', 'EMG', 'rgastroc'}];
% continuous_variable_info = [continuous_variable_info; {'rperolng_normalized_all', 'right Peroneus Longus', 'EMG', 'rperolng'}];

% discrete_variable_info = [discrete_variable_info; {'step_length_all', 'step length'}];
discrete_variable_info = [discrete_variable_info; {'step_width_all', 'step width'}];
% discrete_variable_info = [discrete_variable_info; {'step_times_all', 'step times'}];

%% load data
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

condition_labels = {'stance foot', 'perturbation', 'delay', 'index'};

% for phase-dependent GVS
% conditions_control = ...
%   {
%     'LEFT', 'CONTROL', 'CONTROL', 'CONTROL'; ...
%     'RIGHT', 'CONTROL', 'CONTROL', 'CONTROL'; ...
%   };
% 
% conditions_to_analyze = ...
%   {
%     'LEFT', 'POSITIVE', '0ms', 'TWO'; ...
%     'LEFT', 'POSITIVE', '150ms', 'TWO'; ...
%     'LEFT', 'POSITIVE', '450ms', 'TWO'; ...
%     'LEFT', 'NEGATIVE', '0ms', 'TWO'; ...
%     'LEFT', 'NEGATIVE', '150ms', 'TWO'; ...
%     'LEFT', 'NEGATIVE', '450ms', 'TWO'; ...
%     'RIGHT', 'POSITIVE', '0ms', 'ONE'; ...
%     'RIGHT', 'POSITIVE', '150ms', 'ONE'; ...
%     'RIGHT', 'POSITIVE', '450ms', 'ONE'; ...
%     'RIGHT', 'NEGATIVE', '0ms', 'ONE'; ...
%     'RIGHT', 'NEGATIVE', '150ms', 'ONE'; ...
%     'RIGHT', 'NEGATIVE', '450ms', 'ONE'; ...
%   };

% for vision
% conditions_control = ...
%   {
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL'; ...
%     'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'CONTROL'; ...
%   };

% all
% conditions_to_analyze = ...
%   {
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'ONE'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'ONE'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'ONE'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'ONE'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'TWO'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'TWO'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'TWO'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'TWO'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'THREE'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'THREE'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'THREE'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'THREE'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'FOUR'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'FOUR'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'FOUR'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'FOUR'; ...
%   };

% ACHTUNG! do not use this for now. The list of conditions to analyze should be exactly the same as defined in the
% previous processing step. Rather, save this list, then replace this with conditions_to_plot or something.

% first step left stance
% conditions_to_analyze = ...
%   {
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'ONE'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'ONE'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'TWO'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'TWO'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'THREE'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'THREE'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'FOUR'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'FOUR'; ...
%   };

% % first step right stance
% conditions_to_analyze = ...
%   {
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'ONE'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'ONE'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'TWO'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'TWO'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'THREE'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'THREE'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'FOUR'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'FOUR'; ...
%   };

number_of_conditions_control = size(conditions_control, 1);
number_of_conditions_to_analyze = size(conditions_to_analyze, 1);
    
applicable_control_condition_indices = zeros(number_of_conditions_to_analyze, 1);
for i_condition = 1 : number_of_conditions_to_analyze
    if strcmp(conditions_to_analyze(i_condition, 1), 'STANCE_LEFT')
        applicable_control_condition_indices(i_condition) = 1;
    elseif strcmp(conditions_to_analyze(i_condition, 1), 'STANCE_RIGHT')
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

%% do plots
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



%% plot detailed
if plot_detailed
    % discrete variables
    for i_variable = 1 : size(discrete_variable_info, 1)
        evalstring = ['points_to_plot = ' discrete_variable_info{i_variable, 1} ';'];
        eval(evalstring);
        for i_comparison = 1 : length(comparison_indices);
            figure; axes; hold on;
            
            % ... control plot is missing here
           u  
            this_comparison = comparison_indices{i_comparison};
            for i_condition = 1 : length(this_comparison)
                condition_indicator = conditions_to_analyze_indicators(:, this_comparison(i_condition));
                label_string = conditions_to_analyze{comparison_indices{i_comparison}(i_condition), comparison_to_make};
                histogram ...
                  ( ...
                    points_to_plot(condition_indicator), ...
                    'edgecolor', colors_comparison(i_condition, :), ...
                    'facecolor', lightenColor(colors_comparison(i_condition, :), 0.5), ...
                    'DisplayName', label_string ...
                  )
            end

            % annotate
            if show_legend
                legend('toggle')
            end
            title_string = discrete_variable_info{i_variable, 2};
            for i_label = 1 : length(condition_labels);
                if i_label ~= comparison_to_make
                    title_string = [title_string ' - ' conditions_to_analyze{comparison_indices{i_comparison}(1), i_label}];
                end
            end
            title(title_string); set(gca, 'Fontsize', 12)
        end
    end    
    
    
    
    
    for i_variable = 1 : size(continuous_variable_info, 1)
        evalstring = ['trajectories_to_plot = ' continuous_variable_info{i_variable, 1} ';'];
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


%% plot overview
if plot_overview
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
            
            group_order = conditions_to_analyze(comparison_indices{i_comparison}, comparison_to_make);
            
            figure; axes; hold on;
            box_plot_data = boxplot(points_to_plot_pruned, condition_labels_for_boxplot_pruned, 'grouporder', group_order);
            
            
            

            

            experimental_conditions = group_order;
            
            setBoxPlotColors(gca, box_plot_data, group_order, experimental_conditions, colors_comparison);
            

            % annotate
            title_string = discrete_variable_info{i_variable, 2};
            filename_string = discrete_variable_info{i_variable, 2};
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
%                 label_string = conditionStringToTitle(conditions_to_analyze{comparison_indices{i_comparison}(i_condition), comparison_to_make});
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
            text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), [continuous_variable_info{i_variable, 6} ' $\rightarrow$'] , 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right', 'interpreter', 'LaTeX')
            text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), ['$\leftarrow$ ' continuous_variable_info{i_variable, 7}], 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left', 'interpreter', 'LaTeX')
            
            % save
            if save_figures
                filename = ['../figures/' filename_string '.eps'];
                saveas(gcf, filename, 'epsc2')
                
            end
        end
    end    
    
    
    
    
    
end











