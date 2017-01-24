%     This file is part of the CoBaL code base
%     Copyright (C) 2017 Hendrik Reimann <hendrikreimann@gmail.com>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% plot results

% input
% ... results.mat for each subject

function plotResults(varargin)
    %% parse input
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'subjects', [])
    addParameter(parser, 'dictate_axes', false)
    addParameter(parser, 'show_legend', false)
    parse(parser, varargin{:})
    subjects = parser.Results.subjects;
    dictate_axes = parser.Results.dictate_axes;
    show_legend = parser.Results.show_legend;

    % load study settings
    if exist('studySettings.txt', 'file')
        % study settings not found here, so we're in a subject folder, load from level above
        study_settings = loadSettingsFile('studySettings.txt');
    else
        study_settings = loadSettingsFile(['..' filesep 'studySettings.txt']);
    end

    %% determine subjects
    if isempty(subjects)
        % nothing passed, so load from file
        subject_data_file = ['..' filesep 'subjects.csv'];
        format = '%s';
        fid = fopen(subject_data_file);
        fgetl(fid);
        fgetl(fid);
        data_raw = textscan(fid, format);
        fclose(fid);
        
        % transform to cell
        data_lines = data_raw{1};
        data_cell = {};
        for i_line = 1 : length(data_lines)
            line_split = strsplit(data_lines{i_line}, ',');
            data_cell = [data_cell; line_split]; %#ok<AGROW>
        end
        
        subjects = data_cell(:, 1);
    end
    if ischar(subjects)
        % single string (i.e. char array) was passed, make a cell out of this
        subjects = {subjects};
    end
    
    comparison_indices = determineComparisons(study_settings);
    number_of_comparisons = length(comparison_indices);
    episode_indices = determineEpisodes(study_settings, comparison_indices);
    number_of_episodes = length(episode_indices);
    
    %% collect data from all subjects
    number_of_variables_to_plot = size(study_settings.variables_to_plot, 1);
    condition_stance_foot_list_all = {};
    condition_perturbation_list_all = {};
    condition_delay_list_all = {};
    condition_index_list_all = {};
    condition_experimental_list_all = {};
    variable_data_all = cell(number_of_variables_to_plot, 1);

    for i_subject = 1 : length(subjects)
        % load subject data
        path = strsplit(pwd, filesep);
        if strcmp(path(end), subjects{i_subject})
            % we're already in the subject folder, so just load from here
            data_path = '';
        else
            % we're in the study root, so load from subject folder
            data_path = [subjects{i_subject} filesep];
        end
        load([data_path 'subjectInfo.mat'], 'date', 'subject_id');
        load([data_path 'analysis' filesep date '_' subject_id '_results.mat']);

        % append data from this subject to containers for all subjects
        condition_stance_foot_list_all = [condition_stance_foot_list_all; condition_stance_foot_list_subject]; %#ok<AGROW>
        condition_perturbation_list_all = [condition_perturbation_list_all; condition_perturbation_list_subject]; %#ok<AGROW>
        condition_delay_list_all = [condition_delay_list_all; condition_delay_list_subject]; %#ok<AGROW>
        condition_index_list_all = [condition_index_list_all; condition_index_list_subject]; %#ok<AGROW>
        condition_experimental_list_all = [condition_experimental_list_all; condition_experimental_list_subject]; %#ok<AGROW>
        for i_variable = 1 : number_of_variables_to_plot
            % load and extract data
            this_variable_name = study_settings.variables_to_plot{i_variable, 1};
            index_in_saved_data = find(strcmp(variable_names_subject, this_variable_name), 1, 'first');
            this_variable_data = variable_data_subject{index_in_saved_data}; %#ok<USENS>
            
            % store
            variable_data_all{i_variable} = [variable_data_all{i_variable} this_variable_data];
        end


    end
    
    %% create figures and determine abscissae for each comparison
    comparison_variable_to_axes_index_map = zeros(number_of_comparisons, 1);
    abscissae_cell = cell(number_of_comparisons, number_of_variables_to_plot);
    
    if strcmp(study_settings.plot_mode, 'detailed') || strcmp(study_settings.plot_mode, 'overview')
        % make one figure per comparison and variable
        axes_handles = zeros(number_of_comparisons, number_of_variables_to_plot);
        for i_variable = 1 : number_of_variables_to_plot
            for i_comparison = 1 : number_of_comparisons
                % make figure and axes
                figure; new_axes = axes; hold on;
                
                % store handles and determine abscissa data
                axes_handles(i_comparison, i_variable) = new_axes;
                comparison_variable_to_axes_index_map(i_comparison) = i_comparison;
                if isDiscreteVariable(i_variable, variable_data_all)
                    % abscissae gives the bin edges here
                    if dictate_axes
                        lower_bound = str2double(study_settings.variables_to_plot{i_variable, 5});
                        upper_bound = str2double(study_settings.variables_to_plot{i_variable, 6});
                    else
                        lower_bound = min(variable_data_all{i_variable});
                        upper_bound = max(variable_data_all{i_variable});
                    end
                    abscissae_cell{i_comparison, i_variable} = linspace(lower_bound, upper_bound, study_settings.number_of_bins_in_histogram);
                end
                if isContinuousVariable(i_variable, variable_data_all)
                    abscissae_cell{i_comparison, i_variable} = 1 : 100;
                end
                
                % set axes properties
                if dictate_axes
    %                 set(gca, 'xlim', [time_normalized(1), time_normalized(end)]);
                    set(gca, 'ylim', [str2double(study_settings.variables_to_plot{i_variable, 5}), str2double(study_settings.variables_to_plot{i_variable, 6})]);
                    if strcmp(study_settings.plot_mode, 'overview')
                        set(gca, 'xlim', [-0.5 length(comparison_indices{i_comparison})+0.5]);
                    end
                end

                % determine title
                title_string = study_settings.variables_to_plot{i_variable, 2};
                for i_label = 1 : length(study_settings.condition_labels);
                    if i_label ~= study_settings.comparison_to_make
                        title_string = [title_string ' - ' strrep(study_settings.conditions_to_plot{comparison_indices{i_comparison}(1), i_label}, '_', ' ')]; %#ok<AGROW>
                    end
                end
                title(title_string); set(gca, 'Fontsize', 12)
            end
        end
    end
    if strcmp(study_settings.plot_mode, 'episodes')
        % make one figure per episode and variable
        axes_handles = zeros(number_of_episodes, number_of_variables_to_plot);
        for i_variable = 1 : number_of_variables_to_plot
            for i_episode = 1 : number_of_episodes
                % make figure and axes and store handles
                figure; new_axes = axes; hold on;
                axes_handles(i_episode, i_variable) = new_axes;
                this_episode = episode_indices{i_episode};

                % store handles and determine abscissa data for all comparisons in this episode
                for i_comparison = 1 : length(this_episode)
                    comparison_variable_to_axes_index_map(this_episode(i_comparison)) = i_episode;

                    % determine which step this is
                    this_comparison = this_episode(i_comparison);
                    conditions_in_this_comparison = comparison_indices{this_comparison};
                    example_condition = conditions_in_this_comparison(1);
                    condition_identifier = study_settings.conditions_to_plot(example_condition, :);
                    if strcmp(condition_identifier{4}, 'ONE')
                        abscissae_cell{this_episode(i_comparison), i_variable} = 1 : 100;
                    elseif strcmp(condition_identifier{4}, 'TWO')
                        abscissae_cell{this_episode(i_comparison), i_variable} = 101 : 200;
                    elseif strcmp(condition_identifier{4}, 'THREE')
                        abscissae_cell{this_episode(i_comparison), i_variable} = 201 : 300;
                    elseif strcmp(condition_identifier{4}, 'FOUR')
                        abscissae_cell{this_episode(i_comparison), i_variable} = 301 : 400;
                    end
                end
                
                % set axes properties
                if dictate_axes
    %                 set(gca, 'xlim', [time_normalized(1), time_normalized(end)]);
                    set(gca, 'ylim', [str2double(study_settings.variables_to_plot{i_variable, 5}), str2double(study_settings.variables_to_plot{i_variable, 6})]);
                end

                % determine title
                title_string = study_settings.variables_to_plot{i_variable, 2};
                for i_label = 1 : length(study_settings.condition_labels);
                    if (i_label ~= study_settings.comparison_to_make) && (i_label ~= 4)
                        title_string = [title_string ' - ' strrep(study_settings.conditions_to_plot{comparison_indices{i_comparison}(1), i_label}, '_', ' ')]; %#ok<AGROW>
                    end
                end
                title(title_string); set(gca, 'Fontsize', 12)
                
            end
        end
        
        
    end
    
    %% plot data
    for i_variable = 1 : number_of_variables_to_plot
        data_to_plot = variable_data_all{i_variable, 1};
        for i_comparison = 1 : length(comparison_indices);
            % find correct condition indicator for control
            conditions_this_comparison = comparison_indices{i_comparison};
            top_level_plots = [];
            target_axes_handle = axes_handles(comparison_variable_to_axes_index_map(i_comparison), i_variable);
            target_abscissa = abscissae_cell{i_comparison, i_variable};
            
            % plot control
            if ~isempty(study_settings.conditions_control)
                % determine which control condition applies here
                representant_condition_index = conditions_this_comparison(1);
                this_condition_stance_foot = study_settings.conditions_to_plot(representant_condition_index, 1);
                applicable_control_condition_index = find(strcmp(study_settings.conditions_control, this_condition_stance_foot), 1, 'first');
                applicable_control_condition_labels = study_settings.conditions_control(applicable_control_condition_index, :);
                % NOTE: so far the applicable control conditions is determined only by stance foot. A general solution 
                % is not needed at this time, but can be added later if desired
                
                % extract data for control condition
                stance_foot_indicator = strcmp(condition_stance_foot_list_all, applicable_control_condition_labels{1});
                perturbation_indicator = strcmp(condition_perturbation_list_all, applicable_control_condition_labels{2});
                delay_indicator = strcmp(condition_delay_list_all, applicable_control_condition_labels{3});
                index_indicator = strcmp(condition_index_list_all, applicable_control_condition_labels{4});
                experimental_indicator = strcmp(condition_experimental_list_all, applicable_control_condition_labels{5});
                this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator;
                data_to_plot_this_condition = data_to_plot(:, this_condition_indicator);
                
                if isDiscreteVariable(i_variable, variable_data_all)
                    if strcmp(study_settings.plot_mode, 'detailed')
                        histogram ...
                          ( ...
                            target_axes_handle, ...
                            data_to_plot_this_condition, ...
                            'binEdges', target_abscissa, ...
                            'edgecolor', study_settings.color_control, ...
                            'facecolor', lightenColor(study_settings.color_control, 0.5) ...
                          );
                    end
                end
                if strcmp(study_settings.plot_mode, 'overview')
%                     box_plot_data = boxplot(target_axes_handle, data_to_plot_this_condition, 'widths', 0.8);
                    
%                     fake_data = [randn(100, 1); 3; 3.3];
%                     box_plot_data = boxplot(target_axes_handle, fake_data, 'widths', 0.8);
                    
%                     setBoxAbscissa(box_plot_data, 0);
%                     setBoxColors(box_plot_data, study_settings.color_control);
                    
                    singleBoxPlot(target_axes_handle, 0, data_to_plot_this_condition, study_settings.color_control)
                end
                if isContinuousVariable(i_variable, variable_data_all)
                    if strcmp(study_settings.plot_mode, 'detailed')
                        % individual trajectories
                        plot ...
                          ( ...
                            target_axes_handle, ...
                            target_abscissa, ...
                            data_to_plot_this_condition, ...
                            'HandleVisibility', 'off', ...
                            'color', lightenColor(study_settings.color_control, 0.5) ...
                          );
                        % condition average
                        control_mean_plot = plot ...
                          ( ...
                            target_axes_handle, ...
                            target_abscissa, ...
                            mean(data_to_plot_this_condition, 2), ...
                            'DisplayName', 'CONTROL', ...
                            'linewidth', 5, ...
                            'color', study_settings.color_control ...
                          );
                        top_level_plots = [top_level_plots control_mean_plot]; %#ok<AGROW>
                    end
                    if strcmp(study_settings.plot_mode, 'overview') || strcmp(study_settings.plot_mode, 'episodes')
                        plot_handles = shadedErrorBar ...
                          ( ...
                            target_abscissa, ...
                            mean(data_to_plot_this_condition, 2), ...
                            cinv(data_to_plot_this_condition, 2), ...
                            { ...
                              'color', study_settings.color_control, ...
                              'linewidth', 6 ...
                            }, ...
                            1, ...
                            target_axes_handle ...
                          );
                        top_level_plots = [top_level_plots plot_handles.mainLine]; %#ok<AGROW>
                    end
                end
            end
            
            % plot stimulus
            for i_condition = 1 : length(conditions_this_comparison)
                % find correct condition indicator
                condition_identifier = study_settings.conditions_to_plot(conditions_this_comparison(i_condition), :);
                stance_foot_indicator = strcmp(condition_stance_foot_list_all, condition_identifier{1});
                perturbation_indicator = strcmp(condition_perturbation_list_all, condition_identifier{2});
                delay_indicator = strcmp(condition_delay_list_all, condition_identifier{3});
                index_indicator = strcmp(condition_index_list_all, condition_identifier{4});
                experimental_indicator = strcmp(condition_experimental_list_all, condition_identifier{5});
                this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator;
                data_to_plot_this_condition = data_to_plot(:, this_condition_indicator);
                if isDiscreteVariable(i_variable, variable_data_all)
                    if strcmp(study_settings.plot_mode, 'detailed')
                        histogram ...
                          ( ...
                            target_axes_handle, ...
                            data_to_plot_this_condition, ...
                            target_abscissa, ...
                            'edgecolor', study_settings.colors_comparison(i_condition, :), ...
                            'facecolor', lightenColor(study_settings.colors_comparison(i_condition, :), 0.5) ...
                          );
                    end
                    if strcmp(study_settings.plot_mode, 'overview')
%                         box_plot_data = boxplot(target_axes_handle, data_to_plot_this_condition, 'widths', 0.8);
%                         setBoxAbscissa(box_plot_data, i_condition);
%                         setBoxColors(box_plot_data, study_settings.colors_comparison(i_condition, :));
                        singleBoxPlot(target_axes_handle, i_condition, data_to_plot_this_condition, study_settings.colors_comparison(i_condition, :))
                    end
                end
                if isContinuousVariable(i_variable, variable_data_all)
                    if strcmp(study_settings.plot_mode, 'detailed')
                        plot ...
                          ( ...
                            target_axes_handle, ...
                            target_abscissa, ...
                            data_to_plot_this_condition, ...
                            'HandleVisibility', 'off', ...
                            'color', lightenColor(study_settings.colors_comparison(i_condition, :), 0.5) ...
                          );
                        label_string = strrep(study_settings.conditions_to_plot{comparison_indices{i_comparison}(i_condition), study_settings.comparison_to_make}, '_', ' ');
                        condition_mean_plot = plot ...
                          ( ...
                            target_axes_handle, ...
                            target_abscissa, ...
                            mean(data_to_plot_this_condition, 2), ...
                            'DisplayName', label_string, ...
                            'linewidth', 5, ...
                            'color', study_settings.colors_comparison(i_condition, :) ...
                          );
                        top_level_plots = [top_level_plots condition_mean_plot]; %#ok<AGROW>
                    end
                    if strcmp(study_settings.plot_mode, 'overview') || strcmp(study_settings.plot_mode, 'episodes')
                        plot_handles = shadedErrorBar ...
                          ( ...
                            target_abscissa, ...
                            mean(data_to_plot_this_condition, 2), ...
                            cinv(data_to_plot_this_condition, 2), ...
                            { ...
                              'color', study_settings.colors_comparison(i_condition, :), ...
                              'linewidth', 6 ...
                            }, ...
                            1, ...
                            target_axes_handle ...
                          );
                        top_level_plots = [top_level_plots plot_handles.mainLine]; %#ok<AGROW>
                    end
                end
            end

            % reorder to bring mean plots on top
            for i_plot = 1 : length(top_level_plots)
                uistack(top_level_plots(i_plot), 'top');
            end
            
            % toggle legend
            if show_legend
                legend(target_axes_handle, 'toggle')
            end
        end
    end
    
    
    
end

%% helper functions

function comparison_indices = determineComparisons(study_settings)
    
    % determine comparisons
    number_of_conditions_to_plot = size(study_settings.conditions_to_plot, 1);
    comparison_indices = {};
    conditions_already_compared = [];
    % here, we go through all conditions_to_plot and group up those that go into one comparison, i.e. one figure
    while length(conditions_already_compared) < number_of_conditions_to_plot
        % start with the first available condition
        i_condition = 1;
        while ismember(i_condition, conditions_already_compared)
            i_condition = i_condition + 1;
        end

        this_comparison = i_condition; % this is the first condition in this episode, more will be added
        % search for conditions that differ from this one only in the one we're comparing
        for j_condition = 1 : number_of_conditions_to_plot
            if i_condition ~= j_condition
                % check which conditions labels agree between these two conditions
                comparison_table = zeros(1, length(study_settings.condition_labels)); % this is a table indicating equality between the two conditions in questions
                for i_label = 1 : length(study_settings.condition_labels)
                    comparison_table(i_label) = strcmp(study_settings.conditions_to_plot{i_condition, i_label}, study_settings.conditions_to_plot{j_condition, i_label});
                end

                % look at the relevant entries of the comparison table
                comparison_table_relevant = comparison_table;
                comparison_table_relevant(study_settings.comparison_to_make) = [];
                if all(comparison_table_relevant)
                    this_comparison = [this_comparison, j_condition]; %#ok<AGROW>
                end
            end
        end
        comparison_indices = [comparison_indices; this_comparison]; %#ok<AGROW>
        conditions_already_compared = [conditions_already_compared this_comparison]; %#ok<AGROW>
    end    

end

function episode_indices = determineEpisodes(study_settings, comparison_indices)
    condition_column_index = find(strcmp(study_settings.condition_labels, 'index'));
    condition_column_stancefoot = find(strcmp(study_settings.condition_labels, 'stance foot'));
%     episode_first_stretch_indices = find(strcmp(study_settings.conditions_to_plot(:, 4), 'ONE'));
    episode_indices = {};
    comparisons_already_used = [];
    number_of_comparisons = length(comparison_indices);
    while length(comparisons_already_used) < number_of_comparisons
        % start with the first available comparison
        i_comparison = 1;
        while ismember(i_comparison, comparisons_already_used)
            i_comparison = i_comparison + 1;
        end
        comparison_indices_in_this_episode = i_comparison; % this is the first comparison in this episode, more will be added

        % search for comparisons that differ from this one in only the step number
        base_comparison = comparison_indices{i_comparison};
        example_condition_in_base_comparison = base_comparison(1);
        example_condition_in_base_comparison_labels = study_settings.conditions_to_plot(example_condition_in_base_comparison, :);
        for j_comparison = 1 : number_of_comparisons
            if i_comparison ~= j_comparison
                this_comparison = comparison_indices{j_comparison};
                example_condition_in_this_comparison = this_comparison(1);
                example_condition_in_this_comparison_labels = study_settings.conditions_to_plot(example_condition_in_this_comparison, :);
                % check which conditions labels agree between these two conditions
                comparison_table = zeros(1, length(study_settings.condition_labels)); % this is a table indicating equality between the two conditions in questions
                for i_label = 1 : length(study_settings.condition_labels)
                    comparison_table(i_label) = strcmp(example_condition_in_base_comparison_labels{i_label}, example_condition_in_this_comparison_labels{i_label});
                end

                % look at the relevant entries of the comparison table
                comparison_table_relevant = comparison_table;
                comparison_table_relevant([condition_column_stancefoot condition_column_index study_settings.comparison_to_make]) = [];
                if all(comparison_table_relevant)
                    % check if the stance foot is alternating
                    if strcmp(example_condition_in_base_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                        if strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'TWO') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'THREE') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'FOUR') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        end
                    elseif strcmp(example_condition_in_base_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                        if strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'TWO') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'THREE') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'FOUR') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
        episode_indices = [episode_indices; comparison_indices_in_this_episode]; %#ok<AGROW>
        comparisons_already_used = [comparisons_already_used comparison_indices_in_this_episode]; %#ok<AGROW>

    end   
end

function discrete = isDiscreteVariable(variable_index, variable_data)
    discrete = false;
    if size(variable_data{variable_index}, 1) == 1
        discrete = true;
    end
end

function continuous = isContinuousVariable(variable_index, variable_data)
    continuous = false;
    if size(variable_data{variable_index}, 1) == 100
        continuous = true;
    end
end

function singleBoxPlot(target_axes_handle, abscissa, data, color)
    % set some parameters, these should be name-value pair arguments later
    width = 0.8;
    
    % extract data
    data_median = median(data);
    data_quartile_1 = prctile(data, 25);
    data_quartile_3 = prctile(data, 75);
    
    data_percentile_9 = prctile(data, 9);
    data_percentile_91 = prctile(data, 91);
    data_mean = mean(data);
    data_std = std(data);
    
    data_iqr = iqr(data);
    data_upper_inner_fence = data_quartile_3 + 1.5*data_iqr;
    data_lower_inner_fence = data_quartile_1 - 1.5*data_iqr;
    data_upper_adjacent = max(data(data<data_upper_inner_fence));
    data_lower_adjacent = min(data(data>data_lower_inner_fence));
    outliers = data(data>data_upper_inner_fence | data<data_lower_inner_fence);
    
    % plot
    box_x_data = [abscissa-width/2 abscissa+width/2 abscissa+width/2 abscissa-width/2 abscissa-width/2];
    box_y_data = [data_quartile_1 data_quartile_1 data_quartile_3 data_quartile_3 data_quartile_1];
    patch ...
      ( ...
        box_x_data, ...
        box_y_data, ...
        color, ...
        'parent', target_axes_handle, ...
        'EdgeColor', 'none' ...
      );
    plot(target_axes_handle, abscissa + width*[-0.5 0.5], [data_median data_median], 'color', 'k'); % median
    plot(target_axes_handle, [abscissa abscissa], [data_quartile_3 data_upper_adjacent], 'k--'); % upper range
    plot(target_axes_handle, [abscissa abscissa], [data_lower_adjacent data_quartile_1], 'k--'); % lower range
    plot(target_axes_handle, abscissa+width*[-0.25 0.25], [data_lower_adjacent data_lower_adjacent], 'k-'); % max
    plot(target_axes_handle, abscissa+width*[-0.25 0.25], [data_upper_adjacent data_upper_adjacent], 'k-'); % min
    plot(target_axes_handle, abscissa * ones(size(outliers)), outliers, '+', 'color', [1; 1; 1] * 0.7);
end

% function setBoxColors(box_plot_data, color)
%     % median line
%     set(box_plot_data(6), 'color', 'k')
% 
%     % box
%     patch_handle = ...
%         patch ...
%           ( ...
%             get(box_plot_data(5), 'XData'), ...
%             get(box_plot_data(5), 'YData'), ...
%             color, ...
%             'parent', get(box_plot_data(1), 'parent'), ...
%             'EdgeColor', 'none' ...
%           ); 
%     uistack(patch_handle, 'bottom')
%     
%     % outliers
%     set(box_plot_data(7), 'MarkerEdgeColor', [1; 1; 1] * 0.7);
% 
%     % remove box edges
%     set(box_plot_data(5), 'color', 'none');
% end

% function setBoxAbscissa(box_plot_data, abscissa)
%     for i_plot = 1 : length(box_plot_data)
%         set(box_plot_data(i_plot), 'xdata', get(box_plot_data(i_plot), 'xdata') - 1 + abscissa);
%     end
% 
% end


