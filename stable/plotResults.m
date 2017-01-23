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
        header_string = fgetl(fid);
        unit_string = fgetl(fid);
        data_raw = textscan(fid, format);
        fclose(fid);
        
        % transform to cell
        data_lines = data_raw{1};
        data_cell = {};
        for i_line = 1 : length(data_lines)
            line_split = strsplit(data_lines{i_line}, ',');
            try
                data_cell = [data_cell; line_split];
            end
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
        condition_stance_foot_list_all = [condition_stance_foot_list_all; condition_stance_foot_list_subject];
        condition_perturbation_list_all = [condition_perturbation_list_all; condition_perturbation_list_subject];
        condition_delay_list_all = [condition_delay_list_all; condition_delay_list_subject];
        condition_index_list_all = [condition_index_list_all; condition_index_list_subject];
        condition_experimental_list_all = [condition_experimental_list_all; condition_experimental_list_subject];
        for i_variable = 1 : number_of_variables_to_plot
            % load and extract data
            this_variable_name = study_settings.variables_to_plot{i_variable, 1};
            index_in_saved_data = find(strcmp(variable_names_subject, this_variable_name), 1, 'first');
            this_variable_data = variable_data_subject{index_in_saved_data};
            
            % store
            variable_data_all{i_variable} = [variable_data_all{i_variable} this_variable_data];
        end


    end
    
    %% create figures and determine abscissae for each comparison
    comparison_variable_to_axes_index_map = zeros(number_of_comparisons, 1);
    comparison_to_abscissae_map = cell(number_of_comparisons, 1);
    
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
                comparison_to_abscissae_map{i_comparison} = 1 : 100;
                
                % set axes properties
                if dictate_axes
    %                 set(gca, 'xlim', [time_normalized(1), time_normalized(end)]);
                    set(gca, 'ylim', [str2double(study_settings.variables_to_plot{i_variable, 5}), str2double(study_settings.variables_to_plot{i_variable, 6})]);
                end

                % determine title
                title_string = study_settings.variables_to_plot{i_variable, 2};
                for i_label = 1 : length(study_settings.condition_labels);
                    if i_label ~= study_settings.comparison_to_make
                        title_string = [title_string ' - ' study_settings.conditions_to_plot{comparison_indices{i_comparison}(1), i_label}];
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
                        comparison_to_abscissae_map{this_episode(i_comparison)} = 1 : 100;
                    elseif strcmp(condition_identifier{4}, 'TWO')
                        comparison_to_abscissae_map{this_episode(i_comparison)} = 101 : 200;
                    elseif strcmp(condition_identifier{4}, 'THREE')
                        comparison_to_abscissae_map{this_episode(i_comparison)} = 201 : 300;
                    elseif strcmp(condition_identifier{4}, 'FOUR')
                        comparison_to_abscissae_map{this_episode(i_comparison)} = 301 : 400;
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
                        title_string = [title_string ' - ' strrep(study_settings.conditions_to_plot{comparison_indices{i_comparison}(1), i_label}, '_', ' ')];
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
            target_abscissa = comparison_to_abscissae_map{i_comparison};
            
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
                    top_level_plots = [top_level_plots control_mean_plot];
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
                    top_level_plots = [top_level_plots plot_handles.mainLine];
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
                    top_level_plots = [top_level_plots condition_mean_plot];
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
                    top_level_plots = [top_level_plots plot_handles.mainLine];
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
    number_of_conditions_control = size(study_settings.conditions_control, 1);
    number_of_conditions_to_plot = size(study_settings.conditions_to_plot, 1);
    use_control = ~isempty(study_settings.conditions_control);
    comparison_indices = {};
    control_conditions_for_comparisons = [];
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
                    this_comparison = [this_comparison, j_condition];
                end
            end
        end
        comparison_indices = [comparison_indices; this_comparison];
        conditions_already_compared = [conditions_already_compared this_comparison];
    end    

end

function episode_indices = determineEpisodes(study_settings, comparison_indices)
    condition_column_index = find(strcmp(study_settings.condition_labels, 'index'));
    condition_column_stancefoot = find(strcmp(study_settings.condition_labels, 'stance foot'));
    episode_first_stretch_indices = find(strcmp(study_settings.conditions_to_plot(:, 4), 'ONE'));
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
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison];
                        elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'THREE') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison];
                        elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'FOUR') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison];
                        end
                    elseif strcmp(example_condition_in_base_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                        if strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'TWO') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison];
                        elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'THREE') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison];
                        elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'FOUR') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison];
                        end
                    end
                end
            end
        end
        episode_indices = [episode_indices; comparison_indices_in_this_episode];
        comparisons_already_used = [comparisons_already_used comparison_indices_in_this_episode];

    end   
end
