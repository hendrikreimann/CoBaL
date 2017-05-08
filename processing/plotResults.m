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
    addParameter(parser, 'save', false)
    addParameter(parser, 'format', 'epsc')
    parse(parser, varargin{:})
    subjects = parser.Results.subjects;
    dictate_axes = parser.Results.dictate_axes;
    show_legend = parser.Results.show_legend;

    % load study settings
    study_settings_file = '';
    if exist('studySettings.txt', 'file')
        study_settings_file = 'studySettings.txt';
    end    
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    

    %% determine subjects and data folders
    data_folder_list = determineDataStructure(subjects);
    [comparison_indices, conditions_per_comparison_max] = determineComparisons(study_settings);
    number_of_comparisons = length(comparison_indices);
    episode_indices = determineEpisodes(study_settings, comparison_indices);
    number_of_episodes = length(episode_indices);
    
    %% collect data from all data folders
    variables_to_plot = study_settings.get('variables_to_plot');
    number_of_variables_to_plot = size(variables_to_plot, 1);
    condition_stance_foot_list_all = {};
    condition_perturbation_list_all = {};
    condition_delay_list_all = {};
    condition_index_list_all = {};
    condition_experimental_list_all = {};
    condition_stimulus_list_all = {};
    condition_day_list_all = {};
    origin_trial_list_all = [];
    origin_start_time_list_all = [];
    origin_end_time_list_all = [];
    variable_data_all = cell(number_of_variables_to_plot, 1);
    if study_settings.get('plot_response')
        response_data_all = cell(number_of_variables_to_plot, 1);
    end
    step_time_data = [];
    
    for i_folder = 1 : length(data_folder_list)
        % load data
        data_path = data_folder_list{i_folder};
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
        load([data_path filesep 'analysis' filesep date '_' subject_id '_results.mat']);

        % append data from this subject to containers for all subjects
        condition_stance_foot_list_all = [condition_stance_foot_list_all; condition_stance_foot_list_session]; %#ok<AGROW>
        condition_perturbation_list_all = [condition_perturbation_list_all; condition_perturbation_list_session]; %#ok<AGROW>
        condition_delay_list_all = [condition_delay_list_all; condition_delay_list_session]; %#ok<AGROW>
        condition_index_list_all = [condition_index_list_all; condition_index_list_session]; %#ok<AGROW>
        condition_experimental_list_all = [condition_experimental_list_all; condition_experimental_list_session]; %#ok<AGROW>
        condition_stimulus_list_all = [condition_stimulus_list_all; condition_stimulus_list_session]; %#ok<AGROW>
        condition_day_list_all = [condition_day_list_all; condition_day_list_session]; %#ok<AGROW>
        origin_trial_list_all = [origin_trial_list_all; origin_trial_list_session]; %#ok<AGROW>
        origin_start_time_list_all = [origin_start_time_list_all; origin_start_time_list_session]; %#ok<AGROW>
        origin_end_time_list_all = [origin_end_time_list_all; origin_end_time_list_session]; %#ok<AGROW>
        for i_variable = 1 : number_of_variables_to_plot
            % load and extract data
            this_variable_name = variables_to_plot{i_variable, 1};
            index_in_saved_data = find(strcmp(variable_names_session, this_variable_name), 1, 'first');
            this_variable_data = variable_data_session{index_in_saved_data}; %#ok<USENS>
            if study_settings.get('plot_response')
                this_response_data = response_data_session{index_in_saved_data}; %#ok<USENS>
            end
            
            % store
            variable_data_all{i_variable} = [variable_data_all{i_variable} this_variable_data];
            if study_settings.get('plot_response')
                response_data_all{i_variable} = [response_data_all{i_variable} this_response_data];
            end
        end
        if strcmp(study_settings.get('time_plot_style'), 'scaled_to_comparison_mean') || strcmp(study_settings.get('time_plot_style'), 'scaled_to_condition_mean')
            index_in_saved_data = find(strcmp(variable_names_session, 'step_time'), 1, 'first');
            this_step_time_data = variable_data_session{index_in_saved_data};
            step_time_data = [step_time_data this_step_time_data];
        end
    end
    
    %% create figures and determine abscissae for each comparison
    comparison_variable_to_axes_index_map = zeros(number_of_comparisons, 1);
    abscissae_cell = cell(number_of_comparisons, number_of_variables_to_plot);
    
    plot_mode = study_settings.get('plot_mode');
    variables_to_plot = study_settings.get('variables_to_plot');
    conditions_to_plot = study_settings.get('conditions_to_plot');    
    conditions_control = study_settings.get('conditions_control');
    
    if strcmp(plot_mode, 'detailed') || strcmp(plot_mode, 'overview')
        % make one figure per comparison and variable
        figure_handles = zeros(number_of_comparisons, number_of_variables_to_plot);
        axes_handles = zeros(number_of_comparisons, number_of_variables_to_plot);
        pos_text_handles = zeros(number_of_comparisons, number_of_variables_to_plot);
        neg_text_handles = zeros(number_of_comparisons, number_of_variables_to_plot);
        for i_variable = 1 : number_of_variables_to_plot
            for i_comparison = 1 : number_of_comparisons
                % make figure and axes
                new_figure = figure; new_axes = axes; hold on;
                
                % store handles and determine abscissa data
                figure_handles(i_comparison, i_variable) = new_figure;
                axes_handles(i_comparison, i_variable) = new_axes;
                comparison_variable_to_axes_index_map(i_comparison) = i_comparison;
                    
                if isDiscreteVariable(i_variable, variable_data_all)
                    % abscissae gives the bin edges here
                    if dictate_axes
                        lower_bound = str2double(variables_to_plot{i_variable, 5});
                        upper_bound = str2double(variables_to_plot{i_variable, 6});
                    else
                        lower_bound = min(variable_data_all{i_variable});
                        upper_bound = max(variable_data_all{i_variable});
                    end
                    if strcmp(plot_mode, 'detailed')
%                         abscissae_cell{i_comparison, i_variable}(i_condition, :) = linspace(lower_bound, upper_bound, study_settings.get('number_of_bins_in_histogram'));
                        abscissae_cell{i_comparison, i_variable} = linspace(lower_bound, upper_bound, study_settings.get('number_of_bins_in_histogram'));
                    end
                    if strcmp(plot_mode, 'overview')
                        this_comparison = comparison_indices{i_comparison};
                        abscissae_control = 0;
                        abscissae_stimulus = 1 : length(this_comparison);
                        abscissae = {abscissae_control, abscissae_stimulus};
%                         abscissae_cell{i_comparison, i_variable}(i_condition, :) = abscissae;
                        abscissae_cell{i_comparison, i_variable} = abscissae;
                    end
                end
                if isContinuousVariable(i_variable, variable_data_all)
                    abscissa_unscaled = linspace(0, 100, study_settings.get('number_of_time_steps_normalized'));

                    % scale abscissae
                    conditions_this_comparison = comparison_indices{i_comparison};
                    step_time_means_this_comparison = zeros(size(conditions_this_comparison));
                    for i_condition = 1 : length(conditions_this_comparison)
                        % find correct condition indicator
                        condition_identifier = conditions_to_plot(conditions_this_comparison(i_condition), :);
                        stance_foot_indicator = strcmp(condition_stance_foot_list_all, condition_identifier{1});
                        perturbation_indicator = strcmp(condition_perturbation_list_all, condition_identifier{2});
                        delay_indicator = strcmp(condition_delay_list_all, condition_identifier{3});
                        index_indicator = strcmp(condition_index_list_all, condition_identifier{4});
                        experimental_indicator = strcmp(condition_experimental_list_all, condition_identifier{5});
                        stimulus_indicator = strcmp(condition_stimulus_list_all, condition_identifier{6});
                        day_indicator = strcmp(condition_day_list_all, condition_identifier{7});
                        this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
                        step_time_data_this_condition = step_time_data(:, this_condition_indicator);
                        step_time_means_this_comparison(i_condition) = mean(step_time_data_this_condition);
                    end
                    for i_condition = 1 : length(conditions_this_comparison)
                        if strcmp(study_settings.get('time_plot_style'), 'scaled_to_comparison_mean')
                            abscissa_scaled = abscissa_unscaled * mean(step_time_means_this_comparison) / 100;
                            abscissae_cell{i_comparison, i_variable}(i_condition, :) = abscissa_scaled;
                        elseif strcmp(study_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                            abscissa_scaled = abscissa_unscaled * step_time_means_this_comparison(i_condition) / 100;
                            abscissae_cell{i_comparison, i_variable}(i_condition, :) = abscissa_scaled;
                        else
                            abscissae_cell{i_comparison, i_variable}(i_condition, :) = abscissa_unscaled;
                        end
                    end                    
                end
                
                % set axes properties
                if dictate_axes
    %                 set(gca, 'xlim', [time_normalized(1), time_normalized(end)]);
                    set(gca, 'ylim', [str2double(variables_to_plot{i_variable, 5}), str2double(variables_to_plot{i_variable, 6})]);
                end
                if isDiscreteVariable(i_variable, variable_data_all) && strcmp(plot_mode, 'overview')
                    xtick = abscissae_cell{i_comparison, i_variable}{2};
                    if ~isempty(conditions_control)
                        xtick = [abscissae_cell{i_comparison, i_variable}{1} xtick]; %#ok<AGROW>
                    end
%                     set(gca, 'xlim', [-0.5 length(comparison_indices{i_comparison})+0.5]);
%                     set(gca, 'xtick', 0 : length(comparison_indices{i_comparison}));
                    set(gca, 'xlim', [-0.5 + min(xtick) 0.5 + max(xtick(end))]);
                    set(gca, 'xtick', xtick);
                end
                
                % set axis labels
                if strcmp(study_settings.get('time_plot_style'), 'scaled_to_comparison_mean') || strcmp(study_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                    xlabel('normalized time (s)');
                else
                    xlabel('normalized time (%)');
                end
                ylabel(variables_to_plot{i_variable, 3});
                
                % add text labels
                pos_text_handles(i_comparison, i_variable) = ...
                    text ...
                      ( ...
                        0, ...
                        0, ...
                        [variables_to_plot{i_variable, 7} ' $\rightarrow$'], ...
                        'rotation', 90, ...
                        'Fontsize', 24, ...
                        'horizontalalignment', 'right', ...
                        'interpreter', 'LaTeX', ...
                        'parent', new_axes ...
                      );
                neg_text_handles(i_comparison, i_variable) = ...
                    text ...
                      ( ...
                        0, ...
                        0, ...
                        ['$\leftarrow$ ' variables_to_plot{i_variable, 8}], ...
                        'rotation', 90, ...
                        'Fontsize', 24, ...
                        'horizontalalignment', 'left', ...
                        'interpreter', 'LaTeX', ...
                        'parent', new_axes...
                      );
                
                % determine title
                title_string = variables_to_plot{i_variable, 2};
                filename_string = variables_to_plot{i_variable, 4};
                for i_label = 1 : length(study_settings.get('condition_labels'))
                    if (i_label ~= study_settings.get('comparison_to_make')) ...
                        && (i_label ~= 1) ...
                        && (i_label ~= 3) ...
                        && (i_label ~= 5) ...
                        && (i_label ~= 6) ...
                        && (i_label ~= 7)
                        this_condition_label = strrep(conditions_to_plot{comparison_indices{i_comparison}(1), i_label}, '_', ' ');
                        if i_label ~= study_settings.get('comparison_to_make')
                            title_string = [title_string ' - ' this_condition_label]; %#ok<AGROW>
                            filename_string = [filename_string '_' this_condition_label];
                        end
                    end
                end
                stance_label = conditions_to_plot{comparison_indices{i_comparison}(1), 1};
                if strcmp(stance_label, 'STANCE_RIGHT')
                    title_string = [title_string ' - first step stance leg RIGHT'];
                    filename_string = [filename_string '_stanceR'];
                end
                if strcmp(stance_label, 'STANCE_LEFT')
                    title_string = [title_string ' - first step stance leg LEFT'];
                    filename_string = [filename_string '_stanceL'];
                end
                title(title_string); set(gca, 'Fontsize', 12)
                set(gcf, 'UserData', filename_string)
                
                
            end
        end
    end
    if strcmp(plot_mode, 'episodes')
        abscissae_cell_unscaled = cell(size(abscissae_cell));
        % make one figure per episode and variable
        figure_handles = zeros(number_of_episodes, number_of_variables_to_plot);
        axes_handles = zeros(number_of_episodes, number_of_variables_to_plot);
        pos_text_handles = zeros(number_of_episodes, number_of_variables_to_plot);
        neg_text_handles = zeros(number_of_episodes, number_of_variables_to_plot);

        for i_variable = 1 : number_of_variables_to_plot
            for i_episode = 1 : number_of_episodes
                % make figure and axes and store handles
                new_figure = figure; new_axes = axes; hold on;
                figure_handles(i_episode, i_variable) = new_figure;
                axes_handles(i_episode, i_variable) = new_axes;
                this_episode = episode_indices{i_episode};

                % store handles and determine abscissa data for all comparisons in this episode
                xtick = [];
                for i_comparison = 1 : length(this_episode)
                    comparison_variable_to_axes_index_map(this_episode(i_comparison)) = i_episode;
                    
                    % determine which step this is
                    this_comparison = this_episode(i_comparison);
                    conditions_this_comparison = comparison_indices{this_comparison};
                    example_condition = conditions_this_comparison(1);
                    condition_identifier = conditions_to_plot(example_condition, :);
                    gap_between_steps = 1;
                    if strcmp(condition_identifier{4}, 'ONE')
                        step_index = 1;
                    elseif strcmp(condition_identifier{4}, 'TWO')
                        step_index = 2;
                    elseif strcmp(condition_identifier{4}, 'THREE')
                        step_index = 3;
                    elseif strcmp(condition_identifier{4}, 'FOUR')
                        step_index = 4;
                    end
                    if isDiscreteVariable(i_variable, variable_data_all)
                        this_comparison = comparison_indices{i_comparison};
                        abscissae_control = (conditions_per_comparison_max + gap_between_steps) * step_index;
                        abscissae_stimulus = (1 : length(this_comparison)) + (conditions_per_comparison_max + gap_between_steps) * step_index;
                        abscissae = {abscissae_control, abscissae_stimulus};
                        abscissae_cell{this_episode(i_comparison), i_variable} = abscissae;
                        
                        if ~isempty(conditions_control)
                            xtick = [xtick abscissae{1}];
                        end
                        xtick = [xtick abscissae{2}];
                    end
                    if isContinuousVariable(i_variable, variable_data_all)
%                         abscissae_cell{this_episode(i_comparison), i_variable} = (linspace(0, 100, study_settings.get('number_of_time_steps_normalized'))) + (step_index-1)*100;
                        abscissae_cell_unscaled{this_episode(i_comparison), i_variable} = (linspace(0, 100, study_settings.get('number_of_time_steps_normalized')));
                    end
                end
                
                % set axes properties
                if dictate_axes
                    set(gca, 'ylim', [str2double(variables_to_plot{i_variable, 5}), str2double(variables_to_plot{i_variable, 6})]);
                end
                if isDiscreteVariable(i_variable, variable_data_all)
                    set(gca, 'xlim', [-0.5 + min(xtick) 0.5 + max(xtick(end))]);
                    set(gca, 'xtick', xtick);
                    set(gca, 'XTickLabelRotation', 60);
                end

                % set axis labels
                if isContinuousVariable(i_variable, variable_data_all)
                    if strcmp(study_settings.get('time_plot_style'), 'scaled_to_comparison_mean') || strcmp(study_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                        xlabel('normalized time (s)');
                    else
                        xlabel('normalized time (%)');
                    end
                end
                ylabel(variables_to_plot{i_variable, 3});
                
                % add text labels
                if ~strcmp(variables_to_plot{i_variable, 7}, '~')
                    pos_text_handles(i_episode, i_variable) = ...
                        text ...
                          ( ...
                            0, ...
                            0, ...
                            [variables_to_plot{i_variable, 7} ' $\rightarrow$'], ...
                            'rotation', 90, ...
                            'Fontsize', 24, ...
                            'horizontalalignment', 'right', ...
                            'interpreter', 'LaTeX', ...
                            'parent', new_axes ...
                          );
                end
                if ~strcmp(variables_to_plot{i_variable, 8}, '~')
                    neg_text_handles(i_episode, i_variable) = ...
                        text ...
                          ( ...
                            0, ...
                            0, ...
                            ['$\leftarrow$ ' variables_to_plot{i_variable, 8}], ...
                            'rotation', 90, ...
                            'Fontsize', 24, ...
                            'horizontalalignment', 'left', ...
                            'interpreter', 'LaTeX', ...
                            'parent', new_axes...
                          );
                end
                
                % determine title and filename
                title_string = variables_to_plot{i_variable, 2};
                filename_string = variables_to_plot{i_variable, 4};
                for i_label = 1 : length(study_settings.get('condition_labels'))
                    if (i_label ~= study_settings.get('comparison_to_make')) ...
                        && (i_label ~= 1) ...
                        && (i_label ~= 4) ...
                        && (i_label ~= 5) ...
                        && (i_label ~= 6) ...
                        this_condition_label = strrep(conditions_to_plot{comparison_indices{i_comparison}(1), i_label}, '_', ' ');
                        if ~strcmp(this_condition_label, 'N/A')
                            title_string = [title_string ' - ' this_condition_label]; %#ok<AGROW>
                            filename_string = [filename_string '_' this_condition_label];
                        end
                    end
                end
                first_comparison = this_episode(1);
                conditions_first_comparison = comparison_indices{first_comparison};
                example_condition = conditions_first_comparison(1);
                condition_identifier = conditions_to_plot(example_condition, :);
                if strcmp(condition_identifier{1}, 'STANCE_RIGHT')
                    title_string = [title_string ' - first step stance leg RIGHT'];
                    filename_string = [filename_string '_stanceR'];
                end
                if strcmp(condition_identifier{1}, 'STANCE_LEFT')
                    title_string = [title_string ' - first step stance leg LEFT'];
                    filename_string = [filename_string '_stanceL'];
                end
                
                title(title_string); set(gca, 'Fontsize', 12)
                set(gcf, 'UserData', filename_string)
                
                
            end
        end
        
        % calculate average step times and scale abscissa
        for i_variable = 1 : number_of_variables_to_plot
            for i_episode = 1 : number_of_episodes
                this_episode = episode_indices{i_episode};
                for i_comparison = 1 : length(this_episode)
                    
                    % determine which step this is
                    this_comparison = this_episode(i_comparison);
                    conditions_this_comparison = comparison_indices{this_comparison};
                    step_time_means_this_comparison = zeros(size(conditions_this_comparison));
                    for i_condition = 1 : length(conditions_this_comparison)
                        % find correct condition indicator
                        condition_identifier = conditions_to_plot(conditions_this_comparison(i_condition), :);
                        stance_foot_indicator = strcmp(condition_stance_foot_list_all, condition_identifier{1});
                        perturbation_indicator = strcmp(condition_perturbation_list_all, condition_identifier{2});
                        delay_indicator = strcmp(condition_delay_list_all, condition_identifier{3});
                        index_indicator = strcmp(condition_index_list_all, condition_identifier{4});
                        experimental_indicator = strcmp(condition_experimental_list_all, condition_identifier{5});
                        stimulus_indicator = strcmp(condition_stimulus_list_all, condition_identifier{6});
                        day_indicator = strcmp(condition_day_list_all, condition_identifier{7});
                        this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
                        step_time_data_this_condition = step_time_data(:, this_condition_indicator);
                        
                        % calculate average step time
                        step_time_means_this_comparison(i_condition) = mean(step_time_data_this_condition);
                    end
                    
                    % scale abscissa
                    if isContinuousVariable(i_variable, variable_data_all)
                        abscissa_unscaled = abscissae_cell_unscaled{this_episode(i_comparison), i_variable};
                        for i_condition = 1 : length(conditions_this_comparison)
                            if strcmp(study_settings.get('time_plot_style'), 'scaled_to_comparison_mean')
                                abscissa_scaled = abscissa_unscaled * mean(step_time_means_this_comparison) / 100;
                                abscissae_cell{this_episode(i_comparison), i_variable}(i_condition, :) = abscissa_scaled;
                            elseif strcmp(study_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                                abscissa_scaled = abscissa_unscaled * step_time_means_this_comparison(i_condition) / 100;
                                abscissae_cell{this_episode(i_comparison), i_variable}(i_condition, :) = abscissa_scaled;
                            else
                                abscissae_cell{this_episode(i_comparison), i_variable}(i_condition, :) = abscissa_unscaled;
                            end
                        end                    
                    end                    
                end
            end
        end
        
        % determine abscissa offsets
        for i_step = 2 : 4
            for i_variable = 1 : number_of_variables_to_plot
                if isContinuousVariable(i_variable, variable_data_all)
                    for i_episode = 1 : number_of_episodes
                        this_episode = episode_indices{i_episode};
                        for i_comparison = 1 : length(this_episode)

                            % determine which step this is
                            this_comparison = this_episode(i_comparison);
                            conditions_this_comparison = comparison_indices{this_comparison};
                            for i_condition = 1 : length(conditions_this_comparison)
                                % determine step index
                                condition_identifier = conditions_to_plot(conditions_this_comparison(i_condition), :);
                                if strcmp(condition_identifier{4}, 'ONE')
                                    step_index = 1;
                                elseif strcmp(condition_identifier{4}, 'TWO')
                                    step_index = 2;
                                    previous_step_label = 'ONE';
                                elseif strcmp(condition_identifier{4}, 'THREE')
                                    step_index = 3;
                                    previous_step_label = 'TWO';
                                elseif strcmp(condition_identifier{4}, 'FOUR')
                                    step_index = 4;
                                    previous_step_label = 'THREE';
                                end

                                if step_index == i_step
                                    % find condition index for previous step in same condition
                                    previous_step_condition_index = [];
                                    previous_step_comparison_index = [];
                                    for j_comparison = 1 : length(this_episode)
                                        candidate_comparison = this_episode(j_comparison);
                                        conditions_candidate_comparison = comparison_indices{candidate_comparison};

                                        for j_condition = 1 : length(conditions_candidate_comparison)
                                            candidate_condition_identifier = conditions_to_plot(conditions_candidate_comparison(j_condition), :);


                                            if strcmp(condition_identifier{2}, candidate_condition_identifier{2}) ...
                                            && strcmp(condition_identifier{3}, candidate_condition_identifier{3}) ...
                                            && strcmp(previous_step_label, candidate_condition_identifier{4}) ...
                                            && strcmp(condition_identifier{5}, candidate_condition_identifier{5}) ...
                                            && strcmp(condition_identifier{6}, candidate_condition_identifier{6}) ...
                                            && strcmp(condition_identifier{7}, candidate_condition_identifier{7})
                                                previous_step_comparison_index = j_comparison;
                                                previous_step_condition_index = j_condition;
                                            end
                                        end
                                    end

                                    % find out where abscissa for previous step ends
                                    previous_step_last_data_point = abscissae_cell{this_episode(previous_step_comparison_index), i_variable}(previous_step_condition_index, end);
                                    abscissae_cell{this_episode(i_comparison), i_variable}(i_condition, :) = abscissae_cell{this_episode(i_comparison), i_variable}(i_condition, :) + previous_step_last_data_point;


                                end

                            end
                        end
                    end
                end
            end
        end
        
        % set x-limits
        for i_variable = 1 : number_of_variables_to_plot
            if isContinuousVariable(i_variable, variable_data_all)
                for i_episode = 1 : number_of_episodes
                    % determine time window to show
                    this_episode = episode_indices{i_episode};
                    first_comparison_in_episode_index = this_episode(1);
                    first_step_abscissae = abscissae_cell{first_comparison_in_episode_index, i_variable};
                    episode_start_time = first_step_abscissae(1, 1);
                    last_comparison_in_episode_index = this_episode(end);
                    last_step_abscissae = abscissae_cell{last_comparison_in_episode_index, i_variable};
                    episode_end_time = last_step_abscissae(1, end);

                    % set x-limits accordingly
                    set(axes_handles(i_episode, i_variable), 'xlim', [episode_start_time episode_end_time]);
                end
            end
        end
    end
    
    %% plot data
    for i_variable = 1 : number_of_variables_to_plot
        if study_settings.get('plot_response')
            data_to_plot = response_data_all{i_variable, 1};
        else
            data_to_plot = variable_data_all{i_variable, 1};
        end
        
        colors_comparison = study_settings.get('colors_comparison');
        for i_comparison = 1 : length(comparison_indices)
            % find correct condition indicator for control
            conditions_this_comparison = comparison_indices{i_comparison};
            top_level_plots = [];
            target_axes_handle = axes_handles(comparison_variable_to_axes_index_map(i_comparison), i_variable);
            
            % plot control
            if study_settings.get('plot_control') && ~isempty(conditions_control) && ~study_settings.get('plot_response')
                % determine which control condition applies here
                representant_condition_index = conditions_this_comparison(1);
                applicable_control_condition_index = findApplicableControlConditionIndex(conditions_to_plot(representant_condition_index, :), conditions_control);
                applicable_control_condition_labels = conditions_control(applicable_control_condition_index, :);
                
                % extract data for control condition
                stance_foot_indicator = strcmp(condition_stance_foot_list_all, applicable_control_condition_labels{1});
                perturbation_indicator = strcmp(condition_perturbation_list_all, applicable_control_condition_labels{2});
                delay_indicator = strcmp(condition_delay_list_all, applicable_control_condition_labels{3});
                index_indicator = strcmp(condition_index_list_all, applicable_control_condition_labels{4});
                experimental_indicator = strcmp(condition_experimental_list_all, applicable_control_condition_labels{5});
                stimulus_indicator = strcmp(condition_stimulus_list_all, applicable_control_condition_labels{6});
                day_indicator = strcmp(condition_day_list_all, applicable_control_condition_labels{7});
                this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
                data_to_plot_this_condition = data_to_plot(:, this_condition_indicator);
                
                if ~isempty(data_to_plot_this_condition)
                    if isDiscreteVariable(i_variable, variable_data_all)
                        target_abscissa = abscissae_cell{i_comparison, i_variable};
                        if strcmp(plot_mode, 'detailed')
                            histogram ...
                              ( ...
                                target_axes_handle, ...
                                data_to_plot_this_condition, ...
                                'binEdges', target_abscissa, ...
                                'edgecolor', study_settings.get('color_control'), ...
                                'facecolor', lightenColor(study_settings.get('color_control'), 0.5), ...
                                'DisplayName', 'CONTROL' ...
                              );
                        end
                        if strcmp(plot_mode, 'overview') || strcmp(plot_mode, 'episodes')
                            singleBoxPlot(target_axes_handle, target_abscissa{1}, data_to_plot_this_condition, study_settings.get('color_control'), 'CONTROL', false)
                        end
                    end
                    if isContinuousVariable(i_variable, variable_data_all)
                        target_abscissa = abscissae_cell{i_comparison, i_variable}(i_condition, :);
                        if strcmp(plot_mode, 'detailed')
                            % individual trajectories
                            plot ...
                              ( ...
                                target_axes_handle, ...
                                target_abscissa, ...
                                data_to_plot_this_condition, ...
                                'HandleVisibility', 'off', ...
                                'color', lightenColor(study_settings.get('color_control'), 0.5) ...
                              );
                            % condition average
                            control_mean_plot = plot ...
                              ( ...
                                target_axes_handle, ...
                                target_abscissa, ...
                                mean(data_to_plot_this_condition, 2), ...
                                'DisplayName', 'CONTROL', ...
                                'linewidth', 5, ...
                                'color', study_settings.get('color_control') ...
                              );
                            top_level_plots = [top_level_plots control_mean_plot]; %#ok<AGROW>
                        end
                        if strcmp(plot_mode, 'overview') || strcmp(plot_mode, 'episodes')
                            plot_handles = shadedErrorBar ...
                              ( ...
                                target_abscissa, ...
                                mean(data_to_plot_this_condition, 2), ...
                                cinv(data_to_plot_this_condition, 2), ...
                                { ...
                                  'color', study_settings.get('color_control'), ...
                                  'linewidth', 6 ...
                                }, ...
                                1, ...
                                target_axes_handle ...
                              );
                            set(plot_handles.edge, 'HandleVisibility', 'off');
                            set(plot_handles.patch, 'HandleVisibility', 'off');
                            set(plot_handles.mainLine, 'DisplayName', 'CONTROL');

                            top_level_plots = [top_level_plots plot_handles.mainLine]; %#ok<AGROW>
                        end
                    end
                end
            end
            
            % plot stimulus
            for i_condition = 1 : length(conditions_this_comparison)
                label_string = strrep(conditions_to_plot{comparison_indices{i_comparison}(i_condition), study_settings.get('comparison_to_make')}, '_', ' ');
                
                % find correct condition indicator
                condition_identifier = conditions_to_plot(conditions_this_comparison(i_condition), :);
                stance_foot_indicator = strcmp(condition_stance_foot_list_all, condition_identifier{1});
                perturbation_indicator = strcmp(condition_perturbation_list_all, condition_identifier{2});
                delay_indicator = strcmp(condition_delay_list_all, condition_identifier{3});
                index_indicator = strcmp(condition_index_list_all, condition_identifier{4});
                experimental_indicator = strcmp(condition_experimental_list_all, condition_identifier{5});
                stimulus_indicator = strcmp(condition_stimulus_list_all, condition_identifier{6});
                day_indicator = strcmp(condition_day_list_all, condition_identifier{7});
                this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
                data_to_plot_this_condition = data_to_plot(:, this_condition_indicator);
                origin_trial_list_this_condition = origin_trial_list_all(this_condition_indicator);
                if isDiscreteVariable(i_variable, variable_data_all)
                    target_abscissa = abscissae_cell{i_comparison, i_variable};
                    if strcmp(plot_mode, 'detailed')
                        histogram ...
                          ( ...
                            target_axes_handle, ...
                            data_to_plot_this_condition, ...
                            target_abscissa, ...
                            'edgecolor', colors_comparison(i_condition, :), ...
                            'facecolor', lightenColor(colors_comparison(i_condition, :), 0.5), ...
                            'DisplayName', label_string ...
                          );
                    end
                    if strcmp(plot_mode, 'overview') || strcmp(plot_mode, 'episodes')
                        if ~any(isnan(data_to_plot_this_condition))
                            singleBoxPlot(target_axes_handle, target_abscissa{2}(i_condition), data_to_plot_this_condition, colors_comparison(i_condition, :), label_string, false)
                        end
                    end
                end
                if isContinuousVariable(i_variable, variable_data_all)
                    target_abscissa = abscissae_cell{i_comparison, i_variable}(i_condition, :);
                    if strcmp(plot_mode, 'detailed')
                        for i_stretch = 1 : size(data_to_plot_this_condition, 2)
                            origin_trial_data = ones(size(target_abscissa)) * origin_trial_list_this_condition(i_stretch);
                            plot3 ...
                              ( ...
                                target_axes_handle, ...
                                target_abscissa, ...
                                data_to_plot_this_condition(:, i_stretch), ...
                                origin_trial_data, ...
                                'HandleVisibility', 'off', ...
                                'color', lightenColor(colors_comparison(i_condition, :), 0.5) ...
                              );
                        end
                    end
                    if strcmp(plot_mode, 'overview') || strcmp(plot_mode, 'episodes')
                        plot_handles = shadedErrorBar ...
                          ( ...
                            target_abscissa, ...
                            mean(data_to_plot_this_condition, 2), ...
                            cinv(data_to_plot_this_condition, 2), ...
                            { ...
                              'color', colors_comparison(i_condition, :), ...
                              'linewidth', 6 ...
                            }, ...
                            1, ...
                            target_axes_handle ...
                          );
                        top_level_plots = [top_level_plots plot_handles.mainLine]; %#ok<AGROW>
                        set(plot_handles.edge, 'HandleVisibility', 'off');
                        set(plot_handles.patch, 'HandleVisibility', 'off');
                        set(plot_handles.mainLine, 'DisplayName', label_string);
                        if ~strcmp(condition_identifier{4}, 'ONE')
                            set(plot_handles.mainLine, 'HandleVisibility', 'off');
                        end
                    end
                end
            end

            % reorder to bring mean plots on top
            for i_plot = 1 : length(top_level_plots)
                uistack(top_level_plots(i_plot), 'top');
            end
            
            % toggle legend
            if show_legend && ~(isDiscreteVariable(i_variable, variable_data_all) && (strcmp(plot_mode, 'overview') || strcmp(plot_mode, 'episodes')))
                legend(target_axes_handle, 'show')
            end
        end
    end
    
    %% update label positions
    for i_variable = 1 : number_of_variables_to_plot
        for i_axes = 1 : size(axes_handles, 1);
            these_axes = axes_handles(i_axes, i_variable);
            xlimits = get(these_axes, 'xlim'); ylimits = get(these_axes, 'ylim');
            if pos_text_handles(i_axes, i_variable) ~= 0
                pos_text_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.12;
                pos_text_position_y = ylimits(2);
                set(pos_text_handles(i_axes, i_variable), 'Position', [pos_text_position_x pos_text_position_y]);
            end
            if neg_text_handles(i_axes, i_variable) ~= 0
                neg_text_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.12;
                neg_text_position_y = ylimits(1);
                set(neg_text_handles(i_axes, i_variable), 'Position', [neg_text_position_x neg_text_position_y]);
            end
        end
    end
    
    %% save figures
    if parser.Results.save
        % figure out folders
        if ~exist('figures', 'dir')
            mkdir('figures')
        end
        if ~exist(['figures' filesep 'withLabels'], 'dir')
            mkdir(['figures' filesep 'withLabels'])
        end
        if ~exist(['figures' filesep 'noLabels'], 'dir')
            mkdir(['figures' filesep 'noLabels'])
        end
        for i_figure = 1 : length(figure_handles)
            % save with labels
            legend(axes_handles(i_figure), 'show');
            filename = ['figures' filesep 'withLabels' filesep get(figure_handles(i_figure), 'UserData')];
            saveas(figure_handles(i_figure), filename, parser.Results.format)
            
            % save without labels
%             zero_plot = plot(get(axes_handles(i_figure), 'xlimits'), [0 0], 'color', [0.7 0.7 0.7]);
%             uistack(zero_plot, 'bottom')
%             set(postext, 'visible', 'off');
%             set(negtext, 'visible', 'off');
            set(get(axes_handles(i_figure), 'xaxis'), 'visible', 'off');
            set(get(axes_handles(i_figure), 'yaxis'), 'visible', 'off');
            set(get(axes_handles(i_figure), 'xlabel'), 'visible', 'off');
            set(get(axes_handles(i_figure), 'ylabel'), 'visible', 'off');
            set(get(axes_handles(i_figure), 'title'), 'visible', 'off');
            set(axes_handles(i_figure), 'xticklabel', '');
            set(axes_handles(i_figure), 'yticklabel', '');
            set(axes_handles(i_figure), 'position', [0 0 1 1]);
            legend(axes_handles(i_figure), 'hide');
            filename = ['figures' filesep 'noLabels' filesep get(figure_handles(i_figure), 'UserData')];
                saveas(figure_handles(i_figure), filename, parser.Results.format)

            close(figure_handles(i_figure))            
        end
    end
end

%% helper functions


function [comparison_indices, conditions_per_comparison_max] = determineComparisons(study_settings)
    % initialize
    conditions_to_plot = study_settings.get('conditions_to_plot');
    number_of_conditions_to_plot = size(conditions_to_plot, 1);
    comparison_indices = {};
    conditions_already_compared = [];
    conditions_per_comparison_max = 0;
    
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
                comparison_table = zeros(1, length(study_settings.get('condition_labels'))); % this is a table indicating equality between the two conditions in questions
                for i_label = 1 : length(study_settings.get('condition_labels'))
                    comparison_table(i_label) = strcmp(conditions_to_plot{i_condition, i_label}, conditions_to_plot{j_condition, i_label});
                end

                % look at the relevant entries of the comparison table
                comparison_table_relevant = comparison_table;
                comparison_table_relevant(study_settings.get('comparison_to_make')) = [];
                if all(comparison_table_relevant)
                    this_comparison = [this_comparison, j_condition]; %#ok<AGROW>
                end
            end
        end
        comparison_indices = [comparison_indices; this_comparison]; %#ok<AGROW>
        conditions_already_compared = [conditions_already_compared this_comparison]; %#ok<AGROW>
        
        if conditions_per_comparison_max < length(this_comparison)
            conditions_per_comparison_max = length(this_comparison);
        end
    end    
    if ~isempty(study_settings.get('conditions_control'))
        conditions_per_comparison_max = conditions_per_comparison_max + 1;
    end
end

function episode_indices = determineEpisodes(study_settings, comparison_indices)
    conditions_to_plot = study_settings.get('conditions_to_plot');
    condition_column_index = find(strcmp(study_settings.get('condition_labels'), 'index'));
    condition_column_stancefoot = find(strcmp(study_settings.get('condition_labels'), 'stance foot'));
%     episode_first_stretch_indices = find(strcmp(study_settings.get('conditions_to_plot(:, 4), 'ONE'));
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
        example_condition_in_base_comparison_labels = conditions_to_plot(example_condition_in_base_comparison, :);
        for j_comparison = 1 : number_of_comparisons
            if i_comparison ~= j_comparison
                this_comparison = comparison_indices{j_comparison};
                example_condition_in_this_comparison = this_comparison(1);
                example_condition_in_this_comparison_labels = conditions_to_plot(example_condition_in_this_comparison, :);
                % check which conditions labels agree between these two conditions
                comparison_table = zeros(1, length(study_settings.get('condition_labels'))); % this is a table indicating equality between the two conditions in questions
                for i_label = 1 : length(study_settings.get('condition_labels'))
                    comparison_table(i_label) = strcmp(example_condition_in_base_comparison_labels{i_label}, example_condition_in_this_comparison_labels{i_label});
                end

                % look at the relevant entries of the comparison table
                comparison_table_relevant = comparison_table;
                comparison_table_relevant([condition_column_stancefoot condition_column_index study_settings.get('comparison_to_make')]) = [];
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
                    end
                    if strcmp(example_condition_in_base_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                        if strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'TWO') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'THREE') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'FOUR') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                            comparison_indices_in_this_episode = [comparison_indices_in_this_episode, j_comparison]; %#ok<AGROW>
                        end
                    end
                    if strcmp(example_condition_in_base_comparison_labels(condition_column_stancefoot), 'STANCE_BOTH')
                        if strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'TWO') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
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

function singleBoxPlot(target_axes_handle, abscissa, data, color, label, show_outliers)
    if nargin < 6
        show_outliers = true;
    end

    % set some parameters, these should be name-value pair arguments later
    width = 0.8;
    
    % extract data
    data_median = median(data);
    data_quartile_1 = prctile(data, 25);
    data_quartile_3 = prctile(data, 75);
    data_iqr = iqr(data);
    data_upper_inner_fence = data_quartile_3 + 1.5*data_iqr;
    data_lower_inner_fence = data_quartile_1 - 1.5*data_iqr;
    data_upper_adjacent = max([data(data<=data_upper_inner_fence) -inf]);
    data_lower_adjacent = min([data(data>=data_lower_inner_fence) inf]);
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
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off' ...
      );
    plot(target_axes_handle, abscissa + width*[-0.5 0.5], [data_median data_median], 'color', 'k', 'HandleVisibility', 'off'); % median
    plot(target_axes_handle, [abscissa abscissa], [data_quartile_3 data_upper_adjacent], 'k--', 'HandleVisibility', 'off'); % upper range
    plot(target_axes_handle, [abscissa abscissa], [data_lower_adjacent data_quartile_1], 'k--', 'HandleVisibility', 'off'); % lower range
    plot(target_axes_handle, abscissa+width*[-0.25 0.25], [data_lower_adjacent data_lower_adjacent], 'k-', 'HandleVisibility', 'off'); % max
    plot(target_axes_handle, abscissa+width*[-0.25 0.25], [data_upper_adjacent data_upper_adjacent], 'k-', 'HandleVisibility', 'off'); % min
    if show_outliers
        plot(target_axes_handle, abscissa * ones(size(outliers)), outliers, '+', 'color', [1; 1; 1] * 0.7, 'HandleVisibility', 'off');
    end
    
    % labels
    xtick = get(target_axes_handle, 'xtick');
    xticklabels = get(target_axes_handle, 'xticklabel');
    xticklabels{xtick == abscissa} = label;
    set(target_axes_handle, 'xticklabel', xticklabels);
end




