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
    study_settings = loadSettingsFile(study_settings_file);
    

    %% determine subjects and data folders
    [data_folder_list, subject_list] = determineDataStructure(subjects);
    [comparison_indices, conditions_per_comparison_max] = determineComparisons(study_settings);
    number_of_comparisons = length(comparison_indices);
    episode_indices = determineEpisodes(study_settings, comparison_indices);
    number_of_episodes = length(episode_indices);
    
    %% collect data from all data folders
    number_of_variables_to_plot = size(study_settings.variables_to_plot, 1);
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

    for i_folder = 1 : length(data_folder_list)
        % load data
        data_path = data_folder_list{i_folder};
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
        load([data_path filesep 'analysis' filesep date '_' subject_id '_results.mat']);

        % append data from this subject to containers for all subjects
        condition_stance_foot_list_all = [condition_stance_foot_list_all; condition_stance_foot_list_subject]; %#ok<AGROW>
        condition_perturbation_list_all = [condition_perturbation_list_all; condition_perturbation_list_subject]; %#ok<AGROW>
        condition_delay_list_all = [condition_delay_list_all; condition_delay_list_subject]; %#ok<AGROW>
        condition_index_list_all = [condition_index_list_all; condition_index_list_subject]; %#ok<AGROW>
        condition_experimental_list_all = [condition_experimental_list_all; condition_experimental_list_subject]; %#ok<AGROW>
        condition_stimulus_list_all = [condition_stimulus_list_all; condition_stimulus_list_subject]; %#ok<AGROW>
        condition_day_list_all = [condition_day_list_all; condition_day_list_subject]; %#ok<AGROW>
        origin_trial_list_all = [origin_trial_list_all; origin_trial_list_subject];
        origin_start_time_list_all = [origin_start_time_list_all; origin_start_time_list_subject];
        origin_end_time_list_all = [origin_end_time_list_all; origin_end_time_list_subject];
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
                    if strcmp(study_settings.plot_mode, 'detailed')
                        abscissae_cell{i_comparison, i_variable} = linspace(lower_bound, upper_bound, study_settings.number_of_bins_in_histogram);
                    end
                    if strcmp(study_settings.plot_mode, 'overview')
                        this_comparison = comparison_indices{i_comparison};
                        abscissae_control = 0;
                        abscissae_stimulus = 1 : length(this_comparison);
                        abscissae = {abscissae_control, abscissae_stimulus};
                        abscissae_cell{i_comparison, i_variable} = abscissae;
                    end
                end
                if isContinuousVariable(i_variable, variable_data_all)
                    abscissae_cell{i_comparison, i_variable} = 1 : 100;
                end
                
                % set axes properties
                if dictate_axes
    %                 set(gca, 'xlim', [time_normalized(1), time_normalized(end)]);
                    set(gca, 'ylim', [str2double(study_settings.variables_to_plot{i_variable, 5}), str2double(study_settings.variables_to_plot{i_variable, 6})]);
                end
                if isDiscreteVariable(i_variable, variable_data_all) && strcmp(study_settings.plot_mode, 'overview')
                    xtick = abscissae_cell{i_comparison, i_variable}{2};
                    if ~isempty(study_settings.conditions_control)
                        xtick = [abscissae_cell{i_comparison, i_variable}{1} xtick]; %#ok<AGROW>
                    end
%                     set(gca, 'xlim', [-0.5 length(comparison_indices{i_comparison})+0.5]);
%                     set(gca, 'xtick', 0 : length(comparison_indices{i_comparison}));
                    set(gca, 'xlim', [-0.5 + min(xtick) 0.5 + max(xtick(end))]);
                    set(gca, 'xtick', xtick);
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
                xtick = [];
                for i_comparison = 1 : length(this_episode)
                    comparison_variable_to_axes_index_map(this_episode(i_comparison)) = i_episode;
                    
                    % determine which step this is
                    this_comparison = this_episode(i_comparison);
                    conditions_in_this_comparison = comparison_indices{this_comparison};
                    example_condition = conditions_in_this_comparison(1);
                    condition_identifier = study_settings.conditions_to_plot(example_condition, :);
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
                        
                        if ~isempty(study_settings.conditions_control)
                            xtick = [xtick abscissae{1}];
                        end
                        xtick = [xtick abscissae{2}];
                    end
                    if isContinuousVariable(i_variable, variable_data_all)
                        abscissae_cell{this_episode(i_comparison), i_variable} = (1 : 100) + (step_index-1)*100;
                    end
                end
                
                % set axes properties
                if dictate_axes
    %                 set(gca, 'xlim', [time_normalized(1), time_normalized(end)]);
                    set(gca, 'ylim', [str2double(study_settings.variables_to_plot{i_variable, 5}), str2double(study_settings.variables_to_plot{i_variable, 6})]);
                end
                if isDiscreteVariable(i_variable, variable_data_all)
                    set(gca, 'xlim', [-0.5 + min(xtick) 0.5 + max(xtick(end))]);
                    set(gca, 'xtick', xtick);
                    set(gca, 'XTickLabelRotation', 60);
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
                stimulus_indicator = strcmp(condition_stimulus_list_all, applicable_control_condition_labels{6});
                day_indicator = strcmp(condition_day_list_all, applicable_control_condition_labels{7});
                this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
                data_to_plot_this_condition = data_to_plot(:, this_condition_indicator);
                
                if ~isempty(data_to_plot_this_condition)
                    if isDiscreteVariable(i_variable, variable_data_all)
                        if strcmp(study_settings.plot_mode, 'detailed')
                            histogram ...
                              ( ...
                                target_axes_handle, ...
                                data_to_plot_this_condition, ...
                                'binEdges', target_abscissa, ...
                                'edgecolor', study_settings.color_control, ...
                                'facecolor', lightenColor(study_settings.color_control, 0.5), ...
                                'DisplayName', 'CONTROL' ...
                              );
                        end
                        if strcmp(study_settings.plot_mode, 'overview') || strcmp(study_settings.plot_mode, 'episodes')
                            singleBoxPlot(target_axes_handle, target_abscissa{1}, data_to_plot_this_condition, study_settings.color_control, 'CONTROL')
                        end
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
                label_string = strrep(study_settings.conditions_to_plot{comparison_indices{i_comparison}(i_condition), study_settings.comparison_to_make}, '_', ' ');
                
                % find correct condition indicator
                condition_identifier = study_settings.conditions_to_plot(conditions_this_comparison(i_condition), :);
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
                    if strcmp(study_settings.plot_mode, 'detailed')
                        histogram ...
                          ( ...
                            target_axes_handle, ...
                            data_to_plot_this_condition, ...
                            target_abscissa, ...
                            'edgecolor', study_settings.colors_comparison(i_condition, :), ...
                            'facecolor', lightenColor(study_settings.colors_comparison(i_condition, :), 0.5), ...
                            'DisplayName', label_string ...
                          );
                    end
                    if strcmp(study_settings.plot_mode, 'overview') || strcmp(study_settings.plot_mode, 'episodes')
                        if ~any(isnan(data_to_plot_this_condition))
                            singleBoxPlot(target_axes_handle, target_abscissa{2}(i_condition), data_to_plot_this_condition, study_settings.colors_comparison(i_condition, :), label_string)
                        end
                    end
                end
                if isContinuousVariable(i_variable, variable_data_all)
                    if strcmp(study_settings.plot_mode, 'detailed')
                        for i_stretch = 1 : size(data_to_plot_this_condition, 2)
                            origin_trial_data = ones(size(target_abscissa)) * origin_trial_list_this_condition(i_stretch);
                            plot3 ...
                              ( ...
                                target_axes_handle, ...
                                target_abscissa, ...
                                data_to_plot_this_condition(:, i_stretch), ...
                                origin_trial_data, ...
                                'HandleVisibility', 'off', ...
                                'color', lightenColor(study_settings.colors_comparison(i_condition, :), 0.5) ...
                              );
                        end
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
                        set(plot_handles.edge, 'HandleVisibility', 'off');
                        set(plot_handles.patch, 'HandleVisibility', 'off');
                        set(plot_handles.mainLine, 'DisplayName', label_string);
                        
                    end
                end
            end

            % reorder to bring mean plots on top
            for i_plot = 1 : length(top_level_plots)
                uistack(top_level_plots(i_plot), 'top');
            end
            
            % toggle legend
            if show_legend && ~(isDiscreteVariable(i_variable, variable_data_all) && (strcmp(study_settings.plot_mode, 'overview') || strcmp(study_settings.plot_mode, 'episodes')))
                legend(target_axes_handle, 'show')
            end
        end
    end
    
end

%% helper functions

function [data_folder_list, subject_list] = determineDataStructure(subjects)
    % determine subject list and current folder type
    if exist('subjectSettings.txt', 'file')
        % data folders contain subjectSettings.txt files
        current_folder_type = 'data';
    else
        % we're not in a data folder
        if exist('studySettings.txt', 'file')
            % study folders contain studySettings.txt files
            current_folder_type = 'study';
        else
            % we're neither in a data folder nor in a study folder
            if exist(['..' filesep 'studySettings.txt'], 'file')
                % but one level above is a study folder, so we must be in a subject folder containing data folders
                current_folder_type = 'subject';
            else
                error('Something is wrong with the folder structure. Current folder does not appear to be a study, subject or data folder.')
            end
        end
    end
    
    % determine data folders
    if strcmp(current_folder_type, 'data')
        path_split = strsplit(pwd, filesep);
        subject_list = path_split(end-1);
        data_folder_list = {pwd};
    end
    if strcmp(current_folder_type, 'subject')
        % get list of everything in the current directory
        things_in_current_folder = dir;
        % get a logical vector that tells which is a directory
        dir_flags = [things_in_current_folder.isdir];
        % extract only those that are directories.
        dir_name = {things_in_current_folder.name};
        folder_list = dir_name(dir_flags);
        % remove pointers to upper level directories
        folder_list(1:2) = [];
        
        % store data folders in lists
        number_of_data_folders = length(folder_list);
        path_split = strsplit(pwd, filesep);
        subject = path_split{end};
        subject_list = cell(number_of_data_folders, 1);
        data_folder_list = cell(number_of_data_folders, 1);
        for i_folder = 1 : number_of_data_folders
            subject_list{i_folder} = subject;
            data_folder_list{i_folder} = [pwd filesep folder_list{i_folder}];
        end
    end
    if strcmp(current_folder_type, 'study')
        % if no list was passed, load from subjects.csv
        if isempty(subjects)
            % no list passed, but current folder is a data folder, so plot this data
            subject_data_file = 'subjects.csv';
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
        
        % check each subject folder for data folders
        subject_list = {};
        data_folder_list = {};
        for i_subject = 1 : length(subjects)
            subject_path = [pwd filesep subjects{i_subject}];
            if exist([subject_path filesep 'subjectSettings.txt'], 'file')
                % subject folder is also a data folder
                subject_list = [subject_list; subjects{i_subject}];
                data_folder_list = [data_folder_list; subject_path];
            else
                % subject folder contains data folders
                things_in_subject_folder = dir(subject_path);
                % get a logical vector that tells which is a directory
                dir_flags = [things_in_subject_folder.isdir];
                % extract only those that are directories.
                dir_name = {things_in_subject_folder.name};
                folder_list = dir_name(dir_flags);
                % remove pointers to upper level directories
                folder_list(1:2) = [];

                % store data folders in lists
                number_of_data_folders = length(folder_list);
                path_split = strsplit(subject_path, filesep);
                subject = path_split{end};
                subject_list = cell(number_of_data_folders, 1);
                data_folder_list = cell(number_of_data_folders, 1);
                for i_folder = 1 : number_of_data_folders
                    subject_list{i_folder} = subject;
                    data_folder_list{i_folder} = [subject_path filesep folder_list{i_folder}];
                end                
            end
        end
    end
    
    
    
    
    % old way to determine subjects, delete later
%     if strcmp(path_split(end), subjects{i_subject})
%         % we're already in the subject folder, so just load from here
%         data_path = '';
%     end
%     if true % TODO: change this when determineDataStructure has worked
%         % we're in the study root, so load from subject folder
%         data_path = [subjects{i_subject} filesep];
%     end
%     load([data_path 'subjectInfo.mat'], 'date', 'subject_id');
%     load([data_path 'analysis' filesep date '_' subject_id '_results.mat']);
% 
% 
% 
%     % determine subject list
%     subject_list_determined = false;
%     if ~isempty(subjects)
%         subject_list_determined = true;
%     end
%     if ~subject_list_determined && exist('subjectInfo.mat', 'file')
%         % no list passed, but current folder is a data folder, so plot this data
%         subject_data_file = 'subjects.csv';
%         format = '%s';
%         fid = fopen(subject_data_file);
%         fgetl(fid);
%         fgetl(fid);
%         data_raw = textscan(fid, format);
%         fclose(fid);
% 
%         % transform to cell
%         data_lines = data_raw{1};
%         data_cell = {};
%         for i_line = 1 : length(data_lines)
%             line_split = strsplit(data_lines{i_line}, ',');
%             data_cell = [data_cell; line_split]; %#ok<AGROW>
%         end
% 
%         subjects = data_cell(:, 1);
%     end
%     
%     
%     
%     if ~subject_list_determined && (exist(['..' filesep 'subjects.csv'], 'file') || exist(['..' filesep '..' filesep 'subjects.csv'], 'file'))
%         % no list passed, there's a subject list one or two levels up, load from there
%         if exist(['..' filesep 'subjects.csv'], 'file')
%             subject_data_file = ['..' filesep 'subjects.csv'];
%         end
%         if exist(['..' filesep '..' filesep 'subjects.csv'], 'file')
%             subject_data_file = ['..' filesep '..' filesep 'subjects.csv'];
%         end
%         format = '%s';
%         fid = fopen(subject_data_file);
%         fgetl(fid);
%         fgetl(fid);
%         data_raw = textscan(fid, format);
%         fclose(fid);
% 
%         % transform to cell
%         data_lines = data_raw{1};
%         data_cell = {};
%         for i_line = 1 : length(data_lines)
%             line_split = strsplit(data_lines{i_line}, ',');
%             data_cell = [data_cell; line_split]; %#ok<AGROW>
%         end
% 
%         subjects = data_cell(:, 1);
%     end
%     if ischar(subjects)
%         % single string (i.e. char array) was passed, make a cell out of this
%         subjects = {subjects};
%     end

end

function [comparison_indices, conditions_per_comparison_max] = determineComparisons(study_settings)
    % initialize
    number_of_conditions_to_plot = size(study_settings.conditions_to_plot, 1);
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
        
        if conditions_per_comparison_max < length(this_comparison)
            conditions_per_comparison_max = length(this_comparison);
        end
    end    
    if ~isempty(study_settings.conditions_control)
        conditions_per_comparison_max = conditions_per_comparison_max + 1;
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

function singleBoxPlot(target_axes_handle, abscissa, data, color, label)
    % set some parameters, these should be name-value pair arguments later
    width = 0.8;
    
    % extract data
    data_median = median(data);
    data_quartile_1 = prctile(data, 25);
    data_quartile_3 = prctile(data, 75);
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
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off' ...
      );
    plot(target_axes_handle, abscissa + width*[-0.5 0.5], [data_median data_median], 'color', 'k', 'HandleVisibility', 'off'); % median
    plot(target_axes_handle, [abscissa abscissa], [data_quartile_3 data_upper_adjacent], 'k--', 'HandleVisibility', 'off'); % upper range
    plot(target_axes_handle, [abscissa abscissa], [data_lower_adjacent data_quartile_1], 'k--', 'HandleVisibility', 'off'); % lower range
    plot(target_axes_handle, abscissa+width*[-0.25 0.25], [data_lower_adjacent data_lower_adjacent], 'k-', 'HandleVisibility', 'off'); % max
    plot(target_axes_handle, abscissa+width*[-0.25 0.25], [data_upper_adjacent data_upper_adjacent], 'k-', 'HandleVisibility', 'off'); % min
    plot(target_axes_handle, abscissa * ones(size(outliers)), outliers, '+', 'color', [1; 1; 1] * 0.7, 'HandleVisibility', 'off');
    
    % labels
    xtick = get(target_axes_handle, 'xtick');
    xticklabels = get(target_axes_handle, 'xticklabel');
    xticklabels{xtick == abscissa} = label;
    set(target_axes_handle, 'xticklabel', xticklabels);
end




