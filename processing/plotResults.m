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

% concepts
% condition - a condition is a combination of specific levels for each relevant factor. The code will generate a cell 
%   array "condition_combination_labels", where each column is a factor and each row a condition.

% A comparison is a collection of conditions that are mostly the same, but different in one factor. This code generates 
%   a cell array "comparison_indices", where each entry is an array of condition indices, referring to a row in 
%   condition_combination_labels. Usually there will be one comparison per figure. When plotting episodes of multiple
%   consecutive steps, however, there will be multiple comparisons per figure. To keep track of this, the code generates
%   an array "trajectory_axes_handles" that maps comparisons and variables to axes handles, and a cell array 
%   "abscissae_cell" storing the x-values. For both these arrays, rows = comparisons, columns = variables

% An episode is a list of comparisons
%   (this is a bit of a mess and will be fixed eventually... I promise - HR)

function plotResults(varargin)
    %% parse input
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'subjects', [])
    addParameter(parser, 'dictate_axes', false)
    addParameter(parser, 'show_legend', false)
    addParameter(parser, 'save', false)
    addParameter(parser, 'format', 'tiff')
    addParameter(parser, 'settings', 'plotSettings.txt')
    addParameter(parser, 'spread_method', 'cinv')
    parse(parser, varargin{:})
    subjects = parser.Results.subjects;
    dictate_axes = parser.Results.dictate_axes;
    show_legend = parser.Results.show_legend;
    settings_file = parser.Results.settings;
    spread_method = parser.Results.spread_method;

    % load settings
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
    
    plot_settings_file = '';
    if exist(settings_file, 'file')
        plot_settings_file = settings_file;
    end    
    if exist(['..' filesep settings_file], 'file')
        plot_settings_file = ['..' filesep settings_file];
    end    
    if exist(['..' filesep '..' filesep settings_file], 'file')
        plot_settings_file = ['..' filesep '..' filesep settings_file];
    end
    plot_settings = SettingsCustodian(plot_settings_file);
    show_outliers = plot_settings.get('show_outliers');
    plot_mode = plot_settings.get('plot_mode');
    mark_pushoff = plot_settings.get('mark_pushoff');
    mark_bands = plot_settings.get('mark_bands');
    band_labels = study_settings.get('band_labels');
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');

    %% load data
    % declare variables
    conditions_settings = study_settings.get('conditions');
    condition_to_compare = plot_settings.get('condition_to_compare');
    condition_labels = conditions_settings(:, 1)';
    condition_source_variables = conditions_settings(:, 2)';
    number_of_condition_labels = length(condition_labels);
    
    % load data
    data_folder_list = determineDataStructure(subjects);
    variables_to_plot = plot_settings.get('variables_to_plot');
    number_of_variables_to_plot = size(variables_to_plot, 1);
    if size(variables_to_plot, 2) ~= 7
        disp('Expected entries in variables_to_plot in plot settings file are:')
        disp('variable name, variable_source, variable label, y-axis label, save file string, y-axis lower limit, y-axis upper limit')
        error('Wrong number of columns in variables_to_plot in plot settings file. ');
    end
    condition_data_all = {};
    origin_trial_list_all = [];
    origin_start_time_list_all = [];
    origin_end_time_list_all = [];
    data_all = cell(number_of_variables_to_plot, 1);
    directions = cell(number_of_variables_to_plot, 2);
    
    step_time_data = [];
    pushoff_time_data = [];
    bands_per_stretch = [];
    
    for i_folder = 1 : length(data_folder_list)
        % load data
        data_path = data_folder_list{i_folder};
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
        results_file_name = [data_path filesep 'analysis' filesep makeFileName(date, subject_id, 'results')];
        loaded_data = load(results_file_name);
        number_of_stretches_this_session = length(loaded_data.time_list_session);
        bands_per_stretch_this_session = loaded_data.bands_per_stretch;

        % transform conditions into cell array
        conditions_session = loaded_data.conditions_session;
        condition_array_session = cell(number_of_stretches_this_session, number_of_condition_labels);
        for i_condition = 1 : number_of_condition_labels
            condition_array_session(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
        end
        
        condition_data_all = [condition_data_all; condition_array_session]; %#ok<AGROW>
        origin_trial_list_all = [origin_trial_list_all; loaded_data.origin_trial_list_session]; %#ok<AGROW>
        origin_start_time_list_all = [origin_start_time_list_all; loaded_data.origin_start_time_list_session]; %#ok<AGROW>
        origin_end_time_list_all = [origin_end_time_list_all; loaded_data.origin_end_time_list_session]; %#ok<AGROW>
        
        
        if any(strcmp(variables_to_plot(:, 2), 'stretch'))
            stretch_names_session = loaded_data.stretch_names_session;
            stretch_data_session = loaded_data.stretch_data_session;
            stretch_directions_session = loaded_data.stretch_directions_session;
        end
        if any(strcmp(variables_to_plot(:, 2), 'response'))
            response_names_session = loaded_data.response_names_session;
            response_data_session = loaded_data.response_data_session;
            response_directions_session = loaded_data.response_directions_session;
        end
        if any(strcmp(variables_to_plot(:, 2), 'analysis'))
            analysis_names_session = loaded_data.analysis_names_session;
            analysis_data_session = loaded_data.analysis_data_session;
            analysis_directions_session = loaded_data.analysis_directions_session;
        end
        for i_variable = 1 : number_of_variables_to_plot
            % load and extract data
            this_variable_name = variables_to_plot{i_variable, 1};
            this_variable_source = variables_to_plot{i_variable, 2};
            if strcmp(this_variable_source, 'stretch')
                index_in_saved_data = find(strcmp(stretch_names_session, this_variable_name), 1, 'first');
            end
            if strcmp(this_variable_source, 'response')
                index_in_saved_data = find(strcmp(response_names_session, this_variable_name), 1, 'first');
            end
            if strcmp(this_variable_source, 'analysis')
                index_in_saved_data = find(strcmp(analysis_names_session, this_variable_name), 1, 'first');
            end
            
            if isempty(index_in_saved_data)
                error(['Data not found: ' this_variable_name])
            end
            
            if strcmp(this_variable_source, 'stretch')
                this_variable_data = stretch_data_session{index_in_saved_data};
                this_variable_directions = stretch_directions_session(index_in_saved_data, :);
            end
            if strcmp(this_variable_source, 'response')
                this_variable_data = response_data_session{index_in_saved_data};
                this_variable_directions = response_directions_session(index_in_saved_data, :);
            end
            if strcmp(this_variable_source, 'analysis')
                this_variable_data = analysis_data_session{index_in_saved_data};
                this_variable_directions = analysis_directions_session(index_in_saved_data, :);
            end
            
            if plot_settings.get('convert_to_mm') && (strcmp(this_variable_name,'cop_from_com_x') || strcmp(this_variable_name, 'step_placement_x'))
                this_variable_data = this_variable_data * 1000;
            end
            
            % store
            data_all{i_variable} = [data_all{i_variable} this_variable_data];
            directions(i_variable, :) = this_variable_directions;
        end
        % get time variables
        if any(find(strcmp(loaded_data.stretch_names_session, 'step_time')))
            index_in_saved_data = find(strcmp(loaded_data.stretch_names_session, 'step_time'), 1, 'first');
            this_step_time_data = loaded_data.stretch_data_session{index_in_saved_data};
            step_time_data = [step_time_data this_step_time_data]; %#ok<AGROW>
        end
        if any(find(strcmp(loaded_data.stretch_names_session, 'pushoff_time')))
            index_in_saved_data = find(strcmp(loaded_data.stretch_names_session, 'pushoff_time'), 1, 'first');
            this_pushoff_time_data = loaded_data.stretch_data_session{index_in_saved_data};
            pushoff_time_data = [pushoff_time_data this_pushoff_time_data]; %#ok<AGROW>
        end
        if isempty(bands_per_stretch)
            bands_per_stretch = bands_per_stretch_this_session;
        else
            if bands_per_stretch ~= bands_per_stretch_this_session
               warning('Different sessions have different numbers of bands per stretch') 
            end
        end
    end
    % calculate mean pushoff index
    if mark_pushoff
        pushoff_time_ratio = pushoff_time_data ./ step_time_data;
        mean_pushoff_ratio = mean(pushoff_time_ratio);
        pushoff_index = round(mean_pushoff_ratio * 100);
    end
    
    %% populate condition cell
    labels_to_ignore = plot_settings.get('conditions_to_ignore');
    levels_to_remove = plot_settings.get('levels_to_remove');
    preferred_level_order = plot_settings.get('preferred_level_order');
    [condition_combination_labels, condition_combinations_stimulus, condition_combinations_control] = determineConditionCombinations(condition_data_all, conditions_settings, labels_to_ignore, levels_to_remove);
    condition_combinations_stimulus = sortConditionCombinations(condition_combinations_stimulus, condition_combination_labels, condition_to_compare, preferred_level_order);
    
    %% determine subjects and data folders
    [comparison_indices, conditions_per_comparison_max] = determineComparisons(condition_combinations_stimulus, condition_combination_labels, plot_settings);
    number_of_comparisons = length(comparison_indices);
    if strcmp(plot_mode, 'episodes')
        episodes = determineEpisodes(condition_combinations_stimulus, condition_combination_labels, comparison_indices, plot_settings);
        number_of_episodes = length(episodes);
    end
    
    %% create figures and determine abscissae for each comparison
    comparison_variable_to_axes_index_map = zeros(number_of_comparisons, 1);
    abscissae_cell = cell(number_of_comparisons, number_of_variables_to_plot);
%     comparison_path_to_axes_index_map = zeros(number_of_comparisons, 1);
    
    % time plots
    if strcmp(plot_mode, 'detailed') || strcmp(plot_mode, 'overview')
        % make one figure per comparison and variable
        trajectory_figure_handles = zeros(number_of_comparisons, number_of_variables_to_plot);
        trajectory_axes_handles = zeros(number_of_comparisons, number_of_variables_to_plot);
        pos_text_handles = zeros(number_of_comparisons, number_of_variables_to_plot);
        neg_text_handles = zeros(number_of_comparisons, number_of_variables_to_plot);
        pos_arrow_handles = zeros(number_of_comparisons, number_of_variables_to_plot);
        neg_arrow_handles = zeros(number_of_comparisons, number_of_variables_to_plot);
        step_start_times_cell = cell(number_of_comparisons, number_of_variables_to_plot);
        step_end_times_cell = cell(number_of_comparisons, number_of_variables_to_plot);
        step_pushoff_times_cell = cell(number_of_comparisons, number_of_variables_to_plot);
%         step_stance_foot_cell = cell(number_of_comparisons, number_of_variables_to_plot);
        for i_variable = 1 : number_of_variables_to_plot
            for i_comparison = 1 : number_of_comparisons
                this_comparison = comparison_indices{i_comparison};
                % make figure and axes
                new_figure = figure; new_axes = axes; hold on;
                
                % store handles and determine abscissa data
                trajectory_figure_handles(i_comparison, i_variable) = new_figure;
                trajectory_axes_handles(i_comparison, i_variable) = new_axes;
                comparison_variable_to_axes_index_map(i_comparison) = i_comparison;
                    
                if isDiscreteVariable(i_variable, data_all, bands_per_stretch)
                    % abscissae gives the bin edges here
                    data_to_plot = data_all{i_variable, 1};
                    if dictate_axes
                        lower_bound = str2double(variables_to_plot{i_variable, 6});
                        upper_bound = str2double(variables_to_plot{i_variable, 7});
                    else
                        lower_bound = min(data_to_plot);
                        upper_bound = max(data_to_plot);
                    end
                    if strcmp(plot_mode, 'detailed')
%                         abscissae_cell{i_comparison, i_variable}(i_condition, :) = linspace(lower_bound, upper_bound, plot_settings.get('number_of_bins_in_histogram'));
                        abscissae_cell{i_comparison, i_variable} = linspace(lower_bound, upper_bound, plot_settings.get('number_of_bins_in_histogram'));
                    end
                    if strcmp(plot_mode, 'overview')
                        this_comparison = comparison_indices{i_comparison};
                        number_of_entries = length(this_comparison);
                        if plot_settings.get('plot_control')
                            number_of_entries = number_of_entries + 1;
                        end
                        
                        if plot_settings.get('group_bands_within_conditions')
                            gap_between_conditions = 1;
                            
                            abscissae_stimulus = repmat((1 : bands_per_stretch)', 1, length(this_comparison));
                            shifter = (0:number_of_entries-1) * (bands_per_stretch + gap_between_conditions);
                            abscissae_stimulus = abscissae_stimulus + repmat(shifter, bands_per_stretch, 1);
                            if plot_settings.get('merge_bands')
                                abscissae_stimulus = abscissae_stimulus(1, :);
                            end
                            
                        else
                            gap_between_bands = 1;

%                             ab
                            abscissae_stimulus = repmat((1 : number_of_entries), bands_per_stretch, 1);
                            shifter = (0:bands_per_stretch-1)' * (conditions_per_comparison_max + gap_between_bands);
                            abscissae_stimulus = abscissae_stimulus + repmat(shifter, 1, conditions_per_comparison_max);
                            if plot_settings.get('merge_bands')
                                abscissae_stimulus = abscissae_stimulus(1, :);
                            end
                        end
                        abscissae_cell{i_comparison, i_variable} = abscissae_stimulus;
                        
                    end
                end
                if isContinuousVariable(i_variable, data_all, bands_per_stretch)
                    % scale abscissae
                    conditions_this_comparison = comparison_indices{i_comparison};
                    step_time_means_this_comparison = zeros(bands_per_stretch, size(conditions_this_comparison, 2));
                    for i_condition = 1 : length(conditions_this_comparison)
                        this_condition_combination = condition_combinations_stimulus(conditions_this_comparison(i_condition), :);
                        this_condition_indicator = getConditionIndicator(this_condition_combination, condition_combination_labels, condition_data_all, condition_labels);
                        step_time_data_this_condition = step_time_data(:, this_condition_indicator);
                        step_time_means_this_comparison(:, i_condition) = mean(step_time_data_this_condition, 2);
                    end
                    for i_condition = 1 : length(conditions_this_comparison)
                        if strcmp(plot_settings.get('time_plot_style'), 'scaled_to_comparison_mean')
                            band_scales = mean(step_time_means_this_comparison, 2);
                        elseif strcmp(plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                            band_scales = step_time_means_this_comparison(:, i_condition);
                        else
                            band_scales = ones(bands_per_stretch, 1) * (number_of_time_steps_normalized-1);
                        end
                        [abscissa_scaled, band_limits] = createScaledAbscissa(band_scales, number_of_time_steps_normalized);
                        
                        if strcmp(plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                            time_plot_band_anchor_index = plot_settings.get('time_plot_band_anchor');
                            time_plot_band_anchor_time = band_limits(time_plot_band_anchor_index);
                            abscissa_scaled = abscissa_scaled - time_plot_band_anchor_time;
                        end
                        
                        abscissae_cell{i_comparison, i_variable}(i_condition, :) = abscissa_scaled;
                        

                    end                    
                end
                
                % set axes properties
                if dictate_axes && ~(strcmp(plot_mode, 'detailed') && isDiscreteVariable(i_variable, data_all, bands_per_stretch))
    %                 set(gca, 'xlim', [time_normalized(1), time_normalized(end)]);
                    set(gca, 'ylim', [str2double(variables_to_plot{i_variable, 6}), str2double(variables_to_plot{i_variable, 7})]);
                end
                if isDiscreteVariable(i_variable, data_all, bands_per_stretch) && strcmp(plot_mode, 'overview')
%                     xtick = abscissae_cell{i_comparison, i_variable}{2};
%                     if plot_settings.get('plot_control')
%                         xtick = [abscissae_cell{i_comparison, i_variable}{1} xtick]; %#ok<AGROW>
%                     end
                    xtick = sort(reshape(abscissae_cell{i_comparison, i_variable}, 1, numel(abscissae_cell{i_comparison, i_variable})));
                    set(gca, 'xlim', [-0.5 + min(xtick) 0.5 + max(xtick(end))]);
                    set(gca, 'xtick', xtick);
                end
                
                % set axis labels
                if isContinuousVariable(i_variable, data_all, bands_per_stretch)
                    if strcmp(plot_settings.get('time_plot_style'), 'scaled_to_comparison_mean') || strcmp(plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                        xlabel('normalized time (s)');
                    else
                        xlabel('normalized time (%)');
                    end
                end
                ylabel(variables_to_plot{i_variable, 4});
                
                % add text labels
                pos_text_handles(i_comparison, i_variable) = ...
                    text ...
                      ( ...
                        0, ...
                        0, ...
                        directions{i_variable, 1}, ...
                        'rotation', 90, ...
                        'Fontsize', 24, ...
                        'horizontalalignment', 'right', ...
                        'parent', new_axes ...
                      );                pos_arrow_handles(i_comparison, i_variable) = ...
                    text ...
                      ( ...
                        0, ...
                        0, ...
                        ' $\rightarrow$', ...
                        'rotation', 90, ...
                        'Fontsize', 36, ...
                        'horizontalalignment', 'right', ...
                        'interpreter', 'LaTeX', ...
                        'parent', new_axes ...
                      );
                neg_text_handles(i_comparison, i_variable) = ...
                    text ...
                      ( ...
                        0, ...
                        0, ...
                        directions{i_variable, 2}, ...
                        'rotation', 90, ...
                        'Fontsize', 24, ...
                        'horizontalalignment', 'left', ...
                        'parent', new_axes...
                      );
                neg_arrow_handles(i_comparison, i_variable) = ...
                    text ...
                      ( ...
                        0, ...
                        0, ...
                        '$\leftarrow$ ', ...
                        'rotation', 90, ...
                        'Fontsize', 36, ...
                        'horizontalalignment', 'left', ...
                        'interpreter', 'LaTeX', ...
                        'parent', new_axes ...
                      );                
                % determine title
                title_string = variables_to_plot{i_variable, 3};
                filename_string = variables_to_plot{i_variable, 5};
                
                representative_condition = condition_combinations_stimulus(this_comparison(1), :);
                
                for i_label = 1 : length(representative_condition)
                    if ~(strcmp(condition_combination_labels{i_label}, condition_to_compare))
                        this_string = strrep(representative_condition{i_label}, '_', '');
                        filename_string = [filename_string '_' this_string]; %#ok<AGROW>
                        title_string = [title_string ' - ' this_string]; %#ok<AGROW>
                    end
                end
                
                

                title(title_string); set(gca, 'Fontsize', 12)
                set(gcf, 'UserData', filename_string)
                
                
            end
        end
        
        % set x-limits
        for i_variable = 1 : number_of_variables_to_plot
            if isContinuousVariable(i_variable, data_all, bands_per_stretch)
                for i_comparison = 1 : number_of_comparisons
                    target_abscissae = abscissae_cell{i_comparison, i_variable};
                    xlim = [min(target_abscissae(:, 1)) max(target_abscissae(:, end))];
                    % set x-limits accordingly
                    set(trajectory_axes_handles(i_comparison, i_variable), 'xlim', xlim);
                end
            end
        end

        % determine stance start and end times and stance foot
        for i_variable = 1 : number_of_variables_to_plot
            if isContinuousVariable(i_variable, data_all, bands_per_stretch)
                for i_comparison = 1 : number_of_comparisons
                    
                    step_abscissa = abscissae_cell{i_comparison, i_variable};
                    step_start_times_cell{i_comparison, i_variable} = step_abscissa(1, 1);
                    step_end_times_cell{i_comparison, i_variable} = step_abscissa(1, end);
                    
                    if mark_pushoff
                        step_pushoff_times_cell{i_comparison, i_variable} = step_abscissa(1, pushoff_index);
                    end
                end
            end
        end        
    end
    if strcmp(plot_mode, 'episodes')
        abscissae_cell_unscaled = cell(size(abscissae_cell));
        % make one figure per episode and variable
        trajectory_figure_handles = zeros(number_of_episodes, number_of_variables_to_plot);
        trajectory_axes_handles = zeros(number_of_episodes, number_of_variables_to_plot);
        pos_text_handles = zeros(number_of_episodes, number_of_variables_to_plot);
        neg_text_handles = zeros(number_of_episodes, number_of_variables_to_plot);
        pos_arrow_handles = zeros(number_of_episodes, number_of_variables_to_plot);
        neg_arrow_handles = zeros(number_of_episodes, number_of_variables_to_plot);
        step_start_times_cell = cell(number_of_episodes, number_of_variables_to_plot);
        step_end_times_cell = cell(number_of_episodes, number_of_variables_to_plot);
        step_pushoff_times_cell = cell(number_of_episodes, number_of_variables_to_plot);
%         step_stance_foot_cell = cell(number_of_episodes, number_of_variables_to_plot);

        for i_variable = 1 : number_of_variables_to_plot
            for i_episode = 1 : number_of_episodes
                % make figure and axes and store handles
                new_figure = figure; new_axes = axes; hold on;
                trajectory_figure_handles(i_episode, i_variable) = new_figure;
                trajectory_axes_handles(i_episode, i_variable) = new_axes;
                this_episode = episodes{i_episode};

                % store handles and determine abscissa data for all comparisons in this episode
                xtick = [];
                for i_comparison = 1 : size(this_episode, 2)
                    comparison_variable_to_axes_index_map(this_episode(i_comparison)) = i_episode;
                    
                    % determine which step this is
                    this_comparison = this_episode(i_comparison);
                    conditions_this_comparison = comparison_indices{this_comparison};
                    example_condition_index = conditions_this_comparison(1);
                    condition_identifier = condition_combinations_stimulus(example_condition_index, :);
                    gap_between_steps = 1;
                    if strcmp(condition_identifier{strcmp(condition_combination_labels, 'index')}, 'ONE')
                        step_index = 1;
                    elseif strcmp(condition_identifier{strcmp(condition_combination_labels, 'index')}, 'TWO')
                        step_index = 2;
                    elseif strcmp(condition_identifier{strcmp(condition_combination_labels, 'index')}, 'THREE')
                        step_index = 3;
                    elseif strcmp(condition_identifier{strcmp(condition_combination_labels, 'index')}, 'FOUR')
                        step_index = 4;
                    end
                    if isDiscreteVariable(i_variable, data_all, bands_per_stretch)
                        this_comparison = comparison_indices{i_comparison};
                        abscissae_control = (conditions_per_comparison_max + gap_between_steps) * step_index;
                        abscissae_stimulus = (1 : length(this_comparison)) + (conditions_per_comparison_max + gap_between_steps) * step_index;
                        abscissae = {abscissae_control, abscissae_stimulus};
                        abscissae_cell{this_episode(i_comparison), i_variable} = abscissae;
                        
%                         if ~isempty(conditions_control) && plot_settings.get('plot_control') && strcmp(data_source, 'stretch')
                        if plot_settings.get('plot_control')
                            xtick = [xtick abscissae{1}]; %#ok<AGROW>
                        end
                        xtick = [xtick abscissae{2}]; %#ok<AGROW>



                    end
                    if isContinuousVariable(i_variable, data_all, bands_per_stretch)
                        abscissae_cell_unscaled{this_episode(i_comparison), i_variable} = (linspace(0, 100, study_settings.get('number_of_time_steps_normalized')));
                    end
                end
                
                % set axes properties
                if dictate_axes
                    set(gca, 'ylim', [str2double(variables_to_plot{i_variable, 6}), str2double(variables_to_plot{i_variable, 7})]);
                end
                if isDiscreteVariable(i_variable, data_all, bands_per_stretch)
                    set(gca, 'xlim', [-0.5 + min(xtick) 0.5 + max(xtick(end))]);
                    set(gca, 'xtick', xtick);
                    set(gca, 'XTickLabelRotation', 60);
                end

                % set axis labels
                if isContinuousVariable(i_variable, data_all, bands_per_stretch)
                    if strcmp(plot_settings.get('time_plot_style'), 'scaled_to_comparison_mean') || strcmp(plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                        xlabel('normalized time (s)');
                    else
                        xlabel('normalized time (%)');
                    end
                end
                ylabel(variables_to_plot{i_variable, 4});
                
                % add text labels
                pos_text_handles(i_episode, i_variable) = ...
                    text ...
                      ( ...
                        0, ...
                        0, ...
                        directions{i_variable, 1}, ...
                        'rotation', 90, ...
                        'Fontsize', 18, ...
                        'horizontalalignment', 'right', ...
                        'parent', new_axes ...
                      );
                pos_arrow_handles(i_episode, i_variable) = ...
                    text ...
                      ( ...
                        0, ...
                        0, ...
                        ' $\rightarrow$', ...
                        'rotation', 90, ...
                        'Fontsize', 36, ...
                        'horizontalalignment', 'right', ...
                        'interpreter', 'LaTeX', ...
                        'parent', new_axes ...
                      );
                neg_text_handles(i_episode, i_variable) = ...
                    text ...
                      ( ...
                        0, ...
                        0, ...
                        directions{i_variable, 2}, ...
                        'rotation', 90, ...
                        'Fontsize', 18, ...
                        'horizontalalignment', 'left', ...
                        'parent', new_axes...
                      );
                neg_arrow_handles(i_episode, i_variable) = ...
                    text ...
                      ( ...
                        0, ...
                        0, ...
                        '$\leftarrow$ ', ...
                        'rotation', 90, ...
                        'Fontsize', 36, ...
                        'horizontalalignment', 'left', ...
                        'interpreter', 'LaTeX', ...
                        'parent', new_axes...
                      );
                
                % determine title and filename
                first_comparison = this_episode(1);
                conditions_first_comparison = comparison_indices{first_comparison};
                example_condition_index = conditions_first_comparison(1);
                condition_identifier = condition_combinations_stimulus(example_condition_index, :);
                title_string = variables_to_plot{i_variable, 3};
                filename_string = variables_to_plot{i_variable, 5};
                
                % new
                condition_to_compare = plot_settings.get('condition_to_compare');
                relevant_labels = find(~(strcmp(condition_combination_labels, 'index') | strcmp(condition_combination_labels, condition_to_compare)));
                for i_label = relevant_labels
                    title_string = [title_string ' - ' strrep(condition_identifier{1, i_label}, '_', ' ')]; %#ok<AGROW>
                    filename_string = [filename_string '_' strrep(condition_identifier{1, i_label}, '_', '')]; %#ok<AGROW>
                end
                
                % old
%                 for i_label = 1 : length(condition_labels)
%                     if (i_label ~= plot_settings.get('comparison_to_make')) ...
%                         && (i_label ~= 1) ...
%                         && (i_label ~= 4) ...
%                         && (i_label ~= 5) ...
%                         && (i_label ~= 6)
%                         this_condition_label = strrep(condition_identifier{1, i_label}, '_', ' ');
%                         if ~strcmp(this_condition_label, 'N/A')
%                             title_string = [title_string ' - ' this_condition_label]; %#ok<AGROW>
%                             filename_string = [filename_string '_' this_condition_label];
%                         end
%                     end
%                 end
%                 if strcmp(condition_identifier{1}, 'STANCE_RIGHT')
%                     title_string = [title_string ' - first step stance leg RIGHT'];
%                     filename_string = [filename_string '_stanceR'];
%                 end
%                 if strcmp(condition_identifier{1}, 'STANCE_LEFT')
%                     title_string = [title_string ' - first step stance leg LEFT'];
%                     filename_string = [filename_string '_stanceL'];
%                 end
                
                title(title_string); set(gca, 'Fontsize', 12)
                set(gcf, 'UserData', filename_string)
                
                
            end
        end
        
        % calculate average step times and scale abscissa
        for i_variable = 1 : number_of_variables_to_plot
            for i_episode = 1 : number_of_episodes
                this_episode = episodes{i_episode};
                for i_comparison = 1 : length(this_episode)
                    
                    % determine which step this is
                    this_comparison = this_episode(i_comparison);
                    conditions_this_comparison = comparison_indices{this_comparison};
                    step_time_means_this_comparison = zeros(size(conditions_this_comparison));
                    for i_condition = 1 : length(conditions_this_comparison)
                        % find correct condition indicator
                        condition_identifier = condition_combinations_stimulus(conditions_this_comparison(i_condition), :);
                        this_condition_indicator = getConditionIndicator(condition_identifier, condition_combination_labels, condition_data_all, condition_labels);
                        
%                         stance_foot_indicator = strcmp(condition_stance_foot_list_all, condition_identifier{1});
%                         perturbation_indicator = strcmp(condition_perturbation_list_all, condition_identifier{2});
%                         delay_indicator = strcmp(condition_delay_list_all, condition_identifier{3});
%                         index_indicator = strcmp(condition_index_list_all, condition_identifier{4});
%                         experimental_indicator = strcmp(condition_experimental_list_all, condition_identifier{5});
%                         stimulus_indicator = strcmp(condition_stimulus_list_all, condition_identifier{6});
%                         day_indicator = strcmp(condition_day_list_all, condition_identifier{7});
%                         this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
                        
                        
                        if strcmp(plot_settings.get('time_plot_style'), 'scaled_to_comparison_mean')
                            % calculate average step time
                            step_time_data_this_condition = step_time_data(:, this_condition_indicator);
                            step_time_means_this_comparison(i_condition) = mean(step_time_data_this_condition);
                        end
                        
                    end
                    
                    % scale abscissa
                    if isContinuousVariable(i_variable, data_all, bands_per_stretch)
                        abscissa_unscaled = abscissae_cell_unscaled{this_episode(i_comparison), i_variable};
                        for i_condition = 1 : length(conditions_this_comparison)
                            if strcmp(plot_settings.get('time_plot_style'), 'scaled_to_comparison_mean')
                                abscissa_scaled = abscissa_unscaled * mean(step_time_means_this_comparison) / 100;
                                abscissae_cell{this_episode(i_comparison), i_variable}(i_condition, :) = abscissa_scaled;
                            elseif strcmp(plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
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
                if isContinuousVariable(i_variable, data_all, bands_per_stretch)
                    for i_episode = 1 : number_of_episodes
                        this_episode = episodes{i_episode};
                        for i_comparison = 1 : length(this_episode)

                            % determine which step this is
                            this_comparison = this_episode(i_comparison);
                            conditions_this_comparison = comparison_indices{this_comparison};
                            for i_condition = 1 : length(conditions_this_comparison)
                                % determine step index
                                condition_identifier = condition_combinations_stimulus(conditions_this_comparison(i_condition), :);
                                if strcmp(condition_identifier{strcmp(condition_combination_labels, 'index')}, 'ONE')
                                    step_index = 1;
                                elseif strcmp(condition_identifier{strcmp(condition_combination_labels, 'index')}, 'TWO')
                                    step_index = 2;
                                    previous_step_label = 'ONE';
                                elseif strcmp(condition_identifier{strcmp(condition_combination_labels, 'index')}, 'THREE')
                                    step_index = 3;
                                    previous_step_label = 'TWO';
                                elseif strcmp(condition_identifier{strcmp(condition_combination_labels, 'index')}, 'FOUR')
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
                                            candidate_condition_identifier = condition_combinations_stimulus(conditions_candidate_comparison(j_condition), :);
% got up to here: TODO: figure out how to replace the explicit references to columns here... the goal is to find the 
% condition index for the previous step in the same condition. Maybe I can use the episode map that I made in determineEpisodes

                                            if strcmp(candidate_condition_identifier{strcmp(condition_combination_labels, 'index')}, previous_step_label)
                                                % this condition has the right index, now check if it also matches the factor levels
                                                relevant_labels = ~strcmp(condition_combination_labels, 'index');
                                                match = 1;
                                                for i_label = find(relevant_labels)
                                                    if ~strcmp(condition_identifier{i_label}, candidate_condition_identifier{i_label})
                                                        match = 0;
                                                    end
                                                end
                                                if match
                                                    previous_step_comparison_index = j_comparison;
                                                    previous_step_condition_index = j_condition;
                                                end
                                            end
                                                
%                                             if strcmp(condition_identifier{2}, candidate_condition_identifier{2}) ...
%                                             && strcmp(condition_identifier{3}, candidate_condition_identifier{3}) ...
%                                             && strcmp(previous_step_label, candidate_condition_identifier{4}) ...
%                                             && strcmp(condition_identifier{5}, candidate_condition_identifier{5}) ...
%                                             && strcmp(condition_identifier{6}, candidate_condition_identifier{6}) ...
%                                             && strcmp(condition_identifier{7}, candidate_condition_identifier{7})
%                                                 previous_step_comparison_index = j_comparison;
%                                                 previous_step_condition_index = j_condition;
%                                             end
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
            if isContinuousVariable(i_variable, data_all, bands_per_stretch)
                for i_episode = 1 : number_of_episodes
                    % determine time window to show
                    this_episode = episodes{i_episode};
                    first_comparison_in_episode_index = this_episode(1);
                    first_step_abscissae = abscissae_cell{first_comparison_in_episode_index, i_variable};
                    episode_start_time = first_step_abscissae(1, 1);
                    last_comparison_in_episode_index = this_episode(end);
                    last_step_abscissae = abscissae_cell{last_comparison_in_episode_index, i_variable};
                    episode_end_time = last_step_abscissae(1, end);

                    % set x-limits accordingly
                    set(trajectory_axes_handles(i_episode, i_variable), 'xlim', [episode_start_time episode_end_time]);
                end
            end
        end
        
        % determine stance start and end times and stance foot
        for i_variable = 1 : number_of_variables_to_plot
            if isContinuousVariable(i_variable, data_all, bands_per_stretch)
                for i_episode = 1 : number_of_episodes
                    % get start times and end times for the steps
                    this_episode = episodes{i_episode};
                    step_start_times = zeros(1, length(this_episode));
                    step_end_times = zeros(1, length(this_episode));
                    step_pushoff_times = zeros(1, length(this_episode));
%                     step_stance_foot = zeros(1, length(this_episode));
                    
                    for i_comparison = 1 : length(this_episode)
                        this_comparison = this_episode(i_comparison);
                        step_abscissa = abscissae_cell{this_comparison, i_variable};
                        step_start_times(i_comparison) = step_abscissa(1, 1);
                        if mark_pushoff
                            step_pushoff_times(i_comparison) = step_abscissa(1, pushoff_index);
                        end
                        step_end_times(i_comparison) = step_abscissa(1, end);

                        % determine stance foot
%                         conditions_this_comparison = comparison_indices{this_comparison};
%                         example_condition_index = 1;
%                         condition_identifier = condition_combinations_stimulus(conditions_this_comparison(example_condition_index), :);
%                         if strcmp(condition_identifier{1}, 'STANCE_BOTH')
%                             step_stance_foot(i_comparison) = 0;
%                         end
%                         if strcmp(condition_identifier{1}, 'STANCE_LEFT')
%                             step_stance_foot(i_comparison) = 1;
%                         end
%                         if strcmp(condition_identifier{1}, 'STANCE_RIGHT')
%                             step_stance_foot(i_comparison) = 2;
%                         end
                        
                    end
                    
                    step_start_times_cell{i_episode, i_variable} = step_start_times;
                    step_end_times_cell{i_episode, i_variable} = step_end_times;
                    step_pushoff_times_cell{i_episode, i_variable} = step_pushoff_times;
%                     step_stance_foot_cell{i_episode, i_variable} = step_stance_foot;
                end
            end
        end
        
    end
    
    % path plots - removed for now by if false, check back later
    if false
        if strcmp(plot_mode, 'detailed') || strcmp(plot_mode, 'overview') %#ok<UNRCH>
            % make one figure per comparison and variable
            path_figure_handles = zeros(number_of_comparisons, number_of_variables_to_plot);
            path_axes_handles = zeros(number_of_comparisons, number_of_variables_to_plot);
            for i_path = 1 : number_of_paths_to_plot
                for i_comparison = 1 : number_of_comparisons
                    % make figure and axes
                    new_figure = figure; new_axes = axes; hold on;

                    % store handles and determine abscissa data
                    path_figure_handles(i_comparison, i_path) = new_figure;
                    path_axes_handles(i_comparison, i_path) = new_axes;
                    comparison_path_to_axes_index_map(i_comparison) = i_comparison;

                    % determine title
                    title_string = paths_to_plot{i_path, 2};
                    filename_string = paths_to_plot{i_path, 6};
                    for i_label = 1 : length(study_settings.get('condition_labels'))
                        if (i_label ~= plot_settings.get('comparison_to_make')) ...
                            && (i_label ~= 1) ...
                            && (i_label ~= 3) ...
                            && (i_label ~= 5) ...
                            && (i_label ~= 6) ...
                            && (i_label ~= 7)
                            this_condition_label = strrep(conditions_to_plot{comparison_indices{i_comparison}(1), i_label}, '_', ' ');
                            if i_label ~= plot_settings.get('comparison_to_make')
                                title_string = [title_string ' - ' this_condition_label];
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
    end
    
    %% plot data
    colors_comparison = plot_settings.get('colors_comparison');
    colors_bands = plot_settings.get('colors_bands');
    for i_variable = 1 : number_of_variables_to_plot
        data_to_plot = data_all{i_variable, 1};
        for i_comparison = 1 : length(comparison_indices)
            % find correct condition indicator for control
            conditions_this_comparison = comparison_indices{i_comparison};
            top_level_plots = [];
            target_axes_handle = trajectory_axes_handles(comparison_variable_to_axes_index_map(i_comparison), i_variable);
            
            % plot control
            if plot_settings.get('plot_control')
                % determine which control condition applies here
                representant_condition_index = conditions_this_comparison(1);
                this_condition = condition_combinations_control(representant_condition_index, :);
                
                this_condition_indicator = getConditionIndicator(this_condition, condition_combination_labels, condition_data_all, condition_labels);
                data_to_plot_this_condition = data_to_plot(:, this_condition_indicator);
                origin_indices = find(this_condition_indicator);
                
                if ~isempty(data_to_plot_this_condition)
                    if isDiscreteVariable(i_variable, data_all, bands_per_stretch)
                        if plot_settings.get('merge_bands')
                            data_to_plot_this_condition = reshape(data_to_plot_this_condition, 1, numel(data_to_plot_this_condition));
                        end
                        % ----------------------------------------------------------------------------------------------
                        % start of fix
                        % ----------------------------------------------------------------------------------------------
                        for i_band = 1 : size(data_to_plot_this_condition, 1)
                            if ~isempty(band_labels)
                                label_string_this_band = ['control -' band_labels{i_band}];
                            else
                                label_string_this_band = 'control';
                            end
                            if strcmp(plot_mode, 'episodes')
                                % TODO: copied over from stimulus, fix this later
                                this_cell = abscissae_cell{i_comparison, i_variable};
                                target_abscissa = this_cell{1};
                            else
                                target_abscissa = abscissae_cell{i_comparison, i_variable}(i_band, end);
                            end
                            data_to_plot_this_band = data_to_plot_this_condition(i_band, :);
                            if strcmp(plot_mode, 'detailed')
                                % TODO: copied over from stimulus, fix this later
%                                 histogram ...
%                                   ( ...
%                                     target_axes_handle, ...
%                                     data_to_plot_this_band, ...
%                                     target_abscissa, ...
%                                     'edgecolor', colors_comparison(i_condition, :), ...
%                                     'facecolor', lightenColor(colors_comparison(i_condition, :), 0.5), ...
%                                     'DisplayName', label_string ...
%                                   );
                            end
                            if strcmp(plot_mode, 'overview') || strcmp(plot_mode, 'episodes')
                                if ~any(isnan(data_to_plot_this_band))
                                    if plot_settings.get('group_bands_within_conditions')
                                        this_color = colors_bands(i_band, :);
                                    else
                                        this_color = plot_settings.get('color_control');
                                    end

                                    if strcmp(plot_settings.get('discrete_data_plot_style'), 'box')
                                        singleBoxPlot ...
                                          ( ...
                                            target_axes_handle, ...
                                            target_abscissa, ...
                                            data_to_plot_this_band, ...
                                            this_color, ...
                                            label_string_this_band, ...
                                            show_outliers ...
                                          )
                                    end
                                    if strcmp(plot_settings.get('discrete_data_plot_style'), 'bar')
                                        singleBarPlot ...
                                          ( ...
                                            target_axes_handle, ...
                                            target_abscissa, ...
                                            data_to_plot_this_band, ...
                                            this_color, ...
                                            label_string_this_band ...
                                          )
                                    end
                                    if strcmp(plot_settings.get('discrete_data_plot_style'), 'violin')
                                        singleViolinPlot ...
                                          ( ...
                                            data_to_plot_this_band, ...
                                            'axes', target_axes_handle, ...
                                            'abscissa', target_abscissa, ...
                                            'facecolor', this_color, ...
                                            'plot_mean', false, ...
                                            'plot_median', true, ...
                                            'mediancolor', [0 0 0], ...
                                            'show_outliers', show_outliers, ...
                                            'xlabel', label_string_this_band ...
                                          );
                                    end
                                end
                            end
                        end
                        
                        
                        % ----------------------------------------------------------------------------------------------
                        % HR: before fixing this, leaving this around because this hasn't been tested thoroughly
                        % ----------------------------------------------------------------------------------------------
                        
%                         target_abscissa = abscissae_cell{i_comparison, i_variable};
%                         if strcmp(plot_mode, 'detailed')
%                             histogram ...
%                               ( ...
%                                 target_axes_handle, ...
%                                 data_to_plot_this_condition, ...
%                                 'binEdges', target_abscissa, ...
%                                 'edgecolor', plot_settings.get('color_control'), ...
%                                 'facecolor', lightenColor(plot_settings.get('color_control'), 0.5), ...
%                                 'DisplayName', 'CONTROL' ...
%                               );
%                         end
%                         if strcmp(plot_mode, 'overview') || strcmp(plot_mode, 'episodes')
%                             if strcmp(plot_settings.get('discrete_data_plot_style'), 'box')
%                                 singleBoxPlot ...
%                                   ( ...
%                                     target_axes_handle, ...
%                                     target_abscissa{1}, ...
%                                     data_to_plot_this_condition, ...
%                                     plot_settings.get('color_control'), ...
%                                     'CONTROL', ...
%                                     show_outliers ...
%                                   )
%                             end
%                             if strcmp(plot_settings.get('discrete_data_plot_style'), 'bar')
%                                singleBarPlot ...
%                                    ( ...
%                                      target_axes_handle, ...
%                                      target_abscissa{1}, ...
%                                      data_to_plot_this_condition, ...
%                                      plot_settings.get('color_control'), ...
%                                      'CONTROL' ...
%                                    ) 
%                             end
%                             if strcmp(plot_settings.get('discrete_data_plot_style'), 'violin')
%                                 singleViolinPlot ...
%                                   ( ...
%                                     data_to_plot_this_condition, ...
%                                     'axes', target_axes_handle, ...
%                                     'abscissa', target_abscissa{1}, ...
%                                     'facecolor', plot_settings.get('color_control'), ...
%                                     'plot_mean', false, ...
%                                     'plot_median', true, ...
%                                     'mediancolor', [0 0 0], ...
%                                     'show_outliers', show_outliers, ...
%                                     'xlabel', 'CONTROL' ...
%                                   );
%                             end
%                         end
                        % ----------------------------------------------------------------------------------------------
                        % end
                        % ----------------------------------------------------------------------------------------------
                    end
                    if isContinuousVariable(i_variable, data_all, bands_per_stretch)
                        target_abscissa = abscissae_cell{i_comparison, i_variable}(i_condition, :);
                        if strcmp(plot_mode, 'detailed')
                            % individual trajectories
                            for i_stretch = 1 : size(data_to_plot_this_condition, 2)
                                origin_index_data = ones(size(target_abscissa)) * origin_indices(i_stretch);
%                                 plot ...
%                                   ( ...
%                                     target_axes_handle, ...
%                                     target_abscissa, ...
%                                     data_to_plot_this_condition, ...
                                plot3 ...
                                  ( ...
                                    target_axes_handle, ...
                                    target_abscissa, ...
                                    data_to_plot_this_condition(:, i_stretch), ...
                                    origin_index_data, ... %origin_trial_data, ...
                                    'HandleVisibility', 'off', ...
                                    'color', lightenColor(plot_settings.get('color_control'), 0.5) ...
                                  );
                            end
                            % condition average
                            control_mean_plot = plot ...
                              ( ...
                                target_axes_handle, ...
                                target_abscissa, ...
                                mean(data_to_plot_this_condition, 2), ...
                                'DisplayName', 'CONTROL', ...
                                'linewidth', 5, ...
                                'color', plot_settings.get('color_control') ...
                              );
                            top_level_plots = [top_level_plots control_mean_plot]; %#ok<AGROW>
                        end
                        if strcmp(plot_mode, 'overview') || strcmp(plot_mode, 'episodes')
                            plot_handles = shadedErrorBar ...
                              ( ...
                                target_abscissa, ...
                                mean(data_to_plot_this_condition, 2), ...
                                spread(data_to_plot_this_condition, spread_method), ...
                                { ...
                                  'color', plot_settings.get('color_control'), ...
                                  'linewidth', 6 ...
                                }, ...
                                1, ...
                                target_axes_handle ...
                              );
                            set(plot_handles.edge, 'HandleVisibility', 'off');
                            set(plot_handles.patch, 'HandleVisibility', 'off');
                            set(plot_handles.mainLine, 'DisplayName', 'CONTROL');

                            % can't make this work for both episodes and single, revisit this later
%                             if control_already_labeled
                                set(plot_handles.mainLine, 'HandleVisibility', 'off')
%                             else
%                                 control_already_labeled = true;
%                             end
                            
%                             if strcmp(plot_mode, 'episodes') && ~strcmp(this_condition{strcmp(condition_combination_labels, 'index')}, 'ONE')
%                                 set(plot_handles.mainLine, 'HandleVisibility', 'off');
%                             end

                            top_level_plots = [top_level_plots plot_handles.mainLine]; %#ok<AGROW>
                        end
                    end
                end
            end
            
            % plot stimulus
            for i_condition = 1 : length(conditions_this_comparison)
                this_condition_index = conditions_this_comparison(i_condition);
                this_condition = condition_combinations_stimulus(this_condition_index, :);
                label_string = strrep(this_condition{strcmp(condition_combination_labels, condition_to_compare)}, '_', ' ');
                this_condition_indicator = getConditionIndicator(this_condition, condition_combination_labels, condition_data_all, condition_labels);
                data_to_plot_this_condition = data_to_plot(:, this_condition_indicator);
                
                origin_indices = find(this_condition_indicator);
                if isDiscreteVariable(i_variable, data_all, bands_per_stretch)
                    if plot_settings.get('merge_bands')
                        data_to_plot_this_condition = reshape(data_to_plot_this_condition, 1, numel(data_to_plot_this_condition));
                    end
                    for i_band = 1 : size(data_to_plot_this_condition, 1)
                        if ~isempty(band_labels)
                            label_string_this_band = [label_string '-' band_labels{i_band}];
                        else
                            label_string_this_band = label_string;
                        end
                        if strcmp(plot_mode, 'episodes')
                            this_cell = abscissae_cell{i_comparison, i_variable};
                            target_abscissa = this_cell{2}(i_band, i_condition);
                        else
                            target_abscissa = abscissae_cell{i_comparison, i_variable}(i_band, i_condition);
                        end
                        data_to_plot_this_band = data_to_plot_this_condition(i_band, :);
                        if strcmp(plot_mode, 'detailed')
                            histogram ...
                              ( ...
                                target_axes_handle, ...
                                data_to_plot_this_band, ...
                                target_abscissa, ...
                                'edgecolor', colors_comparison(i_condition, :), ...
                                'facecolor', lightenColor(colors_comparison(i_condition, :), 0.5), ...
                                'DisplayName', label_string ...
                              );
                        end
                        if strcmp(plot_mode, 'overview') || strcmp(plot_mode, 'episodes')
                            if ~any(isnan(data_to_plot_this_band))
                                if plot_settings.get('group_bands_within_conditions')
                                    this_color = colors_bands(i_band, :);
                                else
                                    this_color = colors_comparison(i_condition, :);
                                end
                                
                                if strcmp(plot_settings.get('discrete_data_plot_style'), 'box')
                                    singleBoxPlot ...
                                      ( ...
                                        target_axes_handle, ...
                                        target_abscissa, ...
                                        data_to_plot_this_band, ...
                                        this_color, ...
                                        label_string_this_band, ...
                                        show_outliers ...
                                      )
                                end
                                if strcmp(plot_settings.get('discrete_data_plot_style'), 'bar')
                                    singleBarPlot ...
                                      ( ...
                                        target_axes_handle, ...
                                        target_abscissa, ...
                                        data_to_plot_this_band, ...
                                        this_color, ...
                                        label_string_this_band ...
                                      )
                                end
                                if strcmp(plot_settings.get('discrete_data_plot_style'), 'violin')
                                    singleViolinPlot ...
                                      ( ...
                                        data_to_plot_this_band, ...
                                        'axes', target_axes_handle, ...
                                        'abscissa', target_abscissa, ...
                                        'facecolor', this_color, ...
                                        'plot_mean', false, ...
                                        'plot_median', true, ...
                                        'mediancolor', [0 0 0], ...
                                        'show_outliers', show_outliers, ...
                                        'xlabel', label_string_this_band ...
                                      );
                                end
                            end
                        end
                    end
                end
                if isContinuousVariable(i_variable, data_all, bands_per_stretch)
                    target_abscissa = abscissae_cell{i_comparison, i_variable}(i_condition, :);
                    if strcmp(plot_mode, 'detailed')
                        for i_stretch = 1 : size(data_to_plot_this_condition, 2)
%                             origin_trial_data = ones(size(target_abscissa)) * origin_trial_list_this_condition(i_stretch);
                            origin_index_data = ones(size(target_abscissa)) * origin_indices(i_stretch);
                            plot3 ...
                              ( ...
                                target_axes_handle, ...
                                target_abscissa, ...
                                data_to_plot_this_condition(:, i_stretch), ...
                                origin_index_data, ... %origin_trial_data, ...
                                'HandleVisibility', 'off', ...
                                'color', lightenColor(colors_comparison(i_condition, :), 0.5) ...
                              );
                        end
%                         condition_mean_plot = plot ...
%                           ( ...
%                             target_axes_handle, ...
%                             target_abscissa, ...
%                             mean(data_to_plot_this_condition, 2), ...
%                             'linewidth', 5, ...
%                             'color', colors_comparison(i_condition, :) ...
%                           );                    
                    end
                    if strcmp(plot_mode, 'overview') || strcmp(plot_mode, 'episodes')
                        plot_handles = shadedErrorBar ...
                          ( ...
                            target_abscissa, ...
                            mean(data_to_plot_this_condition, 2), ...
                            spread(data_to_plot_this_condition, spread_method), ...
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
                        
                        if strcmp(plot_mode, 'episodes') && ~strcmp(this_condition{strcmp(condition_combination_labels, 'index')}, 'ONE')
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
            if show_legend && ~(isDiscreteVariable(i_variable, data_all, bands_per_stretch) && (strcmp(plot_mode, 'overview') || strcmp(plot_mode, 'episodes')))
                legend(target_axes_handle, 'show')
            end
        end
    end
    if false % removed for now, look back later
%      for i_path = 1 : number_of_paths_to_plot
%         if plot_settings.get('plot_response')
%             data_to_plot_x = path_response_data_all{i_path, 1};
%             data_to_plot_y = path_response_data_all{i_path, 2};
%         else
%             data_to_plot_x = path_data_all{i_path, 1};
%             data_to_plot_y = path_data_all{i_path, 2};
%         end
%         
%         for i_comparison = 1 : length(comparison_indices)
%             % find correct condition indicator for control
%             conditions_this_comparison = comparison_indices{i_comparison};
%             top_level_plots = [];
%             target_axes_handle = path_axes_handles(comparison_path_to_axes_index_map(i_comparison), i_path);
%             
%             % plot control
%             if plot_settings.get('plot_control') && ~isempty(conditions_control) && ~plot_settings.get('plot_response')
%                 % determine which control condition applies here
%                 representant_condition_index = conditions_this_comparison(1);
%                 applicable_control_condition_index = findApplicableControlConditionIndex(conditions_to_plot(representant_condition_index, :), conditions_control);
%                 applicable_control_condition_labels = conditions_control(applicable_control_condition_index, :);
%                 
%                 % extract data for control condition
%                 stance_foot_indicator = strcmp(condition_stance_foot_list_all, applicable_control_condition_labels{1});
%                 perturbation_indicator = strcmp(condition_perturbation_list_all, applicable_control_condition_labels{2});
%                 delay_indicator = strcmp(condition_delay_list_all, applicable_control_condition_labels{3});
%                 index_indicator = strcmp(condition_index_list_all, applicable_control_condition_labels{4});
%                 experimental_indicator = strcmp(condition_experimental_list_all, applicable_control_condition_labels{5});
%                 stimulus_indicator = strcmp(condition_stimulus_list_all, applicable_control_condition_labels{6});
%                 day_indicator = strcmp(condition_day_list_all, applicable_control_condition_labels{7});
%                 this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
%                 data_to_plot_this_condition_x = data_to_plot_x(:, this_condition_indicator);
%                 data_to_plot_this_condition_y = data_to_plot_y(:, this_condition_indicator);
% %                 origin_indices = find(this_condition_indicator);
%                 
%                 if ~isempty(data_to_plot_this_condition_x)
%                     if strcmp(plot_mode, 'detailed')
%                         % individual trajectories
%                         for i_stretch = 1 : size(data_to_plot_this_condition_x, 2)
% %                             origin_index_data = ones(size(target_abscissa)) * origin_indices(i_stretch);
%                             plot ...
%                               ( ...
%                                 target_axes_handle, ...
%                                 data_to_plot_this_condition_x(:, i_stretch), ...
%                                 data_to_plot_this_condition_y(:, i_stretch), ...
%                                 'HandleVisibility', 'off', ...
%                                 'color', lightenColor(plot_settings.get('color_control'), 0.5) ...
%                               );
%                         end
%                         % condition average
%                         control_mean_plot = plot ...
%                           ( ...
%                             target_axes_handle, ...
%                             mean(data_to_plot_this_condition_x, 2), ...
%                             mean(data_to_plot_this_condition_y, 2), ...
%                             'DisplayName', 'CONTROL', ...
%                             'linewidth', 5, ...
%                             'color', plot_settings.get('color_control') ...
%                           );
%                         top_level_plots = [top_level_plots control_mean_plot]; %#ok<AGROW>
%                     end
%                 end
%             end
%             
%             % plot stimulus
%             for i_condition = 1 : length(conditions_this_comparison)
%                 label_string = strrep(conditions_to_plot{comparison_indices{i_comparison}(i_condition), plot_settings.get('comparison_to_make')}, '_', ' ');
%                 
%                 % find correct condition indicator
%                 condition_identifier = conditions_to_plot(conditions_this_comparison(i_condition), :);
%                 stance_foot_indicator = strcmp(condition_stance_foot_list_all, condition_identifier{1});
%                 perturbation_indicator = strcmp(condition_perturbation_list_all, condition_identifier{2});
%                 delay_indicator = strcmp(condition_delay_list_all, condition_identifier{3});
%                 index_indicator = strcmp(condition_index_list_all, condition_identifier{4});
%                 experimental_indicator = strcmp(condition_experimental_list_all, condition_identifier{5});
%                 stimulus_indicator = strcmp(condition_stimulus_list_all, condition_identifier{6});
%                 day_indicator = strcmp(condition_day_list_all, condition_identifier{7});
%                 this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
%                 data_to_plot_this_condition_x = data_to_plot_x(:, this_condition_indicator);
%                 data_to_plot_this_condition_y = data_to_plot_y(:, this_condition_indicator);
% %                 origin_trial_list_this_condition = origin_trial_list_all(this_condition_indicator);
%                 origin_indices = find(this_condition_indicator);
%                 if strcmp(plot_mode, 'detailed')
%                     for i_stretch = 1 : size(data_to_plot_this_condition_x, 2)
% %                             origin_trial_data = ones(size(target_abscissa)) * origin_trial_list_this_condition(i_stretch);
% %                         origin_index_data = ones(size(target_abscissa)) * origin_indices(i_stretch);
%                         plot ...
%                           ( ...
%                             target_axes_handle, ...
%                             data_to_plot_this_condition_x(:, i_stretch), ...
%                             data_to_plot_this_condition_y(:, i_stretch), ...
%                             'HandleVisibility', 'off', ...
%                             'color', lightenColor(colors_comparison(i_condition, :), 0.5) ...
%                           );
%                     end
%                     condition_mean_plot = plot ...
%                       ( ...
%                         target_axes_handle, ...
%                         mean(data_to_plot_this_condition_x, 2), ...
%                         mean(data_to_plot_this_condition_y, 2), ...
%                         'linewidth', 5, ...
%                         'color', colors_comparison(i_condition, :) ...
%                       );
%                     top_level_plots = [top_level_plots condition_mean_plot]; %#ok<AGROW>
%                 end
%             end
% 
%             % reorder to bring mean plots on top
%             for i_plot = 1 : length(top_level_plots)
%                 uistack(top_level_plots(i_plot), 'top');
%             end
%             
%             % toggle legend
%             if show_legend && ~(isDiscreteVariable(i_variable, data_all, bands_per_stretch) && (strcmp(plot_mode, 'overview') || strcmp(plot_mode, 'episodes')))
%                 legend(target_axes_handle, 'show')
%             end
%         end
    end
    
    %% update label positions
    for i_variable = 1 : number_of_variables_to_plot
        for i_axes = 1 : size(trajectory_axes_handles, 1)
            these_axes = trajectory_axes_handles(i_axes, i_variable);
            xlimits = get(these_axes, 'xlim'); ylimits = get(these_axes, 'ylim');
            if pos_arrow_handles(i_axes, i_variable) ~= 0
                pos_arrow_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.09;
                pos_arrow_position_y = ylimits(2);
                set(pos_arrow_handles(i_axes, i_variable), 'Position', [pos_arrow_position_x pos_arrow_position_y]);
                pos_text_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.14;
                pos_text_position_y = ylimits(2);
                set(pos_text_handles(i_axes, i_variable), 'Position', [pos_text_position_x pos_text_position_y]);
            end
            if neg_arrow_handles(i_axes, i_variable) ~= 0
                neg_arrow_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.09;
                neg_arrow_position_y = ylimits(1);
                set(neg_arrow_handles(i_axes, i_variable), 'Position', [neg_arrow_position_x neg_arrow_position_y]);
                neg_text_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.14;
                neg_text_position_y = ylimits(1);
                set(neg_text_handles(i_axes, i_variable), 'Position', [neg_text_position_x neg_text_position_y]);
            end
            if isDiscreteVariable(i_variable, data_all, bands_per_stretch)
                % rotate labels
                xtick_label_rotation = plot_settings.get('xtick_label_rotation');
                set(these_axes, 'XTickLabelRotation', xtick_label_rotation);
                
            end
            
        end
    end
    
    %% add zero line
    if plot_settings.get('plot_zero')
        for i_variable = 1 : number_of_variables_to_plot
            for i_axes = 1 : size(trajectory_axes_handles, 1)
                these_axes = trajectory_axes_handles(i_axes, i_variable);
                xlimits = get(these_axes, 'xlim');
                zero_plot = plot(these_axes, xlimits, [0 0], 'color', [0.7 0.7 0.7]);
                set(zero_plot, 'HandleVisibility', 'off');
                uistack(zero_plot, 'bottom')
                
            end
        end    
    end
    
    %% shade steps
    if mark_bands
        if strcmp(plot_mode, 'overview')
            for i_comparison = 1 : number_of_comparisons
                for i_variable = 1 : number_of_variables_to_plot
                    if isContinuousVariable(i_variable, data_all, bands_per_stretch)
                        these_axes = trajectory_axes_handles(i_comparison, i_variable);
                        these_abscissae = abscissae_cell{i_comparison, i_variable};
                        ylimits = get(these_axes, 'ylim');

                        for i_band = 2 : 2 : bands_per_stretch
                            % double stance patch
                            double_stance_patch_color = plot_settings.get('stance_double_color');

                            [start_index, end_index] = getBandIndices(i_band, number_of_time_steps_normalized);

                            band_start_times = these_abscissae(:, start_index);
                            band_end_times = these_abscissae(:, end_index);

                            % if these rows are all the same, then we're good and we can mark only a single box.

                            % for testing
                            if length(unique(band_start_times)) == 1 && length(unique(band_end_times)) == 1
                                patch_x = [band_start_times(1) band_end_times(1) band_end_times(1) band_start_times(1)];
                                patch_y = [ylimits(1) ylimits(1) ylimits(2) ylimits(2)];
                                patch_handle = ...
                                    patch ...
                                      ( ...
                                        patch_x, ...
                                        patch_y, ...
                                        double_stance_patch_color, ...
                                        'parent', these_axes, ...
                                        'EdgeColor', 'none', ...
                                        'FaceAlpha', plot_settings.get('stance_alpha'), ...
                                        'HandleVisibility', 'off' ...
                                      ); 
                                uistack(patch_handle, 'bottom')                    
                            end

                            % otherwise we'll have to do something else
                            if length(unique(band_start_times)) > 1 || length(unique(band_end_times)) > 1
                                number_of_conditions = length(band_start_times);
                                y_values = linspace(ylimits(1), ylimits(2), number_of_conditions+1);
                                for i_condition = 1 : number_of_conditions
                                    patch_x = [band_start_times(i_condition) band_end_times(i_condition) band_end_times(i_condition) band_start_times(i_condition)];
                                    patch_y = [y_values(i_condition) y_values(i_condition) y_values(i_condition+1) y_values(i_condition+1)];
                                    patch_handle = ...
                                        patch ...
                                          ( ...
                                            patch_x, ...
                                            patch_y, ...
                                            colors_comparison(i_condition, :), ...
                                            'parent', these_axes, ...
                                            'EdgeColor', 'none', ...
                                            'FaceAlpha', plot_settings.get('stance_alpha'), ...
                                            'HandleVisibility', 'off' ...
                                          ); 
                                    uistack(patch_handle, 'bottom')                    

                                end

                            end


                        end
                    end
                end
            end
        end
    end
    
    if mark_pushoff
        if strcmp(plot_mode, 'overview')
            for i_variable = 1 : number_of_variables_to_plot
                if isContinuousVariable(i_variable, data_all, bands_per_stretch)
                    for i_comparison = 1 : number_of_comparisons
                        these_axes = trajectory_axes_handles(i_comparison, i_variable);
                        ylimits = get(these_axes, 'ylim');

                        step_start_time = step_start_times_cell{i_comparison, i_variable};
                        step_pushoff_time = step_pushoff_times_cell{i_comparison, i_variable};

                        % double stance patch
                        double_stance_patch_color = plot_settings.get('stance_double_color');
                        stretch_start = step_start_time;
                        stretch_end = step_pushoff_time;
                        patch_x = [stretch_start stretch_end stretch_end stretch_start];
                        patch_y = [ylimits(1) ylimits(1) ylimits(2) ylimits(2)];
                        patch_handle = ...
                            patch ...
                              ( ...
                                patch_x, ...
                                patch_y, ...
                                double_stance_patch_color, ...
                                'parent', these_axes, ...
                                'EdgeColor', 'none', ...
                                'FaceAlpha', plot_settings.get('stance_alpha'), ...
                                'HandleVisibility', 'off' ...
                              ); 
                        uistack(patch_handle, 'bottom')                
                    end
                end
            end        
        end
        if strcmp(plot_mode, 'episodes')
            for i_variable = 1 : number_of_variables_to_plot
                for i_episode = 1 : number_of_episodes
                    these_axes = trajectory_axes_handles(i_episode, i_variable);
                    ylimits = get(these_axes, 'ylim');

                    step_start_times = step_start_times_cell{i_episode, i_variable};
%                     step_end_times = step_end_times_cell{i_episode, i_variable};
                    step_pushoff_times = step_pushoff_times_cell{i_episode, i_variable};
    %                 step_stance_foot = step_stance_foot_cell{i_episode, i_variable};

                    for i_step = 1 : length(step_start_times)
                        % double stance patch
                        double_stance_patch_color = plot_settings.get('stance_double_color');
                        stretch_start = step_start_times(i_step);
                        stretch_end = step_pushoff_times(i_step);
                        patch_x = [stretch_start stretch_end stretch_end stretch_start];
                        patch_y = [ylimits(1) ylimits(1) ylimits(2) ylimits(2)];
                        patch_handle = ...
                            patch ...
                              ( ...
                                patch_x, ...
                                patch_y, ...
                                double_stance_patch_color, ...
                                'parent', these_axes, ...
                                'EdgeColor', 'none', ...
                                'FaceAlpha', plot_settings.get('stance_alpha'), ...
                                'HandleVisibility', 'off' ...
                              ); 
                        uistack(patch_handle, 'bottom')

                        % single stance patch
    %                     single_stance_patch_color = [1 1 1] * 0.8;
    %                     if step_stance_foot(i_step) == 0
    %                         single_stance_patch_color = plot_settings.get('stance_double_color');
    %                     end
    %                     if step_stance_foot(i_step) == 1
    %                         single_stance_patch_color = plot_settings.get('stance_left_color');
    %                     end
    %                     if step_stance_foot(i_step) == 2
    %                         single_stance_patch_color = plot_settings.get('stance_right_color');
    %                     end
    %                     stretch_start = step_pushoff_times(i_step);
    %                     stretch_end = step_end_times(i_step);
    %                     patch_x = [stretch_start stretch_end stretch_end stretch_start];
    %                     patch_y = [ylimits(1) ylimits(1) ylimits(2) ylimits(2)];
    %                     patch_handle = ...
    %                         patch ...
    %                           ( ...
    %                             patch_x, ...
    %                             patch_y, ...
    %                             single_stance_patch_color, ...
    %                             'parent', these_axes, ...
    %                             'EdgeColor', 'none', ...
    %                             'FaceAlpha', plot_settings.get('stance_alpha'), ...
    %                             'HandleVisibility', 'off' ...
    %                           ); 
    %                     uistack(patch_handle, 'bottom')
                    end
                end
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
        for i_figure = 1 : numel(trajectory_figure_handles)
            % save with labels
%             legend(axes_handles(i_figure), 'show');
            filename = ['figures' filesep 'withLabels' filesep get(trajectory_figure_handles(i_figure), 'UserData')];
            saveas(trajectory_figure_handles(i_figure), filename, parser.Results.format)
            
            % save without labels
%             set(postext, 'visible', 'off');
%             set(negtext, 'visible', 'off');
            
            set(get(trajectory_axes_handles(i_figure), 'xaxis'), 'visible', 'off');
            set(get(trajectory_axes_handles(i_figure), 'yaxis'), 'visible', 'off');
            set(get(trajectory_axes_handles(i_figure), 'xlabel'), 'visible', 'off');
            set(get(trajectory_axes_handles(i_figure), 'ylabel'), 'visible', 'off');
            set(get(trajectory_axes_handles(i_figure), 'title'), 'visible', 'off');
            set(trajectory_axes_handles(i_figure), 'xticklabel', '');
            set(trajectory_axes_handles(i_figure), 'yticklabel', '');
            set(trajectory_axes_handles(i_figure), 'position', [0 0 1 1]);
            legend(trajectory_axes_handles(i_figure), 'hide');
            filename = ['figures' filesep 'noLabels' filesep get(trajectory_figure_handles(i_figure), 'UserData')];
            saveas(trajectory_figure_handles(i_figure), filename, parser.Results.format);

            close(trajectory_figure_handles(i_figure))            
        end
    end
    
end

%% helper functions

function discrete = isDiscreteVariable(variable_index, variable_data, bands_per_stretch)
    discrete = false;
    if size(variable_data{variable_index}, 1) == bands_per_stretch
        discrete = true;
    end
end

function continuous = isContinuousVariable(variable_index, variable_data, bands_per_stretch)
    continuous = false;
    if size(variable_data{variable_index}, 1) ~= bands_per_stretch
        continuous = true;
    end
end

function s = spread(data, method)
    if strcmp(method, 'cinv')
        s = cinv(data, 2);
    end
    if strcmp(method, 'sem')
        s = std(data, 0, 2) * 1 / sqrt(size(data, 1));
    end
    
end

function [scaled_abscissa, band_limits] = createScaledAbscissa(band_scales, number_of_time_steps_normalized)
    number_of_bands = length(band_scales);
    scaled_abscissa = [];
    band_limits = 0;
    for i_band = 1 : number_of_bands
        scaled_abscissa_this_band = linspace(0, band_scales(i_band), number_of_time_steps_normalized)';
        if i_band > 1
            % start time of this band is end time of the last band, so remove the duplicate point
            scaled_abscissa_this_band = scaled_abscissa_this_band + scaled_abscissa(end);
            scaled_abscissa_this_band = scaled_abscissa_this_band(2:end);
        end
        scaled_abscissa = [scaled_abscissa; scaled_abscissa_this_band]; %#ok<AGROW>
        band_limits = [band_limits scaled_abscissa(end)]; %#ok<AGROW>
    end
    
end


