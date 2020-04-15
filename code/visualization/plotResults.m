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
%   a cell array "comparisons.comparison_indices", where each entry is an array of condition indices, referring to a row in 
%   condition_combination_labels. Usually there will be one comparison per figure. To keep track of this, the code generates
%   an array "figure_data.trajectory_axes_handles" that maps comparisons and variables to axes handles, and a cell array 
%   "abscissae_cell" storing the x-values. For both these arrays, rows = comparisons, columns = variables

function plotResults(varargin)
    settings = determineSettings(varargin{:});

    % load data
    data = loadDataToPlot(settings);
    
    % determine condition combinations to plot and comparisons
    comparisons = createComparisonData(settings, data);
    
    % create figures and determine abscissae for each comparison
    figure_data = createFigureData(settings, data, comparisons);      
    
    % plot data
    figure_data = plotData(settings, data, comparisons, figure_data);
    
    % groom axes, labels etc
    figure_data = groomFigures(settings, data, comparisons, figure_data);
    
    % save and close
    saveFigures(settings, figure_data);
    closeFigures(settings, figure_data);
end

%% helper functions
function settings = determineSettings(varargin)
    settings = struct;
    
    % parse 
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'subjects', [])
    addParameter(parser, 'dictate_axes', false)
    addParameter(parser, 'show_legend', false)
    addParameter(parser, 'save', false)
    addParameter(parser, 'close', false)
    addParameter(parser, 'format', 'jpeg')
    addParameter(parser, 'resolution', '300')
    addParameter(parser, 'settings', 'plot')
    addParameter(parser, 'spread_method', 'cinv')
    parse(parser, varargin{:})
    settings.subjects = parser.Results.subjects;
    settings.dictate_axes = parser.Results.dictate_axes;
    settings.show_legend = parser.Results.show_legend;
    settings.settings_file = parser.Results.settings;
    settings.spread_method = parser.Results.spread_method;
    
    % save format
    settings.save_results = parser.Results.save;
    settings.save_format = ['-d' parser.Results.format];
    settings.save_resolution = ['-r' num2str(parser.Results.resolution)];
    settings.close = parser.Results.close;

    % load settings
    study_settings = loadSettingsFromFile('study');
    plot_settings = loadSettingsFromFile(settings.settings_file);
    
    settings.show_single_data_points = plot_settings.get('show_single_data_points', 1);
    settings.mark_pushoff = plot_settings.get('mark_pushoff', 1);
    settings.mark_bands = plot_settings.get('mark_bands', 1);
    settings.band_labels = study_settings.get('band_labels', 1);
    settings.group_bands_within_conditions = plot_settings.get('group_bands_within_conditions', 1);
    settings.number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    settings.show_average_data = plot_settings.get('show_average_data', 1);
    settings.show_spread_data = plot_settings.get('show_spread_data', 1);
    
    settings.study_settings = study_settings;
    settings.plot_settings = plot_settings;
    
    % colors
    settings.colors_bands = settings.plot_settings.get('colors_bands', 1);
    settings.colors_header = settings.plot_settings.get('colors_header', 1);
    settings.colors = settings.plot_settings.get('colors', 1);
    
    % extract and store settings from files
    settings.conditions_settings = settings.study_settings.get('conditions');
    settings.condition_to_compare = settings.plot_settings.get('condition_to_compare');
    settings.condition_labels = settings.conditions_settings(:, 1)';
    settings.variables_to_plot = settings.plot_settings.get('variables_to_plot');
    
    settings.variables_to_plot_header = settings.plot_settings.get('variables_to_plot_header', true);

    settings.number_of_variables_to_plot = size(settings.variables_to_plot, 1);
end

function data = loadDataToPlot(settings)
    data = struct;

    % declare variables
    condition_source_variables = settings.conditions_settings(:, 2)';
    number_of_condition_labels = length(settings.condition_labels);
    data_source_file = settings.plot_settings.get('data_source', 1);

    file_label = ['results' data_source_file];
    
    % load data
    data_folder_list = determineDataStructure(settings.subjects);
    data.condition_data = {};
    origin_trial_list_all = [];
    origin_start_time_list_all = [];
    origin_end_time_list_all = [];
    data.variable_data = cell(settings.number_of_variables_to_plot, 1);
    data.directions = cell(settings.number_of_variables_to_plot, 2);
    
    data.step_time_data = [];
    pushoff_time_data = [];
    data.bands_per_stretch = [];
    
    for i_folder = 1 : length(data_folder_list)
        % get information
        this_data_folder_path = data_folder_list{i_folder};
        subject_settings = loadSettingsFromFile('subject', this_data_folder_path);
        collection_date = subject_settings.get('collection_date');
        subject_id = subject_settings.get('subject_id');
        
        % find results file
        results_file_candidate_analysis = [this_data_folder_path filesep 'analysis' filesep makeFileName(collection_date, subject_id, file_label) '.mat'];
        results_file_candidate_subject = [this_data_folder_path filesep makeFileName(collection_date, subject_id, file_label) '.mat'];
        results_file_candidate_results = [this_data_folder_path filesep 'results' filesep  makeFileName(collection_date, subject_id, file_label) '.mat'];
        if exist(results_file_candidate_analysis, 'file')
            results_file_name = results_file_candidate_analysis;
        end    
        if exist(results_file_candidate_subject, 'file')
            results_file_name = results_file_candidate_subject;
        end    
        if exist(results_file_candidate_results, 'file')
            results_file_name = results_file_candidate_results;
        end    
        
        % load data
        disp(['loading data from ' results_file_name])
        loaded_data = load(results_file_name);
        number_of_stretches_this_session = length(loaded_data.time_list_session);
        bands_per_stretch_this_session = loaded_data.bands_per_stretch;

        % transform conditions into cell array
        conditions_session = loaded_data.conditions_session;
        condition_array_session = cell(number_of_stretches_this_session, number_of_condition_labels);
        for i_condition = 1 : number_of_condition_labels
            condition_array_session(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
        end
        
        data.condition_data = [data.condition_data; condition_array_session];
        origin_trial_list_all = [origin_trial_list_all; loaded_data.origin_trial_list_session]; %#ok<AGROW>
        origin_start_time_list_all = [origin_start_time_list_all; loaded_data.origin_start_time_list_session]; %#ok<AGROW>
        origin_end_time_list_all = [origin_end_time_list_all; loaded_data.origin_end_time_list_session]; %#ok<AGROW>
        
        % extract data
        for i_variable = 1 : settings.number_of_variables_to_plot
            
            this_variable_name = settings.variables_to_plot{i_variable, strcmp(settings.variables_to_plot_header, 'variable name')};
            this_variable_type = settings.variables_to_plot{i_variable, strcmp(settings.variables_to_plot_header, 'variable type')};
            this_variable_source_index = find(strcmp(loaded_data.([this_variable_type '_names_session']), this_variable_name), 1, 'first');
            if isempty(this_variable_source_index)
                error(['Variable not found: ' this_variable_name])
            end
            this_variable_data = loaded_data.([this_variable_type '_data_session']){this_variable_source_index};
            this_variable_directions = loaded_data.([this_variable_type '_directions_session'])(this_variable_source_index, :);

            if settings.plot_settings.get('convert_to_mm', 1) && (strcmp(this_variable_name,'cop_from_com_x') || strcmp(this_variable_name, 'step_placement_x'))
                this_variable_data = this_variable_data * 1000;
            end
            
            % store
            data.variable_data{i_variable} = [data.variable_data{i_variable} this_variable_data];
            data.directions(i_variable, :) = this_variable_directions;
        end
        
        % get time variables
        if isfield(loaded_data, 'stretch_names_session') && any(find(strcmp(loaded_data.stretch_names_session, 'step_time')))
            index_in_saved_data = find(strcmp(loaded_data.stretch_names_session, 'step_time'), 1, 'first');
            this_step_time_data = loaded_data.stretch_data_session{index_in_saved_data};
            data.step_time_data = [data.step_time_data this_step_time_data];
        end
        if isfield(loaded_data, 'stretch_names_session') && any(find(strcmp(loaded_data.stretch_names_session, 'pushoff_time')))
            index_in_saved_data = find(strcmp(loaded_data.stretch_names_session, 'pushoff_time'), 1, 'first');
            this_pushoff_time_data = loaded_data.stretch_data_session{index_in_saved_data};
            pushoff_time_data = [pushoff_time_data this_pushoff_time_data]; %#ok<AGROW>
        end
        if isempty(data.bands_per_stretch)
            data.bands_per_stretch = bands_per_stretch_this_session;
        else
            if data.bands_per_stretch ~= bands_per_stretch_this_session
               warning('Different sessions have different numbers of bands per stretch') 
            end
        end
    end
    
    % calculate mean pushoff index
    if settings.mark_pushoff
        pushoff_time_ratio = pushoff_time_data ./ data.step_time_data;
        mean_pushoff_ratio = mean(pushoff_time_ratio, 2);
        data.pushoff_index = round(mean_pushoff_ratio * 100);
    end
    
end

function comparisons = createComparisonData(settings, data)
    % create container and extract some settings
    comparisons = struct;
    labels_to_ignore = settings.plot_settings.get('conditions_to_ignore');
    levels_to_remove = settings.plot_settings.get('levels_to_remove');
    preferred_level_order = settings.plot_settings.get('preferred_level_order', 1);
    
    % determine comparisons and auxiliary data
    [comparisons.condition_combination_labels, comparisons.condition_combinations] = determineConditionCombinations(data.condition_data, settings.conditions_settings, labels_to_ignore, levels_to_remove);
    comparisons.condition_combinations = sortConditionCombinations(comparisons.condition_combinations, comparisons.condition_combination_labels, settings.condition_to_compare, preferred_level_order);
    [comparisons.comparison_indices, comparisons.conditions_per_comparison_max] = determineComparisons(comparisons.condition_combinations, comparisons.condition_combination_labels, settings);
    comparisons.number_of_comparisons = length(comparisons.comparison_indices);
    
    % determine colors for the combinations
    comparisons.condition_colors = determineConditionColors(settings, comparisons);
end

function condition_colors = determineConditionColors(settings, comparisons)
    % find unique levels of condition to compare
    levels = unique(comparisons.condition_combinations(:, strcmp(comparisons.condition_combination_labels, settings.condition_to_compare)));
    
    % make default color map
    default_colors = lines(length(levels));
    
    % get colors from settings for this condition
    condition_column = find(strcmp(settings.colors_header, 'condition'));
    level_column = find(strcmp(settings.colors_header, 'level'));
    color_column = find(strcmp(settings.colors_header, 'color'));
    if ~isempty(settings.colors)
        colors_from_settings = settings.colors(strcmp(settings.colors(:, condition_column), settings.condition_to_compare), [level_column color_column]); %#ok<FNDSB>
    else
        colors_from_settings = cell(0, 2);
    end
    
    % go through levels and store default color or the one provided in the settings
    condition_colors = [levels cell(size(levels))];
    for i_level = 1 : length(levels)
        this_level = levels(i_level);
        if any(strcmp(colors_from_settings(:, 1), this_level))
            % use color provided in settings
            condition_colors{i_level, 2} = hex2rgb(colors_from_settings{strcmp(colors_from_settings(:, 1), this_level), 2});
        else
            % use default color
            condition_colors{i_level, 2} = default_colors(i_level, :);
        end
        
    end
end

function figure_data = createFigureData(settings, data, comparisons)
    figure_data.comparison_variable_to_axes_index_map = zeros(comparisons.number_of_comparisons, 1);
    figure_data.abscissae_cell = cell(comparisons.number_of_comparisons, settings.number_of_variables_to_plot);
    
    % time plots
    % make one figure per comparison and variable
    figure_data.trajectory_figure_handles = zeros(comparisons.number_of_comparisons, settings.number_of_variables_to_plot);
    figure_data.trajectory_axes_handles = zeros(comparisons.number_of_comparisons, settings.number_of_variables_to_plot);
    figure_data.pos_text_handles = zeros(comparisons.number_of_comparisons, settings.number_of_variables_to_plot);
    figure_data.neg_text_handles = zeros(comparisons.number_of_comparisons, settings.number_of_variables_to_plot);
    figure_data.pos_arrow_handles = zeros(comparisons.number_of_comparisons, settings.number_of_variables_to_plot);
    figure_data.neg_arrow_handles = zeros(comparisons.number_of_comparisons, settings.number_of_variables_to_plot);
    step_start_times_cell = cell(comparisons.number_of_comparisons, settings.number_of_variables_to_plot);
    step_end_times_cell = cell(comparisons.number_of_comparisons, settings.number_of_variables_to_plot);
    for i_variable = 1 : settings.number_of_variables_to_plot
        for i_comparison = 1 : comparisons.number_of_comparisons
            this_comparison = comparisons.comparison_indices{i_comparison};
            % make figure and axes
            new_figure = figure; new_axes = axes; hold on;

            % store handles and determine abscissa data
            figure_data.trajectory_figure_handles(i_comparison, i_variable) = new_figure;
            figure_data.trajectory_axes_handles(i_comparison, i_variable) = new_axes;
            figure_data.comparison_variable_to_axes_index_map(i_comparison) = i_comparison;

            if isDiscreteVariable(i_variable, data.variable_data, data.bands_per_stretch)
                % abscissa gives the bin edges here
                this_comparison = comparisons.comparison_indices{i_comparison};
                number_of_entries = length(this_comparison);

                if settings.group_bands_within_conditions
                    gap_between_conditions = 1;

                    abscissae_stimulus = repmat((1 : data.bands_per_stretch)', 1, length(this_comparison));
                    shifter = (0:number_of_entries-1) * (data.bands_per_stretch + gap_between_conditions);
                    abscissae_stimulus = abscissae_stimulus + repmat(shifter, data.bands_per_stretch, 1);
                    if settings.plot_settings.get('merge_bands', 1)
                        abscissae_stimulus = abscissae_stimulus(1, :);
                    end

                else
                    gap_between_bands = 1;

%                             ab
                    abscissae_stimulus = repmat((1 : number_of_entries), data.bands_per_stretch, 1);
                    shifter = (0:data.bands_per_stretch-1)' * (comparisons.conditions_per_comparison_max + gap_between_bands);
                    abscissae_stimulus = abscissae_stimulus + repmat(shifter, 1, comparisons.conditions_per_comparison_max);
                    if settings.plot_settings.get('merge_bands', 1)
                        abscissae_stimulus = abscissae_stimulus(1, :);
                    end
                end
                figure_data.abscissae_cell{i_comparison, i_variable} = abscissae_stimulus;
            end
            if isContinuousVariable(i_variable, data.variable_data, data.bands_per_stretch)
                % scale abscissae
                conditions_this_comparison = comparisons.comparison_indices{i_comparison};
                step_time_means_this_comparison = zeros(data.bands_per_stretch, size(conditions_this_comparison, 2));
                for i_condition = 1 : length(conditions_this_comparison)
                    this_condition_combination = comparisons.condition_combinations(conditions_this_comparison(i_condition), :);
                    this_condition_indicator = getConditionIndicator(this_condition_combination, comparisons.condition_combination_labels, data.condition_data, settings.condition_labels);
                    step_time_data_this_condition = data.step_time_data(:, this_condition_indicator);
                    step_time_means_this_comparison(:, i_condition) = mean(step_time_data_this_condition, 2);
                end
                for i_condition = 1 : length(conditions_this_comparison)
                    if strcmp(settings.plot_settings.get('time_plot_style'), 'scaled_to_comparison_mean')
                        band_scales = mean(step_time_means_this_comparison, 2);
                    elseif strcmp(settings.plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                        band_scales = step_time_means_this_comparison(:, i_condition);
                    else
                        band_scales = ones(data.bands_per_stretch, 1) * (settings.number_of_time_steps_normalized-1);
                    end
                    [abscissa_scaled, band_limits] = createScaledAbscissa(band_scales, settings.number_of_time_steps_normalized);

                    if strcmp(settings.plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                        time_plot_band_anchor_index = settings.plot_settings.get('time_plot_band_anchor');
                        time_plot_band_anchor_time = band_limits(time_plot_band_anchor_index);
                        abscissa_scaled = abscissa_scaled - time_plot_band_anchor_time;
                    end

                    figure_data.abscissae_cell{i_comparison, i_variable}(i_condition, :) = abscissa_scaled;


                end                    
            end

            % set axes properties
            if isDiscreteVariable(i_variable, data.variable_data, data.bands_per_stretch)
                xtick = sort(reshape(figure_data.abscissae_cell{i_comparison, i_variable}, 1, numel(figure_data.abscissae_cell{i_comparison, i_variable})));
                set(gca, 'xlim', [-0.5 + min(xtick) 0.5 + max(xtick(end))]);
                set(gca, 'xtick', xtick);
            end

            % set axis labels
            if isContinuousVariable(i_variable, data.variable_data, data.bands_per_stretch)
                if strcmp(settings.plot_settings.get('time_plot_style'), 'scaled_to_comparison_mean') || strcmp(settings.plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                    xlabel('normalized time (s)');
                else
                    xlabel('normalized time (%)');
                end
            end
            this_label = settings.variables_to_plot{i_variable, strcmp(settings.variables_to_plot_header, 'y-axis label')};
            ylabel(this_label);

            % add text labels
            figure_data.pos_text_handles(i_comparison, i_variable) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    data.directions{i_variable, 1}, ...
                    'rotation', 90, ...
                    'Fontsize', 24, ...
                    'horizontalalignment', 'right', ...
                    'parent', new_axes ...
                  );                
              figure_data.pos_arrow_handles(i_comparison, i_variable) = ...
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
            figure_data.neg_text_handles(i_comparison, i_variable) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    data.directions{i_variable, 2}, ...
                    'rotation', 90, ...
                    'Fontsize', 24, ...
                    'horizontalalignment', 'left', ...
                    'parent', new_axes...
                  );
            figure_data.neg_arrow_handles(i_comparison, i_variable) = ...
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
            title_string = settings.variables_to_plot{i_variable, strcmp(settings.variables_to_plot_header, 'variable label')};
            filename_string = settings.variables_to_plot{i_variable, strcmp(settings.variables_to_plot_header, 'save file string')};

            representative_condition = comparisons.condition_combinations(this_comparison(1), :);

            for i_label = 1 : length(representative_condition)
                if ~(strcmp(comparisons.condition_combination_labels{i_label}, settings.condition_to_compare))
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
    for i_variable = 1 : settings.number_of_variables_to_plot
        if isContinuousVariable(i_variable, data.variable_data, data.bands_per_stretch)
            for i_comparison = 1 : comparisons.number_of_comparisons
                target_abscissae = figure_data.abscissae_cell{i_comparison, i_variable};
                xlim = [min(target_abscissae(:, 1)) max(target_abscissae(:, end))];

                % set x-limits accordingly
                set(figure_data.trajectory_axes_handles(i_comparison, i_variable), 'xlim', xlim);
            end
        end
    end

    % determine stance start and end times and stance foot
    for i_variable = 1 : settings.number_of_variables_to_plot
        if isContinuousVariable(i_variable, data.variable_data, data.bands_per_stretch)
            for i_comparison = 1 : comparisons.number_of_comparisons

                step_abscissa = figure_data.abscissae_cell{i_comparison, i_variable};
                step_start_times_cell{i_comparison, i_variable} = step_abscissa(1, 1);
                step_end_times_cell{i_comparison, i_variable} = step_abscissa(1, end);
            end
        end
    end        
end

function figure_data = plotData(settings, data, comparisons, figure_data)
    for i_variable = 1 : settings.number_of_variables_to_plot
        data_to_plot = data.variable_data{i_variable, 1};
        for i_comparison = 1 : length(comparisons.comparison_indices)
            % find correct condition indicator for control
            conditions_this_comparison = comparisons.comparison_indices{i_comparison};
            top_level_plots = [];
            target_axes_handle = figure_data.trajectory_axes_handles(figure_data.comparison_variable_to_axes_index_map(i_comparison), i_variable);
            
            % plot stimulus
            for i_condition = 1 : length(conditions_this_comparison)
                this_condition_index = conditions_this_comparison(i_condition);
                this_condition = comparisons.condition_combinations(this_condition_index, :);
                this_label = this_condition{strcmp(comparisons.condition_combination_labels, settings.condition_to_compare)};
                this_color = comparisons.condition_colors{strcmp(comparisons.condition_colors(:, 1), this_label), 2};
                label_string = strrep(this_label, '_', ' ');
                this_condition_indicator = getConditionIndicator(this_condition, comparisons.condition_combination_labels, data.condition_data, settings.condition_labels);
                data_to_plot_this_condition = data_to_plot(:, this_condition_indicator);
                
                origin_indices = find(this_condition_indicator);
                if isDiscreteVariable(i_variable, data.variable_data, data.bands_per_stretch)
                    if settings.plot_settings.get('merge_bands', 1)
                        data_to_plot_this_condition = reshape(data_to_plot_this_condition, 1, numel(data_to_plot_this_condition));
                    end
                    for i_band = 1 : size(data_to_plot_this_condition, 1)
                        if ~isempty(settings.band_labels)
                            label_string_this_band = [label_string '-' settings.band_labels{i_band}];
                        else
                            label_string_this_band = label_string;
                        end
                        target_abscissa = figure_data.abscissae_cell{i_comparison, i_variable}(i_band, i_condition);
                        data_to_plot_this_band = data_to_plot_this_condition(i_band, :);
                        if ~any(isnan(data_to_plot_this_band))
                            if settings.group_bands_within_conditions
                                % override color
                                colors = copper(size(data_to_plot_this_condition, 1));
                                this_color = colors(i_band, :);
                            end
                            plotDiscreteData ...
                              ( ...
                                data_to_plot_this_band, ...
                                'abscissa', target_abscissa, ...
                                'axes', target_axes_handle, ...         % axes
                                'color', this_color, ...                % color
                                'ShowMean', settings.show_average_data, ...
                                'MeanStyle', 'd', ...
                                'MeanColor', [1 1 1]*0.7, ...
                                'ShowMedian', settings.show_spread_data, ...
                                'MedianStyle', 'line', ...
                                'ShowIndividualData', settings.plot_settings.get('show_individual_discrete_data', 1), ...
                                'ShowSpread', settings.show_spread_data, ...
                                'SpreadStyle', settings.plot_settings.get('discrete_data_plot_style'), ...
                                'label', label_string_this_band ...     % label
                              )
                            
                        end
                    end
                end
                if isContinuousVariable(i_variable, data.variable_data, data.bands_per_stretch)
                    target_abscissa = figure_data.abscissae_cell{i_comparison, i_variable}(i_condition, :);                    
                    if settings.plot_settings.get('show_individual_trajectory_data', 1)
                        for i_stretch = 1 : size(data_to_plot_this_condition, 2)
                            origin_index_data = - ones(size(target_abscissa)) * origin_indices(i_stretch);
                            plot3 ...
                              ( ...
                                target_axes_handle, ...
                                target_abscissa, ...
                                data_to_plot_this_condition(:, i_stretch), ...
                                origin_index_data, ...
                                'linewidth', 1, ...
                                'HandleVisibility', 'off', ...
                                'color', lightenColor(this_color, 0.5) ...
                              );
                        end
                        
                    end
                    if settings.show_average_data
                        average_plot = plot ...
                          ( ...
                            target_abscissa, ...
                            nanmean(data_to_plot_this_condition, 2), ...
                            'parent', target_axes_handle, ...
                            'DisplayName', label_string, ...
                            'color', this_color, ...
                            'linewidth', 6 ...
                          );
                        top_level_plots = [top_level_plots average_plot]; %#ok<AGROW>
                    end
                    if settings.show_spread_data
                        plot_handles = shadedErrorBar ...
                          ( ...
                            target_abscissa, ...
                            nanmean(data_to_plot_this_condition, 2), ...
                            spread(data_to_plot_this_condition, settings.spread_method), ...
                            { ...
                              'color', this_color, ...
                              'linewidth', 6 ...
                            }, ...
                            1, ...
                            target_axes_handle ...
                          );
                        set(plot_handles.patch, 'HandleVisibility', 'off');
                        delete(plot_handles.edge);
                        delete(plot_handles.mainLine);
                    end                   
                end
            end
        end
    end
end

function figure_data = groomFigures(settings, data, comparisons, figure_data)
    % set axis limits
    if settings.dictate_axes
        for i_variable = 1 : settings.number_of_variables_to_plot
            % get x-axis limits from settings
            this_variable_x_lower = settings.variables_to_plot{i_variable, strcmp(settings.variables_to_plot_header, 'x-axis lower limit')};
            this_variable_x_upper = settings.variables_to_plot{i_variable, strcmp(settings.variables_to_plot_header, 'x-axis upper limit')};
            
            % get y-axis limits from settings
            this_variable_y_lower = settings.variables_to_plot{i_variable, strcmp(settings.variables_to_plot_header, 'y-axis lower limit')};
            this_variable_y_upper = settings.variables_to_plot{i_variable, strcmp(settings.variables_to_plot_header, 'y-axis upper limit')};
            
            for i_axes = 1 : size(figure_data.trajectory_axes_handles, 1)
                % get current axes and limits
                these_axes = figure_data.trajectory_axes_handles(i_axes, i_variable);
                xlimits = get(these_axes, 'xlim');
                ylimits = get(these_axes, 'ylim');
                
                % apply new limits if any were set
                if ~strcmp(this_variable_x_lower, '~')
                    xlimits(1) = str2double(this_variable_x_lower);
                end
                if ~strcmp(this_variable_x_upper, '~')
                    xlimits(2) = str2double(this_variable_x_upper);
                end
                if ~strcmp(this_variable_y_lower, '~')
                    ylimits(1) = str2double(this_variable_y_lower);
                end
                if ~strcmp(this_variable_y_upper, '~')
                    ylimits(2) = str2double(this_variable_y_upper);
                end
                set(these_axes, 'xlim', xlimits, 'ylim', ylimits);
                
            end
        end
    end

    % update label positions
    for i_variable = 1 : settings.number_of_variables_to_plot
        for i_axes = 1 : size(figure_data.trajectory_axes_handles, 1)
            these_axes = figure_data.trajectory_axes_handles(i_axes, i_variable);
            xlimits = get(these_axes, 'xlim'); ylimits = get(these_axes, 'ylim');
            if figure_data.pos_arrow_handles(i_axes, i_variable) ~= 0
                pos_arrow_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.09;
                pos_arrow_position_y = ylimits(2);
                set(figure_data.pos_arrow_handles(i_axes, i_variable), 'Position', [pos_arrow_position_x pos_arrow_position_y]);
                pos_text_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.14;
                pos_text_position_y = ylimits(2);
                set(figure_data.pos_text_handles(i_axes, i_variable), 'Position', [pos_text_position_x pos_text_position_y]);
            end
            if figure_data.neg_arrow_handles(i_axes, i_variable) ~= 0
                neg_arrow_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.09;
                neg_arrow_position_y = ylimits(1);
                set(figure_data.neg_arrow_handles(i_axes, i_variable), 'Position', [neg_arrow_position_x neg_arrow_position_y]);
                neg_text_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.14;
                neg_text_position_y = ylimits(1);
                set(figure_data.neg_text_handles(i_axes, i_variable), 'Position', [neg_text_position_x neg_text_position_y]);
            end
            if isDiscreteVariable(i_variable, data.variable_data, data.bands_per_stretch)
                % rotate labels
                xtick_label_rotation = settings.plot_settings.get('xtick_label_rotation', 1);
                set(these_axes, 'XTickLabelRotation', xtick_label_rotation);
                
            end
            
        end
    end

    % toggle legend
    for i_variable = 1 : settings.number_of_variables_to_plot
        for i_axes = 1 : size(figure_data.trajectory_axes_handles, 1)
            these_axes = figure_data.trajectory_axes_handles(i_axes, i_variable);
            if settings.show_legend && ~(isDiscreteVariable(i_variable, data.variable_data, data.bands_per_stretch))
                legend(these_axes, 'show')
            end
        end
    end
    
    % add zero line
    if settings.plot_settings.get('plot_zero', 1)
        for i_variable = 1 : settings.number_of_variables_to_plot
            for i_axes = 1 : size(figure_data.trajectory_axes_handles, 1)
                these_axes = figure_data.trajectory_axes_handles(i_axes, i_variable);
                xlimits = get(these_axes, 'xlim');
                zero_plot = plot(these_axes, xlimits, [0 0], 'color', [0.7 0.7 0.7]);
                set(zero_plot, 'HandleVisibility', 'off');
                uistack(zero_plot, 'bottom')
                
            end
        end    
    end
    
    % mark bands
    if settings.mark_bands
        for i_comparison = 1 : comparisons.number_of_comparisons
            for i_variable = 1 : settings.number_of_variables_to_plot
                if isContinuousVariable(i_variable, data.variable_data, data.bands_per_stretch)
                    these_axes = figure_data.trajectory_axes_handles(i_comparison, i_variable);
                    these_abscissae = figure_data.abscissae_cell{i_comparison, i_variable};
                    ylimits = get(these_axes, 'ylim');

                    if settings.mark_bands == 1
                        bands_to_mark = 2 : 2 : data.bands_per_stretch;
                    end
                    if settings.mark_bands == 2
                        bands_to_mark = 1 : 2 : data.bands_per_stretch;
                    end

                    for i_band = bands_to_mark
                        % double stance patch
                        double_stance_patch_color = settings.plot_settings.get('stance_double_color');

                        [start_index, end_index] = getBandIndices(i_band, settings.number_of_time_steps_normalized);

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
                                    'FaceAlpha', settings.plot_settings.get('stance_alpha'), ...
                                    'HandleVisibility', 'off' ...
                                  ); 
                            uistack(patch_handle, 'bottom')                    
                        end

                        % otherwise we'll have to do something else
                        if length(unique(band_start_times)) > 1 || length(unique(band_end_times)) > 1
                            number_of_conditions = length(band_start_times);
                            y_values = linspace(ylimits(1), ylimits(2), number_of_conditions+1);
                            for i_condition = 1 : number_of_conditions
% 2020-APR-15 HR: this is code to mark the double stance for each condition
% with an individual box, with a light shade of the condition color. I
% don't have a good way to access the condition color here, after
% re-working the way these colors are determined. Since this is used very
% rarely, I'll leave it as gray boxes for now, using 
% double_stance_patch_color, to be fixed if it's actually needed
% 
%                                 this_condition_index = conditions_this_comparison(i_condition);
%                                 this_condition = comparisons.condition_combinations(this_condition_index, :);
%                                 this_label = this_condition{strcmp(comparisons.condition_combination_labels, settings.condition_to_compare)};
%                                 this_color = comparisons.condition_colors{strcmp(comparisons.condition_colors(:, 1), this_label), 2};
                                
                                patch_x = [band_start_times(i_condition) band_end_times(i_condition) band_end_times(i_condition) band_start_times(i_condition)];
                                patch_y = [y_values(i_condition) y_values(i_condition) y_values(i_condition+1) y_values(i_condition+1)];
                                patch_handle = ...
                                    patch ...
                                      ( ...
                                        patch_x, ...
                                        patch_y, ...
                                        double_stance_patch_color, ...
                                        'parent', these_axes, ...
                                        'EdgeColor', 'none', ...
                                        'FaceAlpha', settings.plot_settings.get('stance_alpha'), ...
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
    
    % mark pushoff
    if settings.mark_pushoff
        for i_comparison = 1 : comparisons.number_of_comparisons
            for i_variable = 1 : settings.number_of_variables_to_plot
                if isContinuousVariable(i_variable, data.variable_data, data.bands_per_stretch)
                    these_axes = figure_data.trajectory_axes_handles(i_comparison, i_variable);
                    these_abscissae = figure_data.abscissae_cell{i_comparison, i_variable};
                    ylimits = get(these_axes, 'ylim');
                    
                    for i_band = 1 : data.bands_per_stretch
                        % double stance patch
                        double_stance_patch_color = settings.plot_settings.get('stance_double_color');

                        start_index = getBandIndices(i_band, settings.number_of_time_steps_normalized);
                        pushoff_index_here = start_index + data.pushoff_index(i_band);

                        band_start_times = these_abscissae(:, start_index);
                        band_end_times = these_abscissae(:, pushoff_index_here);

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
                                'FaceAlpha', settings.plot_settings.get('stance_alpha'), ...
                                'HandleVisibility', 'off' ...
                              ); 
                        uistack(patch_handle, 'bottom')                    
                    end
                end
            end
        end
    end
    
end

function saveFigures(settings, figure_data)
    % save figures
    if settings.save_results
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
        for i_figure = 1 : numel(figure_data.trajectory_figure_handles)
            % remove some white space on right side and top
            axes_position = get(figure_data.trajectory_axes_handles(i_figure), 'position');
            axes_position(3) = 1 - axes_position(1) - 0.01;
            axes_position(4) = 1 - axes_position(2) - 0.04;
            set(figure_data.trajectory_axes_handles(i_figure), 'position', axes_position);
            
            % save with labels            
            filename_with = ['figures' filesep 'withLabels' filesep get(figure_data.trajectory_figure_handles(i_figure), 'UserData')];
            print(figure_data.trajectory_figure_handles(i_figure), filename_with, settings.save_format, settings.save_resolution)
            
            % remove text and marks to save data lines only
            set(get(figure_data.trajectory_axes_handles(i_figure), 'xaxis'), 'visible', 'off');
            set(get(figure_data.trajectory_axes_handles(i_figure), 'yaxis'), 'visible', 'off');
            set(get(figure_data.trajectory_axes_handles(i_figure), 'xlabel'), 'visible', 'off');
            set(get(figure_data.trajectory_axes_handles(i_figure), 'ylabel'), 'visible', 'off');
            set(get(figure_data.trajectory_axes_handles(i_figure), 'title'), 'visible', 'off');
            set(figure_data.trajectory_axes_handles(i_figure), 'xticklabel', '');
            set(figure_data.trajectory_axes_handles(i_figure), 'yticklabel', '');
            set(figure_data.trajectory_axes_handles(i_figure), 'position', [0 0 1 1]);
            legend(figure_data.trajectory_axes_handles(i_figure), 'hide');
            filename_without = ['figures' filesep 'noLabels' filesep get(figure_data.trajectory_figure_handles(i_figure), 'UserData')];
            print(figure_data.trajectory_figure_handles(i_figure), filename_without, settings.save_format, settings.save_resolution)
            disp(['Saved as ' filename_with ' and ' filename_without])
            
            % put some marks back
            set(get(figure_data.trajectory_axes_handles(i_figure), 'title'), 'visible', 'on');
            set(figure_data.trajectory_axes_handles(i_figure), 'position', [0.05 0.05 0.9 0.9]);
        end
    end

end

function closeFigures(settings, figure_data)
    % close figures
    if settings.close
        for i_figure = 1 : numel(figure_data.trajectory_figure_handles)
            close(figure_data.trajectory_figure_handles(i_figure))            
        end
    end    

end

function discrete = isDiscreteVariable(variable_index, variable_data, bands_per_stretch)
    discrete = false;
    if size(variable_data{variable_index}, 1) == bands_per_stretch
        discrete = true;
    end
    if size(variable_data{variable_index}, 1) == 1
        discrete = true;
    end
end

function continuous = isContinuousVariable(variable_index, variable_data, bands_per_stretch)
    continuous = false;
    if size(variable_data{variable_index}, 1) > bands_per_stretch
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



