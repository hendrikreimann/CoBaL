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
    
    % make data custodian object and load data
    data_source = settings.plot_settings.get('data_source', 1);
    subjects = settings.subjects;
    data_custodian = StretchDataCustodian(pwd, data_source, subjects);
    
    % transform old variables list to new format
    settings = transformVariablesListToNewFormat(settings, data_custodian);
                
    % determine condition combinations to plot and comparisons
    condition_data = data_custodian.getConditionData();
    comparisons = createComparisonData(settings, condition_data);
    
    % create figures and determine abscissae for each comparison
    figure_data_continuous = createFigures_continuous(settings, comparisons, data_custodian);
    figure_data_discrete = createFigures_discrete(settings, comparisons, data_custodian);
    figure_data_range = createFigures_range(settings, comparisons);
    figure_data_pairs = createFigures_pairs(settings, comparisons);
    
    % plot data
    plotData_continuous(settings, comparisons, data_custodian, figure_data_continuous);
    plotData_discrete(settings, comparisons, data_custodian, figure_data_discrete);
    plotData_range(settings, comparisons, data_custodian, figure_data_range);
    plotData_pairs(settings, comparisons, data_custodian, figure_data_pairs);
    
    % groom axes, labels etc
    groomFigures_continuous(settings, data_custodian, comparisons, figure_data_continuous);
    groomFigures_discrete(settings, figure_data_discrete);
    groomFigures_range(settings, figure_data_range);
    groomFigures_pairs(settings, figure_data_pairs);
    
    % save and close
    saveFigures(settings, figure_data_continuous);
    saveFigures(settings, figure_data_discrete);
    saveFigures(settings, figure_data_range);
    saveFigures(settings, figure_data_pairs);
    closeFigures(settings, figure_data_continuous);
    closeFigures(settings, figure_data_discrete);

end

% helper functions
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
    settings.colors_table = settings.plot_settings.getTable('colors', 1);
    
    % extract and store settings from files
    settings.conditions_settings = settings.study_settings.get('conditions');
    settings.condition_to_compare = settings.plot_settings.get('condition_to_compare');
    settings.condition_labels = settings.conditions_settings(:, 1)';
    settings.variables_to_plot_table = settings.plot_settings.getTable('variables_to_plot', true);
    settings.number_of_variables_to_plot = size(settings.variables_to_plot_table, 1);
    
    settings.variables_to_plot_discrete = settings.plot_settings.get('variables_to_plot_discrete', true);
    settings.variables_to_plot_discrete_header = settings.plot_settings.get('variables_to_plot_discrete_header', true);
    settings.number_of_variables_to_plot_discrete = size(settings.variables_to_plot_discrete, 1);
    
    settings.variables_to_plot_continuous = settings.plot_settings.get('variables_to_plot_continuous', true);
    settings.variables_to_plot_continuous_header = settings.plot_settings.get('variables_to_plot_continuous_header', true);
    settings.number_of_variables_to_plot_continuous = size(settings.variables_to_plot_continuous, 1);
    
    settings.variables_to_plot_range = settings.plot_settings.get('variables_to_plot_range', true);
    settings.variables_to_plot_range_header = settings.plot_settings.get('variables_to_plot_range_header', true);
    settings.number_of_variables_to_plot_range = size(settings.variables_to_plot_range, 1);
    settings.variables_to_plot_pairs = settings.plot_settings.getTable('variables_to_plot_pairs', true);
    settings.number_of_variables_to_plot_pairs = size(settings.variables_to_plot_pairs, 1);
    
end

function settings = transformVariablesListToNewFormat(settings, data_custodian)
    % grab variables
    variables_to_plot_table = settings.variables_to_plot_table;
    if isempty(variables_to_plot_table)
        % old list entry exists, but is empty, we don't have to do anything
        return
    end

    % go through variables and transform one by one
    for i_variable = 1 : length(variables_to_plot_table.variable_name)
        this_variable_name = variables_to_plot_table.variable_name{i_variable};
        this_variable_type = variables_to_plot_table.variable_type{i_variable};
        this_variable_data = data_custodian.getData(this_variable_name, this_variable_type);
        
        this_row = variables_to_plot_table(i_variable, :);
        % check size of this variable
        if size(this_variable_data.variable_data, 1) == data_custodian.bands_per_stretch
            % this is a discrete variable
            settings.variables_to_plot_discrete = addRowToVariableTable(this_row, settings.variables_to_plot_discrete, settings.variables_to_plot_discrete_header);
        else
            % this is a continuous variable
            settings.variables_to_plot_continuous = addRowToVariableTable(this_row, settings.variables_to_plot_continuous, settings.variables_to_plot_continuous_header);
        end
    end
    
    % update length of new format lists
    settings.number_of_variables_to_plot_discrete = size(settings.variables_to_plot_discrete, 1);
    settings.number_of_variables_to_plot_continuous = size(settings.variables_to_plot_continuous, 1);
end

function list = addRowToVariableTable(row_to_add, existing_list, existing_list_header)
    new_row_reformatted = cell(1, length(existing_list_header));
    for i_column = 1 : length(existing_list_header)
        % find this column in old list
        this_label = existing_list_header{i_column};
        if any(strcmp(row_to_add.Properties.VariableNames, this_label))
            this_entry = row_to_add{1, strcmp(row_to_add.Properties.VariableNames, this_label)};
            this_entry = this_entry{1}; % temporary fix to transform from table to cell array
        else
            this_entry = '~';
        end
        new_row_reformatted{i_column} = this_entry;
    end
    list = [existing_list; new_row_reformatted];
end

function comparisons = createComparisonData(settings, condition_data)
    % create container and extract some settings
    comparisons = struct;
    labels_to_ignore = settings.plot_settings.get('conditions_to_ignore');
    levels_to_remove = settings.plot_settings.get('levels_to_remove');
    preferred_level_order = settings.plot_settings.get('preferred_level_order', 1);
    
    % determine comparisons and auxiliary data
    [comparisons.condition_combination_labels, comparisons.condition_combinations] = determineConditionCombinations(condition_data, settings.conditions_settings, labels_to_ignore, levels_to_remove);
    [comparisons.condition_combinations, comparisons.level_order] = sortConditionCombinations(comparisons.condition_combinations, comparisons.condition_combination_labels, settings.condition_to_compare, preferred_level_order);
    
    condition_to_compare = settings.plot_settings.get('condition_to_compare');
    condition_table = settings.conditions_settings;
    [comparisons.comparison_indices, comparisons.conditions_per_comparison_max] = determineComparisons(comparisons.condition_combinations, comparisons.condition_combination_labels, condition_to_compare, condition_table);
    comparisons.number_of_comparisons = length(comparisons.comparison_indices);
    
    % determine colors for the combinations
    comparisons.condition_to_compare = settings.condition_to_compare;
    comparisons.condition_colors = determineConditionColors(settings.colors_table, comparisons, settings.condition_to_compare);
end

function figure_data = createFigures_continuous(settings, comparisons, data_custodian)
    figure_data = createFigureData(settings.variables_to_plot_continuous, comparisons);
    if settings.number_of_variables_to_plot_continuous == 0
        return
    end
    
    step_start_times_cell = cell(comparisons.number_of_comparisons, settings.number_of_variables_to_plot_continuous);
    step_end_times_cell = cell(comparisons.number_of_comparisons, settings.number_of_variables_to_plot_continuous);
    bands_per_stretch = data_custodian.bands_per_stretch;
    condition_data = data_custodian.getConditionData();
    stretch_times = data_custodian.getStretchTimes();
    
    % create figures
    for i_variable = 1 : settings.number_of_variables_to_plot_continuous
        for i_comparison = 1 : comparisons.number_of_comparisons
            % make figure and axes
            new_figure = figure;
            dcm_obj = datacursormode(new_figure);
            set(dcm_obj,'UpdateFcn',{@singlePlotTooltip})
            new_axes = axes; hold on;

            % store handles and determine abscissa data
            figure_data.figure_handles{i_comparison, i_variable} = new_figure;
            figure_data.axes_handles{i_comparison, i_variable} = new_axes;
            figure_data.comparison_variable_to_axes_index_map(i_comparison) = i_comparison;

            % scale abscissae
            conditions_this_comparison = comparisons.comparison_indices{i_comparison};
            step_time_means_this_comparison = zeros(bands_per_stretch, size(conditions_this_comparison, 2));
            for i_condition = 1 : length(conditions_this_comparison)
                this_condition_combination = comparisons.condition_combinations(conditions_this_comparison(i_condition), :);
                this_condition_indicator = getConditionIndicator ...
                  ( ...
                    this_condition_combination, ...
                    comparisons.condition_combination_labels, ...
                    condition_data, ...
                    settings.condition_labels, ...
                    settings.plot_settings.get('levels_to_remove') ...
                  );
                stretch_time_data_this_condition = stretch_times(:, this_condition_indicator);
                stretch_duration_data_this_condition = diff(stretch_time_data_this_condition, 1);
                step_time_means_this_comparison(:, i_condition) = mean(stretch_duration_data_this_condition, 2);
            end
            for i_condition = 1 : length(conditions_this_comparison)
                if strcmp(settings.plot_settings.get('time_plot_style'), 'scaled_to_comparison_mean')
                    band_scales = mean(step_time_means_this_comparison, 2);
                elseif strcmp(settings.plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                    band_scales = step_time_means_this_comparison(:, i_condition);
                else
                    band_scales = ones(size(step_time_means_this_comparison, 1), 1) * (settings.number_of_time_steps_normalized-1);
                end
                [abscissa_scaled, band_limits] = createScaledAbscissa(band_scales, settings.number_of_time_steps_normalized);

                if strcmp(settings.plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                    time_plot_band_anchor_index = settings.plot_settings.get('time_plot_band_anchor');
                    time_plot_band_anchor_time = band_limits(time_plot_band_anchor_index);
                    abscissa_scaled = abscissa_scaled - time_plot_band_anchor_time;
                end

                figure_data.abscissae_cell{i_comparison, i_variable}(i_condition, :) = abscissa_scaled;


            end                    

            % set axis labels
            if strcmp(settings.plot_settings.get('time_plot_style'), 'scaled_to_comparison_mean') || strcmp(settings.plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                xlabel('normalized time (s)');
            else
                xlabel('normalized time (%)');
            end
        end
    end

    % set x-limits
    for i_variable = 1 : settings.number_of_variables_to_plot_continuous
        for i_comparison = 1 : comparisons.number_of_comparisons
            target_abscissae = figure_data.abscissae_cell{i_comparison, i_variable};
            xlim = [min(target_abscissae(:, 1)) max(target_abscissae(:, end))];

            % set x-limits accordingly
            set(figure_data.axes_handles{i_comparison, i_variable}, 'xlim', xlim);
        end
    end

    % determine stance start and end times and stance foot
    for i_variable = 1 : settings.number_of_variables_to_plot_continuous
        for i_comparison = 1 : comparisons.number_of_comparisons

            step_abscissa = figure_data.abscissae_cell{i_comparison, i_variable};
            step_start_times_cell{i_comparison, i_variable} = step_abscissa(1, 1);
            step_end_times_cell{i_comparison, i_variable} = step_abscissa(1, end);
        end
    end        

    % add labels
    figure_data = addLabelsAndData(figure_data, comparisons, settings.variables_to_plot_continuous, settings.variables_to_plot_continuous_header);
    
end

function figure_data = createFigures_discrete(settings, comparisons, data_custodian)
    figure_data = createFigureData(settings.variables_to_plot_discrete, comparisons);
    bands_per_stretch = data_custodian.bands_per_stretch;
    
    for i_variable = 1 : settings.number_of_variables_to_plot_discrete
        for i_comparison = 1 : comparisons.number_of_comparisons
            % make figure and axes
            new_figure = figure; 
            dcm_obj = datacursormode(new_figure);
            set(dcm_obj,'UpdateFcn',{@singlePlotTooltip})
            new_axes = axes; hold on;

            % store handles and determine abscissa data
            figure_data.figure_handles{i_comparison, i_variable} = new_figure;
            figure_data.axes_handles{i_comparison, i_variable} = new_axes;
            figure_data.comparison_variable_to_axes_index_map(i_comparison) = i_comparison;

            % abscissa gives the bin edges here
            conditions_this_comparison = comparisons.comparison_indices{i_comparison};

            if settings.group_bands_within_conditions
                number_of_entries = comparisons.conditions_per_comparison_max;
                gap_between_conditions = 1;

                abscissae_base = repmat((1 : data_custodian.bands_per_stretch)', 1, number_of_entries);
                shifter = (0:number_of_entries-1) * (data_custodian.bands_per_stretch + gap_between_conditions);
                abscissae_base = abscissae_base + repmat(shifter, data_custodian.bands_per_stretch, 1);
                if settings.plot_settings.get('merge_bands', 1)
                    abscissae_base = abscissae_base(1, :);
                end

                % go through and select appropriate one based on label
                abscissae_stimulus = zeros(bands_per_stretch, comparisons.conditions_per_comparison_max) * NaN;
                for i_condition = 1 : length(conditions_this_comparison)
                    this_condition_index = conditions_this_comparison(i_condition);
                    this_condition = comparisons.condition_combinations(this_condition_index, :);
                    this_label = this_condition{strcmp(comparisons.condition_combination_labels, settings.condition_to_compare)};
                    this_label_index_in_level_order = strcmp(this_label, comparisons.level_order);
                    abscissae_stimulus(:, i_condition) = abscissae_base(:, this_label_index_in_level_order);
                end
                
                
            else
                gap_between_bands = 1;
                
                % go through each level individually, find its place in the level order and determine its abscissa
                abscissae_template = ones(1, comparisons.conditions_per_comparison_max) * NaN;
                for i_condition = 1 : length(conditions_this_comparison)
                    this_condition_index = conditions_this_comparison(i_condition);
                    this_condition = comparisons.condition_combinations(this_condition_index, :);
                    this_label = this_condition{strcmp(comparisons.condition_combination_labels, settings.condition_to_compare)};
                    this_label_index_in_level_order = find(strcmp(this_label, comparisons.level_order));
                    abscissae_template(i_condition) = this_label_index_in_level_order;
                end
                
                % remove gaps
                i_index = 1;
                while i_index < max(abscissae_template)
                    while ~ismember(i_index, abscissae_template)
                        % while this index is not present, decrement all indices above this one
                        abscissae_template(abscissae_template > i_index) = abscissae_template(abscissae_template > i_index) - 1;
                    end
                    i_index = i_index + 1;
                end
                
                % distribute
                abscissae_stimulus = repmat(abscissae_template, bands_per_stretch, 1);
                shifter = (0:bands_per_stretch-1)' * (comparisons.conditions_per_comparison_max + gap_between_bands);
                abscissae_stimulus = abscissae_stimulus + repmat(shifter, 1, comparisons.conditions_per_comparison_max);
                if settings.plot_settings.get('merge_bands', 1)
                    abscissae_stimulus = abscissae_stimulus(1, :);
                end
            end
            figure_data.abscissae_cell{i_comparison, i_variable} = abscissae_stimulus;

            % set axes properties
            xtick = reshape(figure_data.abscissae_cell{i_comparison, i_variable}, 1, numel(figure_data.abscissae_cell{i_comparison, i_variable}));
            xtick(isnan(xtick)) = [];
            xtick = sort(xtick);
            set(gca, 'xlim', [-0.5 + min(xtick) 0.5 + max(xtick(end))]);
            set(gca, 'xtick', xtick);
        end
    end
    
    % add labels
    figure_data = addLabelsAndData(figure_data, comparisons, settings.variables_to_plot_discrete, settings.variables_to_plot_discrete_header);
end

function figure_data = createFigures_range(settings, comparisons)
    figure_data = createFigureData(settings.variables_to_plot_range, comparisons);
    
    for i_variable = 1 : settings.number_of_variables_to_plot_range
        for i_comparison = 1 : comparisons.number_of_comparisons
            % make figure and axes
            new_figure = figure;
            dcm_obj = datacursormode(new_figure);
            set(dcm_obj,'UpdateFcn',{@singlePlotTooltip})
            new_axes = axes; hold on;

            % store handles and determine abscissa data
            figure_data.figure_handles{i_comparison, i_variable} = new_figure;
            figure_data.axes_handles{i_comparison, i_variable} = new_axes;
            figure_data.comparison_variable_to_axes_index_map(i_comparison) = i_comparison;

            % abscissa gives the bin edges here
            conditions_this_comparison = comparisons.comparison_indices{i_comparison};

            % go through and select appropriate one based on label
            abscissae_stimulus = zeros(1, comparisons.conditions_per_comparison_max) * NaN;
            for i_condition = 1 : length(conditions_this_comparison)
                this_condition_index = conditions_this_comparison(i_condition);
                this_condition = comparisons.condition_combinations(this_condition_index, :);
                this_label = this_condition{strcmp(comparisons.condition_combination_labels, settings.condition_to_compare)};
                this_label_index_in_level_order = strcmp(this_label, comparisons.level_order);
                abscissae_stimulus(i_condition) = find(this_label_index_in_level_order);
            end
                
                
            figure_data.abscissae_cell{i_comparison, i_variable} = abscissae_stimulus;

            % set axes properties
            xtick = reshape(figure_data.abscissae_cell{i_comparison, i_variable}, 1, numel(figure_data.abscissae_cell{i_comparison, i_variable}));
            xtick(isnan(xtick)) = [];
            xtick = sort(xtick);
            set(gca, 'xlim', [-0.5 + min(xtick) 0.5 + max(xtick(end))]);
            set(gca, 'xtick', xtick);
        end
    end
    
    % add labels
    figure_data = addLabelsAndData(figure_data, comparisons, settings.variables_to_plot_range, settings.variables_to_plot_range_header);
end

function figure_data = createFigures_pairs(settings, comparisons)
    figure_data = createFigureData(settings.variables_to_plot_pairs, comparisons);
%     bands_per_stretch = data_custodian.bands_per_stretch;
    
    for i_variable = 1 : settings.number_of_variables_to_plot_pairs
        for i_comparison = 1 : comparisons.number_of_comparisons
            % make figure and axes
            new_figure = figure; 
            new_axes = axes('Position', [0.1300 0.1500 0.7750 0.7750]); hold on;

            % store handles and determine abscissa data
            figure_data.figure_handles{i_comparison, i_variable} = new_figure;
            figure_data.axes_handles{i_comparison, i_variable} = new_axes;
            figure_data.comparison_variable_to_axes_index_map(i_comparison) = i_comparison;
        end
    end
    
    % add labels
    figure_data = addLabelsAndData_pairs(figure_data, comparisons, settings.variables_to_plot_pairs);
end

function text = singlePlotTooltip(~, event_obj)
    this_line = event_obj.Target;
    origin = this_line.UserData;

    % Customizes text of data tips
    pos = get(event_obj,'Position');
    text = { ...
            ['X: ', num2str(pos(1))], ...
            ['Y: ', num2str(pos(2))] ...
          };
%     if isempty(origin)
%         text = position_text;
%     else
%         origin_txt = { ...
%                        ['Subject: ', origin.subject], ...
%                        ['Trial: ', num2str(origin.trial)], ...
%                        ['Start Time: ', num2str(origin.start_time)], ...
%                        ['End Time: ', num2str(origin.end_time)] ...
%                      };
%         text = [position_text, origin_txt];
%     end
    
    if isfield(origin, 'subject')
        text = [text, ['Subject: ', origin.subject]];
    end
    if isfield(origin, 'trial')
        text = [text, ['Trial: ', num2str(origin.trial)]];
    end
    if isfield(origin, 'start_time')
        text = [text, ['Start Time: ', num2str(origin.start_time)]];
    end
    if isfield(origin, 'start_time')
        text = [text, ['End Time: ', num2str(origin.end_time)]];
    end
    
end

function figure_data = addLabelsAndData(figure_data, comparisons, variables_to_plot, variables_to_plot_header)
    for i_variable = 1 : size(variables_to_plot, 1)
        for i_comparison = 1 : comparisons.number_of_comparisons
            these_axes = figure_data.axes_handles{i_comparison, i_variable};
            this_figure = figure_data.figure_handles{i_comparison, i_variable};
            % add text labels for arrows
            arrow_text = 'TBD';
            figure_data.pos_text_handles_ordinate(i_comparison, i_variable) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    arrow_text, ...
                    'rotation', 90, ...
                    'Fontsize', 24, ...
                    'horizontalalignment', 'right', ...
                    'parent', these_axes ...
                  );                
            figure_data.pos_arrow_handles_ordinate(i_comparison, i_variable) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    ' $\rightarrow$', ...
                    'rotation', 90, ...
                    'Fontsize', 36, ...
                    'horizontalalignment', 'right', ...
                    'interpreter', 'LaTeX', ...
                    'parent', these_axes ...
                  );
            figure_data.neg_text_handles_ordinate(i_comparison, i_variable) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    arrow_text, ...
                    'rotation', 90, ...
                    'Fontsize', 24, ...
                    'horizontalalignment', 'left', ...
                    'parent', these_axes...
                  );
            figure_data.neg_arrow_handles_ordinate(i_comparison, i_variable) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    '$\leftarrow$ ', ...
                    'rotation', 90, ...
                    'Fontsize', 36, ...
                    'horizontalalignment', 'left', ...
                    'interpreter', 'LaTeX', ...
                    'parent', these_axes ...
                  );
              
            % set axis labels
            this_label = variables_to_plot{i_variable, strcmp(variables_to_plot_header, 'y_axis_label')};
            ylabel(these_axes, this_label);
            
            % determine title
            title_string = variables_to_plot{i_variable, strcmp(variables_to_plot_header, 'variable_label')};
            filename_string = variables_to_plot{i_variable, strcmp(variables_to_plot_header, 'save_file_string')};
            this_comparison = comparisons.comparison_indices{i_comparison};
            representative_condition = comparisons.condition_combinations(this_comparison(1), :);

            for i_label = 1 : length(representative_condition)
                this_is_the_condition_to_compare = (strcmp(comparisons.condition_combination_labels{i_label}, comparisons.condition_to_compare));
                this_label_is_the_same_for_all_conditions = length(unique(comparisons.condition_combinations(:, i_label))) == 1;
                
                if ~this_is_the_condition_to_compare && ~this_label_is_the_same_for_all_conditions
                    this_string = strrep(representative_condition{i_label}, '_', '');
                    filename_string = [filename_string '_' this_string]; %#ok<AGROW>
                    title_string = [title_string ' - ' this_string]; %#ok<AGROW>
                end
            end

            title(these_axes, title_string); 
            set(these_axes, 'Fontsize', 12)
            set(this_figure, 'UserData', filename_string)
        end
    end
end

function figure_data = addLabelsAndData_pairs(figure_data, comparisons, variables_to_plot)
    for i_variable = 1 : size(variables_to_plot, 1)
        for i_comparison = 1 : comparisons.number_of_comparisons
            these_axes = figure_data.axes_handles{i_comparison, i_variable};
            this_figure = figure_data.figure_handles{i_comparison, i_variable};
            arrow_text = 'TBD';

            % add text labels for arrows - x axis
            figure_data.pos_text_handles_abscissa(i_comparison, i_variable) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    arrow_text, ...
                    'rotation', 0, ...
                    'Fontsize', 24, ...
                    'horizontalalignment', 'right', ...
                    'parent', these_axes ...
                  );                
            figure_data.pos_arrow_handles_abscissa(i_comparison, i_variable) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    ' $\rightarrow$', ...
                    'rotation', 0, ...
                    'Fontsize', 36, ...
                    'horizontalalignment', 'right', ...
                    'interpreter', 'LaTeX', ...
                    'parent', these_axes ...
                  );
            figure_data.neg_text_handles_abscissa(i_comparison, i_variable) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    arrow_text, ...
                    'rotation', 0, ...
                    'Fontsize', 24, ...
                    'horizontalalignment', 'left', ...
                    'parent', these_axes...
                  );
            figure_data.neg_arrow_handles_abscissa(i_comparison, i_variable) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    '$\leftarrow$ ', ...
                    'rotation', 0, ...
                    'Fontsize', 36, ...
                    'horizontalalignment', 'left', ...
                    'interpreter', 'LaTeX', ...
                    'parent', these_axes ...
                  );
              
            % add text labels for arrows - y axis
            figure_data.pos_text_handles_ordinate(i_comparison, i_variable) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    arrow_text, ...
                    'rotation', 90, ...
                    'Fontsize', 24, ...
                    'horizontalalignment', 'right', ...
                    'parent', these_axes ...
                  );                
            figure_data.pos_arrow_handles_ordinate(i_comparison, i_variable) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    ' $\rightarrow$', ...
                    'rotation', 90, ...
                    'Fontsize', 36, ...
                    'horizontalalignment', 'right', ...
                    'interpreter', 'LaTeX', ...
                    'parent', these_axes ...
                  );
            figure_data.neg_text_handles_ordinate(i_comparison, i_variable) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    arrow_text, ...
                    'rotation', 90, ...
                    'Fontsize', 24, ...
                    'horizontalalignment', 'left', ...
                    'parent', these_axes...
                  );
            figure_data.neg_arrow_handles_ordinate(i_comparison, i_variable) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    '$\leftarrow$ ', ...
                    'rotation', 90, ...
                    'Fontsize', 36, ...
                    'horizontalalignment', 'left', ...
                    'interpreter', 'LaTeX', ...
                    'parent', these_axes ...
                  );
              
            % set axis labels
            this_label_A = variables_to_plot.variable_A_label{i_variable};
            xlabel(these_axes, this_label_A);
            this_label_B = variables_to_plot.variable_B_label{i_variable};
            ylabel(these_axes, this_label_B);
            
            % determine title
            title_string = variables_to_plot.title{i_variable};
            filename_string = variables_to_plot.filename{i_variable};
            this_comparison = comparisons.comparison_indices{i_comparison};
            representative_condition = comparisons.condition_combinations(this_comparison(1), :);

            for i_label = 1 : length(representative_condition)
                this_is_the_condition_to_compare = (strcmp(comparisons.condition_combination_labels{i_label}, comparisons.condition_to_compare));
                this_label_is_the_same_for_all_conditions = length(unique(comparisons.condition_combinations(:, i_label))) == 1;
                
                if ~this_is_the_condition_to_compare && ~this_label_is_the_same_for_all_conditions
                    this_string = strrep(representative_condition{i_label}, '_', '');
                    filename_string = [filename_string '_' this_string]; %#ok<AGROW>
                    title_string = [title_string ' - ' this_string]; %#ok<AGROW>
                end
            end

            title(these_axes, title_string); 
            set(these_axes, 'Fontsize', 12)
            set(this_figure, 'UserData', filename_string)
        end
    end
end

function figure_data = createFigureData(variables_to_plot, comparisons)
    number_of_variables_to_plot = size(variables_to_plot, 1);
    figure_data.comparison_variable_to_axes_index_map = zeros(comparisons.number_of_comparisons, 1);
    figure_data.abscissae_cell = cell(comparisons.number_of_comparisons, number_of_variables_to_plot);
    
    % make one figure per comparison and variable
    figure_data.figure_handles = cell(comparisons.number_of_comparisons, number_of_variables_to_plot);
    figure_data.axes_handles = cell(comparisons.number_of_comparisons, number_of_variables_to_plot);
    figure_data.pos_text_handles_abscissa = zeros(comparisons.number_of_comparisons, number_of_variables_to_plot);
    figure_data.neg_text_handles_abscissa = zeros(comparisons.number_of_comparisons, number_of_variables_to_plot);
    figure_data.pos_text_handles_ordinate = zeros(comparisons.number_of_comparisons, number_of_variables_to_plot);
    figure_data.neg_text_handles_ordinate = zeros(comparisons.number_of_comparisons, number_of_variables_to_plot);
    figure_data.pos_arrow_handles_abscissa = zeros(comparisons.number_of_comparisons, number_of_variables_to_plot);
    figure_data.neg_arrow_handles_abscissa = zeros(comparisons.number_of_comparisons, number_of_variables_to_plot);
    figure_data.pos_arrow_handles_ordinate = zeros(comparisons.number_of_comparisons, number_of_variables_to_plot);
    figure_data.neg_arrow_handles_ordinate = zeros(comparisons.number_of_comparisons, number_of_variables_to_plot);
end

function plotData_continuous(settings, comparisons, data_custodian, figure_data)
    [condition_data, condition_labels] = data_custodian.getConditionData();
    [origin_subjects, origin_trials, origin_start_times, origin_end_times] = data_custodian.getOriginData();
    for i_variable = 1 : settings.number_of_variables_to_plot_continuous
        variable_name = settings.variables_to_plot_continuous(i_variable, strcmp(settings.variables_to_plot_continuous_header, 'variable_name'));
        variable_type = settings.variables_to_plot_continuous(i_variable, strcmp(settings.variables_to_plot_continuous_header, 'variable_type'));
        scale_factor = settings.variables_to_plot_continuous(i_variable, strcmp(settings.variables_to_plot_continuous_header, 'scale_factor'));
        if isempty(scale_factor)
            scale_factor = 1;
        else
            scale_factor = str2double(scale_factor);
        end
        data_to_plot = data_custodian.getData(variable_name, variable_type);
        
        for i_comparison = 1 : length(comparisons.comparison_indices)
            % find correct condition indicator for control
            conditions_this_comparison = comparisons.comparison_indices{i_comparison};
            top_level_plots = [];
            target_axes_handle = figure_data.axes_handles{figure_data.comparison_variable_to_axes_index_map(i_comparison), i_variable};
            
            % plot
            for i_condition = 1 : length(conditions_this_comparison)
                this_condition_index = conditions_this_comparison(i_condition);
                this_condition = comparisons.condition_combinations(this_condition_index, :);
                this_label = this_condition{strcmp(comparisons.condition_combination_labels, settings.condition_to_compare)};
                this_color = comparisons.condition_colors{strcmp(comparisons.condition_colors(:, 1), this_label), 2};
                label_string = strrep(this_label, '_', ' ');

                this_condition_indicator = getConditionIndicator ...
                  ( ...
                    this_condition, ...
                    comparisons.condition_combination_labels, ...
                    condition_data, ...
                    settings.condition_labels, ...
                    settings.plot_settings.get('levels_to_remove') ...
                  );
                data_to_plot_this_condition = data_to_plot.variable_data(:, this_condition_indicator) * scale_factor;
                origin_subjects_this_condition = origin_subjects(this_condition_indicator);
                origin_trials_this_condition = origin_trials(this_condition_indicator);
                origin_start_times_this_condition = origin_start_times(this_condition_indicator);
                origin_end_times_this_condition = origin_end_times(this_condition_indicator);
                
                if settings.plot_settings.get('average_within_subjects', 1)
                    condition_data_this_condition = condition_data(this_condition_indicator, :);
                    [data_to_plot_this_condition, origin_subjects_this_condition] ...
                        = averageWithinSubjects(data_to_plot_this_condition, condition_data_this_condition, condition_labels);
                end                
                
                
                
                target_abscissa = figure_data.abscissae_cell{i_comparison, i_variable}(i_condition, :);                    
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
                if settings.plot_settings.get('show_individual_trajectory_data', 1)
                    for i_stretch = 1 : size(data_to_plot_this_condition, 2)
                        this_origin = struct;
                        this_origin.subject = origin_subjects_this_condition{i_stretch};
                        if ~settings.plot_settings.get('average_within_subjects', 1)
                            this_origin.trial = origin_trials_this_condition(i_stretch);
                            this_origin.start_time = origin_start_times_this_condition(i_stretch);
                            this_origin.end_time = origin_end_times_this_condition(i_stretch);
                        end
                        plot ...
                          ( ...
                            target_axes_handle, ...
                            target_abscissa, ...
                            data_to_plot_this_condition(:, i_stretch), ...
                            'linewidth', 1, ...
                            'HandleVisibility', 'off', ...
                            'color', lightenColor(this_color, 0.5), ...
                            'UserData', this_origin ...
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
                % TODO: order these different kinds of lines properly
                % across conditions, so that lines in one condition are not
                % un-clickable below spread patches from another condition
            end
            
            % update direction labels
            set(figure_data.pos_text_handles_ordinate(i_comparison, i_variable), 'string', data_to_plot.directions{1});
            set(figure_data.neg_text_handles_ordinate(i_comparison, i_variable), 'string', data_to_plot.directions{2});
        end
    end
end

function plotData_discrete(settings, comparisons, data_custodian, figure_data)
    [condition_data, condition_labels] = data_custodian.getConditionData();
    [origin_subjects, origin_trials, origin_start_times, origin_end_times] = data_custodian.getOriginData();
    for i_variable = 1 : settings.number_of_variables_to_plot_discrete
        variable_name = settings.variables_to_plot_discrete(i_variable, strcmp(settings.variables_to_plot_discrete_header, 'variable_name'));
        variable_type = settings.variables_to_plot_discrete(i_variable, strcmp(settings.variables_to_plot_discrete_header, 'variable_type'));
        data_to_plot = data_custodian.getData(variable_name, variable_type);
        
        for i_comparison = 1 : length(comparisons.comparison_indices)
            % find correct condition indicator for control
            conditions_this_comparison = comparisons.comparison_indices{i_comparison};
            target_axes_handle = figure_data.axes_handles{figure_data.comparison_variable_to_axes_index_map(i_comparison), i_variable};
            
            % plot
            for i_condition = 1 : length(conditions_this_comparison)
                this_condition_index = conditions_this_comparison(i_condition);
                this_condition = comparisons.condition_combinations(this_condition_index, :);
                this_label = this_condition{strcmp(comparisons.condition_combination_labels, settings.condition_to_compare)};
                this_color = comparisons.condition_colors{strcmp(comparisons.condition_colors(:, 1), this_label), 2};
                label_string = strrep(this_label, '_', ' ');
                this_condition_indicator = getConditionIndicator ...
                  ( ...
                    this_condition, ...
                    comparisons.condition_combination_labels, ...
                    condition_data, ...
                    settings.condition_labels, ...
                    settings.plot_settings.get('levels_to_remove') ...
                  );
                data_to_plot_this_condition = data_to_plot.variable_data(:, this_condition_indicator);
                origin_subjects_this_condition = origin_subjects(this_condition_indicator);
                origin_trials_this_condition = origin_trials(this_condition_indicator);
                origin_start_times_this_condition = origin_start_times(this_condition_indicator);
                origin_end_times_this_condition = origin_end_times(this_condition_indicator);
                
                
                if settings.plot_settings.get('average_within_subjects', 1)
                    condition_data_this_condition = condition_data(this_condition_indicator, :);

%                     data_to_plot_this_condition = averageWithinSubjects(data_to_plot_this_condition, condition_data_this_condition, condition_labels);
                    [data_to_plot_this_condition, origin_subjects_this_condition] ...
                        = averageWithinSubjects(data_to_plot_this_condition, condition_data_this_condition, condition_labels);
                end                
                if settings.plot_settings.get('merge_bands', 1)
                    data_to_plot_this_condition = reshape(data_to_plot_this_condition, 1, numel(data_to_plot_this_condition));
                end
                
                % build origin struct
                origin_data = [];
                for i_stretch = 1 : size(data_to_plot_this_condition, 2)
                    this_origin = struct;
                    if settings.plot_settings.get('average_within_subjects', 1)
                        this_origin.subject = origin_subjects_this_condition;
                    else
                        this_origin.subject = origin_subjects_this_condition{i_stretch};
                        this_origin.trial = origin_trials_this_condition(i_stretch);
                        this_origin.start_time = origin_start_times_this_condition(i_stretch);
                        this_origin.end_time = origin_end_times_this_condition(i_stretch);
                    end
                    origin_data = [origin_data; this_origin];
                end                
                
                for i_band = 1 : size(data_to_plot_this_condition, 1)
                    if ~isempty(settings.band_labels) && ~settings.plot_settings.get('merge_bands', 1)
                        label_string_this_band = strrep([label_string '-' settings.band_labels{i_band}], '_', ' ');
                    else
                        label_string_this_band = strrep(label_string, '_', ' ');
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
                            'MeanStyle', settings.plot_settings.get('mean_style', 1), ...
                            'MeanMarkerSize', settings.plot_settings.get('mean_marker_size', 1), ...
                            'MeanColor', lightenColor(this_color, 1-settings.plot_settings.get('mean_data_saturation', 1)), ...
                            'ShowMedian', settings.show_spread_data, ...
                            'MedianStyle', 'line', ...
                            'ShowIndividualData', settings.plot_settings.get('show_individual_discrete_data', 1), ...
                            'IndividualDataMarkerStyle', settings.plot_settings.get('individual_data_marker_style', 1), ...
                            'IndividualDataMarkerSize', settings.plot_settings.get('individual_data_marker_size', 1), ...
                            'IndividualDataColor', lightenColor(this_color, 1-settings.plot_settings.get('individual_data_saturation', 1)), ...
                            'IndividualDataJitterStd', settings.plot_settings.get('individual_data_jitter_std', 1), ...
                            'SpreadFaceColor', lightenColor(this_color, 1-settings.plot_settings.get('spread_fill_saturation', 1)), ...
                            'SpreadEdgeColor', lightenColor(this_color, 1-settings.plot_settings.get('spread_edge_saturation', 1)), ...
                            'SpreadEdgeLinewidth', settings.plot_settings.get('discrete_data_spread_edge_linewidth', 1), ...
                            'SpreadWiskColor', lightenColor(this_color, 1-settings.plot_settings.get('spread_edge_saturation', 1)), ...
                            'ShowSpread', settings.show_spread_data, ...
                            'SpreadStyle', settings.plot_settings.get('discrete_data_plot_style'), ...
                            'IndividualDataInfo', origin_data, ...
                            'label', label_string_this_band ...     % label
                          );
%                             'SpreadFaceAlpha', settings.plot_settings.get('discrete_data_spread_face_alpha', 1), ...

                    end
                end
                
                
            end
            
            % update direction labels
            set(figure_data.pos_text_handles_ordinate(i_comparison, i_variable), 'string', data_to_plot.directions{1});
            set(figure_data.neg_text_handles_ordinate(i_comparison, i_variable), 'string', data_to_plot.directions{2});
        end
    end
end

function plotData_range(settings, comparisons, data_custodian, figure_data)
    [condition_data, condition_labels] = data_custodian.getConditionData();
    [origin_subjects, origin_trials, origin_start_times, origin_end_times] = data_custodian.getOriginData();
    for i_variable = 1 : settings.number_of_variables_to_plot_range
        variable_name = settings.variables_to_plot_range(i_variable, strcmp(settings.variables_to_plot_range_header, 'variable_name'));
        variable_type = 'range';
        data_to_plot = data_custodian.getData(variable_name, variable_type);
        
        for i_comparison = 1 : length(comparisons.comparison_indices)
            % find correct condition indicator for control
            conditions_this_comparison = comparisons.comparison_indices{i_comparison};
            target_axes_handle = figure_data.axes_handles{figure_data.comparison_variable_to_axes_index_map(i_comparison), i_variable};
            
            % plot
            for i_condition = 1 : length(conditions_this_comparison)
                this_condition_index = conditions_this_comparison(i_condition);
                this_condition = comparisons.condition_combinations(this_condition_index, :);
                this_label = this_condition{strcmp(comparisons.condition_combination_labels, settings.condition_to_compare)};
                this_color = comparisons.condition_colors{strcmp(comparisons.condition_colors(:, 1), this_label), 2};
                label_string = strrep(this_label, '_', ' ');
                this_condition_indicator = getConditionIndicator ...
                  ( ...
                    this_condition, ...
                    comparisons.condition_combination_labels, ...
                    condition_data, ...
                    settings.condition_labels, ...
                    settings.plot_settings.get('levels_to_remove') ...
                  );
                data_to_plot_this_condition = data_to_plot.variable_data(:, this_condition_indicator);
                origin_subjects_this_condition = origin_subjects(this_condition_indicator);
                origin_trials_this_condition = origin_trials(this_condition_indicator);
                origin_start_times_this_condition = origin_start_times(this_condition_indicator);
                origin_end_times_this_condition = origin_end_times(this_condition_indicator);
                
                
                if settings.plot_settings.get('average_within_subjects', 1)
                    condition_data_this_condition = condition_data(this_condition_indicator, :);
%                     data_to_plot_this_condition = averageWithinSubjects(data_to_plot_this_condition, condition_data_this_condition, condition_labels);
                    
                    [data_to_plot_this_condition, origin_subjects_this_condition] ...
                        = averageWithinSubjects(data_to_plot_this_condition, condition_data_this_condition, condition_labels);
                end
                
                % build origin struct
                origin_data = [];
                for i_stretch = 1 : size(data_to_plot_this_condition, 2)
                    this_origin = struct;
                    this_origin.subject = origin_subjects_this_condition{i_stretch};
                    if ~settings.plot_settings.get('average_within_subjects', 1)
                        this_origin.trial = origin_trials_this_condition(i_stretch);
                        this_origin.start_time = origin_start_times_this_condition(i_stretch);
                        this_origin.end_time = origin_end_times_this_condition(i_stretch);
                    end
                    origin_data = [origin_data; this_origin];
                end                
                
                target_abscissa = figure_data.abscissae_cell{i_comparison, i_variable}(i_condition);
                if ~any(isnan(data_to_plot_this_condition))
                    if settings.group_bands_within_conditions
                        % override color
                        colors = copper(size(data_to_plot_this_condition, 1));
                        this_color = colors(i_band, :);
                    end
                    plotDiscreteData ...
                      ( ...
                        data_to_plot_this_condition, ...
                        'abscissa', target_abscissa, ...
                        'axes', target_axes_handle, ...         % axes
                        'color', this_color, ...                % color
                        'ShowMean', settings.show_average_data, ...
                        'MeanStyle', settings.plot_settings.get('mean_style', 1), ...
                        'MeanMarkerSize', settings.plot_settings.get('mean_marker_size', 1), ...
                        'MeanColor', lightenColor(this_color, 1-settings.plot_settings.get('mean_data_saturation', 1)), ...
                        'ShowMedian', settings.show_spread_data, ...
                        'MedianStyle', 'line', ...
                        'ShowIndividualData', settings.plot_settings.get('show_individual_discrete_data', 1), ...
                        'IndividualDataMarkerStyle', settings.plot_settings.get('individual_data_marker_style', 1), ...
                        'IndividualDataMarkerSize', settings.plot_settings.get('individual_data_marker_size', 1), ...
                        'IndividualDataColor', lightenColor(this_color, 1-settings.plot_settings.get('individual_data_saturation', 1)), ...
                        'IndividualDataJitterStd', settings.plot_settings.get('individual_data_jitter_std', 1), ...
                        'SpreadFaceAlpha', settings.plot_settings.get('discrete_data_spread_face_alpha', 1), ...
                        'SpreadFaceColor', this_color, ...
                        'SpreadEdgeColor', lightenColor(this_color, 1-settings.plot_settings.get('spread_edge_saturation', 1)), ...
                        'SpreadEdgeLinewidth', settings.plot_settings.get('discrete_data_spread_edge_linewidth', 1), ...
                        'SpreadWiskColor', lightenColor(this_color, 1-settings.plot_settings.get('spread_edge_saturation', 1)), ...
                        'ShowSpread', settings.show_spread_data, ...
                        'SpreadStyle', settings.plot_settings.get('discrete_data_plot_style'), ...
                        'IndividualDataInfo', origin_data, ...
                        'label', label_string ...     % label
                      );

                end
                
                
            end
            
            % update direction labels
            set(figure_data.pos_text_handles_ordinate(i_comparison, i_variable), 'string', data_to_plot.directions{1});
            set(figure_data.neg_text_handles_ordinate(i_comparison, i_variable), 'string', data_to_plot.directions{2});
        end
    end
end

function plotData_pairs(settings, comparisons, data_custodian, figure_data)
    [condition_data, condition_labels] = data_custodian.getConditionData();
    [origin_subjects, origin_trials, origin_start_times, origin_end_times] = data_custodian.getOriginData();
    for i_pair = 1 : settings.number_of_variables_to_plot_pairs
        variable_A_name = settings.variables_to_plot_pairs.variable_A_name{i_pair};
        variable_A_type = settings.variables_to_plot_pairs.variable_A_type{i_pair};
        variable_B_name = settings.variables_to_plot_pairs.variable_B_name{i_pair};
        variable_B_type = settings.variables_to_plot_pairs.variable_B_type{i_pair};
        data_to_plot_A = data_custodian.getData(variable_A_name, variable_A_type);
        data_to_plot_B = data_custodian.getData(variable_B_name, variable_B_type);
        start_index = settings.variables_to_plot_pairs.start_index{i_pair};
        end_index = settings.variables_to_plot_pairs.end_index{i_pair};
        if strcmp(start_index, '~')
            start_index = 1;
        else
            start_index = str2double(start_index);
        end
        if strcmp(end_index, '~')
            end_index = size(data_to_plot_A.variable_data, 1);
        else
            end_index = str2double(end_index);
        end
        
        for i_comparison = 1 : length(comparisons.comparison_indices)
            % find correct condition indicator for control
            conditions_this_comparison = comparisons.comparison_indices{i_comparison};
            top_level_plots = [];
            target_axes_handle = figure_data.axes_handles{figure_data.comparison_variable_to_axes_index_map(i_comparison), i_pair};
            
            % plot
            for i_condition = 1 : length(conditions_this_comparison)
                this_condition_index = conditions_this_comparison(i_condition);
                this_condition = comparisons.condition_combinations(this_condition_index, :);
                this_label = this_condition{strcmp(comparisons.condition_combination_labels, settings.condition_to_compare)};
                this_color = comparisons.condition_colors{strcmp(comparisons.condition_colors(:, 1), this_label), 2};
                label_string = strrep(this_label, '_', ' ');

                this_condition_indicator = getConditionIndicator ...
                  ( ...
                    this_condition, ...
                    comparisons.condition_combination_labels, ...
                    condition_data, ...
                    settings.condition_labels, ...
                    settings.plot_settings.get('levels_to_remove') ...
                  );
                data_to_plot_this_condition_A = data_to_plot_A.variable_data(start_index:end_index, this_condition_indicator);
                data_to_plot_this_condition_B = data_to_plot_B.variable_data(start_index:end_index, this_condition_indicator);
                origin_subjects_this_condition = origin_subjects(this_condition_indicator);
                origin_trials_this_condition = origin_trials(this_condition_indicator);
                origin_start_times_this_condition = origin_start_times(this_condition_indicator);
                origin_end_times_this_condition = origin_end_times(this_condition_indicator);
                
                if settings.plot_settings.get('average_within_subjects', 1)
                    error('Averaging within subjects not yet implemented for plotting variable pairs')
%                     condition_data_this_condition = condition_data(this_condition_indicator, :);
%                     [data_to_plot_this_condition, origin_subjects_this_condition] ...
%                         = averageWithinSubjects(data_to_plot_this_condition, condition_data_this_condition, condition_labels);
                end                
                
                if settings.plot_settings.get('show_individual_pairs_data', 1)
                    for i_stretch = 1 : size(data_to_plot_this_condition_A, 2)
                        this_origin = struct;
                        this_origin.subject = origin_subjects_this_condition{i_stretch};
                        if ~settings.plot_settings.get('average_within_subjects', 1)
                            this_origin.trial = origin_trials_this_condition(i_stretch);
                            this_origin.start_time = origin_start_times_this_condition(i_stretch);
                            this_origin.end_time = origin_end_times_this_condition(i_stretch);
                        end
                        plot ...
                          ( ...
                            target_axes_handle, ...
                            data_to_plot_this_condition_A(:, i_stretch), ...
                            data_to_plot_this_condition_B(:, i_stretch), ...
                            'linewidth', 1, ...
                            'HandleVisibility', 'off', ...
                            'color', lightenColor(this_color, 0.5), ...
                            'UserData', this_origin ...
                          );
                    end
                end
                if settings.show_average_data
                    average_plot = plot ...
                      ( ...
                        nanmean(data_to_plot_this_condition_A, 2), ...
                        nanmean(data_to_plot_this_condition_B, 2), ...
                        'parent', target_axes_handle, ...
                        'DisplayName', label_string, ...
                        'color', this_color, ...
                        'linewidth', 6 ...
                      );
                    top_level_plots = [top_level_plots average_plot]; %#ok<AGROW>
                end
                % TODO: order these different kinds of lines properly
                % across conditions, so that lines in one condition are not
                % un-clickable below spread patches from another condition
            end
            
            % update direction labels
            set(figure_data.pos_text_handles_abscissa(i_comparison, i_pair), 'string', data_to_plot_A.directions{1});
            set(figure_data.neg_text_handles_abscissa(i_comparison, i_pair), 'string', data_to_plot_A.directions{2});
            set(figure_data.pos_text_handles_ordinate(i_comparison, i_pair), 'string', data_to_plot_B.directions{1});
            set(figure_data.neg_text_handles_ordinate(i_comparison, i_pair), 'string', data_to_plot_B.directions{2});
        end
    end
end

function groomFigures_continuous(settings, data_custodian, comparisons, figure_data)
    % set axis limits
    dictateAxes(settings, figure_data, settings.variables_to_plot_continuous, settings.variables_to_plot_continuous_header);
    updateLabelPositions(figure_data);
    addZeroLine(settings, figure_data);

    % toggle legend
    if settings.show_legend
        for i_variable = 1 : settings.number_of_variables_to_plot_continuous
            for i_axes = 1 : size(figure_data.axes_handles, 1)
                these_axes = figure_data.axes_handles{i_axes, i_variable};
                legend(these_axes, 'show')
            end
        end
    end
    
    % mark bands
    if settings.mark_bands
        for i_comparison = 1 : comparisons.number_of_comparisons
            for i_variable = 1 : settings.number_of_variables_to_plot_continuous
                these_axes = figure_data.axes_handles{i_comparison, i_variable};
                these_abscissae = figure_data.abscissae_cell{i_comparison, i_variable};
                ylimits = get(these_axes, 'ylim');

                if settings.mark_bands == 1
                    bands_to_mark = 2 : 2 : data_custodian.bands_per_stretch;
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
    
    % mark pushoff
    if settings.mark_pushoff
        % determine pushoff times ad indices
        pushoff_time_data = data_custodian.getData('pushoff_time', 'stretch');
        pushoff_times = pushoff_time_data.variable_data;
        step_time_data = data_custodian.getData('step_time', 'stretch');
        step_times = step_time_data.variable_data;
        pushoff_time_ratio = pushoff_times ./ step_times;
        mean_pushoff_ratio = mean(pushoff_time_ratio, 2);
        pushoff_index = round(mean_pushoff_ratio * 100);
        
        for i_comparison = 1 : comparisons.number_of_comparisons
            for i_variable = 1 : settings.number_of_variables_to_plot_continuous
                these_axes = figure_data.axes_handles{i_comparison, i_variable};
                these_abscissae = figure_data.abscissae_cell{i_comparison, i_variable};
                ylimits = get(these_axes, 'ylim');

                for i_band = 1 : data_custodian.bands_per_stretch
                    % double stance patch
                    double_stance_patch_color = settings.plot_settings.get('stance_double_color');

                    start_index = getBandIndices(i_band, settings.number_of_time_steps_normalized);
                    pushoff_index_here = start_index + pushoff_index(i_band);

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

function groomFigures_discrete(settings, figure_data)
    % set axis limits
    dictateAxes(settings, figure_data, settings.variables_to_plot_discrete, settings.variables_to_plot_discrete_header);
    updateLabelPositions(figure_data);
    addZeroLine(settings, figure_data);
    
    % rotate labels
    for i_variable = 1 : settings.number_of_variables_to_plot_discrete
        for i_axes = 1 : size(figure_data.axes_handles, 1)
            these_axes = figure_data.axes_handles{i_axes, i_variable};
            xtick_label_rotation = settings.plot_settings.get('xtick_label_rotation', 1);
            set(these_axes, 'XTickLabelRotation', xtick_label_rotation);            
        end
    end

    
end

function groomFigures_range(settings, figure_data)
    % set axis limits
    dictateAxes(settings, figure_data, settings.variables_to_plot_range, settings.variables_to_plot_range_header);
    updateLabelPositions(figure_data);
    addZeroLine(settings, figure_data);
    
    % rotate labels
    for i_variable = 1 : settings.number_of_variables_to_plot_range
        for i_axes = 1 : size(figure_data.axes_handles, 1)
            these_axes = figure_data.axes_handles{i_axes, i_variable};
            xtick_label_rotation = settings.plot_settings.get('xtick_label_rotation', 1);
            set(these_axes, 'XTickLabelRotation', xtick_label_rotation);            
        end
    end

    
end

function groomFigures_pairs(settings, figure_data)
    % set axis limits
    dictateAxes(settings, figure_data, settings.variables_to_plot_pairs);
    updateLabelPositions(figure_data);
    addZeroLine(settings, figure_data, 'x');
    addZeroLine(settings, figure_data, 'y');
    
    % rotate labels
    for i_variable = 1 : settings.number_of_variables_to_plot_pairs
        for i_axes = 1 : size(figure_data.axes_handles, 1)
            these_axes = figure_data.axes_handles{i_axes, i_variable};
            xtick_label_rotation = settings.plot_settings.get('xtick_label_rotation', 1);
            set(these_axes, 'XTickLabelRotation', xtick_label_rotation);            
        end
    end

    
end

function dictateAxes(settings, figure_data, variables_to_plot, variables_to_plot_header)
    if settings.dictate_axes
        for i_variable = 1 : size(variables_to_plot, 1)
            if nargin == 3
                % no header, assume variables_to_plot is a table
                this_variable_x_lower = variables_to_plot.x_axis_lower_limit{i_variable};
                this_variable_x_upper = variables_to_plot.x_axis_upper_limit{i_variable};
                this_variable_y_lower = variables_to_plot.y_axis_lower_limit{i_variable};
                this_variable_y_upper = variables_to_plot.y_axis_upper_limit{i_variable};
            else
                % get x-axis limits from settings
                this_variable_x_lower = variables_to_plot{i_variable, strcmp(variables_to_plot_header, 'x_axis_lower_limit')};
                this_variable_x_upper = variables_to_plot{i_variable, strcmp(variables_to_plot_header, 'x_axis_upper_limit')};

                % get y-axis limits from settings
                this_variable_y_lower = variables_to_plot{i_variable, strcmp(variables_to_plot_header, 'y_axis_lower_limit')};
                this_variable_y_upper = variables_to_plot{i_variable, strcmp(variables_to_plot_header, 'y_axis_upper_limit')};
            end

            for i_axes = 1 : size(figure_data.axes_handles, 1)
                % get current axes and limits
                these_axes = figure_data.axes_handles{i_axes, i_variable};
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
end

function updateLabelPositions(figure_data)
    % update label positions
    for i_variable = 1 : size(figure_data.axes_handles, 2)
        for i_axes = 1 : size(figure_data.axes_handles, 1)
            these_axes = figure_data.axes_handles{i_axes, i_variable};
            xlimits = get(these_axes, 'xlim'); ylimits = get(these_axes, 'ylim');
            if figure_data.pos_arrow_handles_abscissa(i_axes, i_variable) ~= 0
                pos_arrow_position_x = xlimits(2);
                pos_arrow_position_y = ylimits(1) - (ylimits(2)-ylimits(1))*0.09;
                set(figure_data.pos_arrow_handles_abscissa(i_axes, i_variable), 'Position', [pos_arrow_position_x pos_arrow_position_y]);
                pos_text_position_x = xlimits(2);
                pos_text_position_y = ylimits(1) - (ylimits(2)-ylimits(1))*0.14;
                set(figure_data.pos_text_handles_abscissa(i_axes, i_variable), 'Position', [pos_text_position_x pos_text_position_y]);
            end
            if figure_data.neg_arrow_handles_abscissa(i_axes, i_variable) ~= 0
                neg_arrow_position_x = xlimits(1);
                neg_arrow_position_y = ylimits(1) - (ylimits(2)-ylimits(1))*0.09;
                set(figure_data.neg_arrow_handles_abscissa(i_axes, i_variable), 'Position', [neg_arrow_position_x neg_arrow_position_y]);
                neg_text_position_x = xlimits(1);
                neg_text_position_y = ylimits(1) - (ylimits(2)-ylimits(1))*0.14;
                set(figure_data.neg_text_handles_abscissa(i_axes, i_variable), 'Position', [neg_text_position_x neg_text_position_y]);
            end            
            if figure_data.pos_arrow_handles_ordinate(i_axes, i_variable) ~= 0
                pos_arrow_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.09;
                pos_arrow_position_y = ylimits(2);
                set(figure_data.pos_arrow_handles_ordinate(i_axes, i_variable), 'Position', [pos_arrow_position_x pos_arrow_position_y]);
                pos_text_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.14;
                pos_text_position_y = ylimits(2);
                set(figure_data.pos_text_handles_ordinate(i_axes, i_variable), 'Position', [pos_text_position_x pos_text_position_y]);
            end
            if figure_data.neg_arrow_handles_ordinate(i_axes, i_variable) ~= 0
                neg_arrow_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.09;
                neg_arrow_position_y = ylimits(1);
                set(figure_data.neg_arrow_handles_ordinate(i_axes, i_variable), 'Position', [neg_arrow_position_x neg_arrow_position_y]);
                neg_text_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.14;
                neg_text_position_y = ylimits(1);
                set(figure_data.neg_text_handles_ordinate(i_axes, i_variable), 'Position', [neg_text_position_x neg_text_position_y]);
            end            
        end
    end
end

function addZeroLine(settings, figure_data, axis)
    if nargin < 3
        axis = 'x';
    end
    if settings.plot_settings.get('plot_zero', 1)
        for i_variable = 1 : size(figure_data.axes_handles, 2)
            for i_axes = 1 : size(figure_data.axes_handles, 1)
                these_axes = figure_data.axes_handles{i_axes, i_variable};
                if ~isempty(these_axes) && axis == 'x'
                    xlimits = get(these_axes, 'xlim');
                    zero_plot = plot(these_axes, xlimits, [0 0], 'Color', settings.plot_settings.get('zero_color', 1), 'linewidth', settings.plot_settings.get('zero_linewidth', 1));
%                     zero_plot = yline(these_axes, 0, 'Color', settings.plot_settings.get('zero_color', 1), 'linewidth', settings.plot_settings.get('zero_linewidth', 1));
                end
                if ~isempty(these_axes) && axis == 'y'
%                     ylimits = get(these_axes, 'ylim');
                    zero_plot = xline(these_axes, 0, 'Color', settings.plot_settings.get('zero_color', 1), 'linewidth', settings.plot_settings.get('zero_linewidth', 1));
                end
                set(zero_plot, 'HandleVisibility', 'off');
                uistack(zero_plot, 'bottom')
                
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
        for i_figure = 1 : numel(figure_data.figure_handles)
            % remove some white space on right side and top
            axes_position = get(figure_data.axes_handles{i_figure}, 'position');
            axes_position(3) = 1 - axes_position(1) - 0.01;
            axes_position(4) = 1 - axes_position(2) - 0.04;
            set(figure_data.axes_handles{i_figure}, 'position', axes_position);
            
            % save with labels
            filename_with = ['figures' filesep 'withLabels' filesep get(figure_data.figure_handles{i_figure}, 'UserData')];
            print(figure_data.figure_handles{i_figure}, filename_with, settings.save_format, settings.save_resolution)
            
            % remove text and marks to save data lines only
            set(get(figure_data.axes_handles{i_figure}, 'xaxis'), 'visible', 'off');
            set(get(figure_data.axes_handles{i_figure}, 'yaxis'), 'visible', 'off');
            set(get(figure_data.axes_handles{i_figure}, 'xlabel'), 'visible', 'off');
            set(get(figure_data.axes_handles{i_figure}, 'ylabel'), 'visible', 'off');
            set(get(figure_data.axes_handles{i_figure}, 'title'), 'visible', 'off');
            set(figure_data.axes_handles{i_figure}, 'xticklabel', '');
            set(figure_data.axes_handles{i_figure}, 'yticklabel', '');
            set(figure_data.axes_handles{i_figure}, 'position', [0 0 1 1]);
            legend(figure_data.axes_handles{i_figure}, 'hide');
            filename_without = ['figures' filesep 'noLabels' filesep get(figure_data.figure_handles{i_figure}, 'UserData')];
            print(figure_data.figure_handles{i_figure}, filename_without, settings.save_format, settings.save_resolution)
            disp(['Saved as ' filename_with ' and ' filename_without])
            
            % put some marks back
            set(get(figure_data.axes_handles{i_figure}, 'title'), 'visible', 'on');
            set(figure_data.axes_handles{i_figure}, 'position', [0.05 0.05 0.9 0.9]);
        end
    end

end

function closeFigures(settings, figure_data)
    % close figures
    if settings.close
        for i_figure = 1 : numel(figure_data.figure_handles)
            close(figure_data.figure_handles{i_figure})            
        end
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

function [data_averaged, subjects_unique] = averageWithinSubjects(data_to_plot_this_condition, condition_data_this_condition, condition_labels)
    % extract subjects
    subjects = condition_data_this_condition(:, strcmp(condition_labels, 'subject'));
    subjects_unique = unique(subjects);
    number_of_unique_subjects = numel(subjects_unique);
    
    % calculate average within each subject
    data_averaged = zeros(size(data_to_plot_this_condition, 1), number_of_unique_subjects);
    for i_subject = 1 : number_of_unique_subjects
        this_subject = subjects_unique{i_subject};
        this_subject_indicator = strcmp(subjects, this_subject);
        this_subject_data = data_to_plot_this_condition(:, this_subject_indicator);
        this_subject_data_mean = mean(this_subject_data, 2);
        data_averaged(:, i_subject) = this_subject_data_mean;
    end
    
    
end











