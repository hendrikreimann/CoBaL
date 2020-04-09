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
%   a cell array "data.comparison_indices", where each entry is an array of condition indices, referring to a row in 
%   condition_combination_labels. Usually there will be one comparison per figure. When plotting episodes of multiple
%   consecutive steps, however, there will be multiple comparisons per figure. To keep track of this, the code generates
%   an array "figure_data.trajectory_axes_handles" that maps comparisons and variables to axes handles, and a cell array 
%   "abscissae_cell" storing the x-values. For both these arrays, rows = comparisons, columns = variables

% An episode is a list of comparisons
%   (this is a bit of a mess and will be fixed eventually... I promise - HR)
%   2020-APR-06 HR: fixing this now, fully removing the concept of episode

function plotResults(varargin)
    settings = determineSettings(varargin{:});

    % load data
    data = loadDataToPlot(settings);
    
    % determine condition combinations to plot and comparisons
    data = createComparisonData(settings, data);
    
    % create figures and determine abscissae for each comparison
    figure_data = createFigureData(settings, data);      
    
    % plot data
    figure_data = plotData(settings, data, figure_data);
    
    % groom axes, labels etc
    figure_data = groomFigures(settings, data, figure_data);
    
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
    addParameter(parser, 'format', 'tiff')
    addParameter(parser, 'settings', 'plot')
    addParameter(parser, 'spread_method', 'cinv')
    parse(parser, varargin{:})
    settings.subjects = parser.Results.subjects;
    settings.dictate_axes = parser.Results.dictate_axes;
    settings.show_legend = parser.Results.show_legend;
    settings.settings_file = parser.Results.settings;
    settings.spread_method = parser.Results.spread_method;
    settings.save_results = parser.Results.save;
    settings.close = parser.Results.close;

    % load settings
    study_settings = loadSettingsFromFile('study');
    plot_settings = loadSettingsFromFile(settings.settings_file);
    
    settings.show_outliers = plot_settings.get('show_outliers');
    settings.show_single_data_points = plot_settings.get('show_single_data_points', 1);
    settings.edge_color = plot_settings.get('edge_color', 1);
    settings.plot_mode = plot_settings.get('plot_mode');
    settings.mark_pushoff = plot_settings.get('mark_pushoff', 1);
    settings.mark_bands = plot_settings.get('mark_bands', 1);
    settings.band_labels = study_settings.get('band_labels', 1);
    settings.group_bands_within_conditions = plot_settings.get('group_bands_within_conditions', 1);
    settings.number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    settings.show_individual_data = plot_settings.get('show_individual_data', 1);
    settings.show_average_data = plot_settings.get('show_average_data', 1);
    settings.show_spread_data = plot_settings.get('show_spread_data', 1);
    
    settings.study_settings = study_settings;
    settings.plot_settings = plot_settings;
    
    % colors
    settings.colors_comparison = settings.plot_settings.get('colors_comparison');
    if size(settings.colors_comparison, 2) == 1
        settings.colors_comparison = hex2rgb(colors_comparison);
    end
    settings.colors_bands = settings.plot_settings.get('colors_bands', 1);
    settings.color_control = settings.plot_settings.get('color_control', 1);
end

function data = loadDataToPlot(settings)
    data = struct;

    % declare variables
    data.conditions_settings = settings.study_settings.get('conditions');
    data.condition_to_compare = settings.plot_settings.get('condition_to_compare');
    data.condition_labels = data.conditions_settings(:, 1)';
    condition_source_variables = data.conditions_settings(:, 2)';
    number_of_condition_labels = length(data.condition_labels);
    data_source = settings.plot_settings.get('data_source', 1);
    file_label = ['results' data_source];
    
    % load data
    data_folder_list = determineDataStructure(settings.subjects);
    data.variables_to_plot = settings.plot_settings.get('variables_to_plot');
    variables_to_plot_header = settings.plot_settings.get('variables_to_plot_header', true);

    data.number_of_variables_to_plot = size(data.variables_to_plot, 1);
    data.condition_data_all = {};
    origin_trial_list_all = [];
    origin_start_time_list_all = [];
    origin_end_time_list_all = [];
    data.data_all = cell(data.number_of_variables_to_plot, 1);
    data.directions = cell(data.number_of_variables_to_plot, 2);
    
    data.step_time_data = [];
    pushoff_time_data = [];
    data.bands_per_stretch = [];
    
    for i_folder = 1 : length(data_folder_list)
        % get information
        this_data_folder_path = data_folder_list{i_folder};
        load([this_data_folder_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
        
        % find results file
        results_file_candidate_analysis = [this_data_folder_path filesep 'analysis' filesep makeFileName(date, subject_id, file_label) '.mat'];
        results_file_candidate_subject = [this_data_folder_path filesep makeFileName(date, subject_id, file_label) '.mat'];
        results_file_candidate_results = [this_data_folder_path filesep 'results' filesep  makeFileName(date, subject_id, file_label) '.mat'];
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
        % HR, 2.10.2019 -- TF added this, but I don't think it's needed. Remove for now
%         if settings.study_settings.get('gather_step_minus_one')
%             bands_per_stretch_this_session = bands_per_stretch_this_session+1;
%         end

        % transform conditions into cell array
        conditions_session = loaded_data.conditions_session;
        condition_array_session = cell(number_of_stretches_this_session, number_of_condition_labels);
        for i_condition = 1 : number_of_condition_labels
            condition_array_session(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
        end
        
        data.condition_data_all = [data.condition_data_all; condition_array_session];
        origin_trial_list_all = [origin_trial_list_all; loaded_data.origin_trial_list_session]; %#ok<AGROW>
        origin_start_time_list_all = [origin_start_time_list_all; loaded_data.origin_start_time_list_session]; %#ok<AGROW>
        origin_end_time_list_all = [origin_end_time_list_all; loaded_data.origin_end_time_list_session]; %#ok<AGROW>
        
        % extract data
        for i_variable = 1 : data.number_of_variables_to_plot
            
            this_variable_name = data.variables_to_plot{i_variable, strcmp(variables_to_plot_header, 'variable name')};
            this_variable_type = data.variables_to_plot{i_variable, strcmp(variables_to_plot_header, 'variable type')};
            this_variable_source_index = find(strcmp(loaded_data.([this_variable_type '_names_session']), this_variable_name), 1, 'first');
            if isempty(this_variable_source_index)
                error(['Variable not found: ' this_variable_name])
            end
            this_variable_data = loaded_data.([this_variable_type '_data_session']){this_variable_source_index};
            this_variable_directions = loaded_data.([this_variable_type '_directions_session'])(this_variable_source_index, :);


            if settings.plot_settings.get('convert_to_mm') && (strcmp(this_variable_name,'cop_from_com_x') || strcmp(this_variable_name, 'step_placement_x'))
                this_variable_data = this_variable_data * 1000;
            end
            
            % store
            data.data_all{i_variable} = [data.data_all{i_variable} this_variable_data];
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

function data = createComparisonData(settings, data)
    labels_to_ignore = settings.plot_settings.get('conditions_to_ignore');
    levels_to_remove = settings.plot_settings.get('levels_to_remove');
    preferred_level_order = settings.plot_settings.get('preferred_level_order', 1);
    [data.condition_combination_labels, data.condition_combinations_stimulus, data.condition_combinations_control] = determineConditionCombinations(data.condition_data_all, data.conditions_settings, labels_to_ignore, levels_to_remove);
    [data.condition_combinations_stimulus, data.condition_combinations_control] = sortConditionCombinations(data.condition_combinations_stimulus, data.condition_combinations_control, data.condition_combination_labels, data.condition_to_compare, preferred_level_order);
    [data.comparison_indices, data.conditions_per_comparison_max] = determineComparisons(data.condition_combinations_stimulus, data.condition_combination_labels, settings.plot_settings);
    data.number_of_comparisons = length(data.comparison_indices);
end

function figure_data = createFigureData(settings, data)
    figure_data.comparison_variable_to_axes_index_map = zeros(data.number_of_comparisons, 1);
    figure_data.abscissae_cell = cell(data.number_of_comparisons, data.number_of_variables_to_plot);
    
    % time plots
    if strcmp(settings.plot_mode, 'detailed') || strcmp(settings.plot_mode, 'overview')
        % make one figure per comparison and variable
        figure_data.trajectory_figure_handles = zeros(data.number_of_comparisons, data.number_of_variables_to_plot);
        figure_data.trajectory_axes_handles = zeros(data.number_of_comparisons, data.number_of_variables_to_plot);
        figure_data.pos_text_handles = zeros(data.number_of_comparisons, data.number_of_variables_to_plot);
        figure_data.neg_text_handles = zeros(data.number_of_comparisons, data.number_of_variables_to_plot);
        figure_data.pos_arrow_handles = zeros(data.number_of_comparisons, data.number_of_variables_to_plot);
        figure_data.neg_arrow_handles = zeros(data.number_of_comparisons, data.number_of_variables_to_plot);
        step_start_times_cell = cell(data.number_of_comparisons, data.number_of_variables_to_plot);
        step_end_times_cell = cell(data.number_of_comparisons, data.number_of_variables_to_plot);
        for i_variable = 1 : data.number_of_variables_to_plot
            for i_comparison = 1 : data.number_of_comparisons
                this_comparison = data.comparison_indices{i_comparison};
                % make figure and axes
                new_figure = figure; new_axes = axes; hold on;
                
                % store handles and determine abscissa data
                figure_data.trajectory_figure_handles(i_comparison, i_variable) = new_figure;
                figure_data.trajectory_axes_handles(i_comparison, i_variable) = new_axes;
                figure_data.comparison_variable_to_axes_index_map(i_comparison) = i_comparison;
                    
                if isDiscreteVariable(i_variable, data.data_all, data.bands_per_stretch)
                    % abscissa gives the bin edges here
                    data_to_plot = data.data_all{i_variable, 1};
                    if settings.dictate_axes
                        lower_bound = str2double(data.variables_to_plot{i_variable, 6});
                        upper_bound = str2double(data.variables_to_plot{i_variable, 7});
                    else
                        lower_bound = min(data_to_plot);
                        upper_bound = max(data_to_plot);
                    end
                    if strcmp(settings.plot_mode, 'detailed')
%                         figure_data.abscissae_cell{i_comparison, i_variable}(i_condition, :) = linspace(lower_bound, upper_bound, settings.plot_settings.get('number_of_bins_in_histogram'));
                        figure_data.abscissae_cell{i_comparison, i_variable} = linspace(lower_bound, upper_bound, settings.plot_settings.get('number_of_bins_in_histogram'));
                    end
                    if strcmp(settings.plot_mode, 'overview')
                        this_comparison = data.comparison_indices{i_comparison};
                        number_of_entries = length(this_comparison);
                        if settings.plot_settings.get('plot_control')
                            number_of_entries = number_of_entries + 1;
                        end
                        
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
                            shifter = (0:data.bands_per_stretch-1)' * (data.conditions_per_comparison_max + gap_between_bands);
                            abscissae_stimulus = abscissae_stimulus + repmat(shifter, 1, data.conditions_per_comparison_max);
                            if settings.plot_settings.get('merge_bands', 1)
                                abscissae_stimulus = abscissae_stimulus(1, :);
                            end
                        end
                        figure_data.abscissae_cell{i_comparison, i_variable} = abscissae_stimulus;
                        
                    end
                end
                if isContinuousVariable(i_variable, data.data_all, data.bands_per_stretch)
                    % scale abscissae
                    conditions_this_comparison = data.comparison_indices{i_comparison};
                    step_time_means_this_comparison = zeros(data.bands_per_stretch, size(conditions_this_comparison, 2));
                    for i_condition = 1 : length(conditions_this_comparison)
                        this_condition_combination = data.condition_combinations_stimulus(conditions_this_comparison(i_condition), :);
                        this_condition_indicator = getConditionIndicator(this_condition_combination, data.condition_combination_labels, data.condition_data_all, data.condition_labels);
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
                if settings.dictate_axes && ~(strcmp(settings.plot_mode, 'detailed') && isDiscreteVariable(i_variable, data.data_all, data.bands_per_stretch))
    %                 set(gca, 'xlim', [time_normalized(1), time_normalized(end)]);
                    set(gca, 'ylim', [str2double(data.variables_to_plot{i_variable, 6}), str2double(data.variables_to_plot{i_variable, 7})]);
                end
                if isDiscreteVariable(i_variable, data.data_all, data.bands_per_stretch) && strcmp(settings.plot_mode, 'overview')
%                     xtick = figure_data.abscissae_cell{i_comparison, i_variable}{2};
%                     if settings.plot_settings.get('plot_control')
%                         xtick = [figure_data.abscissae_cell{i_comparison, i_variable}{1} xtick]; %#ok<AGROW>
%                     end
                    xtick = sort(reshape(figure_data.abscissae_cell{i_comparison, i_variable}, 1, numel(figure_data.abscissae_cell{i_comparison, i_variable})));
                    set(gca, 'xlim', [-0.5 + min(xtick) 0.5 + max(xtick(end))]);
                    set(gca, 'xtick', xtick);
                end
                
                % set axis labels
                if isContinuousVariable(i_variable, data.data_all, data.bands_per_stretch)
                    if strcmp(settings.plot_settings.get('time_plot_style'), 'scaled_to_comparison_mean') || strcmp(settings.plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                        xlabel('normalized time (s)');
                    else
                        xlabel('normalized time (%)');
                    end
                end
                ylabel(data.variables_to_plot{i_variable, 4});
                
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
                title_string = data.variables_to_plot{i_variable, 3};
                filename_string = data.variables_to_plot{i_variable, 5};
                
                representative_condition = data.condition_combinations_stimulus(this_comparison(1), :);
                
                for i_label = 1 : length(representative_condition)
                    if ~(strcmp(data.condition_combination_labels{i_label}, data.condition_to_compare))
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
        for i_variable = 1 : data.number_of_variables_to_plot
            if isContinuousVariable(i_variable, data.data_all, data.bands_per_stretch)
                for i_comparison = 1 : data.number_of_comparisons
                    target_abscissae = figure_data.abscissae_cell{i_comparison, i_variable};
                    
                    % HR, 2.10.2019 -- TF added this, but I don't think it's needed. Remove for now
%                     if settings.plot_settings.get('cutoff_2nd_doublestance')
%                         xlim = [min(target_abscissae(:, 1)) 100 + data.pushoff_index(2)];
%                     else
                        xlim = [min(target_abscissae(:, 1)) max(target_abscissae(:, end))];
%                     end
                    % set x-limits accordingly
                    set(figure_data.trajectory_axes_handles(i_comparison, i_variable), 'xlim', xlim);
                end
            end
        end

        % determine stance start and end times and stance foot
        for i_variable = 1 : data.number_of_variables_to_plot
            if isContinuousVariable(i_variable, data.data_all, data.bands_per_stretch)
                for i_comparison = 1 : data.number_of_comparisons
                    
                    step_abscissa = figure_data.abscissae_cell{i_comparison, i_variable};
                    step_start_times_cell{i_comparison, i_variable} = step_abscissa(1, 1);
                    step_end_times_cell{i_comparison, i_variable} = step_abscissa(1, end);
                    
%                     if mark_pushoff
%                         step_pushoff_times_cell{i_comparison, i_variable} = step_abscissa(1, data.pushoff_index);
%                     end
                end
            end
        end        
    end
end

function figure_data = plotData(settings, data, figure_data)
    for i_variable = 1 : data.number_of_variables_to_plot
        data_to_plot = data.data_all{i_variable, 1};
        for i_comparison = 1 : length(data.comparison_indices)
            % find correct condition indicator for control
            conditions_this_comparison = data.comparison_indices{i_comparison};
            top_level_plots = [];
            target_axes_handle = figure_data.trajectory_axes_handles(figure_data.comparison_variable_to_axes_index_map(i_comparison), i_variable);
            
            % plot control
            if settings.plot_settings.get('plot_control')
                % determine which control condition applies here
                representant_condition_index = conditions_this_comparison(1);
                this_condition = data.condition_combinations_control(representant_condition_index, :);
                
                this_condition_indicator = getConditionIndicator(this_condition, data.condition_combination_labels, data.condition_data_all, data.condition_labels);
                data_to_plot_this_condition = data_to_plot(:, this_condition_indicator);
                origin_indices = find(this_condition_indicator);
                
                if ~isempty(data_to_plot_this_condition)
                    if isDiscreteVariable(i_variable, data.data_all, data.bands_per_stretch)
                        if settings.plot_settings.get('merge_bands', 1)
                            data_to_plot_this_condition = reshape(data_to_plot_this_condition, 1, numel(data_to_plot_this_condition));
                        end
                        for i_band = 1 : size(data_to_plot_this_condition, 1)
                            if ~isempty(settings.band_labels)
                                label_string_this_band = ['control -' settings.band_labels{i_band}];
                            else
                                label_string_this_band = 'control';
                            end
                            target_abscissa = figure_data.abscissae_cell{i_comparison, i_variable}(i_band, end);
                            data_to_plot_this_band = data_to_plot_this_condition(i_band, :);
                            if strcmp(settings.plot_mode, 'overview')
                                if ~any(isnan(data_to_plot_this_band))
                                    if settings.group_bands_within_conditions
                                        this_color = settings.colors_bands(i_band, :);
                                    else
                                        this_color = settings.plot_settings.get('color_control');
                                    end

                                    plotDiscreteData ...
                                      ( ...
                                        data_to_plot_this_band, ...
                                        'SpreadStyle', settings.plot_settings.get('discrete_data_plot_style'), ...
                                        'abscissa', target_abscissa, ...
                                        'axes', target_axes_handle, ...         % axes
                                        'color', this_color, ...                % color
                                        'ShowMean', true, ...
                                        'MeanStyle', 'd', ...
                                        'MeanColor', [1 1 1]*0.7, ...
                                        'ShowMedian', true, ...
                                        'MedianStyle', 'line', ...
                                        'label', label_string_this_band ...    % label
                                      )

                                end
                            end
                        end
                    end
                    if isContinuousVariable(i_variable, data.data_all, data.bands_per_stretch)
                        target_abscissa = figure_data.abscissae_cell{i_comparison, i_variable}(end, :);
                        if settings.show_individual_data
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
                                    'color', lightenColor(settings.color_control, 0.5) ...
                                  );
                            end

                        end
                        if settings.show_average_data
                            average_plot = plot ...
                              ( ...
                                target_abscissa, ...
                                nanmean(data_to_plot_this_condition, 2), ...
                                'parent', target_axes_handle, ...
                                'DisplayName', 'CONTROL', ...
                                'color', settings.color_control, ...
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
                                  'color', settings.color_control, ...
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
            
            % plot stimulus
            for i_condition = 1 : length(conditions_this_comparison)
                this_condition_index = conditions_this_comparison(i_condition);
                this_condition = data.condition_combinations_stimulus(this_condition_index, :);
                label_string = strrep(this_condition{strcmp(data.condition_combination_labels, data.condition_to_compare)}, '_', ' ');
                this_condition_indicator = getConditionIndicator(this_condition, data.condition_combination_labels, data.condition_data_all, data.condition_labels);
                data_to_plot_this_condition = data_to_plot(:, this_condition_indicator);
                
                origin_indices = find(this_condition_indicator);
                if isDiscreteVariable(i_variable, data.data_all, data.bands_per_stretch)
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
                                this_color = settings.colors_bands(i_band, :);
                            else
                                this_color = settings.colors_comparison(i_condition, :);
                            end
                            plotDiscreteData ...
                              ( ...
                                data_to_plot_this_band, ...
                                'SpreadStyle', settings.plot_settings.get('discrete_data_plot_style'), ...
                                'abscissa', target_abscissa, ...
                                'axes', target_axes_handle, ...         % axes
                                'color', this_color, ...                % color
                                'ShowMean', true, ...
                                'MeanStyle', 'd', ...
                                'MeanColor', [1 1 1]*0.7, ...
                                'ShowMedian', true, ...
                                'MedianStyle', 'line', ...
                                'label', label_string_this_band ...    % label
                              )
                            
                        end
                    end
                end
                if isContinuousVariable(i_variable, data.data_all, data.bands_per_stretch)
                    target_abscissa = figure_data.abscissae_cell{i_comparison, i_variable}(i_condition, :);                    
                    if settings.show_individual_data
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
                                'color', lightenColor(settings.colors_comparison(i_condition, :), 0.5) ...
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
                            'color', settings.colors_comparison(i_condition, :), ...
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
                              'color', settings.colors_comparison(i_condition, :), ...
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

function figure_data = groomFigures(settings, data, figure_data)
    % update label positions
    for i_variable = 1 : data.number_of_variables_to_plot
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
            if isDiscreteVariable(i_variable, data.data_all, data.bands_per_stretch)
                % rotate labels
                xtick_label_rotation = settings.plot_settings.get('xtick_label_rotation', 1);
                set(these_axes, 'XTickLabelRotation', xtick_label_rotation);
                
            end
            
        end
    end

    % toggle legend
    for i_variable = 1 : data.number_of_variables_to_plot
        for i_axes = 1 : size(figure_data.trajectory_axes_handles, 1)
            these_axes = figure_data.trajectory_axes_handles(i_axes, i_variable);
            if settings.show_legend && ~(isDiscreteVariable(i_variable, data.data_all, data.bands_per_stretch) && (strcmp(settings.plot_mode, 'overview')))
                legend(these_axes, 'show')
            end
        end
    end
    
    % add zero line
    if settings.plot_settings.get('plot_zero', 1)
        for i_variable = 1 : data.number_of_variables_to_plot
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
        for i_comparison = 1 : data.number_of_comparisons
            for i_variable = 1 : data.number_of_variables_to_plot
                if isContinuousVariable(i_variable, data.data_all, data.bands_per_stretch)
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
                                patch_x = [band_start_times(i_condition) band_end_times(i_condition) band_end_times(i_condition) band_start_times(i_condition)];
                                patch_y = [y_values(i_condition) y_values(i_condition) y_values(i_condition+1) y_values(i_condition+1)];
                                patch_handle = ...
                                    patch ...
                                      ( ...
                                        patch_x, ...
                                        patch_y, ...
                                        settings.colors_comparison(i_condition, :), ...
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
        for i_comparison = 1 : data.number_of_comparisons
            for i_variable = 1 : data.number_of_variables_to_plot
                if isContinuousVariable(i_variable, data.data_all, data.bands_per_stretch)
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
            % save with labels
%             legend(axes_handles(i_figure), 'show');
            filename = ['figures' filesep 'withLabels' filesep get(figure_data.trajectory_figure_handles(i_figure), 'UserData')];
            saveas(figure_data.trajectory_figure_handles(i_figure), filename, parser.Results.format)
            
            % save without labels
%             set(postext, 'visible', 'off');
%             set(negtext, 'visible', 'off');
            
            % remove text and marks to save graphs only
            set(get(figure_data.trajectory_axes_handles(i_figure), 'xaxis'), 'visible', 'off');
            set(get(figure_data.trajectory_axes_handles(i_figure), 'yaxis'), 'visible', 'off');
            set(get(figure_data.trajectory_axes_handles(i_figure), 'xlabel'), 'visible', 'off');
            set(get(figure_data.trajectory_axes_handles(i_figure), 'ylabel'), 'visible', 'off');
            set(get(figure_data.trajectory_axes_handles(i_figure), 'title'), 'visible', 'off');
            set(figure_data.trajectory_axes_handles(i_figure), 'xticklabel', '');
            set(figure_data.trajectory_axes_handles(i_figure), 'yticklabel', '');
            set(figure_data.trajectory_axes_handles(i_figure), 'position', [0 0 1 1]);
            legend(figure_data.trajectory_axes_handles(i_figure), 'hide');
            filename = ['figures' filesep 'noLabels' filesep get(figure_data.trajectory_figure_handles(i_figure), 'UserData')];
            saveas(figure_data.trajectory_figure_handles(i_figure), filename, parser.Results.format);
            
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

% function plotDiscreteData(xdata, ydata, settings, varargin)
%     parser = inputParser;
%     parser.KeepUnmatched = true;
%     addParameter(parser, 'axes', gca);
%     addParameter(parser, 'color', [0.5 0.5 0.5]);
%     addParameter(parser, 'label', '');
%     
%     parse(parser, varargin{:})
%     axes_handle = parser.Results.axes;
%     color = parser.Results.color;
%     label = parser.Results.label;
% 
% 
%     if settings.show_individual_data && ~settings.show_average_data && ~settings.show_spread_data
%         % show individual data only
%         singleScatterPlot ...
%           ( ...
%             ydata, ...
%             'axes', axes_handle, ...
%             'abscissa', xdata, ...
%             'color', lightenColor(color, 0.5), ...
%             'xlabel', label ...
%           );
%     end
%     if ~settings.show_individual_data && settings.show_average_data && ~settings.show_spread_data
%         % show average data only
%         plot ...
%           ( ...
%             axes_handle, ...
%             xdata, ...
%             mean(ydata), ...
%             'o', ...
%             'markerSize', 18, ...
%             'markerFaceColor', color, ...
%             'MarkerEdgeColor', 'none' ...
%           );
%     end
%     if ~settings.show_individual_data && ~settings.show_average_data && settings.show_spread_data
%         % show box only
%         singleBoxPlot ...
%           ( ...
%             ydata, ...
%             'axes', axes_handle, ...
%             'abscissa', xdata, ...
%             'FaceColor', color, ...
%             'xlabel', label, ...
%             'MedianColor', settings.edge_color, ...
%             'PlotMean', false, ...
%             'WiskColor', settings.edge_color, ...
%             'MarkerColor', lightenColor(color, 0.5), ...
%             'ShowOutliers', settings.show_outliers ...
%           )
%     end
%     if settings.show_individual_data && settings.show_average_data && ~settings.show_spread_data
%         % show individual data and mean
%         singleScatterPlot ...
%           ( ...
%             ydata, ...
%             'axes', axes_handle, ...
%             'abscissa', xdata, ...
%             'color', lightenColor(color, 0.5), ...
%             'xlabel', label ...
%           );
%         plot ...
%           ( ...
%             axes_handle, ...
%             xdata, ...
%             mean(ydata), ...
%             'o', ...
%             'markerSize', 18, ...
%             'markerFaceColor', color, ...
%             'MarkerEdgeColor', 'none' ...
%           );
%     end
%     if settings.show_individual_data && ~settings.show_average_data && settings.show_spread_data
%         % show box only
%         singleBoxPlot ...
%           ( ...
%             ydata, ...
%             'axes', axes_handle, ...
%             'abscissa', xdata, ...
%             'FaceColor', color, ...
%             'xlabel', label, ...
%             'MedianColor', settings.edge_color, ...
%             'PlotMean', false, ...
%             'ShowData', true, ...
%             'WiskColor', settings.edge_color, ...
%             'MarkerColor', lightenColor(color, 0.5), ...
%             'ShowOutliers', settings.show_outliers ...
%           )
%     end
%     if ~settings.show_individual_data && settings.show_average_data && settings.show_spread_data
%         % show box with mean
%         singleBoxPlot ...
%           ( ...
%             ydata, ...
%             'axes', axes_handle, ...
%             'abscissa', xdata, ...
%             'FaceColor', color, ...
%             'xlabel', label, ...
%             'MedianColor', settings.edge_color, ...
%             'PlotMean', true, ...
%             'MeanLinewidth', 3, ...
%             'MeanColor', lightenColor(color, 0.2), ...
%             'ShowData', false, ...
%             'WiskColor', settings.edge_color, ...
%             'MarkerColor', lightenColor(color, 0.5), ...
%             'ShowOutliers', settings.show_outliers ...
%           )
%     end
%     if settings.show_individual_data && settings.show_average_data && settings.show_spread_data
%         % show box with mean
%         singleBoxPlot ...
%           ( ...
%             ydata, ...
%             'axes', axes_handle, ...
%             'abscissa', xdata, ...
%             'FaceColor', color, ...
%             'xlabel', label, ...
%             'MedianColor', settings.edge_color, ...
%             'PlotMean', true, ...
%             'MeanLinewidth', 3, ...
%             'MeanColor', lightenColor(color, 0.2), ...
%             'ShowData', true, ...
%             'WiskColor', settings.edge_color, ...
%             'MarkerColor', lightenColor(color, 0.5), ...
%             'ShowOutliers', settings.show_outliers ...
%           )
%     end
% end

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



