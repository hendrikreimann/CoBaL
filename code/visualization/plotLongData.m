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


function plotLongData(varargin)
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
    mark_bands = plot_settings.get('mark_bands');
    
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');

    %% load data
    % declare variables
    conditions_settings = study_settings.get('conditions');
    condition_to_compare = plot_settings.get('condition_to_compare');
    condition_labels = conditions_settings(:, 1)';
    condition_source_variables = conditions_settings(:, 2)';
    number_of_condition_labels = length(condition_labels);
    condition_data_all = {};
    condition_control_data_all = {};
    variables_to_plot = plot_settings.get('long_variables_to_plot');
    plot_mode = plot_settings.get('plot_mode');
    
    % load data
    number_of_variables_to_plot = size(variables_to_plot, 1);
    step_time_data = [];
    data_all = cell(number_of_variables_to_plot, 1);
    data_control_all = cell(number_of_variables_to_plot, 1);
    data_response_all = cell(number_of_variables_to_plot, 1);
    data_folder_list = determineDataStructure(subjects);

    for i_folder = 1 : length(data_folder_list)
        % load data
        data_path = data_folder_list{i_folder};
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
        results_file_name = [data_path filesep 'analysis' filesep makeFileName(date, subject_id, 'longStretchResults')];
        loaded_data = load(results_file_name);
        number_of_stretches_this_session = size(loaded_data.long_stretch_times_session, 2);
        number_of_stretches_control_this_session = size(loaded_data.long_stretch_control_times_session, 2);
        bands_per_stretch = loaded_data.bands_per_long_stretch;

        % transform conditions into cell array
        conditions_session = loaded_data.long_stretch_conditions_session;
        condition_array_session = cell(number_of_stretches_this_session, number_of_condition_labels);
        conditions_control_session = loaded_data.long_stretch_control_conditions_session;
        condition_control_array_session = cell(number_of_stretches_control_this_session, number_of_condition_labels);
        for i_condition = 1 : number_of_condition_labels
            condition_array_session(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
            condition_control_array_session(:, i_condition) = conditions_control_session.(condition_source_variables{i_condition});
        end
        condition_data_all = [condition_data_all; condition_array_session]; %#ok<AGROW>
        condition_control_data_all = [condition_control_data_all; condition_control_array_session]; %#ok<AGROW>
        
        % extract data
        this_step_time_data = loaded_data.step_time_data_session;
        step_time_data = [step_time_data this_step_time_data]; %#ok<AGROW>
        for i_variable = 1 : number_of_variables_to_plot
            this_variable_name = variables_to_plot{i_variable, 1};
            this_variable_source = variables_to_plot{i_variable, 2};
            
            if strcmp(this_variable_source, 'stretch')
                index_in_saved_data = find(strcmp(loaded_data.long_stretch_data_labels_session, this_variable_name), 1, 'first');
                this_variable_data = loaded_data.long_stretch_data_session{index_in_saved_data};
                data_all{i_variable} = [data_all{i_variable} this_variable_data];
            end
            if strcmp(this_variable_source, 'response')
                index_in_saved_data = find(strcmp(loaded_data.long_stretch_data_labels_session, this_variable_name), 1, 'first');
                this_variable_response_data = loaded_data.long_stretch_response_data_session{index_in_saved_data};
                data_all{i_variable} = [data_all{i_variable} this_variable_response_data];
            end
            if strcmp(this_variable_source, 'analysis')
                index_in_saved_data = find(strcmp(loaded_data.analysis_names_session, this_variable_name), 1, 'first');
                this_variable_data = loaded_data.analysis_data_session{index_in_saved_data};
                data_all{i_variable} = [data_all{i_variable} this_variable_data];
            end
            
%             data_all{i_variable} = [data_all{i_variable} this_variable_data];
%             data_response_all{i_variable} = [data_response_all{i_variable} this_variable_response_data];
        end
    end

    
    %% populate condition cell
    labels_to_ignore = plot_settings.get('conditions_to_ignore');
    levels_to_remove = plot_settings.get('levels_to_remove');
    preferred_level_order = plot_settings.get('preferred_level_order');
    [condition_combination_labels, condition_combinations_stimulus, condition_combinations_control] = determineConditionCombinations(condition_data_all, conditions_settings, labels_to_ignore, levels_to_remove);
%     condition_combinations_stimulus = sortConditionCombinations(condition_combinations_stimulus, condition_combination_labels, condition_to_compare, preferred_level_order);
    [condition_combinations_stimulus, condition_combinations_control] = sortConditionCombinations(condition_combinations_stimulus, condition_combinations_control, condition_combination_labels, condition_to_compare, preferred_level_order);
    [comparison_indices, conditions_per_comparison_max] = determineComparisons(condition_combinations_stimulus, condition_combination_labels, plot_settings);
    number_of_comparisons = length(comparison_indices);
    
    %% create figures and determine abscissae for each comparison
    comparison_variable_to_axes_index_map = zeros(number_of_comparisons, 1);
    abscissae_cell = cell(number_of_comparisons, number_of_variables_to_plot);
    comparison_path_to_axes_index_map = zeros(number_of_comparisons, 1);

    
    
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
    for i_variable = 1 : number_of_variables_to_plot
        for i_comparison = 1 : number_of_comparisons
            this_comparison = comparison_indices{i_comparison};
            % make figure and axes
            new_figure = figure; new_axes = axes; hold on;

            % store handles and determine abscissa data
            trajectory_figure_handles(i_comparison, i_variable) = new_figure;
            trajectory_axes_handles(i_comparison, i_variable) = new_axes;
            comparison_variable_to_axes_index_map(i_comparison) = i_comparison;

            % scale abscissae
            conditions_this_comparison = comparison_indices{i_comparison};
            step_time_means_this_comparison = zeros(bands_per_stretch, size(conditions_this_comparison, 2));
            for i_condition = 1 : length(conditions_this_comparison)
                this_condition_combination = condition_combinations_stimulus(conditions_this_comparison(i_condition), :);
                this_condition_indicator = getConditionIndicator(this_condition_combination, condition_combination_labels, condition_data_all, condition_labels);
                step_time_data_this_condition = step_time_data(:, this_condition_indicator);
                step_time_means_this_comparison(:, i_condition) = nanmean(step_time_data_this_condition, 2);
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

            % set axes properties
            if dictate_axes && ~(strcmp(plot_mode, 'detailed') && isDiscreteVariable(i_variable, data_all, bands_per_stretch))
%                 set(gca, 'xlim', [time_normalized(1), time_normalized(end)]);
                set(gca, 'ylim', [str2double(variables_to_plot{i_variable, 6}), str2double(variables_to_plot{i_variable, 7})]);
            end

            % set axis labels
            if strcmp(plot_settings.get('time_plot_style'), 'scaled_to_comparison_mean') || strcmp(plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
                xlabel('normalized time (s)');
            else
                xlabel('normalized time (%)');
            end
            ylabel(variables_to_plot{i_variable, 4});

            % add text labels
            pos_text_handles(i_comparison, i_variable) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    'TBD', ...
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
                    'TBD', ...
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
        for i_comparison = 1 : number_of_comparisons
            target_abscissae = abscissae_cell{i_comparison, i_variable};
            xlim = [nanmin(target_abscissae(:, 1)) nanmax(target_abscissae(:, end))];
            % set x-limits accordingly
            set(trajectory_axes_handles(i_comparison, i_variable), 'xlim', xlim);
        end
    end

    % determine stance start and end times and stance foot
    for i_variable = 1 : number_of_variables_to_plot
        for i_comparison = 1 : number_of_comparisons

            step_abscissa = abscissae_cell{i_comparison, i_variable};
            step_start_times_cell{i_comparison, i_variable} = step_abscissa(1, 1);
            step_end_times_cell{i_comparison, i_variable} = step_abscissa(1, end);

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
            
            % plot stimulus
            for i_condition = 1 : length(conditions_this_comparison)
                this_condition_index = conditions_this_comparison(i_condition);
                this_condition = condition_combinations_stimulus(this_condition_index, :);
                label_string = strrep(this_condition{strcmp(condition_combination_labels, condition_to_compare)}, '_', ' ');
                this_condition_indicator = getConditionIndicator(this_condition, condition_combination_labels, condition_data_all, condition_labels);
                data_to_plot_this_condition = data_to_plot(:, this_condition_indicator);
                
                target_abscissa = abscissae_cell{i_comparison, i_variable}(i_condition, :);
                if strcmp(plot_mode, 'detailed')
                    for i_stretch = 1 : size(data_to_plot_this_condition, 2)
                        plot ...
                          ( ...
                            target_axes_handle, ...
                            target_abscissa, ...
                            data_to_plot_this_condition(:, i_stretch), ...
                            'HandleVisibility', 'off', ...
                            'color', lightenColor(colors_comparison(i_condition, :), 0.5) ...
                          );
                    end
                    condition_mean_plot = plot ...
                      ( ...
                        target_axes_handle, ...
                        target_abscissa, ...
                        mean(data_to_plot_this_condition, 2), ...
                        'linewidth', 5, ...
                        'color', colors_comparison(i_condition, :) ...
                      );                    
                end
                if strcmp(plot_mode, 'overview') || strcmp(plot_mode, 'episodes')
                    plot_handles = shadedErrorBar ...
                      ( ...
                        target_abscissa, ...
                        nanmean(data_to_plot_this_condition, 2), ...
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
    
    %% shade bands
    if mark_bands
        for i_comparison = 1 : number_of_comparisons
            for i_variable = 1 : number_of_variables_to_plot
                these_axes = trajectory_axes_handles(i_comparison, i_variable);
                these_abscissae = abscissae_cell{i_comparison, i_variable};
                ylimits = get(these_axes, 'ylim');

                if mark_bands == 1
                    bands_to_mark = 2 : 2 : bands_per_stretch;
                end
                if mark_bands == 2
                    bands_to_mark = 1 : 2 : bands_per_stretch;
                end

                for i_band = bands_to_mark
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















