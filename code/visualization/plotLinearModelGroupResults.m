%     This file is part of the CoBaL code base
%     Copyright (C) 2021 Hendrik Reimann <hendrikreimann@gmail.com>
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

function plotLinearModelGroupResults(varargin)
    % parse input parameters
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'show_legend', true)
    addParameter(parser, 'save', false)
    parse(parser, varargin{:})
    arguments.show_legend = parser.Results.show_legend;
    arguments.save_results = parser.Results.save;

    % load settings
    study_settings = loadSettingsFromFile('study');
    linear_model_settings = loadSettingsFromFile('linearModel');
    condition_to_compare = linear_model_settings.get('condition_to_compare');
    arguments.separate_figures = linear_model_settings.get('separate_figures');
    
    % load data
    model_data_filename = ['groupResults' filesep 'linearModelResults.mat'];
    model_data = load(model_data_filename);
    model_list = linear_model_settings.getTable('plot_table');
    
    % process
    preferred_level_order = linear_model_settings.get('preferred_level_order', 1);
    for i_model = 1 : size(model_list, 1)
        this_model_label = model_list.label{i_model};
        this_model_index = findModelIndex(model_data.model_results, this_model_label);
        this_model_data = model_data.model_results(this_model_index);
        
        % determine comparisons
        conditions_settings = study_settings.get('conditions');
        condition_labels = conditions_settings(:, 1)';
        condition_data = extractConditionData(this_model_data, condition_labels);
        labels_to_ignore = linear_model_settings.get('conditions_to_ignore', 1);
        levels_to_remove = linear_model_settings.get('levels_to_remove', 1);
        [comparisons.condition_combination_labels, comparisons.condition_combinations] = determineConditionCombinations(condition_data, conditions_settings, labels_to_ignore, levels_to_remove);
        [comparisons.condition_combinations, comparisons.level_order] = sortConditionCombinations(comparisons.condition_combinations, comparisons.condition_combination_labels, condition_to_compare, preferred_level_order);
        [comparisons.comparison_indices, comparisons.conditions_per_comparison_max] = determineComparisons(comparisons.condition_combinations, comparisons.condition_combination_labels, condition_to_compare, conditions_settings);
        comparisons.number_of_comparisons = length(comparisons.comparison_indices);
        comparisons.condition_colors = determineConditionColors(linear_model_settings.getTable('colors', 1), comparisons, condition_to_compare);
        comparisons.condition_to_compare = condition_to_compare;

        % loop through comparisons
        for i_comparison = 1 : comparisons.number_of_comparisons
            createComparisonFigure(this_model_data, comparisons, i_comparison, linear_model_settings, study_settings, arguments);
        end
    end
end

function index = findModelIndex(model_data, requested_label)
    % loop through models to find the requested one
    index = 0;
    for i_model = 1 : size(model_data, 1)
        if strcmp(model_data(i_model).label, requested_label)
            index = i_model;
        end
    end
    
    if index == 0
        error(['Model "' requested_label '" not available'])
    end
end

function comparison_label = createComparisonLabel(row_info, comparison_indices, column_with_condition_to_compare)
    this_row_index = comparison_indices(1);

    % identify columns for which all entries are the same
    number_of_columns = size(row_info, 2);
    number_of_unique_entries = zeros(1, number_of_columns);
    for i_column = 1 : number_of_columns
        this_column = row_info(:, i_column);
        number_of_unique_entries(i_column) = numel(unique(this_column));
    end
    columns_with_relevant_information = logical(number_of_unique_entries - 1);
    
    % remove the column for the condition we're actually comparing
    columns_with_relevant_information(column_with_condition_to_compare) = false;
    
    % pull out relevant row info
    relevant_row_info = row_info(this_row_index, columns_with_relevant_information);
    
    % make label
    if isempty(relevant_row_info)
        comparison_label = '';
    else
        comparison_label = strrep(relevant_row_info{1}, '_', ' ');
        for i_column = 2 : length(relevant_row_info)
            comparison_label = [comparison_label ', ' strrep(relevant_row_info{i_column}, '_', ' ')]; %#ok<AGROW>
        end
    end
end

function condition_data = extractConditionData(model_data, condition_labels)
    condition_data = {};
    for i_condition = 1 : length(condition_labels)
        this_condition_label = condition_labels{i_condition};
        
        if strcmp(model_data.type, 'discrete')
            this_condition_data = model_data.data.(this_condition_label);
        elseif strcmp(model_data.type, 'continuous')
            column_index = strcmp(model_data.headers, this_condition_label);
            this_condition_data = model_data.data(:, column_index);
        end
        condition_data = [condition_data this_condition_data]; %#ok<AGROW>
    end

end

function createComparisonFigure(model_data, comparisons, comparison_to_show, linear_model_settings, study_settings, arguments)
    % get model info from settings
    model_list = linear_model_settings.getTable('models');
    model_list = [model_list; linear_model_settings.getTable('models_second_order')];
    this_model_row = strcmp(model_list.label, model_data.label);
    model_label = model_list.label{this_model_row};
    outcome_label = model_list.outcome_variable_name{this_model_row};
    predictors_label = model_list.predictor_variable_list{this_model_row};
    
    % TODO: find a way to check if this predictor is present
    if linear_model_settings.settingIsPresent(predictors_label)
        % predictors_label is the name of a table in the settings file
        predictors_list = linear_model_settings.get(predictors_label, 1);
    else
        % predictors_label is the name of a variable
        predictors_list = {predictors_label};
    end
    
    number_of_predictors = length(predictors_list);
    number_of_time_points_per_stretch = study_settings.get('number_of_time_steps_normalized');

    % prepare comparison
    comparison_indices = comparisons.comparison_indices{comparison_to_show};
    condition_to_compare = linear_model_settings.get('condition_to_compare');
    if strcmp(model_data.type, 'discrete')
        header = model_data.data.Properties.VariableNames;
        data = table2cell(model_data.data);
    elseif strcmp(model_data.type, 'continuous')
        header = model_data.headers;
        data = model_data.data;
    end
    relevant_column = strcmp(comparisons.condition_combination_labels, condition_to_compare);
    comparison_label = createComparisonLabel(comparisons.condition_combinations, comparison_indices, relevant_column);
    number_of_conditions_in_this_comparison = length(comparison_indices);
    
    % prepare figures
    slope_axes = cell(number_of_predictors, 1);
    variance_axes = cell(3, 1);
    
    % create figure and axes
    if arguments.separate_figures
        title_label = comparison_label;
        slope_figures = cell(number_of_predictors, 1);
        slope_figure_file_names = cell(number_of_predictors, 1);
        variance_figures = cell(3, 1);
        variance_figure_file_names = cell(3, 1);
        x_label = comparisons.condition_to_compare;
        file_labels = cell(size(predictors_list));
        if isempty(comparison_label)
            this_file_suffix = [];
        else
            this_file_suffix = ['-' strrep(comparison_label, ', ', '_')];
        end
    
        % create slope figures
        for i_predictor = 1 : number_of_predictors
            slope_figures{i_predictor} = figure;
            slope_figure_file_names{i_predictor} = [model_label '_slope_' num2str(i_predictor) this_file_suffix];
            slope_axes{i_predictor} = axes;
            set(slope_axes{i_predictor}, 'fontsize', 12);
            set(slope_axes{i_predictor}, 'XTickLabel', {}, 'xtick', [])
            xlabel(x_label);
            this_predictor_label = predictors_list{i_predictor};
            this_slope_label = {['\partial[' strrep(outcome_label, '_', ' ') ']'], ['\partial' '[' strrep(this_predictor_label, '_', ' ') ']']};
            ylabel(this_slope_label);
            hold on;
            title(title_label)
            this_file_label = [outcome_label '_VS_' predictors_list{i_predictor} '_' comparison_label];
            file_labels{i_predictor} = this_file_label;
            set(slope_axes{i_predictor}, 'XTickLabel', {}, 'xtick', []);
        end

        % create R^2 figure
        variance_figures{1} = figure;
        variance_figure_file_names{1} = [model_label '_rsquare' this_file_suffix];
        variance_axes{1} = axes;
        set(variance_axes{1}, 'fontsize', 12);
        hold on;
        ylabel([strrep(model_label, '_', ' ') ': R^2']);
        ylim([0 1]);
        xlabel(variance_axes{1}, x_label);
        title(variance_axes{1}, strrep(title_label, '_', ' '));
        set(variance_axes{1}, 'XTickLabel', {}, 'xtick', []);

        % create ss_difference figure
        variance_figures{2} = figure;
        variance_figure_file_names{2} = [model_label '_SSdiff' this_file_suffix];
        variance_axes{2} = axes;
        set(variance_axes{2}, 'fontsize', 12);
        hold on;
        ylabel([strrep(model_label, '_', ' ') ': SS_{total} - SS_{residual}']);
        xlabel(variance_axes{2}, x_label);
        title(variance_axes{2}, strrep(title_label, '_', ' '));
        set(variance_axes{2}, 'XTickLabel', {}, 'xtick', []);

        % create ss_residual figure
        variance_figures{3} = figure;
        variance_figure_file_names{3} = [model_label '_SSres' this_file_suffix];
        variance_axes{3} = axes;
        set(variance_axes{3}, 'fontsize', 12);
        hold on;
        ylabel([strrep(model_label, '_', ' ') ': SS_{residual}']);
        xlabel(variance_axes{3}, x_label);
        title(variance_axes{3}, strrep(title_label, '_', ' '));
        set(variance_axes{3}, 'XTickLabel', {}, 'xtick', []);

    else
        % combine plots into two figures, one for slopes and one for residual measures
        title_label = ['Outcome: ' outcome_label ' - ' comparison_label];
        file_label = [outcome_label '_VS_' predictors_label '_' comparison_label];
        
        % slope
        figure;
        tiledlayout(number_of_predictors, 1);

        for i_predictor = 1 : number_of_predictors
            slope_axes{i_predictor} = nexttile;
            set(slope_axes{i_predictor}, 'fontsize', 12);
            set(slope_axes{i_predictor}, 'xtick', [])
            ylabel('slope');
            hold on;
            title(strrep(predictors_list{i_predictor}, '_', ' '), 'Units', 'normalized', 'Position', [0.5, 0.9, 0])
            set(slope_axes{i_predictor}, 'XTickLabel', {}, 'xtick', []);
        end
        uicontrol('style', 'text', 'string', title_label, 'units', 'normalized', 'position', [0, 0.95, 1, 0.05], 'fontsize', 16, 'FontWeight', 'bold');

        % residuals
        figure;
        tiledlayout(3, 1);
        variance_axes{1} = nexttile; % was: r_square_axes
        set(variance_axes{1}, 'fontsize', 12);
        hold on;
        ylabel('R^2');
        set(variance_axes{1}, 'XTickLabel', {}, 'xtick', []);
        ylim([0 1]);
        
        variance_axes{2} = nexttile; % was: ss_difference_axes
        set(variance_axes{2}, 'fontsize', 12);
        hold on;
        set(variance_axes{2}, 'XTickLabel', {}, 'xtick', []);
        ylabel('SS_{tot} - SS_{res}');

        variance_axes{3} = nexttile; % was: ss_residual_axes
        set(variance_axes{3}, 'fontsize', 12);
        hold on;
        set(variance_axes{3}, 'XTickLabel', {}, 'xtick', []);
        ylabel('SS_{res}');

        % define parameters depending on model type, discrete vs. continuous predictor
        if strcmp(model_data.type, 'discrete')
            linestyle = 'none';
            linewidth = 1;
            marker = 's';
            marker_size = 15;
            x_label = comparisons.condition_to_compare;
            set(variance_axes{1}, 'XTickLabel', {}, 'xtick', []);
        elseif strcmp(model_data.type, 'continuous')
            linestyle = '-';
            linewidth = 2;
            marker = 'none';
            marker_size = 1;
            x_label = 'time (%)';
        end

        % add labels
        xlabel(variance_axes{1}, x_label); 
        uicontrol('style', 'text', 'string', title_label, 'units', 'normalized', 'position', [0, 0.95, 1, 0.05], 'fontsize', 16, 'FontWeight', 'bold');
    end
    
    % loop through conditions
    for i_condition = 1 : number_of_conditions_in_this_comparison
        this_condition_index = comparison_indices(i_condition);
        
        % this_condition_index is in relation to comparisons.condition_combinations, I need to find the relevant index in model_data
        this_condition = comparisons.condition_combinations(this_condition_index, :);
        levels_to_remove = [];
        this_condition_indicator = getConditionIndicator ...
          ( ...
            this_condition, ...
            comparisons.condition_combination_labels, ...
            data, ...
            header, ...
            levels_to_remove ...
          );
        
        
        this_condition_label = strrep(this_condition{relevant_column}, '_', ' ');
        
        if strcmp(model_data.type, 'discrete')
            abscissa_data = i_condition;
            for i_predictor = 1 : number_of_predictors
                slope_axes{i_predictor}.XTick = [slope_axes{i_predictor}.XTick, i_condition];
                slope_axes{i_predictor}.XTickLabel = [slope_axes{i_predictor}.XTickLabel; this_condition_label];
            end
            for i_variance = 1 : 3
                variance_axes{i_variance}.XTick = [variance_axes{i_variance}.XTick, i_condition];
                variance_axes{i_variance}.XTickLabel = [variance_axes{i_variance}.XTickLabel; this_condition_label];
            end
        elseif strcmp(model_data.type, 'continuous')
            abscissa_data = (1 : number_of_time_points_per_stretch) - 1;
        end
        r_square_column = strcmp(header, 'R_square');
        r_square_data_this_condition_cell = data(this_condition_indicator, r_square_column);
        r_square_data_here = cell2mat(r_square_data_this_condition_cell');
        
        ss_difference_column = strcmp(header, 'SS_difference');
        ss_difference_data_this_condition_cell = data(this_condition_indicator, ss_difference_column);
        ss_difference_data_here = cell2mat(ss_difference_data_this_condition_cell');
        
        ss_residual_column = strcmp(header, 'SS_residual');
        ss_residual_data_this_condition_cell = data(this_condition_indicator, ss_residual_column);
        ss_residual_data_here = cell2mat(ss_residual_data_this_condition_cell');
        
        % get slope data
        slope_data_this_condition = cell(size(predictors_list));
        for i_predictor = 1 : length(predictors_list)
            % extract data
%             this_predictor_header = ['d_' outcome_label '_by_d_' predictors_list{i_predictor}];
            this_predictor_header = [model_label '_slope_' num2str(i_predictor)];
            this_predictor_slope_column = strcmp(header, this_predictor_header);
            this_predictor_this_condition_slope_data_cell = data(this_condition_indicator, this_predictor_slope_column);
            this_predictor_this_condition_slope_data = cell2mat(this_predictor_this_condition_slope_data_cell');
            slope_data_this_condition{i_predictor} = this_predictor_this_condition_slope_data;
            
            % get subject data, for debugging
            subjects_here = data(this_condition_indicator, strcmp(header, 'subject'));
        end        
        
        this_condition = comparisons.condition_combinations(this_condition_index, :);
        this_label = this_condition{strcmp(comparisons.condition_combination_labels, comparisons.condition_to_compare)};
        this_color = comparisons.condition_colors{strcmp(comparisons.condition_colors(:, 1), this_label), 2};
        
        if strcmp(model_data.type, 'discrete')
            % plot slope data
            for i_predictor = 1 : length(predictors_list)
                slope_data_here = slope_data_this_condition{i_predictor};
                axes_here = slope_axes{i_predictor};
                singleBoxPlot ...
                 ( ...
                   slope_data_here, ...
                   'axes', axes_here, ...
                   'abscissa', abscissa_data, ...
                   'width', 0.4, ...
                   'FaceColor', lightenColor(this_color, 0.8), ...
                   'EdgeColor', this_color, ...
                   'ShowData', true, ...
                   'MarkerColor', lightenColor(this_color, 0.5), ...
                   'MedianColor', this_color, ...
                   'MeanColor', this_color, ...
                   'WiskColor', this_color, ...
                   'ShowOutliers', false, ...
                   'PlotMean', true, ...
                   'PlotMedian', true, ...
                   'DataMarkerSize', 18 ...
                 )             
            end
            
            % plot R^2 data
            singleBoxPlot ...
              ( ...
                r_square_data_here, ...
                'axes', variance_axes{1}, ...
                'abscissa', abscissa_data, ...
                'width', 0.4, ...
                'FaceColor', lightenColor(this_color, 0.8), ...
                'EdgeColor', this_color, ...
                'ShowData', true, ...
                'MarkerColor', lightenColor(this_color, 0.5), ...
                'MedianColor', this_color, ...
                'MeanColor', this_color, ...
                'WiskColor', this_color, ...
                'ShowOutliers', false, ...
                'PlotMean', true, ...
                'PlotMedian', true, ...
                'DataMarkerSize', 18 ...
              )
            % plot ss_difference data
            singleBoxPlot ...
              ( ...
                ss_difference_data_here, ...
                'axes', variance_axes{2}, ...
                'abscissa', abscissa_data, ...
                'width', 0.4, ...
                'FaceColor', lightenColor(this_color, 0.8), ...
                'EdgeColor', this_color, ...
                'ShowData', true, ...
                'MarkerColor', lightenColor(this_color, 0.5), ...
                'MedianColor', this_color, ...
                'MeanColor', this_color, ...
                'WiskColor', this_color, ...
                'ShowOutliers', false, ...
                'PlotMean', true, ...
                'PlotMedian', true, ...
                'DataMarkerSize', 18 ...
              )
            % plot ss_residual data
            singleBoxPlot ...
              ( ...
                ss_residual_data_here, ...
                'axes', variance_axes{3}, ...
                'abscissa', abscissa_data, ...
                'width', 0.4, ...
                'FaceColor', lightenColor(this_color, 0.8), ...
                'EdgeColor', this_color, ...
                'ShowData', true, ...
                'MarkerColor', lightenColor(this_color, 0.5), ...
                'MedianColor', this_color, ...
                'MeanColor', this_color, ...
                'WiskColor', this_color, ...
                'ShowOutliers', false, ...
                'PlotMean', true, ...
                'PlotMedian', true, ...
                'DataMarkerSize', 18 ...
              )

            
        elseif strcmp(model_data.type, 'continuous')
            plot ...
              ( ...
                variance_axes{1}, ...
                abscissa_data, ...
                mean(r_square_data_here, 2), ...
                'color', this_color, ...
                'linewidth', 3, ...
                'DisplayName', [this_condition_label ' - data'] ...
              );
            % plot slope data
            for i_predictor = 1 : length(predictors_list)
                slope_data_here = slope_data_this_condition{i_predictor};
                axes_here = slope_axes{i_predictor};
                plot ...
                  ( ...
                    axes_here, ...
                    abscissa_data, ...
                    mean(slope_data_here, 2), ...
                    'color', this_color, ...
                    'linewidth', 3, ...
                    'DisplayName', [this_condition_label ' - data'] ...
                  );
                
                
                
            end            
        end
    end

    % adjust axis-limits
%     if model_data.predictor_variable_data_points_per_stretch == 1
        xlimits = [variance_axes{1}.XTick(1)-0.5, variance_axes{1}.XTick(end)+0.5];
        set(variance_axes{1}, 'xlim', xlimits);
        for i_predictor = 1 : number_of_predictors
            set(slope_axes{i_predictor}, 'xlim', xlimits);
        end
%     end
    ylimits_SS_difference = get(variance_axes{2}, 'ylim');
    ylimits_SS_residual = get(variance_axes{3}, 'ylim');
    ylimits_new = [min([ylimits_SS_difference(1) ylimits_SS_residual(1)]) max([ylimits_SS_difference(2) ylimits_SS_residual(2)])];
    set(variance_axes{2}, 'ylim', ylimits_new);
    set(variance_axes{3}, 'ylim', ylimits_new);
    
    % plot zero line
    if linear_model_settings.get('plot_zero', 1)
        for i_predictor = 1 : number_of_predictors
            xlimits = get(slope_axes{i_predictor}, 'xlim');
            zero_plot = plot(slope_axes{i_predictor}, xlimits, [0 0], 'color', [0.7 0.7 0.7], 'linewidth', 2);
            set(zero_plot, 'HandleVisibility', 'off');
            uistack(zero_plot, 'bottom')
        end
    end

    % save
    if arguments.save_results
        if ~directoryExists('figures')
            mkdir figures
        end
        if ~directoryExists(['figures' filesep 'noLabels'])
            mkdir(['figures' filesep 'noLabels'])
        end
        if ~directoryExists(['figures' filesep 'withLabels'])
            mkdir(['figures' filesep 'withLabels'])
        end
        
        
        if arguments.separate_figures
            for i_predictor = 1 : number_of_predictors
                this_figure = slope_figures{i_predictor};
                these_axes = slope_axes{i_predictor};
                this_filename_with = ['figures' filesep 'withLabels' filesep slope_figure_file_names{i_predictor} '.jpg'];
                print(this_figure, this_filename_with, '-r300', '-djpeg')

                set(get(these_axes, 'xaxis'), 'visible', 'off');
                set(get(these_axes, 'yaxis'), 'visible', 'off');
                set(get(these_axes, 'xlabel'), 'visible', 'off');
                set(get(these_axes, 'ylabel'), 'visible', 'off');
                set(get(these_axes, 'title'), 'visible', 'off');
                set(these_axes, 'xticklabel', '');
                set(these_axes, 'yticklabel', '');
                set(these_axes, 'position', [0 0 1 1]);

                filename_without = ['figures' filesep 'noLabels' filesep slope_figure_file_names{i_predictor} '.jpg'];
                print(this_figure, filename_without, '-r300', '-djpeg')
                close(this_figure);
            end

            for i_figure = 1 : 3
                this_figure = variance_figures{i_figure};
                these_axes = variance_axes{i_figure};
                this_filename_with = ['figures' filesep 'withLabels' filesep variance_figure_file_names{i_figure} '.jpg'];
                print(this_figure, this_filename_with, '-r300', '-djpeg')

                set(get(these_axes, 'xaxis'), 'visible', 'off');
                set(get(these_axes, 'yaxis'), 'visible', 'off');
                set(get(these_axes, 'xlabel'), 'visible', 'off');
                set(get(these_axes, 'ylabel'), 'visible', 'off');
                set(get(these_axes, 'title'), 'visible', 'off');
                set(these_axes, 'xticklabel', '');
                set(these_axes, 'yticklabel', '');
                set(these_axes, 'position', [0 0 1 1]);

                filename_without = ['figures' filesep 'noLabels' filesep variance_figure_file_names{i_figure} '.jpg'];
                print(this_figure, filename_without, '-r300', '-djpeg')
                close(this_figure);
            end
        else
            filename = ['figures' filesep file_label '.jpg'];
            print(gcf, filename, '-r300', '-djpeg')
        end
    end    
    
end



