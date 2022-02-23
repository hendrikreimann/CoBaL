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

function plotLinearModelResults_single(varargin)
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
    subject_settings = loadSettingsFromFile('subject');
    linear_model_settings = loadSettingsFromFile('linearModel');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    condition_to_compare = linear_model_settings.get('condition_to_compare');

    % load data
    model_file_name = ['results' filesep collection_date '_' subject_id '_linearModelsSingle.mat'];
    model_data = load(model_file_name);
    model_list = linear_model_settings.getTable('plot_table_single');
    
    preferred_level_order = linear_model_settings.get('preferred_level_order', 1);
    for i_model = 1 : size(model_list, 1)
        this_model = model_list(i_model, :);
        
        % find the requested model in the data
        this_model_index = findModelIndex(model_data, this_model);
        this_model_data = model_data.linear_model_results{this_model_index, strcmp(model_data.linear_model_results_header, 'results')};
        
        % create figures
        relevant_column = strcmp(this_model_data.row_info_headers, condition_to_compare);
        
        % determine comparisons
        condition_data = this_model_data.row_info;
        conditions_settings = study_settings.get('conditions');
        labels_to_ignore = {};
        levels_to_remove = {};
        [comparisons.condition_combination_labels, comparisons.condition_combinations] = determineConditionCombinations(condition_data, conditions_settings, labels_to_ignore, levels_to_remove);
        [comparisons.condition_combinations, comparisons.level_order] = sortConditionCombinations(comparisons.condition_combinations, comparisons.condition_combination_labels, condition_to_compare, preferred_level_order);
        [comparisons.comparison_indices, comparisons.conditions_per_comparison_max] = determineComparisons(comparisons.condition_combinations, comparisons.condition_combination_labels, condition_to_compare, conditions_settings);
        comparisons.number_of_comparisons = length(comparisons.comparison_indices);
        comparisons.condition_colors = determineConditionColors(linear_model_settings.getTable('colors', 1), comparisons, condition_to_compare);
        comparisons.condition_to_compare = condition_to_compare;
        
        % loop through comparisons
        for i_comparison = 1 : comparisons.number_of_comparisons
            if this_model_data.predictor_variable_data_points_per_stretch == 1
                % discrete predictor
                createComparisonFigure_discrete(this_model_data, comparisons, i_comparison, relevant_column, linear_model_settings, arguments);
            else
                % continuous predictor
                createComparisonFigure_continuous(this_model_data, comparisons, i_comparison, relevant_column, linear_model_settings, arguments);
            end
        end
    end

end

function index = findModelIndex(data, model)
    % initialize
    index = [];
    number_of_columns = size(model, 2);
    number_of_data_rows = size(data.linear_model_results, 1);
    
    % loop through columns
    column_names_requested = model.Properties.VariableNames;
    column_names_in_data = data.linear_model_results_header;
    column_matches = false(number_of_data_rows , number_of_columns);
    for i_column = 1 : number_of_columns
        % get column name
        this_column_name = column_names_requested{i_column};
        this_column_value = model{1, i_column}{1};
        
        % find matching column in data
        this_column_index_in_data = strcmp(this_column_name, column_names_in_data);
        
        if any(this_column_index_in_data)
            this_column_names_in_data = data.linear_model_results(:, this_column_index_in_data);
            this_column_match = strcmp(this_column_names_in_data, this_column_value);
            column_matches(:, i_column) = this_column_match;
        else
            return
        end
    end
    
    index = find(all(column_matches, 2), 1);

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
    comparison_label = strrep(relevant_row_info{1}, '_', ' ');
    for i_column = 2 : length(relevant_row_info)
        comparison_label = [comparison_label ', ' strrep(relevant_row_info{i_column}, '_', ' ')]; %#ok<AGROW>
    end
end

function createComparisonFigure_discrete(model_data, comparisons, comparison_to_show, relevant_column, linear_model_settings, arguments)
    comparison_indices = comparisons.comparison_indices{comparison_to_show};
    model_label = [model_data.outcome ' vs. ' model_data.predictor];
    comparison_label = createComparisonLabel(comparisons.condition_combinations, comparison_indices, relevant_column);
    number_of_conditions_in_this_comparison = length(comparison_indices);
    title_label = [model_label ', ' comparison_label];
    file_label = [model_data.outcome '_VS_' model_data.predictor '_' comparison_label];
    predictor_label = strrep(model_data.predictor, '_', ' ');
    outcome_label = strrep(model_data.outcome, '_', ' ');
    
    % create figure and axes
    colors = comparisons.condition_colors;
    figure;

    slope_axes = axes('position', [0.08 0.08, 0.4 0.38], 'fontsize', 12);
    hold on;
    xlabel('time (%)'); ylabel('slope');
    xlim([0.5 number_of_conditions_in_this_comparison+0.5]);
    uicontrol('style', 'text', 'string', title_label, 'units', 'normalized', 'position', [0, 0.95, 1, 0.05], 'fontsize', 16, 'FontWeight', 'bold');

    r_square_axes = axes('position', [0.08 0.5, 0.4 0.38], 'fontsize', 12); 
    hold on;
    ylabel('R^2');
    set(r_square_axes, 'xtick', [])
    xlim([0.5 number_of_conditions_in_this_comparison+0.5]);
    ylim([0 1]);
    set(r_square_axes, 'xtick', [])
    if arguments.show_legend
        legend('Location', 'best')
    end

    data_axes = axes('position', [0.58 0.08, 0.4 0.8], 'fontsize', 12); 
    hold on;
    xlabel(predictor_label);
    ylabel(outcome_label);

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
            model_data.row_info, ...
            model_data.row_info_headers, ...
            levels_to_remove ...
          );
        
        this_condition_label = strrep(model_data.row_info{this_condition_indicator, relevant_column}, '_', ' ');
        
        predictor_data_here = model_data.data.predictor{this_condition_indicator};
        outcome_data_here = model_data.data.outcome{this_condition_indicator};
        r_square_data_here = model_data.R_square{this_condition_indicator};
        slope_data_here = model_data.slope{this_condition_indicator};
        slope_cinv_here = model_data.slope_confidence_interval_width{this_condition_indicator};
        offset_data_here = model_data.offset{this_condition_indicator};

        this_condition = comparisons.condition_combinations(this_condition_index, :);
        this_label = this_condition{strcmp(comparisons.condition_combination_labels, comparisons.condition_to_compare)};
        this_color = comparisons.condition_colors{strcmp(comparisons.condition_colors(:, 1), this_label), 2};
        
        % plot data
        plot ...
          ( ...
            data_axes, ...
            predictor_data_here, ...
            outcome_data_here, ...
            'o', ...
            'DisplayName', [this_condition_label ' - data'], ...
            'color', 'none', ...
            'MarkerFaceColor', lightenColor(this_color, 0.5), ...
            'MarkerSize', 5 ...
          );

        % plot linear model
        predictor_interval = [min(predictor_data_here) max(predictor_data_here)];
        plot ...
          ( ...
            data_axes, ...
            predictor_interval, ...
            predictor_interval * slope_data_here + offset_data_here, ...
            '-', ...
            'DisplayName', [this_condition_label ' - model'], ...
            'color', this_color, ...
            'linewidth', 2 ...
          );

        % plot slope overview
        plot ...
          ( ...
            slope_axes, ...
            i_condition * [1 1], ...
            slope_data_here + [-1 1]*slope_cinv_here, ...
            '-', ...
            'color', this_color, ...
            'linewidth', 5 ...
          );
        plot ...
          ( ...
            slope_axes, ...
            i_condition, ...
            slope_data_here, ...
            's', ...
            'color', 'none', ...
            'MarkerFaceColor', this_color, ...
            'MarkerSize', 15 ...
          );

        % plot R^2 overview
        plot ...
          ( ...
            r_square_axes, ...
            i_condition, ...
            r_square_data_here, ...
            's', ...
            'color', 'none', ...
            'MarkerFaceColor', this_color, ...
            'DisplayName', this_condition_label, ...
            'MarkerSize', 15 ...
          );              
    end

    % plot zero line
    if linear_model_settings.get('plot_zero', 1)
        xlimits = get(slope_axes, 'xlim');
        zero_plot = plot(slope_axes, xlimits, [0 0], 'color', [0.7 0.7 0.7], 'linewidth', 2);
        set(zero_plot, 'HandleVisibility', 'off');
        uistack(zero_plot, 'bottom')
    end

    % save
    if arguments.save_results
        if ~directoryExists('figures')
            mkdir figures
        end
        filename = ['figures' filesep file_label '.jpg'];
        print(gcf, filename, '-r300', '-djpeg')
    end    
    
end

function createComparisonFigure_continuous(model_data, comparisons, comparison_to_show, relevant_column, linear_model_settings, arguments)
    comparison_indices = comparisons.comparison_indices{comparison_to_show};
    model_label = [model_data.outcome ' vs. ' model_data.predictor];
    comparison_label = createComparisonLabel(comparisons.condition_combinations, comparison_indices, relevant_column);
    number_of_conditions_in_this_comparison = length(comparison_indices);
    title_label = [model_label ', ' comparison_label];
    file_label = [model_data.outcome '_VS_' model_data.predictor '_' comparison_label];
    
    % create figure and axes
    colors = lines(number_of_conditions_in_this_comparison);

    figure;
    r_square_axes = axes('position', [0.08 0.5, 0.84 0.4], 'fontsize', 12); 
    hold on;
    ylabel('R^2');
    set(r_square_axes, 'xtick', [])
    ylim([0 1]);
    if arguments.show_legend
        legend('Location', 'best')
    end

    slope_axes = axes('position', [0.08 0.08, 0.84 0.4], 'fontsize', 12);
    hold on;
    xlabel('time (%)'); ylabel('slope');
    
    uicontrol('style', 'text', 'string', title_label, 'units', 'normalized', 'position', [0, 0.95, 1, 0.05], 'fontsize', 16, 'FontWeight', 'bold');

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
            model_data.row_info, ...
            model_data.row_info_headers, ...
            levels_to_remove ...
          );
        
        
        this_condition_label = strrep(model_data.row_info{this_condition_indicator, relevant_column}, '_', ' ');
        
        r_square_data_here = model_data.R_square{this_condition_indicator};
        slope_data_here = model_data.slope{this_condition_indicator};
        slope_cinv_here = model_data.slope_confidence_interval_width{this_condition_indicator};

        this_condition = comparisons.condition_combinations(this_condition_index, :);
        this_label = this_condition{strcmp(comparisons.condition_combination_labels, comparisons.condition_to_compare)};
        this_color = comparisons.condition_colors{strcmp(comparisons.condition_colors(:, 1), this_label), 2};
        
        plot(r_square_axes, 0:100, r_square_data_here, 'linewidth', 2, 'DisplayName', this_condition_label, 'color', this_color);
        shadedErrorBar(0:100, slope_data_here, slope_cinv_here, {'linewidth', 2, 'color', this_color}, 1, slope_axes);
    end

    % plot zero line
    if linear_model_settings.get('plot_zero', 1)
        xlimits = get(slope_axes, 'xlim');
        zero_plot = plot(slope_axes, xlimits, [0 0], 'color', [0.7 0.7 0.7], 'linewidth', 2);
        set(zero_plot, 'HandleVisibility', 'off');
        uistack(zero_plot, 'bottom')
    end

    % save
    if arguments.save_results
        if ~directoryExists('figures')
            mkdir figures
        end
        filename = ['figures' filesep file_label '.jpg'];
        print(gcf, filename, '-r300', '-djpeg')
    end    
    
end




