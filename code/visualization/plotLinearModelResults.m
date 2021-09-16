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

function plotLinearModelResults(varargin)
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
    model_file_name = ['results' filesep collection_date '_' subject_id '_linearModels.mat'];
    model_data = load(model_file_name);
    model_list = linear_model_settings.getTable('plot_table');
    
    preferred_level_order = linear_model_settings.get('preferred_level_order', 1);
    for i_model = 1 : size(model_list, 1)
        this_model = model_list(i_model, :);
        
        % find the requested model in the data
        this_model_index = findModelIndex(model_data, this_model, linear_model_settings);
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
            createComparisonFigure(this_model_data, comparisons, i_comparison, relevant_column, linear_model_settings, arguments);
        end
    end

end

function index = findModelIndex(data, model, settings)
    % extract info
    requested_outcome_variable = model.outcome_variable_name{1};
    requested_predictor_variable_list_name = model.predictor_variable_list{1};
    
    % get list of predictor variables
    requested_predictor_variable_list = settings.get(requested_predictor_variable_list_name, 1);
    if isempty(requested_predictor_variable_list)
        % failed to load list from settings, assume this is already a variable name
        requested_predictor_variable_list = {requested_predictor_variable_list_name};
    end
%     requested_predictor_variable_list = settings.get(requested_predictor_variable_list_name);

    predictor_column = strcmp(data.linear_model_results_header, 'predictor_variables');
    outcome_column = strcmp(data.linear_model_results_header, 'outcome_variable');
    
    
    % loop through models to find the requested one
    index = 0;
    for i_model = 1 : size(data.linear_model_results, 1)
        this_model_outcome_variable = data.linear_model_results{i_model, outcome_column};
        this_model_predictor_variables = data.linear_model_results{i_model, predictor_column};
        
        if isequal(this_model_outcome_variable, requested_outcome_variable) ...
            && isequal(this_model_predictor_variables, requested_predictor_variable_list)
            index = i_model;
        end
    end
    
    if index == 0
        error(['Model not available'])
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
    comparison_label = strrep(relevant_row_info{1}, '_', ' ');
    for i_column = 2 : length(relevant_row_info)
        comparison_label = [comparison_label ', ' strrep(relevant_row_info{i_column}, '_', ' ')]; %#ok<AGROW>
    end
end

function createComparisonFigure(model_data, comparisons, comparison_to_show, relevant_column, linear_model_settings, arguments)
    comparison_indices = comparisons.comparison_indices{comparison_to_show};
    outcome_label = model_data.names.outcome;
    comparison_label = createComparisonLabel(comparisons.condition_combinations, comparison_indices, relevant_column);
    number_of_conditions_in_this_comparison = length(comparison_indices);
    number_of_predictors = length(model_data.names.predictors);
    title_label = ['Outcome: ' outcome_label ', ' comparison_label];
    file_label = [model_data.names.outcome '_VS_' model_data.names.predictors_label '_' comparison_label];
    
    % create figure and axes
    figure;
    tiledlayout(number_of_predictors + 1, 1);

    for i_predictor = 1 : number_of_predictors
        slope_axes(i_predictor) = nexttile; %#ok<AGROW>
        set(slope_axes(i_predictor), 'fontsize', 12);
        set(slope_axes(i_predictor), 'xtick', [])
        ylabel('slope');
        hold on;
        title(strrep(model_data.names.predictors{i_predictor}, '_', ' '), 'Units', 'normalized', 'Position', [0.5, 0.9, 0])
    end
    
    r_square_axes = nexttile;
    set(r_square_axes, 'fontsize', 12);
    hold on;
    ylabel('$R^2$', 'Interpreter', 'latex');
    ylim([0 1]);
    if arguments.show_legend
        legend('Location', 'southeastoutside')
    end
    
    % define parameters depending on model type, discrete vs. continuous predictor
    if model_data.predictor_variable_data_points_per_stretch == 1
        linestyle = 'none';
        linewidth = 1;
        marker = 's';
        marker_size = 15;
        x_label = comparisons.condition_to_compare;
        set(r_square_axes, 'XTickLabel', {}, 'xtick', []);
    else
        linestyle = '-';
        linewidth = 2;
        marker = 'none';
        marker_size = 1;
        x_label = 'time (\%)';
    end
    
    % add labels
    xlabel(r_square_axes, x_label); 
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
        if model_data.predictor_variable_data_points_per_stretch == 1
            abscissa_data = i_condition;
            r_square_axes.XTick = [r_square_axes.XTick, i_condition];
            r_square_axes.XTickLabel = [r_square_axes.XTickLabel; this_condition_label];
        else
            abscissa_data = (1 : model_data.predictor_variable_data_points_per_stretch) - 1;
        end
        r_square_data_here = model_data.R_square{this_condition_indicator};
        slope_data_here = model_data.slope{this_condition_indicator};

        this_condition = comparisons.condition_combinations(this_condition_index, :);
        this_label = this_condition{strcmp(comparisons.condition_combination_labels, comparisons.condition_to_compare)};
        this_color = comparisons.condition_colors{strcmp(comparisons.condition_colors(:, 1), this_label), 2};
        
        plot ...
          ( ...
            r_square_axes, ...
            abscissa_data, r_square_data_here, ...
            'LineStyle', linestyle, ...
            'marker', marker, ...
            'linewidth', linewidth, ...
            'DisplayName', this_condition_label, ...
            'color', this_color, ...
            'MarkerFaceColor', this_color, ...
            'MarkerSize', marker_size, ...
            'DisplayName', this_condition_label ...
          );
        
        for i_predictor = 1 : number_of_predictors
            plot ...
              ( ...
                slope_axes(i_predictor), ...
                abscissa_data, slope_data_here(:, i_predictor), ...
                'LineStyle', linestyle, ...
                'marker', marker, ...
                'linewidth', linewidth, ...
                'DisplayName', this_condition_label, ...
                'color', this_color, ...
                'MarkerFaceColor', this_color, ...
                'MarkerSize', marker_size, ...
                'DisplayName', this_condition_label ...
              );
        end
    end

    % adjust x-limits
    if model_data.predictor_variable_data_points_per_stretch == 1
        xlimits = [r_square_axes.XTick(1)-0.5, r_square_axes.XTick(end)+0.5];
        set(r_square_axes, 'xlim', xlimits);
        for i_predictor = 1 : number_of_predictors
            set(slope_axes(i_predictor), 'xlim', xlimits);
        end
    end
    
    % plot zero line
    if linear_model_settings.get('plot_zero', 1)
        for i_predictor = 1 : number_of_predictors
            xlimits = get(slope_axes(i_predictor), 'xlim');
            zero_plot = plot(slope_axes(i_predictor), xlimits, [0 0], 'color', [0.7 0.7 0.7], 'linewidth', 2);
            set(zero_plot, 'HandleVisibility', 'off');
            uistack(zero_plot, 'bottom')
        end
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




