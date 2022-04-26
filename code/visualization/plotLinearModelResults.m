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
    addParameter(parser, 'settings', 'linearModel')
    addParameter(parser, 'show_legend', true)
    addParameter(parser, 'save', false)
    parse(parser, varargin{:})
    arguments.show_legend = parser.Results.show_legend;
    arguments.save_results = parser.Results.save;
    arguments.settings = parser.Results.settings;

    % load settings
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    linear_model_settings = loadSettingsFromFile(arguments.settings);
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    condition_to_compare = linear_model_settings.get('condition_to_compare');

    % load data
    model_file_name = ['results' filesep collection_date '_' subject_id '_linearModels.mat'];
    model_data = load(model_file_name);
    model_list = linear_model_settings.getTable('plot_table');
    
    preferred_level_order = linear_model_settings.get('preferred_level_order', 1);
    for i_model = 1 : size(model_list, 1)
        
        % find the requested model in the data
        this_model_label = model_list.label{i_model};
        this_model_index = findModelIndex(model_data, this_model_label);
        this_model_data = model_data.linear_model_results{this_model_index, strcmp(model_data.linear_model_results_header, 'results')};
        
        % create figures
        relevant_column = strcmp(this_model_data.row_info_headers, condition_to_compare);
        
        % determine comparisons
        condition_data = this_model_data.row_info;
        conditions_settings = study_settings.get('conditions');
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
            createComparisonFigure(this_model_data, comparisons, i_comparison, relevant_column, linear_model_settings, arguments);
        end
    end

end

function index = findModelIndex_old(data, model, settings)
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

function index = findModelIndex(model_data, requested_label)
    % loop through models to find the requested one
    label_column = strcmp(model_data.linear_model_results_header, 'label');
    index = strcmp(model_data.linear_model_results(:, label_column), requested_label);
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
    
    % create summary figure and axes for slopes and R^2
    summary_figure = figure;
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
    ylabel('R^2');
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
        x_label = 'time (%)';
    end
    
    % add labels
    xlabel(r_square_axes, x_label); 
    uicontrol('style', 'text', 'string', title_label, 'units', 'normalized', 'position', [0, 0.95, 1, 0.05], 'fontsize', 16, 'FontWeight', 'bold');
    
    % create figure for detailed data
    if linear_model_settings.get('show_details', 1)
        details_figure = figure('SizeChangedFcn', @sizeChanged, 'visible', 'off');
        tiledlayout(ceil(number_of_predictors/2), 2);
        
        predictor_directions = model_data.directions.predictor;
        outcome_directions = model_data.directions.outcome;
        
        pos_text_handles_abscissa = zeros(number_of_predictors, 1);
        pos_arrow_handles_abscissa = zeros(number_of_predictors, 1);
        neg_text_handles_abscissa = zeros(number_of_predictors, 1);
        neg_arrow_handles_abscissa = zeros(number_of_predictors, 1);
        pos_text_handles_ordinate = zeros(number_of_predictors, 1);
        pos_arrow_handles_ordinate = zeros(number_of_predictors, 1);
        neg_text_handles_ordinate = zeros(number_of_predictors, 1);
        neg_arrow_handles_ordinate = zeros(number_of_predictors, 1);

        for i_predictor = 1 : number_of_predictors
            details_axes(i_predictor) = nexttile; %#ok<AGROW>
            set(details_axes(i_predictor), 'fontsize', 12);
            xlabel(strrep(model_data.names.predictors{i_predictor}, '_', ' '));
            ylabel(strrep(outcome_label, '_', ' '));
            hold on;
            title(strrep(model_data.names.predictors{i_predictor}, '_', ' '))
            
            %
            fontsize_text = 12;
            fontsize_arrow = 18;
            
            % add text labels for arrows - x axis
            pos_text_handles_abscissa(i_predictor) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    predictor_directions{i_predictor, 1}, ...
                    'rotation', 0, ...
                    'Fontsize', fontsize_text, ...
                    'horizontalalignment', 'right', ...
                    'parent', details_axes(i_predictor) ...
                  );                
            pos_arrow_handles_abscissa(i_predictor) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    ' $\rightarrow$', ...
                    'rotation', 0, ...
                    'Fontsize', fontsize_arrow, ...
                    'horizontalalignment', 'right', ...
                    'interpreter', 'LaTeX', ...
                    'parent', details_axes(i_predictor) ...
                  );
            neg_text_handles_abscissa(i_predictor) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    predictor_directions{i_predictor, 2}, ...
                    'rotation', 0, ...
                    'Fontsize', fontsize_text, ...
                    'horizontalalignment', 'left', ...
                    'parent', details_axes(i_predictor)...
                  );
            neg_arrow_handles_abscissa(i_predictor) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    '$\leftarrow$ ', ...
                    'rotation', 0, ...
                    'Fontsize', fontsize_arrow, ...
                    'horizontalalignment', 'left', ...
                    'interpreter', 'LaTeX', ...
                    'parent', details_axes(i_predictor) ...
                  );
              
            % add text labels for arrows - y axis
            pos_text_handles_ordinate(i_predictor) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    outcome_directions{1}, ...
                    'rotation', 90, ...
                    'Fontsize', fontsize_text, ...
                    'horizontalalignment', 'right', ...
                    'parent', details_axes(i_predictor) ...
                  );                
            pos_arrow_handles_ordinate(i_predictor) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    ' $\rightarrow$', ...
                    'rotation', 90, ...
                    'Fontsize', fontsize_arrow, ...
                    'horizontalalignment', 'right', ...
                    'interpreter', 'LaTeX', ...
                    'parent', details_axes(i_predictor) ...
                  );
            neg_text_handles_ordinate(i_predictor) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    outcome_directions{2}, ...
                    'rotation', 90, ...
                    'Fontsize', fontsize_text, ...
                    'horizontalalignment', 'left', ...
                    'parent', details_axes(i_predictor)...
                  );
            neg_arrow_handles_ordinate(i_predictor) = ...
                text ...
                  ( ...
                    0, ...
                    0, ...
                    '$\leftarrow$ ', ...
                    'rotation', 90, ...
                    'Fontsize', fontsize_arrow, ...
                    'horizontalalignment', 'left', ...
                    'interpreter', 'LaTeX', ...
                    'parent', details_axes(i_predictor) ...
                  );
                     
        end
        
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
            model_data.row_info, ...
            model_data.row_info_headers, ...
            levels_to_remove ...
          );
        
        this_condition_label = strrep(model_data.row_info{find(this_condition_indicator, 1), relevant_column}, '_', ' ');
        
        % get data for r_square and slope
        r_square_data_here = model_data.R_square{this_condition_indicator};
        slope_data_here = model_data.slope{this_condition_indicator};
        
        if model_data.predictor_variable_data_points_per_stretch == 1
            abscissa_data = i_condition;
            r_square_axes.XTick = [r_square_axes.XTick, i_condition];
            r_square_axes.XTickLabel = [r_square_axes.XTickLabel; this_condition_label];
            
            if linear_model_settings.get('show_details', 1)
                % make sure we're dealing with a single trial
                if sum(this_condition_indicator) > 1
                    error('Detailed data can only be shown for single trials, please do not pool across conditions.')
                end
                
                % get predictor and outcome data
                predictor_data_here = model_data.data.predictors(this_condition_indicator, :);
                outcome_data_here = model_data.data.outcome{this_condition_indicator, :};
                
                % create linear model data
                predictor_center = model_data.predictor_offsets{this_condition_indicator}';
                outcome_center = model_data.outcome_offsets{this_condition_indicator}';
                predictor_min = zeros(number_of_predictors, 1);
                predictor_max = zeros(number_of_predictors, 1);
                for i_predictor = 1 : number_of_predictors
                    predictor_min(i_predictor) = min(predictor_data_here{i_predictor});
                    predictor_max(i_predictor) = max(predictor_data_here{i_predictor});
                end
                predictor_data_for_plotting = [predictor_min predictor_max];
                predictor_data_for_plotting_centered = predictor_data_for_plotting - predictor_center;
                outcome_data_for_plotting = outcome_center + slope_data_here * predictor_data_for_plotting_centered;
            end
        else
            abscissa_data = (1 : model_data.predictor_variable_data_points_per_stretch) - 1;
        end

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
            % plot slopes
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
          
            % plot details
            if linear_model_settings.get('show_details', 1)
                % plot individual data
                this_predictor_data = predictor_data_here{i_predictor};
                plot ...
                  ( ...
                    details_axes(i_predictor), ...
                    this_predictor_data, outcome_data_here, ...
                    'LineStyle', 'none', ...
                    'marker', 'o', ...
                    'DisplayName', this_condition_label, ...
                    'MarkerEdgeColor', 'none', ...
                    'MarkerFaceColor', lightenColor(this_color, 0.5), ...
                    'MarkerSize', 6, ...
                    'DisplayName', this_condition_label ...
                  );
                
                % plot linear model
                plot ...
                  ( ...
                    details_axes(i_predictor), ...
                    predictor_data_for_plotting(i_predictor, :), outcome_data_for_plotting, ...
                    'LineStyle', '-', ...
                    'DisplayName', this_condition_label, ...
                    'Color', this_color, ...
                    'linewidth', 6, ...
                    'DisplayName', this_condition_label ...
                  );
            end
        end
    end

    % show figure
    if linear_model_settings.get('show_details', 1)
        set(details_figure, 'visible', 'on')
    end
    
    % adjust x-limits
    if model_data.predictor_variable_data_points_per_stretch == 1
        xlimit_conditions = [r_square_axes.XTick(1)-0.5, r_square_axes.XTick(end)+0.5];
        set(r_square_axes, 'xlim', xlimit_conditions);
        for i_predictor = 1 : number_of_predictors
            set(slope_axes(i_predictor), 'xlim', xlimit_conditions);
            if linear_model_settings.get('show_details', 1)
                sizeChanged();
  
            end
        
        end
    end
    
    % plot zero line
    if linear_model_settings.get('plot_zero', 1)
        for i_predictor = 1 : number_of_predictors
            xlimits_predictor = get(slope_axes(i_predictor), 'xlim');
            zero_plot = plot(slope_axes(i_predictor), xlimits_predictor, [0 0], 'color', [0.7 0.7 0.7], 'linewidth', 2);
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
    
    function sizeChanged(varargin)
        for j_predictor = 1 : number_of_predictors
            these_axes = details_axes(j_predictor);
            xlimits = get(these_axes, 'xlim'); ylimits = get(these_axes, 'ylim');

            pos_arrow_position_x = xlimits(2);
            pos_arrow_position_y = ylimits(1) - (ylimits(2)-ylimits(1))*0.09;
            set(pos_arrow_handles_abscissa(j_predictor), 'Position', [pos_arrow_position_x pos_arrow_position_y]);
            pos_text_position_x = xlimits(2);
            pos_text_position_y = ylimits(1) - (ylimits(2)-ylimits(1))*0.12;
            set(pos_text_handles_abscissa(j_predictor), 'Position', [pos_text_position_x pos_text_position_y]);

            neg_arrow_position_x = xlimits(1);
            neg_arrow_position_y = ylimits(1) - (ylimits(2)-ylimits(1))*0.09;
            set(neg_arrow_handles_abscissa(j_predictor), 'Position', [neg_arrow_position_x neg_arrow_position_y]);
            neg_text_position_x = xlimits(1);
            neg_text_position_y = ylimits(1) - (ylimits(2)-ylimits(1))*0.12;
            set(neg_text_handles_abscissa(j_predictor), 'Position', [neg_text_position_x neg_text_position_y]);

            pos_arrow_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.09;
            pos_arrow_position_y = ylimits(2);
            set(pos_arrow_handles_ordinate(j_predictor), 'Position', [pos_arrow_position_x pos_arrow_position_y]);
            pos_text_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.12;
            pos_text_position_y = ylimits(2);
            set(pos_text_handles_ordinate(j_predictor), 'Position', [pos_text_position_x pos_text_position_y]);

            neg_arrow_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.09;
            neg_arrow_position_y = ylimits(1);
            set(neg_arrow_handles_ordinate(j_predictor), 'Position', [neg_arrow_position_x neg_arrow_position_y]);
            neg_text_position_x = xlimits(1) - (xlimits(2)-xlimits(1))*0.12;
            neg_text_position_y = ylimits(1);
            set(neg_text_handles_ordinate(j_predictor), 'Position', [neg_text_position_x neg_text_position_y]);     
        end
    end
    
end




