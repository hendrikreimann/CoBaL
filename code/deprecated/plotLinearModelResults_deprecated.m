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
    show_legend = parser.Results.show_legend;
    save_results = parser.Results.save;

    % load settings
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    linear_model_settings = loadSettingsFromFile('linearModel');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    condition_to_compare = linear_model_settings.get('condition_to_compare');

    % load data
    model_file_name = ['results' filesep collection_date '_' subject_id '_linearModelsFromSameBand.mat'];
    model_data = load(model_file_name);
    results_file_name = ['results' filesep collection_date '_' subject_id '_results.mat'];
%     results_data = load(results_file_name);
%     number_of_bands = results_data.bands_per_stretch;

    % plot models with continuous predictor
    model_list_continuous_predictor = linear_model_settings.getTable('plot_table_continuous_predictor');
    for i_model = 1 : size(model_list_continuous_predictor, 1)
        this_model = model_list_continuous_predictor(i_model, :);
        
        % find the requested model in the data
        this_model_index = findModelIndex(model_data, this_model);
        this_model_data = model_data.linear_model_results{this_model_index, strcmp(model_data.linear_model_results_header, 'results')};
        
        % create figures
        number_of_rows = size(this_model_data.row_info, 2);
        number_of_columns = size(this_model_data.conditions, 2);
        colors = lines(number_of_columns);
        relevant_row = strcmp(this_model_data.condition_labels, condition_to_compare);
        for i_row = 1 : number_of_rows
            this_row_label = this_model_data.row_info{i_row};
            model_label = [this_model.outcome_variable_name{1} ' vs. ' this_model.predictor_variable_name{1}];
            row_label = [ ' - ' this_row_label ' - '];
            condition_labels = this_model_data.conditions(relevant_row, :);
            
            % create figure
            figure;
            r_square_axes = axes('position', [0.08 0.5, 0.84 0.4]); 
            hold on;
            ylabel('R^2');
            set(r_square_axes, 'xtick', [])
            ylim([0 1]);
            legend('Location', 'best')
            
            slope_axes = axes('position', [0.08 0.08, 0.84 0.4]);
            hold on;
            xlabel('time (%)'); ylabel('slope');
            uicontrol('style', 'text', 'string', model_label, 'units', 'normalized', 'position', [0, 0.95, 1, 0.05], 'fontsize', 16, 'FontWeight', 'bold');
            uicontrol('style', 'text', 'string', row_label, 'units', 'normalized', 'position', [0, 0.9, 1, 0.05], 'fontsize', 16, 'FontWeight', 'bold');
            
            % loop through columns
            for i_column = 1 : number_of_columns
                this_condition_label = strrep(condition_labels{i_column}, '_', ' ');
                r_square_data_here = this_model_data.R_square{i_row, i_column};
                slope_data_here = this_model_data.slope{i_row, i_column};
                slope_cinv_here = this_model_data.slope_confidence_interval_width{i_row, i_column};
                
                this_color = colors(i_column, :);
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
            if save_results
                if ~directoryExists('figure')
                    mkdir figure
                end
                filename = ['figures' filesep model_label this_row_label '.pdf'];
                print(gcf, filename, '-r300', '-dpdf')
            end
        end
        
    end

    % plot models with discrete predictor
    model_list_discrete_predictor = linear_model_settings.getTable('plot_table_discrete_predictor');
    for i_model = 1 : size(model_list_discrete_predictor, 1)
        this_model = model_list_discrete_predictor(i_model, :);
        
        % find the requested model in the data
        this_model_index = findModelIndex(model_data, this_model);
        this_model_data = model_data.linear_model_results{this_model_index, strcmp(model_data.linear_model_results_header, 'results')};
        predictor_label = strrep(this_model.predictor_variable_name{1}, '_', ' ');
        outcome_label = strrep(this_model.outcome_variable_name{1}, '_', ' ');
        
        % create figures
        number_of_rows = size(this_model_data.row_info, 2);
        number_of_columns = size(this_model_data.conditions, 2);
        colors = lines(number_of_columns);
        relevant_row = strcmp(this_model_data.condition_labels, condition_to_compare);
        for i_row = 1 : number_of_rows
            this_row_label = this_model_data.row_info{i_row};
            model_label = [outcome_label ' vs. ' predictor_label];
            row_label = [ ' - ' this_row_label ' - '];
            condition_labels = this_model_data.conditions(relevant_row, :);
            
            % create figure and axes
            figure;
            
            slope_axes = axes('position', [0.08 0.08, 0.4 0.38]);
            hold on;
            xlabel('time (%)'); ylabel('slope');
            xlim([0.5 number_of_columns+0.5]);
            uicontrol('style', 'text', 'string', model_label, 'units', 'normalized', 'position', [0, 0.95, 1, 0.05], 'fontsize', 16, 'FontWeight', 'bold');
            uicontrol('style', 'text', 'string', row_label, 'units', 'normalized', 'position', [0, 0.9, 1, 0.05], 'fontsize', 16, 'FontWeight', 'bold');
            
            r_square_axes = axes('position', [0.08 0.5, 0.4 0.38]); 
            hold on;
            ylabel('R^2');
            set(r_square_axes, 'xtick', [])
            xlim([0.5 number_of_columns+0.5]);
            ylim([0 1]);
            set(r_square_axes, 'xtick', [])
            legend('Location', 'best')
            
            data_axes = axes('position', [0.58 0.08, 0.4 0.8]); 
            hold on;
            xlabel(predictor_label);
            ylabel(outcome_label);
            
            % loop through columns
            for i_column = 1 : number_of_columns
                this_condition_label = strrep(condition_labels{i_column}, '_', ' ');
                predictor_data_here = this_model_data.data.predictor{i_row, i_column};
                outcome_data_here = this_model_data.data.outcome{i_row, i_column};
                r_square_data_here = this_model_data.R_square{i_row, i_column};
                slope_data_here = this_model_data.slope{i_row, i_column};
                slope_cinv_here = this_model_data.slope_confidence_interval_width{i_row, i_column};
                offset_data_here = this_model_data.offset{i_row, i_column};
                
                this_color = colors(i_column, :);
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
                    i_column * [1 1], ...
                    slope_data_here + [-1 1]*slope_cinv_here, ...
                    '-', ...
                    'color', this_color, ...
                    'linewidth', 5 ...
                  );
                plot ...
                  ( ...
                    slope_axes, ...
                    i_column, ...
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
                    i_column, ...
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
            if save_results
                if ~directoryExists('figure')
                    mkdir figure
                end
                filename = ['figures' filesep model_label this_row_label '.pdf'];
                print(gcf, filename, '-r300', '-dpdf')
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









