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

function plotLinearModelResults
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

    % loop through models
    model_list = linear_model_settings.getTable('plot_table');
    for i_model = 1 : size(model_list, 1)
        this_model = model_list(i_model, :);
        
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
            
            filename = ['figures' filesep model_label this_row_label '.pdf'];
            print(gcf, filename, '-r300', '-dpdf')
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









