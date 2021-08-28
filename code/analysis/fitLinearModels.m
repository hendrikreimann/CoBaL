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


function [linear_model_results, linear_model_results_header] = fitLinearModels(data, linear_model_settings, models_to_fit, condition_data_all, condition_labels)
    % get table with analysis information
    model_table = linear_model_settings.getTable(models_to_fit);
    number_of_models = height(model_table);
    [condition_combinations_unique, condition_indicators] = getUniqueConditionInformation(condition_data_all, condition_labels);
    number_of_condition_combinations = size(condition_combinations_unique, 1);
    
    % analyze
    linear_model_results_header = {'results', 'outcome_variable', 'predictor_variables'};
    linear_model_results = cell(number_of_models, 3);
    for i_model = 1 : number_of_models
        predictor_variable_list_name = model_table.predictor_variable_list{i_model};
        predictor_variable_list = linear_model_settings.get(predictor_variable_list_name);
        number_of_predictor_variables = length(predictor_variable_list);
        predictor_variable_data = cell(number_of_predictor_variables, 1);
        predictor_variable_directions = cell(number_of_predictor_variables, 2);
        predictor_variable_data_points_per_stretch_all = zeros(number_of_predictor_variables, 1);
        for i_predictor = 1 : number_of_predictor_variables
            this_predictor_variable_name = predictor_variable_list{i_predictor};
            predictor_variable_data{i_predictor} = data.variable_data{strcmp(data.variable_names, this_predictor_variable_name)};
            predictor_variable_directions(i_predictor, :) = data.variable_directions(strcmp(data.variable_names, this_predictor_variable_name), :);
            predictor_variable_data_points_per_stretch_all(i_predictor) = size(predictor_variable_data{i_predictor}, 1);
        end
        if numel(unique(predictor_variable_data_points_per_stretch_all)) > 1
            error('Predictor variables have different numbers of data points');
        end
        predictor_variable_data_points_per_stretch = predictor_variable_data_points_per_stretch_all(1);
        outcome_variable_name = model_table.outcome_variable_name{i_model};
        outcome_variable_data = data.variable_data{strcmp(data.variable_names, outcome_variable_name)};
        outcome_variable_directions = data.variable_directions(strcmp(data.variable_names, outcome_variable_name), :);
        
        % make containers for results
        results_this_model = struct;
        results_this_model.names.predictors_label = predictor_variable_list_name;
        results_this_model.names.predictors = predictor_variable_list;
        results_this_model.names.outcome = outcome_variable_name;
        results_this_model.data.predictors = cell(number_of_condition_combinations, number_of_predictor_variables);
        results_this_model.data.outcome = cell(number_of_condition_combinations, 1);
        results_this_model.directions.predictor = predictor_variable_directions;
        results_this_model.directions.outcome = outcome_variable_directions;
        results_this_model.R_square = cell(number_of_condition_combinations, 1);
        results_this_model.slope = cell(number_of_condition_combinations, 1);
        results_this_model.row_info = condition_combinations_unique;
        results_this_model.row_info_headers = condition_labels;
        results_this_model.predictor_variable_data_points_per_stretch = predictor_variable_data_points_per_stretch;
        
        % loop through conditions
        for i_condition = 1 : number_of_condition_combinations
            % extract data from this condition
            this_condition_indicator = condition_indicators(:, i_condition);
            number_of_data_points_this_condition = sum(this_condition_indicator);
        
            predictor_variable_data_this_condition = cell(number_of_predictor_variables, 1);
            for i_variable = 1 : number_of_predictor_variables
                this_variable_data_this_condition = predictor_variable_data{i_variable}(:, this_condition_indicator);
                % subtract mean
                this_variable_data_this_condition_mean_free = zeros(size(this_variable_data_this_condition));
                for i_time = 1 : predictor_variable_data_points_per_stretch
                    data_this_time_point = this_variable_data_this_condition(i_time, :);
                    data_this_time_point_mean_free = data_this_time_point - mean(data_this_time_point);
                    this_variable_data_this_condition_mean_free(i_time, :) = data_this_time_point_mean_free;
                end
                predictor_variable_data_this_condition{i_variable} = this_variable_data_this_condition_mean_free;
            end
            outcome_variable_data_this_condition = outcome_variable_data(:, this_condition_indicator) - mean(outcome_variable_data(:, this_condition_indicator));
        
            % create containers
            R_square_table_here = zeros(predictor_variable_data_points_per_stretch, 1) * NaN;
            slope_table_here = zeros(predictor_variable_data_points_per_stretch, number_of_predictor_variables) * NaN;
            
            % loop through stretches in predictor variable
            for i_time = 1 : predictor_variable_data_points_per_stretch
                % extract data for predictor and outcome
                predictor = zeros(number_of_predictor_variables, number_of_data_points_this_condition);
                for i_variable = 1 : number_of_predictor_variables
                    predictor(i_variable, :) = predictor_variable_data_this_condition{i_variable}(i_time, :);
                end
                predictor = predictor';
                outcome = outcome_variable_data_this_condition';
                if ~any(any(isnan(predictor))) || any(isnan(outcome))
                    % fit
                    X = 
                    
                    jacobian = predictor\outcome;
                    
                    % calculate R^2
                    prediction = predictor*jacobian;
                    SS_residual = sum((outcome - prediction).^2);
                    SS_total = sum((outcome - mean(outcome)).^2);
                    R_square = 1 - SS_residual/SS_total;
                    
                    % store
                    R_square_table_here(i_time) = R_square;
                    slope_table_here(i_time, :) = jacobian;
                end
            end
            
            % store
            results_this_model.data.predictors(i_condition, :) = predictor_variable_data_this_condition;
            results_this_model.data.predictors(i_condition, :) = predictor_variable_data_this_condition;
            results_this_model.data.outcome{i_condition} = outcome_variable_data_this_condition;
            results_this_model.R_square{i_condition} = R_square_table_here;
            results_this_model.slope{i_condition} = slope_table_here;
        end
        linear_model_results{i_model, 1} = results_this_model;
        linear_model_results{i_model, 2} = outcome_variable_name;
        linear_model_results{i_model, 3} = predictor_variable_list;
        disp(['Finished model ' num2str(i_model) ' of ' num2str(number_of_models) ', predictors: ' predictor_variable_list_name ', outcome: ' outcome_variable_name])
        
    end



end