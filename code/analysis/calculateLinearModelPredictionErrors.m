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


function [data, error_data] = calculateLinearModelPredictionErrors(data, model_results, model_results_header, linear_model_settings, prediction_error_table, condition_data_all, condition_labels)

    number_of_models = height(prediction_error_table);
    [~, condition_indicators] = getUniqueConditionInformation(condition_data_all, condition_labels);
    
    error_data = struct;
    error_data.variable_data = {};
    error_data.variable_names = {};
    error_data.variable_directions = {};

    % calculate errors
    for i_model = 1 : number_of_models
        % get variable data for predictors and outcome for this model
        predictor_variable_list_name = prediction_error_table.model_predictor_variable_list{i_model};
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
        predictor_variable_data_points_per_stretch = max(predictor_variable_data_points_per_stretch_all);
        outcome_variable_name = prediction_error_table.model_outcome_variable{i_model};
        outcome_variable_data = data.variable_data{strcmp(data.variable_names, outcome_variable_name)};
        outcome_variable_directions = data.variable_directions(strcmp(data.variable_names, outcome_variable_name), :);
        
        % find corresponding model in model results
        outcome_variable_matches = find(strcmp(model_results(:, strcmp(model_results_header, 'outcome_variable')), outcome_variable_name));
        this_model_index_in_results = [];
        for j_model = 1 : length(outcome_variable_matches)
            this_candidate_index = outcome_variable_matches(j_model);
            this_candidate_predictor_variables = model_results{this_candidate_index, strcmp(model_results_header, 'predictor_variables')};
            if isequal(predictor_variable_list, this_candidate_predictor_variables)
                this_model_index_in_results = this_candidate_index;
            end            
        end
        this_model_results = model_results{this_model_index_in_results, strcmp(model_results_header, 'results')};

        % loop through stretches and calculate error
        new_variable_data = zeros(predictor_variable_data_points_per_stretch, data.number_of_stretches);
        new_variable_name = prediction_error_table.new_variable_name{i_model};
        new_variable_directions = outcome_variable_directions;
        
        for i_stretch = 1 : data.number_of_stretches
            % get data for this stretch
            this_condition_indicator = condition_indicators(i_stretch, :);
            
            % loop through time
            for i_time = 1 : predictor_variable_data_points_per_stretch
                % get predictor and outcome data
                predictor_data_here = zeros(1, number_of_predictor_variables);
                for i_variable = 1 : number_of_predictor_variables
                    predictor_data_here(i_variable) = predictor_variable_data{i_variable}(i_time, i_stretch);
                end
                
                % get offsets and slope
                jacobian = this_model_results.slope{this_condition_indicator}(i_time, :)';
                predictor_offset = this_model_results.predictor_offsets{this_condition_indicator}(i_time, :);
                outcome_offset = this_model_results.outcome_offsets{this_condition_indicator}(i_time, :);
                
                % calculate error
                prediction = predictor_data_here*jacobian - predictor_offset*jacobian + outcome_offset;
                outcome_data_here = outcome_variable_data(i_stretch);
                error_here = outcome_data_here - prediction;
                
                % store
                new_variable_data(i_time, i_stretch) = error_here;
            end
        end
        
        % store
        data.variable_data = [data.variable_data; new_variable_data];
        data.variable_names = [data.variable_names; new_variable_name];
        data.variable_directions = [data.variable_directions; new_variable_directions];

        error_data.variable_data = [error_data.variable_data; new_variable_data];
        error_data.variable_names = [error_data.variable_names; new_variable_name];
        error_data.variable_directions = [error_data.variable_directions; new_variable_directions];

        disp(['Calculating prediction errors for model ' num2str(i_model) ' of ' num2str(number_of_models) ', predictors: ' predictor_variable_list_name ', outcome: ' outcome_variable_name])
    end


end