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


% this function loads the linear model variables from results/..._linearModelVariables.mat and fits linear models for 
% each pair of predictor and outcome variables listed in linearModelVariables.mat

function processLinearModels(varargin)
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    linear_model_settings = loadSettingsFromFile('linearModel');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');

    % load settings and existing results
    data_file_name = ['results' filesep makeFileName(collection_date, subject_id, 'linearModelVariables')];
    data = load(data_file_name);
    
    % make condition data tables
    conditions.settings_table = study_settings.get('conditions');
    conditions.factor_labels = conditions.settings_table(:, 1)';
    conditions.source_variables = conditions.settings_table(:, 2)';
    conditions.number_of_factor_labels = length(conditions.factor_labels);
    conditions.conditions_session = data.conditions;
    
    % get list of conditions
    conditions_settings = study_settings.get('conditions');
    condition_labels = conditions_settings(:, 1)';
    condition_source_variables = conditions_settings(:, 2)';
    number_of_condition_labels = length(condition_labels);
    conditions_session = conditions.conditions_session;
    number_of_stretches = data.number_of_stretches; %#ok<*USENS>
    condition_data_all = cell(number_of_stretches, number_of_condition_labels);
    for i_condition = 1 : number_of_condition_labels
        condition_data_all(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
    end
    
    % fit first-order models
    [linear_model_results, linear_model_results_header] = fitLinearModels(data, linear_model_settings, 'models', condition_data_all, condition_labels);
    
    % calculate prediction errors
    prediction_error_table = linear_model_settings.getTable('model_prediction_error');
    [data, error_data] = calculateLinearModelPredictionErrors(data, linear_model_results, linear_model_results_header, linear_model_settings, prediction_error_table, condition_data_all, condition_labels);
    
    % save errors
    stretch_data_file_name = ['results' filesep makeFileName(collection_date, subject_id, 'results')];
    stretch_data = load(stretch_data_file_name);
    for i_error = 1 : size(error_data.variable_data, 1)
        new_data = struct;
        new_data.data = error_data.variable_data{1, :};
        new_data.directions = error_data.variable_directions(i_error, :);
        new_data.name = error_data.variable_names{i_error, :};
        
        stretch_data = addOrReplaceResultsData(stretch_data, new_data, 'range');
    end
    save(stretch_data_file_name, '-struct', 'stretch_data');    

    % fit second-order models
    linear_model_results_2 = fitLinearModels(data, linear_model_settings, 'models_second_order', condition_data_all, condition_labels);
    
    % merge
    linear_model_results = [linear_model_results; linear_model_results_2];
    
    % save
    results_file_name = ['results' filesep collection_date '_' subject_id '_linearModels.mat'];
    save(results_file_name, 'linear_model_results', 'linear_model_results_header');
end    
    
    









