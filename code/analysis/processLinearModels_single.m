
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

function processLinearModels_single(varargin)
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
    [condition_combinations_unique, condition_indicators] = getUniqueConditionInformation(condition_data_all, condition_labels);
    number_of_condition_combinations = size(condition_combinations_unique, 1);
    
    % get table with analysis information
    model_table = linear_model_settings.getTable('models_single');
    number_of_variables = height(model_table);
    
    % analyze
    linear_model_results_header = {'results', 'predictor_variable_name', 'outcome_variable_name'};
    linear_model_results = cell(number_of_variables, 3);
    for i_model = 1 : number_of_variables
        predictor_variable_name = model_table.predictor_variable_name{i_model};
        predictor_variable_data = data.variable_data{strcmp(data.variable_names, predictor_variable_name)};
        predictor_variable_directions = data.variable_directions(strcmp(data.variable_names, predictor_variable_name), :);
        outcome_variable_name = model_table.outcome_variable_name{i_model};
        outcome_variable_data = data.variable_data{strcmp(data.variable_names, outcome_variable_name)};
        outcome_variable_directions = data.variable_directions(strcmp(data.variable_names, outcome_variable_name), :);
        predictor_variable_data_points_per_stretch = size(predictor_variable_data, 1);
        
        % make containers for results
        results_this_model = struct;
        results_this_model.predictor = predictor_variable_name;
        results_this_model.outcome = outcome_variable_name;
        results_this_model.data.predictor = cell(number_of_condition_combinations, 1);
        results_this_model.data.outcome = cell(number_of_condition_combinations, 1);
        results_this_model.directions.predictor = predictor_variable_directions;
        results_this_model.directions.outcome = outcome_variable_directions;
        results_this_model.R_square = cell(number_of_condition_combinations, 1);
        results_this_model.slope = cell(number_of_condition_combinations, 1);
        results_this_model.slope_confidence_interval_width = cell(number_of_condition_combinations, 1);
        results_this_model.offset{i_condition} = cell(number_of_condition_combinations, 1);
        results_this_model.offset_confidence_interval_width{i_condition} = cell(number_of_condition_combinations, 1);
        results_this_model.correlation_p = cell(number_of_condition_combinations, 1);
        results_this_model.row_info = condition_combinations_unique;
        results_this_model.row_info_headers = condition_labels;
        results_this_model.predictor_variable_data_points_per_stretch = predictor_variable_data_points_per_stretch;
        
        % loop through conditions
        for i_condition = 1 : number_of_condition_combinations
            % extract data from this condition
%             this_condition_labels = condition_combinations_unique(i_condition, :); % probably not needed
            this_condition_indicator = condition_indicators(:, i_condition);
        
            predictor_variable_data_this_condition = predictor_variable_data(:, this_condition_indicator);
            outcome_variable_data_this_condition = outcome_variable_data(:, this_condition_indicator);
        
            % create containers
            R_square_table_here = zeros(predictor_variable_data_points_per_stretch, 1) * NaN;
            slope_table_here = zeros(predictor_variable_data_points_per_stretch, 1) * NaN;
            slope_confidence_interval_width_table_here = zeros(predictor_variable_data_points_per_stretch, 1) * NaN;
            offset_table_here = zeros(predictor_variable_data_points_per_stretch, 1) * NaN;
            offset_confidence_interval_width_table_here = zeros(predictor_variable_data_points_per_stretch, 1) * NaN;
            correlation_p_table_here = zeros(predictor_variable_data_points_per_stretch, 1) * NaN;
            
            % loop through stretches in predictor variable
            for i_stretch = 1 : predictor_variable_data_points_per_stretch
                predictor = predictor_variable_data_this_condition(i_stretch, :)';
                outcome = outcome_variable_data_this_condition';
                if ~any(isnan(predictor)) || any(isnan(outcome))
                    % fit
                    [fit_object, fit_stats] = fit(predictor, outcome, 'poly1');
                    [~, correlation_p] = corrcoef(predictor, outcome);

                    % store
                    R_square_table_here(i_stretch) = fit_stats.rsquare;
                    slope_table_here(i_stretch) = fit_object.p1;
                    offset_table_here(i_stretch) = fit_object.p2;
                    confidence_intervals = confint(fit_object);
                    slope_confidence_interval_width_table_here(i_stretch) = range(confidence_intervals(:, 1));
                    offset_confidence_interval_width_table_here(i_stretch) = range(confidence_intervals(:, 2));
                    correlation_p_table_here(i_stretch) = correlation_p(1, 2);
                    
%                     % plot
%                     figure; 
%                     plot(fit_object, predictor', outcome')
%                     side_label = strrep(this_condition_labels{2}, '_', ' '); % TODO: hard-coded for now, solve this properly later
%                     metronome_label = strrep(this_condition_labels{3}, '_', ' '); % TODO: hard-coded for now, solve this properly later
%                     title_string = [metronome_label ' - ' side_label ' - R^2 = ' num2str(fit_stats.rsquare) ', slope = ' num2str(fit_object.p1)];
%                     title(title_string);
% 
%                     xlabel(strrep(predictor_variable_name, '_', ' '));
%                     ylabel(strrep(outcome_variable_name, '_', ' '));
%                     
                end
            end
            
            % store
            results_this_model.data.predictor{i_condition} = predictor_variable_data_this_condition;
            results_this_model.data.outcome{i_condition} = outcome_variable_data_this_condition;
            results_this_model.R_square{i_condition} = R_square_table_here;
            results_this_model.slope{i_condition} = slope_table_here;
            results_this_model.slope_confidence_interval_width{i_condition} = slope_confidence_interval_width_table_here;
            results_this_model.offset{i_condition} = offset_table_here;
            results_this_model.offset_confidence_interval_width{i_condition} = offset_confidence_interval_width_table_here;
            results_this_model.correlation_p{i_condition} = correlation_p_table_here;


%             % calculate variance and standard deviation
%             std_predictor = std(predictor);
%             std_outcome = std(outcome);
%             std_residual = std(residual);
% 
%             var_predictor = var(predictor);
%             var_outcome = var(outcome);
%             var_residual = var(residual);

%             disp([relevant_condition_label ' - ' variable_B_name ' STD: ' num2str(std_outcome)])
%             disp([relevant_condition_label ' - ' variable_B_name ' Variance: ' num2str(var_outcome)])
            
        end
        linear_model_results{i_model, 1} = results_this_model;
        linear_model_results{i_model, 2} = predictor_variable_name;
        linear_model_results{i_model, 3} = outcome_variable_name;
        disp(['Finished model ' num2str(i_model) ' of ' num2str(number_of_variables) ', predictor: ' predictor_variable_name ', outcome: ' outcome_variable_name])
        
    end
    
    % save
    results_file_name = ['results' filesep collection_date '_' subject_id '_linearModels.mat'];
    save(results_file_name, 'linear_model_results', 'linear_model_results_header');
end    
    
    









