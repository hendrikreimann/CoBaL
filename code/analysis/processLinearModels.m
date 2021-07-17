
%     This file is part of the CoBaL code base
%     Copyright (C) 2017 Hendrik Reimann <hendrikreimann@gmail.com>
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

% This function uses the previously calculated stretch variables to process analysis variables. For all stretch variables,
% the response is calculated, i.e. the difference from the control mean.

function processLinearModels(varargin)
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    linear_model_settings = loadSettingsFromFile('linearModel');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');

    % load settings and existing results
    results_file_name = ['results' filesep makeFileName(collection_date, subject_id, 'results')];
    data = load(results_file_name);
    
    % make condition data tables
    conditions.settings_table = study_settings.get('conditions');
    conditions.factor_labels = conditions.settings_table(:, 1)';
    conditions.source_variables = conditions.settings_table(:, 2)';
    conditions.number_of_factor_labels = length(conditions.factor_labels);
    conditions.conditions_session = data.conditions_session;
    
    % get table with analysis information
    analysis_table = getAnalysisTable(linear_model_settings);
    
    
    % analyze
    for i_row = 1 : height(analysis_table)
        this_action = analysis_table{i_row, 'style'};
        this_settings_table_name = analysis_table{i_row, 'settings_table'};
        this_settings_table = linear_model_settings.getTable(this_settings_table_name{1}, 1);
        
        
        if strcmp(this_action, 'relate two variables from same band')
            relateTwoVariablesFromSameBand(this_settings_table, study_settings, subject_settings, data, conditions);
        end        
        
    end
    
end    
    
    
function relateTwoVariablesFromSameBand(variable_table, study_settings, subject_settings, data, conditions)
    % get list of conditions
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    conditions_settings = study_settings.get('conditions');
    condition_labels = conditions_settings(:, 1)';
    condition_source_variables = conditions_settings(:, 2)';
    number_of_condition_labels = length(condition_labels);
    conditions_session = conditions.conditions_session;
    number_of_stretches = size(data.time_list_session, 1); %#ok<*USENS>
    condition_data_all = cell(number_of_stretches, number_of_condition_labels);
    for i_condition = 1 : number_of_condition_labels
        condition_data_all(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
    end
%     labels_to_ignore = study_settings.get('conditions_to_ignore');
%     levels_to_remove = study_settings.get('levels_to_remove', 1);
    [condition_combinations_unique, condition_indicators] = getUniqueConditionInformation(condition_data_all, condition_labels);
    
    % prepare
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    number_of_time_steps_per_continuous_stretch = (number_of_time_steps_normalized-1) * data.bands_per_stretch + 1;
    number_of_bands = data.bands_per_stretch;
    number_of_condition_combinations = size(condition_combinations_unique, 1);
    number_of_variables = size(variable_table, 1);
    
    % loop through variables
    linear_model_results_header = {'results', 'predictor_variable_name', 'predictor_variable_type', 'outcome_variable_name', 'outcome_variable_type'};
    linear_model_results = cell(number_of_variables, 5);
    disp('Processing linear models to relate two variables from same band')
    for i_variable = 1 : number_of_variables
        % get data
        variable_A_name = variable_table.variable_A_name{i_variable};
        variable_A_type = variable_table.variable_A_type{i_variable};
        variable_B_name = variable_table.variable_B_name{i_variable};
        variable_B_type = variable_table.variable_B_type{i_variable};
        
        % pick data depending on source specification
        variable_A_data_source = data.([variable_A_type '_data_session']);
        variable_A_names_source = data.([variable_A_type '_names_session']);
        variable_B_data_source = data.([variable_B_type '_data_session']);
        variable_B_names_source = data.([variable_B_type '_names_session']);
        
        % extract
        variable_A_data_all = variable_A_data_source{strcmp(variable_A_names_source, variable_A_name)};
        if size(variable_A_data_all, 1) == number_of_bands
            data_points_per_band_A = 1;
        elseif size(variable_A_data_all, 1) == number_of_time_steps_per_continuous_stretch
            data_points_per_band_A = number_of_time_steps_normalized;
        else
            error(['Variable ' variable_A_name ' has unexpected size'])
        end
        variable_B_data_all = variable_B_data_source{strcmp(variable_B_names_source, variable_B_name)};
        if size(variable_B_data_all, 1) == number_of_bands
            data_points_per_band_B = 1;
        elseif size(variable_B_data_all, 1) == number_of_time_steps_per_continuous_stretch
            data_points_per_band_B = number_of_time_steps_normalized;
        else
            error(['Variable ' variable_B_name ' has unexpected size'])
        end
        
        % make containers for global results
        results_this_variable = struct;
        results_this_variable.predictor = variable_A_name;
        results_this_variable.outcome = variable_B_name;
        results_this_variable.predictor_type = variable_A_type;
        results_this_variable.outcome_type = variable_B_type;
        results_this_variable.data.predictor = cell(number_of_bands, number_of_condition_combinations);;
        results_this_variable.data.outcome = cell(number_of_bands, number_of_condition_combinations);;
        results_this_variable.R_square = cell(number_of_bands, number_of_condition_combinations);
        results_this_variable.slope = cell(number_of_bands, number_of_condition_combinations);
        results_this_variable.slope_confidence_interval_width = cell(number_of_bands, number_of_condition_combinations);
        results_this_variable.correlation_p = cell(number_of_bands, number_of_condition_combinations);
        results_this_variable.row_info = study_settings.get('band_labels');
        results_this_variable.column_info = 'conditions, see the corresponding column in the conditions field for details';
        results_this_variable.conditions = condition_combinations_unique';
        results_this_variable.condition_labels = condition_labels';
        
        % loop through conditions
        for i_condition = 1 : number_of_condition_combinations
            % extract data from this condition
            this_condition_labels = condition_combinations_unique(i_condition, :);
            this_condition_indicator = condition_indicators(:, i_condition);
%             variable_A_data_this_condition = variable_A_data_all(:, this_condition_indicator);
%             variable_B_data_this_condition = variable_B_data_all(:, this_condition_indicator);
            
            % loop through bands
            for i_band = 1 : data.bands_per_stretch
                % get data for this band
                if data_points_per_band_A == 1
                    this_band_indices_A = i_band;
                else
                    [start_index, end_index] = getBandIndices(i_band, number_of_time_steps_normalized);
                    this_band_indices_A = start_index : end_index;
                end
                if data_points_per_band_B == 1
                    this_band_indices_B = i_band;
                else
                    this_band_indices_B = getBandIndices(i_band, number_of_time_steps_normalized);
                end
                variable_A_data_this_band = variable_A_data_all(this_band_indices_A, this_condition_indicator);
                variable_B_data_this_band = variable_B_data_all(this_band_indices_B, this_condition_indicator);
                
                % create containers
                R_square_table_here = zeros(data_points_per_band_A, data_points_per_band_B) * NaN;
                slope_table_here = zeros(data_points_per_band_A, data_points_per_band_B) * NaN;
                slope_confidence_interval_width_table_here = zeros(data_points_per_band_A, data_points_per_band_B) * NaN;
                offset_table_here = zeros(data_points_per_band_A, data_points_per_band_B) * NaN;
                offset_confidence_interval_width_table_here = zeros(data_points_per_band_A, data_points_per_band_B) * NaN;
                correlation_p_table_here = zeros(data_points_per_band_A, data_points_per_band_B) * NaN;
                
                results_this_variable.data.predictor{i_band, i_condition} = variable_A_data_this_band;
                results_this_variable.data.outcome{i_band, i_condition} = variable_B_data_this_band;
                
                for i_A = 1 : size(variable_A_data_this_band, 1)
                    for i_B = 1 : size(variable_B_data_this_band, 1)
                        predictor = variable_A_data_this_band(i_A, :)';
                        outcome = variable_B_data_this_band(i_B, :)';
                        if any(isnan(predictor)) || any(isnan(outcome))
                            
                        else
                            % fit
                            [fit_object, fit_stats] = fit(predictor, outcome, 'poly1');
                            [~, correlation_p] = corrcoef(predictor, outcome);
                            
                            % store
                            R_square_table_here(i_A, i_B) = fit_stats.rsquare;
                            slope_table_here(i_A, i_B) = fit_object.p1;
                            offset_table_here(i_A, i_B) = fit_object.p2;
                            confidence_intervals = confint(fit_object);
                            slope_confidence_interval_width_table_here(i_A, i_B) = range(confidence_intervals(:, 1));
                            offset_confidence_interval_width_table_here(i_A, i_B) = range(confidence_intervals(:, 2));
                            correlation_p_table_here(i_A, i_B) = correlation_p(1, 2);
                        end                        
                    end
                end
                
                % store
                results_this_variable.R_square{i_band, i_condition} = R_square_table_here;
                results_this_variable.slope{i_band, i_condition} = slope_table_here;
                results_this_variable.slope_confidence_interval_width{i_band, i_condition} = slope_confidence_interval_width_table_here;
                results_this_variable.offset{i_band, i_condition} = offset_table_here;
                results_this_variable.offset_confidence_interval_width{i_band, i_condition} = offset_confidence_interval_width_table_here;
                results_this_variable.correlation_p{i_band, i_condition} = correlation_p_table_here;

%                 % plot
%                 figure; 
%                 plot(fit_object, variable_A_data_this_condition(i_band, :)', variable_B_data_this_condition(i_band, :)')
%                 relevant_condition_label = strrep(this_condition_labels{2}, '_', ' '); % TODO: hard-coded for now, solve this properly later
%                 title_string = [relevant_condition_label ' - ' band_labels{i_band} ' - R^2 = ' num2str(R_square) ', slope = ' num2str(fit_object.p1) ', p = ' num2str(P(1, 2))];
%                 title(title_string);
% 
%                 xlabel(strrep(variable_A_name, '_', ' '));
%                 ylabel(strrep(variable_B_name, '_', ' '));
%                 
%                 % calculate variance and standard deviation
%                 std_predictor = std(predictor);
%                 std_outcome = std(outcome);
%                 std_residual = std(residual);
%                 
%                 var_predictor = var(predictor);
%                 var_outcome = var(outcome);
%                 var_residual = var(residual);
%                 
%                 disp([relevant_condition_label ' - ' variable_B_name ' STD: ' num2str(std_outcome)])
%                 disp([relevant_condition_label ' - ' variable_B_name ' Variance: ' num2str(var_outcome)])
            end
            
        end
        linear_model_results{i_variable, 1} = results_this_variable;
        linear_model_results{i_variable, 2} = variable_A_name;
        linear_model_results{i_variable, 3} = variable_A_type;
        linear_model_results{i_variable, 4} = variable_B_name;
        linear_model_results{i_variable, 5} = variable_B_type;
        disp(['Finished model ' num2str(i_variable) ' of ' num2str(number_of_variables) ', predictor: ' variable_A_name ', outcome: ' variable_B_name])
    end
    
    % save
    results_file_name = ['results' filesep collection_date '_' subject_id '_linearModelsFromSameBand.mat'];
    save(results_file_name, 'linear_model_results', 'linear_model_results_header');
end









