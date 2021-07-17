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

% this function loads the results file containing stretch and analysis variables. Each stretch is assumed to be a gait
% cycle, starting with either the left foot or the right foot, as specified in the first_stance_leg condition. For each 
% stretch, we calculate all variables specified in linearModelSettings.txt
function extractLinearModelVariables(varargin)
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    linear_model_settings = loadSettingsFromFile('linearModel');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');

    % load settings and existing results
    variable_table = linear_model_settings.getTable('variables', 1);
    results_file_name = ['results' filesep makeFileName(collection_date, subject_id, 'results')];
    data = load(results_file_name);
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');
    number_of_variables = size(variable_table, 1);
    
    % make condition data tables
    conditions.settings_table = study_settings.get('conditions');
    conditions.factor_labels = conditions.settings_table(:, 1)';
    conditions.source_variables = conditions.settings_table(:, 2)';
    conditions.number_of_factor_labels = length(conditions.factor_labels);
    conditions.conditions_session = data.conditions_session;

    % go through variables
    variable_data = cell(number_of_variables, 1);
    variable_names = cell(number_of_variables, 1);
    variable_directions = cell(number_of_variables, 2);
    for i_variable = 1 : number_of_variables
        % extract information from table
        variable_name = variable_table.variable_name{i_variable};
        source_variable_name = variable_table.source_variable_name{i_variable};
        source_variable_type = variable_table.source_variable_type{i_variable};
        source_band = str2double(variable_table.source_band{i_variable});
        
        % get data
        variable_data_source = data.([source_variable_type '_data_session']);
        variable_names_source = data.([source_variable_type '_names_session']);
        variable_directions_source = data.([source_variable_type '_directions_session']);
        variable_data_all = variable_data_source{strcmp(variable_names_source, source_variable_name)};
        
        % get data from specified band
        if size(variable_data_all, 1) == data.bands_per_stretch
            % this is a discrete variable, take the one data point at the specified band
            this_variable_data = variable_data_all(source_band, :);
        else
            % this is a continuous variable, take the data points corresponding to the specified band
            [band_start_index, band_end_index] = getBandIndices(source_band, number_of_time_steps_normalized);
            band_indices = band_start_index : band_end_index;
            this_variable_data = variable_data_all(band_indices, :);
        end
        
        % store
        variable_data{i_variable} = this_variable_data;
        variable_names{i_variable} = variable_name;
        variable_directions(i_variable, :) = variable_directions_source(strcmp(variable_names_source, source_variable_name), :);
    end
    
    % save
    data_to_save = struct;
    data_to_save.variable_data = variable_data;
    data_to_save.variable_names = variable_names;
    data_to_save.variable_directions = variable_directions;
    data_to_save.conditions = data.conditions_session;
    
    results_file_name = ['results' filesep makeFileName(collection_date, subject_id, 'linearModelVariables.mat')];
    save(results_file_name, '-struct', 'data_to_save');
end











