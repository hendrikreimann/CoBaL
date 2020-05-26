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


% TODO: turn this into a function. For now, hard-code an input
% source_label = 'ucmAcrossEvents';
source_label = 'varianceAcrossEvents';
% source_label = 'displacement';

%% load data
loaded_data = load(['results_' source_label '.mat']);

%% load settings
if ~exist('studySettings.txt', 'file')
    error('No studySettings.txt file found. This function should be run from a study folder')
end    
study_settings = SettingsCustodian('studySettings.txt');
variables_to_export_discrete = study_settings.get('variables_to_export_discrete');
variables_to_export_discrete_header = study_settings.get('variables_to_export_discrete_header');

% remove variables from this list that don't match the source type
source_column = strcmp(variables_to_export_discrete_header, 'source file');
source_match = strcmp(variables_to_export_discrete(:, source_column), source_label);
variables_to_export_discrete(~source_match, :) = [];

number_of_variables_to_export_discrete = size(variables_to_export_discrete, 1);
number_of_data_points = size(loaded_data.time_list, 1);
number_of_bands = size(loaded_data.step_time_data, 1);

%% gather univariate data for each band
data_cell = {};
data_header = {};
for i_variable = 1 : number_of_variables_to_export_discrete
    % load and assemble data
    this_variable_name = variables_to_export_discrete{i_variable, strcmp(variables_to_export_discrete_header, 'export_variable_name')};
    this_variable_source = variables_to_export_discrete{i_variable, strcmp(variables_to_export_discrete_header, 'source_variable')};
    this_variable_band = str2double(variables_to_export_discrete{i_variable, strcmp(variables_to_export_discrete_header, 'source_band')});
    variable_index_in_data_cell = find(strcmp(this_variable_source, loaded_data.variable_names));
    this_variable_data_all_bands = loaded_data.variable_data{variable_index_in_data_cell}';
    this_variable_data = this_variable_data_all_bands(:, this_variable_band);
    
    % transform to strings
    data_strings = num2str(this_variable_data);
    this_variable_data_cell = strtrim(cellstr(data_strings));

    % store in appropriate location in data cell and header
    data_cell = [data_cell this_variable_data_cell];
    data_header = [data_header this_variable_name];
end

%% gather condition cell
conditions_settings = study_settings.get('conditions');
condition_labels = conditions_settings(:, 1)';
condition_source_variables = conditions_settings(:, 2)';
number_of_conditions = length(condition_labels);
condition_cell = cell(number_of_data_points, number_of_conditions);
condition_header = cell(1, number_of_conditions);
for i_condition = 1 : number_of_conditions
    this_condition_data = loaded_data.conditions.(condition_source_variables{i_condition});
    condition_cell(:, i_condition) = this_condition_data;
    condition_header(:, i_condition) = condition_labels(i_condition);
end

%% gather origin cell
if study_settings.get('export_origin_information', 1)
    origin_trial_number_data_cell = strtrim(cellstr(num2str(loaded_data.origin_trial_number_data)));
    origin_stretch_start_time_data_cell = strtrim(cellstr(num2str(loaded_data.origin_stretch_start_time_data)));
    origin_stretch_end_time_data_cell = strtrim(cellstr(num2str(loaded_data.origin_stretch_end_time_data)));
    origin_cell = [loaded_data.origin_session_folder_data origin_trial_number_data_cell origin_stretch_start_time_data_cell origin_stretch_end_time_data_cell];
    origin_header = {'origin folder', 'origin trial number', 'stretch start time within trial', 'stretch end time within trial'};
else
    origin_cell = {};
    origin_header = [];
end

%% join cells
header_cell = [condition_header, origin_header, data_header];
body_cell = [condition_cell, origin_cell, data_cell];

%% remove levels
levels_to_remove = study_settings.get('levels_to_remove_for_export', 1);
for i_level = 1 : size(levels_to_remove, 1)
    % determine rows to remove
    this_condition_label = levels_to_remove{i_level, 1};
    this_level_label = levels_to_remove{i_level, 2};
    relevant_column = strcmp(header_cell, this_condition_label);
    rows_to_remove = strcmp(body_cell(:, relevant_column), this_level_label);
    
    % remove
    body_cell(rows_to_remove, :) = [];
end

%% save
save_file = ['results_' source_label '.csv'];
export_cell = ...
  [ ...
    header_cell; ...
    body_cell ...
  ];
cell2csv(save_file, export_cell);



