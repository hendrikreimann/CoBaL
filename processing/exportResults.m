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

% function exportResults
    
load results.mat

% load settings
if ~exist('studySettings.txt', 'file')
    error('No studySettings.txt file found. This function should be run from a study folder')
end    
study_settings = SettingsCustodian('studySettings.txt');
variables_to_export = study_settings.get('variables_to_export');
band_labels = study_settings.get('band_labels');
number_of_variables_to_export = size(variables_to_export, 1);
number_of_data_points = size(time_list, 1);
number_of_bands = length(band_labels);

% make data export cell
data_cell = cell(number_of_data_points, number_of_variables_to_export * number_of_bands);
data_header = cell(1, number_of_variables_to_export * number_of_bands);
for i_variable = 1 : number_of_variables_to_export
    % load and assemble data
    this_variable_name = variables_to_export{i_variable};
    variable_index_in_data_cell = find(strcmp(this_variable_name, variable_names));
    this_variable_data = variable_data{variable_index_in_data_cell}';
    
    % export
    for i_band = 1 : number_of_bands
        % extract band data and transform to strings
        data_strings = num2str(this_variable_data(:, i_band));
        this_variable_data_cell = strtrim(cellstr(data_strings));
        
        % store in appropriate location in data cell
        this_band_column = (i_variable-1)*3 + i_band;
        data_cell(:, this_band_column) = this_variable_data_cell;
        
        % make and store label
        data_header{this_band_column} = [this_variable_name '_' band_labels{i_band}];
    end
end
% time_data_strings = num2str(time_list);
% time_data_cell = strtrim(cellstr(time_data_strings));

% gather condition cell
conditions_settings = study_settings.get('conditions');
condition_labels = conditions_settings(:, 1)';
condition_source_variables = conditions_settings(:, 2)';
number_of_conditions = length(condition_labels);
condition_cell = cell(number_of_data_points, number_of_conditions);
condition_header = cell(1, number_of_conditions);
for i_condition = 1 : number_of_conditions
    this_condition_data = conditions.(condition_source_variables{i_condition});
    condition_cell(:, i_condition) = this_condition_data;
    condition_header(:, i_condition) = condition_labels(i_condition);
end

% gather origin cell
origin_trial_number_data_cell = strtrim(cellstr(num2str(origin_trial_number_data)));
origin_stretch_start_time_data_cell = strtrim(cellstr(num2str(origin_stretch_start_time_data)));
origin_stretch_end_time_data_cell = strtrim(cellstr(num2str(origin_stretch_end_time_data)));
origin_cell = [origin_session_folder_data origin_trial_number_data_cell origin_stretch_start_time_data_cell origin_stretch_end_time_data_cell];
origin_header = {'origin folder', 'origin trial number', 'stretch start time within trial', 'stretch end time within trial'};

% join to export cell
export_cell = ...
  [ ...
    condition_header, origin_header, data_header; ...
    condition_cell, origin_cell, data_cell ...
  ];

% save as .csv
% export_cell = [header; export_cell];
cell2csv('results.csv', export_cell);




