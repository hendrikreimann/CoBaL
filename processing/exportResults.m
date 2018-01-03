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
    
% load results.mat

% load settings
if ~exist('studySettings.txt', 'file')
    error('No studySettings.txt file found. This function should be run from a study folder')
end    
study_settings = SettingsCustodian('studySettings.txt');
variables_to_export = study_settings.get('variables_to_export');
number_of_variables_to_export = size(variables_to_export, 1);
number_of_data_points = length(subject_list);

% make data export cell
data_cell = cell(number_of_data_points, number_of_variables_to_export);
data_array = zeros(number_of_data_points, number_of_variables_to_export) * NaN;
for i_variable = 1 : number_of_variables_to_export
    % load and assemble data
    this_variable_name = variables_to_export{i_variable};
    variable_index = find(strcmp(this_variable_name, variable_names));
    this_variable_data = variable_data{variable_index}';
    
    % export
    this_variable_data_strings = num2str(this_variable_data);
    this_variable_data_cell = strtrim(cellstr(this_variable_data_strings));
    data_cell(:, i_variable) = this_variable_data_cell;
    data_array(:, i_variable) = this_variable_data;
end
time_data_strings = num2str(time_list);
time_data_cell = strtrim(cellstr(time_data_strings));

% gather export cell
header = {'subject', 'stance_foot', 'perturbation', 'stimulus', 'step_index', 'direction', 'time_category', 'time'};
export_cell = [subject_list, condition_stance_foot_list, condition_perturbation_list, condition_stimulus_list, condition_index_list, condition_direction_list, time_category, time_data_cell];
if exist('condition_group_list', 'var')
    header = [header, 'group'];
    export_cell = [export_cell, condition_group_list];
end
header = [header, variables_to_export(:, 1)'];
export_cell = [export_cell, data_cell];

% remove control steps
control_steps = strcmp(export_cell(:, 3), 'CONTROL');
export_cell(control_steps, :) = [];
data_array(control_steps, :) = [];

% remove steps 2-4
TWO_steps = strcmp(export_cell(:, strcmp(header, 'step_index')), 'TWO');
export_cell(TWO_steps, :) = [];
data_array(TWO_steps, :) = [];
THREE_steps = strcmp(export_cell(:, strcmp(header, 'step_index')), 'THREE');
export_cell(THREE_steps, :) = [];
data_array(THREE_steps, :) = [];
FOUR_steps = strcmp(export_cell(:, strcmp(header, 'step_index')), 'FOUR');
export_cell(FOUR_steps, :) = [];
data_array(FOUR_steps, :) = [];

% save as .csv
export_cell = [header; export_cell];
cell2csv('results.csv', export_cell);




