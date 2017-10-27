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


variables_to_export = {'step_placement_x', 'stimulus_response_x'};
number_of_variables_to_export = length(variables_to_export);
number_of_data_points = length(subject_list);

% find entries for left stim, to negate
left_stim_data = strcmp(condition_perturbation_list, 'ILLUSION_LEFT');

% make data export cell
data_cell = cell(number_of_data_points, number_of_variables_to_export);
data_array = zeros(number_of_data_points, number_of_variables_to_export) * NaN;
for i_variable = 1 : number_of_variables_to_export
    variable_index = find(strcmp(variables_to_export{i_variable}, variable_names));
    this_variable_data = response_data{variable_index}';
    
    % negate left stim data
    this_variable_data(left_stim_data) = -this_variable_data(left_stim_data);
    
    % export
    this_variable_data_strings = num2str(this_variable_data);
    this_variable_data_cell = strtrim(cellstr(this_variable_data_strings));
    data_cell(:, i_variable) = this_variable_data_cell;
    data_array(:, i_variable) = this_variable_data;
end

% process time
time_category = cell(size(time_list));
time_category_borders = [0 600 1200 1800 2400];
% time_category_borders = [0 1200 2400];
for i_point = 1 : number_of_data_points
    if time_category_borders(1) < time_list(i_point) && time_list(i_point) < time_category_borders(2)
        time_category{i_point} = 'TIME_ONE';
    end
    if time_category_borders(2) < time_list(i_point) && time_list(i_point) < time_category_borders(3)
        time_category{i_point} = 'TIME_TWO';
    end
    if time_category_borders(3) < time_list(i_point) && time_list(i_point) < time_category_borders(4)
        time_category{i_point} = 'TIME_THREE';
    end
    if time_category_borders(4) < time_list(i_point) && time_list(i_point) < time_category_borders(5)
        time_category{i_point} = 'TIME_FOUR';
    end
end
time_data_strings = num2str(time_list);
time_data_cell = strtrim(cellstr(time_data_strings));


% create stimulus label - TOWARDS or AWAY from stance foot
condition_stimulus_list = cell(size(condition_perturbation_list));
for i_point = 1 : number_of_data_points
    if strcmp(condition_stance_foot_list{i_point}, 'STANCE_LEFT') && strcmp(condition_perturbation_list{i_point}, 'ILLUSION_LEFT')
        condition_stimulus_list{i_point} = 'TOWARDS';
    end
    if strcmp(condition_stance_foot_list{i_point}, 'STANCE_LEFT') && strcmp(condition_perturbation_list{i_point}, 'ILLUSION_RIGHT')
        condition_stimulus_list{i_point} = 'AWAY';
    end
    if strcmp(condition_stance_foot_list{i_point}, 'STANCE_RIGHT') && strcmp(condition_perturbation_list{i_point}, 'ILLUSION_LEFT')
        condition_stimulus_list{i_point} = 'AWAY';
    end
    if strcmp(condition_stance_foot_list{i_point}, 'STANCE_RIGHT') && strcmp(condition_perturbation_list{i_point}, 'ILLUSION_RIGHT')
        condition_stimulus_list{i_point} = 'TOWARDS';
    end
end

% gather export cell
export_cell = [subject_list, condition_stance_foot_list, condition_perturbation_list, condition_stimulus_list, condition_index_list, time_category, time_data_cell, data_cell];

% remove control steps
control_steps = strcmp(export_cell(:, 3), 'CONTROL');
export_cell(control_steps, :) = [];
data_array(control_steps, :) = [];

% remove steps 2-4
ONE_steps = strcmp(export_cell(:, 4), 'ONE');
not_ONE_steps = ~ONE_steps;
export_cell(not_ONE_steps, :) = [];
data_array(not_ONE_steps, :) = [];

% save as .csv
header = [{'subject', 'stance_foot', 'perturbation', 'stimulus', 'step_index', 'time_category', 'time'}, variables_to_export];
export_cell = [header; export_cell];
cell2csv('results.csv', export_cell);




