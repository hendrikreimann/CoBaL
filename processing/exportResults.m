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
variables_to_export_means = study_settings.get('variables_to_export_means');
band_labels = study_settings.get('band_labels');
number_of_variables_to_export = size(variables_to_export, 1);
number_of_variables_to_export_means = size(variables_to_export_means, 1);
number_of_data_points = size(time_list, 1);
number_of_bands = length(band_labels);

% export univariate data for each band
data_cell = cell(number_of_data_points, number_of_variables_to_export * number_of_bands);
data_header = cell(1, number_of_variables_to_export * number_of_bands);
for i_variable = 1 : number_of_variables_to_export
    % load and assemble data
    this_variable_name = variables_to_export{i_variable};
    variable_index_in_data_cell = find(strcmp(this_variable_name, variable_names));
    this_variable_data = variable_data{variable_index_in_data_cell}';
    
    % export data
    for i_band = 1 : number_of_bands
        % extract band data and transform to strings
        data_strings = num2str(this_variable_data(:, i_band));
        this_variable_data_cell = strtrim(cellstr(data_strings));
        
        % store in appropriate location in data cell
        this_band_column = (i_variable-1)*number_of_bands + i_band;
        data_cell(:, this_band_column) = this_variable_data_cell;
        
        % make and store label
        if number_of_bands > 1
            data_header{this_band_column} = [this_variable_name '_' band_labels{i_band}];
        else
            data_header{this_band_column} = this_variable_name;
        end
    end
end

% gather trajectory data
source_data_cell_for_means = cell(number_of_variables_to_export_means, 1);
for i_variable = 1 : number_of_variables_to_export_means
    % load and assemble data
    this_variable_name = variables_to_export_means{i_variable};
    variable_index_in_data_cell = find(strcmp(this_variable_name, variable_names));
    source_data_cell_for_means{i_variable} = variable_data{variable_index_in_data_cell}';
end
number_of_time_points_per_stretch = size(source_data_cell_for_means{1}, 2);

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

% join cells
header_cell = [condition_header, origin_header, data_header];
body_cell = [condition_cell, origin_cell, data_cell];

% remove levels
levels_to_remove = study_settings.get('levels_to_remove_for_export');
for i_level = 1 : size(levels_to_remove, 1)
    this_condition_label = levels_to_remove{i_level, 1};
    this_level_label = levels_to_remove{i_level, 2};
    relevant_column = strcmp(header_cell, this_condition_label);
    rows_to_remove = strcmp(body_cell(:, relevant_column), this_level_label);
    
    % remove from univariate data
    body_cell(rows_to_remove, :) = [];
    
    % remove from trajectory data
    for i_variable = 1 : number_of_variables_to_export_means
        source_data_cell_for_means{i_variable}(rows_to_remove, :) = [];
    end
    
    % remove from condition and origin data
    condition_cell(rows_to_remove, :) = [];
    origin_cell(rows_to_remove, :) = [];
end



% calculate means across repetitions for each condition
[unique_condition_combination_labels, unique_condition_combination_indicators] = getUniqueConditionInformation(condition_cell, condition_header);
number_of_condition_combinations = size(unique_condition_combination_labels, 1);
mean_data_cell = cell(number_of_variables_to_export_means, 1);
for i_variable = 1 : number_of_variables_to_export_means
    this_variable_mean_data = zeros(number_of_time_points_per_stretch, number_of_condition_combinations);
    
    source_data_this_variable = source_data_cell_for_means{i_variable};
    for i_combination = 1 : number_of_condition_combinations
        indicator_this_combination = unique_condition_combination_indicators(:, i_combination);
        data_this_variable_this_combination = source_data_this_variable(indicator_this_combination, :);
        this_variable_mean_data(:, i_combination) = mean(data_this_variable_this_combination);
    end
    mean_data_cell{i_variable} = this_variable_mean_data;
end

% repackage means for export
time_point_strings = num2str((1 : number_of_time_points_per_stretch)');
time_point_cell = strtrim(cellstr(time_point_strings));
mean_header_cell = [condition_header time_point_cell'];
for i_variable = 1 : number_of_variables_to_export_means
    % get data in shape
    this_variable_mean_data = mean_data_cell{i_variable};
    this_variable_mean_data_flat = reshape(this_variable_mean_data, numel(this_variable_mean_data), 1);
    this_variable_mean_data_strings = num2str(this_variable_mean_data_flat);
    this_variable_mean_data_cell_flat = cellstr(this_variable_mean_data_strings);
    this_variable_mean_data_cell = strtrim(reshape(this_variable_mean_data_cell_flat, number_of_time_points_per_stretch, size(this_variable_mean_data, 2)))';
    
    mean_body_cell = [unique_condition_combination_labels this_variable_mean_data_cell];
    export_cell = ...
      [ ...
        mean_header_cell; ...
        mean_body_cell ...
      ];
    save_file_name = ['results_' variables_to_export_means{i_variable} '.csv'];
    cell2csv(save_file_name, export_cell);
end




% save as .csv
export_cell = ...
  [ ...
    header_cell; ...
    body_cell ...
  ];
cell2csv('results.csv', export_cell);





