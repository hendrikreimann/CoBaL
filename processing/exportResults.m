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
variables_to_export_trajectories = study_settings.get('variables_to_export_trajectories');
band_labels = study_settings.get('band_labels');
number_of_variables_to_export = size(variables_to_export, 1);
number_of_variables_to_export_means = size(variables_to_export_means, 1);
number_of_variables_to_export_trajectories = size(variables_to_export_trajectories, 1);
number_of_data_points = size(time_list, 1);
number_of_bands = size(step_time_data, 1);
number_of_time_points_normalized = study_settings.get('number_of_time_steps_normalized');

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

% resample normalized time to get equidistant samples
step_time_means = mean(step_time_data, 2);
time_normalized = createScaledAbscissa(step_time_means, number_of_time_points_normalized);
time_rescaled = linspace(time_normalized(1), time_normalized(end), length(time_normalized))';

% gather and resample trajectory data
trajectory_data_cell_wide = cell(number_of_variables_to_export_means, 1);
for i_variable = 1 : number_of_variables_to_export_means
    % load and assemble data
    this_variable_name = variables_to_export_means{i_variable};
    variable_index_in_data_cell = find(strcmp(this_variable_name, variable_names));
    this_variable_data = variable_data{variable_index_in_data_cell}';
    this_variable_data_resampled = spline(time_normalized, this_variable_data, time_rescaled);
    trajectory_data_cell_wide{i_variable} = this_variable_data_resampled;
end
number_of_time_points_per_stretch = size(trajectory_data_cell_wide{1}, 2);

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
        trajectory_data_cell_wide{i_variable}(rows_to_remove, :) = [];
    end
    
    % remove from condition and origin data
    condition_cell(rows_to_remove, :) = [];
    origin_cell(rows_to_remove, :) = [];
end


% calculate means across repetitions for each condition
[unique_condition_combination_labels, unique_condition_combination_indicators] = getUniqueConditionInformation(condition_cell, condition_header);
number_of_condition_combinations = size(unique_condition_combination_labels, 1);
mean_data_cell_wide = cell(number_of_variables_to_export_means, 1);
for i_variable = 1 : number_of_variables_to_export_means
    this_variable_mean_data = zeros(number_of_time_points_per_stretch, number_of_condition_combinations);
    
    source_data_this_variable = trajectory_data_cell_wide{i_variable};
    for i_combination = 1 : number_of_condition_combinations
        indicator_this_combination = unique_condition_combination_indicators(:, i_combination);
        data_this_variable_this_combination = source_data_this_variable(indicator_this_combination, :);
        this_variable_mean_data(:, i_combination) = mean(data_this_variable_this_combination);
    end
    mean_data_cell_wide{i_variable} = this_variable_mean_data;
end

% repackage means for export
mean_header_cell = [condition_header 'time' variables_to_export_means'];
time_point_strings = num2str(time_rescaled);
time_point_cell_single = strtrim(cellstr(time_point_strings));
time_point_cell_long = repmat(time_point_cell_single, number_of_condition_combinations, 1);
condition_cell_long = cell(number_of_condition_combinations * number_of_time_points_per_stretch, 4);
for i_combination = 1 : number_of_condition_combinations
    this_combination = unique_condition_combination_labels(i_combination, :);
    condition_cell_long((i_combination-1)*number_of_time_points_per_stretch+1 : i_combination*number_of_time_points_per_stretch, :) ...
        = repmat(this_combination, number_of_time_points_per_stretch, 1);
end

mean_data_cell = cell(number_of_condition_combinations * number_of_time_points_per_stretch, number_of_variables_to_export_means);
for i_variable = 1 : number_of_variables_to_export_means
    % get data in shape
    this_variable_mean_data = mean_data_cell_wide{i_variable};
    this_variable_mean_data_flat = reshape(this_variable_mean_data, numel(this_variable_mean_data), 1);
    this_variable_mean_data_strings = num2str(this_variable_mean_data_flat);
    this_variable_mean_data_cell = strtrim(cellstr(this_variable_mean_data_strings));
    mean_data_cell(:, i_variable) = this_variable_mean_data_cell;
end
mean_body_cell = [condition_cell_long time_point_cell_long mean_data_cell];



% repackage trajectories for export - separate export for each subject
unique_subjects = table2cell(unique(cell2table(condition_cell(:, strcmp(condition_header, 'subject'))), 'rows'));
trajectory_header_cell = [condition_header 'time' variables_to_export_means'];
for i_subject = 1 : length(unique_subjects)
    % get indicators for this subject's data
    this_subject_label = unique_subjects{i_subject};
    this_subject_indicator = strcmp(condition_cell(:, strcmp(condition_header, 'subject')), this_subject_label);
    
    % create long containers for this subject's data
    number_of_stretches_this_subject = sum(this_subject_indicator);
    trajectory_data_cell_this_subject = cell(number_of_stretches_this_subject * number_of_time_points_per_stretch, number_of_variables_to_export_trajectories);

    % redistribute this subject's data into long containers, stretch by stretch
    time_point_cell_this_subject = repmat(time_point_cell_single, number_of_stretches_this_subject, 1);
    condition_cell_this_subject = condition_cell(this_subject_indicator, :);
    condition_cell_this_subject_long = cell(length(time_point_cell_this_subject), size(condition_cell_this_subject, 2));
    for i_stretch = 1 : number_of_stretches_this_subject
        this_combination = condition_cell_this_subject(i_stretch, :);
        condition_cell_this_subject_long((i_stretch-1)*number_of_time_points_per_stretch+1 : i_stretch*number_of_time_points_per_stretch, :) ...
            = repmat(this_combination, number_of_time_points_per_stretch, 1);
    end
    for i_variable = 1 : number_of_variables_to_export_trajectories
        % get data in shape
        this_variable_trajectory_data = trajectory_data_cell_wide{i_variable};
        this_variable_this_subject_trajectory_data = this_variable_trajectory_data(this_subject_indicator, :);
        
        this_variable_this_subject_trajectory_data_flat = reshape(this_variable_this_subject_trajectory_data', numel(this_variable_this_subject_trajectory_data), 1);
        this_variable_this_subject_trajectory_data_strings = num2str(this_variable_this_subject_trajectory_data_flat);
        this_variable_this_subject_trajectory_data_cell = strtrim(cellstr(this_variable_this_subject_trajectory_data_strings));
        trajectory_data_cell_this_subject(:, i_variable) = this_variable_this_subject_trajectory_data_cell;
    end
    trajectory_body_cell = [condition_cell_this_subject_long time_point_cell_this_subject trajectory_data_cell_this_subject];
    
    % save as .csv
    export_cell = ...
      [ ...
        trajectory_header_cell; ...
        trajectory_body_cell ...
      ];
    save_file_name = ['results_trajectories_' this_subject_label '.csv'];
    cell2csv(save_file_name, export_cell);
    disp(['Exported trajectory data for subject ' this_subject_label]);
end

% save as .csv
export_cell = ...
  [ ...
    header_cell; ...
    body_cell ...
  ];
cell2csv('results.csv', export_cell);

export_cell = ...
  [ ...
    mean_header_cell; ...
    mean_body_cell ...
  ];
save_file_name = 'results_trajectory_means.csv';
cell2csv(save_file_name, export_cell);




