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

function importApdmCsv (varargin)
    source_dir = 'CSV';
    
    % load settings
    subject_settings = loadSettingsFromFile('subject');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');

    % load file containing list of walking files
    file_name = 'Walk_trials.csv';
    file_table = readtable([source_dir filesep file_name]);
    number_of_files = size(file_table, 1);
    gait_speed_averages = zeros(number_of_files, 1);
    condition_labels = cell(number_of_files, 1);
    for i_file = 1 : number_of_files
        this_file_time = file_table.RecordDate_Time{i_file};
        this_file_note = file_table.TrialNotes{i_file};
        
        % figure out trial type and number
        trial_number = 1;
        trial_type = 'TBD';
        if strcmp(this_file_note, 'Preferred walking speed')
            trial_type = 'Preferred';
        end
        if strcmp(this_file_note, 'Preferred walking speed + cognitive task')
            trial_type = 'CogTask';
        end
        if strcmp(this_file_note, 'Fast walking speed')
            trial_type = 'Fast';
        end
        
        % load data
        this_file_name = [this_file_time '_Walk_Trial.csv'];
        % this_file_name = ['Walk_trials.csv'];
        data_table = readtable([source_dir filesep this_file_name]);
        
        % extract data
        data_column_start = 6; % was 6, hardcoded, assuming this is fixed for all APDM csv files we'll try to import
        label_column = 1;
        labels = data_table{:, label_column};

        left_touchdown_times_1 = data_table{strcmp(labels, 'Gait - Lower Limb - Gait Cycle L (s)'), data_column_start : end}';
        right_touchdown_times_1 = data_table{strcmp(labels, 'Gait - Lower Limb - Gait Cycle R (s)'), data_column_start : end}';
        left_midswing_times = data_table{strcmp(labels, 'Gait - Lower Limb - Midswing L (s)'), data_column_start : end}';
        right_midswing_times = data_table{strcmp(labels, 'Gait - Lower Limb - Midswing R (s)'), data_column_start : end}';
        left_pushoff_times = data_table{strcmp(labels, 'Gait - Lower Limb - Toe Off L (s)'), data_column_start : end}';
        right_pushoff_times = data_table{strcmp(labels, 'Gait - Lower Limb - Toe Off R (s)'), data_column_start : end}';
        left_touchdown_times_2 = data_table{strcmp(labels, 'Gait - Lower Limb - Initial Contact Next L (s)'), data_column_start : end}';
        right_touchdown_times_2 = data_table{strcmp(labels, 'Gait - Lower Limb - Initial Contact Next R (s)'), data_column_start : end}';
        turn_start_times = data_table{strcmp(labels, 'Turns - Turn (s)'), data_column_start : end}';
        turn_start_times(isnan(turn_start_times)) = [];
        
        left_touchdown_times(isnan(left_touchdown_times)) = [];
        right_touchdown_times(isnan(right_touchdown_times)) = [];
        left_midswing_times(isnan(left_midswing_times)) = [];
        right_midswing_times(isnan(right_midswing_times)) = [];
        left_pushoff_times(isnan(left_pushoff_times)) = [];
        right_pushoff_times(isnan(right_pushoff_times)) = [];

        % merge touchdowns that start and end the gait cycles identified by APDM
        left_touchdown_times = sort(unique([left_touchdown_times_1; left_touchdown_times_2]));
        right_touchdown_times = sort(unique([right_touchdown_times_1; right_touchdown_times_2]));
        
        variables_to_save = struct;

        event_data = ...
          { ...
            left_pushoff_times; ...
            left_touchdown_times; ...
            left_midswing_times; ...
            right_pushoff_times; ...
            right_touchdown_times; ...
            right_midswing_times; ...
            turn_start_times; ...
          };
        event_labels = ...
          { ...
            'left_pushoff'; ...
            'left_touchdown'; ...
            'left_midswing'; ...
            'right_pushoff'; ...
            'right_touchdown'; ...
            'right_midswing'; ...
            'turn_start'; ...
          };

        % add new variables to be saved
        variables_to_save.event_data = event_data;
        variables_to_save.event_labels = event_labels;

        step_events_file_name = ['analysis' filesep makeFileName(collection_date, subject_id, trial_type, trial_number, 'events.mat')];
        saveDataToFile(step_events_file_name, variables_to_save);

        disp(['Finding Step Events: type ' trial_type ', Trial ' num2str(trial_number) ' imported, saved as ' step_events_file_name]);

        
        % extract mean gait speed
        mean_data_column = 4;
        gait_speed_average_left = data_table{strcmp(labels, 'Gait - Lower Limb - Gait Speed L (m/s)'), mean_data_column}';
        gait_speed_average_right = data_table{strcmp(labels, 'Gait - Lower Limb - Gait Speed R (m/s)'), mean_data_column}';
        gait_speed_average_this_file = (gait_speed_average_left + gait_speed_average_right) * 0.5;

        gait_speed_averages(i_file) = gait_speed_average_this_file;
        condition_labels{i_file} = trial_type;

    end
    
    % save gait speed results
    if ~directoryExists('results')
        mkdir results
    end
    filename = ['results' filesep makeFileName(collection_date, subject_id, 'gaitSpeed.mat')];
    save(filename, 'gait_speed_averages', 'condition_labels')
end