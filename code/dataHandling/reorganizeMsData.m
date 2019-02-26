%     This file is part of the QuietStanceUcm code base
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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% input


function reorganizeMsData(varargin)
    % set some parameters
    millimeter_to_meter = 1e-3;
    milliseconds_to_seconds = 1e-3;
    
    % set and prepare folders
    source_folder = 'original';
    if ~directoryExists('processed')
        mkdir('processed')
    end
    if ~directoryExists('analysis')
        mkdir('analysis')
    end
    
    % load settings
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    
    % get list of files to reorganize
    clear file_name_list;
    data_dir = dir([source_folder filesep '*.mat']);
    [file_name_list{1:length(data_dir)}] = deal(data_dir.name);
    number_of_files = length(file_name_list);
    
    % go through files and reorganize
    for i_file = 1 : number_of_files
        data_file_name = file_name_list{i_file};
        load([source_folder filesep data_file_name]);
        
        % general items
        date = reformatDate(Date);
        file_name_split = strsplit(data_file_name, '.');
        file_name_split = strsplit(file_name_split{1}, '_');
        trial_type = file_name_split{1};
        subject_id = file_name_split{2};
        trial_number = zeroPrefixedIntegerString(str2num(file_name_split{3}(2:end)), 3);
        
        % condition
%         addTrialToConditionFile(date, subject_id, str2num(trial_number));
        condition_label = determineConditionLabel(str2num(trial_number));
        
        % marker data
        number_of_markers = length(pos_labels);
        number_of_time_steps = length(pos_time);
        marker_trajectories = zeros(number_of_time_steps, number_of_markers*3);
        marker_trajectories(:, 1 : 3 : end) = pos_X * millimeter_to_meter;
        marker_trajectories(:, 2 : 3 : end) = pos_Y * millimeter_to_meter;
        marker_trajectories(:, 3 : 3 : end) = pos_Z * millimeter_to_meter;
        time_mocap = pos_time * milliseconds_to_seconds;
        sampling_rate_mocap = round(median(diff(time_mocap))^(-1));
        
        % triplicate labels
        number_of_markers = length(pos_labels);
        marker_labels = cell(3, number_of_markers);
        for i_marker = 1 : length(marker_labels)
            marker_labels{1, i_marker} = [pos_labels{i_marker} '_x'];
            marker_labels{2, i_marker} = [pos_labels{i_marker} '_y'];
            marker_labels{3, i_marker} = [pos_labels{i_marker} '_z'];
        end
        marker_labels = reshape(marker_labels, 1, number_of_markers*3);
        
        % make directions
        % NOTE: this defines directions and makes assumptions, make sure everything is right here
        number_of_marker_trajectories = size(marker_trajectories, 2);
        marker_directions = cell(2, number_of_marker_trajectories);
        [marker_directions{1, 1 : 3 : number_of_marker_trajectories}] = deal('forward');
        [marker_directions{2, 1 : 3 : number_of_marker_trajectories}] = deal('backward');
        [marker_directions{1, 2 : 3 : number_of_marker_trajectories}] = deal('left');
        [marker_directions{2, 2 : 3 : number_of_marker_trajectories}] = deal('right');
        [marker_directions{1, 3 : 3 : number_of_marker_trajectories}] = deal('up');
        [marker_directions{2, 3 : 3 : number_of_marker_trajectories}] = deal('down');
        
        
        % filter
        if study_settings.get('filter_marker_data')
            filter_order = study_settings.get('marker_data_filter_order');
            cutoff_frequency = study_settings.get('marker_data_cutoff_frequency'); % in Hz
            [b_marker, a_marker] = butter(filter_order, cutoff_frequency/(sampling_rate_mocap/2));
            marker_trajectories = nanfiltfilt(b_marker, a_marker, marker_trajectories);
        end

        % save
        save_folder = 'processed';
        save_file_name = makeFileName(date, subject_id, condition_label, trial_number, 'markerTrajectories.mat');
        save ...
          ( ...
            [save_folder filesep save_file_name], ...
            'marker_trajectories', ...
            'time_mocap', ...
            'sampling_rate_mocap', ...
            'marker_labels',  ...
            'marker_directions' ...
          );
        addAvailableData ...
          ( ...
            'marker_trajectories', ...
            'time_mocap', ...
            'sampling_rate_mocap', ...
            '_marker_labels', ...
            '_marker_directions', ...
            save_folder, ...
            save_file_name ...
          );
        
        
        
        disp(['Reorganizing data, trial type ' trial_type ', number ' trial_number ' completed, saved as ' save_folder filesep save_file_name]);
        
        

        
    end





end

function new_string = reformatDate(old_string)
    % extract
    parts = strsplit(old_string, '-');
    old_day_string = parts{1};
    old_month_string = parts{2};
    old_year_string = parts{3};

    % reformat
    new_day_string = old_day_string;
    if strcmp(old_month_string, 'Jan')
        new_month_string = '01';
    end
    if strcmp(old_month_string, 'Feb')
        new_month_string = '02';
    end
    if strcmp(old_month_string, 'Mar')
        new_month_string = '03';
    end
    if strcmp(old_month_string, 'Apr')
        new_month_string = '04';
    end
    if strcmp(old_month_string, 'May')
        new_month_string = '05';
    end
    if strcmp(old_month_string, 'Jun')
        new_month_string = '06';
    end
    if strcmp(old_month_string, 'Jul')
        new_month_string = '07';
    end
    if strcmp(old_month_string, 'Aug')
        new_month_string = '08';
    end
    if strcmp(old_month_string, 'Sep')
        new_month_string = '09';
    end
    if strcmp(old_month_string, 'Oct')
        new_month_string = '10';
    end
    if strcmp(old_month_string, 'Nov')
        new_month_string = '11';
    end
    if strcmp(old_month_string, 'Dec')
        new_month_string = '12';
    end
    new_year_string = ['20' old_year_string];
    
    new_string = [new_year_string new_month_string new_day_string];
end

function new_string = reformatCondition(old_string, trial_number)
    if strcmp(old_string(1:2), 'IL')
        block_1 = 31:35;
        block_2 = 36:40;
        block_3 = 41:45;
        block_4 = 46:50;
        block_5 = 51:55;
        if ismember(trial_number, block_1)
            new_string = ['continuous_B1'];
        end
        if ismember(trial_number, block_2)
            new_string = ['continuous_B2'];
        end
        if ismember(trial_number, block_3)
            new_string = ['continuous_B3'];
        end
        if ismember(trial_number, block_4)
            new_string = ['continuous_B4'];
        end
        if ismember(trial_number, block_5)
            new_string = ['continuous_B5'];
        end
    else
        new_string = old_string;
    end
end

function condition_label = determineConditionLabel(trial_number)
    % determine condition string
    if trial_number == 2
        condition_label = 'calibration';
    end
    if trial_number == 3
        condition_label = 'quietEO';
    end
    if trial_number == 4
        condition_label = 'quietEC';
    end
    if ismember(trial_number, 5:9)
        condition_label = 'ramp012';
    end
    if ismember(trial_number, 10:14)
        condition_label = 'ramp036';
    end
    if ismember(trial_number, 15:19)
        condition_label = 'ramp060';
    end
    if ismember(trial_number, 20:24)
        condition_label = 'ramp084';
    end
    if ismember(trial_number, 25:29)
        condition_label = 'ramp120';
    end

    if ismember(trial_number, 31:35)
        condition_label = 'continuous01';
    end
    if ismember(trial_number, 36:40)
        condition_label = 'continuous02';
    end
    if ismember(trial_number, 41:45)
        condition_label = 'continuous03';
    end
    if ismember(trial_number, 46:50)
        condition_label = 'continuous04';
    end
    if ismember(trial_number, 51:55)
        condition_label = 'continuous05';
    end
end

% function addTrialToConditionFile(date, subject_id, trial_number)
% 
%     
% 
% 
%     conditions_file_name = makeFileName(date, subject_id, 'conditions.csv');
%     
%     
%     
%     
%     condition = reformatCondition(condition_string, trial_number);
%     
%     if ~exist(conditions_file_name, 'file')
%         header_line = 'trial,condition\n';
%         file_id = fopen(conditions_file_name, 'w');
%         fprintf(file_id, header_line);
%         fclose(file_id);
%     end
%     
%     % load conditions file
%     file_id = fopen(conditions_file_name, 'r');
%     header_line = fgetl(file_id);
%     text_cell = {};
%     text_line = fgetl(file_id);
%     while ischar(text_line)
%         text_cell = [text_cell; text_line]; %#ok<AGROW>
%         text_line = fgetl(file_id);
%     end
%     fclose(file_id);
%     
%     % transform to arrays
%     header = strsplit(strrep(header_line, ' ', ''), ',');
%     condition_header = header(2 : end);
%     number_of_trials = size(text_cell, 1);
%     number_of_conditions = size(header, 2) - 1;
%     trials_from_condition_file_list = zeros(number_of_trials, 1) * NaN;
%     condition_cell = cell(number_of_trials, number_of_conditions);
%     for i_trial = 1 : number_of_trials
%         text_line = text_cell{i_trial};
%         line_split = strsplit(text_line, ',');
%         trials_from_condition_file_list(i_trial) = str2num(line_split{1});
%         condition_cell(i_trial, :) = line_split(2:end);
%     end
%     
%     % check if current trial is already there
%     condition_column = find(strcmp(condition_header, 'condition'));
%     condition_fit = strcmp(condition_cell(:, condition_column), condition);
%     trial_fit = (trials_from_condition_file_list == trial_number);
%     complete_fit = trial_fit & condition_fit;
%     if ~any(complete_fit)
%         % this trial isn't listed yet, add it
%         condition_line_cell = cell(1, number_of_conditions);
%         condition_line_cell{condition_column} = condition;
%         condition_line = num2str(trial_number);
%         for i_condition = 1 : number_of_conditions
%             condition_line = [condition_line ',' condition_line_cell{i_condition}];
%         end
%         
%         % write to file
%         header_line = [condition_line '\n'];
%         file_id = fopen(conditions_file_name, 'a');
%         fprintf(file_id, header_line);
%         fclose(file_id);
%     end
%     
%     
%     
% end









