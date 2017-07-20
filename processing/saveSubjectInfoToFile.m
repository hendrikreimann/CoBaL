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

function saveSubjectInfoToFile
    
    % get subject code
    current_path = pwd;
    path_split = strsplit(current_path, filesep);

    % open subject list from root and extract subject data
    subject_data_file = '';
    if exist(['..' filesep 'subjects.csv'], 'file')
        subject_data_file = ['..' filesep 'subjects.csv'];
        subject_code = path_split{end};
    end    
    if exist(['..' filesep '..' filesep 'subjects.csv'], 'file')
        subject_data_file = ['..' filesep '..' filesep 'subjects.csv'];
        subject_code = path_split{end-1};
    end
    if isempty(subject_data_file)
        disp('Failed to load "subjects.csv".')
        return
    end
    
    format = '%s';
    fid = fopen(subject_data_file);

    header_string = fgetl(fid);
    unit_string = fgetl(fid); %#ok<NASGU>
    data_raw = textscan(fid, format);
    fclose(fid);

    variables_to_save = struct;
    
    
    % find header info
    header = strsplit(header_string, ',');
    
    % find line for this subject
    data_lines = data_raw{1};
    data_cell = {};
    for i_line = 1 : size(data_lines, 1)
        line_split = strsplit(data_lines{i_line}, ',');
        data_cell = [data_cell; line_split]; %#ok<AGROW>
    end
    subject_row = find(strcmp(data_cell(:, strcmp(header, 'ID')), subject_code));
    
    % extract data
%     gender = data_cell{subject_row, strcmp(header, 'gender')}; %#ok<NASGU>
%     height = str2num(data_cell{subject_row, strcmp(header, 'height')}); %#ok<ST2NM,NASGU>
%     weight = str2num(data_cell{subject_row, strcmp(header, 'weight')}); %#ok<ST2NM,NASGU>
%     knee_width = str2num(data_cell{subject_row, strcmp(header, 'knee width')}); %#ok<ST2NM,NASGU>
%     ankle_width = str2num(data_cell{subject_row, strcmp(header, 'ankle width')}); %#ok<ST2NM,NASGU>
%     elbow_width = str2num(data_cell{subject_row, strcmp(header, 'elbow width')}); %#ok<ST2NM,NASGU>
    
    for i_column = 2 : length(header)
        variable_name = strrep(header{i_column}, ' ', '_');
        variable_value = data_cell{subject_row, i_column};
        if all(ismember(variable_value, '0123456789-.'))
            variable_value = str2num(variable_value); %#ok<ST2NM,NASGU>
        end
        evalstring = ['variables_to_save.' variable_name ' = variable_value;'];
        eval(evalstring);
    end
    
    % find entries mapping EMG headers to muscle codes
    emg_sensor_map = {};
    for i_column = 1 : length(header)
        if length(header{i_column}) >=3 && strcmp(header{i_column}(1:3), 'EMG')
            muscle_code = data_cell{subject_row, i_column};
            emg_sensor_map = [emg_sensor_map, {header{i_column}; muscle_code}]; %#ok<AGROW>
        end
    end
    variables_to_save.emg_sensor_map = emg_sensor_map;
    
    % get parameters
    data_dir = dir(['raw' filesep '*.mat']);
    clear file_name_list;
    [file_name_list{1:length(data_dir)}] = deal(data_dir.name);
    sample_file_name = file_name_list{1};
    [date, subject_id] = getFileParameters(sample_file_name);
    variables_to_save.date = date;
    variables_to_save.subject_id = subject_id;

    % get list of conditions
    condition_list = {};
    trial_number_list = {};
    for i_file = 1 : length(file_name_list)
        data_file_name = file_name_list{i_file};
        [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(data_file_name); %#ok<ASGLU>
        
        if ~any(strcmp(condition_list, trial_type))
            condition_list = [condition_list; trial_type]; %#ok<AGROW>
            trial_number_list = [trial_number_list; trial_number]; %#ok<AGROW>
        else
            % add current trial to trial number list
            condition_index = find(strcmp(condition_list, trial_type), 1);
            trial_number_list{condition_index} = [trial_number_list{condition_index} trial_number]; %#ok<AGROW>
            % remove duplicates
            trial_number_list{condition_index} = unique(trial_number_list{condition_index}); %#ok<AGROW>
        end
    end
    variables_to_save.condition_list = condition_list;

    
    % if we have a conditions file, remove the trials not listed there
    conditions_file_name = [];
    if exist('conditions.csv', 'file')
        conditions_file_name = 'conditions.csv';
    end
    if exist(makeFileName(date, subject_id, 'conditions.csv'), 'file')
        conditions_file_name = makeFileName(date, subject_id, 'conditions.csv');
    end
    if ~isempty(conditions_file_name)
        % read file
        fileID = fopen(conditions_file_name, 'r');
        header_line = fgetl(fileID);
        text_cell = {};
        text_line = fgetl(fileID);
        while ischar(text_line)
            text_cell = [text_cell; text_line]; %#ok<AGROW>
            text_line = fgetl(fileID);
        end
        fclose(fileID);

        % transform to arrays
        header = strsplit(strrep(header_line, ' ', ''), ',');
        condition_header = header(2 : end);
        number_of_trials = size(text_cell, 1);
        number_of_conditions = size(header, 2) - 1;
        trials_from_condition_file_list = zeros(number_of_trials, 1) * NaN;
        condition_cell = cell(number_of_trials, number_of_conditions);
        for i_trial = 1 : number_of_trials
            text_line = text_cell{i_trial};
            line_split = strsplit(text_line, ',');
            trials_from_condition_file_list(i_trial) = str2num(line_split{1});
            condition_cell(i_trial, :) = line_split(2:end);
        end
        
        walking_condition_index = find(strcmp(condition_list, 'walking'));
        if ~isempty(walking_condition_index)
            trial_number_list_walking = trial_number_list{walking_condition_index};
            matches = ismember(trial_number_list_walking, trials_from_condition_file_list);
            matching_indices = find(matches);
            trial_number_list_walking_pruned = trial_number_list_walking(matching_indices);
            trial_number_list{walking_condition_index} = trial_number_list_walking_pruned;
        end
        
    end
    variables_to_save.trial_number_list = trial_number_list;
    
    % save
    save_file_name = 'subjectInfo.mat';
    save(save_file_name, '-struct', 'variables_to_save');
    
end