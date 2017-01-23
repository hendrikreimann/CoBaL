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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% compare the kinematic tree against the kinematic chain

% function saveSubjectInfoToFile(gender, height, weight)
function saveSubjectInfoToFile
    
    % get subject code
    current_path = pwd;
    path_split = strsplit(current_path, filesep);
    subject_code = path_split{end};

    % try to open subject list from root and extract subject data
    subject_data_file = ['..' filesep 'subjects.csv'];
    if ~exist(subject_data_file, 'file')
        fileID = fopen(subject_data_file, 'w');
        fprintf(fileID,'ID,gender,height,weight,knee width,ankle width,elbow width,EMG1,EMG2,EMG3,EMG4,EMG5,EMG6,EMG7,EMG8\n');
        fprintf(fileID,',,m,kg,m,m,m\n');
        fprintf(fileID,'XYZ,0,0,0,0,0,0,LGLUTMED,LDELTANT,LGASTROC,LPEROLNG,RGLUTMED,RDELTANT,RGASTROC,RPEROLNG\n');
        fclose(fileID);
        
        disp('Failed to load "subjects.csv", a sample file has been created. Please edit it with your subject information or copy the correct file.')
        return
    end
    
    format = '%s';
    fid = fopen(subject_data_file);

    header_string = fgetl(fid);
    unit_string = fgetl(fid);
    data_raw = textscan(fid, format);
    fclose(fid);

    % find header info
    header = strsplit(header_string, ',');
    id_col = find(strcmp(header, 'ID'));
    gender_col = find(strcmp(header, 'gender'));
    height_col = find(strcmp(header, 'height'));
    weight_col = find(strcmp(header, 'weight'));
    knee_width_col = find(strcmp(header, 'knee width'));
    ankle_width_col = find(strcmp(header, 'ankle width'));
    elbow_width_col = find(strcmp(header, 'elbow width'));
    
    % find line for this subject
    data_lines = data_raw{1};
    data_cell = {};
    for i_line = 1 : length(data_lines)
        line_split = strsplit(data_lines{i_line}, ',');
        try
            data_cell = [data_cell; line_split];
        end
    end
    subject_row = find(strcmp(data_cell(:, id_col), subject_code));
    
    % extract data
    gender = data_cell{subject_row, gender_col};
    height = str2num(data_cell{subject_row, height_col});
    weight = str2num(data_cell{subject_row, weight_col});
    knee_width = str2num(data_cell{subject_row, knee_width_col});
    ankle_width = str2num(data_cell{subject_row, ankle_width_col});
    elbow_width = str2num(data_cell{subject_row, elbow_width_col});
    
    % find entries mapping EMG headers to muscle codes
    emg_sensor_map = {};
    for i_column = 1 : length(header)
        if length(header{i_column}) >=3 && strcmp(header{i_column}(1:3), 'EMG')
            muscle_code = data_cell{subject_row, i_column};
            emg_sensor_map = [emg_sensor_map, {header{i_column}; muscle_code}];
        end
    end
    
    % get parameters
    data_dir = dir(['raw' filesep '*.mat']);
    clear file_name_list;
    [file_name_list{1:length(data_dir)}] = deal(data_dir.name);
    sample_file_name = file_name_list{1};
    [date, subject_id] = getFileParameters(sample_file_name);

    % get list of conditions
    condition_list = {};
    trial_number_list = {};
    for i_file = 1 : length(file_name_list)
        data_file_name = file_name_list{i_file};
        [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(data_file_name);
        
        if ~any(strcmp(condition_list, trial_type))
            condition_list = [condition_list; trial_type];
            trial_number_list = [trial_number_list; trial_number];
        else
            % add current trial to trial number list
            condition_index = find(strcmp(condition_list, trial_type), 1);
            trial_number_list{condition_index} = [trial_number_list{condition_index} trial_number];
            % remove duplicates
            trial_number_list{condition_index} = unique(trial_number_list{condition_index});
        end
    end
    
    % save to file
    save ...
      ( ...
        'subjectInfo', ...
        'height', ...
        'weight', ...
        'gender', ...
        'knee_width', ...
        'ankle_width', ...
        'elbow_width', ...
        'date', ...
        'subject_id', ...
        'emg_sensor_map', ...
        'condition_list', ...
        'trial_number_list' ...
      );
end