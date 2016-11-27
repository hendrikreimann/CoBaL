% function saveSubjectInfoToFile(gender, height, weight)
function saveSubjectInfoToFile
    
%     if nargin < 1
%         gender = 'unknown';
%     end
%     if nargin < 2
%         height = -1;
%     end
%     if nargin < 3
%         weight = -1;
%     end
    
    % get subject code
    current_path = pwd;
    path_split = strsplit(current_path, filesep);
    subject_code = path_split{end};

    % try to open subject list from root and extract subject data
    subject_data_file = ['..' filesep 'subjects.csv'];
    
    
    number_of_header_lines = 1;
    format = '%s';
    fid = fopen(subject_data_file);

    header_string = fgetl(fid);
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
        'condition_list', ...
        'trial_number_list' ...
      );
end