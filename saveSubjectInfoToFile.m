function saveSubjectInfoToFile(height, weight, gender)
    
    if nargin < 1
        height = -1;
    end
    if nargin < 2
        weight = -1;
    end
    if nargin < 3
        gender = 'unknown';
    end
    
    data_dir = dir('*_markerTrajectories.mat');
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
            trial_number_list{condition_index} = [trial_number_list{condition_index}; trial_number];
        end
    end

    
    
    save ...
      ( ...
        'subjectInfo', ...
        'height', ...
        'weight', ...
        'gender', ...
        'date', ...
        'subject_id', ...
        'condition_list', ...
        'trial_number_list' ...
      );
end