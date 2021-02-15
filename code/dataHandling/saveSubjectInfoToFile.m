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

function saveSubjectInfoToFile(varargin)
    variables_to_save = struct;

    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'screen_folder', 'analysis')
    parse(parser, varargin{:})
    screen_folder = parser.Results.screen_folder;
    
    % get subject settings
    subject_settings = SettingsCustodian('subjectSettings.txt');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');

    trial_types_to_ignore = subject_settings.get('trial_types_to_ignore', 1);
   
    % get parameters from settings file
    parameters_cell = subject_settings.get('subject_info', 1);
    for i_parameter = 1 : size(parameters_cell, 1)
        this_parameter_value = parameters_cell{i_parameter, 2};
        if all(isstrprop(this_parameter_value, 'digit'))
            this_parameter_value = str2double(this_parameter_value);
        end
        variables_to_save.(parameters_cell{i_parameter, 1}) = this_parameter_value;
    end

    % get parameters
    data_dir = dir([screen_folder filesep '*.mat']);
    clear file_name_list;
    [file_name_list{1:length(data_dir)}] = deal(data_dir.name);
%     sample_file_name = file_name_list{1};
%     [date, subject_id] = getFileParameters(sample_file_name);
%     variables_to_save.date = date;
%     variables_to_save.subject_id = subject_id;

    % get list of conditions
    trial_type_list = {};
    trial_number_list = {};
    for i_file = 1 : length(file_name_list)
        data_file_name = file_name_list{i_file};
        % is this a data file?
        if sum(data_file_name=='_') == 4
            [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(data_file_name); %#ok<ASGLU>

            if ~any(strcmp(trial_type_list, trial_type))
                trial_type_list = [trial_type_list; trial_type]; %#ok<AGROW>
                trial_number_list = [trial_number_list; trial_number]; %#ok<AGROW>
            else
                % add current trial to trial number list
                condition_index = find(strcmp(trial_type_list, trial_type), 1);
                trial_number_list{condition_index} = [trial_number_list{condition_index} trial_number]; %#ok<AGROW>
                % remove duplicates
                trial_number_list{condition_index} = unique(trial_number_list{condition_index}); %#ok<AGROW>
            end
        end
    end
    
    % remove trials listed in subjectSettings.txt
    trials_to_exclude = subject_settings.get('trials_to_exclude', true);
    for i_trial = 1 : size(trials_to_exclude, 1)
        condition_index = find(strcmp(trial_type_list, trials_to_exclude{i_trial, 1}));
        if ~isempty(condition_index)
            trial_number_list_this_condition = trial_number_list{condition_index};
            matches = ismember(trial_number_list_this_condition, str2num(trials_to_exclude{i_trial, 2})); %#ok<ST2NM>
            trial_number_list_this_condition_pruned = trial_number_list_this_condition(~matches);
            trial_number_list{condition_index} = trial_number_list_this_condition_pruned;
        end
        
    end
    
    % remove conditions with no trials left
    empty_types = cellfun(@isempty, trial_number_list);
    trial_type_list(empty_types) = [];
    trial_number_list(empty_types) = [];
    
    % if we have a conditions file, remove the trials not listed there
    conditions_file_name = [];
    if exist('conditions.csv', 'file')
        conditions_file_name = 'conditions.csv';
    end
    if exist(makeFileName(collection_date, subject_id, 'conditions.csv'), 'file')
        conditions_file_name = makeFileName(collection_date, subject_id, 'conditions.csv');
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
        
        walking_condition_index = find(strcmp(trial_type_list, 'walking'));
        if ~isempty(walking_condition_index)
            trial_number_list_walking = trial_number_list{walking_condition_index};
            matches = ismember(trial_number_list_walking, trials_from_condition_file_list);
            matching_indices = find(matches);
            trial_number_list_walking_pruned = trial_number_list_walking(matching_indices);
            trial_number_list{walking_condition_index} = trial_number_list_walking_pruned;
        end
        
    end
    
    % save
    variables_to_save.condition_list = trial_type_list;
    variables_to_save.trial_number_list = trial_number_list;
    variables_to_save.trial_types_to_ignore = trial_types_to_ignore; %#ok<STRNU>
    if subject_settings.get('most_affected_side', 1)
        variables_to_save.affected_side = subject_settings.get('most_affected_side');
    end
    save_file_name = 'subjectInfo.mat';
    save(save_file_name, '-struct', 'variables_to_save');
    
end