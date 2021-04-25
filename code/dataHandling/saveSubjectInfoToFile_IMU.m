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

function saveSubjectInfoToFile_IMU(varargin)
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
    trials_to_exclude = 0;
  
    
  

    
    
    % save
    variables_to_save.condition_list = trial_type_list;
    variables_to_save.trial_number_list = trial_number_list;
    variables_to_save.trial_types_to_ignore = trial_types_to_ignore; %#ok<STRNU>
    variables_to_save.subject_ID = subject_id;
    variables_to_save.collection_date = collection_date;
    save_file_name = 'subjectInfo.mat';
    save(save_file_name, '-struct', 'variables_to_save');
    
end