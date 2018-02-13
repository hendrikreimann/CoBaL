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



function processGroupInformation(varargin)
    load('subjectInfo.mat', 'date', 'subject_id');
    % load settings and existing results
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'results')];
    loaded_data = load(results_file_name);
    
    number_of_stretch_variables = length(loaded_data.stretch_names_session);
    number_of_stretches = size(loaded_data.stretch_data_session{1}, 2); %#ok<*USENS>
    
    % make group condition list
    group_assignment = study_settings.get('group_assignment');
    if isempty(group_assignment)
        error('Trying to assign groups, but no group_assignment information found in studySettings.txt')
    end
    condition_group_list_session = cell(number_of_stretch_variables, 1);
    for i_stretch = 1 : number_of_stretches
        groups_this_subject = group_assignment(strcmp(group_assignment(:, 1), subject_id), 2:3);
        if isempty(groups_this_subject)
            error(['Trying to assign groups, but subject "' subject_id '" not found in group_assignment in studySettings.txt'])
        end
        % TO DO: figure out a way to decipher for Early-Late comparison
        if strcmp(loaded_data.conditions_session.condition_startfoot_list{i_stretch}, 'STANCE_LEFT')
            condition_group_list_session{i_stretch} = groups_this_subject{1};
        elseif strcmp(loaded_data.conditions_session.condition_startfoot_list{i_stretch}, 'STANCE_RIGHT')
            condition_group_list_session{i_stretch} = groups_this_subject{2};
        else
            error('Can only assign a group for STANCE_LEFT or STANCE_RIGHT')
        end
    end

    % save
    variables_to_save = loaded_data;
    variables_to_save.conditions_session.condition_group_list = condition_group_list_session;
    save(results_file_name, '-struct', 'variables_to_save');
    
end