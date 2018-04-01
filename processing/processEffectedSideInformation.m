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
    load('subjectInfo.mat', 'date', 'subject_id', 'most_affected');
    
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
    
    condition_effectedSide_list_session = cell(number_of_stretch_variables, 1);
    for i_stretch = 1 : number_of_stretches
        if most_affected == 'L'
            if strcmp(loaded_data.conditions_session.condition_startfoot_list{i_stretch}, 'STANCE_LEFT')
                condition_effectedSide_list_session{i_stretch} = 'Effected'
            else
                condition_effectedSide_list_session{i_stretch} = 'Non-Effected'
            end
        elseif most_affected == 'R'
            if strcmp(loaded_data.conditions_session.condition_startfoot_list{i_stretch}, 'STANCE_RIGHT')
                condition_effectedSide_list_session{i_stretch} = 'Effected'
            else
                condition_effectedSide_list_session{i_stretch} = 'Non-Effected'
            end
        end
        
    end
        % save
    variables_to_save = loaded_data;
    variables_to_save.conditions_session.condition_effectedSide_list = condition_effectedSide_list_session;
    save(results_file_name, '-struct', 'variables_to_save');
end