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



function addGroupsWithinPopulation(varargin)
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'subjects', [])
    parse(parser, varargin{:})
    subjects = parser.Results.subjects;

 % load settings
    if ~exist('studySettings.txt', 'file')
        error('No studySettings.txt file found. This function should be run from a study folder')
    end    
    study_settings = SettingsCustodian('studySettings_TF.txt');
    
    data_folder_list = determineDataStructure(subjects);
    
    condition_stance_foot_list = {};
    condition_perturbation_list = {};
    condition_delay_list = {};
    condition_index_list = {};
    condition_experimental_list = {};
    condition_stimulus_list = {};
    condition_day_list = {};
    subject_list = {};
    
     for i_folder = 1 : length(data_folder_list)
        % load data
        data_path = data_folder_list{i_folder};
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
        disp(['Collecting from ' subject_id]);
        load([data_path filesep 'analysis' filesep date '_' subject_id '_results.mat']);

        % append data from this subject to containers for all subjects
        condition_stance_foot_list = [condition_stance_foot_list; condition_stance_foot_list_session]; %#ok<AGROW>
        condition_perturbation_list = [condition_perturbation_list; condition_perturbation_list_session]; %#ok<AGROW>
        condition_delay_list = [condition_delay_list; condition_delay_list_session]; %#ok<AGROW>
        condition_index_list = [condition_index_list; condition_index_list_session]; %#ok<AGROW>
        condition_experimental_list = [condition_experimental_list; condition_experimental_list_session]; %#ok<AGROW>
        condition_stimulus_list = [condition_stimulus_list; condition_stimulus_list_session]; %#ok<AGROW>
        condition_day_list = [condition_day_list; condition_day_list_session]; %#ok<AGROW>
        [subject_list{length(subject_list)+(1 : length(condition_stance_foot_list_session))}] = deal(subject_id); %#ok<AGROW>
     end
    
     condition_group_assignment = study_settings.get('condition_group_assignment');
     condition_group_list = cell(size(condition_stimulus_list));
     condition_stimulus_list = condition_perturbation_list;
    for i_stretch = 1 : length(condition_stimulus_list)
         % find the row that the subject of this stretch matches in the foot
        group_assignment_subject_row = find(strcmp(subject_list{i_stretch}, condition_group_assignment(:,1)));
        
           
        group_assignment_foot_col = find(strcmp(condition_stance_foot_list{i_stretch}, condition_group_assignment(group_assignment_subject_row,:)));
        if isempty(group_assignment_foot_col)
            group_assignment_foot_col = find(strcmp(condition_group_assignment(group_assignment_subject_row,:), 'BOTH')); 
        end
        
        if group_assignment_foot_col == 2
            condition_group_list{i_stretch} = 'EARLY';
        elseif group_assignment_foot_col == 3
            condition_group_list{i_stretch} = 'LATE';           
        elseif group_assignment_foot_col == 4
            condition_group_list{i_stretch} = 'NO';
        else
            error('Something wrong with the group matching: No match found')
        end
    end
    
    variables_to_save.condition_stance_foot_list = condition_stance_foot_list;
    variables_to_save.condition_perturbation_list = condition_perturbation_list;
    variables_to_save.condition_stimulus_list = condition_perturbation_list;
    variables_to_save.condition_delay_list = condition_delay_list;
    variables_to_save.condition_index_list = condition_index_list;
    variables_to_save.condition_experimental_list = condition_experimental_list;
    variables_to_save.condition_stimulus_list = condition_stimulus_list;
    variables_to_save.condition_day_list = condition_day_list;
    variables_to_save.condition_group_list = condition_group_list;
    save('results', '-struct', 'variables_to_save');
end