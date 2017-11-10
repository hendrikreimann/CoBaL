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

% consolidate data from all subjects

function collectPopulationResults(varargin)

    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'subjects', [])
    parse(parser, varargin{:})
    subjects = parser.Results.subjects;

    % load settings
    if ~exist('studySettings.txt', 'file')
        error('No studySettings.txt file found. This function should be run from a study folder')
    end    
    study_settings = SettingsCustodian('studySettings.txt');
    
    data_folder_list = determineDataStructure(subjects);

    %% collect data from all data folders
    variables_to_collect = study_settings.get('variables_to_collect');
    
    number_of_variables_to_collect = size(variables_to_collect, 1);
    condition_stance_foot_list = {};
    condition_perturbation_list = {};
    condition_direction_list = {};
    condition_delay_list = {};
    condition_index_list = {};
    condition_experimental_list = {};
    condition_stimulus_list = {};
    condition_day_list = {};
    condition_group_list = {};
    subject_list = {};
    origin_trial_list = [];
    origin_start_time_list = [];
    origin_end_time_list = [];
    time_list = [];
    variable_data = cell(number_of_variables_to_collect, 1);
    step_time_data = [];
    pushoff_time_data = [];
   
    for i_folder = 1 : length(data_folder_list)
        % load data
        data_path = data_folder_list{i_folder};
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
        disp(['Collecting from ' subject_id]);
        load([data_path filesep 'analysis' filesep date '_' subject_id '_results.mat']);

        % append data from this subject to containers for all subjects
        condition_stance_foot_list = [condition_stance_foot_list; condition_stance_foot_list_session]; %#ok<AGROW>
        condition_perturbation_list = [condition_perturbation_list; condition_perturbation_list_session]; %#ok<AGROW>
        condition_direction_list = [condition_direction_list; condition_direction_list_session]; %#ok<AGROW>
        condition_delay_list = [condition_delay_list; condition_delay_list_session]; %#ok<AGROW>
        condition_index_list = [condition_index_list; condition_index_list_session]; %#ok<AGROW>
        condition_experimental_list = [condition_experimental_list; condition_experimental_list_session]; %#ok<AGROW>
        condition_stimulus_list = [condition_stimulus_list; condition_stimulus_list_session]; %#ok<AGROW>
        condition_day_list = [condition_day_list; condition_day_list_session]; %#ok<AGROW>
        condition_group_list = [condition_group_list; condition_group_list_session]; %#ok<AGROW>
        [subject_list{length(subject_list)+(1 : length(condition_stance_foot_list_session))}] = deal(subject_id); %#ok<AGROW>
        origin_trial_list = [origin_trial_list; origin_trial_list_session]; %#ok<AGROW>
        origin_start_time_list = [origin_start_time_list; origin_start_time_list_session]; %#ok<AGROW>
        origin_end_time_list = [origin_end_time_list; origin_end_time_list_session]; %#ok<AGROW>
        time_list = [time_list; time_list_session]; %#ok<AGROW>
        
        for i_variable = 1 : number_of_variables_to_collect
            % load and extract data
            this_variable_name = variables_to_collect{i_variable, 1};
            this_variable_source_type = variables_to_collect{i_variable, 2};
            if strcmp(this_variable_source_type, 'response')
                this_variable_source_index = find(strcmp(response_names_session, this_variable_name), 1, 'first');
                this_variable_data = response_data_session{this_variable_source_index}; %#ok<USENS>
            end
            if strcmp(this_variable_source_type, 'analysis')
                this_variable_source_index = find(strcmp(analysis_names_session, this_variable_name), 1, 'first');
                this_variable_data = analysis_data_session{this_variable_source_index}; %#ok<USENS>
            end
            
            % store
            variable_data{i_variable} = [variable_data{i_variable} this_variable_data];
        end
        step_time_source_index = find(strcmp(response_names_session, 'step_time'), 1, 'first');
        step_time_data = [step_time_data stretch_data_session{step_time_source_index}];
        pushoff_time_source_index = find(strcmp(response_names_session, 'pushoff_time'), 1, 'first');
        pushoff_time_data = [pushoff_time_data stretch_data_session{pushoff_time_source_index}];

    end
    subject_list = subject_list';
    
    
    %% make time categorical variable
    time_category = cell(size(time_list));
    time_category_borders = study_settings.get('time_category_borders');
    for i_point = 1 : length(time_list)
        if time_category_borders(1) < time_list(i_point) && time_list(i_point) < time_category_borders(2)
            time_category{i_point} = 'TIME_ONE';
        end
        if time_category_borders(2) < time_list(i_point) && time_list(i_point) < time_category_borders(3)
            time_category{i_point} = 'TIME_TWO';
        end
        if time_category_borders(3) < time_list(i_point) && time_list(i_point) < time_category_borders(4)
            time_category{i_point} = 'TIME_THREE';
        end
        if time_category_borders(4) < time_list(i_point) && time_list(i_point) < time_category_borders(5)
            time_category{i_point} = 'TIME_FOUR';
        end
    end
    

    
    
    %% save
    variables_to_save = struct;
    variables_to_save.variable_data = variable_data;
    variables_to_save.variable_names = variables_to_collect(:, 1);
    variables_to_save.step_time_data = step_time_data;
    variables_to_save.pushoff_time_data = pushoff_time_data;
    variables_to_save.condition_stance_foot_list = condition_stance_foot_list;
    variables_to_save.condition_perturbation_list = condition_perturbation_list;
    variables_to_save.condition_direction_list = condition_direction_list;
    variables_to_save.condition_stimulus_list = condition_perturbation_list;
    variables_to_save.condition_delay_list = condition_delay_list;
    variables_to_save.condition_index_list = condition_index_list;
    variables_to_save.condition_experimental_list = condition_experimental_list;
    variables_to_save.condition_stimulus_list = condition_stimulus_list;
    variables_to_save.condition_day_list = condition_day_list;
    variables_to_save.condition_group_list = condition_group_list;
    variables_to_save.time_list = time_list;
    variables_to_save.time_category = time_category;
    variables_to_save.subject_list = subject_list; %#ok<STRNU>
    save('results', '-struct', 'variables_to_save');





end