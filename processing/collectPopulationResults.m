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
    variables_to_analyze = study_settings.get('variables_to_analyze');
    number_of_variables_to_analyze = size(variables_to_analyze, 1);
    condition_stance_foot_list = {};
    condition_perturbation_list = {};
    condition_delay_list = {};
    condition_index_list = {};
    condition_experimental_list = {};
    condition_stimulus_list = {};
    condition_day_list = {};
    subject_list = {};
    origin_trial_list = [];
    origin_start_time_list = [];
    origin_end_time_list = [];
    variable_data = cell(number_of_variables_to_analyze, 1);
    response_data = cell(number_of_variables_to_analyze, 1);
    
    for i_folder = 1 : length(data_folder_list)
        % load data
        data_path = data_folder_list{i_folder};
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
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
        origin_trial_list = [origin_trial_list; origin_trial_list_session]; %#ok<AGROW>
        origin_start_time_list = [origin_start_time_list; origin_start_time_list_session]; %#ok<AGROW>
        origin_end_time_list = [origin_end_time_list; origin_end_time_list_session]; %#ok<AGROW>
        for i_variable = 1 : number_of_variables_to_analyze
            % load and extract data
            this_variable_name = variables_to_analyze{i_variable, 1};
            index_in_saved_data = find(strcmp(variable_names_session, this_variable_name), 1, 'first');
            this_variable_data = variable_data_session{index_in_saved_data}; %#ok<USENS>
            this_response_data = response_data_session{index_in_saved_data}; %#ok<USENS>
            
            % store
            variable_data{i_variable} = [variable_data{i_variable} this_variable_data];
            response_data{i_variable} = [response_data{i_variable} this_response_data];
        end
    end
    subject_list = subject_list';
    
    % save data
    variables_to_save = struct;
    variables_to_save.variable_data = variable_data;
    variables_to_save.response_data = response_data;
    variables_to_save.variable_names = variables_to_analyze;
    variables_to_save.condition_stance_foot_list = condition_stance_foot_list;
    variables_to_save.condition_perturbation_list = condition_perturbation_list;
    variables_to_save.condition_delay_list = condition_delay_list;
    variables_to_save.condition_index_list = condition_index_list;
    variables_to_save.condition_experimental_list = condition_experimental_list;
    variables_to_save.condition_stimulus_list = condition_stimulus_list;
    variables_to_save.condition_day_list = condition_day_list;
    variables_to_save.subject_list = subject_list;
    save('results', '-struct', 'variables_to_save');





end