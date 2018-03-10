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
    conditions_settings = study_settings.get('conditions');
    condition_labels = conditions_settings(:, 1)';
    condition_source_variables = conditions_settings(:, 2)';
    number_of_condition_labels = length(condition_labels);
    conditions = struct;
    for i_condition = 1 : number_of_condition_labels
        conditions.(condition_source_variables{i_condition}) = {};
    end

    origin_trial_list = [];
    origin_start_time_list = [];
    origin_end_time_list = [];
    time_list = [];
    variable_data = cell(number_of_variables_to_collect, 1);
    step_time_data = [];
%     pushoff_time_data = [];
    origin_trial_number_data = [];
    origin_stretch_start_time_data = [];
    origin_stretch_end_time_data = [];
    origin_session_folder_data = {};
    variables_to_save = struct;
   
    for i_folder = 1 : length(data_folder_list)
        % load data
        data_path = data_folder_list{i_folder};
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
        disp(['Collecting from ' subject_id]);
        load([data_path filesep 'analysis' filesep date '_' subject_id '_results.mat']); %#ok<LOAD>

        % append data from this subject to containers for all subjects
        for i_condition = 1 : number_of_condition_labels
            conditions.(condition_source_variables{i_condition}) = [conditions.(condition_source_variables{i_condition}); conditions_session.(condition_source_variables{i_condition}) ];
        end
        origin_trial_list = [origin_trial_list; origin_trial_list_session]; %#ok<AGROW>
        origin_start_time_list = [origin_start_time_list; origin_start_time_list_session]; %#ok<AGROW>
        origin_end_time_list = [origin_end_time_list; origin_end_time_list_session]; %#ok<AGROW>
        time_list = [time_list; time_list_session]; %#ok<AGROW>
        
        for i_variable = 1 : number_of_variables_to_collect
            % load and extract data
            this_variable_name = variables_to_collect{i_variable, 1};
            this_variable_source_type = variables_to_collect{i_variable, 2};
            if strcmp(this_variable_source_type, 'stretch')
                this_variable_source_index = find(strcmp(stretch_names_session, this_variable_name), 1, 'first');
                if isempty(this_variable_source_index)
                    error(['Variable not found: ' this_variable_name])
                end
                this_variable_data = stretch_data_session{this_variable_source_index}; %#ok<IDISVAR,USENS>
            end
            if strcmp(this_variable_source_type, 'response')
                this_variable_source_index = find(strcmp(response_names_session, this_variable_name), 1, 'first');
                if isempty(this_variable_source_index)
                    error(['Variable not found: ' this_variable_name])
                end
                this_variable_data = response_data_session{this_variable_source_index}; %#ok<IDISVAR,USENS>
            end
            if strcmp(this_variable_source_type, 'analysis')
                this_variable_source_index = find(strcmp(analysis_names_session, this_variable_name), 1, 'first');
                if isempty(this_variable_source_index)
                    error(['Variable not found: ' this_variable_name])
                end
                this_variable_data = analysis_data_session{this_variable_source_index}; %#ok<IDISVAR,USENS>
            end
            
            % store
            variable_data{i_variable} = [variable_data{i_variable} this_variable_data];
        end
        step_time_source_index = find(strcmp(stretch_names_session, 'step_time'), 1, 'first');
        step_time_data = [step_time_data stretch_data_session{step_time_source_index}]; %#ok<AGROW>
        
        origin_trial_number_data = [origin_trial_number_data; origin_trial_list_session]; %#ok<AGROW>
        origin_stretch_start_time_data = [origin_stretch_start_time_data; origin_start_time_list_session]; %#ok<AGROW>
        origin_stretch_end_time_data = [origin_stretch_end_time_data; origin_end_time_list_session]; %#ok<AGROW>
        session_folder_list = cell(size(origin_trial_list_session));
        data_path_split = strsplit(data_path, filesep);
        session_folder = data_path_split{end};
        [session_folder_list{:}] = deal(session_folder);
        origin_session_folder_data = [origin_session_folder_data; session_folder_list]; %#ok<AGROW>

        
        
        % HR: removed, because the whole pushoff thing shouldn't have to be treated in a special way
        % this might generate problems for legacy formats in Vision and GVS, resolve then
%         pushoff_time_source_index = find(strcmp(response_names_session, 'pushoff_time'), 1, 'first');
%         pushoff_time_data = [pushoff_time_data stretch_data_session{pushoff_time_source_index}];

    end
    
    
    
    %% save
    variables_to_save.variable_data = variable_data;
    variables_to_save.variable_names = variables_to_collect(:, 1);
    variables_to_save.step_time_data = step_time_data;
    variables_to_save.conditions = conditions;
    variables_to_save.time_list = time_list;
    variables_to_save.origin_session_folder_data = origin_session_folder_data;
    variables_to_save.origin_trial_number_data = origin_trial_number_data;
    variables_to_save.origin_stretch_start_time_data = origin_stretch_start_time_data;
    variables_to_save.origin_stretch_end_time_data = origin_stretch_end_time_data; %#ok<STRNU>
    save('results', '-struct', 'variables_to_save');





end