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
    variables_to_collect_long = study_settings.get('variables_to_collect_long');
    number_of_variables_to_collect = size(variables_to_collect, 1);
    number_of_variables_to_collect_long = size(variables_to_collect_long, 1);
    conditions_settings = study_settings.get('conditions');
    condition_labels = conditions_settings(:, 1)';
    condition_source_variables = conditions_settings(:, 2)';
    number_of_condition_labels = length(condition_labels);
    conditions = struct;
    conditions_long = struct;
    for i_condition = 1 : number_of_condition_labels
        conditions.(condition_source_variables{i_condition}) = {};
        conditions_long.(condition_source_variables{i_condition}) = {};
    end

%     origin_trial_list = [];
%     origin_start_time_list = [];
%     origin_end_time_list = [];
    time_list = [];
    variable_data = cell(number_of_variables_to_collect, 1);
    step_time_data = [];
    origin_trial_number_data = [];
    origin_stretch_start_time_data = [];
    origin_stretch_end_time_data = [];
    origin_session_folder_data = {};
    variables_to_save = struct;
    
    variable_data_long = cell(number_of_variables_to_collect_long, 1);
    step_time_data_long = [];
    variables_to_save_long = struct;
   
    for i_folder = 1 : length(data_folder_list)
        % load data
        data_path = data_folder_list{i_folder};
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
        disp(['Collecting from ' subject_id]);
        results_data = load([data_path filesep 'analysis' filesep date '_' subject_id '_results.mat']);
        if ~isempty(variables_to_collect_long)
            results_data_long = load([data_path filesep 'analysis' filesep date '_' subject_id '_longStretchResults.mat']);
        end

        % extract and store results
        for i_condition = 1 : number_of_condition_labels
            conditions.(condition_source_variables{i_condition}) = [conditions.(condition_source_variables{i_condition}); results_data.conditions_session.(condition_source_variables{i_condition}) ];
        end
%         origin_trial_list = [origin_trial_list; results_data.origin_trial_list_session]; %#ok<AGROW>
%         origin_start_time_list = [origin_start_time_list; results_data.origin_start_time_list_session]; %#ok<AGROW>
%         origin_end_time_list = [origin_end_time_list; results_data.origin_end_time_list_session]; %#ok<AGROW>
        time_list = [time_list; results_data.time_list_session]; %#ok<AGROW>
        
        for i_variable = 1 : number_of_variables_to_collect
            % load and extract data
            new_variable_name = variables_to_collect{i_variable, 1};
            source_variable_name = variables_to_collect{i_variable, 2};
            this_variable_source_type = variables_to_collect{i_variable, 3};
            if strcmp(this_variable_source_type, 'stretch')
                this_variable_source_index = find(strcmp(results_data.stretch_names_session, source_variable_name), 1, 'first');
                if isempty(this_variable_source_index)
                    error(['Variable not found: ' source_variable_name])
                end
                this_variable_data = results_data.stretch_data_session{this_variable_source_index};
            end
            if strcmp(this_variable_source_type, 'response')
                this_variable_source_index = find(strcmp(results_data.response_names_session, source_variable_name), 1, 'first');
                if isempty(this_variable_source_index)
                    error(['Variable not found: ' source_variable_name])
                end
                this_variable_data = results_data.response_data_session{this_variable_source_index};
            end
            if strcmp(this_variable_source_type, 'analysis')
                this_variable_source_index = find(strcmp(results_data.analysis_names_session, source_variable_name), 1, 'first');
                if isempty(this_variable_source_index)
                    error(['Variable not found: ' source_variable_name])
                end
                this_variable_data = results_data.analysis_data_session{this_variable_source_index};
            end
            
            % store
            variable_data{i_variable} = [variable_data{i_variable} this_variable_data];
        end
        step_time_source_index = find(strcmp(results_data.stretch_names_session, 'step_time'), 1, 'first');
        step_time_data = [step_time_data results_data.stretch_data_session{step_time_source_index}]; %#ok<AGROW>
        
        origin_trial_number_data = [origin_trial_number_data; results_data.origin_trial_list_session]; %#ok<AGROW>
        origin_stretch_start_time_data = [origin_stretch_start_time_data; results_data.origin_start_time_list_session]; %#ok<AGROW>
        origin_stretch_end_time_data = [origin_stretch_end_time_data; results_data.origin_end_time_list_session]; %#ok<AGROW>
        session_folder_list = cell(size(results_data.origin_trial_list_session));
        data_path_split = strsplit(data_path, filesep);
        session_folder = data_path_split{end};
        [session_folder_list{:}] = deal(session_folder);
        origin_session_folder_data = [origin_session_folder_data; session_folder_list]; %#ok<AGROW>

        % extract and store long stretches
        if ~isempty(variables_to_collect_long)
            for i_condition = 1 : number_of_condition_labels
                conditions_long.(condition_source_variables{i_condition}) = [conditions_long.(condition_source_variables{i_condition}); results_data_long.long_stretch_conditions_session.(condition_source_variables{i_condition}) ];
            end
        
            for i_variable = 1 : number_of_variables_to_collect_long
                % load and extract data
                new_variable_name = variables_to_collect_long{i_variable, 1};
                source_variable_name = variables_to_collect_long{i_variable, 2};
                this_variable_source_type = variables_to_collect_long{i_variable, 3};
                if strcmp(this_variable_source_type, 'stretch')
                    this_variable_source_index = find(strcmp(results_data_long.long_stretch_data_labels_session, source_variable_name), 1, 'first');
                    if isempty(this_variable_source_index)
                        error(['Variable not found: ' source_variable_name])
                    end
                    this_variable_data = results_data_long.long_stretch_data_session{this_variable_source_index};
                end
                if strcmp(this_variable_source_type, 'response')
                    this_variable_source_index = find(strcmp(results_data_long.long_stretch_data_labels_session, source_variable_name), 1, 'first');
                    if isempty(this_variable_source_index)
                        error(['Variable not found: ' source_variable_name])
                    end
                    this_variable_data = results_data_long.long_stretch_response_data_session{this_variable_source_index};
                end
                if strcmp(this_variable_source_type, 'analysis')
                    this_variable_source_index = find(strcmp(results_data_long.analysis_names_session, source_variable_name), 1, 'first');
                    if isempty(this_variable_source_index)
                        error(['Variable not found: ' source_variable_name])
                    end
                    this_variable_data = results_data_long.analysis_data_session{this_variable_source_index};
                end

                % store
                variable_data_long{i_variable} = [variable_data_long{i_variable} this_variable_data];
            end
            step_time_data_long = [step_time_data_long results_data_long.step_time_data_session]; %#ok<AGROW>
        end

        
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
    
    if ~isempty(variables_to_collect_long)
        variables_to_save.variable_data_long = variable_data_long;
        variables_to_save.variable_names_long = variables_to_collect_long(:, 1);
        variables_to_save.step_time_data_long = step_time_data_long;
        variables_to_save.conditions_long = conditions_long;
    end
    
    save('results', '-struct', 'variables_to_save');





end