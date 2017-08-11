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
    study_settings_file = '';
    if ~exist('studySettings.txt', 'file')
        error('No studySettings.txt file found. This function should be run from a study folder')
    end    
    study_settings = SettingsCustodian('studySettings.txt');
    
    data_folder_list = determineDataStructure(subjects);

    %% collect data from all data folders
    variables_to_analyze = study_settings.get('variables_to_analyze');
    number_of_variables_to_analyze = size(variables_to_analyze, 1);
    condition_stance_foot_list_all = {};
    condition_perturbation_list_all = {};
    condition_delay_list_all = {};
    condition_index_list_all = {};
    condition_experimental_list_all = {};
    condition_stimulus_list_all = {};
    condition_day_list_all = {};
    subject_list_all = {};
    origin_trial_list_all = [];
    origin_start_time_list_all = [];
    origin_end_time_list_all = [];
    variable_data_all = cell(number_of_variables_to_analyze, 1);
    if plot_settings.get('plot_response')
        response_data_all = cell(number_of_variables_to_analyze, 1);
    end
    step_time_data = [];
    
    for i_folder = 1 : length(data_folder_list)
        % load data
        data_path = data_folder_list{i_folder};
        load([data_path filesep 'subjectInfo.mat'], 'date', 'subject_id');
        load([data_path filesep 'analysis' filesep date '_' subject_id '_results.mat']);

        % append data from this subject to containers for all subjects
        condition_stance_foot_list_all = [condition_stance_foot_list_all; condition_stance_foot_list_session]; %#ok<AGROW>
        condition_perturbation_list_all = [condition_perturbation_list_all; condition_perturbation_list_session]; %#ok<AGROW>
        condition_delay_list_all = [condition_delay_list_all; condition_delay_list_session]; %#ok<AGROW>
        condition_index_list_all = [condition_index_list_all; condition_index_list_session]; %#ok<AGROW>
        condition_experimental_list_all = [condition_experimental_list_all; condition_experimental_list_session]; %#ok<AGROW>
        condition_stimulus_list_all = [condition_stimulus_list_all; condition_stimulus_list_session]; %#ok<AGROW>
        condition_day_list_all = [condition_day_list_all; condition_day_list_session]; %#ok<AGROW>
        origin_trial_list_all = [origin_trial_list_all; origin_trial_list_session]; %#ok<AGROW>
        origin_start_time_list_all = [origin_start_time_list_all; origin_start_time_list_session]; %#ok<AGROW>
        origin_end_time_list_all = [origin_end_time_list_all; origin_end_time_list_session]; %#ok<AGROW>
        for i_variable = 1 : number_of_variables_to_analyze
            % load and extract data
            this_variable_name = variables_to_analyze{i_variable, 1};
            index_in_saved_data = find(strcmp(variable_names_session, this_variable_name), 1, 'first');
            this_variable_data = variable_data_session{index_in_saved_data}; %#ok<USENS>
            if plot_settings.get('plot_response')
                this_response_data = response_data_session{index_in_saved_data}; %#ok<USENS>
            end
            
            % store
            variable_data_all{i_variable} = [variable_data_all{i_variable} this_variable_data];
            if plot_settings.get('plot_response')
                response_data_all{i_variable} = [response_data_all{i_variable} this_response_data];
            end
        end
        if strcmp(plot_settings.get('time_plot_style'), 'scaled_to_comparison_mean') || strcmp(plot_settings.get('time_plot_style'), 'scaled_to_condition_mean')
            index_in_saved_data = find(strcmp(variable_names_session, 'step_time'), 1, 'first');
            this_step_time_data = variable_data_session{index_in_saved_data};
            step_time_data = [step_time_data this_step_time_data];
        end
    end
    
%     % save data for quick stats
%     variables_to_save = struct;
%     variables_to_save.variable_data_all = variable_data_all;
%     if plot_settings.get('plot_response')
%         variables_to_save.response_data_all = response_data_all;
%     end
%     variables_to_save.variable_names = variables_to_plot;
%     variables_to_save.condition_stance_foot_list_all = condition_stance_foot_list_all;
%     variables_to_save.condition_perturbation_list_all = condition_perturbation_list_all;
%     variables_to_save.condition_delay_list_all = condition_delay_list_all;
%     variables_to_save.condition_index_list_all = condition_index_list_all;
%     variables_to_save.condition_experimental_list_all = condition_experimental_list_all;
%     variables_to_save.condition_stimulus_list_all = condition_stimulus_list_all;
%     variables_to_save.condition_day_list_all = condition_day_list_all;
%     save('results', '-struct', 'variables_to_save');





end