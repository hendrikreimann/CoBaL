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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% input

% analyze the data

% input
% relevantDataStretches.mat

function analyzeData(varargin)
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    study_settings = loadSettingsFile(['..' filesep 'studySettings.txt']);
    load('subjectInfo.mat', 'date', 'subject_id');
%     data_custodian = WalkingDataCustodian(date, subject_id, study_settings.variables_to_analyze);
    data_custodian = WalkingDataCustodian();
    number_of_required_variables = length(data_custodian.required_variable_names);
    
    % make containers to hold the data
    data_total = cell(number_of_required_variables, 1);
    condition_stance_foot_list_total = {};
    condition_perturbation_list_total = {};
    condition_delay_list_total = {};
    condition_index_list_total = {};
    condition_experimental_list_total = {};
    
    % analyze and store data
    for i_condition = 1 : length(condition_list)
        condition = condition_list{i_condition};
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            % load and prepare data
            data_trial = cell(number_of_required_variables, 1);
            data_custodian.prepareData(condition, i_trial);
            load(['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'relevantDataStretches')]);
            number_of_stretches = length(condition_stance_foot_list);
            
            % extract and normalize data from stretches
            for i_stretch = 1 : number_of_stretches
                % time
                this_stretch_start_time = stretch_start_times(i_stretch);
                this_stretch_end_time = stretch_end_times(i_stretch);
                
                for i_variable = 1 : number_of_required_variables
                    time_normalized_data = data_custodian.getTimeNormalizedData(data_custodian.required_variable_names{i_variable}, this_stretch_start_time, this_stretch_end_time);
                    data_trial{i_variable} = [data_trial{i_variable} time_normalized_data];
                end
                
            end
            
            % append the data and condition lists from this trial to the total lists
            for i_variable = 1 : number_of_required_variables
                data_total{i_variable} = [data_total{i_variable} data_trial{i_variable}];
            end
            condition_stance_foot_list_total = [condition_stance_foot_list_total; condition_stance_foot_list];
            condition_perturbation_list_total = [condition_perturbation_list_total; condition_perturbation_list];
            condition_delay_list_total = [condition_delay_list_total; condition_delay_list];
            condition_index_list_total = [condition_index_list_total; condition_index_list];
            condition_experimental_list_total = [condition_experimental_list_total; condition_experimental_list];
        end
    end
    
    % save data
    condition_stance_foot_list = condition_stance_foot_list_total;
    condition_perturbation_list = condition_perturbation_list_total;
    condition_delay_list = condition_delay_list_total;
    condition_index_list = condition_index_list_total;
    condition_experimental_list = condition_experimental_list_total;
    variable_data = data_total;
    variable_names = data_custodian.required_variable_names;
    
    results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'results')];
    save ...
      ( ...
        results_file_name, ...
        'variable_data', ...
        'variable_names', ...
        'condition_stance_foot_list', ...
        'condition_perturbation_list', ...
        'condition_delay_list', ...
        'condition_index_list', ...
        'condition_experimental_list' ...
      )
end

          

