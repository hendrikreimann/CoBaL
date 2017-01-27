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
    load('subjectInfo.mat', 'date', 'subject_id');
    study_settings = loadSettingsFile(['..' filesep 'studySettings.txt']);
    data_custodian = WalkingDataCustodian();
    number_of_stretch_variables = length(data_custodian.stretch_variable_names);
    
    % make containers to hold the data
    data_subject = cell(number_of_stretch_variables, 1);
    condition_stance_foot_list_subject = {};
    condition_perturbation_list_subject = {};
    condition_delay_list_subject = {};
    condition_index_list_subject = {};
    condition_experimental_list_subject = {};
    condition_stimulus_list_subject = {};
    condition_day_list_subject = {};
    
    % analyze and store data
    for i_condition = 1 : length(condition_list)
        condition = condition_list{i_condition};
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            % load and prepare data
            data_custodian.prepareBasicVariables(condition, i_trial);
            load(['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'relevantDataStretches')]);
            data_trial = data_custodian.calculateStretchVariables(stretch_start_times, stretch_end_times, condition_stance_foot_list_trial);
            
            % append the data and condition lists from this trial to the total lists
            for i_variable = 1 : number_of_stretch_variables
                data_subject{i_variable} = [data_subject{i_variable} data_trial{i_variable}];
            end
            condition_stance_foot_list_subject = [condition_stance_foot_list_subject; condition_stance_foot_list_trial]; %#ok<AGROW>
            condition_perturbation_list_subject = [condition_perturbation_list_subject; condition_perturbation_list_trial]; %#ok<AGROW>
            condition_delay_list_subject = [condition_delay_list_subject; condition_delay_list_trial]; %#ok<AGROW>
            condition_index_list_subject = [condition_index_list_subject; condition_index_list_trial]; %#ok<AGROW>
            condition_experimental_list_subject = [condition_experimental_list_subject; condition_experimental_list_trial]; %#ok<AGROW>
            condition_stimulus_list_subject = [condition_stimulus_list_subject; condition_stimulus_list_trial]; %#ok<AGROW>
            condition_day_list_subject = [condition_day_list_subject; condition_day_list_trial]; %#ok<AGROW>
        end
    end
    
    % calculate some subject-level data and report
    number_of_stretches_subject = length(condition_stance_foot_list_subject);
    
    % extract indicators for control
    number_of_conditions_control = size(study_settings.conditions_control, 1);
    conditions_control_indicators = false(number_of_stretches_subject, number_of_conditions_control);
    for i_condition = 1 : number_of_conditions_control
        stance_foot_indicator = strcmp(condition_stance_foot_list_subject, study_settings.conditions_control(i_condition, 1));
        perturbation_indicator = strcmp(condition_perturbation_list_subject, study_settings.conditions_control(i_condition, 2));
        delay_indicator = strcmp(condition_delay_list_subject, study_settings.conditions_control(i_condition, 3));
        index_indicator = strcmp(condition_index_list_subject, study_settings.conditions_control(i_condition, 4));
        experimental_indicator = strcmp(condition_experimental_list_subject, study_settings.conditions_control(i_condition, 5));
        stimulus_indicator = strcmp(condition_stimulus_list_subject, study_settings.conditions_control(i_condition, 6));
        day_indicator = strcmp(condition_day_list_subject, study_settings.conditions_control(i_condition, 7));

        this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
        conditions_control_indicators(:, i_condition) = this_condition_indicator;
    end
    
    % report control
    trials_per_condition_control = sum(conditions_control_indicators)';
    conditions_control_with_number = study_settings.conditions_control;
    for i_condition = 1 : number_of_conditions_control
        conditions_control_with_number{i_condition, size(study_settings.conditions_control, 2)+1} = num2str(trials_per_condition_control(i_condition));
    end
    conditions_control_with_labels = [study_settings.condition_labels 'number of stretches'; conditions_control_with_number];
    disp('Control conditions:')
    disp(conditions_control_with_labels);

    % extract indicators for conditions to analyze
    number_of_conditions_to_analyze = size(study_settings.conditions_to_analyze, 1);
    conditions_to_analyze_indicators = false(number_of_stretches_subject, number_of_conditions_to_analyze);
    for i_condition = 1 : number_of_conditions_to_analyze
        stance_foot_indicator = strcmp(condition_stance_foot_list_subject, study_settings.conditions_to_analyze(i_condition, 1));
        perturbation_indicator = strcmp(condition_perturbation_list_subject, study_settings.conditions_to_analyze(i_condition, 2));
        delay_indicator = strcmp(condition_delay_list_subject, study_settings.conditions_to_analyze(i_condition, 3));
        index_indicator = strcmp(condition_index_list_subject, study_settings.conditions_to_analyze(i_condition, 4));
        experimental_indicator = strcmp(condition_experimental_list_subject, study_settings.conditions_to_analyze(i_condition, 5));
        stimulus_indicator = strcmp(condition_stimulus_list_subject, study_settings.conditions_to_analyze(i_condition, 6));
        day_indicator = strcmp(condition_day_list_subject, study_settings.conditions_to_analyze(i_condition, 7));

        this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
        conditions_to_analyze_indicators(:, i_condition) = this_condition_indicator;
    end
    
    % report conditions to analyze
    trials_per_condition_to_analyze = sum(conditions_to_analyze_indicators)';
    conditions_to_analyze_with_number = study_settings.conditions_to_analyze;
    for i_condition = 1 : number_of_conditions_to_analyze
        conditions_to_analyze_with_number{i_condition, size(study_settings.conditions_to_analyze, 2)+1} = num2str(trials_per_condition_to_analyze(i_condition));
    end
    conditions_to_analyze_with_labels = [study_settings.condition_labels 'number of stretches'; conditions_to_analyze_with_number];
    disp('Conditions to analyze:')
    disp(conditions_to_analyze_with_labels);
    
    disp(['Number of control stretches: ' num2str(sum(trials_per_condition_control))]);
    disp(['Number of stimulus stretches: ' num2str(sum(trials_per_condition_to_analyze))]);
    disp(['Number of unassigned stretches: ' num2str(number_of_stretches_subject - sum(trials_per_condition_control) - sum(trials_per_condition_to_analyze))]);
    
    % save data
    variable_data_subject = data_subject; %#ok<NASGU>
    variable_names_subject = data_custodian.stretch_variable_names; %#ok<NASGU>
    
    results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'results')];
    save ...
      ( ...
        results_file_name, ...
        'variable_data_subject', ...
        'variable_names_subject', ...
        'condition_stance_foot_list_subject', ...
        'condition_perturbation_list_subject', ...
        'condition_delay_list_subject', ...
        'condition_index_list_subject', ...
        'condition_experimental_list_subject', ...
        'condition_stimulus_list_subject', ...
        'condition_day_list_subject' ...
      )
end

          

