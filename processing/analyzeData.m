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

% analyze the data

% input
% relevantDataStretches.mat

function analyzeData(varargin)
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    load('subjectInfo.mat', 'date', 'subject_id');
    % load settings
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    data_custodian = WalkingDataCustodian();
    number_of_stretch_variables = length(data_custodian.stretch_variable_names);
    
    % make containers to hold the data
    data_session = cell(number_of_stretch_variables, 1);
    condition_stance_foot_list_session = {};
    condition_perturbation_list_session = {};
    condition_delay_list_session = {};
    condition_index_list_session = {};
    condition_experimental_list_session = {};
    condition_stimulus_list_session = {};
    condition_day_list_session = {};
    
    % make containers to store origin information for the stretches
    origin_trial_list_session = [];
    origin_start_time_list_session = [];
    origin_end_time_list_session = [];
    time_list_session = [];
    
    %% analyze and store data
    for i_type = 1 : length(condition_list)
        condition = condition_list{i_type};
        trials_to_process = trial_number_list{i_type};
        for i_trial = trials_to_process
            disp(['i_trial = ' num2str(i_trial)])
            % load and prepare data
            data_custodian.prepareBasicVariables(condition, i_trial);
            load(['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'relevantDataStretches')]);
            data_trial = data_custodian.calculateStretchVariables(stretch_start_times, stretch_end_times, stretch_pushoff_times, condition_stance_foot_list_trial, condition_experimental_list_trial);
            
            % append the data and condition lists from this trial to the total lists
            for i_variable = 1 : number_of_stretch_variables
                data_session{i_variable} = [data_session{i_variable} data_trial{i_variable}];
            end
            condition_stance_foot_list_session = [condition_stance_foot_list_session; condition_stance_foot_list_trial]; %#ok<AGROW>
            condition_perturbation_list_session = [condition_perturbation_list_session; condition_perturbation_list_trial]; %#ok<AGROW>
            condition_delay_list_session = [condition_delay_list_session; condition_delay_list_trial]; %#ok<AGROW>
            condition_index_list_session = [condition_index_list_session; condition_index_list_trial]; %#ok<AGROW>
            condition_experimental_list_session = [condition_experimental_list_session; condition_experimental_list_trial]; %#ok<AGROW>
            condition_stimulus_list_session = [condition_stimulus_list_session; condition_stimulus_list_trial]; %#ok<AGROW>
            condition_day_list_session = [condition_day_list_session; condition_day_list_trial]; %#ok<AGROW>
            
            origin_trial_list_session = [origin_trial_list_session; ones(size(stretch_start_times)) * i_trial]; %#ok<AGROW>
            origin_start_time_list_session = [origin_start_time_list_session; stretch_start_times]; %#ok<AGROW>
            origin_end_time_list_session = [origin_end_time_list_session; stretch_end_times]; %#ok<AGROW>
            time_list = ones(size(stretch_start_times)) * (i_trial - 1) * study_settings.get('trial_length') + stretch_start_times;
            time_list_session = [time_list_session; time_list]; %#ok<AGROW>
        end
        disp(['Finished condition "' condition '".'])
    end
    
    %% calculate some subject-level data and report
    number_of_stretches_session = length(condition_stance_foot_list_session);
    
    % extract indicators for control
    conditions_control = study_settings.get('conditions_control');
    number_of_conditions_control = size(study_settings.get('conditions_control'), 1);
    conditions_control_indicators = false(number_of_stretches_session, number_of_conditions_control);
    for i_condition = 1 : number_of_conditions_control
        stance_foot_indicator = strcmp(condition_stance_foot_list_session, conditions_control(i_condition, 1));
        perturbation_indicator = strcmp(condition_perturbation_list_session, conditions_control(i_condition, 2));
        delay_indicator = strcmp(condition_delay_list_session, conditions_control(i_condition, 3));
        index_indicator = strcmp(condition_index_list_session, conditions_control(i_condition, 4));
        experimental_indicator = strcmp(condition_experimental_list_session, conditions_control(i_condition, 5));
        stimulus_indicator = strcmp(condition_stimulus_list_session, conditions_control(i_condition, 6));
        day_indicator = strcmp(condition_day_list_session, conditions_control(i_condition, 7));

        this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
        conditions_control_indicators(:, i_condition) = this_condition_indicator;
    end
    
    % extract indicators for conditions to analyze
    conditions_to_analyze = study_settings.get('conditions_to_analyze');
    number_of_conditions_to_analyze = size(conditions_to_analyze, 1);
    conditions_to_analyze_indicators = false(number_of_stretches_session, number_of_conditions_to_analyze);
    for i_condition = 1 : number_of_conditions_to_analyze
        stance_foot_indicator = strcmp(condition_stance_foot_list_session, conditions_to_analyze(i_condition, 1));
        if study_settings.get('analyze_total_response')
            perturbation_indicator = strcmp(condition_perturbation_list_all,'ILLUSION_RIGHT') | strcmp(condition_perturbation_list_all,'ILLUSION_LEFT');
        else
            perturbation_indicator = strcmp(condition_perturbation_list_session, conditions_to_analyze(i_condition, 2)); 
        end
        delay_indicator = strcmp(condition_delay_list_session, conditions_to_analyze(i_condition, 3));
        index_indicator = strcmp(condition_index_list_session, conditions_to_analyze(i_condition, 4));
        experimental_indicator = strcmp(condition_experimental_list_session, conditions_to_analyze(i_condition, 5));
        stimulus_indicator = strcmp(condition_stimulus_list_session, conditions_to_analyze(i_condition, 6));
        day_indicator = strcmp(condition_day_list_session, conditions_to_analyze(i_condition, 7));
        
        this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
        conditions_to_analyze_indicators(:, i_condition) = this_condition_indicator;
    end
    
    % check the unassigned stretches
    assigned_stretch_indicator = sum([conditions_to_analyze_indicators conditions_control_indicators], 2);
    unassigned_stretch_indicator = ~assigned_stretch_indicator;
    unassigned_stretch_indices = find(unassigned_stretch_indicator);
    unassigned_stretch_labels_session = cell(length(unassigned_stretch_indices), length(study_settings.get('condition_labels')));
    for i_stretch = 1 : length(unassigned_stretch_indices)
        unassigned_stretch_labels_session(i_stretch, 1) = condition_stance_foot_list_session(unassigned_stretch_indices(i_stretch));
        unassigned_stretch_labels_session(i_stretch, 2) = condition_perturbation_list_session(unassigned_stretch_indices(i_stretch));
        unassigned_stretch_labels_session(i_stretch, 3) = condition_delay_list_session(unassigned_stretch_indices(i_stretch));
        unassigned_stretch_labels_session(i_stretch, 4) = condition_index_list_session(unassigned_stretch_indices(i_stretch));
        unassigned_stretch_labels_session(i_stretch, 5) = condition_experimental_list_session(unassigned_stretch_indices(i_stretch));
        unassigned_stretch_labels_session(i_stretch, 6) = condition_stimulus_list_session(unassigned_stretch_indices(i_stretch));
        unassigned_stretch_labels_session(i_stretch, 7) = condition_day_list_session(unassigned_stretch_indices(i_stretch));
    end
    wd = unassigned_stretch_labels_session;
    [~, idx] = unique(strcat(wd(:,1), wd(:,2), wd(:,3), wd(:,4), wd(:,5), wd(:,6), wd(:,7)));
    unassigned_stretch_labels = wd(idx,:);

    % report control
    trials_per_condition_control = sum(conditions_control_indicators)';
    conditions_control_with_number = conditions_control;
    for i_condition = 1 : number_of_conditions_control
        conditions_control_with_number{i_condition, size(conditions_control, 2)+1} = num2str(trials_per_condition_control(i_condition));
    end
    conditions_control_with_labels = [study_settings.get('condition_labels') 'number of stretches'; conditions_control_with_number];
    disp('Control conditions:')
    disp(conditions_control_with_labels);

    % report conditions to analyze
    trials_per_condition_to_analyze = sum(conditions_to_analyze_indicators)';
    conditions_to_analyze_with_number = conditions_to_analyze;
    for i_condition = 1 : number_of_conditions_to_analyze
        conditions_to_analyze_with_number{i_condition, size(conditions_to_analyze, 2)+1} = num2str(trials_per_condition_to_analyze(i_condition));
    end
    conditions_to_analyze_with_labels = [study_settings.get('condition_labels') 'number of stretches'; conditions_to_analyze_with_number];
    disp('Conditions to analyze:')
    disp(conditions_to_analyze_with_labels);
    
    disp('Stretches with these conditions were found but not analyzed:')
    disp([study_settings.get('condition_labels'); unassigned_stretch_labels])
    
    disp(['Number of control stretches: ' num2str(sum(trials_per_condition_control))]);
    disp(['Number of stimulus stretches: ' num2str(sum(trials_per_condition_to_analyze))]);
    disp(['Number of un-analyzed stretches: ' num2str(number_of_stretches_session - sum(trials_per_condition_control) - sum(trials_per_condition_to_analyze))]);
    
    variable_data_session = data_session; %#ok<NASGU>
    variable_names_session = data_custodian.stretch_variable_names;
    
    %% calculate response (i.e. difference from control mean)
    response_data_session = {};
    if ~isempty(conditions_control)
        % prepare container
        response_data_session = cell(size(data_session));
        for i_variable = 1 : number_of_stretch_variables
            response_data_session{i_variable} = zeros(size(data_session{i_variable}));
        end        
        
        % new code - go stretch by stretch
        number_of_stretches = size(variable_data_session{1}, 2);
        for i_stretch = 1 : number_of_stretches
            % determine conditions and applicable control
            this_stretch_condition_string = ...
              { ...
                condition_stance_foot_list_session{i_stretch}, ...
                condition_perturbation_list_session{i_stretch}, ...
                condition_delay_list_session{i_stretch}, ...
                condition_index_list_session{i_stretch}, ...
                condition_experimental_list_session{i_stretch}, ...
                condition_stimulus_list_session{i_stretch}, ...
                condition_day_list_session{i_stretch} ...
              };
            applicable_control_condition_index = findApplicableControlConditionIndex(this_stretch_condition_string, conditions_control);
            applicable_control_condition_labels = conditions_control(applicable_control_condition_index, :);
            
            % determine indicator for control
            stance_foot_indicator = strcmp(condition_stance_foot_list_session, applicable_control_condition_labels{1});
            perturbation_indicator = strcmp(condition_perturbation_list_session, applicable_control_condition_labels{2});
            delay_indicator = strcmp(condition_delay_list_session, applicable_control_condition_labels{3});
            index_indicator = strcmp(condition_index_list_session, applicable_control_condition_labels{4});
            experimental_indicator = strcmp(condition_experimental_list_session, applicable_control_condition_labels{5});
            stimulus_indicator = strcmp(condition_stimulus_list_session, applicable_control_condition_labels{6});
            day_indicator = strcmp(condition_day_list_session, applicable_control_condition_labels{7});
            this_condition_control_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
            
            % calculate responses
            for i_variable = 1 : number_of_stretch_variables
                % calculate control mean
                data_this_variable = data_session{i_variable};
                this_condition_control_data = data_this_variable(:, this_condition_control_indicator);
                this_condition_control_mean = mean(this_condition_control_data, 2);
                
                % calculate response
                response_data_session{i_variable}(:, i_stretch) = data_session{i_variable}(:, i_stretch) - this_condition_control_mean;
            end
            
        end

    end
    
    %% integrate response
    integrated_data_session = cell(size(response_data_session));
    step_time_index_in_saved_data = find(strcmp(data_custodian.stretch_variable_names, 'step_time'), 1, 'first');
    this_time_data = data_session{step_time_index_in_saved_data};
    for i_variable = 1 : number_of_stretch_variables
        % find response data for this variable
        this_variable_name = variable_names_session{i_variable, 1};
        data_index_in_saved_data = find(strcmp(variable_names_session, this_variable_name), 1, 'first');
        this_response_data = response_data_session{data_index_in_saved_data};
        
        % integrate
        if size(this_response_data, 1) > 1
            this_response_data_integrated = zeros(size(this_response_data));
            for i_stretch = 1 : size(this_response_data, 2)
                step_time = this_time_data(i_stretch);
                time = linspace(0, step_time, study_settings.get('number_of_time_steps_normalized'));
                variable_data = this_response_data(:, i_stretch);
                this_response_data_integrated_this_stretch = cumtrapz(time, variable_data);
                this_response_data_integrated(:, i_stretch) = this_response_data_integrated_this_stretch;
            end
        else
            this_response_data_integrated = NaN;
        end
        % store integrated data
        integrated_data_session{data_index_in_saved_data} = this_response_data_integrated;
    end
    
    
    %% save data
    results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'results')];
    save ...
      ( ...
        results_file_name, ...
        'variable_data_session', ...
        'response_data_session', ...
        'integrated_data_session', ...
        'variable_names_session', ...
        'condition_stance_foot_list_session', ...
        'condition_perturbation_list_session', ...
        'condition_delay_list_session', ...
        'condition_index_list_session', ...
        'condition_experimental_list_session', ...
        'condition_stimulus_list_session', ...
        'condition_day_list_session', ...
        'origin_trial_list_session', ...
        'origin_start_time_list_session', ...
        'origin_end_time_list_session', ...
        'time_list_session' ...
      )
end

          

