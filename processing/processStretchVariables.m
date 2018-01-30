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

% process stretches

% input
% relevantDataStretches.mat

function processStretchVariables(varargin)
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
    
    % find relevant conditions
    conditions_settings = study_settings.get('conditions');
    condition_labels = conditions_settings(:, 1)';
    condition_source_variables = conditions_settings(:, 2)';
    number_of_condition_labels = length(condition_labels);
    
    % make containers to hold the data
    data_session = cell(number_of_stretch_variables, 1);
    conditions_session = struct;
    for i_condition = 1 : number_of_condition_labels
        conditions_session.(condition_source_variables{i_condition}) = {};
    end
    
    % make containers to store origin information for the stretches
    origin_trial_list_session = [];
    origin_start_time_list_session = [];
    origin_end_time_list_session = [];
    time_list_session = [];
    bands_per_stretch_session = [];
    
    %% analyze and store data
    for i_type = 1 : length(condition_list)
        condition = condition_list{i_type};
        trials_to_process = trial_number_list{i_type};
        for i_trial = trials_to_process
            % load and prepare data
            data_custodian.prepareBasicVariables(condition, i_trial);
            load(['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'relevantDataStretches')]);
            number_of_stretches_this_trial = size(stretch_times, 1);
            bands_per_stretch_session = [bands_per_stretch_session; bands_per_stretch];
%             data_trial = data_custodian.calculateStretchVariables(stretch_start_times, stretch_end_times, stretch_pushoff_times, time_normalization_markers, stance_foot_data, conditions_trial.condition_experimental_list);
            data_trial = data_custodian.calculateStretchVariables(stretch_times, stance_foot_data, conditions_trial.condition_experimental_list, stretch_pushoff_times);
            
            % append the data and condition lists from this trial to the total lists
            for i_variable = 1 : number_of_stretch_variables
                data_session{i_variable} = [data_session{i_variable} data_trial{i_variable}];
            end
            for i_condition = 1 : number_of_condition_labels
                conditions_session.(condition_source_variables{i_condition}) = [conditions_session.(condition_source_variables{i_condition}); conditions_trial.(condition_source_variables{i_condition}) ];
            end
            origin_trial_list_session = [origin_trial_list_session; ones(number_of_stretches_this_trial, 1) * i_trial]; %#ok<AGROW>
            origin_start_time_list_session = [origin_start_time_list_session; stretch_times(:, 1)]; %#ok<AGROW>
            origin_end_time_list_session = [origin_end_time_list_session; stretch_times(:, end)]; %#ok<AGROW>
            time_list = ones(number_of_stretches_this_trial, 1) * (i_trial - 1) * study_settings.get('trial_length') + stretch_times(:, 1);
            time_list_session = [time_list_session; time_list]; %#ok<AGROW>
            disp(['Processing stretch variables: condition ' condition_list{i_type} ', Trial ' num2str(i_trial) ' completed']);
        end
    end
    
    %% calculate some subject-level data and report
    number_of_stretches_session = size(data_session{1}, 2);
    
    % make condition data tables
    condition_data_all = cell(number_of_stretches_session, number_of_condition_labels);
    for i_condition = 1 : number_of_condition_labels
        condition_data_all(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
    end
    labels_to_ignore = study_settings.get('conditions_to_ignore');
    levels_to_remove = study_settings.get('levels_to_remove');
    [condition_combination_labels, condition_combinations_stimulus, condition_combinations_control] = determineConditionCombinations(condition_data_all, conditions_settings, labels_to_ignore, levels_to_remove);

    % extract indicators for control
    condition_combinations_control_unique = table2cell(unique(cell2table(condition_combinations_control), 'rows'));
    number_of_conditions_control = size(condition_combinations_control_unique, 1);
    conditions_control_indicators = true(number_of_stretches_session, number_of_conditions_control);
    for i_condition = 1 : number_of_conditions_control
        for i_label = 1 : length(condition_combination_labels)
            this_label = condition_combination_labels{i_label};
            this_label_list = condition_data_all(:, strcmp(conditions_settings(:, 1), this_label));
            this_label_indicator = strcmp(this_label_list, condition_combinations_control_unique(i_condition, i_label));
            conditions_control_indicators(:, i_condition) = conditions_control_indicators(:, i_condition) .* this_label_indicator;
        end
    end
    
    % report control
    trials_per_condition_control = sum(conditions_control_indicators)';
    conditions_control_with_number = condition_combinations_control_unique;
    for i_condition = 1 : number_of_conditions_control
        conditions_control_with_number{i_condition, size(condition_combinations_control_unique, 2)+1} = num2str(trials_per_condition_control(i_condition));
    end
    conditions_control_with_labels = [condition_combination_labels 'number of stretches'; conditions_control_with_number];
    disp('Control conditions:')
    disp(conditions_control_with_labels);
    
%     % extract indicators for control - old
%     conditions_control = study_settings.get('conditions_control');
%     number_of_conditions_control = size(study_settings.get('conditions_control'), 1);
%     conditions_control_indicators = false(number_of_stretches_session, number_of_conditions_control);
%     for i_condition = 1 : number_of_conditions_control
%         stance_foot_indicator = strcmp(condition_stance_foot_list_session, conditions_control(i_condition, 1));
%         perturbation_indicator = strcmp(condition_perturbation_list_session, conditions_control(i_condition, 2));
%         delay_indicator = strcmp(condition_delay_list_session, conditions_control(i_condition, 3));
%         index_indicator = strcmp(condition_index_list_session, conditions_control(i_condition, 4));
%         experimental_indicator = strcmp(condition_experimental_list_session, conditions_control(i_condition, 5));
%         stimulus_indicator = strcmp(condition_stimulus_list_session, conditions_control(i_condition, 6));
%         day_indicator = strcmp(condition_day_list_session, conditions_control(i_condition, 7));
% 
%         this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
%         conditions_control_indicators(:, i_condition) = this_condition_indicator;
%     end
    
    % extract indicators for stimulus
    number_of_conditions_stimulus = size(condition_combinations_stimulus, 1);
    conditions_stimulus_indicators = true(number_of_stretches_session, number_of_conditions_stimulus);
    for i_condition = 1 : number_of_conditions_stimulus
        for i_label = 1 : length(condition_combination_labels)
            this_label = condition_combination_labels{i_label};
            this_label_list = condition_data_all(:, strcmp(conditions_settings(:, 1), this_label));
            this_label_indicator = strcmp(this_label_list, condition_combinations_stimulus(i_condition, i_label));
            conditions_stimulus_indicators(:, i_condition) = conditions_stimulus_indicators(:, i_condition) .* this_label_indicator;
        end
    end

    % report stimulus
    trials_per_condition_stimulus = sum(conditions_stimulus_indicators)';
    conditions_stimulus_with_number = condition_combinations_stimulus;
    for i_condition = 1 : number_of_conditions_stimulus
        conditions_stimulus_with_number{i_condition, size(condition_combinations_stimulus, 2)+1} = num2str(trials_per_condition_stimulus(i_condition));
    end
    conditions_stimulus_with_labels = [condition_combination_labels 'number of stretches'; conditions_stimulus_with_number];
    disp('Stimulus conditions:')
    disp(conditions_stimulus_with_labels);
    
%     % extract indicators for conditions to analyze - old
%     conditions_to_analyze = study_settings.get('conditions_to_analyze');
%     number_of_conditions_to_analyze = size(conditions_to_analyze, 1);
%     conditions_to_analyze_indicators = false(number_of_stretches_session, number_of_conditions_to_analyze);
%     for i_condition = 1 : number_of_conditions_to_analyze
%         stance_foot_indicator = strcmp(condition_stance_foot_list_session, conditions_to_analyze(i_condition, 1));
%         if study_settings.get('analyze_total_response')
%             perturbation_indicator = strcmp(condition_perturbation_list_all,'ILLUSION_RIGHT') | strcmp(condition_perturbation_list_all,'ILLUSION_LEFT');
%         else
%             perturbation_indicator = strcmp(condition_perturbation_list_session, conditions_to_analyze(i_condition, 2)); 
%         end
%         delay_indicator = strcmp(condition_delay_list_session, conditions_to_analyze(i_condition, 3));
%         index_indicator = strcmp(condition_index_list_session, conditions_to_analyze(i_condition, 4));
%         experimental_indicator = strcmp(condition_experimental_list_session, conditions_to_analyze(i_condition, 5));
%         stimulus_indicator = strcmp(condition_stimulus_list_session, conditions_to_analyze(i_condition, 6));
%         day_indicator = strcmp(condition_day_list_session, conditions_to_analyze(i_condition, 7));
%         
%         this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
%         conditions_to_analyze_indicators(:, i_condition) = this_condition_indicator;
%     end
    
    % check the unassigned stretches
    assigned_stretch_indicator = sum([conditions_stimulus_indicators conditions_control_indicators], 2);
    unassigned_stretch_indicator = ~assigned_stretch_indicator;
    unassigned_stretch_indices = find(unassigned_stretch_indicator);
    if ~isempty(unassigned_stretch_indices)
        disp(['There were ' num2str(length(unassigned_stretch_indices)) ' unassigned stretches. Code to display more details needs to be updated (l. 206f).'])
    end
%     unassigned_stretch_labels_session = cell(length(unassigned_stretch_indices), length(study_settings.get('condition_labels')));
%     for i_stretch = 1 : length(unassigned_stretch_indices)
%         unassigned_stretch_labels_session(i_stretch, 1) = condition_stance_foot_list_session(unassigned_stretch_indices(i_stretch));
%         unassigned_stretch_labels_session(i_stretch, 2) = condition_perturbation_list_session(unassigned_stretch_indices(i_stretch));
%         unassigned_stretch_labels_session(i_stretch, 3) = condition_delay_list_session(unassigned_stretch_indices(i_stretch));
%         unassigned_stretch_labels_session(i_stretch, 4) = condition_index_list_session(unassigned_stretch_indices(i_stretch));
%         unassigned_stretch_labels_session(i_stretch, 5) = condition_experimental_list_session(unassigned_stretch_indices(i_stretch));
%         unassigned_stretch_labels_session(i_stretch, 6) = condition_stimulus_list_session(unassigned_stretch_indices(i_stretch));
%         unassigned_stretch_labels_session(i_stretch, 7) = condition_day_list_session(unassigned_stretch_indices(i_stretch));
%     end
%     wd = unassigned_stretch_labels_session;
%     [~, idx] = unique(strcat(wd(:,1), wd(:,2), wd(:,3), wd(:,4), wd(:,5), wd(:,6), wd(:,7)));
%     unassigned_stretch_labels = wd(idx,:);

%     % report control - old
%     trials_per_condition_control = sum(conditions_control_indicators)';
%     conditions_control_with_number = conditions_control;
%     for i_condition = 1 : number_of_conditions_control
%         conditions_control_with_number{i_condition, size(conditions_control, 2)+1} = num2str(trials_per_condition_control(i_condition));
%     end
%     conditions_control_with_labels = [study_settings.get('condition_labels') 'number of stretches'; conditions_control_with_number];
%     disp('Control conditions:')
%     disp(conditions_control_with_labels);
% 
%     % report conditions to analyze - old
%     trials_per_condition_to_analyze = sum(conditions_to_analyze_indicators)';
%     conditions_to_analyze_with_number = conditions_to_analyze;
%     for i_condition = 1 : number_of_conditions_to_analyze
%         conditions_to_analyze_with_number{i_condition, size(conditions_to_analyze, 2)+1} = num2str(trials_per_condition_to_analyze(i_condition));
%     end
%     conditions_to_analyze_with_labels = [study_settings.get('condition_labels') 'number of stretches'; conditions_to_analyze_with_number];
%     disp('Conditions to analyze:')
%     disp(conditions_to_analyze_with_labels);
    
%     disp('Stretches with these conditions were found but not analyzed:')
%     disp([study_settings.get('condition_labels'); unassigned_stretch_labels])
    
    disp(['Number of control stretches: ' num2str(sum(trials_per_condition_control))]);
    disp(['Number of stimulus stretches: ' num2str(sum(trials_per_condition_stimulus))]);
    disp(['Number of un-analyzed stretches: ' num2str(number_of_stretches_session - sum(trials_per_condition_control) - sum(trials_per_condition_stimulus))]);
    
    stretch_data_session = data_session; %#ok<NASGU>
    stretch_names_session = data_custodian.stretch_variable_names; %#ok<NASGU>
    
    bands_per_stretch = median(bands_per_stretch_session);
    if any(bands_per_stretch_session ~= bands_per_stretch)
        warning('Different trials have different numbers of bands per stretch')
    end
    
    %% save data
    results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'results')];
    save ...
      ( ...
        results_file_name, ...
        'conditions_session', ...
        'stretch_data_session', ...
        'stretch_names_session', ...
        'bands_per_stretch', ...
        'origin_trial_list_session', ...
        'origin_start_time_list_session', ...
        'origin_end_time_list_session', ...
        'time_list_session' ...
      )
end

          

