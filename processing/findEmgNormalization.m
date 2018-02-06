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

function findEmgNormalization(varargin)
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    
    
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
    subject_settings = SettingsCustodian('subjectSettings.txt');
    emg_variable_names = subject_settings.get('emg_labels')';
    data_custodian = WalkingDataCustodian(emg_variable_names);
    number_of_stretch_variables = length(emg_variable_names);
    
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
%     condition_stance_foot_list_subject = {};
%     condition_perturbation_list_subject = {};
%     condition_delay_list_subject = {};
%     condition_index_list_subject = {};
%     condition_experimental_list_subject = {};
%     condition_stimulus_list_subject = {};
%     condition_day_list_subject = {};
    
    
    % analyze and store data
    for i_type = 1 : length(condition_list)
        condition = condition_list{i_type};
        trials_to_process = trial_number_list{i_type};
        for i_trial = trials_to_process
            % TODO: check whether this trial contains any stretches relevant for EMG normalization
            disp(['Finding EMG normalization: condition ' condition_list{i_type} ', Trial ' num2str(i_trial) ' completed']);
            
            % load and prepare data
            data_custodian.prepareBasicVariables(condition, i_trial, [{'emg_trajectories'}; emg_variable_names]);
            
            load(['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'relevantDataStretches')]);
%             data_trial = data_custodian.calculateStretchVariables(stretch_start_times, stretch_end_times, stretch_pushoff_times, conditions_trial.condition_stance_foot_list, conditions_trial.condition_experimental_list, emg_variable_names);
            data_trial = data_custodian.calculateStretchVariables(stretch_times, stance_foot_data, conditions_trial.condition_experimental_list, stretch_pushoff_times, emg_variable_names);
            
            % append the data and condition lists from this trial to the total lists
            for i_variable = 1 : number_of_stretch_variables
                data_session{i_variable} = [data_session{i_variable} data_trial{i_variable}];
            end
            for i_condition = 1 : number_of_condition_labels
                conditions_session.(condition_source_variables{i_condition}) = [conditions_session.(condition_source_variables{i_condition}); conditions_trial.(condition_source_variables{i_condition}) ];
            end
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
    [condition_combination_labels, condition_combinations_stimulus, condition_combinations_control, condition_combinations_emg_unique] = determineConditionCombinations(condition_data_all, conditions_settings, labels_to_ignore, levels_to_remove);

    % extract indicators for emg
    number_of_conditions_emg = size(condition_combinations_emg_unique, 1);
    conditions_emg_indicators = true(number_of_stretches_session, number_of_conditions_emg);
    for i_condition = 1 : number_of_conditions_emg
        for i_label = 1 : length(condition_combination_labels)
            this_label = condition_combination_labels{i_label};
            this_label_list = condition_data_all(:, strcmp(conditions_settings(:, 1), this_label));
            this_label_indicator = strcmp(this_label_list, condition_combinations_emg_unique(i_condition, i_label));
            conditions_emg_indicators(:, i_condition) = conditions_emg_indicators(:, i_condition) .* this_label_indicator;
        end
    end    
    
    % report emg
    trials_per_condition_emg = sum(conditions_emg_indicators)';
    conditions_emg_with_number = condition_combinations_emg_unique;
    for i_condition = 1 : number_of_conditions_emg
        conditions_emg_with_number{i_condition, size(condition_combinations_emg_unique, 2)+1} = num2str(trials_per_condition_emg(i_condition));
    end
    conditions_emg_with_labels = [condition_combination_labels 'number of stretches'; conditions_emg_with_number];
    disp('EMG normalization conditions:')
    disp(conditions_emg_with_labels);
    
%     % calculate some subject-level data and report
%     number_of_stretches_subject = length(condition_stance_foot_list_subject);
%     
%     % extract indicators for emg normalization
%     conditions_emg_normalization = study_settings.get('conditions_emg_normalization');
%     number_of_conditions_emg_normalization = size(conditions_emg_normalization, 1);
%     conditions_emg_normalization_indicators = false(number_of_stretches_subject, number_of_conditions_emg_normalization);
%     for i_condition = 1 : number_of_conditions_emg_normalization
%         stance_foot_indicator = strcmp(condition_stance_foot_list_subject, conditions_emg_normalization(i_condition, 1));
%         perturbation_indicator = strcmp(condition_perturbation_list_subject, conditions_emg_normalization(i_condition, 2));
%         delay_indicator = strcmp(condition_delay_list_subject, conditions_emg_normalization(i_condition, 3));
%         index_indicator = strcmp(condition_index_list_subject, conditions_emg_normalization(i_condition, 4));
%         experimental_indicator = strcmp(condition_experimental_list_subject, conditions_emg_normalization(i_condition, 5));
%         stimulus_indicator = strcmp(condition_stimulus_list_subject, conditions_emg_normalization(i_condition, 6));
%         day_indicator = strcmp(condition_day_list_subject, conditions_emg_normalization(i_condition, 7));
% 
%         this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator & stimulus_indicator & day_indicator;
%         conditions_emg_normalization_indicators(:, i_condition) = this_condition_indicator;
%     end

%     % report conditions for emg normalization
%     trials_per_condition_emg_normalization = sum(conditions_emg_normalization_indicators)';
%     conditions_emg_normalization_with_number = conditions_emg_normalization;
%     for i_condition = 1 : number_of_conditions_emg_normalization
%         conditions_emg_normalization_with_number{i_condition, size(conditions_emg_normalization, 2)+1} = num2str(trials_per_condition_emg_normalization(i_condition));
%     end
%     conditions_emg_normalization_with_labels = [study_settings.get('condition_labels') 'number of stretches'; conditions_emg_normalization_with_number];
%     disp('EMG normalization conditions:')
%     disp(conditions_emg_normalization_with_labels);

    % average across stretches
    emg_normalization_values = zeros(length(emg_variable_names), 1) * NaN;
    for i_variable = 1 : number_of_stretch_variables
        data_this_variable = data_session{i_variable};
        condition_averages_this_variable = zeros(1, number_of_conditions_emg) * NaN;
        
        if ~(strcmp(emg_variable_names{i_variable}, 'empty') || strcmp(emg_variable_names{i_variable}, '~'))
            for i_condition = 1 : number_of_conditions_emg
                this_condition_indicator = conditions_emg_indicators(:, i_condition);
                data_this_condition = data_this_variable(:, this_condition_indicator);
                average_this_condition = median(data_this_condition, 2); % average across trials
                condition_averages_this_variable(i_condition) = mean(average_this_condition); % average across time

                if visualize
                    figure; hold on;
                    title(emg_variable_names{i_variable});
                    plot(data_this_condition);
                    plot(average_this_condition, 'linewidth', 5);
                end
            end
            emg_normalization_values(i_variable) = mean(condition_averages_this_variable);
        end
    end
    
    % save data
    results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'emgNormalization')];
    save ...
      ( ...
        results_file_name, ...
        'emg_normalization_values', ...
        'emg_variable_names' ...
      )
end

          

