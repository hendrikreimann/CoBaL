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

function processEmgNormalization(varargin)
    %% prepare
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    
    load('subjectInfo.mat', 'date', 'subject_id');
    % load settings   
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    
    emg_import_map_header = subject_settings.get('emg_import_map_header');
    emg_import_map = subject_settings.get('emg_import_map');
    emg_labels = emg_import_map(:, strcmp(emg_import_map_header, 'label_in_cobal'));
    emg_labels_normalization = emg_labels;
    for i_channel = 1 : length(emg_labels_normalization)
        emg_labels_normalization{i_channel} = ['emg:' emg_labels_normalization{i_channel}];
    end

    %% collect data for normalization
    data_custodian = WalkingDataCustodian(emg_labels_normalization);
    number_of_stretch_variables = length(data_custodian.stretch_variable_names);
    
    % find relevant conditions
    conditions_settings = study_settings.get('conditions');
    condition_labels = conditions_settings(:, 1)';
    condition_source_variables = conditions_settings(:, 2)';
    number_of_condition_labels = length(condition_labels);
   
    condition_relevant_for_analysis = study_settings.get('condition_relevant_for_analysis');
    condition_relevant_name = conditions_settings{strcmp(conditions_settings(:, 1), condition_relevant_for_analysis), 2};
    
    % make containers to hold the data
    data_session = cell(number_of_stretch_variables, 1);
    conditions_session = struct;
    for i_condition = 1 : number_of_condition_labels
        conditions_session.(condition_source_variables{i_condition}) = {};
    end
    
    % analyze and store data
    for i_type = 1 : length(condition_list)
        this_type = condition_list{i_type};
        trials_to_process = trial_number_list{i_type};
        for i_trial = trials_to_process
            % TODO: check whether this trial contains any stretches relevant for EMG normalization
            disp(['Finding EMG normalization: condition ' condition_list{i_type} ', Trial ' num2str(i_trial) ' completed']);
            
            % load and prepare data
            data_custodian.prepareBasicVariables(this_type, i_trial, [{'emg_trajectories'}; emg_labels_normalization]);
            
            load(['analysis' filesep makeFileName(date, subject_id, this_type, i_trial, 'relevantDataStretches')], 'stretch_times', 'stance_foot_data', 'conditions_trial');
            data_trial = data_custodian.calculateStretchVariables(stretch_times, stance_foot_data, conditions_trial.(condition_relevant_name), emg_labels_normalization);
            
            % append the data and condition lists from this trial to the total lists
            for i_variable = 1 : number_of_stretch_variables
                data_session{i_variable} = [data_session{i_variable} data_trial{i_variable}];
            end
            for i_condition = 1 : number_of_condition_labels
                conditions_session.(condition_source_variables{i_condition}) = [conditions_session.(condition_source_variables{i_condition}); conditions_trial.(condition_source_variables{i_condition}) ];
            end
        end
    end
    
    %% calculate normalization values
    number_of_stretches_session = size(data_session{1}, 2);
    variable_names = data_custodian.stretch_variable_names;
    
    % make condition data tables
    condition_data_all = cell(number_of_stretches_session, number_of_condition_labels);
    for i_condition = 1 : number_of_condition_labels
        condition_data_all(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
    end
    labels_to_ignore = study_settings.get('conditions_to_ignore');
    levels_to_remove = study_settings.get('levels_to_remove', 1);
    [condition_combination_labels, ~, condition_combinations_emg] = determineConditionCombinations(condition_data_all, conditions_settings, labels_to_ignore, levels_to_remove);

    % extract indicators for emg
    number_of_conditions_emg = size(condition_combinations_emg, 1);
    conditions_emg_indicators = true(number_of_stretches_session, number_of_conditions_emg);
    for i_condition = 1 : number_of_conditions_emg
        for i_label = 1 : length(condition_combination_labels)
            this_label = condition_combination_labels{i_label};
            this_label_list = condition_data_all(:, strcmp(conditions_settings(:, 1), this_label));
            this_label_indicator = strcmp(this_label_list, condition_combinations_emg(i_condition, i_label));
            conditions_emg_indicators(:, i_condition) = conditions_emg_indicators(:, i_condition) .* this_label_indicator;
        end
    end    
    
    % report emg
    trials_per_condition_emg = sum(conditions_emg_indicators)';
    conditions_emg_with_number = condition_combinations_emg;
    for i_condition = 1 : number_of_conditions_emg
        conditions_emg_with_number{i_condition, size(condition_combinations_emg, 2)+1} = num2str(trials_per_condition_emg(i_condition));
    end
    conditions_emg_with_labels = [condition_combination_labels 'number of stretches'; conditions_emg_with_number];
    disp('EMG normalization conditions:')
    disp(conditions_emg_with_labels);
    
    % average across stretches
    emg_normalization_values = zeros(length(variable_names), 1) * NaN;
    for i_variable = 1 : number_of_stretch_variables
        data_this_variable = data_session{i_variable};
        condition_averages_this_variable = zeros(1, number_of_conditions_emg) * NaN;
        
        for i_condition = 1 : number_of_conditions_emg
            this_condition_indicator = conditions_emg_indicators(:, i_condition);
            data_this_condition = data_this_variable(:, this_condition_indicator);
            average_this_condition = nanmedian(data_this_condition, 2); % average across trials
            condition_averages_this_variable(i_condition) = mean(average_this_condition); % average across time

            if visualize
                figure; hold on;
                title(variable_names{i_variable});
                plot(data_this_condition);
                plot(average_this_condition, 'linewidth', 5);
            end
        end
        emg_normalization_values(i_variable) = mean(condition_averages_this_variable);
    end
    
    %% scale EMG data
    data_to_remove_header = subject_settings.get('data_to_remove_header', 1);
    data_to_remove = subject_settings.get('data_to_remove', 1);
    for i_type = 1 : length(condition_list)
        this_type = condition_list{i_type};
        trials_to_process = trial_number_list{i_type};
        this_trial_type_to_remove_rows = strcmp(data_to_remove(:, strcmp(data_to_remove_header, 'trial_type')), this_type);
        for i_trial = trials_to_process
            this_trial_number_to_remove_rows = strcmp(data_to_remove(:, strcmp(data_to_remove_header, 'trial_number')), num2str(i_trial));
            if isempty(data_to_remove)
                data_to_remove_this_trial = [];
            else
                data_to_remove_this_trial = data_to_remove(this_trial_number_to_remove_rows & this_trial_type_to_remove_rows, 1);
            end
            
            % load
            [emg_trajectories, time_emg, sampling_rate_emg, emg_labels_trial, emg_directions] = loadData(date, subject_id, this_type, i_trial, 'emg_trajectories');
            
            % scale
            emg_scaled_trajectories = zeros(size(emg_trajectories)) * NaN;
            for i_channel = 1 : size(emg_scaled_trajectories, 2)
                if i_channel <= length(emg_labels_trial) && any(strcmp(emg_labels_normalization, ['emg:' emg_labels_trial{i_channel}]))
                    this_label = emg_labels_trial{i_channel};
                    this_channel_weight = 1 / emg_normalization_values(strcmp(emg_labels_normalization, ['emg:' emg_labels_trial{i_channel}]));
                else
                    this_label = '';
                    this_channel_weight = NaN;
                end
                if any(strcmp(data_to_remove_this_trial, ['emg:' this_label]))
                    this_channel_weight = NaN;
                end
                emg_scaled_trajectories(:, i_channel) = emg_trajectories(:, i_channel) * this_channel_weight;
            end
            emg_labels = emg_labels_trial;
            
            % save
            save_folder = 'processed';
            emg_normalization_labels = variable_names;
            save_file_name = makeFileName(date, subject_id, this_type, i_trial, 'emgScaledTrajectories.mat');
            save ...
              ( ...
                [save_folder filesep save_file_name], ...
                'emg_scaled_trajectories', ...
                'time_emg', ...
                'sampling_rate_emg', ...
                'emg_labels', ...
                'emg_directions', ...
                'emg_normalization_values', ...
                'emg_normalization_labels' ...
              );
            addAvailableData ...
              ( ...
                'emg_scaled_trajectories', ...
                'time_emg', ...
                'sampling_rate_emg', ...
                '_emg_labels', ...
                '_emg_directions', ...
                save_folder, ...
                save_file_name ...
              );
            
            disp(['normalized EMG data for trial ' num2str(i_trial) ' of type ' this_type ' and saved as ' save_file_name])
        end
    end
    
    
end

          

