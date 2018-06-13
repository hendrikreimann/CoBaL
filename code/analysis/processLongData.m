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

function processLongData(varargin)
    % parse arguments
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});

    %% prepare
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
    number_of_time_steps_normalized = study_settings.get('number_of_time_steps_normalized');

    % create containers
    long_stretch_times_session = [];
    long_stretch_control_times_session = [];
    long_stretch_data_labels_session = {'mpsis_x', 'mpsis_x_vel'};
    long_stretch_conditions_session = struct;
    long_stretch_control_conditions_session = struct;
    number_of_variables = length(long_stretch_data_labels_session);
    long_stretch_data_session = cell(number_of_variables, 1);
    long_stretch_control_data_session = cell(number_of_variables, 1);
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            %% load data
            load(['analysis' filesep makeFileName(date, subject_id, condition_list{i_condition}, i_trial, 'events')]);
            load(['analysis' filesep makeFileName(date, subject_id, condition_list{i_condition}, i_trial, 'relevantDataStretches')]);
            
            % create condition container
            long_stretch_conditions_trial = struct;
            long_stretch_control_conditions_trial = struct;
            condition_labels = fieldnames(conditions_trial);
            for i_field = 1 : numel(condition_labels)
                long_stretch_conditions_trial.(condition_labels{i_field}) = {};
                long_stretch_control_conditions_trial.(condition_labels{i_field}) = {};
            end
            if isempty(fieldnames(long_stretch_conditions_session))
                for i_field = 1 : numel(condition_labels)
                    long_stretch_conditions_session.(condition_labels{i_field}) = {};
                    long_stretch_control_conditions_session.(condition_labels{i_field}) = {};
                end
            end            
            
            %% extract times
            number_of_strides = study_settings.get('strides_to_collect_long');
            long_stretch_times_trial = [];
            long_stretch_control_times_trial = [];
            for i_stretch = 1 : size(stretch_times, 1)
                if strcmp(conditions_trial.stimulus_list{i_stretch}, 'STIM_NONE')
                    % find heelstrike that started this stretch
                    if strcmp(conditions_trial.trigger_foot_list{i_stretch}, 'TRIGGER_RIGHT')
                        trigger_foot_heelstrikes = event_data{strcmp(event_labels, 'right_touchdown')};
                        contra_foot_heelstrikes = event_data{strcmp(event_labels, 'left_touchdown')};
                    elseif strcmp(conditions_trial.trigger_foot_list{i_stretch}, 'TRIGGER_LEFT')
                        trigger_foot_heelstrikes = event_data{strcmp(event_labels, 'left_touchdown')};
                        contra_foot_heelstrikes = event_data{strcmp(event_labels, 'right_touchdown')};
                    else
                        error('Trigger should be "TRIGGER_LEFT" or "TRIGGER_RIGHT"');
                    end

                    % find first indices
                    trigger_first_heelstrike_time = stretch_times(i_stretch, 1);
                    [~, trigger_first_heelstrike_index] = min(abs(trigger_foot_heelstrikes - trigger_first_heelstrike_time));
                    contra_first_heelstrike_time = min(contra_foot_heelstrikes(contra_foot_heelstrikes >= trigger_first_heelstrike_time));
                    [~, contra_first_heelstrike_index] = min(abs(contra_foot_heelstrikes - contra_first_heelstrike_time));

                    % find following indices
                    trigger_foot_indices = trigger_first_heelstrike_index : trigger_first_heelstrike_index + 1;
                    contra_foot_indices = contra_first_heelstrike_index;
                    trigger_foot_times = trigger_foot_heelstrikes(trigger_foot_indices);
                    contra_foot_times = contra_foot_heelstrikes(contra_foot_indices);

                    % find times
                    this_long_stretch_times = zeros(3, 1);
                    this_long_stretch_times([1 3]) = trigger_foot_times;
                    this_long_stretch_times(2) = contra_foot_times;
                    if ~issorted(this_long_stretch_times)
                        error('Something went wrong with determining the stretch times here');
                    end
                    long_stretch_control_times_trial = [long_stretch_control_times_trial this_long_stretch_times];
                    
                    % copy conditions
                    for i_field = 1 : numel(condition_labels)
                        long_stretch_control_conditions_trial.(condition_labels{i_field}) = ...
                          [ ...
                            long_stretch_control_conditions_trial.(condition_labels{i_field}); ...
                            conditions_trial.(condition_labels{i_field}){i_stretch} ...
                          ];
                    end
                end
                
                if ~strcmp(conditions_trial.stimulus_list{i_stretch}, 'STIM_NONE')
                    % find heelstrike that started this stretch
                    if strcmp(conditions_trial.trigger_foot_list{i_stretch}, 'TRIGGER_RIGHT')
                        trigger_foot_heelstrikes = event_data{strcmp(event_labels, 'right_touchdown')};
                        contra_foot_heelstrikes = event_data{strcmp(event_labels, 'left_touchdown')};
                    elseif strcmp(conditions_trial.trigger_foot_list{i_stretch}, 'TRIGGER_LEFT')
                        trigger_foot_heelstrikes = event_data{strcmp(event_labels, 'left_touchdown')};
                        contra_foot_heelstrikes = event_data{strcmp(event_labels, 'right_touchdown')};
                    else
                        error('Trigger should be "TRIGGER_LEFT" or "TRIGGER_RIGHT"');
                    end

                    % find first indices
                    trigger_first_heelstrike_time = stretch_times(i_stretch, 1);
                    [~, trigger_first_heelstrike_index] = min(abs(trigger_foot_heelstrikes - trigger_first_heelstrike_time));
                    contra_first_heelstrike_time = min(contra_foot_heelstrikes(contra_foot_heelstrikes >= trigger_first_heelstrike_time));
                    [~, contra_first_heelstrike_index] = min(abs(contra_foot_heelstrikes - contra_first_heelstrike_time));

                    % find following indices
                    trigger_foot_indices = trigger_first_heelstrike_index : trigger_first_heelstrike_index + number_of_strides;
                    contra_foot_indices = contra_first_heelstrike_index : contra_first_heelstrike_index + number_of_strides - 1;
                    if length(trigger_foot_heelstrikes) >= trigger_foot_indices(end) && length(trigger_foot_heelstrikes) >= contra_foot_indices(end)
                    
                        trigger_foot_times = trigger_foot_heelstrikes(trigger_foot_indices);
                        contra_foot_times = contra_foot_heelstrikes(contra_foot_indices);

                        % find times
                        this_long_stretch_times = zeros(number_of_strides*2+1, 1);
                        this_long_stretch_times(1 : 2 : number_of_strides*2+1) = trigger_foot_times;
                        this_long_stretch_times(2 : 2 : number_of_strides*2) = contra_foot_times;
                        if ~issorted(this_long_stretch_times)
                            error('Something went wrong with determining the stretch times here');
                        end
                        long_stretch_times_trial = [long_stretch_times_trial this_long_stretch_times];

                        % copy conditions
                        for i_field = 1 : numel(condition_labels)
                            long_stretch_conditions_trial.(condition_labels{i_field}) = ...
                              [ ...
                                long_stretch_conditions_trial.(condition_labels{i_field}); ...
                                conditions_trial.(condition_labels{i_field}){i_stretch} ...
                              ];
                        end
                    end
                end
                
            end
            
            %% extract data
            [lpsi_x_data_base, time, sampling_rate, labels, directions, success] = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'marker_trajectories:LPSI_x'); %#ok<ASGLU>
            [rpsi_x_data_base, time, sampling_rate, labels, directions, success] = loadData(date, subject_id, condition_list{i_condition}, i_trial, 'marker_trajectories:RPSI_x'); %#ok<ASGLU>
            mpsis_x_data_base = (lpsi_x_data_base + rpsi_x_data_base) * 0.5;
            mpsis_x_vel_data_base = deriveByTime(mpsis_x_data_base, time);
            
            
%             if study_settings.get('filter_marker_data')
%                 filter_order = study_settings.get('filter_order_com_vel');
%                 cutoff_frequency = study_settings.get('filter_cutoff_com_vel'); % in Hz
%                 [b_marker, a_marker] = butter(filter_order, cutoff_frequency/(sampling_rate/2));
%                 
%                 mpsis_x_vel_data_filtered = nanfiltfilt(b_marker, a_marker, mpsis_x_vel_data_base);
% 
%                 figure; hold on
%                 plot(mpsis_x_vel_data_base)
%                 plot(mpsis_x_vel_data_filtered);
%             end
            
            
            number_of_long_stretches_trial = size(long_stretch_times_trial, 2);
            long_stretch_data_trial = cell(number_of_variables, 1);
            long_stretch_data_trial{1} = zeros(number_of_strides * 2 * (number_of_time_steps_normalized-1) + 1, number_of_long_stretches_trial);
            long_stretch_data_trial{2} = zeros(number_of_strides * 2 * (number_of_time_steps_normalized-1) + 1, number_of_long_stretches_trial);
            for i_stretch = 1 : number_of_long_stretches_trial
                this_stretch_times = long_stretch_times_trial(:, i_stretch);
                long_stretch_data_trial{1}(:, i_stretch) = getTimeNormalizedData(time, mpsis_x_data_base, this_stretch_times, number_of_time_steps_normalized);
                long_stretch_data_trial{2}(:, i_stretch) = getTimeNormalizedData(time, mpsis_x_vel_data_base, this_stretch_times, number_of_time_steps_normalized);
            end
            
            number_of_long_control_stretches_trial = size(long_stretch_control_times_trial, 2);
            long_stretch_control_data_trial = cell(number_of_variables, 1);
            long_stretch_control_data_trial{1} = zeros(2 * (number_of_time_steps_normalized-1) + 1, number_of_long_control_stretches_trial);
            long_stretch_control_data_trial{2} = zeros(2 * (number_of_time_steps_normalized-1) + 1, number_of_long_control_stretches_trial);
            for i_stretch = 1 : number_of_long_control_stretches_trial
                this_stretch_times = long_stretch_control_times_trial(:, i_stretch);
                long_stretch_control_data_trial{1}(:, i_stretch) = getTimeNormalizedData(time, mpsis_x_data_base, this_stretch_times, number_of_time_steps_normalized);
                long_stretch_control_data_trial{2}(:, i_stretch) = getTimeNormalizedData(time, mpsis_x_vel_data_base, this_stretch_times, number_of_time_steps_normalized);
            end
            
            %% store
            long_stretch_times_session = [long_stretch_times_session long_stretch_times_trial];
            long_stretch_control_times_session = [long_stretch_control_times_session long_stretch_control_times_trial];
            for i_variable = 1 : number_of_variables
                long_stretch_data_session{i_variable} = [long_stretch_data_session{i_variable}  long_stretch_data_trial{i_variable}];
                long_stretch_control_data_session{i_variable} = [long_stretch_control_data_session{i_variable} long_stretch_control_data_trial{i_variable}];
            end
            for i_field = 1 : numel(condition_labels)
                long_stretch_conditions_session.(condition_labels{i_field}) = ...
                  [ ...
                    long_stretch_conditions_session.(condition_labels{i_field}), ...
                    long_stretch_conditions_trial.(condition_labels{i_field})' ...
                  ];
                long_stretch_control_conditions_session.(condition_labels{i_field}) = ...
                  [ ...
                    long_stretch_control_conditions_session.(condition_labels{i_field}), ...
                    long_stretch_control_conditions_trial.(condition_labels{i_field})' ...
                  ];
            end
            
            disp(['Processing long-term data: condition ' condition_list{i_condition} ', Trial ' num2str(i_trial) ' completed']);                
        end
    end
    step_time_data_session = diff(long_stretch_times_session);
    
    %% process response data
    
    % create container
    long_stretch_response_data_session = cell(size(long_stretch_data_session));
    for i_variable = 1 : length(long_stretch_data_labels_session)
        long_stretch_response_data_session{i_variable} = zeros(size(long_stretch_data_session{i_variable}));
    end
    
    bands_per_long_stretch = size(long_stretch_times_session, 1) - 1; %#ok<NASGU>
    conditions_settings = study_settings.get('conditions');
    condition_labels = conditions_settings(:, 1)';
    number_of_condition_labels = length(condition_labels);
    condition_source_variables = conditions_settings(:, 2)';
    number_of_stretches_this_session = size(long_stretch_times_session, 2);
    number_of_stretches_control_this_session = size(long_stretch_control_times_session, 2);
    
    % transform conditions into cell array
    conditions_session = long_stretch_conditions_session;
    condition_array_session = cell(number_of_stretches_this_session, number_of_condition_labels);
    conditions_control_session = long_stretch_control_conditions_session;
    condition_control_array_session = cell(number_of_stretches_control_this_session, number_of_condition_labels);
    for i_condition = 1 : number_of_condition_labels
        condition_array_session(:, i_condition) = conditions_session.(condition_source_variables{i_condition});
        condition_control_array_session(:, i_condition) = conditions_control_session.(condition_source_variables{i_condition});
    end
    
    % extract data
    labels_to_ignore = study_settings.get('conditions_to_ignore');
    levels_to_remove = study_settings.get('levels_to_remove');
    [condition_combination_labels, ~, condition_combinations_control] = determineConditionCombinations(condition_array_session, conditions_settings, labels_to_ignore, levels_to_remove);
    condition_combinations_control_unique = table2cell(unique(cell2table(condition_combinations_control), 'rows'));

    % go stretch by stretch
    for i_stretch = 1 : number_of_stretches_this_session
        % extract this stretches relevant conditions
        this_stretch_condition_string = cell(1, length(condition_combination_labels));
        for i_label = 1 : length(condition_combination_labels)
            this_label = condition_combination_labels{i_label};
            this_label_source_variable = condition_source_variables{strcmp(condition_labels, this_label)};
            this_label_condition_list = conditions_session.(this_label_source_variable);
            this_stretch_condition_string{i_label} = this_label_condition_list{i_stretch};
        end

        % determine applicable control condition index
        if strcmp(study_settings.get('experimental_paradigm'), 'GVS_old')
            if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_LEFT')
                applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_LEFT'));
            end
            if strcmp(this_stretch_condition_string{strcmp(condition_combination_labels, 'trigger_foot')}, 'TRIGGER_RIGHT')
                applicable_control_condition_index = find(strcmp(condition_combinations_control_unique(:, strcmp(condition_combination_labels, 'trigger_foot')), 'TRIGGER_RIGHT'));
            end
        end

        % determine indicator for control
        control_condition_indicator = true(number_of_stretches_control_this_session, 1);
        for i_label = 1 : length(condition_combination_labels)
            this_label = condition_combination_labels{i_label};
            this_label_list = condition_control_array_session(:, strcmp(conditions_settings(:, 1), this_label));
            this_label_indicator = strcmp(this_label_list, condition_combinations_control_unique(applicable_control_condition_index, i_label));
            control_condition_indicator = control_condition_indicator .* this_label_indicator;
        end        
        control_condition_indicator = logical(control_condition_indicator);            

        % calculate responses
        for i_variable = 1 : number_of_variables
            % calculate control mean
            this_variable_data = long_stretch_data_session{i_variable}(:, i_stretch);
            this_variable_control_data = long_stretch_control_data_session{i_variable}(:, control_condition_indicator);
            this_variable_control_mean = mean(this_variable_control_data, 2);
            this_variable_control_mean_long = this_variable_control_mean;
            this_variable_control_mean_circle_glue = (this_variable_control_mean_long(1) + this_variable_control_mean_long(end)) / 2;
            this_variable_control_mean(1) = this_variable_control_mean_circle_glue;
            this_variable_control_mean(end) = this_variable_control_mean_circle_glue;
            for i_band = 2 : (bands_per_long_stretch / 2)
                this_variable_control_mean_long = [this_variable_control_mean_long; this_variable_control_mean(2:end)]; %#ok<AGROW>
            end
            
            % smooth connection points
            this_variable_control_mean_long_filtered = this_variable_control_mean_long;
            this_variable_control_mean_long_replaced = this_variable_control_mean_long;

            filter_order = 6;
            cutoff_frequency = 10; % in Hz
            sampling_rate = (0.6 / 100)^(-1);
            [b_filter, a_filter] = butter(filter_order, cutoff_frequency/(sampling_rate/2), 'low');
            smooth_range_width = 10;
            sinusoid_original = (cos(linspace(0, 2*pi, smooth_range_width*2+1)') + 1) * 0.5;
            sinusoid_filtered = 1 - sinusoid_original;
            for i_band = 2 : 2 : bands_per_long_stretch-1
                merge_index = (number_of_time_steps_normalized-1) * i_band + 1;
                merge_range = merge_index-smooth_range_width : merge_index+smooth_range_width;

%                 this_variable_control_mean_range = this_variable_control_mean_long(merge_range);
%                 this_variable_control_mean_range_filtered = filtfilt(b_filter, a_filter, this_variable_control_mean_range);
%                 this_variable_control_mean_range_replace = this_variable_control_mean_range .* sinusoid_original + this_variable_control_mean_range_filtered .* sinusoid_filtered;
%                 this_variable_control_mean_long_filtered(merge_range) = this_variable_control_mean_range_filtered;
%                 this_variable_control_mean_long_replaced(merge_range) = this_variable_control_mean_range_replace;
                
                this_variable_control_mean_long_gap = this_variable_control_mean_long;
                this_variable_control_mean_long_gap(merge_range) = NaN;
                warning('off', 'MATLAB:chckxy:IgnoreNaN')
                this_variable_control_mean_long_splined = spline((1:length(this_variable_control_mean_long_gap))', this_variable_control_mean_long_gap, (1:length(this_variable_control_mean_long_gap))');
                warning('on', 'MATLAB:chckxy:IgnoreNaN')

                this_variable_control_mean_long_replaced(merge_range) = this_variable_control_mean_long_splined(merge_range);

            end

%             figure; hold on
%             plot(this_variable_control_mean_long)
%             plot(this_variable_control_mean_long_replaced)
            
            % calculate maximal splining error
            maximal_smoothing_error = max(abs(this_variable_control_mean_long_replaced - this_variable_control_mean_long));
%             disp(['maximal smoothing error: ' num2str(maximal_smoothing_error)])
            
            % store
            long_stretch_response_data_session{i_variable}(:, i_stretch) = this_variable_data - this_variable_control_mean_long_replaced;
            
        end

    end

%     % plot some things
%     figure; hold on;
%     plot(this_variable_control_data)
%     plot(this_variable_control_mean, 'linewidth', 3)
%     
%     
%     % plot some things
%     figure; hold on;
%     plot(long_stretch_data_session{1})
%     plot(this_variable_control_mean_long_replaced, 'linewidth', 3)
    
    %% save
    results_file_name = ['analysis' filesep makeFileName(date, subject_id, 'longStretchResults')];
    save ...
      ( ...
        results_file_name, ...
        'step_time_data_session', ...
        'long_stretch_times_session', ...
        'long_stretch_data_session', ...
        'long_stretch_response_data_session', ...
        'long_stretch_conditions_session', ...
        'long_stretch_control_times_session', ...
        'long_stretch_control_data_session', ...
        'long_stretch_control_conditions_session', ...
        'long_stretch_data_labels_session', ...
        'bands_per_long_stretch' ...
      )
    
end

function data_normalized = getTimeNormalizedData(variable_time, variable_data, stretch_times, number_of_time_steps_normalized)
    % extract data
    stretch_time_indices = zeros(size(stretch_times));
    for i_stretch_time = 1 : length(stretch_times)
        [~, time_index] = min(abs(variable_time - stretch_times(i_stretch_time)));
        stretch_time_indices(i_stretch_time) = time_index;
    end

    time_extracted = variable_time(stretch_time_indices(1) : stretch_time_indices(end));
    data_extracted = variable_data(stretch_time_indices(1) : stretch_time_indices(end));
    stretch_time_indices_local = stretch_time_indices - stretch_time_indices(1) + 1;

    % normalize data in time
    if ~isempty(time_extracted) && ~any(isnan(data_extracted))
        % create normalized time
        number_of_stretches = length(stretch_time_indices_local) - 1;
        time_normalized = [];
        for i_stretch = 1 : number_of_stretches
            time_normalized_this_stretch = linspace(time_extracted(stretch_time_indices_local(i_stretch)), time_extracted(stretch_time_indices_local(i_stretch+1)), number_of_time_steps_normalized)';
            if i_stretch > 1
                % start time of this stretch is end time of the last stretch, so remove the duplicate point
                time_normalized_this_stretch = time_normalized_this_stretch(2:end);
            end
            time_normalized = [time_normalized; time_normalized_this_stretch]; %#ok<AGROW>
        end

        % time-normalize data
        data_normalized = spline(time_extracted, data_extracted, time_normalized);
    else
        number_of_stretches = length(stretch_time_indices_local) - 1;
        time_normalized = [];
        for i_stretch = 1 : number_of_stretches
            time_normalized_this_stretch = linspace(time_extracted(stretch_time_indices_local(i_stretch)), time_extracted(stretch_time_indices_local(i_stretch+1)), this.number_of_time_steps_normalized)';
            if i_stretch > 1
                % start time of this stretch is end time of the last stretch, so remove the duplicate point
                time_normalized_this_stretch = time_normalized_this_stretch(2:end);
            end
            time_normalized = [time_normalized; time_normalized_this_stretch]; %#ok<AGROW>
        end

        % time-normalize data
        data_normalized = time_normalized * NaN;
    end
end
