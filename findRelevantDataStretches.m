% This script finds the stretches of data relevant to the experimental paradigm.
function findRelevantDataStretches(varargin)
    % parse arguments
    parser = inputParser;
    parser.KeepUnmatched = true;
    
    default_stimulus_type = 'none';
    valid_stimulus_types = {'none', 'gvs', 'visual'};
    check_stimulus_type = @(x) any(validatestring(x,valid_stimulus_types));
    addParameter(parser, 'stimulus', default_stimulus_type, check_stimulus_type)
    
    parse(parser, varargin{:})
    stimulus_type = parser.Results.stimulus;
    
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});

    %% Choose Data Type
    emg_present                     = 1;
%     stimulus_type = 'gvs';
    % stimulus_type = 'visual';
%     stimulus_type = 'none';

    remove_crossover_steps          = 0;
    visualize                       = 0;

    
    % stretches_to_analyze = 'single stance';
    stretches_to_analyze = 'full stride';

    %% prepare
    wait_times = [0 0.150 0.450];
    wait_time_labels = {'0ms', '150ms', '450ms'};

%     trials_to_process = 1 : 20;
    % trials_to_process = 5;

%     total_positive_steps = [];
%     total_negative_steps = [];
%     total_perturbation_condition_list = [];
%     total_step_condition_list = [];

%     swing_foot_fz_zero_threshold = 20; % threshold for counting a vertical force reading as zero, in Nm
    swing_foot_nonzero_stretch_length_threshold_time = 0.05; % the force plate should register non-zero for at less than this long
    time_to_nearest_heelstrike_before_trigger_threshold = 0.10; % a heelstrike should happen less than this long after a trigger
    time_to_nearest_heelstrike_after_trigger_threshold = 0.3; % a heelstrike should happen less than this long after a trigger
    
    
    duration_prelude = 0.4; % duration that should have zero vertical force prior to a heelstrike, to identify crossover

    colors = parula(6);


    %% extract data
    load('subjectInfo.mat', 'date', 'subject_id');
    stretch_length_indices_forceplate = [];
    condition_stance_foot_list = {};
    condition_perturbation_list = {};
    condition_delay_list = {};
    condition_index_list = {};
    stim_start_time_relative_to_stretch = [];

    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
        
            %% load data
            condition = condition_list{i_condition};

            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
            
            % XXX this should be replaced with a system where the subject info file stores the available data for each trial
            labview_file_name = ['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'labviewData')];
            if exist([labview_file_name '.mat'], 'file')
                load(labview_file_name);
            end                
            forceplate_file_name = ['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'forceplateTrajectories')];
            if exist([forceplate_file_name '.mat'], 'file')
                load(forceplate_file_name);
                if isnan(sampling_rate_forceplate)
                    sampling_rate_forceplate = 1 / median(diff(time_forceplate));
                end
                number_of_prelude_time_steps = ceil(duration_prelude * sampling_rate_forceplate);
                swing_foot_nonzero_stretch_length_threshold_indices = ceil(swing_foot_nonzero_stretch_length_threshold_time * sampling_rate_forceplate);
                force_plate_data_available = 1;
            else
                force_plate_data_available = 0;
            end
            load(['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'stepEvents')]);
            if emg_present
                load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'emgTrajectories')]);
            end

            % determine illusion
            if ~strcmp(stimulus_type, 'none')
                illusion_trajectory = zeros(size(time_labviewData)); % 1 = RIGHT, -1 = LEFT
                if strcmp(stimulus_type, 'gvs')
                    % use GVS_out_trajectory
                    for i_time = 1 : length(time_labviewData)
                        if GVS_out_trajectory(i_time) > 0
                            % anode is on the right, cathode is on the left, illusory fall is towards the cathode, i.e. LEFT
                            illusion_trajectory(i_time) = -1;
                        end
                        if GVS_out_trajectory(i_time) < 0
                            % anode is on the left, cathode is on the right, illusory fall is towards the cathode, i.e. RIGHT
                            illusion_trajectory(i_time) = 1;
                        end
                    end
                elseif strcmp(stimulus_type, 'visual')
                    for i_time = 1 : length(time_labviewData)
                        if visual_scene_ml_translation__trajectory(i_time) > 0
                            % angle change is positive, horizon rotates counter-clockwise, illusion is to the RIGHT
                            illusion_trajectory(i_time) = 1;
                        end
                        if visual_scene_ml_translation__trajectory(i_time) < 0
                            % angle change is negative, horizon rotates clockwise, illusion is to the LEFT
                            illusion_trajectory(i_time) = -1;
                        end
                    end
                end
            end
            
            % get CoP trajectories
            if force_plate_data_available
                left_copx_trajectory = left_forceplate_cop_Acw(:, 1);
                right_copx_trajectory = right_forceplate_cop_Acw(:, 1);
                left_fz_trajectory = left_forceplate_wrench_Acw(:, 3);
                right_fz_trajectory = right_forceplate_wrench_Acw(:, 3);
            end
            
            %% find triggers
            %
            % Find the triggering events that indicate a stretch of interest. For perturbation experiments, this is the onset of
            % a perturbation. For unperturbed walking, this is any heelstrike.
            % The result is trigger_indices_labview.
            %
            if strcmp(stimulus_type, 'none')
                % use all touchdown events as triggers
                trigger_times = [left_touchdown_times; right_touchdown_times];
            else
                % find the time steps where the stimulus state crosses a threshold
                stimulus_threshold = 0.5;
                trigger_indices_labview = find(diff(sign(stimulus_state_trajectory - stimulus_threshold)) > 0) + 1;

                %
                epsilon = 1e-5;
                stim_start_indices_labview = find(diff(sign(abs(illusion_trajectory) - epsilon)) > 0) + 1;
                trigger_indices_labview = trigger_indices_labview(1 : length(stim_start_indices_labview)); % in case a stim is triggered, but not recorded
                
                trigger_times = time_labviewData(trigger_indices_labview);
            end
            % calculate indices
            trigger_indices_mocap = zeros(size(trigger_times));
            for i_index = 1 : length(trigger_times)
                [~, index_mocap] = min(abs(time_mocap - trigger_times(i_index)));
                trigger_indices_mocap(i_index) = index_mocap;
            end

            % visualize triggers
            if visualize
                if force_plate_data_available
                    left_cop_x_trajectory_relevant = left_copx_trajectory; left_cop_x_trajectory_relevant(left_cop_x_trajectory_relevant==0) = NaN;
                    right_cop_x_trajectory_relevant = right_copx_trajectory; right_cop_x_trajectory_relevant(right_cop_x_trajectory_relevant==0) = NaN;
                end
                figure; axes; hold on
%                 plot(time_labviewData, stimulus_state_trajectory*0.02);
                if force_plate_data_available
                    plot(time_forceplate, left_cop_x_trajectory_relevant, 'linewidth', 2, 'Displayname', 'cop left');
                    plot(time_forceplate, right_cop_x_trajectory_relevant, 'linewidth', 2, 'Displayname', 'cop right');
                end
%                 plot(time_forceplate(left_touchdown_indices_forceplate), zeros(size(left_touchdown_indices_forceplate)), 'o')
%                 plot(time_forceplate(right_touchdown_indices_forceplate), zeros(size(right_touchdown_indices_forceplate)), 'o')
                plot(time_mocap(trigger_indices_mocap), zeros(size(trigger_indices_mocap)), 'x', 'Displayname', 'heelstrikes')
                plot(time_labviewData, illusion_trajectory, 'Displayname', 'illusion')
%                 legend('stimulus state', 'left cop', 'right cop', 'left touchdown', 'right touchdown', 'trigger', 'stim start')
                legend('toggle')
            end

            %% extract event data
            %
            % For each trigger, determine the conditions and the relevant step events.
            % The result is 
            % stretch_start_indices_forceplate, stretch_end_indices_forceplate
            % stretch_start_times, stretch_end_times
            % condition_stance_foot_list, condition_perturbation_list, condition_delay_list, condition_index_list

            number_of_triggers = length(trigger_indices_mocap);
            removal_flags = zeros(number_of_triggers, 1);
            if strcmp(stimulus_type, 'none')
                stretch_start_times = zeros(number_of_triggers, 1);
                stretch_end_times = zeros(number_of_triggers, 1);
                closest_heelstrike_distance_times = zeros(number_of_triggers, 1);
                condition_stance_foot_list = cell(number_of_triggers, 1);
                condition_perturbation_list = cell(number_of_triggers, 1);
                condition_delay_list = cell(number_of_triggers, 1);
                condition_index_list = cell(number_of_triggers, 1);
                condition_experimental_list = cell(number_of_triggers, 1);
                
                for i_trigger = 1 : number_of_triggers
                    condition_perturbation_list{i_trigger, 1} = 'CONTROL';
                    condition_delay_list{i_trigger, 1} = 'CONTROL';
                    condition_index_list{i_trigger, 1} = 'CONTROL';
                    condition_experimental_list{i_trigger, 1} = condition;
                
                    % find out which heelstrike triggered
                    % XXX change this to use the interval, same as in the stimulus case
                    [distance_to_trigger_left_time, index_left] = min(abs(left_touchdown_times - trigger_times(i_trigger)));
                    [distance_to_trigger_right_time, index_right] = min(abs(right_touchdown_times - trigger_times(i_trigger)));
                    closest_heelstrike_distance_time = min([distance_to_trigger_left_time distance_to_trigger_right_time]);
                    
                    if distance_to_trigger_left_time < distance_to_trigger_right_time
                        condition_stance_foot_list{i_trigger, 1} = 'STANCE_LEFT';
                        stretch_start_times(i_trigger, 1) = trigger_times(i_trigger);
                        stretch_end_time_index = find(right_touchdown_times > trigger_times(i_trigger), 1, 'first');
                        if isempty(stretch_end_time_index)
                            stretch_end_times(i_trigger, 1) = -1;
                            removal_flags(i_trigger) = 1;
                        else
                            stretch_end_times(i_trigger, 1) = right_touchdown_times(stretch_end_time_index);
                        end
                    end    
                    if distance_to_trigger_right_time < distance_to_trigger_left_time
                        condition_stance_foot_list{i_trigger, 1} = 'STANCE_RIGHT';
                        stretch_start_times(i_trigger, 1) = trigger_times(i_trigger);
                        stretch_end_time_index = find(left_touchdown_times > trigger_times(i_trigger), 1, 'first');
                        if isempty(stretch_end_time_index)
                            stretch_end_times(i_trigger, 1) = -1;
                            removal_flags(i_trigger) = 1;
                        else
                            stretch_end_times(i_trigger, 1) = left_touchdown_times(stretch_end_time_index);
                        end
                    end                        
                    condition_experimental_list{i_trigger, 1} = condition;
                   
                end
                
                % remove flagged triggers
                unflagged_indices = ~removal_flags;
                stretch_start_times = stretch_start_times(unflagged_indices, :);
                stretch_end_times = stretch_end_times(unflagged_indices, :);
                condition_experimental_list = condition_experimental_list(unflagged_indices, :);
                condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
                closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :); 
                
                
                if visualize
                    for i_trigger = 1 : length(stretch_start_times)
                        if strcmp(condition_stance_foot_list(i_trigger), 'STANCE_RIGHT')
                            stretch_indicator_height = 0.01;
                        else
                            stretch_indicator_height = -0.01;
                        end
                            
                        plot([stretch_start_times(i_trigger) stretch_end_times(i_trigger)], [1 1]*stretch_indicator_height, 'linewidth', 3);
                    end
                end
                
                % check step times and flag outliers
                number_of_stretches = length(stretch_start_times);
                stretch_times = stretch_end_times - stretch_start_times;
                stretch_time_outlier_limits = median(stretch_times) * [0.8 1.2];
                
                removal_flags = zeros(number_of_stretches, 1);
                removal_flags(stretch_times < stretch_time_outlier_limits(1)) = 1;
                removal_flags(stretch_times > stretch_time_outlier_limits(2)) = 1;
                
                % check data availability and flag stretches with gaps
                for i_stretch = 1 : number_of_stretches
                    [~, start_index_mocap] = min(abs(time_mocap - stretch_start_times(i_stretch)));
                    [~, end_index_mocap] = min(abs(time_mocap - stretch_end_times(i_stretch)));
                    if any(any(isnan(marker_trajectories(start_index_mocap : end_index_mocap, :))))
                        removal_flags(i_stretch) = 1;
                    end
                end
                
                
                % remove flagged triggers
                unflagged_indices = ~removal_flags;
                stretch_start_times = stretch_start_times(unflagged_indices, :);
                stretch_end_times = stretch_end_times(unflagged_indices, :);
                condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
                condition_experimental_list = condition_experimental_list(unflagged_indices, :);
                closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :);
                
                % save data
                data_stretches_file_name = ['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'relevantDataStretches')];
                save ...
                  ( ...
                    data_stretches_file_name, ...
                    'condition_stance_foot_list', ...
                    'condition_experimental_list', ...
                    'stretch_start_times', ...
                    'stretch_end_times' ...
                  )

                disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' data_stretches_file_name]);                
                
                
                
                
            end  
            if ~strcmp(stimulus_type, 'none')
                % for each trigger, extract conditions and relevant step events
                number_of_triggers = length(trigger_indices_mocap);
                removal_flags = zeros(number_of_triggers, 1);
                stim_start_indices_labview = [stim_start_indices_labview zeros(number_of_triggers, 5)]; %#ok<AGROW>
                stretch_start_times = zeros(number_of_triggers, 6);
                stretch_end_times = zeros(number_of_triggers, 6);
                closest_heelstrike_distance_times = zeros(number_of_triggers, 1);
                condition_stance_foot_list = cell(number_of_triggers, 6);
                condition_perturbation_list = cell(number_of_triggers, 6);
                condition_delay_list = cell(number_of_triggers, 6);
                condition_index_list = cell(number_of_triggers, 6);
                condition_experimental_list = cell(number_of_triggers, 6);
                
                for i_trigger = 1 : number_of_triggers
                   % perturbation condition
                    if illusion_trajectory(stim_start_indices_labview(i_trigger)) > 0
                        condition_perturbation_list{i_trigger, 1} = 'ILLUSION_RIGHT';
                        condition_perturbation_list{i_trigger, 2} = 'ILLUSION_RIGHT';
                        condition_perturbation_list{i_trigger, 3} = 'ILLUSION_RIGHT';
                        condition_perturbation_list{i_trigger, 4} = 'ILLUSION_RIGHT';
                    elseif illusion_trajectory(stim_start_indices_labview(i_trigger)) < 0
                        condition_perturbation_list{i_trigger, 1} = 'ILLUSION_LEFT';
                        condition_perturbation_list{i_trigger, 2} = 'ILLUSION_LEFT';
                        condition_perturbation_list{i_trigger, 3} = 'ILLUSION_LEFT';
                        condition_perturbation_list{i_trigger, 4} = 'ILLUSION_LEFT';
                    else
        %                 disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_labviewData(trigger_indices_labview(i_trigger))) ' - no stim']);
                    end
                    condition_perturbation_list{i_trigger, 5} = 'CONTROL';
                    condition_perturbation_list{i_trigger, 6} = 'CONTROL';
                
                    % delay condition
                    wait_time_stim = time_labviewData(stim_start_indices_labview(i_trigger)) - time_labviewData(trigger_indices_labview(i_trigger));
                    [~, wait_condition_index] = min(abs(wait_times - wait_time_stim));
                    condition_delay_list{i_trigger, 1} = wait_time_labels{wait_condition_index};
                    condition_delay_list{i_trigger, 2} = wait_time_labels{wait_condition_index};
                    condition_delay_list{i_trigger, 3} = wait_time_labels{wait_condition_index};
                    condition_delay_list{i_trigger, 4} = wait_time_labels{wait_condition_index};
                    condition_delay_list{i_trigger, 5} = 'CONTROL';
                    condition_delay_list{i_trigger, 6} = 'CONTROL';
                    
                    % experimental condition
                    condition_experimental_list{i_trigger, 1} = condition;
                    condition_experimental_list{i_trigger, 2} = condition;
                    condition_experimental_list{i_trigger, 3} = condition;
                    condition_experimental_list{i_trigger, 4} = condition;
                    condition_experimental_list{i_trigger, 5} = condition;
                    condition_experimental_list{i_trigger, 6} = condition;
                
                    % get closest heelstrike
                    [~, index_left] = min(abs(left_touchdown_times - trigger_times(i_trigger)));
                    [~, index_right] = min(abs(right_touchdown_times - trigger_times(i_trigger)));
                    
                    % is the closest left heelstrike within the acceptable interval?
                    closest_left_heelstrike = left_touchdown_times(index_left);
                    time_difference_left = closest_left_heelstrike - trigger_times(i_trigger); % where does the closest left heelstrike lie relative to the trigger?
                    if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_left && time_difference_left < time_to_nearest_heelstrike_after_trigger_threshold
                    	% left heelstrike is acceptable
                        left_heelstrike_acceptable = true;
                    else
                        left_heelstrike_acceptable = false;
                    end
                    
                    % is the closest right heelstrike within the acceptable interval?
                    closest_right_heelstrike = right_touchdown_times(index_right);
                    time_difference_right = closest_right_heelstrike - trigger_times(i_trigger); % where does the closest right heelstrike lie relative to the trigger?
                    if -time_to_nearest_heelstrike_before_trigger_threshold < time_difference_right && time_difference_right < time_to_nearest_heelstrike_after_trigger_threshold
                    	% right heelstrike is acceptable
                        right_heelstrike_acceptable = true;
                    else
                        right_heelstrike_acceptable = false;
                    end
                    
                    % accept the acceptable one
                    if left_heelstrike_acceptable && ~right_heelstrike_acceptable
                        % triggered by left heelstrike
                        trigger_foot = 'left';
                        closest_heelstrike_distance_times(i_trigger) = time_difference_left;
                    elseif ~left_heelstrike_acceptable && right_heelstrike_acceptable
                        % triggered by right heelstrike
                        trigger_foot = 'right';
                        closest_heelstrike_distance_times(i_trigger) = time_difference_right;
                    elseif left_heelstrike_acceptable && right_heelstrike_acceptable
                        removal_flags(i_trigger) = 1;
                    elseif ~left_heelstrike_acceptable && ~right_heelstrike_acceptable
                        removal_flags(i_trigger) = 1;
                    end
                    
%                     closest_heelstrike_distance_time = min([distance_to_trigger_left_time distance_to_trigger_right_time]);

                    % confirm that the trigger happened less than 100ms (or whatever time's set) before the actual heelstrike
%                     if closest_heelstrike_distance_time > time_to_nearest_heelstrike_before_trigger_threshold
%                         removal_flags(i_trigger) = 1;
%         %                 disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_labviewData(trigger_indices_labview(i_trigger))) ' - trigger not in sync with heelstrike']);
%                     end

                    if strcmp(trigger_foot, 'left')
                        % check whether time offset is positive or negative
%                         if left_touchdown_times(index_left) - trigger_times(i_trigger) >= 0
%                             closest_heelstrike_distance_times(i_trigger) = closest_heelstrike_distance_time;
%                         else
%                             closest_heelstrike_distance_times(i_trigger) = -closest_heelstrike_distance_time;
%                         end
                        if index_left == 1 || length(left_touchdown_times) < index_left + 1 || removal_flags(i_trigger) == 1
                            % data doesn't include previous or next step
                            removal_flags(i_trigger) = 1;
                            right_foot_heelstrike_minus_1 = NaN;
                            right_foot_heelstrike_0 = NaN;
                            right_foot_heelstrike_plus_1 = NaN;
                            right_foot_heelstrike_plus_2 = NaN;
                            left_foot_heelstrike_minus_1 = NaN;
                            left_foot_heelstrike_0 = NaN;
                            left_foot_heelstrike_plus_1 = NaN;
                            left_foot_heelstrike_plus_2 = NaN;
                        else
                            left_foot_heelstrike_minus_1    = left_touchdown_times(index_left-1);
                            left_foot_heelstrike_0          = left_touchdown_times(index_left);
                            left_foot_heelstrike_plus_1     = left_touchdown_times(index_left+1);
                            left_foot_heelstrike_plus_2     = left_touchdown_times(index_left+2);
                            left_foot_pushoff_minus_1       = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_minus_1));
                            left_foot_pushoff_0             = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_0));
                            left_foot_pushoff_plus_1        = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_plus_1));
                            left_foot_pushoff_plus_2        = min(left_pushoff_times(left_pushoff_times >= left_foot_heelstrike_plus_2));

                            right_foot_heelstrike_minus_1   = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_minus_1));
                            right_foot_heelstrike_0         = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_0));
                            right_foot_heelstrike_plus_1    = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_plus_1));
                            right_foot_heelstrike_plus_2    = min(right_touchdown_times(right_touchdown_times >= left_foot_heelstrike_plus_2));
                            right_foot_pushoff_minus_1      = max(right_pushoff_times(right_pushoff_times <= left_foot_pushoff_minus_1));
                            right_foot_pushoff_0            = max(right_pushoff_times(right_pushoff_times <= left_foot_pushoff_0));
                            right_foot_pushoff_plus_1       = max(right_pushoff_times(right_pushoff_times <= left_foot_pushoff_plus_1));
                            right_foot_pushoff_plus_2       = max(right_pushoff_times(right_pushoff_times <= left_foot_pushoff_plus_2));

                            % notify if events are not sorted properly
                            if ~issorted ...
                                  ( ...
                                    [ ...
                                      left_foot_heelstrike_minus_1 ...
                                      right_foot_pushoff_minus_1 ...
                                      right_foot_heelstrike_minus_1 ...
                                      left_foot_pushoff_minus_1 ...
                                      left_foot_heelstrike_0 ...
                                      right_foot_pushoff_0 ...
                                      right_foot_heelstrike_0 ...
                                      left_foot_pushoff_0 ...
                                      left_foot_heelstrike_plus_1 ...
                                      right_foot_pushoff_plus_1 ...
                                      right_foot_heelstrike_plus_1 ...
                                      left_foot_pushoff_plus_1 ...
                                      left_foot_heelstrike_plus_2 ...
                                      right_foot_pushoff_plus_2 ...
                                      right_foot_heelstrike_plus_2 ...
                                      left_foot_pushoff_plus_2 ...
                                    ] ...
                                  );
                                disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_labviewData(trigger_indices_labview(i_trigger)))]);
                            end
                            % check check
                            if visualize
                                plot([left_foot_heelstrike_minus_1 left_foot_heelstrike_0 left_foot_heelstrike_plus_1 left_foot_heelstrike_plus_2], [0 0 0 0]-0.01, 'v', 'linewidth', 3);
                                plot([left_foot_pushoff_minus_1  left_foot_pushoff_0 left_foot_pushoff_plus_1 left_foot_pushoff_plus_2], [0 0 0 0]-0.01, '^', 'linewidth', 3);
                                plot([right_foot_heelstrike_minus_1 right_foot_heelstrike_0 right_foot_heelstrike_plus_1 right_foot_heelstrike_plus_2], [0 0 0 0]+0.01, 'v', 'linewidth', 3);
                                plot([right_foot_pushoff_minus_1  right_foot_pushoff_0 right_foot_pushoff_plus_1 right_foot_pushoff_plus_2], [0 0 0 0]+0.01, '^', 'linewidth', 3);
                                % note: this can crash if one of thse events is empty, because we are plotting before we
                                % have checked that
                            end                
                        end
                    elseif strcmp(trigger_foot, 'right')
                        % check whether time offset is positive or negative
%                         if right_touchdown_times(index_right) - trigger_times(i_trigger) > 0
%                             closest_heelstrike_distance_times(i_trigger) = closest_heelstrike_distance_time;
%                         else
%                             closest_heelstrike_distance_times(i_trigger) = -closest_heelstrike_distance_time;
%                         end
                        if index_right == 1 || length(right_touchdown_times) < index_right + 1 || removal_flags(i_trigger) == 1
                            % data doesn't include previous or next step
                            removal_flags(i_trigger) = 1;
                            left_foot_heelstrike_minus_1 = NaN;
                            left_foot_heelstrike_0 = NaN;
                            left_foot_heelstrike_plus_1 = NaN;
                            left_foot_heelstrike_plus_2 = NaN;
                            right_foot_heelstrike_minus_1 = NaN;
                            right_foot_heelstrike_0 = NaN;
                            right_foot_heelstrike_plus_1 = NaN;
                            right_foot_heelstrike_plus_2 = NaN;
                        else
                            right_foot_heelstrike_minus_1       = right_touchdown_times(index_right-1);
                            right_foot_heelstrike_0             = right_touchdown_times(index_right);
                            right_foot_heelstrike_plus_1        = right_touchdown_times(index_right+1);
                            right_foot_heelstrike_plus_2        = right_touchdown_times(index_right+2);
                            right_foot_pushoff_minus_1          = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_minus_1));
                            right_foot_pushoff_0                = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_0));
                            right_foot_pushoff_plus_1           = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_plus_1));
                            right_foot_pushoff_plus_2           = min(right_pushoff_times(right_pushoff_times >= right_foot_heelstrike_plus_2));

                            left_foot_heelstrike_minus_1        = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_minus_1));
                            left_foot_heelstrike_0              = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_0));
                            left_foot_heelstrike_plus_1         = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_plus_1));
                            left_foot_heelstrike_plus_2         = min(left_touchdown_times(left_touchdown_times >= right_foot_heelstrike_plus_2));
                            left_foot_pushoff_minus_1           = max(left_pushoff_times(left_pushoff_times <= right_foot_pushoff_minus_1));
                            left_foot_pushoff_0                 = max(left_pushoff_times(left_pushoff_times <= right_foot_pushoff_0));
                            left_foot_pushoff_plus_1            = max(left_pushoff_times(left_pushoff_times <= right_foot_pushoff_plus_1));
                            left_foot_pushoff_plus_2            = max(left_pushoff_times(left_pushoff_times <= right_foot_pushoff_plus_2));


                            % notify if events are not sorted properly
                            if ~issorted ...
                                  ( ...
                                    [ ...
                                      right_foot_heelstrike_minus_1 ...
                                      left_foot_pushoff_minus_1 ...
                                      left_foot_heelstrike_minus_1 ...
                                      right_foot_pushoff_minus_1 ...
                                      right_foot_heelstrike_0 ...
                                      left_foot_pushoff_0 ...
                                      left_foot_heelstrike_0 ...
                                      right_foot_pushoff_0 ...
                                      right_foot_heelstrike_plus_1 ...
                                      left_foot_pushoff_plus_1 ...
                                      left_foot_heelstrike_plus_1 ...
                                      right_foot_pushoff_plus_1 ...
                                      right_foot_heelstrike_plus_2 ...
                                      left_foot_pushoff_plus_2 ...
                                      left_foot_heelstrike_plus_2 ...
                                      right_foot_pushoff_plus_2 ...
                                    ] ...
                                  );
                              disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_labviewData(trigger_indices_labview(i_trigger)))]);

                            end

                            if visualize
                                plot([left_foot_heelstrike_minus_1 left_foot_heelstrike_0 left_foot_heelstrike_plus_1 left_foot_heelstrike_plus_2], [0 0 0 0]-0.01, 'v', 'linewidth', 3);
                                plot([left_foot_pushoff_minus_1  left_foot_pushoff_0 left_foot_pushoff_plus_1 left_foot_pushoff_plus_2], [0 0 0 0]-0.01, '^', 'linewidth', 3);
                                plot([right_foot_heelstrike_minus_1 right_foot_heelstrike_0 right_foot_heelstrike_plus_1 right_foot_heelstrike_plus_2], [0 0 0 0]+0.01, 'v', 'linewidth', 3);
                                plot([right_foot_pushoff_minus_1  right_foot_pushoff_0 right_foot_pushoff_plus_1 right_foot_pushoff_plus_2], [0 0 0 0]+0.01, '^', 'linewidth', 3);
                            end            
                        end            
                    else
                        trigger_foot = 'unclear';
                        disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_labviewData(trigger_indices_labview(i_trigger))) ' - trigger exactly between two heelstrikes']);
                        removal_flags(i_trigger) = 1;
                    end


                    % flag for removal if not all events are present
                    if any( ...
                            [ ...
                              isempty(left_foot_heelstrike_minus_1) ...
                              isempty(left_foot_heelstrike_0) ...
                              isempty(left_foot_heelstrike_plus_1) ...
                              isempty(left_foot_heelstrike_plus_2) ...
                              isempty(right_foot_heelstrike_minus_1) ...
                              isempty(right_foot_heelstrike_0) ...
                              isempty(right_foot_heelstrike_plus_1) ...
                              isempty(right_foot_heelstrike_plus_2) ...
                            ] ...
                          ) || removal_flags(i_trigger) == 1
                        % not all events are present
                        removal_flags(i_trigger) = 1;
                        left_foot_heelstrike_minus_1 = NaN;
                        left_foot_heelstrike_0 = NaN;
                        left_foot_heelstrike_plus_1 = NaN;
                        left_foot_heelstrike_plus_2 = NaN;
                        right_foot_heelstrike_minus_1 = NaN;
                        right_foot_heelstrike_0 = NaN;
                        right_foot_heelstrike_plus_1 = NaN;
                        right_foot_heelstrike_plus_2 = NaN;

                    else % all events are present
                        % mark relevant event delimiters depending on wait time and triggering foot
                        if strcmp(trigger_foot, 'right')
                            % triggered by right foot heelstrike
                            % triggering step
                            condition_index_list{i_trigger, 1} = 'ONE';
                            condition_stance_foot_list{i_trigger, 1} = 'STANCE_RIGHT';
                            % post stimulus steps
                            condition_index_list{i_trigger, 2} = 'TWO';
                            condition_stance_foot_list{i_trigger, 2} = 'STANCE_LEFT';
                            condition_index_list{i_trigger, 3} = 'THREE';
                            condition_stance_foot_list{i_trigger, 3} = 'STANCE_RIGHT';
                            condition_index_list{i_trigger, 4} = 'FOUR';
                            condition_stance_foot_list{i_trigger, 4} = 'STANCE_LEFT';
                            % use two previous steps as control
                            condition_index_list{i_trigger, 5} = 'CONTROL';
                            condition_stance_foot_list{i_trigger, 5} = 'STANCE_RIGHT';
                            condition_index_list{i_trigger, 6} = 'CONTROL';
                            condition_stance_foot_list{i_trigger, 6} = 'STANCE_LEFT';

                            if strcmp(stretches_to_analyze, 'single stance')
                                % this step
                                stretch_start_times(i_trigger, 1) = left_foot_pushoff_0;
                                stretch_end_times(i_trigger, 1) = left_foot_heelstrike_0;
                                % post-stimulus steps
                                stretch_start_times(i_trigger, 2) = right_foot_pushoff_0;
                                stretch_end_times(i_trigger, 2) = right_foot_heelstrike_plus_1;

                                stretch_start_times(i_trigger, 3) = left_foot_pushoff_plus_1;
                                stretch_end_times(i_trigger, 3) = left_foot_heelstrike_plus_1;
                                stretch_start_times(i_trigger, 4) = right_foot_pushoff_plus_1;
                                stretch_end_times(i_trigger, 4) = right_foot_heelstrike_plus_2;

                                % use two previous steps as control
                                stretch_start_times(i_trigger, 5) = left_foot_pushoff_minus_1;
                                stretch_end_times(i_trigger, 5) = left_foot_heelstrike_minus_1;
                                stretch_start_times(i_trigger, 6) = right_foot_pushoff_minus_1;
                                stretch_end_times(i_trigger, 6) = right_foot_heelstrike_0;

                            elseif strcmp(stretches_to_analyze, 'full stride')
                                % this step
                                stretch_start_times(i_trigger, 1) = right_foot_heelstrike_0;
                                stretch_end_times(i_trigger, 1) = left_foot_heelstrike_0;
                                % post-stimulus steps
                                stretch_start_times(i_trigger, 2) = left_foot_heelstrike_0;
                                stretch_end_times(i_trigger, 2) = right_foot_heelstrike_plus_1;

                                stretch_start_times(i_trigger, 3) = right_foot_heelstrike_plus_1;
                                stretch_end_times(i_trigger, 3) = left_foot_heelstrike_plus_1;
                                stretch_start_times(i_trigger, 4) = left_foot_heelstrike_plus_1;
                                stretch_end_times(i_trigger, 4) = right_foot_heelstrike_plus_2;
                                % use two previous steps as control
                                stretch_start_times(i_trigger, 5) = right_foot_heelstrike_minus_1;
                                stretch_end_times(i_trigger, 5) = left_foot_heelstrike_minus_1;
                                stretch_start_times(i_trigger, 6) = left_foot_heelstrike_minus_1;
                                stretch_end_times(i_trigger, 6) = right_foot_heelstrike_0;
                            end

                            % visualize
                            if visualize
                                plot([stretch_start_times(i_trigger, 1) stretch_end_times(i_trigger, 1)], [-1 1]*-0.01, 'color', 'r', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 2) stretch_end_times(i_trigger, 2)], [-1 1]*+0.01, 'color', 'g', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 3) stretch_end_times(i_trigger, 3)], [-1 1]*-0.01, 'color', 'b', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 4) stretch_end_times(i_trigger, 4)], [-1 1]*+0.01, 'color', 'm', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 5) stretch_end_times(i_trigger, 5)], [-1 1]*-0.01, 'color', 'y', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 6) stretch_end_times(i_trigger, 6)], [-1 1]*+0.01, 'color', 'k', 'linewidth', 3);
                            end

                        elseif strcmp(trigger_foot, 'left')
                            % triggered by left foot heelstrike
                            % this step
                            condition_index_list{i_trigger, 1} = 'ONE';
                            condition_stance_foot_list{i_trigger, 1} = 'STANCE_LEFT';
                            % post stimulus steps
                            condition_index_list{i_trigger, 2} = 'TWO';
                            condition_stance_foot_list{i_trigger, 2} = 'STANCE_RIGHT';
                            condition_index_list{i_trigger, 3} = 'THREE';
                            condition_stance_foot_list{i_trigger, 3} = 'STANCE_LEFT';
                            condition_index_list{i_trigger, 4} = 'FOUR';
                            condition_stance_foot_list{i_trigger, 4} = 'STANCE_RIGHT';
                            % use two previous steps as control
                            condition_index_list{i_trigger, 5} = 'CONTROL';
                            condition_stance_foot_list{i_trigger, 5} = 'STANCE_LEFT';
                            condition_index_list{i_trigger, 6} = 'CONTROL';
                            condition_stance_foot_list{i_trigger, 6} = 'STANCE_RIGHT';
                            if strcmp(stretches_to_analyze, 'single stance')
                                % this step
                                stretch_start_times(i_trigger, 1) = right_foot_pushoff_0;
                                stretch_end_times(i_trigger, 1) = right_foot_heelstrike_0;
                                % post stimulus steps
                                stretch_start_times(i_trigger, 2) = left_foot_pushoff_0;
                                stretch_end_times(i_trigger, 2) = left_foot_heelstrike_plus_1;
                                stretch_start_times(i_trigger, 3) = right_foot_pushoff_plus_1;
                                stretch_end_times(i_trigger, 3) = right_foot_heelstrike_plus_1;
                                stretch_start_times(i_trigger, 4) = left_foot_pushoff_plus_1;
                                stretch_end_times(i_trigger, 4) = left_foot_heelstrike_plus_2;
                                % use two previous steps as control
                                stretch_start_times(i_trigger, 3) = right_foot_pushoff_minus_1;
                                stretch_end_times(i_trigger, 3) = right_foot_heelstrike_minus_1;
                                stretch_start_times(i_trigger, 4) = left_foot_pushoff_minus_1;
                                stretch_end_times(i_trigger, 4) = left_foot_heelstrike_0;

                            elseif strcmp(stretches_to_analyze, 'full stride')
                                % this step
                                stretch_start_times(i_trigger, 1) = left_foot_heelstrike_0;
                                stretch_end_times(i_trigger, 1) = right_foot_heelstrike_0;
                                % post stimulus steps
                                stretch_start_times(i_trigger, 2) = right_foot_heelstrike_0;
                                stretch_end_times(i_trigger, 2) = left_foot_heelstrike_plus_1;
                                stretch_start_times(i_trigger, 3) = left_foot_heelstrike_plus_1;
                                stretch_end_times(i_trigger, 3) = right_foot_heelstrike_plus_1;
                                stretch_start_times(i_trigger, 4) = right_foot_heelstrike_plus_1;
                                stretch_end_times(i_trigger, 4) = left_foot_heelstrike_plus_2;
                                % use two previous steps as control
                                stretch_start_times(i_trigger, 5) = left_foot_heelstrike_minus_1;
                                stretch_end_times(i_trigger, 5) = right_foot_heelstrike_minus_1;
                                stretch_start_times(i_trigger, 6) = right_foot_heelstrike_minus_1;
                                stretch_end_times(i_trigger, 6) = left_foot_heelstrike_0;
                            end


                            % visualize
                            if visualize
                                plot([stretch_start_times(i_trigger, 1) stretch_end_times(i_trigger, 1)], [-1 1]*0.01, 'color', 'r', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 2) stretch_end_times(i_trigger, 2)], [-1 1]*-0.01, 'color', 'g', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 3) stretch_end_times(i_trigger, 3)], [-1 1]*0.01, 'color', 'b', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 4) stretch_end_times(i_trigger, 4)], [-1 1]*-0.01, 'color', 'm', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 5) stretch_end_times(i_trigger, 5)], [-1 1]*0.01, 'color', 'y', 'linewidth', 3);
                                plot([stretch_start_times(i_trigger, 6) stretch_end_times(i_trigger, 6)], [-1 1]*-0.01, 'color', 'k', 'linewidth', 3);
                            end
                        else
                            removal_flags(i_trigger) = 1;
                            condition_stance_foot_list{i_trigger, 1} = 'UNCLEAR';
                            condition_stance_foot_list{i_trigger, 2} = 'UNCLEAR';
                            condition_stance_foot_list{i_trigger, 3} = 'UNCLEAR';
                            condition_stance_foot_list{i_trigger, 4} = 'UNCLEAR';
                            condition_stance_foot_list{i_trigger, 5} = 'UNCLEAR';
                            condition_stance_foot_list{i_trigger, 6} = 'UNCLEAR';
                        end
                    end
                end

                % remove flagged triggers
                unflagged_indices = ~removal_flags;
                stim_start_indices_labview = stim_start_indices_labview(unflagged_indices, :);
                stretch_start_times = stretch_start_times(unflagged_indices, :);
                stretch_end_times = stretch_end_times(unflagged_indices, :);
                condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
                condition_perturbation_list = condition_perturbation_list(unflagged_indices, :);
                condition_delay_list = condition_delay_list(unflagged_indices, :);
                condition_index_list = condition_index_list(unflagged_indices, :);
                condition_experimental_list = condition_experimental_list(unflagged_indices, :); % XXX needs to be updated
                closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :);
                
                % reorder
                stretch_start_times = ...
                  [ ...
                    stretch_start_times(:, 1); ...
                    stretch_start_times(:, 2); ...
                    stretch_start_times(:, 3); ...
                    stretch_start_times(:, 4); ...
                    stretch_start_times(:, 5); ...
                    stretch_start_times(:, 6); ...
                  ];
                stretch_end_times = ...
                  [ ...
                    stretch_end_times(:, 1); ...
                    stretch_end_times(:, 2); ...
                    stretch_end_times(:, 3); ...
                    stretch_end_times(:, 4); ...
                    stretch_end_times(:, 5); ...
                    stretch_end_times(:, 6); ...
                  ];
                condition_stance_foot_list = ...
                    [ ...
                      condition_stance_foot_list(:, 1); ...
                      condition_stance_foot_list(:, 2); ...
                      condition_stance_foot_list(:, 3); ...
                      condition_stance_foot_list(:, 4); ...
                      condition_stance_foot_list(:, 5); ...
                      condition_stance_foot_list(:, 6); ...
                    ];
                condition_perturbation_list = ...
                    [ ...
                      condition_perturbation_list(:, 1); ...
                      condition_perturbation_list(:, 2); ...
                      condition_perturbation_list(:, 3); ...
                      condition_perturbation_list(:, 4); ...
                      condition_perturbation_list(:, 5); ...
                      condition_perturbation_list(:, 6); ...
                    ];
                condition_delay_list = ...
                    [ ...
                      condition_delay_list(:, 1); ...
                      condition_delay_list(:, 2); ...
                      condition_delay_list(:, 3); ...
                      condition_delay_list(:, 4); ...
                      condition_delay_list(:, 5); ...
                      condition_delay_list(:, 6); ...
                    ];
                condition_index_list = ...
                    [ ...
                      condition_index_list(:, 1); ...
                      condition_index_list(:, 2); ...
                      condition_index_list(:, 3); ...
                      condition_index_list(:, 4); ...
                      condition_index_list(:, 5); ...
                      condition_index_list(:, 6); ...
                    ];
                condition_experimental_list = ...
                    [ ...
                      condition_experimental_list(:, 1); ...
                      condition_experimental_list(:, 2); ...
                      condition_experimental_list(:, 3); ...
                      condition_experimental_list(:, 4); ...
                      condition_experimental_list(:, 5); ...
                      condition_experimental_list(:, 6); ...
                    ];
                
                
                % we now have a neatly ordered list of stretches which we can prune
                % check step times and flag outliers
                number_of_stretches = length(stretch_start_times);
                stretch_times = stretch_end_times - stretch_start_times;
                stretch_time_outlier_limits = median(stretch_times) * [0.8 1.2];

                removal_flags = zeros(number_of_stretches, 1);
                removal_flags(stretch_times < stretch_time_outlier_limits(1)) = 1;
                removal_flags(stretch_times > stretch_time_outlier_limits(2)) = 1;

                % check data availability and flag stretches with gaps
%                 for i_stretch = 1 : number_of_stretches
%                     [~, start_index_mocap] = min(abs(time_mocap - stretch_start_times(i_stretch)));
%                     [~, end_index_mocap] = min(abs(time_mocap - stretch_end_times(i_stretch)));
%                     if any(any(isnan(marker_trajectories(start_index_mocap : end_index_mocap, :))))
%                         removal_flags(i_stretch) = 1;
%                     end
%                 end
% XXX removed for now. Only those stretches with gaps in markers that actually interest me should be flagged.
% information about markers of interest should go to the experiment file

                % remove flagged triggers
                unflagged_indices = ~removal_flags;
                stretch_start_times = stretch_start_times(unflagged_indices, :);
                stretch_end_times = stretch_end_times(unflagged_indices, :);
                condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
                condition_perturbation_list = condition_perturbation_list(unflagged_indices, :);
                condition_delay_list = condition_delay_list(unflagged_indices, :);
                condition_index_list = condition_index_list(unflagged_indices, :);
                condition_experimental_list = condition_experimental_list(unflagged_indices, :);
                                
                % save data
                data_stretches_file_name = ['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'relevantDataStretches')];
                save ...
                  ( ...
                    data_stretches_file_name, ...
                    'condition_stance_foot_list', ...
                    'condition_perturbation_list', ...
                    'condition_delay_list', ...
                    'condition_index_list', ...
                    'condition_experimental_list', ...
                    'stretch_start_times', ...
                    'stretch_end_times', ...
                    'closest_heelstrike_distance_times' ...
                  )

                disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' data_stretches_file_name]);                
                
                
            end            
            

            


%             % extract and time-normalize data
%             number_of_stretches = length(stretch_start_indices_forceplate);
%             stim_start_time_relative_to_stretch = zeros(number_of_stretches, 1);
%             number_of_indices_per_step = zeros(number_of_stretches, 1);
%             start_indices_mocap = zeros(number_of_stretches, 1);
%             end_indices_mocap = zeros(number_of_stretches, 1);
%             start_indices_forceplate = zeros(number_of_stretches, 1);
%             end_indices_forceplate = zeros(number_of_stretches, 1);
%             start_indices_labview = zeros(number_of_stretches, 1);
%             end_indices_labview = zeros(number_of_stretches, 1);
%             start_indices_emg = zeros(number_of_stretches, 1);
%             end_indices_emg = zeros(number_of_stretches, 1);
%             removal_flags = zeros(number_of_stretches, 1);
%             for i_stretch = 1 : number_of_stretches
% 
%                 start_index_forceplate = stretch_start_indices_forceplate(i_stretch);
%                 end_index_forceplate = stretch_end_indices_forceplate(i_stretch);
%                 stretch_start_time = time_forceplate(start_index_forceplate);
%                 stretch_end_time = time_forceplate(end_index_forceplate);
%                     if remove_crossover_steps
%                         % determine if this step should be removed due to crossover of the contralateral foot onto the ipsilateral forceplate
%                         if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')
%             %                 stance_foot_cop_stretch = right_copx_trajectory(start_index_forceplate : end_index_forceplate);
%             %                 touchdown_index_forceplate = stretch_start_index_forceplate + find(stance_foot_cop_stretch~=0, 1, 'first') - 1;
%                             stance_foot_fz_step = right_fz_trajectory(start_index_forceplate : end_index_forceplate);
%                             swing_foot_fz_step = left_fz_trajectory(start_index_forceplate : end_index_forceplate);
%             %                 stance_foot_fz_prelude = right_fz_trajectory(stretch_start_index_forceplate-number_of_prelude_time_steps : touchdown_index_forceplate-1);
%                         elseif strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')
%             %                 stance_foot_cop_stretch = left_copx_trajectory(start_index_forceplate : end_index_forceplate);
%             %                 touchdown_index_forceplate = stretch_start_index_forceplate + find(stance_foot_cop_stretch~=0, 1, 'first') - 1;
%                             stance_foot_fz_step = left_fz_trajectory(start_index_forceplate : end_index_forceplate);
%                             swing_foot_fz_step = right_fz_trajectory(start_index_forceplate : end_index_forceplate);
%             %                 stance_foot_fz_prelude = left_fz_trajectory(touchdown_index_forceplate-number_of_prelude_time_steps : touchdown_index_forceplate-1);
%                         end
%             %             touchdown_time_forceplate = time_forceplate(touchdown_index_forceplate);
%                         time_extracted_forceplate = time_forceplate(start_index_forceplate : end_index_forceplate);
% 
%                         % check if this stretch is usable
%                         number_of_indices_per_step(i_stretch) = length(time_extracted_forceplate);
%                         number_of_swing_foot_fz_zeros_stretch = sum(abs(swing_foot_fz_step) < swing_foot_fz_zero_threshold);
%                         number_of_swing_foot_fz_nonzeros_stretch = sum(abs(swing_foot_fz_step) >= swing_foot_fz_zero_threshold);
% 
%             %             number_of_swing_foot_fz_zeros_prelude = sum(-stance_foot_fz_prelude < swing_foot_fz_zero_threshold);
% 
%                         if number_of_swing_foot_fz_nonzeros_stretch > swing_foot_nonzero_stretch_length_threshold_indices
%                             removal_flags(i_stretch) = 1;
%                             disp(['excluding stretch index ' num2str(i_stretch) ' - cross over at time ' num2str(stretch_start_time)]);
%                         end
% 
% 
% 
% 
% 
% 
%             %             if number_of_swing_foot_fz_zeros_prelude < swing_foot_zero_stretch_length_threshold_indices
%             %                 removal_flags(i_stretch) = 1;
%             %                 disp(['excluding stretch index ' num2str(i_stretch) ' - cross over at time ' num2str(time_forceplate(stretch_start_index_forceplate))]);
%             %             end
%             %             if isempty(touchdown_index_forceplate)
%             %                 touchdown_index_forceplate = stretch_start_index_forceplate;
%             %                 touchdown_time_forceplate = time_forceplate(touchdown_index_forceplate);
%             %                 removal_flags(i_stretch) = 1;
%             %             end
%                     end


        %         plot(stance_foot_fz_step, 'displayname', num2str(i_stretch));
        %         plot(swing_foot_fz_step, 'displayname', num2str(i_stretch));
        %         plot(stance_foot_cop_stretch, 'displayname', num2str(i_stretch));


%                 % time
%                 [~, start_index_mocap] = min(abs(time_mocap - stretch_start_time));
%                 [~, end_index_mocap] = min(abs(time_mocap - stretch_end_time));
%                 [~, start_index_labview] = min(abs(time_labviewData - stretch_start_time));
%                 [~, end_index_labview] = min(abs(time_labviewData - stretch_end_time));
%                 if emg_present
%                     [~, start_index_emg] = min(abs(time_emg - stretch_start_time));
%                     [~, end_index_emg] = min(abs(time_emg - stretch_end_time));
%                 end
% 
%                 % store
%                 start_indices_mocap(i_stretch) = start_index_mocap;
%                 end_indices_mocap(i_stretch) = end_index_mocap;
%                 start_indices_forceplate(i_stretch) = start_index_forceplate;
%                 end_indices_forceplate(i_stretch) = end_index_forceplate;
%                 start_indices_labview(i_stretch) = start_index_labview;
%                 end_indices_labview(i_stretch) = end_index_labview;
%                 if emg_present
%                     start_indices_emg(i_stretch) = start_index_emg;
%                     end_indices_emg(i_stretch) = end_index_emg;
%                 end
%                 if ~strcmp(stimulus_type, 'none')
%                     if stim_start_indices_labview(i_stretch) == 0
%                         stim_start_time_relative_to_stretch(i_stretch) = NaN;
%                     else
%                         stim_start_time = time_labviewData(stim_start_indices_labview(i_stretch));
%                         stim_start_time_relative_to_stretch(i_stretch) = stim_start_time - stretch_start_time;
%                     end
%                 end
% 
% 
%                 % extract force plate data
%                 left_cop_x_extracted_stretch = left_copx_trajectory(start_index_forceplate : end_index_forceplate);
%                 right_cop_x_extracted_stretch = right_copx_trajectory(start_index_forceplate : end_index_forceplate);
% 
%                 if remove_crossover_steps && strcmp(stretches_to_analyze, 'single stance')
%                     % check for data correctness - one force plate data stretch should have no zeros
%                     if any(left_cop_x_extracted_stretch==0) & any(right_cop_x_extracted_stretch==0)
%                         removal_flags(i_stretch) = 1;
%                     end
%                 end
%             end
% 
%             % remove flagged stretches
%             unflagged_indices = ~removal_flags;
%             condition_stance_foot_list = condition_stance_foot_list(unflagged_indices);
%             condition_perturbation_list = condition_perturbation_list(unflagged_indices);
%             condition_delay_list = condition_delay_list(unflagged_indices);
%             condition_index_list = condition_index_list(unflagged_indices);
%             start_indices_mocap = start_indices_mocap(unflagged_indices);
%             end_indices_mocap = end_indices_mocap(unflagged_indices);
%             start_indices_forceplate = start_indices_forceplate(unflagged_indices);
%             end_indices_forceplate = end_indices_forceplate(unflagged_indices);
%             start_indices_labview = start_indices_labview(unflagged_indices);
%             end_indices_labview = end_indices_labview(unflagged_indices);
%             start_indices_emg = start_indices_emg(unflagged_indices);
%             end_indices_emg = end_indices_emg(unflagged_indices);
%             stim_start_time_relative_to_stretch = stim_start_time_relative_to_stretch(unflagged_indices);
% 
% 
%             % prepare stretches for further processing
%             data_points_to_process_mocap = [];
%             number_of_padding_steps = 30;
%             for i_stretch = 1 : size(start_indices_mocap)
%                 % get start and end indices and apply padding
%                 start_index_mocap_stretch = start_indices_mocap(i_stretch) - number_of_padding_steps;
%                 end_index_mocap_stretch = end_indices_mocap(i_stretch) + number_of_padding_steps;
% 
%                 % crop if necessary
%                 if start_index_mocap_stretch < 1
%                     start_index_mocap_stretch = 1;
%                 end
%                 if end_index_mocap_stretch > size(marker_trajectories, 1);
%                     end_index_mocap_stretch = size(marker_trajectories, 1);
%                 end
% 
%                 time_steps_stretch = start_index_mocap_stretch : end_index_mocap_stretch;
%                 data_points_to_process_mocap = [data_points_to_process_mocap time_steps_stretch];
%             end
% 
%             %% save data
%             file_name = makeFileName(date, subject_id, 'walking', i_trial, 'relevantDataStretches');
%             save ...
%               ( ...
%                 file_name, ...
%                 'condition_stance_foot_list', ...
%                 'condition_perturbation_list', ...
%                 'condition_delay_list', ...
%                 'condition_index_list', ...
%                 'start_indices_mocap', ...
%                 'end_indices_mocap', ...
%                 'start_indices_forceplate', ...
%                 'end_indices_forceplate', ...
%                 'start_indices_labview', ...
%                 'end_indices_labview', ...
%                 'start_indices_emg', ...
%                 'end_indices_emg', ...
%                 'stim_start_time_relative_to_stretch', ...
%                 'data_points_to_process_mocap' ...
%               )
% 
%             disp(['Trial ' num2str(i_trial) ' completed']);

%         end
    end

end





















