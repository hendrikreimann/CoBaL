% find the relevant stretches of data after the stimuli


%% Choose Data Type
emg_present                     = 1;
stimulus_type = 'gvs';
% stimulus_type = 'visual';
% stimulus_type = 'none';

ignore_crossover                = 0;

visualize_triggers              = 0;

%% prepare
wait_times = [0 0.150 0.450];
wait_time_labels = {'0ms', '150ms', '450ms'};
load subjectInfo.mat;

% trials_to_process = 1 : 23;
% trials_to_process = [2:12 14:15 18:20];
% trials_to_process = 12 : 23;
trials_to_process = 2:21;
% trials_to_process = 1;
trials_to_process = 3 : 43;
% trials_to_process = 1 : 20;
% trials_to_process = 3;

total_positive_steps = [];
total_negative_steps = [];
total_polarity_condition_list = [];
total_step_condition_list = [];

swing_foot_fz_zero_threshold = 20; % threshold for counting a vertical force reading as zero, in Nm
swing_foot_zero_stretch_length_threshold_time = 0.2; % the force plate should register zero for at least this long
duration_until_nearest_future_heelstrike_threshold = 0.20; % a heelstrike should happen less than this long after a trigger
duration_prelude = 0.4; % duration that should have zero vertical force prior to a heelstrike, to identify crossover


%% extract data
stretch_length_indices_forceplate = [];
condition_stance_foot_list = {};
condition_polarity_list = {};
condition_delay_list = {};
stim_start_time_relative_to_stretch = [];

for i_trial = trials_to_process
    % load data
    load(makeFileName(date, subject_id, 'walking', i_trial, 'markerTrajectories'));
    load(makeFileName(date, subject_id, 'walking', i_trial, 'labviewTrajectories'));
    load(makeFileName(date, subject_id, 'walking', i_trial, 'forceplateTrajectories'));
    load(makeFileName(date, subject_id, 'walking', i_trial, 'stepEvents'));
    if emg_present
        load(makeFileName(date, subject_id, 'walking', i_trial, 'emgTrajectories'));
    end

    if isnan(sampling_rate_forceplate)
        sampling_rate_forceplate = 1 / median(diff(time_forceplate));
    end
    number_of_prelude_time_steps = ceil(duration_prelude * sampling_rate_forceplate);
    swing_foot_zero_stretch_length_threshold_indices = ceil(swing_foot_zero_stretch_length_threshold_time * 1/sampling_rate_forceplate);

    if strcmp(stimulus_type, 'gvs')
        stim_sent_trajectory = gvs_out_trajectory;
    elseif strcmp(stimulus_type, 'visual')
        stim_sent_trajectory = visual_shift_ml_trajectory;
    end

    left_copx_trajectory = left_forceplate_cop_Acw(:, 1);
    right_copx_trajectory = right_forceplate_cop_Acw(:, 1);
    left_fz_trajectory = left_forceplate_wrench_Acw(:, 3);
    right_fz_trajectory = right_forceplate_wrench_Acw(:, 3);

    % find trigger
    if strcmp(stimulus_type, 'none')
        trigger_indices_labview = [left_touchdown_indices_labview right_touchdown_indices_labview];
    else
        % -- find the first index where the heel marker crosses the threshold
        % and is then recognized as heel strike
        trigger_indices_labview = find(diff(sign(stimulus_foot_state - 0.5)) > 0) + 1;
        epsilon = 1e-5;
        stim_start_indices_labview = find(diff(sign(abs(stim_sent_trajectory) - epsilon)) > 0) + 1;
        trigger_indices_labview = trigger_indices_labview(1 : length(stim_start_indices_labview)); % in case a stim is triggered, but not recorded
    end

    % visualize triggers
    if visualize_triggers
        left_cop_x_trajectory_relevant = left_copx_trajectory; left_cop_x_trajectory_relevant(left_cop_x_trajectory_relevant==0) = NaN;
        right_cop_x_trajectory_relevant = right_copx_trajectory; right_cop_x_trajectory_relevant(right_cop_x_trajectory_relevant==0) = NaN;
        figure; axes; hold on
        plot(time_labview, stimulus_foot_state);
        plot(time_forceplate, left_cop_x_trajectory_relevant, 'linewidth', 2);
        plot(time_forceplate, right_cop_x_trajectory_relevant, 'linewidth', 2);
        plot(time_forceplate(left_touchdown_indices_forceplate), zeros(size(left_touchdown_indices_forceplate)), 'o')
        plot(time_forceplate(right_touchdown_indices_forceplate), zeros(size(right_touchdown_indices_forceplate)), 'o')
        plot(time_labview(trigger_indices_labview), zeros(size(trigger_indices_labview)), 'x')
%         plot(time_labview(stim_start_indices_labview), zeros(size(stim_start_indices_labview)), 'x')
        legend('stimulus state', 'left cop', 'right cop', 'left touchdown', 'right touchdown', 'trigger', 'stim start')
    end

    % for each trigger, extract conditions and relevant step events
    number_of_triggers = length(trigger_indices_labview);
    removal_flags = zeros(number_of_triggers, 1);
    if strcmp(stimulus_type, 'none')
        stretch_start_indices_forceplate = zeros(number_of_triggers, 1);
        stretch_end_indices_forceplate = zeros(number_of_triggers, 1);
        closest_heelstrike_distance_times = zeros(number_of_triggers, 1);
        condition_stance_foot_list = cell(number_of_triggers, 1);
        condition_polarity_list = cell(number_of_triggers, 1);
        condition_delay_list = cell(number_of_triggers, 1);
    else
        stim_start_indices_labview = [stim_start_indices_labview zeros(number_of_triggers, 1)]; %#ok<AGROW>
        stretch_start_indices_forceplate = zeros(number_of_triggers, 2);
        stretch_end_indices_forceplate = zeros(number_of_triggers, 2);
        closest_heelstrike_distance_times = zeros(number_of_triggers, 2);
        condition_stance_foot_list = cell(number_of_triggers, 2);
        condition_polarity_list = cell(number_of_triggers, 2);
        condition_delay_list = cell(number_of_triggers, 2);
    end

    for i_trigger = 1 : number_of_triggers
        % determine condition
        if strcmp(stimulus_type, 'none')
            condition_polarity_list{i_trigger, 1} = 'CONTROL';
            condition_delay_list{i_trigger, 1} = 'CONTROL';
        else
            % polarity condition
            if stim_sent_trajectory(stim_start_indices_labview(i_trigger)) > 0
                condition_polarity_list{i_trigger, 1} = 'POSITIVE';
            elseif stim_sent_trajectory(stim_start_indices_labview(i_trigger)) < 0
                condition_polarity_list{i_trigger, 1} = 'NEGATIVE';
            else
                disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_labview(trigger_indices_labview(i_trigger))) ' - no stim']);
            end
            condition_polarity_list{i_trigger, 2} = 'CONTROL';
            
            % delay
            wait_time_stim = time_labview(stim_start_indices_labview(i_trigger)) - time_labview(trigger_indices_labview(i_trigger));
            [~, wait_condition_index] = min(abs(wait_times - wait_time_stim));
            condition_delay_list{i_trigger, 1} = wait_time_labels{wait_condition_index};
            condition_delay_list{i_trigger, 2} = 'CONTROL';
        end
        

        % find out which heelstrike triggered
        [distance_to_trigger_left_time, index_left] = min(abs(time_forceplate(left_touchdown_indices_forceplate) - time_labview(trigger_indices_labview(i_trigger))));
        [distance_to_trigger_right_time, index_right] = min(abs(time_forceplate(right_touchdown_indices_forceplate) - time_labview(trigger_indices_labview(i_trigger))));
        closest_heelstrike_distance_time = min([distance_to_trigger_left_time distance_to_trigger_right_time]);
        
        % confirm that the trigger happened less than 100ms (or whatever time's set) before the actual heelstrike
        if ~strcmp(stimulus_type, 'none')
            if closest_heelstrike_distance_time > duration_until_nearest_future_heelstrike_threshold
                removal_flags(i_trigger) = 1;
                disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_labview(trigger_indices_labview(i_trigger))) ' - trigger not in sync with heelstrike']);
            end
        end
        
        if distance_to_trigger_left_time < distance_to_trigger_right_time
            % triggered by left heelstrike
            trigger_foot = 'left';
            % check whether time offset is positive or negative
            if left_touchdown_indices_forceplate(index_left) - trigger_indices_labview(i_trigger) > 0
                closest_heelstrike_distance_times(i_trigger) = closest_heelstrike_distance_time;
            else
                closest_heelstrike_distance_times(i_trigger) = -closest_heelstrike_distance_time;
            end
            if index_left == 1 || length(left_touchdown_indices_forceplate) < index_left + 1
                % data doesn't include previous or next step
                removal_flags(i_trigger) = 1;
                last_right_foot_heelstrike = NaN;
                this_right_foot_heelstrike = NaN;
                next_right_foot_heelstrike = NaN;
                last_left_foot_heelstrike = NaN;
                this_left_foot_heelstrike = NaN;
                next_left_foot_heelstrike = NaN;
            else
                last_left_foot_heelstrike = left_touchdown_indices_forceplate(index_left-1);
                this_left_foot_heelstrike = left_touchdown_indices_forceplate(index_left);
                next_left_foot_heelstrike = left_touchdown_indices_forceplate(index_left+1);
                
                last_right_foot_heelstrike = min(right_touchdown_indices_forceplate(right_touchdown_indices_forceplate >= last_left_foot_heelstrike));
                this_right_foot_heelstrike = min(right_touchdown_indices_forceplate(right_touchdown_indices_forceplate >= this_left_foot_heelstrike));
                next_right_foot_heelstrike = min(right_touchdown_indices_forceplate(right_touchdown_indices_forceplate >= next_left_foot_heelstrike));
                
%                 % check check
%                 figure; hold on
%                 plot([last_right_foot_heelstrike this_right_foot_heelstrike next_right_foot_heelstrike], [0 0 0], 'x');
%                 plot([last_left_foot_heelstrike this_left_foot_heelstrike next_left_foot_heelstrike], [0 0 0], 'x');
                
            end            
        elseif distance_to_trigger_right_time < distance_to_trigger_left_time
            % triggered by right heelstrike
            trigger_foot = 'right';
            % check whether time offset is positive or negative
            if right_touchdown_indices_forceplate(index_right) - trigger_indices_labview(i_trigger) > 0
                closest_heelstrike_distance_times(i_trigger) = closest_heelstrike_distance_time;
            else
                closest_heelstrike_distance_times(i_trigger) = -closest_heelstrike_distance_time;
            end
            if index_right == 1 || length(right_touchdown_indices_forceplate) < index_right + 1
                % data doesn't include previous or next step
                removal_flags(i_trigger) = 1;
                last_left_foot_heelstrike = NaN;
                this_left_foot_heelstrike = NaN;
                next_left_foot_heelstrike = NaN;
                last_right_foot_heelstrike = NaN;
                this_right_foot_heelstrike = NaN;
                next_right_foot_heelstrike = NaN;
            else
                last_right_foot_heelstrike = right_touchdown_indices_forceplate(index_right-1);
                this_right_foot_heelstrike = right_touchdown_indices_forceplate(index_right);
                next_right_foot_heelstrike = right_touchdown_indices_forceplate(index_right+1);
                
                last_left_foot_heelstrike = min(left_touchdown_indices_forceplate(left_touchdown_indices_forceplate >= last_right_foot_heelstrike));
                this_left_foot_heelstrike = min(left_touchdown_indices_forceplate(left_touchdown_indices_forceplate >= this_right_foot_heelstrike));
                next_left_foot_heelstrike = min(left_touchdown_indices_forceplate(left_touchdown_indices_forceplate >= next_right_foot_heelstrike));
                
                
%                 % check check
%                 figure; hold on
%                 plot([last_right_foot_heelstrike this_right_foot_heelstrike next_right_foot_heelstrike], [0 0 0], 'x');
%                 plot([last_left_foot_heelstrike this_left_foot_heelstrike next_left_foot_heelstrike], [0 0 0], 'x');
                
            end            
        else
            trigger_foot = 'unclear';
            disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_labview(trigger_indices_labview(i_trigger))) ' - trigger exactly between two heelstrikes']);
            removal_flags(i_trigger) = 1;
        end
        
        % flag for removal if not all events are present
        if any( ...
                [ ...
                  isempty(last_left_foot_heelstrike) ...
                  isempty(this_left_foot_heelstrike) ...
                  isempty(next_left_foot_heelstrike) ...
                  isempty(last_right_foot_heelstrike) ...
                  isempty(this_right_foot_heelstrike) ...
                  isempty(next_right_foot_heelstrike) ...
                ] ...
              )
            removal_flags(i_trigger) = 1;
            last_left_foot_heelstrike = NaN;
            this_left_foot_heelstrike = NaN;
            next_left_foot_heelstrike = NaN;
            last_right_foot_heelstrike = NaN;
            this_right_foot_heelstrike = NaN;
            next_right_foot_heelstrike = NaN;
            
        else
            % mark relevant event delimiters depending on wait time and triggering foot
            if strcmp(trigger_foot, 'right')
                % triggered by right foot heelstrike
                if strcmp(stimulus_type, 'none')
                    condition_stance_foot_list{i_trigger, 1} = 'RIGHT';
                    stretch_start_indices_forceplate(i_trigger, 1) = this_right_foot_heelstrike;
                    stretch_end_indices_forceplate(i_trigger, 1) = this_left_foot_heelstrike;
                else
                    if strcmp(condition_delay_list{i_trigger, 1}, '0ms') || strcmp(condition_delay_list{i_trigger, 1}, '150ms')
                        % this step is of interest
                        condition_stance_foot_list{i_trigger, 1} = 'RIGHT';
                        condition_stance_foot_list{i_trigger, 2} = 'RIGHT';
                        stretch_start_indices_forceplate(i_trigger, 1) = this_right_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 1) = this_left_foot_heelstrike;
                        stretch_start_indices_forceplate(i_trigger, 2) = last_right_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 2) = last_left_foot_heelstrike;
                    elseif strcmp(condition_delay_list{i_trigger, 1}, '450ms')
                        % next step is of interest
                        condition_stance_foot_list{i_trigger, 1} = 'LEFT';
                        condition_stance_foot_list{i_trigger, 2} = 'LEFT';
                        stretch_start_indices_forceplate(i_trigger, 1) = this_left_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 1) = next_right_foot_heelstrike;
                        stretch_start_indices_forceplate(i_trigger, 2) = last_left_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 2) = this_right_foot_heelstrike;
                    else
                        disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_forceplate(trigger_indices_labview(i_trigger))) ' - delay not defined']);
                    end
                end
            elseif strcmp(trigger_foot, 'left')
                % triggered by left foot heelstrike
                if strcmp(stimulus_type, 'none')
                    condition_stance_foot_list{i_trigger, 1} = 'LEFT';
                    stretch_start_indices_forceplate(i_trigger, 1) = this_left_foot_heelstrike;
                    stretch_end_indices_forceplate(i_trigger, 1) = this_right_foot_heelstrike;
                else
                    if strcmp(condition_delay_list{i_trigger, 1}, '0ms') || strcmp(condition_delay_list{i_trigger, 1}, '150ms')
                        % this step is of interest
                        condition_stance_foot_list{i_trigger, 1} = 'LEFT';
                        condition_stance_foot_list{i_trigger, 2} = 'LEFT';
                        stretch_start_indices_forceplate(i_trigger, 1) = this_left_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 1) = this_right_foot_heelstrike;
                        stretch_start_indices_forceplate(i_trigger, 2) = last_left_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 2) = last_right_foot_heelstrike;
                    elseif strcmp(condition_delay_list{i_trigger, 1}, '450ms')
                        % next step is of interest
                        condition_stance_foot_list{i_trigger, 1} = 'RIGHT';
                        condition_stance_foot_list{i_trigger, 2} = 'RIGHT';
                        stretch_start_indices_forceplate(i_trigger, 1) = this_right_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 1) = next_left_foot_heelstrike;
                        stretch_start_indices_forceplate(i_trigger, 2) = last_right_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 2) = this_left_foot_heelstrike;
                    else
                        disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_forceplate(trigger_indices_labview(i_trigger))) ' - delay not defined']);
                    end
                end
            else
                removal_flags(i_trigger) = 1;
                condition_stance_foot_list{i_trigger, 1} = 'UNCLEAR';
                condition_stance_foot_list{i_trigger, 2} = 'UNCLEAR';
            end
        end
    end

    % remove flagged triggers
    unflagged_indices = ~removal_flags;
    if ~strcmp(stimulus_type, 'none')
        stim_start_indices_labview = stim_start_indices_labview(unflagged_indices, :);
    end
    stretch_start_indices_forceplate = stretch_start_indices_forceplate(unflagged_indices, :);
    stretch_end_indices_forceplate = stretch_end_indices_forceplate(unflagged_indices, :);
    condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
    condition_polarity_list = condition_polarity_list(unflagged_indices, :);
    condition_delay_list = condition_delay_list(unflagged_indices, :);
    closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :);

    % form stretches
    if ~strcmp(stimulus_type, 'none')
        stim_start_indices_labview = [stim_start_indices_labview(:, 1); stim_start_indices_labview(:, 2)];
        stretch_start_indices_forceplate = [stretch_start_indices_forceplate(:, 1); stretch_start_indices_forceplate(:, 2)];
        stretch_end_indices_forceplate = [stretch_end_indices_forceplate(:, 1); stretch_end_indices_forceplate(:, 2)];
        condition_stance_foot_list = [condition_stance_foot_list(:, 1); condition_stance_foot_list(:, 2)];
        condition_polarity_list = [condition_polarity_list(:, 1); condition_polarity_list(:, 2)];
        condition_delay_list = [condition_delay_list(:, 1); condition_delay_list(:, 2)];
    end

    % check step times and remove outliers
    if strcmp(stimulus_type, 'none')
        number_of_stretches = length(stretch_start_indices_forceplate);
        stretch_times = time_forceplate(stretch_end_indices_forceplate) - time_forceplate(stretch_start_indices_forceplate);
        stretch_time_outlier_limits = median(stretch_times) * [0.8 1.2];
        removal_flags = zeros(number_of_stretches, 1);
        removal_flags(stretch_times < stretch_time_outlier_limits(1)) = 1;
        removal_flags(stretch_times > stretch_time_outlier_limits(2)) = 1;
        unflagged_indices = ~removal_flags;
        if ~strcmp(stimulus_type, 'none')
            stim_start_indices_labview = stim_start_indices_labview(unflagged_indices, :);
        end
        stretch_start_indices_forceplate = stretch_start_indices_forceplate(unflagged_indices, :);
        stretch_end_indices_forceplate = stretch_end_indices_forceplate(unflagged_indices, :);
        condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
        condition_polarity_list = condition_polarity_list(unflagged_indices, :);
        condition_delay_list = condition_delay_list(unflagged_indices, :);
        closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :);
    end
    
    % extract and time-normalize data
    number_of_stretches = length(stretch_start_indices_forceplate);
    stim_start_time_relative_to_stretch = zeros(number_of_stretches, 1);
    number_of_indices_per_step = zeros(number_of_stretches, 1);
    start_indices_mocap = zeros(number_of_stretches, 1);
    end_indices_mocap = zeros(number_of_stretches, 1);
    start_indices_forceplate = zeros(number_of_stretches, 1);
    end_indices_forceplate = zeros(number_of_stretches, 1);
    start_indices_labview = zeros(number_of_stretches, 1);
    end_indices_labview = zeros(number_of_stretches, 1);
    start_indices_emg = zeros(number_of_stretches, 1);
    end_indices_emg = zeros(number_of_stretches, 1);
    removal_flags = zeros(number_of_stretches, 1);
    for i_stretch = 1 : number_of_stretches
        % we have the delimiters for the stretches of interest, based on touchdown and pushoff identification from mocap
        % now look at the force plate data and decide whether this stretch is usable or not due to cross-over
        % if the stretch is usable, find the touchdown and pushoff events from force plate loading

        stretch_start_index_forceplate = stretch_start_indices_forceplate(i_stretch);
        stretch_end_index_forceplate = stretch_end_indices_forceplate(i_stretch);
        

        % find start of CoP data
        if strcmp(condition_stance_foot_list{i_stretch}, 'RIGHT')
            stance_foot_cop_stretch = right_copx_trajectory(stretch_start_index_forceplate : stretch_end_index_forceplate);
            touchdown_index_forceplate = stretch_start_index_forceplate + find(stance_foot_cop_stretch~=0, 1, 'first') - 1;
            stance_foot_fz_step = right_fz_trajectory(touchdown_index_forceplate : stretch_end_index_forceplate);
            swing_foot_fz_step = left_fz_trajectory(touchdown_index_forceplate : stretch_end_index_forceplate);
            stance_foot_fz_prelude = right_fz_trajectory(touchdown_index_forceplate-number_of_prelude_time_steps : touchdown_index_forceplate-1);
        elseif strcmp(condition_stance_foot_list{i_stretch}, 'LEFT')
            stance_foot_cop_stretch = left_copx_trajectory(stretch_start_index_forceplate : stretch_end_index_forceplate);
            touchdown_index_forceplate = stretch_start_index_forceplate + find(stance_foot_cop_stretch~=0, 1, 'first') - 1;
            stance_foot_fz_step = left_fz_trajectory(touchdown_index_forceplate : stretch_end_index_forceplate);
            swing_foot_fz_step = right_fz_trajectory(touchdown_index_forceplate : stretch_end_index_forceplate);
            stance_foot_fz_prelude = left_fz_trajectory(touchdown_index_forceplate-number_of_prelude_time_steps : touchdown_index_forceplate-1);
        end
        touchdown_time_forceplate = time_forceplate(touchdown_index_forceplate);
        time_extracted_forceplate = time_forceplate(touchdown_index_forceplate : stretch_end_index_forceplate);

        % check if this stretch is usable
        number_of_indices_per_step(i_stretch) = length(time_extracted_forceplate);
        number_of_swing_foot_fz_zeros_stretch = sum(-swing_foot_fz_step < swing_foot_fz_zero_threshold);
        number_of_swing_foot_fz_zeros_prelude = sum(-stance_foot_fz_prelude < swing_foot_fz_zero_threshold);
        
%         plot(stance_foot_fz_step, 'displayname', num2str(i_stretch));
%         plot(swing_foot_fz_step, 'displayname', num2str(i_stretch));
%         plot(stance_foot_cop_stretch, 'displayname', num2str(i_stretch));
        
% XXX removed these three conditionals for the balance beam data. This is a rough hack, make this right later!
        if number_of_swing_foot_fz_zeros_stretch < swing_foot_zero_stretch_length_threshold_indices
            removal_flags(i_stretch) = 1;
            disp(['excluding stretch index ' num2str(i_stretch) ' - cross over at time ' num2str(time_forceplate(stretch_start_index_forceplate))]);
        end
        if number_of_swing_foot_fz_zeros_prelude < swing_foot_zero_stretch_length_threshold_indices
            removal_flags(i_stretch) = 1;
            disp(['excluding stretch index ' num2str(i_stretch) ' - cross over at time ' num2str(time_forceplate(stretch_start_index_forceplate))]);
        end
        if isempty(touchdown_index_forceplate)
            touchdown_index_forceplate = stretch_start_index_forceplate;
            touchdown_time_forceplate = time_forceplate(touchdown_index_forceplate);
            removal_flags(i_stretch) = 1;
        end

        % time
        if ignore_crossover
            start_index_forceplate = stretch_start_index_forceplate;
        else
            start_index_forceplate = touchdown_index_forceplate;
        end
        end_index_forceplate = stretch_end_indices_forceplate(i_stretch);
        [~, start_index_mocap] = min(abs(time_mocap - time_forceplate(start_index_forceplate)));
        [~, end_index_mocap] = min(abs(time_mocap - time_forceplate(end_index_forceplate)));
        [~, start_index_labview] = min(abs(time_labview - time_forceplate(start_index_forceplate)));
        [~, end_index_labview] = min(abs(time_labview - time_forceplate(end_index_forceplate)));
        if emg_present
            [~, start_index_emg] = min(abs(time_emg - time_forceplate(start_index_forceplate)));
            [~, end_index_emg] = min(abs(time_emg - time_forceplate(end_index_forceplate)));
        end
        
        % store
        start_indices_mocap(i_stretch) = start_index_mocap;
        end_indices_mocap(i_stretch) = end_index_mocap;
        start_indices_forceplate(i_stretch) = start_index_forceplate;
        end_indices_forceplate(i_stretch) = end_index_forceplate;
        start_indices_labview(i_stretch) = start_index_labview;
        end_indices_labview(i_stretch) = end_index_labview;
        if emg_present
            start_indices_emg(i_stretch) = start_index_emg;
            end_indices_emg(i_stretch) = end_index_emg;
        end
        if ~strcmp(stimulus_type, 'none')
            if stim_start_indices_labview(i_stretch) == 0
                stim_start_time_relative_to_stretch(i_stretch) = NaN;
            else
                stim_start_time = time_forceplate(stim_start_indices_labview(i_stretch));
                stim_start_time_relative_to_stretch(i_stretch) = stim_start_time - touchdown_time_forceplate;
            end
        end
        

        % extract force plate data
        left_cop_x_extracted_stretch = left_copx_trajectory(start_index_forceplate : end_index_forceplate);
        right_cop_x_extracted_stretch = right_copx_trajectory(start_index_forceplate : end_index_forceplate);

        % check for data correctness - one force plate data stretch should have no zeros
        if any(left_cop_x_extracted_stretch==0) & any(right_cop_x_extracted_stretch==0)
            removal_flags(i_stretch) = 1;
        end
    end

    % remove flagged stretches
    unflagged_indices = ~removal_flags;
    condition_stance_foot_list = condition_stance_foot_list(unflagged_indices);
    condition_polarity_list = condition_polarity_list(unflagged_indices);
    condition_delay_list = condition_delay_list(unflagged_indices);
    start_indices_mocap = start_indices_mocap(unflagged_indices);
    end_indices_mocap = end_indices_mocap(unflagged_indices);
    start_indices_forceplate = start_indices_forceplate(unflagged_indices);
    end_indices_forceplate = end_indices_forceplate(unflagged_indices);
    start_indices_emg = start_indices_emg(unflagged_indices);
    end_indices_emg = end_indices_emg(unflagged_indices);
    stim_start_time_relative_to_stretch = stim_start_time_relative_to_stretch(unflagged_indices);


    % prepare stretches for further processing
    data_points_to_process_mocap = [];
    number_of_padding_steps = 30;
    for i_stretch = 1 : size(start_indices_mocap)
        % get start and end indices and apply padding
        start_index_mocap_stretch = start_indices_mocap(i_stretch) - number_of_padding_steps;
        end_index_mocap_stretch = end_indices_mocap(i_stretch) + number_of_padding_steps;
        
        % crop if necessary
        if start_index_mocap_stretch < 1
            start_index_mocap_stretch = 1;
        end
        if end_index_mocap_stretch > size(marker_trajectories, 1);
            end_index_mocap_stretch = size(marker_trajectories, 1);
        end
        
        time_steps_stretch = start_index_mocap_stretch : end_index_mocap_stretch;
        data_points_to_process_mocap = [data_points_to_process_mocap time_steps_stretch];
    end
    
    %% save data
    file_name = makeFileName(date, subject_id, 'walking', i_trial, 'relevantDataStretches');
    save ...
      ( ...
        file_name, ...
        'condition_stance_foot_list', ...
        'condition_polarity_list', ...
        'condition_delay_list', ...
        'start_indices_mocap', ...
        'end_indices_mocap', ...
        'start_indices_forceplate', ...
        'end_indices_forceplate', ...
        'start_indices_labview', ...
        'end_indices_labview', ...
        'start_indices_emg', ...
        'end_indices_emg', ...
        'stim_start_time_relative_to_stretch', ...
        'data_points_to_process_mocap' ...
      )

    disp(['Trial ' num2str(i_trial) ' completed']);

end

























