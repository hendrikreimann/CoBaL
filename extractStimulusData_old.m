
%% Choose Data Type
visual_perturbation             = 1;
phase_dependent                 = 0;
emg_present                     = 0;

%% Choose Analysis Processes
extract_data                    = 0;
calculate_responses             = 0;
calculate_stats                 = 0;
calculate_response_extrema      = 0;

view_totals_and_removals        = 0;
visualize_triggers              = 0;
visualize_steps_during_extract  = 0;
do_cop_plots_absolute_right     = 0;
do_cop_plots_absolute_left      = 0;
do_heel_plots_absolute_right    = 0;
do_heel_plots_absolute_left     = 0;
do_emg_plots_absolute           = 0;
do_cop_plots_response_right     = 0;
do_cop_plots_response_left      = 0;
do_heel_plots_response_right    = 0;
do_heel_plots_response_left     = 0;
do_side_comparison_plots        = 0;
do_response_extrema_plots       = 0;
do_step_response_plot           = 0;
do_stim_response_plot           = 0;
do_stim_start_time_histograms   = 0;

%% prepare
wait_times = [0 0.150 0.450];
wait_time_labels = {'0ms', '150ms', '450ms'};
load subjectInfo.mat;

% trials_to_process = 1 : 23;
trials_to_process = [1:7 9:23];
% trials_to_process = 12 : 23;

total_positive_steps = [];
total_negative_steps = [];
total_polarity_condition_list = [];
total_step_condition_list = [];

number_of_time_steps_normalized = 256;
swing_foot_fz_zero_threshold = 20; % threshold for counting a vertical force reading as zero, in Nm
swing_foot_zero_stretch_length_threshold = 30; % the number of zero indices in the vertical swing foot force has to be larger than this number
duration_until_nearest_future_heelstrike_threshold = 0.1; % a heelstrike should happen less than this long after a trigger

lasi_marker = 24;
rasi_marker = 25;
lpsi_marker = 26;
rpsi_marker = 27;
left_heel_marker = 34;
left_toes_marker = 35;
right_heel_marker = 42;
right_toes_marker = 43;

lasi_marker_indices = reshape([(lasi_marker - 1) * 3 + 1; (lasi_marker - 1) * 3 + 2; (lasi_marker - 1) * 3 + 3], 1, length(lasi_marker)*3);
rasi_marker_indices = reshape([(rasi_marker - 1) * 3 + 1; (rasi_marker - 1) * 3 + 2; (rasi_marker - 1) * 3 + 3], 1, length(rasi_marker)*3);
lpsi_marker_indices = reshape([(lpsi_marker - 1) * 3 + 1; (lpsi_marker - 1) * 3 + 2; (lpsi_marker - 1) * 3 + 3], 1, length(lpsi_marker)*3);
rpsi_marker_indices = reshape([(rpsi_marker - 1) * 3 + 1; (rpsi_marker - 1) * 3 + 2; (rpsi_marker - 1) * 3 + 3], 1, length(rpsi_marker)*3);
left_heel_marker_indices = reshape([(left_heel_marker - 1) * 3 + 1; (left_heel_marker - 1) * 3 + 2; (left_heel_marker - 1) * 3 + 3], 1, length(left_heel_marker)*3);
left_toes_marker_indices = reshape([(left_toes_marker - 1) * 3 + 1; (left_toes_marker - 1) * 3 + 2; (left_toes_marker - 1) * 3 + 3], 1, length(left_toes_marker)*3);
right_heel_marker_indices = reshape([(right_heel_marker - 1) * 3 + 1; (right_heel_marker - 1) * 3 + 2; (right_heel_marker - 1) * 3 + 3], 1, length(right_heel_marker)*3);
right_toes_marker_indices = reshape([(right_toes_marker - 1) * 3 + 1; (right_toes_marker - 1) * 3 + 2; (right_toes_marker - 1) * 3 + 3], 1, length(right_toes_marker)*3);

if emg_present == 1
    left_glutmed_sensor_index = 7;
    left_tibiant_sensor_index = 8;
    left_perolng_sensor_index = 9;
    right_glutmed_sensor_index = 1;
    right_tibiant_sensor_index = 2;
    right_perolng_sensor_index = 3;
end

color_left_control = [0.3 0.1 1];
color_left_positive = [1 0.3 .1] * 0.7;
color_left_negative = [0.3 1 0.1] * 0.7;
color_right_control = [0.3 0.1 1];
color_right_positive = [1 0.3 .1] * 0.7;
color_right_negative = [0.3 1 0.1] * 0.7;

%% extract data
if extract_data
    stretch_length_indices_forceplate = [];
    condition_stance_foot_list = {};
    condition_polarity_list = {};
    condition_delay_list = {};
    left_cop_x_normalized = [];
    right_cop_x_normalized = [];
    lasi_x_pos_normalized = [];
    rasi_x_pos_normalized = [];
    lpsi_x_pos_normalized = [];
    rpsi_x_pos_normalized = [];
    lasi_x_vel_normalized = [];
    rasi_x_vel_normalized = [];
    lpsi_x_vel_normalized = [];
    rpsi_x_vel_normalized = [];
    pelvis_x_pos_normalized = [];
    pelvis_x_vel_normalized = [];
    left_heel_x_pos_normalized = [];
    right_heel_x_pos_normalized = [];
    left_heel_y_pos_normalized = [];
    right_heel_y_pos_normalized = [];
    
    if emg_present
        left_glutmed_emg_normalized = [];
        left_tibiant_emg_normalized = [];
        left_perolng_emg_normalized = [];
        right_glutmed_emg_normalized = [];
        right_tibiant_emg_normalized = [];
        right_perolng_emg_normalized = [];
    end
    
    step_times = [];
    stim_start_time_relative_to_stretch = [];
    
    for i_trial = trials_to_process
        % load data
        load(makeFileName(date, subject_id, 'walking', i_trial, 'forcePlateData'));
        load(makeFileName(date, subject_id, 'walking', i_trial, 'markerTrajectories'));
        if emg_present
            load(makeFileName(date, subject_id, 'walking', i_trial, 'emgTrajectories'));
        end
        load(makeFileName(date, subject_id, 'walking', i_trial, 'stepEvents'));
        
        if visual_perturbation
            stim_sent_trajectory = visual_shift_ml_trajectory;
        end
        
        left_copx_trajectory = left_force_plate_cop_Acw(:, 1);
        right_copx_trajectory = right_force_plate_cop_Acw(:, 1);
        left_fz_trajectory = left_force_plate_wrench_Acw(:, 3);
        right_fz_trajectory = right_force_plate_wrench_Acw(:, 3);
        lasi_x_pos_trajectory = marker_trajectories(:, lasi_marker_indices(1));
        rasi_x_pos_trajectory = marker_trajectories(:, rasi_marker_indices(1));
        lpsi_x_pos_trajectory = marker_trajectories(:, lpsi_marker_indices(1));
        rpsi_x_pos_trajectory = marker_trajectories(:, rpsi_marker_indices(1));
        lasi_x_vel_trajectory = deriveByTime(lasi_x_pos_trajectory, sampling_rate_mocap^(-1));
        rasi_x_vel_trajectory = deriveByTime(rasi_x_pos_trajectory, sampling_rate_mocap^(-1));
        lpsi_x_vel_trajectory = deriveByTime(lpsi_x_pos_trajectory, sampling_rate_mocap^(-1));
        rpsi_x_vel_trajectory = deriveByTime(rpsi_x_pos_trajectory, sampling_rate_mocap^(-1));
        left_heel_x_pos_trajectory = marker_trajectories(:, left_heel_marker_indices(1));
        right_heel_x_pos_trajectory = marker_trajectories(:, right_heel_marker_indices(1));
        left_heel_y_pos_trajectory = marker_trajectories(:, left_heel_marker_indices(2));
        right_heel_y_pos_trajectory = marker_trajectories(:, right_heel_marker_indices(2));
        
        if emg_present
            left_glutmed_emg = emg_trajectories(:, left_glutmed_sensor_index);
            left_tibiant_emg = emg_trajectories(:, left_tibiant_sensor_index);
            left_perolng_emg = emg_trajectories(:, left_perolng_sensor_index);
            right_glutmed_emg = emg_trajectories(:, right_glutmed_sensor_index);
            right_tibiant_emg = emg_trajectories(:, right_tibiant_sensor_index);
            right_perolng_emg = emg_trajectories(:, right_perolng_sensor_index);
        end
        
        % find trigger
        trigger_indices_forceplate_trial = find(diff(sign(stimulus_foot_state - 0.5)) > 0) + 1;
        epsilon = 1e-5;
        stim_start_indices_forceplate_trial = find(diff(sign(abs(stim_sent_trajectory) - epsilon)) > 0) + 1;
        trigger_indices_forceplate_trial = trigger_indices_forceplate_trial(1 : length(stim_start_indices_forceplate_trial)); % in case a stim is triggered, but not recorded
        % there might be a problem here where the stim start index is identified as the last time step before the stim
        
        % visualize triggers
        if visualize_triggers
            left_cop_x_trajectory_relevant = left_copx_trajectory; left_cop_x_trajectory_relevant(left_cop_x_trajectory_relevant==0) = NaN;
            right_cop_x_trajectory_relevant = right_copx_trajectory; right_cop_x_trajectory_relevant(right_cop_x_trajectory_relevant==0) = NaN;
            figure; axes; hold on
            plot(time_force_plate, stimulus_foot_state);
            plot(time_force_plate, left_cop_x_trajectory_relevant, 'linewidth', 2);
            plot(time_force_plate, right_cop_x_trajectory_relevant, 'linewidth', 2);
            plot(time_force_plate(left_touchdown_indices_force_plate), zeros(size(left_touchdown_indices_force_plate)), 'o')
            plot(time_force_plate(right_touchdown_indices_force_plate), zeros(size(right_touchdown_indices_force_plate)), 'o')
            plot(time_force_plate(trigger_indices_forceplate_trial), zeros(size(trigger_indices_forceplate_trial)), 'x')
            plot(time_force_plate(stim_start_indices_forceplate_trial), zeros(size(stim_start_indices_forceplate_trial)), 'x')
            legend('stimulus state', 'left cop', 'right cop', 'left touchdown', 'right touchdown', 'trigger', 'stim start')
        end
        
        % for each trigger, extract conditions and relevant step events
        number_of_triggers = length(trigger_indices_forceplate_trial);
        removal_flags = zeros(number_of_triggers, 1);
        stim_start_indices_forceplate_trial = [stim_start_indices_forceplate_trial zeros(number_of_triggers, 1)]; %#ok<AGROW>
        stretch_start_indices_forceplate_trial = zeros(number_of_triggers, 2);
        stretch_end_indices_forceplate_trial = zeros(number_of_triggers, 2);
        duration_until_nearest_future_heelstrike_trial = zeros(number_of_triggers, 2);
        condition_stance_foot_list_trial = cell(number_of_triggers, 2);
        condition_polarity_list_trial = cell(number_of_triggers, 2);
        condition_delay_list_trial = cell(number_of_triggers, 2);

        %% Check Each Stim Time: Proper Trigger and FP Data
        for i_trigger = 1 : number_of_triggers
            % polarity condition
            if stim_sent_trajectory(stim_start_indices_forceplate_trial(i_trigger)) > 0
                condition_polarity_list_trial{i_trigger, 1} = 'POSITIVE';
            elseif stim_sent_trajectory(stim_start_indices_forceplate_trial(i_trigger)) < 0
                condition_polarity_list_trial{i_trigger, 1} = 'NEGATIVE';
            else
                disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate_trial(i_trigger))) ' - no stim']);
            end
            condition_polarity_list_trial{i_trigger, 2} = 'CONTROL';
            
            % delay
            wait_time_stim = time_force_plate(stim_start_indices_forceplate_trial(i_trigger)) - time_force_plate(trigger_indices_forceplate_trial(i_trigger));
            [~, wait_condition_index] = min(abs(wait_times - wait_time_stim));
            condition_delay_list_trial{i_trigger, 1} = wait_time_labels{wait_condition_index};
            condition_delay_list_trial{i_trigger, 2} = 'CONTROL';
            
            %% check which foot is the stance foot for time period of interest
            [last_left_foot_heelstrike, index_left] = max(left_touchdown_indices_force_plate(left_touchdown_indices_force_plate <= trigger_indices_forceplate_trial(i_trigger)));
            if isempty(index_left)
                removal_flags(i_trigger) = 1;
            else
                if length(left_touchdown_indices_force_plate) < index_left + 2
                    removal_flags(i_trigger) = 1;
                else
                    this_left_foot_heelstrike = left_touchdown_indices_force_plate(index_left+1);
                    next_left_foot_heelstrike = left_touchdown_indices_force_plate(index_left+2);
                end
            end
            [last_right_foot_heelstrike, index_right] = max(right_touchdown_indices_force_plate(right_touchdown_indices_force_plate <= trigger_indices_forceplate_trial(i_trigger)));
            if isempty(index_right)
                removal_flags(i_trigger) = 1;
            else
                if length(right_touchdown_indices_force_plate) < index_right + 2
                    removal_flags(i_trigger) = 1;
                else
                    this_right_foot_heelstrike = right_touchdown_indices_force_plate(index_right+1);
                    next_right_foot_heelstrike = right_touchdown_indices_force_plate(index_right+2);
                end
            end
            
            % confirm that the trigger happened less than 100ms before the actual heelstrike
            nearest_future_heelstrike_index = min([this_left_foot_heelstrike this_right_foot_heelstrike]);
            nearest_future_heelstrike_time = time_force_plate(nearest_future_heelstrike_index);
            duration_until_nearest_future_heelstrike_trial(i_trigger, 1) = nearest_future_heelstrike_time - time_force_plate(trigger_indices_forceplate_trial(i_trigger));
            if duration_until_nearest_future_heelstrike_trial(i_trigger, 1) > duration_until_nearest_future_heelstrike_threshold
                removal_flags(i_trigger) = 1;
                disp(['XXXXXXX Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate_trial(i_trigger))) ' - trigger not in sync with heelstrike']);
            end
            
            % mark relevant event delimiters depending on wait time and triggering foot
            if left_contact_indicators_force_plate(trigger_indices_forceplate_trial(i_trigger)) == 1 && right_contact_indicators_force_plate(trigger_indices_forceplate_trial(i_trigger)) == 0
                % triggered by right foot heelstrike
                if strcmp(condition_delay_list_trial{i_trigger, 1}, '0ms') || strcmp(condition_delay_list_trial{i_trigger, 1}, '150ms')
                    % this step is of interest
                    condition_stance_foot_list_trial{i_trigger, 1} = 'RIGHT';
                    condition_stance_foot_list_trial{i_trigger, 2} = 'RIGHT';
                    stretch_start_indices_forceplate_trial(i_trigger, 1) = this_right_foot_heelstrike;
                    stretch_end_indices_forceplate_trial(i_trigger, 1) = this_left_foot_heelstrike;
                    stretch_start_indices_forceplate_trial(i_trigger, 2) = last_right_foot_heelstrike;
                    stretch_end_indices_forceplate_trial(i_trigger, 2) = last_left_foot_heelstrike;
                elseif strcmp(condition_delay_list_trial{i_trigger, 1}, '450ms')
                    % next step is of interest
                    condition_stance_foot_list_trial{i_trigger, 1} = 'LEFT';
                    condition_stance_foot_list_trial{i_trigger, 2} = 'LEFT';
                    stretch_start_indices_forceplate_trial(i_trigger, 1) = this_left_foot_heelstrike;
                    stretch_end_indices_forceplate_trial(i_trigger, 1) = next_right_foot_heelstrike;
                    stretch_start_indices_forceplate_trial(i_trigger, 2) = last_left_foot_heelstrike;
                    stretch_end_indices_forceplate_trial(i_trigger, 2) = this_right_foot_heelstrike;
                else
                    disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate_trial(i_trigger))) ' - delay not defined']);
                end
            elseif left_contact_indicators_force_plate(trigger_indices_forceplate_trial(i_trigger)) == 0 && right_contact_indicators_force_plate(trigger_indices_forceplate_trial(i_trigger)) == 1
                % triggered by left foot heelstrike
                if strcmp(condition_delay_list_trial{i_trigger, 1}, '0ms') || strcmp(condition_delay_list_trial{i_trigger, 1}, '150ms')
                    % this step is of interest
                    condition_stance_foot_list_trial{i_trigger, 1} = 'LEFT';
                    condition_stance_foot_list_trial{i_trigger, 2} = 'LEFT';
                    stretch_start_indices_forceplate_trial(i_trigger, 1) = this_left_foot_heelstrike;
                    stretch_end_indices_forceplate_trial(i_trigger, 1) = this_right_foot_heelstrike;
                    stretch_start_indices_forceplate_trial(i_trigger, 2) = last_left_foot_heelstrike;
                    stretch_end_indices_forceplate_trial(i_trigger, 2) = last_right_foot_heelstrike;
                elseif strcmp(condition_delay_list_trial{i_trigger, 1}, '450ms')
                    % next step is of interest
                    condition_stance_foot_list_trial{i_trigger, 1} = 'RIGHT';
                    condition_stance_foot_list_trial{i_trigger, 2} = 'RIGHT';
                    stretch_start_indices_forceplate_trial(i_trigger, 1) = this_right_foot_heelstrike;
                    stretch_end_indices_forceplate_trial(i_trigger, 1) = next_left_foot_heelstrike;
                    stretch_start_indices_forceplate_trial(i_trigger, 2) = last_right_foot_heelstrike;
                    stretch_end_indices_forceplate_trial(i_trigger, 2) = this_left_foot_heelstrike;
                else
                    disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate_trial(i_trigger))) ' - delay not defined']);
                end
            else
                removal_flags(i_trigger) = 1;
                condition_stance_foot_list_trial{i_trigger, 1} = 'UNCLEAR';
                condition_stance_foot_list_trial{i_trigger, 2} = 'UNCLEAR';
                disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate_trial(i_trigger))) ' - stance foot unclear']);
            end
            
        end

        % remove flagged stims
        unflagged_indices = ~removal_flags;
        stim_start_indices_forceplate_trial = stim_start_indices_forceplate_trial(unflagged_indices, :);
        stretch_start_indices_forceplate_trial = stretch_start_indices_forceplate_trial(unflagged_indices, :);
        stretch_end_indices_forceplate_trial = stretch_end_indices_forceplate_trial(unflagged_indices, :);
        condition_stance_foot_list_trial = condition_stance_foot_list_trial(unflagged_indices, :);
        condition_polarity_list_trial = condition_polarity_list_trial(unflagged_indices, :);
        condition_delay_list_trial = condition_delay_list_trial(unflagged_indices, :);
        duration_until_nearest_future_heelstrike_trial = duration_until_nearest_future_heelstrike_trial(unflagged_indices, :);

        % form stretches
        stim_start_indices_forceplate_trial = [stim_start_indices_forceplate_trial(:, 1); stim_start_indices_forceplate_trial(:, 2)];
        stretch_start_indices_forceplate_trial = [stretch_start_indices_forceplate_trial(:, 1); stretch_start_indices_forceplate_trial(:, 2)];
        stretch_end_indices_forceplate_trial = [stretch_end_indices_forceplate_trial(:, 1); stretch_end_indices_forceplate_trial(:, 2)];
        condition_stance_foot_list_trial = [condition_stance_foot_list_trial(:, 1); condition_stance_foot_list_trial(:, 2)];
        condition_polarity_list_trial = [condition_polarity_list_trial(:, 1); condition_polarity_list_trial(:, 2)];
        condition_delay_list_trial = [condition_delay_list_trial(:, 1); condition_delay_list_trial(:, 2)];

        % extract and time-normalize data
        number_of_stretches_trial = length(stretch_start_indices_forceplate_trial);
        step_times_trial = zeros(1, number_of_stretches_trial);
        number_of_time_steps_extracted = zeros(1, number_of_stretches_trial);

        left_cop_x_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        right_cop_x_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        
        lasi_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        rasi_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        lpsi_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        rpsi_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        lasi_x_vel_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        rasi_x_vel_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        lpsi_x_vel_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        rpsi_x_vel_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        pelvis_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        pelvis_x_vel_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        left_heel_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        right_heel_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        left_heel_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        right_heel_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        
        if emg_present == 1
            left_glutmed_emg_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            left_tibiant_emg_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            left_perolng_emg_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            right_glutmed_emg_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            right_tibiant_emg_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            right_perolng_emg_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        end
        
        stim_start_time_relative_to_stretch_trial = zeros(1, number_of_stretches_trial);
        number_of_indices_per_step_trial = zeros(1, number_of_stretches_trial);
        removal_flags = zeros(number_of_stretches_trial, 1);
        for i_stretch = 1 : number_of_stretches_trial
            % we have the delimiters for the stretches of interest, based on touchdown and pushoff identification from mocap
            % now look at the force plate data and decide whether this stretch is usable or not due to cross-over
            % if the stretch is usable, find the touchdown and pushoff events from force plate loading
            
            stretch_start_index_forceplate = stretch_start_indices_forceplate_trial(i_stretch);
            stretch_end_index_forceplate = stretch_end_indices_forceplate_trial(i_stretch);
            
            % find start of CoP data
            if strcmp(condition_stance_foot_list_trial{i_stretch}, 'RIGHT')
                stance_foot_cop_stretch = right_copx_trajectory(stretch_start_index_forceplate : stretch_end_index_forceplate);
                touchdown_index_forceplate = stretch_start_index_forceplate + find(stance_foot_cop_stretch~=0, 1, 'first') - 1;
                stance_foot_fz_step = right_fz_trajectory(touchdown_index_forceplate : stretch_end_index_forceplate);
                swing_foot_fz_step = left_fz_trajectory(touchdown_index_forceplate : stretch_end_index_forceplate);
            elseif strcmp(condition_stance_foot_list_trial{i_stretch}, 'LEFT')
                stance_foot_cop_stretch = left_copx_trajectory(stretch_start_index_forceplate : stretch_end_index_forceplate);
                touchdown_index_forceplate = stretch_start_index_forceplate + find(stance_foot_cop_stretch~=0, 1, 'first') - 1;
                stance_foot_fz_step = left_fz_trajectory(touchdown_index_forceplate : stretch_end_index_forceplate);
                swing_foot_fz_step = right_fz_trajectory(touchdown_index_forceplate : stretch_end_index_forceplate);
            end
            touchdown_time_forceplate = time_force_plate(touchdown_index_forceplate);
            time_extracted_forceplate = time_force_plate(touchdown_index_forceplate : stretch_end_index_forceplate);
            
            % check if this stretch is usable
            number_of_indices_per_step_trial(i_stretch) = length(time_extracted_forceplate);
            number_of_swing_foot_fz_zeros_stretch = sum(-swing_foot_fz_step < swing_foot_fz_zero_threshold);
            if number_of_swing_foot_fz_zeros_stretch < swing_foot_zero_stretch_length_threshold
                removal_flags(i_stretch) = 1;
                disp(['excluding stretch index ' num2str(i_stretch) ' - cross over']);
            end
            if isempty(touchdown_index_forceplate)
                touchdown_index_forceplate = stretch_start_index_forceplate;
                removal_flags(i_stretch) = 1;
            end

%             figure; axes; hold on; title('left')
%             plot(time_extracted_forceplate, stance_foot_fz_step)
%             plot(time_extracted_forceplate, swing_foot_fz_step)
% %             plot(time_extracted_forceplate, 
%             plot(touchdown_time_forceplate, 0, 'o')
%             legend('stance cop', 'swing cop', 'touchdown')
            
            
            
            
            
            
            
            
            
            
            
            
            % time
            start_index_force_plate = touchdown_index_forceplate;
            end_index_force_plate = stretch_end_indices_forceplate_trial(i_stretch);
            [~, start_index_mocap] = min(abs(time_mocap - time_force_plate(start_index_force_plate)));
            [~, end_index_mocap] = min(abs(time_mocap - time_force_plate(end_index_force_plate)));
            
            if emg_present
                [~, start_index_emg] = min(abs(time_emg - time_force_plate(start_index_force_plate)));
                [~, end_index_emg] = min(abs(time_emg - time_force_plate(end_index_force_plate)));
            end
            
            step_times_trial(i_stretch) = time_force_plate(end_index_force_plate) - time_force_plate(start_index_force_plate);
            
            if stim_start_indices_forceplate_trial(i_stretch) == 0
                stim_start_time_relative_to_stretch_trial(i_stretch) = NaN;
            else
                stim_start_time = time_force_plate(stim_start_indices_forceplate_trial(i_stretch));
                stim_start_time_relative_to_stretch_trial(i_stretch) = stim_start_time - touchdown_time_forceplate;
            end
            


            % extract force plate data and normalize time
            left_cop_x_extracted_stretch = left_copx_trajectory(start_index_force_plate : end_index_force_plate);
            right_cop_x_extracted_stretch = right_copx_trajectory(start_index_force_plate : end_index_force_plate);
            time_extracted_forceplate = time_force_plate(start_index_force_plate : end_index_force_plate);
            time_normalized_forceplate = linspace(time_extracted_forceplate(1), time_extracted_forceplate(end), number_of_time_steps_normalized);
            left_cop_x_normalized_stretch = spline(time_extracted_forceplate, left_cop_x_extracted_stretch, time_normalized_forceplate);
            right_cop_x_normalized_stretch = spline(time_extracted_forceplate, right_cop_x_extracted_stretch, time_normalized_forceplate);

            % check for data correctness - one force plate data stretch should have no zeros
            if any(left_cop_x_extracted_stretch==0) & any(right_cop_x_extracted_stretch==0)
                removal_flags(i_stretch) = 1;
            end

            % extract mocap data
            time_extracted_mocap = time_mocap(start_index_mocap : end_index_mocap);
            lasi_x_pos_extracted_stretch = lasi_x_pos_trajectory(start_index_mocap : end_index_mocap);
            rasi_x_pos_extracted_stretch = rasi_x_pos_trajectory(start_index_mocap : end_index_mocap);
            lpsi_x_pos_extracted_stretch = lpsi_x_pos_trajectory(start_index_mocap : end_index_mocap);
            rpsi_x_pos_extracted_stretch = rpsi_x_pos_trajectory(start_index_mocap : end_index_mocap);
            lasi_x_vel_extracted_stretch = lasi_x_vel_trajectory(start_index_mocap : end_index_mocap);
            rasi_x_vel_extracted_stretch = rasi_x_vel_trajectory(start_index_mocap : end_index_mocap);
            lpsi_x_vel_extracted_stretch = lpsi_x_vel_trajectory(start_index_mocap : end_index_mocap);
            rpsi_x_vel_extracted_stretch = rpsi_x_vel_trajectory(start_index_mocap : end_index_mocap);
            left_heel_x_pos_extracted_stretch = left_heel_x_pos_trajectory(start_index_mocap : end_index_mocap);
            right_heel_x_pos_extracted_stretch = right_heel_x_pos_trajectory(start_index_mocap : end_index_mocap);
            left_heel_y_pos_extracted_stretch = left_heel_y_pos_trajectory(start_index_mocap : end_index_mocap);
            right_heel_y_pos_extracted_stretch = right_heel_y_pos_trajectory(start_index_mocap : end_index_mocap);
            
            % normalize mocap data in time
            time_normalized_mocap = linspace(time_extracted_mocap(1), time_extracted_mocap(end), number_of_time_steps_normalized);
%             if any(any(isnan([lasi_x_pos_extracted_stretch; rasi_x_pos_extracted_stretch; lpsi_x_pos_extracted_stretch; rpsi_x_pos_extracted_stretch])))
%                 removal_flags(i_stretch) = 1;
%             else
            lasi_x_pos_normalized_stretch = spline(time_extracted_mocap, lasi_x_pos_extracted_stretch, time_normalized_mocap);
            rasi_x_pos_normalized_stretch = spline(time_extracted_mocap, rasi_x_pos_extracted_stretch, time_normalized_mocap);
            lpsi_x_pos_normalized_stretch = spline(time_extracted_mocap, lpsi_x_pos_extracted_stretch, time_normalized_mocap);
            rpsi_x_pos_normalized_stretch = spline(time_extracted_mocap, rpsi_x_pos_extracted_stretch, time_normalized_mocap);
            lasi_x_vel_normalized_stretch = spline(time_extracted_mocap, lasi_x_vel_extracted_stretch, time_normalized_mocap);
            rasi_x_vel_normalized_stretch = spline(time_extracted_mocap, rasi_x_vel_extracted_stretch, time_normalized_mocap);
            lpsi_x_vel_normalized_stretch = spline(time_extracted_mocap, lpsi_x_vel_extracted_stretch, time_normalized_mocap);
            rpsi_x_vel_normalized_stretch = spline(time_extracted_mocap, rpsi_x_vel_extracted_stretch, time_normalized_mocap);
            left_heel_x_pos_normalized_stretch = spline(time_extracted_mocap, left_heel_x_pos_extracted_stretch, time_normalized_mocap);
            right_heel_x_pos_normalized_stretch = spline(time_extracted_mocap, right_heel_x_pos_extracted_stretch, time_normalized_mocap);
            left_heel_y_pos_normalized_stretch = spline(time_extracted_mocap, left_heel_y_pos_extracted_stretch, time_normalized_mocap);
            right_heel_y_pos_normalized_stretch = spline(time_extracted_mocap, right_heel_y_pos_extracted_stretch, time_normalized_mocap);
%             end
            % normalize in space to stance foot heel at start
            if strcmp(condition_stance_foot_list_trial{i_stretch}, 'RIGHT')
                stance_foot_heel_x_initial = right_heel_x_pos_extracted_stretch(1);
                stance_foot_heel_y_initial = right_heel_y_pos_extracted_stretch(1);
            elseif strcmp(condition_stance_foot_list_trial{i_stretch}, 'LEFT')
                stance_foot_heel_x_initial = left_heel_x_pos_extracted_stretch(1);
                stance_foot_heel_y_initial = left_heel_y_pos_extracted_stretch(1);
            end
            
            if emg_present
                % extract emg data
                time_extracted_emg = time_emg(start_index_emg : end_index_emg);
                left_glutmed_emg_extracted_stretch = left_glutmed_emg(start_index_emg : end_index_emg);
                left_tibiant_emg_extracted_stretch = left_tibiant_emg(start_index_emg : end_index_emg);
                left_perolng_emg_extracted_stretch = left_perolng_emg(start_index_emg : end_index_emg);
                right_glutmed_emg_extracted_stretch = right_glutmed_emg(start_index_emg : end_index_emg);
                right_tibiant_emg_extracted_stretch = right_tibiant_emg(start_index_emg : end_index_emg);
                right_perolng_emg_extracted_stretch = right_perolng_emg(start_index_emg : end_index_emg);

                % normalize emg data in time
                time_normalized_emg = linspace(time_extracted_emg(1), time_extracted_emg(end), number_of_time_steps_normalized);
                left_glutmed_emg_normalized_stretch = spline(time_extracted_emg, left_glutmed_emg_extracted_stretch, time_normalized_emg);
                left_tibiant_emg_normalized_stretch = spline(time_extracted_emg, left_tibiant_emg_extracted_stretch, time_normalized_emg);
                left_perolng_emg_normalized_stretch = spline(time_extracted_emg, left_perolng_emg_extracted_stretch, time_normalized_emg);
                right_glutmed_emg_normalized_stretch = spline(time_extracted_emg, right_glutmed_emg_extracted_stretch, time_normalized_emg);
                right_tibiant_emg_normalized_stretch = spline(time_extracted_emg, right_tibiant_emg_extracted_stretch, time_normalized_emg);
                right_perolng_emg_normalized_stretch = spline(time_extracted_emg, right_perolng_emg_extracted_stretch, time_normalized_emg);
            end
            
            % store
            % XXX ACHTUNG! this mingles normalization with storing. Separate that
            left_cop_x_normalized_trial(:, i_stretch) = left_cop_x_normalized_stretch - stance_foot_heel_x_initial;
            right_cop_x_normalized_trial(:, i_stretch) = right_cop_x_normalized_stretch - stance_foot_heel_x_initial;
            lasi_x_pos_normalized_trial(:, i_stretch) = lasi_x_pos_normalized_stretch - stance_foot_heel_x_initial;
            rasi_x_pos_normalized_trial(:, i_stretch) = rasi_x_pos_normalized_stretch - stance_foot_heel_x_initial;
            lpsi_x_pos_normalized_trial(:, i_stretch) = lpsi_x_pos_normalized_stretch - stance_foot_heel_x_initial;
            rpsi_x_pos_normalized_trial(:, i_stretch) = rpsi_x_pos_normalized_stretch - stance_foot_heel_x_initial;
            lasi_x_vel_normalized_trial(:, i_stretch) = lasi_x_vel_normalized_stretch;
            rasi_x_vel_normalized_trial(:, i_stretch) = rasi_x_vel_normalized_stretch;
            lpsi_x_vel_normalized_trial(:, i_stretch) = lpsi_x_vel_normalized_stretch;
            rpsi_x_vel_normalized_trial(:, i_stretch) = rpsi_x_vel_normalized_stretch;
            pelvis_x_pos_normalized_trial(:, i_stretch) = mean([lasi_x_pos_normalized_stretch; rasi_x_pos_normalized_stretch; lpsi_x_pos_normalized_stretch; rpsi_x_pos_normalized_stretch]) - stance_foot_heel_x_initial;
            pelvis_x_vel_normalized_trial(:, i_stretch) = mean([lasi_x_vel_normalized_stretch; rasi_x_vel_normalized_stretch; lpsi_x_vel_normalized_stretch; rpsi_x_vel_normalized_stretch]);
            left_heel_x_pos_normalized_trial(:, i_stretch) = left_heel_x_pos_normalized_stretch - stance_foot_heel_x_initial;
            right_heel_x_pos_normalized_trial(:, i_stretch) = right_heel_x_pos_normalized_stretch - stance_foot_heel_x_initial;
            left_heel_y_pos_normalized_trial(:, i_stretch) = left_heel_y_pos_normalized_stretch - stance_foot_heel_y_initial;
            right_heel_y_pos_normalized_trial(:, i_stretch) = right_heel_y_pos_normalized_stretch - stance_foot_heel_y_initial;
            
            if emg_present
                left_glutmed_emg_normalized_trial(:, i_stretch) = left_glutmed_emg_normalized_stretch;
                left_tibiant_emg_normalized_trial(:, i_stretch) = left_tibiant_emg_normalized_stretch;
                left_perolng_emg_normalized_trial(:, i_stretch) = left_perolng_emg_normalized_stretch;
                right_glutmed_emg_normalized_trial(:, i_stretch) = right_glutmed_emg_normalized_stretch;
                right_tibiant_emg_normalized_trial(:, i_stretch) = right_tibiant_emg_normalized_stretch;
                right_perolng_emg_normalized_trial(:, i_stretch) = right_perolng_emg_normalized_stretch;
            end
            
%             % visualize
%             figure; axes; hold on; title('left')
%             plot(time_extracted_forceplate, left_cop_x_extracted_stretch)
%             plot(time_extracted_forceplate, right_cop_x_extracted_stretch)
%             plot(time_normalized_forceplate, left_cop_x_normalized_stretch)
%             plot(time_normalized_forceplate, right_cop_x_normalized_stretch)
%             plot(time_extracted_emg, left_glutmed_emg_extracted_stretch)
%             plot(time_normalized_emg, left_glutmed_emg_normalized_stretch)
%             legend('left cop extracted', 'right cop extracted', 'left cop normalized', 'right cop normalized', 'left glut med emg extracted', 'left glut med emg normalized')
        end
        
        %% Concat for Totals..
        total_polarity_condition_list = [total_polarity_condition_list; condition_polarity_list_trial];
        total_step_condition_list = [total_step_condition_list; condition_stance_foot_list_trial];
        
        %% 2nd removal flags...??
        % remove flagged stretches
        unflagged_indices = ~removal_flags;
        stretch_start_indices_forceplate_trial = stretch_start_indices_forceplate_trial(unflagged_indices);
        stretch_end_indices_forceplate_trial = stretch_end_indices_forceplate_trial(unflagged_indices);
        condition_stance_foot_list_trial = condition_stance_foot_list_trial(unflagged_indices);
        condition_polarity_list_trial = condition_polarity_list_trial(unflagged_indices);
        condition_delay_list_trial = condition_delay_list_trial(unflagged_indices);
        left_cop_x_normalized_trial = left_cop_x_normalized_trial(:, unflagged_indices);
        right_cop_x_normalized_trial = right_cop_x_normalized_trial(:, unflagged_indices);
        lasi_x_pos_normalized_trial = lasi_x_pos_normalized_trial(:, unflagged_indices);
        rasi_x_pos_normalized_trial = rasi_x_pos_normalized_trial(:, unflagged_indices);
        lpsi_x_pos_normalized_trial = lpsi_x_pos_normalized_trial(:, unflagged_indices);
        rpsi_x_pos_normalized_trial = rpsi_x_pos_normalized_trial(:, unflagged_indices);
        lasi_x_vel_normalized_trial = lasi_x_vel_normalized_trial(:, unflagged_indices);
        rasi_x_vel_normalized_trial = rasi_x_vel_normalized_trial(:, unflagged_indices);
        lpsi_x_vel_normalized_trial = lpsi_x_vel_normalized_trial(:, unflagged_indices);
        rpsi_x_vel_normalized_trial = rpsi_x_vel_normalized_trial(:, unflagged_indices);
        pelvis_x_pos_normalized_trial = pelvis_x_pos_normalized_trial(:, unflagged_indices);
        pelvis_x_vel_normalized_trial = pelvis_x_vel_normalized_trial(:, unflagged_indices);
        left_heel_x_pos_normalized_trial = left_heel_x_pos_normalized_trial(:, unflagged_indices);
        right_heel_x_pos_normalized_trial = right_heel_x_pos_normalized_trial(:, unflagged_indices);
        left_heel_y_pos_normalized_trial = left_heel_y_pos_normalized_trial(:, unflagged_indices);
        right_heel_y_pos_normalized_trial = right_heel_y_pos_normalized_trial(:, unflagged_indices);
        
        if emg_present
            left_glutmed_emg_normalized_trial = left_glutmed_emg_normalized_trial(:, unflagged_indices);
            left_tibiant_emg_normalized_trial = left_tibiant_emg_normalized_trial(:, unflagged_indices);
            left_perolng_emg_normalized_trial = left_perolng_emg_normalized_trial(:, unflagged_indices);
            right_glutmed_emg_normalized_trial = right_glutmed_emg_normalized_trial(:, unflagged_indices);
            right_tibiant_emg_normalized_trial = right_tibiant_emg_normalized_trial(:, unflagged_indices);
            right_perolng_emg_normalized_trial = right_perolng_emg_normalized_trial(:, unflagged_indices);
        end
        
        step_times_trial = step_times_trial(unflagged_indices);
        stim_start_time_relative_to_stretch_trial = stim_start_time_relative_to_stretch_trial(unflagged_indices);
        
        % add data to lists
%         stretch_length_indices_forceplate = [stretch_length_indices_forceplate; stretch_length_indices_forceplate_trial];
        condition_stance_foot_list = [condition_stance_foot_list; condition_stance_foot_list_trial];
        condition_polarity_list = [condition_polarity_list; condition_polarity_list_trial];
        condition_delay_list = [condition_delay_list; condition_delay_list_trial];
        left_cop_x_normalized = [left_cop_x_normalized left_cop_x_normalized_trial];
        right_cop_x_normalized = [right_cop_x_normalized right_cop_x_normalized_trial];
        lasi_x_pos_normalized = [lasi_x_pos_normalized lasi_x_pos_normalized_trial];
        rasi_x_pos_normalized = [rasi_x_pos_normalized rasi_x_pos_normalized_trial];
        lpsi_x_pos_normalized = [lpsi_x_pos_normalized lpsi_x_pos_normalized_trial];
        rpsi_x_pos_normalized = [rpsi_x_pos_normalized rpsi_x_pos_normalized_trial];
        lasi_x_vel_normalized = [lasi_x_vel_normalized lasi_x_vel_normalized_trial];
        rasi_x_vel_normalized = [rasi_x_vel_normalized rasi_x_vel_normalized_trial];
        lpsi_x_vel_normalized = [lpsi_x_vel_normalized lpsi_x_vel_normalized_trial];
        rpsi_x_vel_normalized = [rpsi_x_vel_normalized rpsi_x_vel_normalized_trial];
        pelvis_x_pos_normalized = [pelvis_x_pos_normalized pelvis_x_pos_normalized_trial];
        pelvis_x_vel_normalized = [pelvis_x_vel_normalized pelvis_x_vel_normalized_trial];
        left_heel_x_pos_normalized = [left_heel_x_pos_normalized left_heel_x_pos_normalized_trial];
        right_heel_x_pos_normalized = [right_heel_x_pos_normalized right_heel_x_pos_normalized_trial];
        left_heel_y_pos_normalized = [left_heel_y_pos_normalized left_heel_y_pos_normalized_trial];
        right_heel_y_pos_normalized = [right_heel_y_pos_normalized right_heel_y_pos_normalized_trial];
        
        if emg_present
            left_glutmed_emg_normalized = [left_glutmed_emg_normalized left_glutmed_emg_normalized_trial];
            left_tibiant_emg_normalized = [left_tibiant_emg_normalized left_tibiant_emg_normalized_trial];
            left_perolng_emg_normalized = [left_perolng_emg_normalized left_perolng_emg_normalized_trial];
            right_glutmed_emg_normalized = [right_glutmed_emg_normalized right_glutmed_emg_normalized_trial];
            right_tibiant_emg_normalized = [right_tibiant_emg_normalized right_tibiant_emg_normalized_trial];
            right_perolng_emg_normalized = [right_perolng_emg_normalized right_perolng_emg_normalized_trial];
        end
        
        step_times = [step_times step_times_trial];
        stim_start_time_relative_to_stretch = [stim_start_time_relative_to_stretch stim_start_time_relative_to_stretch_trial];

        % visualize steps
        if visualize_steps_during_extract
            figure; check_axes = axes; hold on
            plot(time_force_plate, heel_strike_count);
            plot(time_force_plate(trigger_indices_forceplate_trial), heel_strike_count(trigger_indices_forceplate_trial), 'o', 'linewidth', 1, 'markersize', 10);
            legend('heel strike count', 'triggered heel strike?')
            
            % left events
            figure; axes_left = axes; hold on; title('Left'); xlim([0, recordingTime]);
            plot(time_force_plate, fzl_trajectory)
            plot(time_mocap(left_touchdown_indices_mocap), zeros(size(left_touchdown_indices_mocap)), 'v', 'linewidth', 1, 'markersize', 10);
            plot(time_mocap(left_pushoff_indices_mocap), zeros(size(left_pushoff_indices_mocap)), '^', 'linewidth', 1, 'markersize', 10);
            plot(time_force_plate(trigger_indices_forceplate_trial), fzl_trajectory(trigger_indices_forceplate_trial), 'o', 'linewidth', 1, 'markersize', 10);
            plot(time_force_plate(stretch_start_indices_forceplate_trial), fzl_trajectory(stretch_start_indices_forceplate_trial), '<', 'linewidth', 1, 'markersize', 10);
            plot(time_force_plate(stretch_end_indices_forceplate_trial), fzl_trajectory(stretch_end_indices_forceplate_trial), '>', 'linewidth', 1, 'markersize', 10);
            legend('fzl', 'left touchdown','left pushoff','triggered','strech start?','stretch end')
            
            % right events
            figure; axes_right = axes; hold on; title('Right'); xlim([0, recordingTime]);
            plot(time_force_plate, fzr_trajectory)
            plot(time_mocap(right_touchdown_indices_mocap), zeros(size(right_touchdown_indices_mocap)), 'v', 'linewidth', 1, 'markersize', 10);
            plot(time_mocap(right_pushoff_indices_mocap), zeros(size(right_pushoff_indices_mocap)), '^', 'linewidth', 1, 'markersize', 10);
            plot(time_force_plate(trigger_indices_forceplate_trial), fzr_trajectory(trigger_indices_forceplate_trial), 'o', 'linewidth', 1, 'markersize', 10);
            plot(time_force_plate(stretch_start_indices_forceplate_trial), fzr_trajectory(stretch_start_indices_forceplate_trial), '<', 'linewidth', 1, 'markersize', 10);
            plot(time_force_plate(stretch_end_indices_forceplate_trial), fzr_trajectory(stretch_end_indices_forceplate_trial), '>', 'linewidth', 1, 'markersize', 10);
            legend('fzr','right touchdown','right pushoff', 'triggered', ' stretch start','stretch end')
            
            linkaxes([axes_left, axes_right check_axes], 'x')
            distFig('rows', 3)
        end

        disp(['Trial ' num2str(i_trial) ' completed']);
        
    end

    median_number_of_indices_per_step = median(stretch_length_indices_forceplate);
    step_time_mean = mean(step_times);
    time_normalized = linspace(0, step_time_mean, number_of_time_steps_normalized);
end 

if view_totals_and_removals
      % Subject Totals and Removals
            total_steps_bytrial{i_trial+1} = removal_flags;
            total_steps = length(find(strcmp(total_polarity_condition_list,'POSITIVE'))) + ...
                length(find(strcmp(total_polarity_condition_list,'NEGATIVE'))) + ...
                length(find(strcmp(total_polarity_condition_list,'CONTROL')));
            total_steps_analyzed = length(find(strcmp(condition_polarity_list,'POSITIVE'))) + ...
                length(find(strcmp(condition_polarity_list,'NEGATIVE'))) + ...
                length(find(strcmp(condition_polarity_list,'CONTROL')));
            total_positive_steps = length(find(strcmp(total_polarity_condition_list,'POSITIVE')));
            total_positive_steps_analyzed = length(find(strcmp(condition_polarity_list,'POSITIVE'))); 
            total_negative_steps = length(find(strcmp(total_polarity_condition_list,'NEGATIVE')));
            total_negative_steps_analyzed = length(find(strcmp(condition_polarity_list,'NEGATIVE')));
            total_control_steps = length(find(strcmp(total_polarity_condition_list,'CONTROL')));
            total_control_steps_analyzed = length(find(strcmp(condition_polarity_list,'CONTROL')));
            total_left_steps =  length(find(strcmp(total_step_condition_list, 'LEFT'))); 
            total_left_steps_analyzed = length(find(strcmp(condition_stance_foot_list, 'LEFT')));
            total_right_steps = length(find(strcmp(total_step_condition_list, 'RIGHT')));
            total_right_steps_analyzed = length(find(strcmp(condition_stance_foot_list, 'RIGHT')));
            
            bar_data = [total_steps total_steps_analyzed ; total_left_steps total_left_steps_analyzed; ...
                total_right_steps total_right_steps_analyzed ; ...
                total_control_steps total_control_steps_analyzed; ...
                total_positive_steps total_positive_steps_analyzed; ...
                total_negative_steps total_negative_steps_analyzed];
            labels = {'','Total Events','Left Steps Triggered','Right Steps Triggered','Control Perturbations', 'Positive Perturbations', 'Negative Perturbations',''};
            figure; hold on
            x = [0:2:10];
            bar(x,bar_data)
            legend('Number of Events Collected',  'Number of Events Used', 'Fontsize', 12)
            set(gca,'xticklabel',labels, 'Fontsize',12)
            
            total_left_steps =  length(find(strcmp(total_step_condition_list, 'LEFT'))); 
            total_left_steps_analyzed = length(find(strcmp(condition_stance_foot_list, 'LEFT')));
            total_right_steps = length(find(strcmp(total_step_condition_list, 'RIGHT')));
            total_right_steps_analyzed = length(find(strcmp(condition_stance_foot_list, 'RIGHT')));
            
%             positive_analyzed_steps_concat  = 
%             negative_analyzed_steps_concat  = 
%             control_analyzed_steps_concat   = 
%             positive_analyzed_steps_removed = 
%             negative_analyzed_steps_removed =
%             control_analyzed_steps_removed =

end

%% calculate responses
if calculate_responses
    % extract conditions
    conditions_left_control = strcmp(condition_stance_foot_list, 'LEFT') & strcmp(condition_polarity_list, 'CONTROL');
    conditions_right_control = strcmp(condition_stance_foot_list, 'RIGHT') & strcmp(condition_polarity_list, 'CONTROL');
    conditions_left_positive_0ms = strcmp(condition_stance_foot_list, 'LEFT') & strcmp(condition_polarity_list, 'POSITIVE') & strcmp(condition_delay_list, '0ms');
    conditions_left_negative_0ms = strcmp(condition_stance_foot_list, 'LEFT') & strcmp(condition_polarity_list, 'NEGATIVE') & strcmp(condition_delay_list, '0ms');
    conditions_right_positive_0ms = strcmp(condition_stance_foot_list, 'RIGHT') & strcmp(condition_polarity_list, 'POSITIVE') & strcmp(condition_delay_list, '0ms');
    conditions_right_negative_0ms = strcmp(condition_stance_foot_list, 'RIGHT') & strcmp(condition_polarity_list, 'NEGATIVE') & strcmp(condition_delay_list, '0ms');
    conditions_left_positive_150ms = strcmp(condition_stance_foot_list, 'LEFT') & strcmp(condition_polarity_list, 'POSITIVE') & strcmp(condition_delay_list, '150ms');
    conditions_left_negative_150ms = strcmp(condition_stance_foot_list, 'LEFT') & strcmp(condition_polarity_list, 'NEGATIVE') & strcmp(condition_delay_list, '150ms');
    conditions_right_positive_150ms = strcmp(condition_stance_foot_list, 'RIGHT') & strcmp(condition_polarity_list, 'POSITIVE') & strcmp(condition_delay_list, '150ms');
    conditions_right_negative_150ms = strcmp(condition_stance_foot_list, 'RIGHT') & strcmp(condition_polarity_list, 'NEGATIVE') & strcmp(condition_delay_list, '150ms');
    conditions_left_positive_450ms = strcmp(condition_stance_foot_list, 'LEFT') & strcmp(condition_polarity_list, 'POSITIVE') & strcmp(condition_delay_list, '450ms');
    conditions_left_negative_450ms = strcmp(condition_stance_foot_list, 'LEFT') & strcmp(condition_polarity_list, 'NEGATIVE') & strcmp(condition_delay_list, '450ms');
    conditions_right_positive_450ms = strcmp(condition_stance_foot_list, 'RIGHT') & strcmp(condition_polarity_list, 'POSITIVE') & strcmp(condition_delay_list, '450ms');
    conditions_right_negative_450ms = strcmp(condition_stance_foot_list, 'RIGHT') & strcmp(condition_polarity_list, 'NEGATIVE') & strcmp(condition_delay_list, '450ms');
    
    conditions_left_0ms = conditions_left_positive_0ms | conditions_left_negative_0ms;
    conditions_right_0ms = conditions_right_positive_0ms | conditions_right_negative_0ms;
    conditions_left_150ms = conditions_left_positive_150ms | conditions_left_negative_150ms;
    conditions_right_150ms = conditions_right_positive_150ms | conditions_right_negative_150ms;
    conditions_left_450ms = conditions_left_positive_450ms | conditions_left_negative_450ms;
    conditions_right_450ms = conditions_right_positive_450ms | conditions_right_negative_450ms;
    
    % calculate control means
    left_cop_x_mean_left_control = mean(left_cop_x_normalized(:, conditions_left_control), 2);
    right_cop_x_mean_right_control = mean(right_cop_x_normalized(:, conditions_right_control), 2);
    left_heel_x_pos_mean_left_control = mean(left_heel_x_pos_normalized(:, conditions_left_control), 2);
    left_heel_x_pos_mean_right_control = mean(left_heel_x_pos_normalized(:, conditions_right_control), 2);
    right_heel_x_pos_mean_left_control = mean(right_heel_x_pos_normalized(:, conditions_left_control), 2);
    right_heel_x_pos_mean_right_control = mean(right_heel_x_pos_normalized(:, conditions_right_control), 2);
    left_heel_y_pos_mean_left_control = mean(left_heel_y_pos_normalized(:, conditions_left_control), 2);
    left_heel_y_pos_mean_right_control = mean(left_heel_y_pos_normalized(:, conditions_right_control), 2);
    right_heel_y_pos_mean_left_control = mean(right_heel_y_pos_normalized(:, conditions_left_control), 2);
    right_heel_y_pos_mean_right_control = mean(right_heel_y_pos_normalized(:, conditions_right_control), 2);
    
    if emg_present
        left_glutmed_emg_mean_left_control = mean(left_glutmed_emg_normalized(:, conditions_left_control), 2);
        left_tibiant_emg_mean_left_control = mean(left_tibiant_emg_normalized(:, conditions_left_control), 2);
        left_perolng_emg_mean_left_control = mean(left_perolng_emg_normalized(:, conditions_left_control), 2);
        right_glutmed_emg_mean_left_control = mean(right_glutmed_emg_normalized(:, conditions_left_control), 2);
        right_tibiant_emg_mean_left_control = mean(right_tibiant_emg_normalized(:, conditions_left_control), 2);
        right_perolng_emg_mean_left_control = mean(right_perolng_emg_normalized(:, conditions_left_control), 2);
        left_glutmed_emg_mean_right_control = mean(left_glutmed_emg_normalized(:, conditions_right_control), 2);
        left_tibiant_emg_mean_right_control = mean(left_tibiant_emg_normalized(:, conditions_right_control), 2);
        left_perolng_emg_mean_right_control = mean(left_perolng_emg_normalized(:, conditions_right_control), 2);
        right_glutmed_emg_mean_right_control = mean(right_glutmed_emg_normalized(:, conditions_right_control), 2);
        right_tibiant_emg_mean_right_control = mean(right_tibiant_emg_normalized(:, conditions_right_control), 2);
        right_perolng_emg_mean_right_control = mean(right_perolng_emg_normalized(:, conditions_right_control), 2);
    end
    
    % CoP response
    number_of_stretches = length(step_times);
    left_cop_x_response = left_cop_x_normalized - repmat(left_cop_x_mean_left_control, 1, number_of_stretches);
    right_cop_x_response = right_cop_x_normalized - repmat(right_cop_x_mean_right_control, 1, number_of_stretches);
    
    % heel response
    left_heel_x_pos_response = left_heel_x_pos_normalized - repmat(left_heel_x_pos_mean_right_control, 1, number_of_stretches);
    right_heel_x_pos_response = right_heel_x_pos_normalized - repmat(right_heel_x_pos_mean_left_control, 1, number_of_stretches);
    left_heel_y_pos_response = left_heel_y_pos_normalized - repmat(left_heel_y_pos_mean_right_control, 1, number_of_stretches);
    right_heel_y_pos_response = right_heel_y_pos_normalized - repmat(right_heel_y_pos_mean_left_control, 1, number_of_stretches);
    
    % calculate midstance
    [~, midstance_index_left] = min(abs(left_heel_y_pos_mean_left_control - right_heel_y_pos_mean_left_control));
    [~, midstance_index_right] = min(abs(left_heel_y_pos_mean_right_control - right_heel_y_pos_mean_right_control));
    midstance_index = floor(mean([midstance_index_left, midstance_index_right]));
    
    % calculate step responses
    [ ...
      Jacobian_left_positive, ...
      correlation_c_left_positive, ...
      correlation_p_left_positive, ...
      step_response_left_positive, ...
      stim_response_left_positive ...
    ] ...
    = ...
    calculateStepResponse ...
      ( ...
        right_heel_x_pos_normalized(:, conditions_left_positive_0ms), ...
        right_heel_x_pos_normalized(:, conditions_left_control), ...
        pelvis_x_pos_normalized(:, conditions_left_positive_0ms), ...
        pelvis_x_pos_normalized(:, conditions_left_control), ...
        pelvis_x_vel_normalized(:, conditions_left_positive_0ms), ...
        pelvis_x_vel_normalized(:, conditions_left_control), ...
        midstance_index ...
      );    
    [ ...
      Jacobian_left_negative, ...
      correlation_c_left_negative, ...
      correlation_p_left_negative, ...
      step_response_left_negative, ...
      stim_response_left_negative ...
    ] ...
    = ...
    calculateStepResponse ...
      ( ...
        right_heel_x_pos_normalized(:, conditions_left_negative_0ms), ...
        right_heel_x_pos_normalized(:, conditions_left_control), ...
        pelvis_x_pos_normalized(:, conditions_left_negative_0ms), ...
        pelvis_x_pos_normalized(:, conditions_left_control), ...
        pelvis_x_vel_normalized(:, conditions_left_negative_0ms), ...
        pelvis_x_vel_normalized(:, conditions_left_control), ...
        midstance_index ...
      );    

    % calculate step responses
    [ ...
      Jacobian_right_positive, ...
      correlation_c_right_positive, ...
      correlation_p_right_positive, ...
      step_response_right_positive, ...
      stim_response_right_positive ...
    ] ...
    = ...
    calculateStepResponse ...
      ( ...
        left_heel_x_pos_normalized(:, conditions_right_positive_0ms), ...
        left_heel_x_pos_normalized(:, conditions_right_control), ...
        pelvis_x_pos_normalized(:, conditions_right_positive_0ms), ...
        pelvis_x_pos_normalized(:, conditions_right_control), ...
        pelvis_x_vel_normalized(:, conditions_right_positive_0ms), ...
        pelvis_x_vel_normalized(:, conditions_right_control), ...
        midstance_index ...
      );    
    [ ...
      Jacobian_right_negative, ...
      correlation_c_right_negative, ...
      correlation_p_right_negative, ...
      step_response_right_negative, ...
      stim_response_right_negative ...
    ] ...
    = ...
    calculateStepResponse ...
      ( ...
        left_heel_x_pos_normalized(:, conditions_right_negative_0ms), ...
        left_heel_x_pos_normalized(:, conditions_right_control), ...
        pelvis_x_pos_normalized(:, conditions_right_negative_0ms), ...
        pelvis_x_pos_normalized(:, conditions_right_control), ...
        pelvis_x_vel_normalized(:, conditions_right_negative_0ms), ...
        pelvis_x_vel_normalized(:, conditions_right_control), ...
        midstance_index ...
      );    
    
    
    

    
end

%% calculate stats
if calculate_stats
    % absolute data
    left_cop_x_civ_left_control = tinv(0.975, sum(conditions_left_control)-1) * std(left_cop_x_normalized(:, conditions_left_control), 1, 2)/sqrt(sum(conditions_left_control));
    left_cop_x_mean_left_positive_0ms = mean(left_cop_x_normalized(:, conditions_left_positive_0ms), 2);
    left_cop_x_civ_left_positive_0ms = tinv(0.975, sum(conditions_left_positive_0ms)-1) * std(left_cop_x_normalized(:, conditions_left_positive_0ms), 1, 2)/sqrt(sum(conditions_left_positive_0ms));
    left_cop_x_mean_left_negative_0ms = mean(left_cop_x_normalized(:, conditions_left_negative_0ms), 2);
    left_cop_x_civ_left_negative_0ms = tinv(0.975, sum(conditions_left_negative_0ms)-1) * std(left_cop_x_normalized(:, conditions_left_negative_0ms), 1, 2)/sqrt(sum(conditions_left_negative_0ms));
    left_cop_x_mean_left_positive_150ms = mean(left_cop_x_normalized(:, conditions_left_positive_150ms), 2);
    left_cop_x_civ_left_positive_150ms = tinv(0.975, sum(conditions_left_positive_150ms)-1) * std(left_cop_x_normalized(:, conditions_left_positive_150ms), 1, 2)/sqrt(sum(conditions_left_positive_150ms));
    left_cop_x_mean_left_negative_150ms = mean(left_cop_x_normalized(:, conditions_left_negative_150ms), 2);
    left_cop_x_civ_left_negative_150ms = tinv(0.975, sum(conditions_left_negative_150ms)-1) * std(left_cop_x_normalized(:, conditions_left_negative_150ms), 1, 2)/sqrt(sum(conditions_left_negative_150ms));
    left_cop_x_mean_left_positive_450ms = mean(left_cop_x_normalized(:, conditions_left_positive_450ms), 2);
    left_cop_x_civ_left_positive_450ms = tinv(0.975, sum(conditions_left_positive_450ms)-1) * std(left_cop_x_normalized(:, conditions_left_positive_450ms), 1, 2)/sqrt(sum(conditions_left_positive_450ms));
    left_cop_x_mean_left_negative_450ms = mean(left_cop_x_normalized(:, conditions_left_negative_450ms), 2);
    left_cop_x_civ_left_negative_450ms = tinv(0.975, sum(conditions_left_negative_450ms)-1) * std(left_cop_x_normalized(:, conditions_left_negative_450ms), 1, 2)/sqrt(sum(conditions_left_negative_450ms));

    right_cop_x_civ_right_control = tinv(0.975, sum(conditions_right_control)-1) * std(right_cop_x_normalized(:, conditions_right_control), 1, 2)/sqrt(sum(conditions_right_control));
    right_cop_x_mean_right_positive_0ms = mean(right_cop_x_normalized(:, conditions_right_positive_0ms), 2);
    right_cop_x_civ_right_positive_0ms = tinv(0.975, sum(conditions_right_positive_0ms)-1) * std(right_cop_x_normalized(:, conditions_right_positive_0ms), 1, 2)/sqrt(sum(conditions_right_positive_0ms));
    right_cop_x_mean_right_negative_0ms = mean(right_cop_x_normalized(:, conditions_right_negative_0ms), 2);
    right_cop_x_civ_right_negative_0ms = tinv(0.975, sum(conditions_right_negative_0ms)-1) * std(right_cop_x_normalized(:, conditions_right_negative_0ms), 1, 2)/sqrt(sum(conditions_right_negative_0ms));
    right_cop_x_mean_right_positive_150ms = mean(right_cop_x_normalized(:, conditions_right_positive_150ms), 2);
    right_cop_x_civ_right_positive_150ms = tinv(0.975, sum(conditions_right_positive_150ms)-1) * std(right_cop_x_normalized(:, conditions_right_positive_150ms), 1, 2)/sqrt(sum(conditions_right_positive_150ms));
    right_cop_x_mean_right_negative_150ms = mean(right_cop_x_normalized(:, conditions_right_negative_150ms), 2);
    right_cop_x_civ_right_negative_150ms = tinv(0.975, sum(conditions_right_negative_150ms)-1) * std(right_cop_x_normalized(:, conditions_right_negative_150ms), 1, 2)/sqrt(sum(conditions_right_negative_150ms));
    right_cop_x_mean_right_positive_450ms = mean(right_cop_x_normalized(:, conditions_right_positive_450ms), 2);
    right_cop_x_civ_right_positive_450ms = tinv(0.975, sum(conditions_right_positive_450ms)-1) * std(right_cop_x_normalized(:, conditions_right_positive_450ms), 1, 2)/sqrt(sum(conditions_right_positive_450ms));
    right_cop_x_mean_right_negative_450ms = mean(right_cop_x_normalized(:, conditions_right_negative_450ms), 2);
    right_cop_x_civ_right_negative_450ms = tinv(0.975, sum(conditions_right_negative_450ms)-1) * std(right_cop_x_normalized(:, conditions_right_negative_450ms), 1, 2)/sqrt(sum(conditions_right_negative_450ms));

    left_heel_x_pos_civ_right_control = tinv(0.975, sum(conditions_right_control)-1) * std(left_heel_x_pos_normalized(:, conditions_right_control), 1, 2)/sqrt(sum(conditions_right_control));
    left_heel_x_pos_mean_right_positive_0ms = mean(left_heel_x_pos_normalized(:, conditions_right_positive_0ms), 2);
    left_heel_x_pos_civ_right_positive_0ms = tinv(0.975, sum(conditions_right_positive_0ms)-1) * std(left_heel_x_pos_normalized(:, conditions_right_positive_0ms), 1, 2)/sqrt(sum(conditions_right_positive_0ms));
    left_heel_x_pos_mean_right_negative_0ms = mean(left_heel_x_pos_normalized(:, conditions_right_negative_0ms), 2);
    left_heel_x_pos_civ_right_negative_0ms = tinv(0.975, sum(conditions_right_negative_0ms)-1) * std(left_heel_x_pos_normalized(:, conditions_right_negative_0ms), 1, 2)/sqrt(sum(conditions_right_negative_0ms));
    left_heel_x_pos_mean_right_positive_150ms = mean(left_heel_x_pos_normalized(:, conditions_right_positive_150ms), 2);
    left_heel_x_pos_civ_right_positive_150ms = tinv(0.975, sum(conditions_right_positive_150ms)-1) * std(left_heel_x_pos_normalized(:, conditions_right_positive_150ms), 1, 2)/sqrt(sum(conditions_right_positive_150ms));
    left_heel_x_pos_mean_right_negative_150ms = mean(left_heel_x_pos_normalized(:, conditions_right_negative_150ms), 2);
    left_heel_x_pos_civ_right_negative_150ms = tinv(0.975, sum(conditions_right_negative_150ms)-1) * std(left_heel_x_pos_normalized(:, conditions_right_negative_150ms), 1, 2)/sqrt(sum(conditions_right_negative_150ms));
    left_heel_x_pos_mean_right_positive_450ms = mean(left_heel_x_pos_normalized(:, conditions_right_positive_450ms), 2);
    left_heel_x_pos_civ_right_positive_450ms = tinv(0.975, sum(conditions_right_positive_450ms)-1) * std(left_heel_x_pos_normalized(:, conditions_right_positive_450ms), 1, 2)/sqrt(sum(conditions_right_positive_450ms));
    left_heel_x_pos_mean_right_negative_450ms = mean(left_heel_x_pos_normalized(:, conditions_right_negative_450ms), 2);
    left_heel_x_pos_civ_right_negative_450ms = tinv(0.975, sum(conditions_right_negative_450ms)-1) * std(left_heel_x_pos_normalized(:, conditions_right_negative_450ms), 1, 2)/sqrt(sum(conditions_right_negative_450ms));

    right_heel_x_pos_civ_left_control = tinv(0.975, sum(conditions_left_control)-1) * std(right_heel_x_pos_normalized(:, conditions_left_control), 1, 2)/sqrt(sum(conditions_left_control));
    right_heel_x_pos_mean_left_positive_0ms = mean(right_heel_x_pos_normalized(:, conditions_left_positive_0ms), 2);
    right_heel_x_pos_civ_left_positive_0ms = tinv(0.975, sum(conditions_left_positive_0ms)-1) * std(right_heel_x_pos_normalized(:, conditions_left_positive_0ms), 1, 2)/sqrt(sum(conditions_left_positive_0ms));
    right_heel_x_pos_mean_left_negative_0ms = mean(right_heel_x_pos_normalized(:, conditions_left_negative_0ms), 2);
    right_heel_x_pos_civ_left_negative_0ms = tinv(0.975, sum(conditions_left_negative_0ms)-1) * std(right_heel_x_pos_normalized(:, conditions_left_negative_0ms), 1, 2)/sqrt(sum(conditions_left_negative_0ms));
    right_heel_x_pos_mean_left_positive_150ms = mean(right_heel_x_pos_normalized(:, conditions_left_positive_150ms), 2);
    right_heel_x_pos_civ_left_positive_150ms = tinv(0.975, sum(conditions_left_positive_150ms)-1) * std(right_heel_x_pos_normalized(:, conditions_left_positive_150ms), 1, 2)/sqrt(sum(conditions_left_positive_150ms));
    right_heel_x_pos_mean_left_negative_150ms = mean(right_heel_x_pos_normalized(:, conditions_left_negative_150ms), 2);
    right_heel_x_pos_civ_left_negative_150ms = tinv(0.975, sum(conditions_left_negative_150ms)-1) * std(right_heel_x_pos_normalized(:, conditions_left_negative_150ms), 1, 2)/sqrt(sum(conditions_left_negative_150ms));
    right_heel_x_pos_mean_left_positive_450ms = mean(right_heel_x_pos_normalized(:, conditions_left_positive_450ms), 2);
    right_heel_x_pos_civ_left_positive_450ms = tinv(0.975, sum(conditions_left_positive_450ms)-1) * std(right_heel_x_pos_normalized(:, conditions_left_positive_450ms), 1, 2)/sqrt(sum(conditions_left_positive_450ms));
    right_heel_x_pos_mean_left_negative_450ms = mean(right_heel_x_pos_normalized(:, conditions_left_negative_450ms), 2);
    right_heel_x_pos_civ_left_negative_450ms = tinv(0.975, sum(conditions_left_negative_450ms)-1) * std(right_heel_x_pos_normalized(:, conditions_left_negative_450ms), 1, 2)/sqrt(sum(conditions_left_negative_450ms));
    
    if emg_present
        left_glutmed_emg_civ_left_control = tinv(0.975, sum(conditions_left_control)-1) * std(left_glutmed_emg_normalized(:, conditions_left_control), 1, 2)/sqrt(sum(conditions_left_control));
        left_glutmed_emg_mean_left_positive_0ms = mean(left_glutmed_emg_normalized(:, conditions_left_positive_0ms), 2);
        left_glutmed_emg_civ_left_positive_0ms = tinv(0.975, sum(conditions_left_positive_0ms)-1) * std(left_glutmed_emg_normalized(:, conditions_left_positive_0ms), 1, 2)/sqrt(sum(conditions_left_positive_0ms));
        left_glutmed_emg_mean_left_negative_0ms = mean(left_glutmed_emg_normalized(:, conditions_left_negative_0ms), 2);
        left_glutmed_emg_civ_left_negative_0ms = tinv(0.975, sum(conditions_left_negative_0ms)-1) * std(left_glutmed_emg_normalized(:, conditions_left_negative_0ms), 1, 2)/sqrt(sum(conditions_left_negative_0ms));
        left_glutmed_emg_mean_left_positive_150ms = mean(left_glutmed_emg_normalized(:, conditions_left_positive_150ms), 2);
        left_glutmed_emg_civ_left_positive_150ms = tinv(0.975, sum(conditions_left_positive_150ms)-1) * std(left_glutmed_emg_normalized(:, conditions_left_positive_150ms), 1, 2)/sqrt(sum(conditions_left_positive_150ms));
        left_glutmed_emg_mean_left_negative_150ms = mean(left_glutmed_emg_normalized(:, conditions_left_negative_150ms), 2);
        left_glutmed_emg_civ_left_negative_150ms = tinv(0.975, sum(conditions_left_negative_150ms)-1) * std(left_glutmed_emg_normalized(:, conditions_left_negative_150ms), 1, 2)/sqrt(sum(conditions_left_negative_150ms));
        left_glutmed_emg_mean_left_positive_450ms = mean(left_glutmed_emg_normalized(:, conditions_left_positive_450ms), 2);
        left_glutmed_emg_civ_left_positive_450ms = tinv(0.975, sum(conditions_left_positive_450ms)-1) * std(left_glutmed_emg_normalized(:, conditions_left_positive_450ms), 1, 2)/sqrt(sum(conditions_left_positive_450ms));
        left_glutmed_emg_mean_left_negative_450ms = mean(left_glutmed_emg_normalized(:, conditions_left_negative_450ms), 2);
        left_glutmed_emg_civ_left_negative_450ms = tinv(0.975, sum(conditions_left_negative_450ms)-1) * std(left_glutmed_emg_normalized(:, conditions_left_negative_450ms), 1, 2)/sqrt(sum(conditions_left_negative_450ms));

        left_glutmed_emg_civ_right_control = tinv(0.975, sum(conditions_right_control)-1) * std(left_glutmed_emg_normalized(:, conditions_right_control), 1, 2)/sqrt(sum(conditions_right_control));
        left_glutmed_emg_mean_right_positive_0ms = mean(left_glutmed_emg_normalized(:, conditions_right_positive_0ms), 2);
        left_glutmed_emg_civ_right_positive_0ms = tinv(0.975, sum(conditions_right_positive_0ms)-1) * std(left_glutmed_emg_normalized(:, conditions_right_positive_0ms), 1, 2)/sqrt(sum(conditions_right_positive_0ms));
        left_glutmed_emg_mean_right_negative_0ms = mean(left_glutmed_emg_normalized(:, conditions_right_negative_0ms), 2);
        left_glutmed_emg_civ_right_negative_0ms = tinv(0.975, sum(conditions_right_negative_0ms)-1) * std(left_glutmed_emg_normalized(:, conditions_right_negative_0ms), 1, 2)/sqrt(sum(conditions_right_negative_0ms));
        left_glutmed_emg_mean_right_positive_150ms = mean(left_glutmed_emg_normalized(:, conditions_right_positive_150ms), 2);
        left_glutmed_emg_civ_right_positive_150ms = tinv(0.975, sum(conditions_right_positive_150ms)-1) * std(left_glutmed_emg_normalized(:, conditions_right_positive_150ms), 1, 2)/sqrt(sum(conditions_right_positive_150ms));
        left_glutmed_emg_mean_right_negative_150ms = mean(left_glutmed_emg_normalized(:, conditions_right_negative_150ms), 2);
        left_glutmed_emg_civ_right_negative_150ms = tinv(0.975, sum(conditions_right_negative_150ms)-1) * std(left_glutmed_emg_normalized(:, conditions_right_negative_150ms), 1, 2)/sqrt(sum(conditions_right_negative_150ms));
        left_glutmed_emg_mean_right_positive_450ms = mean(left_glutmed_emg_normalized(:, conditions_right_positive_450ms), 2);
        left_glutmed_emg_civ_right_positive_450ms = tinv(0.975, sum(conditions_right_positive_450ms)-1) * std(left_glutmed_emg_normalized(:, conditions_right_positive_450ms), 1, 2)/sqrt(sum(conditions_right_positive_450ms));
        left_glutmed_emg_mean_right_negative_450ms = mean(left_glutmed_emg_normalized(:, conditions_right_negative_450ms), 2);
        left_glutmed_emg_civ_right_negative_450ms = tinv(0.975, sum(conditions_right_negative_450ms)-1) * std(left_glutmed_emg_normalized(:, conditions_right_negative_450ms), 1, 2)/sqrt(sum(conditions_right_negative_450ms));

        left_tibiant_emg_civ_left_control = tinv(0.975, sum(conditions_left_control)-1) * std(left_tibiant_emg_normalized(:, conditions_left_control), 1, 2)/sqrt(sum(conditions_left_control));
        left_tibiant_emg_mean_left_positive_0ms = mean(left_tibiant_emg_normalized(:, conditions_left_positive_0ms), 2);
        left_tibiant_emg_civ_left_positive_0ms = tinv(0.975, sum(conditions_left_positive_0ms)-1) * std(left_tibiant_emg_normalized(:, conditions_left_positive_0ms), 1, 2)/sqrt(sum(conditions_left_positive_0ms));
        left_tibiant_emg_mean_left_negative_0ms = mean(left_tibiant_emg_normalized(:, conditions_left_negative_0ms), 2);
        left_tibiant_emg_civ_left_negative_0ms = tinv(0.975, sum(conditions_left_negative_0ms)-1) * std(left_tibiant_emg_normalized(:, conditions_left_negative_0ms), 1, 2)/sqrt(sum(conditions_left_negative_0ms));
        left_tibiant_emg_mean_left_positive_150ms = mean(left_tibiant_emg_normalized(:, conditions_left_positive_150ms), 2);
        left_tibiant_emg_civ_left_positive_150ms = tinv(0.975, sum(conditions_left_positive_150ms)-1) * std(left_tibiant_emg_normalized(:, conditions_left_positive_150ms), 1, 2)/sqrt(sum(conditions_left_positive_150ms));
        left_tibiant_emg_mean_left_negative_150ms = mean(left_tibiant_emg_normalized(:, conditions_left_negative_150ms), 2);
        left_tibiant_emg_civ_left_negative_150ms = tinv(0.975, sum(conditions_left_negative_150ms)-1) * std(left_tibiant_emg_normalized(:, conditions_left_negative_150ms), 1, 2)/sqrt(sum(conditions_left_negative_150ms));
        left_tibiant_emg_mean_left_positive_450ms = mean(left_tibiant_emg_normalized(:, conditions_left_positive_450ms), 2);
        left_tibiant_emg_civ_left_positive_450ms = tinv(0.975, sum(conditions_left_positive_450ms)-1) * std(left_tibiant_emg_normalized(:, conditions_left_positive_450ms), 1, 2)/sqrt(sum(conditions_left_positive_450ms));
        left_tibiant_emg_mean_left_negative_450ms = mean(left_tibiant_emg_normalized(:, conditions_left_negative_450ms), 2);
        left_tibiant_emg_civ_left_negative_450ms = tinv(0.975, sum(conditions_left_negative_450ms)-1) * std(left_tibiant_emg_normalized(:, conditions_left_negative_450ms), 1, 2)/sqrt(sum(conditions_left_negative_450ms));

        left_tibiant_emg_civ_right_control = tinv(0.975, sum(conditions_right_control)-1) * std(left_tibiant_emg_normalized(:, conditions_right_control), 1, 2)/sqrt(sum(conditions_right_control));
        left_tibiant_emg_mean_right_positive_0ms = mean(left_tibiant_emg_normalized(:, conditions_right_positive_0ms), 2);
        left_tibiant_emg_civ_right_positive_0ms = tinv(0.975, sum(conditions_right_positive_0ms)-1) * std(left_tibiant_emg_normalized(:, conditions_right_positive_0ms), 1, 2)/sqrt(sum(conditions_right_positive_0ms));
        left_tibiant_emg_mean_right_negative_0ms = mean(left_tibiant_emg_normalized(:, conditions_right_negative_0ms), 2);
        left_tibiant_emg_civ_right_negative_0ms = tinv(0.975, sum(conditions_right_negative_0ms)-1) * std(left_tibiant_emg_normalized(:, conditions_right_negative_0ms), 1, 2)/sqrt(sum(conditions_right_negative_0ms));
        left_tibiant_emg_mean_right_positive_150ms = mean(left_tibiant_emg_normalized(:, conditions_right_positive_150ms), 2);
        left_tibiant_emg_civ_right_positive_150ms = tinv(0.975, sum(conditions_right_positive_150ms)-1) * std(left_tibiant_emg_normalized(:, conditions_right_positive_150ms), 1, 2)/sqrt(sum(conditions_right_positive_150ms));
        left_tibiant_emg_mean_right_negative_150ms = mean(left_tibiant_emg_normalized(:, conditions_right_negative_150ms), 2);
        left_tibiant_emg_civ_right_negative_150ms = tinv(0.975, sum(conditions_right_negative_150ms)-1) * std(left_tibiant_emg_normalized(:, conditions_right_negative_150ms), 1, 2)/sqrt(sum(conditions_right_negative_150ms));
        left_tibiant_emg_mean_right_positive_450ms = mean(left_tibiant_emg_normalized(:, conditions_right_positive_450ms), 2);
        left_tibiant_emg_civ_right_positive_450ms = tinv(0.975, sum(conditions_right_positive_450ms)-1) * std(left_tibiant_emg_normalized(:, conditions_right_positive_450ms), 1, 2)/sqrt(sum(conditions_right_positive_450ms));
        left_tibiant_emg_mean_right_negative_450ms = mean(left_tibiant_emg_normalized(:, conditions_right_negative_450ms), 2);
        left_tibiant_emg_civ_right_negative_450ms = tinv(0.975, sum(conditions_right_negative_450ms)-1) * std(left_tibiant_emg_normalized(:, conditions_right_negative_450ms), 1, 2)/sqrt(sum(conditions_right_negative_450ms));

        left_perolng_emg_civ_left_control = tinv(0.975, sum(conditions_left_control)-1) * std(left_perolng_emg_normalized(:, conditions_left_control), 1, 2)/sqrt(sum(conditions_left_control));
        left_perolng_emg_mean_left_positive_0ms = mean(left_perolng_emg_normalized(:, conditions_left_positive_0ms), 2);
        left_perolng_emg_civ_left_positive_0ms = tinv(0.975, sum(conditions_left_positive_0ms)-1) * std(left_perolng_emg_normalized(:, conditions_left_positive_0ms), 1, 2)/sqrt(sum(conditions_left_positive_0ms));
        left_perolng_emg_mean_left_negative_0ms = mean(left_perolng_emg_normalized(:, conditions_left_negative_0ms), 2);
        left_perolng_emg_civ_left_negative_0ms = tinv(0.975, sum(conditions_left_negative_0ms)-1) * std(left_perolng_emg_normalized(:, conditions_left_negative_0ms), 1, 2)/sqrt(sum(conditions_left_negative_0ms));
        left_perolng_emg_mean_left_positive_150ms = mean(left_perolng_emg_normalized(:, conditions_left_positive_150ms), 2);
        left_perolng_emg_civ_left_positive_150ms = tinv(0.975, sum(conditions_left_positive_150ms)-1) * std(left_perolng_emg_normalized(:, conditions_left_positive_150ms), 1, 2)/sqrt(sum(conditions_left_positive_150ms));
        left_perolng_emg_mean_left_negative_150ms = mean(left_perolng_emg_normalized(:, conditions_left_negative_150ms), 2);
        left_perolng_emg_civ_left_negative_150ms = tinv(0.975, sum(conditions_left_negative_150ms)-1) * std(left_perolng_emg_normalized(:, conditions_left_negative_150ms), 1, 2)/sqrt(sum(conditions_left_negative_150ms));
        left_perolng_emg_mean_left_positive_450ms = mean(left_perolng_emg_normalized(:, conditions_left_positive_450ms), 2);
        left_perolng_emg_civ_left_positive_450ms = tinv(0.975, sum(conditions_left_positive_450ms)-1) * std(left_perolng_emg_normalized(:, conditions_left_positive_450ms), 1, 2)/sqrt(sum(conditions_left_positive_450ms));
        left_perolng_emg_mean_left_negative_450ms = mean(left_perolng_emg_normalized(:, conditions_left_negative_450ms), 2);
        left_perolng_emg_civ_left_negative_450ms = tinv(0.975, sum(conditions_left_negative_450ms)-1) * std(left_perolng_emg_normalized(:, conditions_left_negative_450ms), 1, 2)/sqrt(sum(conditions_left_negative_450ms));

        left_perolng_emg_civ_right_control = tinv(0.975, sum(conditions_right_control)-1) * std(left_perolng_emg_normalized(:, conditions_right_control), 1, 2)/sqrt(sum(conditions_right_control));
        left_perolng_emg_mean_right_positive_0ms = mean(left_perolng_emg_normalized(:, conditions_right_positive_0ms), 2);
        left_perolng_emg_civ_right_positive_0ms = tinv(0.975, sum(conditions_right_positive_0ms)-1) * std(left_perolng_emg_normalized(:, conditions_right_positive_0ms), 1, 2)/sqrt(sum(conditions_right_positive_0ms));
        left_perolng_emg_mean_right_negative_0ms = mean(left_perolng_emg_normalized(:, conditions_right_negative_0ms), 2);
        left_perolng_emg_civ_right_negative_0ms = tinv(0.975, sum(conditions_right_negative_0ms)-1) * std(left_perolng_emg_normalized(:, conditions_right_negative_0ms), 1, 2)/sqrt(sum(conditions_right_negative_0ms));
        left_perolng_emg_mean_right_positive_150ms = mean(left_perolng_emg_normalized(:, conditions_right_positive_150ms), 2);
        left_perolng_emg_civ_right_positive_150ms = tinv(0.975, sum(conditions_right_positive_150ms)-1) * std(left_perolng_emg_normalized(:, conditions_right_positive_150ms), 1, 2)/sqrt(sum(conditions_right_positive_150ms));
        left_perolng_emg_mean_right_negative_150ms = mean(left_perolng_emg_normalized(:, conditions_right_negative_150ms), 2);
        left_perolng_emg_civ_right_negative_150ms = tinv(0.975, sum(conditions_right_negative_150ms)-1) * std(left_perolng_emg_normalized(:, conditions_right_negative_150ms), 1, 2)/sqrt(sum(conditions_right_negative_150ms));
        left_perolng_emg_mean_right_positive_450ms = mean(left_perolng_emg_normalized(:, conditions_right_positive_450ms), 2);
        left_perolng_emg_civ_right_positive_450ms = tinv(0.975, sum(conditions_right_positive_450ms)-1) * std(left_perolng_emg_normalized(:, conditions_right_positive_450ms), 1, 2)/sqrt(sum(conditions_right_positive_450ms));
        left_perolng_emg_mean_right_negative_450ms = mean(left_perolng_emg_normalized(:, conditions_right_negative_450ms), 2);
        left_perolng_emg_civ_right_negative_450ms = tinv(0.975, sum(conditions_right_negative_450ms)-1) * std(left_perolng_emg_normalized(:, conditions_right_negative_450ms), 1, 2)/sqrt(sum(conditions_right_negative_450ms));

        right_glutmed_emg_civ_left_control = tinv(0.975, sum(conditions_left_control)-1) * std(right_glutmed_emg_normalized(:, conditions_left_control), 1, 2)/sqrt(sum(conditions_left_control));
        right_glutmed_emg_mean_left_positive_0ms = mean(right_glutmed_emg_normalized(:, conditions_left_positive_0ms), 2);
        right_glutmed_emg_civ_left_positive_0ms = tinv(0.975, sum(conditions_left_positive_0ms)-1) * std(right_glutmed_emg_normalized(:, conditions_left_positive_0ms), 1, 2)/sqrt(sum(conditions_left_positive_0ms));
        right_glutmed_emg_mean_left_negative_0ms = mean(right_glutmed_emg_normalized(:, conditions_left_negative_0ms), 2);
        right_glutmed_emg_civ_left_negative_0ms = tinv(0.975, sum(conditions_left_negative_0ms)-1) * std(right_glutmed_emg_normalized(:, conditions_left_negative_0ms), 1, 2)/sqrt(sum(conditions_left_negative_0ms));
        right_glutmed_emg_mean_left_positive_150ms = mean(right_glutmed_emg_normalized(:, conditions_left_positive_150ms), 2);
        right_glutmed_emg_civ_left_positive_150ms = tinv(0.975, sum(conditions_left_positive_150ms)-1) * std(right_glutmed_emg_normalized(:, conditions_left_positive_150ms), 1, 2)/sqrt(sum(conditions_left_positive_150ms));
        right_glutmed_emg_mean_left_negative_150ms = mean(right_glutmed_emg_normalized(:, conditions_left_negative_150ms), 2);
        right_glutmed_emg_civ_left_negative_150ms = tinv(0.975, sum(conditions_left_negative_150ms)-1) * std(right_glutmed_emg_normalized(:, conditions_left_negative_150ms), 1, 2)/sqrt(sum(conditions_left_negative_150ms));
        right_glutmed_emg_mean_left_positive_450ms = mean(right_glutmed_emg_normalized(:, conditions_left_positive_450ms), 2);
        right_glutmed_emg_civ_left_positive_450ms = tinv(0.975, sum(conditions_left_positive_450ms)-1) * std(right_glutmed_emg_normalized(:, conditions_left_positive_450ms), 1, 2)/sqrt(sum(conditions_left_positive_450ms));
        right_glutmed_emg_mean_left_negative_450ms = mean(right_glutmed_emg_normalized(:, conditions_left_negative_450ms), 2);
        right_glutmed_emg_civ_left_negative_450ms = tinv(0.975, sum(conditions_left_negative_450ms)-1) * std(right_glutmed_emg_normalized(:, conditions_left_negative_450ms), 1, 2)/sqrt(sum(conditions_left_negative_450ms));

        right_glutmed_emg_civ_right_control = tinv(0.975, sum(conditions_right_control)-1) * std(right_glutmed_emg_normalized(:, conditions_right_control), 1, 2)/sqrt(sum(conditions_right_control));
        right_glutmed_emg_mean_right_positive_0ms = mean(right_glutmed_emg_normalized(:, conditions_right_positive_0ms), 2);
        right_glutmed_emg_civ_right_positive_0ms = tinv(0.975, sum(conditions_right_positive_0ms)-1) * std(right_glutmed_emg_normalized(:, conditions_right_positive_0ms), 1, 2)/sqrt(sum(conditions_right_positive_0ms));
        right_glutmed_emg_mean_right_negative_0ms = mean(right_glutmed_emg_normalized(:, conditions_right_negative_0ms), 2);
        right_glutmed_emg_civ_right_negative_0ms = tinv(0.975, sum(conditions_right_negative_0ms)-1) * std(right_glutmed_emg_normalized(:, conditions_right_negative_0ms), 1, 2)/sqrt(sum(conditions_right_negative_0ms));
        right_glutmed_emg_mean_right_positive_150ms = mean(right_glutmed_emg_normalized(:, conditions_right_positive_150ms), 2);
        right_glutmed_emg_civ_right_positive_150ms = tinv(0.975, sum(conditions_right_positive_150ms)-1) * std(right_glutmed_emg_normalized(:, conditions_right_positive_150ms), 1, 2)/sqrt(sum(conditions_right_positive_150ms));
        right_glutmed_emg_mean_right_negative_150ms = mean(right_glutmed_emg_normalized(:, conditions_right_negative_150ms), 2);
        right_glutmed_emg_civ_right_negative_150ms = tinv(0.975, sum(conditions_right_negative_150ms)-1) * std(right_glutmed_emg_normalized(:, conditions_right_negative_150ms), 1, 2)/sqrt(sum(conditions_right_negative_150ms));
        right_glutmed_emg_mean_right_positive_450ms = mean(right_glutmed_emg_normalized(:, conditions_right_positive_450ms), 2);
        right_glutmed_emg_civ_right_positive_450ms = tinv(0.975, sum(conditions_right_positive_450ms)-1) * std(right_glutmed_emg_normalized(:, conditions_right_positive_450ms), 1, 2)/sqrt(sum(conditions_right_positive_450ms));
        right_glutmed_emg_mean_right_negative_450ms = mean(right_glutmed_emg_normalized(:, conditions_right_negative_450ms), 2);
        right_glutmed_emg_civ_right_negative_450ms = tinv(0.975, sum(conditions_right_negative_450ms)-1) * std(right_glutmed_emg_normalized(:, conditions_right_negative_450ms), 1, 2)/sqrt(sum(conditions_right_negative_450ms));

        right_tibiant_emg_civ_left_control = tinv(0.975, sum(conditions_left_control)-1) * std(right_tibiant_emg_normalized(:, conditions_left_control), 1, 2)/sqrt(sum(conditions_left_control));
        right_tibiant_emg_mean_left_positive_0ms = mean(right_tibiant_emg_normalized(:, conditions_left_positive_0ms), 2);
        right_tibiant_emg_civ_left_positive_0ms = tinv(0.975, sum(conditions_left_positive_0ms)-1) * std(right_tibiant_emg_normalized(:, conditions_left_positive_0ms), 1, 2)/sqrt(sum(conditions_left_positive_0ms));
        right_tibiant_emg_mean_left_negative_0ms = mean(right_tibiant_emg_normalized(:, conditions_left_negative_0ms), 2);
        right_tibiant_emg_civ_left_negative_0ms = tinv(0.975, sum(conditions_left_negative_0ms)-1) * std(right_tibiant_emg_normalized(:, conditions_left_negative_0ms), 1, 2)/sqrt(sum(conditions_left_negative_0ms));
        right_tibiant_emg_mean_left_positive_150ms = mean(right_tibiant_emg_normalized(:, conditions_left_positive_150ms), 2);
        right_tibiant_emg_civ_left_positive_150ms = tinv(0.975, sum(conditions_left_positive_150ms)-1) * std(right_tibiant_emg_normalized(:, conditions_left_positive_150ms), 1, 2)/sqrt(sum(conditions_left_positive_150ms));
        right_tibiant_emg_mean_left_negative_150ms = mean(right_tibiant_emg_normalized(:, conditions_left_negative_150ms), 2);
        right_tibiant_emg_civ_left_negative_150ms = tinv(0.975, sum(conditions_left_negative_150ms)-1) * std(right_tibiant_emg_normalized(:, conditions_left_negative_150ms), 1, 2)/sqrt(sum(conditions_left_negative_150ms));
        right_tibiant_emg_mean_left_positive_450ms = mean(right_tibiant_emg_normalized(:, conditions_left_positive_450ms), 2);
        right_tibiant_emg_civ_left_positive_450ms = tinv(0.975, sum(conditions_left_positive_450ms)-1) * std(right_tibiant_emg_normalized(:, conditions_left_positive_450ms), 1, 2)/sqrt(sum(conditions_left_positive_450ms));
        right_tibiant_emg_mean_left_negative_450ms = mean(right_tibiant_emg_normalized(:, conditions_left_negative_450ms), 2);
        right_tibiant_emg_civ_left_negative_450ms = tinv(0.975, sum(conditions_left_negative_450ms)-1) * std(right_tibiant_emg_normalized(:, conditions_left_negative_450ms), 1, 2)/sqrt(sum(conditions_left_negative_450ms));

        right_tibiant_emg_civ_right_control = tinv(0.975, sum(conditions_right_control)-1) * std(right_tibiant_emg_normalized(:, conditions_right_control), 1, 2)/sqrt(sum(conditions_right_control));
        right_tibiant_emg_mean_right_positive_0ms = mean(right_tibiant_emg_normalized(:, conditions_right_positive_0ms), 2);
        right_tibiant_emg_civ_right_positive_0ms = tinv(0.975, sum(conditions_right_positive_0ms)-1) * std(right_tibiant_emg_normalized(:, conditions_right_positive_0ms), 1, 2)/sqrt(sum(conditions_right_positive_0ms));
        right_tibiant_emg_mean_right_negative_0ms = mean(right_tibiant_emg_normalized(:, conditions_right_negative_0ms), 2);
        right_tibiant_emg_civ_right_negative_0ms = tinv(0.975, sum(conditions_right_negative_0ms)-1) * std(right_tibiant_emg_normalized(:, conditions_right_negative_0ms), 1, 2)/sqrt(sum(conditions_right_negative_0ms));
        right_tibiant_emg_mean_right_positive_150ms = mean(right_tibiant_emg_normalized(:, conditions_right_positive_150ms), 2);
        right_tibiant_emg_civ_right_positive_150ms = tinv(0.975, sum(conditions_right_positive_150ms)-1) * std(right_tibiant_emg_normalized(:, conditions_right_positive_150ms), 1, 2)/sqrt(sum(conditions_right_positive_150ms));
        right_tibiant_emg_mean_right_negative_150ms = mean(right_tibiant_emg_normalized(:, conditions_right_negative_150ms), 2);
        right_tibiant_emg_civ_right_negative_150ms = tinv(0.975, sum(conditions_right_negative_150ms)-1) * std(right_tibiant_emg_normalized(:, conditions_right_negative_150ms), 1, 2)/sqrt(sum(conditions_right_negative_150ms));
        right_tibiant_emg_mean_right_positive_450ms = mean(right_tibiant_emg_normalized(:, conditions_right_positive_450ms), 2);
        right_tibiant_emg_civ_right_positive_450ms = tinv(0.975, sum(conditions_right_positive_450ms)-1) * std(right_tibiant_emg_normalized(:, conditions_right_positive_450ms), 1, 2)/sqrt(sum(conditions_right_positive_450ms));
        right_tibiant_emg_mean_right_negative_450ms = mean(right_tibiant_emg_normalized(:, conditions_right_negative_450ms), 2);
        right_tibiant_emg_civ_right_negative_450ms = tinv(0.975, sum(conditions_right_negative_450ms)-1) * std(right_tibiant_emg_normalized(:, conditions_right_negative_450ms), 1, 2)/sqrt(sum(conditions_right_negative_450ms));

        right_perolng_emg_civ_left_control = tinv(0.975, sum(conditions_left_control)-1) * std(right_perolng_emg_normalized(:, conditions_left_control), 1, 2)/sqrt(sum(conditions_left_control));
        right_perolng_emg_mean_left_positive_0ms = mean(right_perolng_emg_normalized(:, conditions_left_positive_0ms), 2);
        right_perolng_emg_civ_left_positive_0ms = tinv(0.975, sum(conditions_left_positive_0ms)-1) * std(right_perolng_emg_normalized(:, conditions_left_positive_0ms), 1, 2)/sqrt(sum(conditions_left_positive_0ms));
        right_perolng_emg_mean_left_negative_0ms = mean(right_perolng_emg_normalized(:, conditions_left_negative_0ms), 2);
        right_perolng_emg_civ_left_negative_0ms = tinv(0.975, sum(conditions_left_negative_0ms)-1) * std(right_perolng_emg_normalized(:, conditions_left_negative_0ms), 1, 2)/sqrt(sum(conditions_left_negative_0ms));
        right_perolng_emg_mean_left_positive_150ms = mean(right_perolng_emg_normalized(:, conditions_left_positive_150ms), 2);
        right_perolng_emg_civ_left_positive_150ms = tinv(0.975, sum(conditions_left_positive_150ms)-1) * std(right_perolng_emg_normalized(:, conditions_left_positive_150ms), 1, 2)/sqrt(sum(conditions_left_positive_150ms));
        right_perolng_emg_mean_left_negative_150ms = mean(right_perolng_emg_normalized(:, conditions_left_negative_150ms), 2);
        right_perolng_emg_civ_left_negative_150ms = tinv(0.975, sum(conditions_left_negative_150ms)-1) * std(right_perolng_emg_normalized(:, conditions_left_negative_150ms), 1, 2)/sqrt(sum(conditions_left_negative_150ms));
        right_perolng_emg_mean_left_positive_450ms = mean(right_perolng_emg_normalized(:, conditions_left_positive_450ms), 2);
        right_perolng_emg_civ_left_positive_450ms = tinv(0.975, sum(conditions_left_positive_450ms)-1) * std(right_perolng_emg_normalized(:, conditions_left_positive_450ms), 1, 2)/sqrt(sum(conditions_left_positive_450ms));
        right_perolng_emg_mean_left_negative_450ms = mean(right_perolng_emg_normalized(:, conditions_left_negative_450ms), 2);
        right_perolng_emg_civ_left_negative_450ms = tinv(0.975, sum(conditions_left_negative_450ms)-1) * std(right_perolng_emg_normalized(:, conditions_left_negative_450ms), 1, 2)/sqrt(sum(conditions_left_negative_450ms));

        right_perolng_emg_civ_right_control = tinv(0.975, sum(conditions_right_control)-1) * std(right_perolng_emg_normalized(:, conditions_right_control), 1, 2)/sqrt(sum(conditions_right_control));
        right_perolng_emg_mean_right_positive_0ms = mean(right_perolng_emg_normalized(:, conditions_right_positive_0ms), 2);
        right_perolng_emg_civ_right_positive_0ms = tinv(0.975, sum(conditions_right_positive_0ms)-1) * std(right_perolng_emg_normalized(:, conditions_right_positive_0ms), 1, 2)/sqrt(sum(conditions_right_positive_0ms));
        right_perolng_emg_mean_right_negative_0ms = mean(right_perolng_emg_normalized(:, conditions_right_negative_0ms), 2);
        right_perolng_emg_civ_right_negative_0ms = tinv(0.975, sum(conditions_right_negative_0ms)-1) * std(right_perolng_emg_normalized(:, conditions_right_negative_0ms), 1, 2)/sqrt(sum(conditions_right_negative_0ms));
        right_perolng_emg_mean_right_positive_150ms = mean(right_perolng_emg_normalized(:, conditions_right_positive_150ms), 2);
        right_perolng_emg_civ_right_positive_150ms = tinv(0.975, sum(conditions_right_positive_150ms)-1) * std(right_perolng_emg_normalized(:, conditions_right_positive_150ms), 1, 2)/sqrt(sum(conditions_right_positive_150ms));
        right_perolng_emg_mean_right_negative_150ms = mean(right_perolng_emg_normalized(:, conditions_right_negative_150ms), 2);
        right_perolng_emg_civ_right_negative_150ms = tinv(0.975, sum(conditions_right_negative_150ms)-1) * std(right_perolng_emg_normalized(:, conditions_right_negative_150ms), 1, 2)/sqrt(sum(conditions_right_negative_150ms));
        right_perolng_emg_mean_right_positive_450ms = mean(right_perolng_emg_normalized(:, conditions_right_positive_450ms), 2);
        right_perolng_emg_civ_right_positive_450ms = tinv(0.975, sum(conditions_right_positive_450ms)-1) * std(right_perolng_emg_normalized(:, conditions_right_positive_450ms), 1, 2)/sqrt(sum(conditions_right_positive_450ms));
        right_perolng_emg_mean_right_negative_450ms = mean(right_perolng_emg_normalized(:, conditions_right_negative_450ms), 2);
        right_perolng_emg_civ_right_negative_450ms = tinv(0.975, sum(conditions_right_negative_450ms)-1) * std(right_perolng_emg_normalized(:, conditions_right_negative_450ms), 1, 2)/sqrt(sum(conditions_right_negative_450ms));
    end
    
    % responses
    left_cop_x_response_mean_left_positive_0ms = mean(left_cop_x_response(:, conditions_left_positive_0ms), 2);
    left_cop_x_response_civ_left_positive_0ms = tinv(0.975, sum(conditions_left_positive_0ms)-1) * std(left_cop_x_response(:, conditions_left_positive_0ms), 1, 2)/sqrt(sum(conditions_left_positive_0ms));
    left_cop_x_response_mean_left_negative_0ms = mean(left_cop_x_response(:, conditions_left_negative_0ms), 2);
    left_cop_x_response_civ_left_negative_0ms = tinv(0.975, sum(conditions_left_negative_0ms)-1) * std(left_cop_x_response(:, conditions_left_negative_0ms), 1, 2)/sqrt(sum(conditions_left_negative_0ms));
    left_cop_x_response_mean_left_positive_150ms = mean(left_cop_x_response(:, conditions_left_positive_150ms), 2);
    left_cop_x_response_civ_left_positive_150ms = tinv(0.975, sum(conditions_left_positive_150ms)-1) * std(left_cop_x_response(:, conditions_left_positive_150ms), 1, 2)/sqrt(sum(conditions_left_positive_150ms));
    left_cop_x_response_mean_left_negative_150ms = mean(left_cop_x_response(:, conditions_left_negative_150ms), 2);
    left_cop_x_response_civ_left_negative_150ms = tinv(0.975, sum(conditions_left_negative_150ms)-1) * std(left_cop_x_response(:, conditions_left_negative_150ms), 1, 2)/sqrt(sum(conditions_left_negative_150ms));
    left_cop_x_response_mean_left_positive_450ms = mean(left_cop_x_response(:, conditions_left_positive_450ms), 2);
    left_cop_x_response_civ_left_positive_450ms = tinv(0.975, sum(conditions_left_positive_450ms)-1) * std(left_cop_x_response(:, conditions_left_positive_450ms), 1, 2)/sqrt(sum(conditions_left_positive_450ms));
    left_cop_x_response_mean_left_negative_450ms = mean(left_cop_x_response(:, conditions_left_negative_450ms), 2);
    left_cop_x_response_civ_left_negative_450ms = tinv(0.975, sum(conditions_left_negative_450ms)-1) * std(left_cop_x_response(:, conditions_left_negative_450ms), 1, 2)/sqrt(sum(conditions_left_negative_450ms));
    
    right_cop_x_response_mean_right_positive_0ms = mean(right_cop_x_response(:, conditions_right_positive_0ms), 2);
    right_cop_x_response_civ_right_positive_0ms = tinv(0.975, sum(conditions_right_positive_0ms)-1) * std(right_cop_x_response(:, conditions_right_positive_0ms), 1, 2)/sqrt(sum(conditions_right_positive_0ms));
    right_cop_x_response_mean_right_negative_0ms = mean(right_cop_x_response(:, conditions_right_negative_0ms), 2);
    right_cop_x_response_civ_right_negative_0ms = tinv(0.975, sum(conditions_right_negative_0ms)-1) * std(right_cop_x_response(:, conditions_right_negative_0ms), 1, 2)/sqrt(sum(conditions_right_negative_0ms));
    right_cop_x_response_mean_right_positive_150ms = mean(right_cop_x_response(:, conditions_right_positive_150ms), 2);
    right_cop_x_response_civ_right_positive_150ms = tinv(0.975, sum(conditions_right_positive_150ms)-1) * std(right_cop_x_response(:, conditions_right_positive_150ms), 1, 2)/sqrt(sum(conditions_right_positive_150ms));
    right_cop_x_response_mean_right_negative_150ms = mean(right_cop_x_response(:, conditions_right_negative_150ms), 2);
    right_cop_x_response_civ_right_negative_150ms = tinv(0.975, sum(conditions_right_negative_150ms)-1) * std(right_cop_x_response(:, conditions_right_negative_150ms), 1, 2)/sqrt(sum(conditions_right_negative_150ms));
    right_cop_x_response_mean_right_positive_450ms = mean(right_cop_x_response(:, conditions_right_positive_450ms), 2);
    right_cop_x_response_civ_right_positive_450ms = tinv(0.975, sum(conditions_right_positive_450ms)-1) * std(right_cop_x_response(:, conditions_right_positive_450ms), 1, 2)/sqrt(sum(conditions_right_positive_450ms));
    right_cop_x_response_mean_right_negative_450ms = mean(right_cop_x_response(:, conditions_right_negative_450ms), 2);
    right_cop_x_response_civ_right_negative_450ms = tinv(0.975, sum(conditions_right_negative_450ms)-1) * std(right_cop_x_response(:, conditions_right_negative_450ms), 1, 2)/sqrt(sum(conditions_right_negative_450ms));
    
    left_heel_x_pos_response_mean_right_positive_0ms = mean(left_heel_x_pos_response(:, conditions_right_positive_0ms), 2);
    left_heel_x_pos_response_civ_right_positive_0ms = tinv(0.975, sum(conditions_right_positive_0ms)-1) * std(left_heel_x_pos_response(:, conditions_right_positive_0ms), 1, 2)/sqrt(sum(conditions_right_positive_0ms));
    left_heel_x_pos_response_mean_right_negative_0ms = mean(left_heel_x_pos_response(:, conditions_right_negative_0ms), 2);
    left_heel_x_pos_response_civ_right_negative_0ms = tinv(0.975, sum(conditions_right_negative_0ms)-1) * std(left_heel_x_pos_response(:, conditions_right_negative_0ms), 1, 2)/sqrt(sum(conditions_right_negative_0ms));
    left_heel_x_pos_response_mean_right_positive_150ms = mean(left_heel_x_pos_response(:, conditions_right_positive_150ms), 2);
    left_heel_x_pos_response_civ_right_positive_150ms = tinv(0.975, sum(conditions_right_positive_150ms)-1) * std(left_heel_x_pos_response(:, conditions_right_positive_150ms), 1, 2)/sqrt(sum(conditions_right_positive_150ms));
    left_heel_x_pos_response_mean_right_negative_150ms = mean(left_heel_x_pos_response(:, conditions_right_negative_150ms), 2);
    left_heel_x_pos_response_civ_right_negative_150ms = tinv(0.975, sum(conditions_right_negative_150ms)-1) * std(left_heel_x_pos_response(:, conditions_right_negative_150ms), 1, 2)/sqrt(sum(conditions_right_negative_150ms));
    left_heel_x_pos_response_mean_right_positive_450ms = mean(left_heel_x_pos_response(:, conditions_right_positive_450ms), 2);
    left_heel_x_pos_response_civ_right_positive_450ms = tinv(0.975, sum(conditions_right_positive_450ms)-1) * std(left_heel_x_pos_response(:, conditions_right_positive_450ms), 1, 2)/sqrt(sum(conditions_right_positive_450ms));
    left_heel_x_pos_response_mean_right_negative_450ms = mean(left_heel_x_pos_response(:, conditions_right_negative_450ms), 2);
    left_heel_x_pos_response_civ_right_negative_450ms = tinv(0.975, sum(conditions_right_negative_450ms)-1) * std(left_heel_x_pos_response(:, conditions_right_negative_450ms), 1, 2)/sqrt(sum(conditions_right_negative_450ms));
    
    right_heel_x_pos_response_mean_left_positive_0ms = mean(right_heel_x_pos_response(:, conditions_left_positive_0ms), 2);
    right_heel_x_pos_response_civ_left_positive_0ms = tinv(0.975, sum(conditions_left_positive_0ms)-1) * std(right_heel_x_pos_response(:, conditions_left_positive_0ms), 1, 2)/sqrt(sum(conditions_left_positive_0ms));
    right_heel_x_pos_response_mean_left_negative_0ms = mean(right_heel_x_pos_response(:, conditions_left_negative_0ms), 2);
    right_heel_x_pos_response_civ_left_negative_0ms = tinv(0.975, sum(conditions_left_negative_0ms)-1) * std(right_heel_x_pos_response(:, conditions_left_negative_0ms), 1, 2)/sqrt(sum(conditions_left_negative_0ms));
    right_heel_x_pos_response_mean_left_positive_150ms = mean(right_heel_x_pos_response(:, conditions_left_positive_150ms), 2);
    right_heel_x_pos_response_civ_left_positive_150ms = tinv(0.975, sum(conditions_left_positive_150ms)-1) * std(right_heel_x_pos_response(:, conditions_left_positive_150ms), 1, 2)/sqrt(sum(conditions_left_positive_150ms));
    right_heel_x_pos_response_mean_left_negative_150ms = mean(right_heel_x_pos_response(:, conditions_left_negative_150ms), 2);
    right_heel_x_pos_response_civ_left_negative_150ms = tinv(0.975, sum(conditions_left_negative_150ms)-1) * std(right_heel_x_pos_response(:, conditions_left_negative_150ms), 1, 2)/sqrt(sum(conditions_left_negative_150ms));
    right_heel_x_pos_response_mean_left_positive_450ms = mean(right_heel_x_pos_response(:, conditions_left_positive_450ms), 2);
    right_heel_x_pos_response_civ_left_positive_450ms = tinv(0.975, sum(conditions_left_positive_450ms)-1) * std(right_heel_x_pos_response(:, conditions_left_positive_450ms), 1, 2)/sqrt(sum(conditions_left_positive_450ms));
    right_heel_x_pos_response_mean_left_negative_450ms = mean(right_heel_x_pos_response(:, conditions_left_negative_450ms), 2);
    right_heel_x_pos_response_civ_left_negative_450ms = tinv(0.975, sum(conditions_left_negative_450ms)-1) * std(right_heel_x_pos_response(:, conditions_left_negative_450ms), 1, 2)/sqrt(sum(conditions_left_negative_450ms));
    
    left_heel_y_pos_response_mean_right_positive_0ms = mean(left_heel_y_pos_response(:, conditions_right_positive_0ms), 2);
    left_heel_y_pos_response_civ_right_positive_0ms = tinv(0.975, sum(conditions_right_positive_0ms)-1) * std(left_heel_y_pos_response(:, conditions_right_positive_0ms), 1, 2)/sqrt(sum(conditions_right_positive_0ms));
    left_heel_y_pos_response_mean_right_negative_0ms = mean(left_heel_y_pos_response(:, conditions_right_negative_0ms), 2);
    left_heel_y_pos_response_civ_right_negative_0ms = tinv(0.975, sum(conditions_right_negative_0ms)-1) * std(left_heel_y_pos_response(:, conditions_right_negative_0ms), 1, 2)/sqrt(sum(conditions_right_negative_0ms));
    left_heel_y_pos_response_mean_right_positive_150ms = mean(left_heel_y_pos_response(:, conditions_right_positive_150ms), 2);
    left_heel_y_pos_response_civ_right_positive_150ms = tinv(0.975, sum(conditions_right_positive_150ms)-1) * std(left_heel_y_pos_response(:, conditions_right_positive_150ms), 1, 2)/sqrt(sum(conditions_right_positive_150ms));
    left_heel_y_pos_response_mean_right_negative_150ms = mean(left_heel_y_pos_response(:, conditions_right_negative_150ms), 2);
    left_heel_y_pos_response_civ_right_negative_150ms = tinv(0.975, sum(conditions_right_negative_150ms)-1) * std(left_heel_y_pos_response(:, conditions_right_negative_150ms), 1, 2)/sqrt(sum(conditions_right_negative_150ms));
    left_heel_y_pos_response_mean_right_positive_450ms = mean(left_heel_y_pos_response(:, conditions_right_positive_450ms), 2);
    left_heel_y_pos_response_civ_right_positive_450ms = tinv(0.975, sum(conditions_right_positive_450ms)-1) * std(left_heel_y_pos_response(:, conditions_right_positive_450ms), 1, 2)/sqrt(sum(conditions_right_positive_450ms));
    left_heel_y_pos_response_mean_right_negative_450ms = mean(left_heel_y_pos_response(:, conditions_right_negative_450ms), 2);
    left_heel_y_pos_response_civ_right_negative_450ms = tinv(0.975, sum(conditions_right_negative_450ms)-1) * std(left_heel_y_pos_response(:, conditions_right_negative_450ms), 1, 2)/sqrt(sum(conditions_right_negative_450ms));
    
    right_heel_y_pos_response_mean_left_positive_0ms = mean(right_heel_y_pos_response(:, conditions_left_positive_0ms), 2);
    right_heel_y_pos_response_civ_left_positive_0ms = tinv(0.975, sum(conditions_left_positive_0ms)-1) * std(right_heel_y_pos_response(:, conditions_left_positive_0ms), 1, 2)/sqrt(sum(conditions_left_positive_0ms));
    right_heel_y_pos_response_mean_left_negative_0ms = mean(right_heel_y_pos_response(:, conditions_left_negative_0ms), 2);
    right_heel_y_pos_response_civ_left_negative_0ms = tinv(0.975, sum(conditions_left_negative_0ms)-1) * std(right_heel_y_pos_response(:, conditions_left_negative_0ms), 1, 2)/sqrt(sum(conditions_left_negative_0ms));
    right_heel_y_pos_response_mean_left_positive_150ms = mean(right_heel_y_pos_response(:, conditions_left_positive_150ms), 2);
    right_heel_y_pos_response_civ_left_positive_150ms = tinv(0.975, sum(conditions_left_positive_150ms)-1) * std(right_heel_y_pos_response(:, conditions_left_positive_150ms), 1, 2)/sqrt(sum(conditions_left_positive_150ms));
    right_heel_y_pos_response_mean_left_negative_150ms = mean(right_heel_y_pos_response(:, conditions_left_negative_150ms), 2);
    right_heel_y_pos_response_civ_left_negative_150ms = tinv(0.975, sum(conditions_left_negative_150ms)-1) * std(right_heel_y_pos_response(:, conditions_left_negative_150ms), 1, 2)/sqrt(sum(conditions_left_negative_150ms));
    right_heel_y_pos_response_mean_left_positive_450ms = mean(right_heel_y_pos_response(:, conditions_left_positive_450ms), 2);
    right_heel_y_pos_response_civ_left_positive_450ms = tinv(0.975, sum(conditions_left_positive_450ms)-1) * std(right_heel_y_pos_response(:, conditions_left_positive_450ms), 1, 2)/sqrt(sum(conditions_left_positive_450ms));
    right_heel_y_pos_response_mean_left_negative_450ms = mean(right_heel_y_pos_response(:, conditions_left_negative_450ms), 2);
    right_heel_y_pos_response_civ_left_negative_450ms = tinv(0.975, sum(conditions_left_negative_450ms)-1) * std(right_heel_y_pos_response(:, conditions_left_negative_450ms), 1, 2)/sqrt(sum(conditions_left_negative_450ms));
    
    % step and stim response
    step_response_mean_left_positive = mean(step_response_left_positive);
    step_response_civ_left_positive = tinv(0.975, length(step_response_left_positive)-1) * std(step_response_left_positive)/sqrt(length(step_response_left_positive));
    step_response_mean_left_negative = mean(step_response_left_negative);
    step_response_civ_left_negative = tinv(0.975, length(step_response_left_negative)-1) * std(step_response_left_negative)/sqrt(length(step_response_left_negative));
    step_response_mean_right_positive = mean(step_response_right_positive);
    step_response_civ_right_positive = tinv(0.975, length(step_response_right_positive)-1) * std(step_response_right_positive)/sqrt(length(step_response_right_positive));
    step_response_mean_right_negative = mean(step_response_right_negative);
    step_response_civ_right_negative = tinv(0.975, length(step_response_right_negative)-1) * std(step_response_right_negative)/sqrt(length(step_response_right_negative));
    stim_response_mean_left_positive = mean(stim_response_left_positive);
    stim_response_civ_left_positive = tinv(0.975, length(stim_response_left_positive)-1) * std(stim_response_left_positive)/sqrt(length(stim_response_left_positive));
    stim_response_mean_left_negative = mean(stim_response_left_negative);
    stim_response_civ_left_negative = tinv(0.975, length(stim_response_left_negative)-1) * std(stim_response_left_negative)/sqrt(length(stim_response_left_negative));
    stim_response_mean_right_positive = mean(stim_response_right_positive);
    stim_response_civ_right_positive = tinv(0.975, length(stim_response_right_positive)-1) * std(stim_response_right_positive)/sqrt(length(stim_response_right_positive));
    stim_response_mean_right_negative = mean(stim_response_right_negative);
    stim_response_civ_right_negative = tinv(0.975, length(stim_response_right_negative)-1) * std(stim_response_right_negative)/sqrt(length(stim_response_right_negative));
    
    % stim start times
    stim_start_times_mean_left_0ms = mean(stim_start_time_relative_to_stretch(conditions_left_0ms));
    stim_start_times_mean_right_0ms = mean(stim_start_time_relative_to_stretch(conditions_right_0ms));
    stim_start_times_mean_left_150ms = mean(stim_start_time_relative_to_stretch(conditions_left_150ms));
    stim_start_times_mean_right_150ms = mean(stim_start_time_relative_to_stretch(conditions_right_150ms));
    stim_start_times_mean_left_450ms = mean(stim_start_time_relative_to_stretch(conditions_left_450ms));
    stim_start_times_mean_right_450ms = mean(stim_start_time_relative_to_stretch(conditions_right_450ms));
    
end

%% calculate response maxima
if calculate_response_extrema
    cop_response_extrema = zeros(number_of_stretches, 1);
    for i_stretch = 1 : number_of_stretches
        if conditions_left_positive_0ms(i_stretch)
            response = left_cop_x_response(:, i_stretch);
            extremum = min(response);
            cop_response_extrema(i_stretch) = extremum;
        elseif conditions_left_negative_0ms(i_stretch)
            response = left_cop_x_response(:, i_stretch);
            extremum = max(response);
            cop_response_extrema(i_stretch) = extremum;
        elseif conditions_right_positive_0ms(i_stretch)
            response = right_cop_x_response(:, i_stretch);
            extremum = min(response);
            cop_response_extrema(i_stretch) = extremum;
        elseif conditions_right_negative_0ms(i_stretch)
            response = right_cop_x_response(:, i_stretch);
            extremum = max(response);
            cop_response_extrema(i_stretch) = extremum;
        end
    end
end

%% do_cop_plots_absolute
if do_cop_plots_absolute_right
    % figure; axes; hold on
    % plot(time_normalized, left_cop_x_normalized(:, conditions_left_control))

    % figure; axes; hold on
    % plot(time_normalized, right_cop_x_normalized(:, conditions_right_control))

    % return

    figure; axes; hold on; title('right foot medial-lateral CoP, 0ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, right_cop_x_mean_right_control, right_cop_x_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, right_cop_x_mean_right_positive_0ms, right_cop_x_civ_right_positive_0ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, right_cop_x_mean_right_negative_0ms, right_cop_x_civ_right_negative_0ms, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean]); 
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    
    if phase_dependent
        figure; axes; hold on; title('right foot medial-lateral CoP, 150ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, right_cop_x_mean_right_control, right_cop_x_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, right_cop_x_mean_right_positive_150ms, right_cop_x_civ_right_positive_150ms, {'color', color_right_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, right_cop_x_mean_right_negative_150ms, right_cop_x_civ_right_negative_150ms, {'color', color_right_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')

        figure; axes; hold on; title('left foot medial-lateral CoP, 450ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, left_cop_x_mean_left_control, left_cop_x_civ_left_control, {'color', color_left_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, left_cop_x_mean_left_positive_450ms, left_cop_x_civ_left_positive_450ms, {'color', color_left_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, left_cop_x_mean_left_negative_450ms, left_cop_x_civ_left_negative_450ms, {'color', color_left_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHEAST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    end
%     % single ones
%     figure; axes; hold on; title('left foot medial-lateral CoP, 450ms')
%     plot(time_normalized, left_cop_x_normalized(:, conditions_left_control), 'color', color_left_control, 'linewidth', 1);
%     plot(time_normalized, left_cop_x_normalized(:, conditions_left_positive_450ms), 'color', color_left_positive, 'linewidth', 1);
%     plot(time_normalized, left_cop_x_normalized(:, conditions_left_negative_450ms), 'color', color_left_negative, 'linewidth', 1);
%     xlabel('time')
% 
%     figure; axes; hold on; title('right foot medial-lateral CoP, 0ms')
%     plot(time_normalized, right_cop_x_normalized(:, conditions_right_control), 'color', color_right_control, 'linewidth', 1);
%     plot(time_normalized, right_cop_x_normalized(:, conditions_right_positive_0ms), 'color', color_left_positive, 'linewidth', 1);
%     plot(time_normalized, right_cop_x_normalized(:, conditions_right_negative_0ms), 'color', color_left_negative, 'linewidth', 1);
%     
%     figure; axes; hold on; title('right foot medial-lateral CoP, 150ms')
%     plot(time_normalized, right_cop_x_normalized(:, conditions_right_control), 'color', color_right_control, 'linewidth', 1);
%     plot(time_normalized, right_cop_x_normalized(:, conditions_right_positive_150ms), 'color', color_left_positive, 'linewidth', 1);
%     plot(time_normalized, right_cop_x_normalized(:, conditions_right_negative_150ms), 'color', color_left_negative, 'linewidth', 1);
%     xlabel('time')
end

if do_cop_plots_absolute_left
    figure; axes; hold on; title('left foot medial-lateral CoP, 0ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, left_cop_x_mean_left_control, left_cop_x_civ_left_control, {'color', color_left_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, left_cop_x_mean_left_positive_0ms, left_cop_x_civ_left_positive_0ms, {'color', color_left_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, left_cop_x_mean_left_negative_0ms, left_cop_x_civ_left_negative_0ms, {'color', color_left_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean]); 
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
end

%% do_heel_plots_absolute
if do_heel_plots_absolute_left
    figure; axes; hold on; title('left foot medial-lateral heel, 0ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, left_heel_x_pos_mean_right_control, left_heel_x_pos_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, left_heel_x_pos_mean_right_positive_0ms, left_heel_x_pos_civ_right_positive_0ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, left_heel_x_pos_mean_right_negative_0ms, left_heel_x_pos_civ_right_negative_0ms, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    
    if phase_dependent
        figure; axes; hold on; title('left foot medial-lateral heel, 150ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, left_heel_x_pos_mean_right_control, left_heel_x_pos_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, left_heel_x_pos_mean_right_positive_150ms, left_heel_x_pos_civ_right_positive_150ms, {'color', color_right_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, left_heel_x_pos_mean_right_negative_150ms, left_heel_x_pos_civ_right_negative_150ms, {'color', color_right_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')

        figure; axes; hold on; title('right foot medial-lateral heel, 450ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, right_heel_x_pos_mean_left_control, right_heel_x_pos_civ_left_control, {'color', color_left_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, right_heel_x_pos_mean_left_positive_450ms, right_heel_x_pos_civ_left_positive_450ms, {'color', color_left_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, right_heel_x_pos_mean_left_negative_450ms, right_heel_x_pos_civ_left_negative_450ms, {'color', color_left_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    end
%     % shaded error bars
%     figure; axes; hold on; title('left foot medial-lateral heel marker')
%     shadedErrorBar(time_normalized, left_heel_x_mean_right_control, left_heel_x_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
%     shadedErrorBar(time_normalized, left_heel_x_mean_right_positive, left_heel_x_civ_right_positive, {'color', color_right_positive, 'linewidth', 5}, 1);
%     shadedErrorBar(time_normalized, left_heel_x_mean_right_negative, left_heel_x_civ_right_negative, {'color', color_right_negative, 'linewidth', 5}, 1);
%     xlabel('time')
% 
%     figure; axes; hold on; title('right foot medial-lateral heel marker')
%     shadedErrorBar(time_normalized, right_heel_x_mean_left_control, right_heel_x_civ_left_control, {'color', color_left_control, 'linewidth', 5}, 1);
%     shadedErrorBar(time_normalized, right_heel_x_mean_left_positive, right_heel_x_civ_left_positive, {'color', color_left_positive, 'linewidth', 5}, 1);
%     shadedErrorBar(time_normalized, right_heel_x_mean_left_negative, right_heel_x_civ_left_negative, {'color', color_left_negative, 'linewidth', 5}, 1);
%     xlabel('time')
%     
%     % x
%     figure; axes; hold on; title('left stance foot - x')
%     plot(time_normalized, left_heel_x_mean_right_control, 'linewidth', 5);
%     plot(time_normalized, left_heel_x_mean_left_control, 'linewidth', 5);
%     
%     figure; axes; hold on; title('right stance foot - x')
%     plot(time_normalized, right_heel_x_mean_right_control, 'linewidth', 5);
%     plot(time_normalized, right_heel_x_mean_left_control, 'linewidth', 5);
    
%     % y
%     figure; axes; hold on; title('left stance foot - y')
%     plot(time_normalized, left_heel_y_pos_mean_right_control, 'linewidth', 5);
%     plot(time_normalized, left_heel_y_pos_mean_left_control, 'linewidth', 5);
%     
%     figure; axes; hold on; title('right stance foot - y')
%     plot(time_normalized, right_heel_y_pos_mean_right_control, 'linewidth', 5);
%     plot(time_normalized, right_heel_y_pos_mean_left_control, 'linewidth', 5);
    
end

if do_heel_plots_absolute_right
    figure; axes; hold on; title('right foot medial-lateral heel, 0ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, right_heel_x_pos_mean_left_control, right_heel_x_pos_civ_left_control, {'color', color_right_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, right_heel_x_pos_mean_left_positive_0ms, right_heel_x_pos_civ_left_positive_0ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, right_heel_x_pos_mean_left_negative_0ms, right_heel_x_pos_civ_left_negative_0ms, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
end

%% do_emg_plots_absolute
if do_emg_plots_absolute
    % left gluteus medius
    left_glutmed_right_0ms_figure = figure; axes; hold on; title('left gluteus medius EMG, right foot stance, 0ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, left_glutmed_emg_mean_right_control, left_glutmed_emg_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, left_glutmed_emg_mean_right_positive_0ms, left_glutmed_emg_civ_right_positive_0ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, left_glutmed_emg_mean_right_negative_0ms, left_glutmed_emg_civ_right_negative_0ms, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    
    if phase_dependent
        left_glutmed_right_150ms_figure = figure; axes; hold on; title('left gluteus medius EMG, right foot stance, 150ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, left_glutmed_emg_mean_right_control, left_glutmed_emg_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, left_glutmed_emg_mean_right_positive_150ms, left_glutmed_emg_civ_right_positive_150ms, {'color', color_right_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, left_glutmed_emg_mean_right_negative_150ms, left_glutmed_emg_civ_right_negative_150ms, {'color', color_right_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')

        left_glutmed_left_450ms_figure = figure; axes; hold on; title('left gluteus medius EMG, left foot stance, 450ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, left_glutmed_emg_mean_left_control, left_glutmed_emg_civ_left_control, {'color', color_right_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, left_glutmed_emg_mean_left_positive_450ms, left_glutmed_emg_civ_left_positive_450ms, {'color', color_right_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, left_glutmed_emg_mean_left_negative_450ms, left_glutmed_emg_civ_left_negative_450ms, {'color', color_right_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
    end
    
    % left tibialis anterior
    left_tibiant_right_0ms_figure = figure; axes; hold on; title('left tibialis anterior EMG, right foot stance, 0ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, left_tibiant_emg_mean_right_control, left_tibiant_emg_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, left_tibiant_emg_mean_right_positive_0ms, left_tibiant_emg_civ_right_positive_0ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, left_tibiant_emg_mean_right_negative_0ms, left_tibiant_emg_civ_right_negative_0ms, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    
    if phase_dependent
        left_tibiant_right_150ms_figure = figure; axes; hold on; title('left tibialis anterior EMG, right foot stance, 150ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, left_tibiant_emg_mean_right_control, left_tibiant_emg_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, left_tibiant_emg_mean_right_positive_150ms, left_tibiant_emg_civ_right_positive_150ms, {'color', color_right_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, left_tibiant_emg_mean_right_negative_150ms, left_tibiant_emg_civ_right_negative_150ms, {'color', color_right_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')

        left_tibiant_left_450ms_figure = figure; axes; hold on; title('left tibialis anterior EMG, left foot stance, 450ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, left_tibiant_emg_mean_left_control, left_tibiant_emg_civ_left_control, {'color', color_right_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, left_tibiant_emg_mean_left_positive_450ms, left_tibiant_emg_civ_left_positive_450ms, {'color', color_right_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, left_tibiant_emg_mean_left_negative_450ms, left_tibiant_emg_civ_left_negative_450ms, {'color', color_right_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
    end
    
    % left peroneus longus
    left_perolng_right_0ms_figure = figure; axes; hold on; title('left peroneus longus EMG, right foot stance, 0ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, left_perolng_emg_mean_right_control, left_perolng_emg_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, left_perolng_emg_mean_right_positive_0ms, left_perolng_emg_civ_right_positive_0ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, left_perolng_emg_mean_right_negative_0ms, left_perolng_emg_civ_right_negative_0ms, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    
    if phase_dependent
        left_perolng_right_150ms_figure = figure; axes; hold on; title('left peroneus longus EMG, right foot stance, 150ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, left_perolng_emg_mean_right_control, left_perolng_emg_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, left_perolng_emg_mean_right_positive_150ms, left_perolng_emg_civ_right_positive_150ms, {'color', color_right_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, left_perolng_emg_mean_right_negative_150ms, left_perolng_emg_civ_right_negative_150ms, {'color', color_right_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')

        left_perolng_left_450ms_figure = figure; axes; hold on; title('left peroneus longus EMG, left foot stance, 450ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, left_perolng_emg_mean_left_control, left_perolng_emg_civ_left_control, {'color', color_right_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, left_perolng_emg_mean_left_positive_450ms, left_perolng_emg_civ_left_positive_450ms, {'color', color_right_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, left_perolng_emg_mean_left_negative_450ms, left_perolng_emg_civ_left_negative_450ms, {'color', color_right_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
    end
    
    % right gluteus medius
    right_glutmed_left_0ms_figure = figure; axes; hold on; title('right gluteus medius EMG, right foot stance, 0ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, right_glutmed_emg_mean_right_control, right_glutmed_emg_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, right_glutmed_emg_mean_right_positive_0ms, right_glutmed_emg_civ_right_positive_0ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, right_glutmed_emg_mean_right_negative_0ms, right_glutmed_emg_civ_right_negative_0ms, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    
    if phase_dependent
        right_glutmed_left_150ms_figure = figure; axes; hold on; title('right gluteus medius EMG, right foot stance, 150ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, right_glutmed_emg_mean_right_control, right_glutmed_emg_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, right_glutmed_emg_mean_right_positive_150ms, right_glutmed_emg_civ_right_positive_150ms, {'color', color_right_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, right_glutmed_emg_mean_right_negative_150ms, right_glutmed_emg_civ_right_negative_150ms, {'color', color_right_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')

        right_glutmed_right_450ms_figure = figure; axes; hold on; title('right gluteus medius EMG, right foot stance, 450ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, right_glutmed_emg_mean_left_control, right_glutmed_emg_civ_left_control, {'color', color_right_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, right_glutmed_emg_mean_left_positive_450ms, right_glutmed_emg_civ_left_positive_450ms, {'color', color_right_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, right_glutmed_emg_mean_left_negative_450ms, right_glutmed_emg_civ_left_negative_450ms, {'color', color_right_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
    end
    
    % right tibialis anterior
    right_tibiant_left_0ms_figure = figure; axes; hold on; title('right tibialis anterior EMG, right foot stance, 0ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, right_tibiant_emg_mean_right_control, right_tibiant_emg_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, right_tibiant_emg_mean_right_positive_0ms, right_tibiant_emg_civ_right_positive_0ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, right_tibiant_emg_mean_right_negative_0ms, right_tibiant_emg_civ_right_negative_0ms, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    
    if phase_dependent
        right_tibiant_left_150ms_figure = figure; axes; hold on; title('right tibialis anterior EMG, right foot stance, 150ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, right_tibiant_emg_mean_right_control, right_tibiant_emg_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, right_tibiant_emg_mean_right_positive_150ms, right_tibiant_emg_civ_right_positive_150ms, {'color', color_right_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, right_tibiant_emg_mean_right_negative_150ms, right_tibiant_emg_civ_right_negative_150ms, {'color', color_right_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')

        right_tibiant_right_450ms_figure = figure; axes; hold on; title('right tibialis anterior EMG, right foot stance, 450ms'); set(gca, 'Fontsize', 12)
        control_plot = shadedErrorBar(time_normalized, right_tibiant_emg_mean_left_control, right_tibiant_emg_civ_left_control, {'color', color_right_control, 'linewidth', 5}, 1);
        positive_plot = shadedErrorBar(time_normalized, right_tibiant_emg_mean_left_positive_450ms, right_tibiant_emg_civ_left_positive_450ms, {'color', color_right_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, right_tibiant_emg_mean_left_negative_450ms, right_tibiant_emg_civ_left_negative_450ms, {'color', color_right_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
    end
    
    % right peroneus longus
    right_perolng_left_0ms_figure = figure; axes; hold on; title('right peroneus longus EMG, right foot stance, 0ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, right_perolng_emg_mean_right_control, right_perolng_emg_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, right_perolng_emg_mean_right_positive_0ms, right_perolng_emg_civ_right_positive_0ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, right_perolng_emg_mean_right_negative_0ms, right_perolng_emg_civ_right_negative_0ms, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    
    if phase_dependent
    right_perolng_left_150ms_figure = figure; axes; hold on; title('right peroneus longus EMG, right foot stance, 150ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, right_perolng_emg_mean_right_control, right_perolng_emg_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, right_perolng_emg_mean_right_positive_150ms, right_perolng_emg_civ_right_positive_150ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, right_perolng_emg_mean_right_negative_150ms, right_perolng_emg_civ_right_negative_150ms, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    
    right_perolng_right_450ms_figure = figure; axes; hold on; title('right peroneus longus EMG, right foot stance, 450ms'); set(gca, 'Fontsize', 12)
    control_plot = shadedErrorBar(time_normalized, right_perolng_emg_mean_left_control, right_perolng_emg_civ_left_control, {'color', color_right_control, 'linewidth', 5}, 1);
    positive_plot = shadedErrorBar(time_normalized, right_perolng_emg_mean_left_positive_450ms, right_perolng_emg_civ_left_positive_450ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, right_perolng_emg_mean_left_negative_450ms, right_perolng_emg_civ_left_negative_450ms, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
    this_legend = legend([control_plot.mainLine positive_plot.mainLine negative_plot.mainLine], 'control', 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')    
    end
    
%     saveas(left_glutmed_right_0ms_figure, 'left_glutmed_right_0ms_figure', 'eps2c');
%     saveas(left_glutmed_right_150ms_figure, 'left_glutmed_right_150ms_figure', 'eps2c');
%     saveas(left_glutmed_left_450ms_figure, 'left_glutmed_right_450ms_figure', 'eps2c');
%     saveas(left_tibiant_right_0ms_figure, 'left_tibiant_right_0ms_figure', 'eps2c');
%     saveas(left_tibiant_right_150ms_figure, 'left_tibiant_right_150ms_figure', 'eps2c');
%     saveas(left_tibiant_left_450ms_figure, 'left_tibiant_right_450ms_figure', 'eps2c');
%     saveas(left_perolng_right_0ms_figure, 'left_perolng_right_0ms_figure', 'eps2c');
%     saveas(left_perolng_right_150ms_figure, 'left_perolng_right_150ms_figure', 'eps2c');
%     saveas(left_perolng_left_450ms_figure, 'left_perolng_right_450ms_figure', 'eps2c');
%     saveas(right_glutmed_left_0ms_figure, 'right_glutmed_left_0ms_figure', 'eps2c');
%     saveas(right_glutmed_left_150ms_figure, 'right_glutmed_left_150ms_figure', 'eps2c');
%     saveas(right_glutmed_right_450ms_figure, 'right_glutmed_left_450ms_figure', 'eps2c');
%     saveas(right_tibiant_left_0ms_figure, 'right_tibiant_left_0ms_figure', 'eps2c');
%     saveas(right_tibiant_left_150ms_figure, 'right_tibiant_left_150ms_figure', 'eps2c');
%     saveas(right_tibiant_right_450ms_figure, 'right_tibiant_left_450ms_figure', 'eps2c');
%     saveas(right_perolng_left_0ms_figure, 'right_perolng_left_0ms_figure', 'eps2c');
%     saveas(right_perolng_left_150ms_figure, 'right_perolng_left_150ms_figure', 'eps2c');
%     saveas(right_perolng_right_450ms_figure, 'right_perolng_left_450ms_figure', 'eps2c');
    
%     % single ones
%     figure; axes; hold on; title('left gluteus medius EMG, 0ms')
%     plot(time_normalized, left_glutmed_emg_normalized(:, conditions_right_control), 'color', color_right_control, 'linewidth', 1);
%     plot(time_normalized, left_glutmed_emg_normalized(:, conditions_right_positive_0ms), 'color', color_left_positive, 'linewidth', 1);
%     plot(time_normalized, left_glutmed_emg_normalized(:, conditions_right_negative_0ms), 'color', color_left_negative, 'linewidth', 1);
%     
%     figure; axes; hold on; title('left gluteus medius EMG, 0ms')
%     plot(time_normalized, left_tibiant_emg_normalized(:, conditions_right_control), 'color', color_right_control, 'linewidth', 1);
%     plot(time_normalized, left_tibiant_emg_normalized(:, conditions_right_positive_0ms), 'color', color_left_positive, 'linewidth', 1);
%     plot(time_normalized, left_tibiant_emg_normalized(:, conditions_right_negative_0ms), 'color', color_left_negative, 'linewidth', 1);
%     
%     figure; axes; hold on; title('left gluteus medius EMG, 0ms')
%     plot(time_normalized, left_perolng_emg_normalized(:, conditions_right_control), 'color', color_right_control, 'linewidth', 1);
%     plot(time_normalized, left_perolng_emg_normalized(:, conditions_right_positive_0ms), 'color', color_left_positive, 'linewidth', 1);
%     plot(time_normalized, left_perolng_emg_normalized(:, conditions_right_negative_0ms), 'color', color_left_negative, 'linewidth', 1);
    
end

%% do_cop_plots_response
if do_cop_plots_response_right
    figure; axes; hold on; title('right foot medial-lateral CoP response, 0ms'); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, right_cop_x_response_mean_right_positive_0ms, right_cop_x_response_civ_right_positive_0ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, right_cop_x_response_mean_right_negative_0ms, right_cop_x_response_civ_right_positive_0ms, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean], 'ylim', [-0.01 0.01]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    
    if phase_dependent
        figure; axes; hold on; title('right foot medial-lateral CoP response, 150ms'); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, right_cop_x_response_mean_right_positive_150ms, right_cop_x_response_civ_right_positive_150ms, {'color', color_right_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, right_cop_x_response_mean_right_negative_150ms, right_cop_x_response_civ_right_positive_150ms, {'color', color_right_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean], 'ylim', [-0.01 0.01]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')

        figure; axes; hold on; title('left foot medial-lateral CoP response, 450ms'); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, left_cop_x_response_mean_left_positive_450ms, left_cop_x_response_civ_left_positive_450ms, {'color', color_left_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, left_cop_x_response_mean_left_negative_450ms, left_cop_x_response_civ_left_positive_450ms, {'color', color_left_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean], 'ylim', [-0.01 0.01]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    end
    
%     figure; axes; hold on; title('right foot medial-lateral CoP - response')
%     shadedErrorBar(time_normalized, right_cop_x_response_mean_right_positive, right_cop_x_response_civ_right_positive, {'color', color_right_positive, 'linewidth', 5}, 1);
%     shadedErrorBar(time_normalized, right_cop_x_response_mean_right_negative, right_cop_x_response_civ_right_negative, {'color', color_right_negative, 'linewidth', 5}, 1);
%     xlabel('time')
end

if do_cop_plots_response_left
    figure; axes; hold on; title('left foot medial-lateral CoP response, 0ms'); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, left_cop_x_response_mean_left_positive_0ms, left_cop_x_response_civ_left_positive_0ms, {'color', color_left_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, left_cop_x_response_mean_left_negative_0ms, left_cop_x_response_civ_left_positive_0ms, {'color', color_left_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean], 'ylim', [-0.01 0.01]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
end

%% do_heel_plots_response
if do_heel_plots_response_left
    figure; axes; hold on; title('right foot medial-lateral heel response, 0ms'); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, left_heel_x_pos_response_mean_right_positive_0ms, left_heel_x_pos_response_civ_right_positive_0ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, left_heel_x_pos_response_mean_right_negative_0ms, left_heel_x_pos_response_civ_right_positive_0ms, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    
    if phase_dependent
        figure; axes; hold on; title('right foot medial-lateral heel response, 150ms'); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, left_heel_x_pos_response_mean_right_positive_150ms, left_heel_x_pos_response_civ_right_positive_150ms, {'color', color_right_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, left_heel_x_pos_response_mean_right_negative_150ms, left_heel_x_pos_response_civ_right_positive_150ms, {'color', color_right_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')

        figure; axes; hold on; title('left foot medial-lateral heel response, 450ms'); set(gca, 'Fontsize', 12)
        positive_plot = shadedErrorBar(time_normalized, right_heel_x_pos_response_mean_left_positive_450ms, right_heel_x_pos_response_civ_left_positive_450ms, {'color', color_left_positive, 'linewidth', 5}, 1);
        negative_plot = shadedErrorBar(time_normalized, right_heel_x_pos_response_mean_left_negative_450ms, right_heel_x_pos_response_civ_left_positive_450ms, {'color', color_left_negative, 'linewidth', 5}, 1);
        xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
        this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
        set(this_legend, 'Location', 'NORTHWEST')
        xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
        text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    end
end

if do_heel_plots_response_right
    figure; axes; hold on; title('right foot medial-lateral heel response, 0ms'); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, right_heel_x_pos_response_mean_left_positive_0ms, right_heel_x_pos_response_civ_left_positive_0ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, right_heel_x_pos_response_mean_left_negative_0ms, right_heel_x_pos_response_civ_left_positive_0ms, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
    
    figure; axes; hold on; title('left foot medial-lateral heel response, 450ms'); set(gca, 'Fontsize', 12)
    positive_plot = shadedErrorBar(time_normalized, right_heel_x_pos_response_mean_left_positive_450ms, right_heel_x_pos_response_civ_left_positive_450ms, {'color', color_left_positive, 'linewidth', 5}, 1);
    negative_plot = shadedErrorBar(time_normalized, right_heel_x_pos_response_mean_left_negative_450ms, right_heel_x_pos_response_civ_left_positive_450ms, {'color', color_left_negative, 'linewidth', 5}, 1);
    xlabel('time'); set(gca, 'xlim', [0, step_time_mean]);
    this_legend = legend([positive_plot.mainLine negative_plot.mainLine], 'illusion left', 'illusion right');
    set(this_legend, 'Location', 'NORTHWEST')
    xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), 'right \rightarrow', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right')
    text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), '\leftarrow left', 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left')
end

%% do_side_comparison_plots
if do_side_comparison_plots
%     figure; axes; hold on; title('medial-lateral CoP')
%     shadedErrorBar(time_normalized, left_cop_x_mean_left_control - left_cop_x_mean_left_control(1), left_cop_x_civ_left_control, {'color', color_left_control, 'linewidth', 5}, 1);
%     shadedErrorBar(time_normalized, left_cop_x_mean_left_positive - left_cop_x_mean_left_positive(1), left_cop_x_civ_left_positive, {'color', color_left_positive, 'linewidth', 5}, 1);
%     shadedErrorBar(time_normalized, left_cop_x_mean_left_negative - left_cop_x_mean_left_negative(1), left_cop_x_civ_left_negative, {'color', color_left_negative, 'linewidth', 5}, 1);
%     plot(time_normalized, -right_cop_x_mean_right_control + right_cop_x_mean_right_control(1), 'color', color_left_control, 'linewidth', 5);
%     plot(time_normalized, -right_cop_x_mean_right_positive + right_cop_x_mean_right_positive(1), 'color', color_left_positive, 'linewidth', 5);
%     plot(time_normalized, -right_cop_x_mean_right_negative + right_cop_x_mean_right_negative(1), 'color', color_left_negative, 'linewidth', 5);
%     xlabel('time')
%     
%     figure; axes; hold on; title('left foot medial-lateral heel marker')
%     shadedErrorBar(time_normalized, left_heel_x_mean_right_control, left_heel_x_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
%     shadedErrorBar(time_normalized, left_heel_x_mean_right_positive, left_heel_x_civ_right_positive, {'color', color_right_positive, 'linewidth', 5}, 1);
%     shadedErrorBar(time_normalized, left_heel_x_mean_right_negative, left_heel_x_civ_right_negative, {'color', color_right_negative, 'linewidth', 5}, 1);
%     plot(time_normalized, -right_heel_x_mean_left_control, 'color', color_right_control, 'linewidth', 5);
%     plot(time_normalized, -right_heel_x_mean_left_positive, 'color', color_right_positive, 'linewidth', 5);
%     plot(time_normalized, -right_heel_x_mean_left_negative, 'color', color_right_negative, 'linewidth', 5);
%     xlabel('time')
    
    figure; axes; hold on; title('left foot medial-lateral CoP - response')
    shadedErrorBar(time_normalized, left_cop_x_response_mean_left_positive_0ms, left_cop_x_response_civ_left_positive_0ms, {'color', color_left_positive, 'linewidth', 5}, 1);
    shadedErrorBar(time_normalized, left_cop_x_response_mean_left_negative_0ms, left_cop_x_response_civ_left_negative_0ms, {'color', color_left_negative, 'linewidth', 5}, 1);
    shadedErrorBar(time_normalized, -right_cop_x_response_mean_right_positive_0ms, right_cop_x_response_civ_right_positive_0ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    shadedErrorBar(time_normalized, -right_cop_x_response_mean_right_negative_0ms, right_cop_x_response_civ_right_negative_0ms, {'color', color_right_negative, 'linewidth', 5}, 1);
%     plot(time_normalized, -right_cop_x_response_mean_right_positive, 'color', color_right_positive, 'linewidth', 5);
%     plot(time_normalized, -right_cop_x_response_mean_right_negative, 'color', color_right_negative, 'linewidth', 5);
    xlabel('time')
    
    figure; axes; hold on; title('left foot medial-lateral heel marker')
    shadedErrorBar(time_normalized, left_heel_x_pos_response_mean_right_positive_0ms, left_heel_x_pos_response_civ_right_positive_0ms, {'color', color_right_positive, 'linewidth', 5}, 1);
    shadedErrorBar(time_normalized, left_heel_x_pos_response_mean_right_negative_0ms, left_heel_x_pos_response_civ_right_negative_0ms, {'color', color_right_negative, 'linewidth', 5}, 1);
    shadedErrorBar(time_normalized, -right_heel_x_pos_response_mean_left_positive_0ms, right_heel_x_pos_response_civ_left_positive_0ms, {'color', color_left_positive, 'linewidth', 5}, 1);
    shadedErrorBar(time_normalized, -right_heel_x_pos_response_mean_left_negative_0ms, right_heel_x_pos_response_civ_left_negative_0ms, {'color', color_left_negative, 'linewidth', 5}, 1);
    xlabel('time')
    
end

%% do_response_extrema_plots
if do_response_extrema_plots
    figure; axes; hold on; title('response extrema')
    plot(cop_response_extrema(conditions_left_positive_0ms), 'x-', 'linewidth', 2);
    plot(cop_response_extrema(conditions_left_negative_0ms), 'x-', 'linewidth', 2);
    plot(cop_response_extrema(conditions_right_positive_0ms), 'x-', 'linewidth', 2);
    plot(cop_response_extrema(conditions_right_negative_0ms), 'x-', 'linewidth', 2);
    legend('left positive', 'left negative', 'right positive', 'right negative')
end

%% do step response plots
if do_step_response_plot
    figure; axes; hold on; title('step response');
    plot(step_response_mean_left_positive + [-1 1]*step_response_civ_left_positive, [1 1]-0.1, 'color', color_left_positive, 'linewidth', 15);
    plot(step_response_mean_left_negative + [-1 1]*step_response_civ_left_negative, [2 2]-0.1, 'color', color_left_negative, 'linewidth', 15);
    plot(step_response_mean_right_positive + [-1 1]*step_response_civ_right_positive, [3 3]-0.1, 'color', color_right_positive, 'linewidth', 15);
    plot(step_response_mean_right_negative + [-1 1]*step_response_civ_right_negative, [4 4]-0.1, 'color', color_right_negative, 'linewidth', 15);
    set(gca, 'ytick', 1:4, 'yticklabels', {'left, positive', 'left, negative', 'right, positive', 'right, negative'});
    set(gca, 'ylim', [0.5 4.5], 'xlim', [-0.015 0.015]);
    
    
end

%% do stim response plots
if do_stim_response_plot
%     figure; axes; hold on; title('stim response');
    plot(stim_response_mean_left_positive + [-1 1]*stim_response_civ_left_positive, [1 1]+0.1, 'color', color_left_positive, 'linewidth', 15);
    plot(stim_response_mean_left_negative + [-1 1]*stim_response_civ_left_negative, [2 2]+0.1, 'color', color_left_negative, 'linewidth', 15);
    plot(stim_response_mean_right_positive + [-1 1]*stim_response_civ_right_positive, [3 3]+0.1, 'color', color_right_positive, 'linewidth', 15);
    plot(stim_response_mean_right_negative + [-1 1]*stim_response_civ_right_negative, [4 4]+0.1, 'color', color_right_negative, 'linewidth', 15);
    set(gca, 'ytick', 1:4, 'yticklabels', {'left, positive', 'left, negative', 'right, positive', 'right, negative'});
    set(gca, 'ylim', [0.5 4.5], 'xlim', [-0.015 0.015]);
end

%% do_stim_start_time_histograms
if do_stim_start_time_histograms
    stim_start_times_right_0ms = stim_start_time_relative_to_stretch(conditions_right_0ms);
    stim_start_times_right_150ms = stim_start_time_relative_to_stretch(conditions_right_150ms);
    stim_start_times_left_450ms = stim_start_time_relative_to_stretch(conditions_left_450ms);
    
    figure; axes; hold on; title('stim start times');
    plot(find(conditions_right_0ms), stim_start_times_right_0ms, 'x')
    plot(find(conditions_right_150ms), stim_start_times_right_150ms, 'x')
    plot(find(conditions_left_450ms), stim_start_times_left_450ms, 'x')
    legend('0ms', '150ms', '450ms')
    
    
    figure; axes; hold on; title('stim start times');
    histogram(stim_start_time_relative_to_stretch(conditions_right_0ms))
    histogram(stim_start_time_relative_to_stretch(conditions_right_150ms))
    histogram(stim_start_time_relative_to_stretch(conditions_left_450ms))
    legend('0ms', '150ms', '450ms')
end


% save('MCA_DATA',
distFig






















