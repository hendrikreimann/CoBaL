% find the relevant stretches of data after the stimuli


%% Choose Data Type
visual_perturbation             = 0;
emg_present                     = 0;
% stimulus_type = 'gvs';
stimulus_type = 'visual';

view_totals_and_removals       =  1;
visualize_triggers              = 0;

%% prepare
wait_times = [0 0.150 0.450];
wait_time_labels = {'0ms', '150ms', '450ms'};
load subjectInfo.mat;

% trials_to_process = 1 : 23;
trials_to_process = [1:7 9:23];
% trials_to_process = 12 : 23;
% trials_to_process = [1:21];
% trials_to_process = 2;
% trials_to_process = 12 : 23;
trials_to_process = 1 : 21;
% trials_to_process = 4;

total_positive_steps = [];
total_negative_steps = [];
total_polarity_condition_list = [];
total_step_condition_list = [];

number_of_time_steps_normalized = 256;
swing_foot_fz_zero_threshold = 20; % threshold for counting a vertical force reading as zero, in Nm
swing_foot_zero_stretch_length_threshold = 20; % the number of zero indices in the vertical swing foot force has to be larger than this number
duration_until_nearest_future_heelstrike_threshold = 0.05; % a heelstrike should happen less than this long after a trigger

% left_heel_marker = 34;
% left_toes_marker = 35;
% right_heel_marker = 42;
% right_toes_marker = 43;
left_heel_marker = 33;
left_toes_marker = 34;
right_heel_marker = 42;
right_toes_marker = 43;

left_heel_marker_indices = reshape([(left_heel_marker - 1) * 3 + 1; (left_heel_marker - 1) * 3 + 2; (left_heel_marker - 1) * 3 + 3], 1, length(left_heel_marker)*3);
left_toes_marker_indices = reshape([(left_toes_marker - 1) * 3 + 1; (left_toes_marker - 1) * 3 + 2; (left_toes_marker - 1) * 3 + 3], 1, length(left_toes_marker)*3);
right_heel_marker_indices = reshape([(right_heel_marker - 1) * 3 + 1; (right_heel_marker - 1) * 3 + 2; (right_heel_marker - 1) * 3 + 3], 1, length(right_heel_marker)*3);
right_toes_marker_indices = reshape([(right_toes_marker - 1) * 3 + 1; (right_toes_marker - 1) * 3 + 2; (right_toes_marker - 1) * 3 + 3], 1, length(right_toes_marker)*3);


%% extract data
stretch_length_indices_forceplate = [];
condition_stance_foot_list = {};
condition_polarity_list = {};
condition_delay_list = {};
stim_start_time_relative_to_stretch = [];

for i_trial = trials_to_process
    % load data
    load(makeFileName(date, subject_id, 'walking', i_trial, 'forcePlateData'));
    load(makeFileName(date, subject_id, 'walking', i_trial, 'markerTrajectories'));
    load(makeFileName(date, subject_id, 'walking', i_trial, 'stepEvents'));
    if emg_present
        load(makeFileName(date, subject_id, 'walking', i_trial, 'emgTrajectories'));
    end
    load(makeFileName(date, subject_id, 'walking', i_trial, 'stepEvents'));

    if strcmp(stimulus_type, 'gvs')
        stim_sent_trajectory = gvs_out_trajectory;
    elseif strcmp(stimulus_type, 'visual')
        stim_sent_trajectory = visual_shift_ml_trajectory;
    end

    left_copx_trajectory = left_force_plate_cop_Acw(:, 1);
    right_copx_trajectory = right_force_plate_cop_Acw(:, 1);
    left_fz_trajectory = left_force_plate_wrench_Acw(:, 3);
    right_fz_trajectory = right_force_plate_wrench_Acw(:, 3);

    % find trigger
    % -- find the first index where the heel marker crosses the threshold
    % and is then recognized as heel strike
    trigger_indices_forceplate = find(diff(sign(stimulus_foot_state - 0.5)) > 0) + 1;
    epsilon = 1e-5;
    stim_start_indices_forceplate = find(diff(sign(abs(stim_sent_trajectory) - epsilon)) > 0) + 1;
    trigger_indices_forceplate = trigger_indices_forceplate(1 : length(stim_start_indices_forceplate)); % in case a stim is triggered, but not recorded

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
        plot(time_force_plate(trigger_indices_forceplate), zeros(size(trigger_indices_forceplate)), 'x')
        plot(time_force_plate(stim_start_indices_forceplate), zeros(size(stim_start_indices_forceplate)), 'x')
        legend('stimulus state', 'left cop', 'right cop', 'left touchdown', 'right touchdown', 'trigger', 'stim start')
    end

    % for each trigger, extract conditions and relevant step events
    number_of_triggers = length(trigger_indices_forceplate);
    removal_flags = zeros(number_of_triggers, 1);
    stim_start_indices_forceplate = [stim_start_indices_forceplate zeros(number_of_triggers, 1)]; %#ok<AGROW>
    stretch_start_indices_forceplate = zeros(number_of_triggers, 2);
    stretch_end_indices_forceplate = zeros(number_of_triggers, 2);
    closest_heelstrike_distance_times = zeros(number_of_triggers, 2);
    condition_stance_foot_list = cell(number_of_triggers, 2);
    condition_polarity_list = cell(number_of_triggers, 2);
    condition_delay_list = cell(number_of_triggers, 2);

    for i_trigger = 1 : number_of_triggers
        % polarity condition
        if stim_sent_trajectory(stim_start_indices_forceplate(i_trigger)) > 0
            condition_polarity_list{i_trigger, 1} = 'POSITIVE';
        elseif stim_sent_trajectory(stim_start_indices_forceplate(i_trigger)) < 0
            condition_polarity_list{i_trigger, 1} = 'NEGATIVE';
        else
            disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate(i_trigger))) ' - no stim']);
        end
        condition_polarity_list{i_trigger, 2} = 'CONTROL';

        % delay
        wait_time_stim = time_force_plate(stim_start_indices_forceplate(i_trigger)) - time_force_plate(trigger_indices_forceplate(i_trigger));
        [~, wait_condition_index] = min(abs(wait_times - wait_time_stim));
        condition_delay_list{i_trigger, 1} = wait_time_labels{wait_condition_index};
        condition_delay_list{i_trigger, 2} = 'CONTROL';

        % find out which heelstrike triggered
        [distance_to_trigger_left_indices, index_left] = min(abs(left_touchdown_indices_force_plate - trigger_indices_forceplate(i_trigger)));
        [distance_to_trigger_right_indices, index_right] = min(abs(right_touchdown_indices_force_plate - trigger_indices_forceplate(i_trigger)));
        
        % confirm that the trigger happened less than 100ms before the actual heelstrike
%         nearest_future_heelstrike_index = min([this_left_foot_heelstrike this_right_foot_heelstrike]);
%         nearest_future_heelstrike_time = time_force_plate(nearest_future_heelstrike_index);
%         duration_until_nearest_future_heelstrike(i_trigger, 1) = nearest_future_heelstrike_time - time_force_plate(trigger_indices_forceplate(i_trigger));
%         if duration_until_nearest_future_heelstrike(i_trigger, 1) > duration_until_nearest_future_heelstrike_threshold
%             removal_flags(i_trigger) = 1;
%             disp(['XXXXXXX Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate(i_trigger))) ' - trigger not in sync with heelstrike']);
%         end
        closest_heelstrike_distance_indices = min([distance_to_trigger_left_indices distance_to_trigger_right_indices]);
        closest_heelstrike_distance_time = time_mocap(closest_heelstrike_distance_indices+1);
        if closest_heelstrike_distance_time > duration_until_nearest_future_heelstrike_threshold
            removal_flags(i_trigger) = 1;
            disp(['XXXXXXX Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate(i_trigger))) ' - trigger not in sync with heelstrike']);
        end
        
        if distance_to_trigger_left_indices < distance_to_trigger_right_indices
            % triggered by left heelstrike
            trigger_foot = 'left';
            % check whether time offset is positive or negative
            if left_touchdown_indices_force_plate(index_left) - trigger_indices_forceplate(i_trigger) > 0
                closest_heelstrike_distance_times(i_trigger) = closest_heelstrike_distance_time;
            else
                closest_heelstrike_distance_times(i_trigger) = -closest_heelstrike_distance_time;
            end
            if index_left == 1 || length(left_touchdown_indices_force_plate) < index_left + 1
                % data doesn't include previous or next step
                removal_flags(i_trigger) = 1;
                last_right_foot_heelstrike = NaN;
                this_right_foot_heelstrike = NaN;
                next_right_foot_heelstrike = NaN;
                last_left_foot_heelstrike = NaN;
                this_left_foot_heelstrike = NaN;
                next_left_foot_heelstrike = NaN;
            else
                last_left_foot_heelstrike = left_touchdown_indices_force_plate(index_left-1);
                this_left_foot_heelstrike = left_touchdown_indices_force_plate(index_left);
                next_left_foot_heelstrike = left_touchdown_indices_force_plate(index_left+1);
                
                last_right_foot_heelstrike = min(right_touchdown_indices_force_plate(right_touchdown_indices_force_plate >= last_left_foot_heelstrike));
                this_right_foot_heelstrike = min(right_touchdown_indices_force_plate(right_touchdown_indices_force_plate >= this_left_foot_heelstrike));
                next_right_foot_heelstrike = min(right_touchdown_indices_force_plate(right_touchdown_indices_force_plate >= next_left_foot_heelstrike));
                
%                 % check check
%                 figure; hold on
%                 plot([last_right_foot_heelstrike this_right_foot_heelstrike next_right_foot_heelstrike], [0 0 0], 'x');
%                 plot([last_left_foot_heelstrike this_left_foot_heelstrike next_left_foot_heelstrike], [0 0 0], 'x');
                
            end            
        elseif distance_to_trigger_right_indices < distance_to_trigger_left_indices
            % triggered by right heelstrike
            trigger_foot = 'right';
            % check whether time offset is positive or negative
            if right_touchdown_indices_force_plate(index_right) - trigger_indices_forceplate(i_trigger) > 0
                closest_heelstrike_distance_times(i_trigger) = closest_heelstrike_distance_time;
            else
                closest_heelstrike_distance_times(i_trigger) = -closest_heelstrike_distance_time;
            end
            if index_right == 1 || length(right_touchdown_indices_force_plate) < index_right + 1
                % data doesn't include previous or next step
                removal_flags(i_trigger) = 1;
                last_left_foot_heelstrike = NaN;
                this_left_foot_heelstrike = NaN;
                next_left_foot_heelstrike = NaN;
                last_right_foot_heelstrike = NaN;
                this_right_foot_heelstrike = NaN;
                next_right_foot_heelstrike = NaN;
            else
                last_right_foot_heelstrike = right_touchdown_indices_force_plate(index_right-1);
                this_right_foot_heelstrike = right_touchdown_indices_force_plate(index_right);
                next_right_foot_heelstrike = right_touchdown_indices_force_plate(index_right+1);
                
                last_left_foot_heelstrike = min(left_touchdown_indices_force_plate(left_touchdown_indices_force_plate >= last_right_foot_heelstrike));
                this_left_foot_heelstrike = min(left_touchdown_indices_force_plate(left_touchdown_indices_force_plate >= this_right_foot_heelstrike));
                next_left_foot_heelstrike = min(left_touchdown_indices_force_plate(left_touchdown_indices_force_plate >= next_right_foot_heelstrike));
                
%                 % check check
%                 figure; hold on
%                 plot([last_right_foot_heelstrike this_right_foot_heelstrike next_right_foot_heelstrike], [0 0 0], 'x');
%                 plot([last_left_foot_heelstrike this_left_foot_heelstrike next_left_foot_heelstrike], [0 0 0], 'x');
                
            end            
        else
            trigger_foot = 'unclear';
            disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate(i_trigger))) ' - trigger exactly between two heelstrikes']);
            removal_flags(i_trigger) = 1;
        end
        
        
%         % stance foot for time period of interest
%         [last_left_foot_heelstrike, index_left] = max(left_touchdown_indices_force_plate(left_touchdown_indices_force_plate <= trigger_indices_forceplate(i_trigger)));
%         if isempty(index_left)
%             removal_flags(i_trigger) = 1;
%         else
%             if length(left_touchdown_indices_force_plate) < index_left + 2
%                 removal_flags(i_trigger) = 1;
%             else
%                 this_left_foot_heelstrike = left_touchdown_indices_force_plate(index_left+1);
%                 next_left_foot_heelstrike = left_touchdown_indices_force_plate(index_left+2);
%             end
%         end
%         [last_right_foot_heelstrike, index_right] = max(right_touchdown_indices_force_plate(right_touchdown_indices_force_plate <= trigger_indices_forceplate(i_trigger)));
%         if isempty(index_right)
%             removal_flags(i_trigger) = 1;
%         else
%             if length(right_touchdown_indices_force_plate) < index_right + 2
%                 removal_flags(i_trigger) = 1;
%             else
%                 this_right_foot_heelstrike = right_touchdown_indices_force_plate(index_right+1);
%                 next_right_foot_heelstrike = right_touchdown_indices_force_plate(index_right+2);
%             end
%         end


        % mark relevant event delimiters depending on wait time and triggering foot
%         if left_contact_indicators_force_plate(trigger_indices_forceplate(i_trigger)) == 1 && right_contact_indicators_force_plate(trigger_indices_forceplate(i_trigger)) == 0
        if strcmp(trigger_foot, 'right')
            % triggered by right foot heelstrike
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
                disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate(i_trigger))) ' - delay not defined']);
            end
%         elseif left_contact_indicators_force_plate(trigger_indices_forceplate(i_trigger)) == 0 && right_contact_indicators_force_plate(trigger_indices_forceplate(i_trigger)) == 1
        elseif strcmp(trigger_foot, 'left')
            % triggered by left foot heelstrike
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
                disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate(i_trigger))) ' - delay not defined']);
            end
        else
            removal_flags(i_trigger) = 1;
            condition_stance_foot_list{i_trigger, 1} = 'UNCLEAR';
            condition_stance_foot_list{i_trigger, 2} = 'UNCLEAR';
%             disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate(i_trigger))) ' - stance foot unclear']);
        end

    end

    % remove flagged triggers
    unflagged_indices = ~removal_flags;
    stim_start_indices_forceplate = stim_start_indices_forceplate(unflagged_indices, :);
    stretch_start_indices_forceplate = stretch_start_indices_forceplate(unflagged_indices, :);
    stretch_end_indices_forceplate = stretch_end_indices_forceplate(unflagged_indices, :);
    condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
    condition_polarity_list = condition_polarity_list(unflagged_indices, :);
    condition_delay_list = condition_delay_list(unflagged_indices, :);
    closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :);

    % form stretches
    stim_start_indices_forceplate = [stim_start_indices_forceplate(:, 1); stim_start_indices_forceplate(:, 2)];
    stretch_start_indices_forceplate = [stretch_start_indices_forceplate(:, 1); stretch_start_indices_forceplate(:, 2)];
    stretch_end_indices_forceplate = [stretch_end_indices_forceplate(:, 1); stretch_end_indices_forceplate(:, 2)];
    condition_stance_foot_list = [condition_stance_foot_list(:, 1); condition_stance_foot_list(:, 2)];
    condition_polarity_list = [condition_polarity_list(:, 1); condition_polarity_list(:, 2)];
    condition_delay_list = [condition_delay_list(:, 1); condition_delay_list(:, 2)];

    % extract and time-normalize data
    number_of_stretches = length(stretch_start_indices_forceplate);

    stim_start_time_relative_to_stretch = zeros(number_of_stretches, 1);
    number_of_indices_per_step = zeros(number_of_stretches, 1);
    start_indices_mocap = zeros(number_of_stretches, 1);
    end_indices_mocap = zeros(number_of_stretches, 1);
    start_indices_forceplate = zeros(number_of_stretches, 1);
    end_indices_forceplate = zeros(number_of_stretches, 1);
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
        elseif strcmp(condition_stance_foot_list{i_stretch}, 'LEFT')
            stance_foot_cop_stretch = left_copx_trajectory(stretch_start_index_forceplate : stretch_end_index_forceplate);
            touchdown_index_forceplate = stretch_start_index_forceplate + find(stance_foot_cop_stretch~=0, 1, 'first') - 1;
            stance_foot_fz_step = left_fz_trajectory(touchdown_index_forceplate : stretch_end_index_forceplate);
            swing_foot_fz_step = right_fz_trajectory(touchdown_index_forceplate : stretch_end_index_forceplate);
        end
        touchdown_time_forceplate = time_force_plate(touchdown_index_forceplate);
        time_extracted_forceplate = time_force_plate(touchdown_index_forceplate : stretch_end_index_forceplate);

        % check if this stretch is usable
        number_of_indices_per_step(i_stretch) = length(time_extracted_forceplate);
        number_of_swing_foot_fz_zeros_stretch = sum(-swing_foot_fz_step < swing_foot_fz_zero_threshold);
        if number_of_swing_foot_fz_zeros_stretch < swing_foot_zero_stretch_length_threshold
            removal_flags(i_stretch) = 1;
            disp(['excluding stretch index ' num2str(i_stretch) ' - cross over at time ' num2str(time_force_plate(stretch_start_index_forceplate))]);
        end
        if isempty(touchdown_index_forceplate)
            touchdown_index_forceplate = stretch_start_index_forceplate;
            touchdown_time_forceplate = time_force_plate(touchdown_index_forceplate);
            removal_flags(i_stretch) = 1;
        end

        % time
        start_index_force_plate = touchdown_index_forceplate;
        end_index_force_plate = stretch_end_indices_forceplate(i_stretch);
        [~, start_index_mocap] = min(abs(time_mocap - time_force_plate(start_index_force_plate)));
        [~, end_index_mocap] = min(abs(time_mocap - time_force_plate(end_index_force_plate)));
        if emg_present
            [~, start_index_emg] = min(abs(time_emg - time_force_plate(start_index_force_plate)));
            [~, end_index_emg] = min(abs(time_emg - time_force_plate(end_index_force_plate)));
        end
        
        % store
        start_indices_mocap(i_stretch) = start_index_mocap;
        end_indices_mocap(i_stretch) = end_index_mocap;
        start_indices_forceplate(i_stretch) = start_index_force_plate;
        end_indices_forceplate(i_stretch) = end_index_force_plate;
        if emg_present
            start_indices_emg(i_stretch) = start_index_emg;
            end_indices_emg(i_stretch) = end_index_emg;
        end
        
        if stim_start_indices_forceplate(i_stretch) == 0
            stim_start_time_relative_to_stretch(i_stretch) = NaN;
        else
            stim_start_time = time_force_plate(stim_start_indices_forceplate(i_stretch));
            stim_start_time_relative_to_stretch(i_stretch) = stim_start_time - touchdown_time_forceplate;
        end

        % extract force plate data
        left_cop_x_extracted_stretch = left_copx_trajectory(start_index_force_plate : end_index_force_plate);
        right_cop_x_extracted_stretch = right_copx_trajectory(start_index_force_plate : end_index_force_plate);

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
        'start_indices_emg', ...
        'end_indices_emg', ...
        'stim_start_time_relative_to_stretch', ...
        'data_points_to_process_mocap' ...
      )

    disp(['Trial ' num2str(i_trial) ' completed']);

end

























