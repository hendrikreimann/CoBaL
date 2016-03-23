
visualize_steps                 = 0;
do_cop_plots                    = 1;
do_heel_plots                   = 1;



wait_times = [0 0.150 0.450];
wait_time_labels = {'0ms', '150ms', '450ms'};
load subjectInfo.mat;

trials_to_process_control = 1;
trials_to_process_stim = 6 : 14;
trials_to_process_stim = 6;


% trials_to_process = [trials_to_process_control trials_to_process_stim];
trials_to_process = trials_to_process_stim;
number_of_time_steps_normalized = 64;

left_heel_marker = 34;
left_toes_marker = 35;
right_heel_marker = 42;
right_toes_marker = 43;

left_heel_marker_indices = reshape([(left_heel_marker - 1) * 3 + 1; (left_heel_marker - 1) * 3 + 2; (left_heel_marker - 1) * 3 + 3], 1, length(left_heel_marker)*3);
left_toes_marker_indices = reshape([(left_toes_marker - 1) * 3 + 1; (left_toes_marker - 1) * 3 + 2; (left_toes_marker - 1) * 3 + 3], 1, length(left_toes_marker)*3);
right_heel_marker_indices = reshape([(right_heel_marker - 1) * 3 + 1; (right_heel_marker - 1) * 3 + 2; (right_heel_marker - 1) * 3 + 3], 1, length(right_heel_marker)*3);
right_toes_marker_indices = reshape([(right_toes_marker - 1) * 3 + 1; (right_toes_marker - 1) * 3 + 2; (right_toes_marker - 1) * 3 + 3], 1, length(right_toes_marker)*3);


condition_stance_foot_list = {};
condition_polarity_list = {};
condition_delay_list = {};
left_cop_x_normalized = [];
right_cop_x_normalized = [];
left_heel_marker_x_normalized = [];
right_heel_marker_x_normalized = [];
step_times = [];
for i_trial = trials_to_process
% for i_trial = 2
    % load data
    load(makeFileName(date, subject_id, 'walking', i_trial, 'forcePlateData'));
    load(makeFileName(date, subject_id, 'walking', i_trial, 'markerTrajectories'));
    load(makeFileName(date, subject_id, 'walking', i_trial, 'stepEvents'));
    
    left_heel_marker_x_trajectory = marker_trajectories(:, left_heel_marker_indices(1));
    left_toes_marker_x_trajectory = marker_trajectories(:, left_toes_marker_indices(1));
    right_heel_marker_x_trajectory = marker_trajectories(:, right_heel_marker_indices(1));
    right_toes_marker_x_trajectory = marker_trajectories(:, right_toes_marker_indices(1));
    left_cop_x_trajectory = left_force_plate_cop_Acw(:, 1);
    right_cop_x_trajectory = right_force_plate_cop_Acw(:, 1);
    
    % find trigger
    if ismember(i_trial, trials_to_process_control)
        trigger_indices_forceplate_trial = find(diff(abs(heel_strike_count)) ~= 0);
        stim_start_indices_forceplate_trial = find(diff(abs(heel_strike_count)) ~= 0) + 1;
    elseif ismember(i_trial, trials_to_process_stim)
        trigger_indices_forceplate_trial = find(diff(sign(stimulus_foot_state - 0.5)) > 0) + 1;
        stim_start_indices_forceplate_trial = find(diff(abs(stim_sent_trajectory)) > 0) + 1;
    else
        disp(['Trial ' num2str(i_trial) ': something went wrong, this trial is neither stim nor control']);
    end
    
    % for each trigger, extract conditions and relevant step events
    number_of_triggers = length(trigger_indices_forceplate_trial);
    stim_step_start_indices_forceplate_trial = zeros(size(trigger_indices_forceplate_trial));
    stim_step_end_indices_forceplate_trial = zeros(size(trigger_indices_forceplate_trial));
    control_step_start_indices_forceplate_trial = zeros(size(trigger_indices_forceplate_trial));
    control_step_end_indices_forceplate_trial = zeros(size(trigger_indices_forceplate_trial));
    condition_stance_foot_list_trial = cell(number_of_triggers, 2);
    condition_polarity_list_trial = cell(number_of_triggers, 2);
    condition_delay_list_trial = cell(number_of_triggers, 2);
    for i_trigger = 1 : number_of_triggers
        % polarity
        if ismember(i_trial, trials_to_process_control)
            condition_polarity_list_trial{i_trigger, 1} = 'CONTROL';
        else
            if stim_sent_trajectory(stim_start_indices_forceplate_trial(i_trigger)) > 0
                condition_polarity_list_trial{i_trigger, 1} = 'POSITIVE';
            elseif stim_sent_trajectory(stim_start_indices_forceplate_trial(i_trigger)) < 0
                condition_polarity_list_trial{i_trigger, 1} = 'NEGATIVE';
            else
                disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate_trial(i_trigger))) ' - no stim']);
            end
            condition_polarity_list_trial{i_trigger, 2} = 'CONTROL';
        end
        % delay
        wait_time_stim = time_force_plate(stim_start_indices_forceplate_trial(i_trigger)) - time_force_plate(trigger_indices_forceplate_trial(i_trigger));
        [~, wait_condition_index] = min(abs(wait_times - wait_time_stim));
        condition_delay_list_trial{i_trigger, 1} = wait_time_labels{wait_condition_index};
        condition_delay_list_trial{i_trigger, 2} = 'CONTROL';
        
        % stance foot
        [last_left_foot_heelstrike, index_left] = max(left_touchdown_indices_force_plate(left_touchdown_indices_force_plate <= trigger_indices_forceplate_trial(i_trigger)));
        [last_right_foot_heelstrike, index_right] = max(right_touchdown_indices_force_plate(right_touchdown_indices_force_plate <= trigger_indices_forceplate_trial(i_trigger)));
        this_left_foot_heelstrike = left_touchdown_indices_force_plate(index_left+1);
        this_right_foot_heelstrike = right_touchdown_indices_force_plate(index_right+1);
        next_left_foot_heelstrike = left_touchdown_indices_force_plate(index_left+2);
        next_right_foot_heelstrike = right_touchdown_indices_force_plate(index_right+2);
        
        if left_contact_indicators_force_plate(trigger_indices_forceplate_trial(i_trigger)) == 1 && right_contact_indicators_force_plate(trigger_indices_forceplate_trial(i_trigger)) == 0
            % triggered by right foot heelstrike
            if strcmp(condition_delay_list_trial{i_trigger, 1}, '0ms') || strcmp(condition_delay_list_trial{i_trigger, 1}, '150ms')
                % this step is of interest
                condition_stance_foot_list_trial{i_trigger, 1} = 'RIGHT';
                condition_stance_foot_list_trial{i_trigger, 2} = 'RIGHT';
                stim_step_start_indices_forceplate_trial(i_trigger) = this_right_foot_heelstrike;
                stim_step_end_indices_forceplate_trial(i_trigger) = this_left_foot_heelstrike;
                control_step_start_indices_forceplate_trial(i_trigger) = last_right_foot_heelstrike;
                control_step_end_indices_forceplate_trial(i_trigger) = last_left_foot_heelstrike;
            elseif strcmp(condition_delay_list_trial{i_trigger, 1}, '450ms')
                % next step is of interest
                condition_stance_foot_list_trial{i_trigger, 1} = 'LEFT';
                condition_stance_foot_list_trial{i_trigger, 2} = 'LEFT';
                stim_step_start_indices_forceplate_trial(i_trigger) = next_left_foot_heelstrike;
                stim_step_end_indices_forceplate_trial(i_trigger) = next_right_foot_heelstrike;
                control_step_start_indices_forceplate_trial(i_trigger) = last_left_foot_heelstrike;
                control_step_end_indices_forceplate_trial(i_trigger) = this_right_foot_heelstrike;
            else
                disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate_trial(i_trigger))) ' - delay not defined']);
            end
        elseif left_contact_indicators_force_plate(trigger_indices_forceplate_trial(i_trigger)) == 0 && right_contact_indicators_force_plate(trigger_indices_forceplate_trial(i_trigger)) == 1
            % triggered by left foot heelstrike
            if strcmp(condition_delay_list_trial{i_trigger, 1}, '0ms') || strcmp(condition_delay_list_trial{i_trigger, 1}, '150ms')
                % this step is of interest
                condition_stance_foot_list_trial{i_trigger, 1} = 'LEFT';
                condition_stance_foot_list_trial{i_trigger, 2} = 'LEFT';
                stim_step_start_indices_forceplate_trial(i_trigger) = this_left_foot_heelstrike;
                stim_step_end_indices_forceplate_trial(i_trigger) = this_right_foot_heelstrike;
                control_step_start_indices_forceplate_trial(i_trigger) = last_left_foot_heelstrike;
                control_step_end_indices_forceplate_trial(i_trigger) = last_right_foot_heelstrike;
            elseif strcmp(condition_delay_list_trial{i_trigger, 1}, '450ms')
                % next step is of interest
                condition_stance_foot_list_trial{i_trigger, 1} = 'RIGHT';
                condition_stance_foot_list_trial{i_trigger, 2} = 'RIGHT';
                stim_step_start_indices_forceplate_trial(i_trigger) = next_right_foot_heelstrike;
                stim_step_end_indices_forceplate_trial(i_trigger) = next_left_foot_heelstrike;
                control_step_start_indices_forceplate_trial(i_trigger) = last_right_foot_heelstrike;
                control_step_end_indices_forceplate_trial(i_trigger) = this_left_foot_heelstrike;
            else
                disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate_trial(i_trigger))) ' - delay not defined']);
            end
            condition_stance_foot_list_trial{i_trigger, 1} = 'UNCLEAR';
            condition_stance_foot_list_trial{i_trigger, 2} = 'UNCLEAR';
            disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate_trial(i_trigger))) ' - stance foot unclear']);
        end
        
        
        
        % stance foot for time period of interest
%         if left_contact_indicators_force_plate(trigger_indices_forceplate_trial(i_trigger)) == 1 && right_contact_indicators_force_plate(trigger_indices_forceplate_trial(i_trigger)) == 0
%             % right foot was in swing, so right heelstrike triggered the stim
%             condition_stance_foot_list_trial{i_trigger, 1} = 'RIGHT';
%             condition_stance_foot_list_trial{i_trigger, 2} = 'RIGHT';
%             this_right_foot_heelstrike = min(right_touchdown_indices_force_plate(right_touchdown_indices_force_plate > trigger_indices_forceplate_trial(i_trigger)));
%             this_left_foot_heelstrike = min(left_touchdown_indices_force_plate(left_touchdown_indices_force_plate > trigger_indices_forceplate_trial(i_trigger)));
%             last_right_foot_heelstrike = max(right_touchdown_indices_force_plate(right_touchdown_indices_force_plate <= trigger_indices_forceplate_trial(i_trigger)));
%             last_left_foot_heelstrike = max(left_touchdown_indices_force_plate(left_touchdown_indices_force_plate <= trigger_indices_forceplate_trial(i_trigger)));
%             if ~(isempty(this_right_foot_heelstrike) || isempty(this_left_foot_heelstrike))
%                 stim_step_start_indices_forceplate_trial(i_trigger) = this_right_foot_heelstrike;
%                 stim_step_end_indices_forceplate_trial(i_trigger) = this_left_foot_heelstrike;
%                 control_step_start_indices_forceplate_trial(i_trigger) = last_right_foot_heelstrike;
%                 control_step_end_indices_forceplate_trial(i_trigger) = last_left_foot_heelstrike;
%             else
%                 % data not complete, flag for removal
%                 stim_step_start_indices_forceplate_trial(i_trigger) = NaN;
%                 stim_step_end_indices_forceplate_trial(i_trigger) = NaN;
%                 control_step_start_indices_forceplate_trial(i_trigger) = NaN;
%                 control_step_end_indices_forceplate_trial(i_trigger) = NaN;
%             end
%         elseif left_contact_indicators_force_plate(trigger_indices_forceplate_trial(i_trigger)) == 0 && right_contact_indicators_force_plate(trigger_indices_forceplate_trial(i_trigger)) == 1
%             condition_stance_foot_list_trial{i_trigger, 1} = 'LEFT';
%             condition_stance_foot_list_trial{i_trigger, 2} = 'LEFT';
%             this_left_foot_heelstrike = min(left_touchdown_indices_force_plate(left_touchdown_indices_force_plate > trigger_indices_forceplate_trial(i_trigger)));
%             this_right_foot_heelstrike = min(right_touchdown_indices_force_plate(right_touchdown_indices_force_plate > trigger_indices_forceplate_trial(i_trigger)));
%             last_left_foot_heelstrike = max(left_touchdown_indices_force_plate(left_touchdown_indices_force_plate <= trigger_indices_forceplate_trial(i_trigger)));
%             last_right_foot_heelstrike = max(right_touchdown_indices_force_plate(right_touchdown_indices_force_plate <= trigger_indices_forceplate_trial(i_trigger)));
%             if ~(isempty(this_right_foot_heelstrike) || isempty(this_left_foot_heelstrike))
%                 stim_step_start_indices_forceplate_trial(i_trigger) = this_left_foot_heelstrike;
%                 stim_step_end_indices_forceplate_trial(i_trigger) = this_right_foot_heelstrike;
%                 control_step_start_indices_forceplate_trial(i_trigger) = last_left_foot_heelstrike;
%                 control_step_end_indices_forceplate_trial(i_trigger) = last_right_foot_heelstrike;
%             else
%                 % data not complete, flag for removal
%                 stim_step_start_indices_forceplate_trial(i_trigger) = NaN;
%                 stim_step_end_indices_forceplate_trial(i_trigger) = NaN;
%                 control_step_start_indices_forceplate_trial(i_trigger) = NaN;
%                 control_step_end_indices_forceplate_trial(i_trigger) = NaN;
%             end
%         else
%             condition_stance_foot_list_trial{i_trigger, 1} = 'UNCLEAR';
%             condition_stance_foot_list_trial{i_trigger, 2} = 'UNCLEAR';
%             disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_force_plate(trigger_indices_forceplate_trial(i_trigger))) ' - stance foot unclear']);
%         end
    end
    
    % remove flagged stims
    unflagged_indices = ...
      ~( ...
         isnan(stim_step_start_indices_forceplate_trial(:, 1)) ...
         | isnan(stim_step_end_indices_forceplate_trial(:, 1)) ...
         | strcmp(condition_stance_foot_list_trial(:, 1), 'UNCLEAR') ...
       );
    stim_step_start_indices_forceplate_trial = stim_step_start_indices_forceplate_trial(unflagged_indices);
    stim_step_end_indices_forceplate_trial = stim_step_end_indices_forceplate_trial(unflagged_indices);
    control_step_start_indices_forceplate_trial = control_step_start_indices_forceplate_trial(unflagged_indices);
    control_step_end_indices_forceplate_trial = control_step_end_indices_forceplate_trial(unflagged_indices);
    condition_stance_foot_list_trial = condition_stance_foot_list_trial(unflagged_indices, :);
    condition_polarity_list_trial = condition_polarity_list_trial(unflagged_indices, :);
    condition_delay_list_trial = condition_delay_list_trial(unflagged_indices, :);
    
    % form stretches
    stretch_start_indices_forceplate_trial = [stim_step_start_indices_forceplate_trial; control_step_start_indices_forceplate_trial];
    stretch_end_indices_forceplate_trial = [stim_step_end_indices_forceplate_trial; control_step_end_indices_forceplate_trial];
    condition_stance_foot_list_trial = [condition_stance_foot_list_trial(:, 1); condition_stance_foot_list_trial(:, 2)];
    condition_polarity_list_trial = [condition_polarity_list_trial(:, 1); condition_polarity_list_trial(:, 2)];
    condition_delay_list_trial = [condition_delay_list_trial(:, 1); condition_delay_list_trial(:, 2)];
    
    % extract and time-normalize data
    number_of_stretches = length(stretch_start_indices_forceplate_trial);
    step_times_trial = zeros(1, number_of_stretches);
    number_of_time_steps_extracted = zeros(1, number_of_stretches);
    left_cop_x_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches);
    right_cop_x_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches);
    left_heel_marker_x_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches);
    right_heel_marker_x_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches);
    removal_flags = zeros(number_of_stretches, 1);
    for i_stretch = 1 : number_of_stretches
        % time
        start_index_force_plate = stretch_start_indices_forceplate_trial(i_stretch) + 1;
        end_index_force_plate = stretch_end_indices_forceplate_trial(i_stretch);
        [~, start_index_mocap] = min(abs(time_mocap - time_force_plate(start_index_force_plate)));
        [~, end_index_mocap] = min(abs(time_mocap - time_force_plate(end_index_force_plate)));
        step_times_trial(i_stretch) = time_force_plate(end_index_force_plate) - time_force_plate(start_index_force_plate);
        
        % normalize force plate data
        left_cop_x_extracted_stim = left_cop_x_trajectory(start_index_force_plate : end_index_force_plate);
        right_cop_x_extracted_stim = right_cop_x_trajectory(start_index_force_plate : end_index_force_plate);
        time_extracted = time_force_plate(start_index_force_plate : end_index_force_plate);
        time_normalized = linspace(time_extracted(1), time_extracted(end), number_of_time_steps_normalized);
        left_cop_x_normalized_stim = spline(time_extracted, left_cop_x_extracted_stim, time_normalized);
        right_cop_x_normalized_stim = spline(time_extracted, right_cop_x_extracted_stim, time_normalized);
        
        left_cop_x_normalized_stim = left_cop_x_normalized_stim - left_cop_x_normalized_stim(1);
        right_cop_x_normalized_stim = right_cop_x_normalized_stim - right_cop_x_normalized_stim(1);
        
        left_cop_x_normalized_trial(:, i_stretch) = left_cop_x_normalized_stim;
        right_cop_x_normalized_trial(:, i_stretch) = right_cop_x_normalized_stim;
        
        % check for data - one force plate data stretch should have no zeros
        if any(left_cop_x_extracted_stim==0) & any(right_cop_x_extracted_stim==0)
            removal_flags(i_stretch) = 1;
        end
        
        % mocap data
        left_heel_marker_x_extracted_stim = left_heel_marker_x_trajectory(start_index_mocap : end_index_mocap);
        right_heel_marker_x_extracted_stim = right_heel_marker_x_trajectory(start_index_mocap : end_index_mocap);
        time_extracted = time_mocap(start_index_mocap : end_index_mocap);
        time_normalized = linspace(time_extracted(1), time_extracted(end), number_of_time_steps_normalized);
        left_heel_marker_x_normalized_stim = spline(time_extracted, left_heel_marker_x_extracted_stim, time_normalized);
        right_heel_marker_x_normalized_stim = spline(time_extracted, right_heel_marker_x_extracted_stim, time_normalized);
                
%         left_heel_marker_x_normalized_stim = left_heel_marker_x_normalized_stim - left_heel_marker_x_normalized_stim(1);
%         right_heel_marker_x_normalized_stim = right_heel_marker_x_normalized_stim - right_heel_marker_x_normalized_stim(1);
        
        left_heel_marker_x_normalized_trial(:, i_stretch) = left_heel_marker_x_normalized_stim;
        right_heel_marker_x_normalized_trial(:, i_stretch) = right_heel_marker_x_normalized_stim;
    end
    
    % remove flagged stretches
    unflagged_indices = ~removal_flags;
    stretch_start_indices_forceplate_trial = stretch_start_indices_forceplate_trial(unflagged_indices);
    stretch_end_indices_forceplate_trial = stretch_end_indices_forceplate_trial(unflagged_indices);
    condition_stance_foot_list_trial = condition_stance_foot_list_trial(unflagged_indices);
    condition_polarity_list_trial = condition_polarity_list_trial(unflagged_indices);
    left_cop_x_normalized_trial = left_cop_x_normalized_trial(:, unflagged_indices);
    right_cop_x_normalized_trial = right_cop_x_normalized_trial(:, unflagged_indices);
    left_heel_marker_x_normalized_trial = left_heel_marker_x_normalized_trial(:, unflagged_indices);
    right_heel_marker_x_normalized_trial = right_heel_marker_x_normalized_trial(:, unflagged_indices);
    
    
    % add data to lists
    condition_stance_foot_list = [condition_stance_foot_list; condition_stance_foot_list_trial];
    condition_polarity_list = [condition_polarity_list; condition_polarity_list_trial];
    left_cop_x_normalized = [left_cop_x_normalized left_cop_x_normalized_trial];
    right_cop_x_normalized = [right_cop_x_normalized right_cop_x_normalized_trial];
    left_heel_marker_x_normalized = [left_heel_marker_x_normalized left_heel_marker_x_normalized_trial];
    right_heel_marker_x_normalized = [right_heel_marker_x_normalized right_heel_marker_x_normalized_trial];
    step_times = [step_times step_times_trial];
    
    
    
    
    
    
    %% visualize
    if visualize_steps
        figure; check_axes = axes; hold on
        plot(time_force_plate, heel_strike_count);
        plot(time_force_plate(trigger_indices_forceplate_trial), heel_strike_count(trigger_indices_forceplate_trial), 'o', 'linewidth', 1, 'markersize', 10);

        % left events
        figure; axes_left = axes; hold on; title('Left'); xlim([0, recordingTime]);
        plot(time_force_plate, fzl_trajectory)
        plot(time_mocap(left_touchdown_indices_mocap), zeros(size(left_touchdown_indices_mocap)), 'v', 'linewidth', 1, 'markersize', 10);
        plot(time_mocap(left_pushoff_indices_mocap), zeros(size(left_pushoff_indices_mocap)), '^', 'linewidth', 1, 'markersize', 10);
        plot(time_force_plate(trigger_indices_forceplate_trial), fzl_trajectory(trigger_indices_forceplate_trial), 'o', 'linewidth', 1, 'markersize', 10);
        plot(time_force_plate(stim_step_start_indices_forceplate_trial), fzl_trajectory(stim_step_start_indices_forceplate_trial), '<', 'linewidth', 1, 'markersize', 10);
        plot(time_force_plate(stim_step_end_indices_forceplate_trial), fzl_trajectory(stim_step_end_indices_forceplate_trial), '>', 'linewidth', 1, 'markersize', 10);

        % right events
        figure; axes_right = axes; hold on; title('Right'); xlim([0, recordingTime]);
        plot(time_force_plate, fzr_trajectory)
        plot(time_mocap(right_touchdown_indices_mocap), zeros(size(right_touchdown_indices_mocap)), 'v', 'linewidth', 1, 'markersize', 10);
        plot(time_mocap(right_pushoff_indices_mocap), zeros(size(right_pushoff_indices_mocap)), '^', 'linewidth', 1, 'markersize', 10);
        plot(time_force_plate(trigger_indices_forceplate_trial), fzr_trajectory(trigger_indices_forceplate_trial), 'o', 'linewidth', 1, 'markersize', 10);
        plot(time_force_plate(stim_step_start_indices_forceplate_trial), fzr_trajectory(stim_step_start_indices_forceplate_trial), '<', 'linewidth', 1, 'markersize', 10);
        plot(time_force_plate(stim_step_end_indices_forceplate_trial), fzr_trajectory(stim_step_end_indices_forceplate_trial), '>', 'linewidth', 1, 'markersize', 10);

        linkaxes([axes_left, axes_right check_axes], 'x')
        distFig('rows', 3)
    end
    
    
end

%% extract conditions
step_time_mean = mean(step_times);
time_normalized = linspace(0, step_time_mean, number_of_time_steps_normalized);
conditions_left_control = strcmp(condition_stance_foot_list, 'LEFT') & strcmp(condition_polarity_list, 'CONTROL');
conditions_left_positive = strcmp(condition_stance_foot_list, 'LEFT') & strcmp(condition_polarity_list, 'POSITIVE');
conditions_left_negative = strcmp(condition_stance_foot_list, 'LEFT') & strcmp(condition_polarity_list, 'NEGATIVE');
conditions_right_control = strcmp(condition_stance_foot_list, 'RIGHT') & strcmp(condition_polarity_list, 'CONTROL');
conditions_right_positive = strcmp(condition_stance_foot_list, 'RIGHT') & strcmp(condition_polarity_list, 'POSITIVE');
conditions_right_negative = strcmp(condition_stance_foot_list, 'RIGHT') & strcmp(condition_polarity_list, 'NEGATIVE');

%% calculate means and confidence intervals
left_cop_x_mean_left_control = mean(left_cop_x_normalized(:, conditions_left_control), 2);
left_cop_x_civ_left_control = tinv(0.975, sum(conditions_left_control)-1) * std(left_cop_x_normalized(:, conditions_left_control), 1, 2)/sqrt(sum(conditions_left_control));                      % Confidence Intervals    
left_cop_x_mean_left_positive = mean(left_cop_x_normalized(:, conditions_left_positive), 2);
left_cop_x_civ_left_positive = tinv(0.975, sum(conditions_left_positive)-1) * std(left_cop_x_normalized(:, conditions_left_positive), 1, 2)/sqrt(sum(conditions_left_positive));                      % Confidence Intervals    
left_cop_x_mean_left_negative = mean(left_cop_x_normalized(:, conditions_left_negative), 2);
left_cop_x_civ_left_negative = tinv(0.975, sum(conditions_left_negative)-1) * std(left_cop_x_normalized(:, conditions_left_negative), 1, 2)/sqrt(sum(conditions_left_negative));                      % Confidence Intervals    

right_cop_x_mean_right_control = mean(right_cop_x_normalized(:, conditions_right_control), 2);
right_cop_x_civ_right_control = tinv(0.975, sum(conditions_right_control)-1) * std(right_cop_x_normalized(:, conditions_right_control), 1, 2)/sqrt(sum(conditions_right_control));                      % Confidence Intervals    
right_cop_x_mean_right_positive = mean(right_cop_x_normalized(:, conditions_right_positive), 2);
right_cop_x_civ_right_positive = tinv(0.975, sum(conditions_right_positive)-1) * std(right_cop_x_normalized(:, conditions_right_positive), 1, 2)/sqrt(sum(conditions_right_positive));                      % Confidence Intervals    
right_cop_x_mean_right_negative = mean(right_cop_x_normalized(:, conditions_right_negative), 2);
right_cop_x_civ_right_negative = tinv(0.975, sum(conditions_right_negative)-1) * std(right_cop_x_normalized(:, conditions_right_negative), 1, 2)/sqrt(sum(conditions_right_negative));                      % Confidence Intervals    

left_heel_marker_x_mean_right_control = mean(left_heel_marker_x_normalized(:, conditions_right_control), 2);
left_heel_marker_x_civ_right_control = tinv(0.975, sum(conditions_right_control)-1) * std(left_heel_marker_x_normalized(:, conditions_right_control), 1, 2)/sqrt(sum(conditions_right_control));                      % Confidence Intervals    
left_heel_marker_x_mean_right_positive = mean(left_heel_marker_x_normalized(:, conditions_right_positive), 2);
left_heel_marker_x_civ_right_positive = tinv(0.975, sum(conditions_right_positive)-1) * std(left_heel_marker_x_normalized(:, conditions_right_positive), 1, 2)/sqrt(sum(conditions_right_positive));                      % Confidence Intervals    
left_heel_marker_x_mean_right_negative = mean(left_heel_marker_x_normalized(:, conditions_right_negative), 2);
left_heel_marker_x_civ_right_negative = tinv(0.975, sum(conditions_right_negative)-1) * std(left_heel_marker_x_normalized(:, conditions_right_negative), 1, 2)/sqrt(sum(conditions_right_negative));                      % Confidence Intervals    

right_heel_marker_x_mean_left_control = mean(right_heel_marker_x_normalized(:, conditions_left_control), 2);
right_heel_marker_x_civ_left_control = tinv(0.975, sum(conditions_left_control)-1) * std(right_heel_marker_x_normalized(:, conditions_left_control), 1, 2)/sqrt(sum(conditions_left_control));                      % Confidence Intervals    
right_heel_marker_x_mean_left_positive = mean(right_heel_marker_x_normalized(:, conditions_left_positive), 2);
right_heel_marker_x_civ_left_positive = tinv(0.975, sum(conditions_left_positive)-1) * std(right_heel_marker_x_normalized(:, conditions_left_positive), 1, 2)/sqrt(sum(conditions_left_positive));                      % Confidence Intervals    
right_heel_marker_x_mean_left_negative = mean(right_heel_marker_x_normalized(:, conditions_left_negative), 2);
right_heel_marker_x_civ_left_negative = tinv(0.975, sum(conditions_left_negative)-1) * std(right_heel_marker_x_normalized(:, conditions_left_negative), 1, 2)/sqrt(sum(conditions_left_negative));                      % Confidence Intervals    


%% plots
color_left_control = [0.3 0.1 1];
color_left_positive = [1 0.3 .1] * 0.7;
color_left_negative = [0.3 1 0.1] * 0.7;
color_right_control = [0.3 0.1 1];
color_right_positive = [1 0.3 .1] * 0.7;
color_right_negative = [0.3 1 0.1] * 0.7;

%% do_cop_plots
if do_cop_plots
    % figure; axes; hold on
    % plot(time_normalized, left_cop_x_normalized(:, conditions_left_control))

    % figure; axes; hold on
    % plot(time_normalized, right_cop_x_normalized(:, conditions_right_control))

    % return
    %
    figure; axes; hold on; title('left foot medial-lateral CoP')
    shadedErrorBar(time_normalized, left_cop_x_mean_left_control, left_cop_x_civ_left_control, {'color', color_left_control, 'linewidth', 5}, 1);
    shadedErrorBar(time_normalized, left_cop_x_mean_left_positive, left_cop_x_civ_left_positive, {'color', color_left_positive, 'linewidth', 5}, 1);
    shadedErrorBar(time_normalized, left_cop_x_mean_left_negative, left_cop_x_civ_left_negative, {'color', color_left_negative, 'linewidth', 5}, 1);
    xlabel('time')

    figure; axes; hold on; title('right foot medial-lateral CoP')
    shadedErrorBar(time_normalized, right_cop_x_mean_right_control, right_cop_x_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
    shadedErrorBar(time_normalized, right_cop_x_mean_right_positive, right_cop_x_civ_right_positive, {'color', color_right_positive, 'linewidth', 5}, 1);
    shadedErrorBar(time_normalized, right_cop_x_mean_right_negative, right_cop_x_civ_right_negative, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time')
end

%% do_heel_plots
if do_heel_plots
    
%     figure; axes; hold on; title('left foot medial-lateral heel marker')
%     plot(time_normalized, left_heel_marker_x_normalized(:, conditions_right_control), 'color', color_right_control);
%     plot(time_normalized, left_heel_marker_x_normalized(:, conditions_right_positive), 'color', color_right_positive);
%     plot(time_normalized, left_heel_marker_x_normalized(:, conditions_right_negative), 'color', color_right_negative);
%     xlabel('time')
% 
%     figure; axes; hold on; title('right foot medial-lateral heel marker')
%     plot(time_normalized, right_heel_marker_x_normalized(:, conditions_left_control), 'color', color_left_control);
%     plot(time_normalized, right_heel_marker_x_normalized(:, conditions_left_positive), 'color', color_left_positive);
%     plot(time_normalized, right_heel_marker_x_normalized(:, conditions_left_negative), 'color', color_left_negative);
%     xlabel('time')
    
    % shaded error bars
    figure; axes; hold on; title('left foot medial-lateral heel marker')
    shadedErrorBar(time_normalized, left_heel_marker_x_mean_right_control, left_heel_marker_x_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
    shadedErrorBar(time_normalized, left_heel_marker_x_mean_right_positive, left_heel_marker_x_civ_right_positive, {'color', color_right_positive, 'linewidth', 5}, 1);
    shadedErrorBar(time_normalized, left_heel_marker_x_mean_right_negative, left_heel_marker_x_civ_right_negative, {'color', color_right_negative, 'linewidth', 5}, 1);
    xlabel('time')

    figure; axes; hold on; title('right foot medial-lateral heel marker')
    shadedErrorBar(time_normalized, right_heel_marker_x_mean_left_control, right_heel_marker_x_civ_left_control, {'color', color_left_control, 'linewidth', 5}, 1);
    shadedErrorBar(time_normalized, right_heel_marker_x_mean_left_positive, right_heel_marker_x_civ_left_positive, {'color', color_left_positive, 'linewidth', 5}, 1);
    shadedErrorBar(time_normalized, right_heel_marker_x_mean_left_negative, right_heel_marker_x_civ_left_negative, {'color', color_left_negative, 'linewidth', 5}, 1);
    xlabel('time')
end

return
% figure; axes; hold on; title('left foot medial-lateral heel marker')
% plot(time_normalized, left_heel_marker_x_mean_left_control, 'color', color_left_control, 'linewidth', 5);
% plot(time_normalized, left_heel_marker_x_mean_left_positive, 'color', color_left_positive, 'linewidth', 5);
% plot(time_normalized, left_heel_marker_x_mean_left_negative, 'color', color_left_negative, 'linewidth', 5);
% xlabel('time')
% 
% figure; axes; hold on; title('right foot medial-lateral heel marker')
% plot(time_normalized, right_heel_marker_x_mean_right_control, 'color', color_right_control, 'linewidth', 5);
% plot(time_normalized, right_heel_marker_x_mean_right_positive, 'color', color_right_positive, 'linewidth', 5);
% plot(time_normalized, right_heel_marker_x_mean_right_negative, 'color', color_right_negative, 'linewidth', 5);
% xlabel('time')

%% 
% plot CoP and heel marker in same figure
figure; axes; hold on; title('left foot medial-lateral CoP and heel marker')
shadedErrorBar(time_normalized, left_cop_x_mean_left_control, left_cop_x_civ_left_control, {'color', color_left_control, 'linewidth', 5}, 1);
shadedErrorBar(time_normalized, left_heel_marker_x_mean_left_control, left_heel_marker_x_civ_left_control, {'color', color_left_control, 'linewidth', 5}, 1);
xlabel('time')

% plot CoP and heel marker in same figure
figure; axes; hold on; title('right foot medial-lateral CoP and heel marker')
shadedErrorBar(time_normalized, right_cop_x_mean_right_control, right_cop_x_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
shadedErrorBar(time_normalized, right_heel_marker_x_mean_right_control, right_heel_marker_x_civ_right_control, {'color', color_right_control, 'linewidth', 5}, 1);
xlabel('time')

% plot CoP and heel marker in same figure
figure; axes; hold on; title('left foot medial-lateral CoP and heel marker')
shadedErrorBar(time_normalized, left_cop_x_mean_left_positive, left_cop_x_civ_left_positive, {'color', color_left_positive, 'linewidth', 5}, 1);
shadedErrorBar(time_normalized, left_heel_marker_x_mean_left_positive, left_heel_marker_x_civ_left_positive, {'color', color_left_positive, 'linewidth', 5}, 1);
xlabel('time')

% plot CoP and heel marker in same figure
figure; axes; hold on; title('right foot medial-lateral CoP and heel marker')
shadedErrorBar(time_normalized, right_cop_x_mean_right_positive, right_cop_x_civ_right_positive, {'color', color_right_positive, 'linewidth', 5}, 1);
shadedErrorBar(time_normalized, right_heel_marker_x_mean_right_positive, right_heel_marker_x_civ_right_positive, {'color', color_right_positive, 'linewidth', 5}, 1);
xlabel('time')

% plot CoP and heel marker in same figure
figure; axes; hold on; title('left foot medial-lateral CoP and heel marker')
shadedErrorBar(time_normalized, left_cop_x_mean_left_negative, left_cop_x_civ_left_negative, {'color', color_left_negative, 'linewidth', 5}, 1);
shadedErrorBar(time_normalized, left_heel_marker_x_mean_left_negative, left_heel_marker_x_civ_left_negative, {'color', color_left_negative, 'linewidth', 5}, 1);
xlabel('time')

% plot CoP and heel marker in same figure
figure; axes; hold on; title('right foot medial-lateral CoP and heel marker')
shadedErrorBar(time_normalized, right_cop_x_mean_right_negative, right_cop_x_civ_right_negative, {'color', color_right_negative, 'linewidth', 5}, 1);
shadedErrorBar(time_normalized, right_heel_marker_x_mean_right_negative, right_heel_marker_x_civ_right_negative, {'color', color_right_negative, 'linewidth', 5}, 1);
xlabel('time')














