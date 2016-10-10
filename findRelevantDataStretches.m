% This script finds the stretches of data relevant to the experimental paradigm.


%% Choose Data Type
emg_present                     = 1;
% stimulus_type = 'gvs';
stimulus_type = 'visual';
% stimulus_type = 'none';

remove_crossover_steps          = 0;
visualize                       = 1;

% stretches_to_analyze = 'single stance';
stretches_to_analyze = 'full stride';

%% prepare
wait_times = [0 0.150 0.450];
wait_time_labels = {'0ms', '150ms', '450ms'};
load subjectInfo.mat;

trials_to_process = 1 : 20;
trials_to_process = 5;

total_positive_steps = [];
total_negative_steps = [];
total_perturbation_condition_list = [];
total_step_condition_list = [];

swing_foot_fz_zero_threshold = 20; % threshold for counting a vertical force reading as zero, in Nm
swing_foot_nonzero_stretch_length_threshold_time = 0.05; % the force plate should register non-zero for at less than this long
time_to_nearest_future_heelstrike_threshold = 0.30; % a heelstrike should happen less than this long after a trigger
duration_prelude = 0.4; % duration that should have zero vertical force prior to a heelstrike, to identify crossover


%% extract data
stretch_length_indices_forceplate = [];
condition_stance_foot_list = {};
condition_perturbation_list = {};
condition_delay_list = {};
condition_index_list = {};
stim_start_time_relative_to_stretch = [];

for i_trial = trials_to_process
    %% load data
    %
    % Load the trial data from disk. 
    %
    
    load(makeFileName(date, subject_id, 'walking', i_trial, 'markerTrajectories'));
    load(makeFileName(date, subject_id, 'walking', i_trial, 'labviewTrajectories'));
    load(makeFileName(date, subject_id, 'walking', i_trial, 'forceplateTrajectories'));
    load(makeFileName(date, subject_id, 'walking', i_trial, 'stepEvents'));
    if emg_present
        load(makeFileName(date, subject_id, 'walking', i_trial, 'emgTrajectories'));
    end

    % estimate forceplate sampling rate if not provided
    if isnan(sampling_rate_forceplate)
        sampling_rate_forceplate = 1 / median(diff(time_forceplate));
    end
    number_of_prelude_time_steps = ceil(duration_prelude * sampling_rate_forceplate);
    swing_foot_nonzero_stretch_length_threshold_indices = ceil(swing_foot_nonzero_stretch_length_threshold_time * sampling_rate_forceplate);

    % determine illusion
    illusion_trajectory = zeros(size(time_labview)); % 1 = RIGHT, -1 = LEFT
    if strcmp(stimulus_type, 'gvs')
        % use GVS_out_trajectory
        for i_time = 1 : length(time_labview)
            if GVS_out_trajectory(i_time) > 0
                % anode is on the right, cathode is on the left, illusory fall is towards the cathode, i.e. LEFT
                illusion_trajectory(i_time) = -1;
            end
            if GVS_out_trajectory(i_time) < 0
                % anode is on the left, cathode is on the right, illusory fall is towards the cathode, i.e. RIGHT
                illusion_trajectory(i_time) = 1;
            end
        end
%         stim_sent_trajectory = GVS_out_trajectory;
    elseif strcmp(stimulus_type, 'visual')
        for i_time = 1 : length(time_labview)
            if visual_scene_ml_translation__trajectory(i_time) > 0
                % angle change is positive, horizon rotates counter-clockwise, illusion is to the RIGHT
                illusion_trajectory(i_time) = 1;
            end
            if visual_scene_ml_translation__trajectory(i_time) < 0
                % angle change is negative, horizon rotates clockwise, illusion is to the LEFT
                illusion_trajectory(i_time) = -1;
            end
        end
%         stim_sent_trajectory = visual_scene_ml_translation__trajectory;
    end
    
    % get CoP trajectories
    left_copx_trajectory = left_forceplate_cop_Acw(:, 1);
    right_copx_trajectory = right_forceplate_cop_Acw(:, 1);
    left_fz_trajectory = left_forceplate_wrench_Acw(:, 3);
    right_fz_trajectory = right_forceplate_wrench_Acw(:, 3);

    %% find triggers
    %
    % Find the triggering events that indicate a stretch of interest. For perturbation experiments, this is the onset of
    % a perturbation. For unperturbed walking, this is any heelstrike.
    % The result is trigger_indices_labview.
    %
    if strcmp(stimulus_type, 'none')
        % use all touchdown events as triggers
        trigger_indices_labview = [left_touchdown_indices_labview right_touchdown_indices_labview];
    else
        % find the time steps where the stimulus state crosses a threshold
        stimulus_threshold = 0.5;
        trigger_indices_labview = find(diff(sign(stimulus_state_trajectory - stimulus_threshold)) > 0) + 1;
        
        %
        epsilon = 1e-5;
        stim_start_indices_labview = find(diff(sign(abs(illusion_trajectory) - epsilon)) > 0) + 1;
        trigger_indices_labview = trigger_indices_labview(1 : length(stim_start_indices_labview)); % in case a stim is triggered, but not recorded
    end

    % visualize triggers
    if visualize
        left_cop_x_trajectory_relevant = left_copx_trajectory; left_cop_x_trajectory_relevant(left_cop_x_trajectory_relevant==0) = NaN;
        right_cop_x_trajectory_relevant = right_copx_trajectory; right_cop_x_trajectory_relevant(right_cop_x_trajectory_relevant==0) = NaN;
        figure; axes; hold on
        plot(time_labview, stimulus_state_trajectory*0.02);
        plot(time_forceplate, left_cop_x_trajectory_relevant, 'linewidth', 2);
        plot(time_forceplate, right_cop_x_trajectory_relevant, 'linewidth', 2);
        plot(time_forceplate(left_touchdown_indices_forceplate), zeros(size(left_touchdown_indices_forceplate)), 'o')
        plot(time_forceplate(right_touchdown_indices_forceplate), zeros(size(right_touchdown_indices_forceplate)), 'o')
        plot(time_labview(trigger_indices_labview), zeros(size(trigger_indices_labview)), 'x')
        plot(time_labview(stim_start_indices_labview), zeros(size(stim_start_indices_labview)), 'x')
        legend('stimulus state', 'left cop', 'right cop', 'left touchdown', 'right touchdown', 'trigger', 'stim start')
    end

    %% extract event data
    %
    % For each trigger, determine the conditions and the relevant step events.
    % The result is 
    % stretch_start_indices_forceplate, stretch_end_indices_forceplate
    % condition_stance_foot_list, condition_perturbation_list, condition_delay_list, condition_index_list
    
    % for each trigger, extract conditions and relevant step events
    number_of_triggers = length(trigger_indices_labview);
    removal_flags = zeros(number_of_triggers, 1);
    if strcmp(stimulus_type, 'none')
        stretch_start_indices_forceplate = zeros(number_of_triggers, 1);
        stretch_end_indices_forceplate = zeros(number_of_triggers, 1);
        closest_heelstrike_distance_times = zeros(number_of_triggers, 1);
        condition_stance_foot_list = cell(number_of_triggers, 1);
        condition_perturbation_list = cell(number_of_triggers, 1);
        condition_delay_list = cell(number_of_triggers, 1);
        condition_index_list = cell(number_of_triggers, 1);
    else
        stim_start_indices_labview = [stim_start_indices_labview zeros(number_of_triggers, 1) zeros(number_of_triggers, 1) zeros(number_of_triggers, 1)]; %#ok<AGROW>
        stretch_start_indices_forceplate = zeros(number_of_triggers, 4);
        stretch_end_indices_forceplate = zeros(number_of_triggers, 4);
        closest_heelstrike_distance_times = zeros(number_of_triggers, 4);
        condition_stance_foot_list = cell(number_of_triggers, 4);
        condition_perturbation_list = cell(number_of_triggers, 4);
        condition_delay_list = cell(number_of_triggers, 4);
        condition_index_list = cell(number_of_triggers, 4);
    end

    for i_trigger = 1 : number_of_triggers
        % determine condition
        if strcmp(stimulus_type, 'none')
            condition_perturbation_list{i_trigger, 1} = 'CONTROL';
            condition_delay_list{i_trigger, 1} = 'CONTROL';
            condition_index_list{i_trigger, 1} = 'CONTROL';
        else
            % perturbation condition
            if illusion_trajectory(stim_start_indices_labview(i_trigger)) > 0
                condition_perturbation_list{i_trigger, 1} = 'ILLUSION_RIGHT';
                condition_perturbation_list{i_trigger, 2} = 'ILLUSION_RIGHT';
            elseif illusion_trajectory(stim_start_indices_labview(i_trigger)) < 0
                condition_perturbation_list{i_trigger, 1} = 'ILLUSION_LEFT';
                condition_perturbation_list{i_trigger, 2} = 'ILLUSION_LEFT';
            else
%                 disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_labview(trigger_indices_labview(i_trigger))) ' - no stim']);
            end
            condition_perturbation_list{i_trigger, 3} = 'CONTROL';
            condition_perturbation_list{i_trigger, 4} = 'CONTROL';
            
            % delay
            wait_time_stim = time_labview(stim_start_indices_labview(i_trigger)) - time_labview(trigger_indices_labview(i_trigger));
            [~, wait_condition_index] = min(abs(wait_times - wait_time_stim));
            condition_delay_list{i_trigger, 1} = wait_time_labels{wait_condition_index};
            condition_delay_list{i_trigger, 2} = wait_time_labels{wait_condition_index};
            condition_delay_list{i_trigger, 3} = 'CONTROL';
            condition_delay_list{i_trigger, 4} = 'CONTROL';
        end
        

        % find out which heelstrike triggered
        [distance_to_trigger_left_time, index_left] = min(abs(time_forceplate(left_touchdown_indices_forceplate) - time_labview(trigger_indices_labview(i_trigger))));
        [distance_to_trigger_right_time, index_right] = min(abs(time_forceplate(right_touchdown_indices_forceplate) - time_labview(trigger_indices_labview(i_trigger))));
        closest_heelstrike_distance_time = min([distance_to_trigger_left_time distance_to_trigger_right_time]);
        
        % confirm that the trigger happened less than 100ms (or whatever time's set) before the actual heelstrike
        if ~strcmp(stimulus_type, 'none')
            if closest_heelstrike_distance_time > time_to_nearest_future_heelstrike_threshold
                removal_flags(i_trigger) = 1;
%                 disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_labview(trigger_indices_labview(i_trigger))) ' - trigger not in sync with heelstrike']);
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
            if index_left == 1 || length(left_touchdown_indices_forceplate) < index_left + 1 || removal_flags(i_trigger) == 1
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
                last_left_foot_pushoff = min(left_pushoff_indices_forceplate(left_pushoff_indices_forceplate >= last_left_foot_heelstrike));
                this_left_foot_pushoff = min(left_pushoff_indices_forceplate(left_pushoff_indices_forceplate >= this_left_foot_heelstrike));
                next_left_foot_pushoff = min(left_pushoff_indices_forceplate(left_pushoff_indices_forceplate >= next_left_foot_heelstrike));
                
                last_right_foot_heelstrike = min(right_touchdown_indices_forceplate(right_touchdown_indices_forceplate >= last_left_foot_heelstrike));
                this_right_foot_heelstrike = min(right_touchdown_indices_forceplate(right_touchdown_indices_forceplate >= this_left_foot_heelstrike));
                next_right_foot_heelstrike = min(right_touchdown_indices_forceplate(right_touchdown_indices_forceplate >= next_left_foot_heelstrike));
                last_right_foot_pushoff = max(right_pushoff_indices_forceplate(right_pushoff_indices_forceplate <= last_left_foot_pushoff));
                this_right_foot_pushoff = max(right_pushoff_indices_forceplate(right_pushoff_indices_forceplate <= this_left_foot_pushoff));
                next_right_foot_pushoff = max(right_pushoff_indices_forceplate(right_pushoff_indices_forceplate <= next_left_foot_pushoff));
                
                % notify if events are not sorted properly
                if ~issorted ...
                      ( ...
                        [ ...
                          last_left_foot_heelstrike ...
                          last_right_foot_pushoff ...
                          last_right_foot_heelstrike ...
                          last_left_foot_pushoff ...
                          this_left_foot_heelstrike ...
                          this_right_foot_pushoff ...
                          this_right_foot_heelstrike ...
                          this_left_foot_pushoff ...
                          next_left_foot_heelstrike ...
                          next_right_foot_pushoff ...
                          next_right_foot_heelstrike ...
                          next_left_foot_pushoff ...
                        ] ...
                      );
                  disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_labview(trigger_indices_labview(i_trigger)))]);
                    
                end
                
                
                % check check
                if visualize
                    plot(time_forceplate([last_left_foot_heelstrike this_left_foot_heelstrike next_left_foot_heelstrike]), [0 0 0]-0.01, 'v', 'linewidth', 3);
                    plot(time_forceplate([last_left_foot_pushoff  this_left_foot_pushoff next_left_foot_pushoff]), [0 0 0]-0.01, '^', 'linewidth', 3);
                    plot(time_forceplate([last_right_foot_heelstrike this_right_foot_heelstrike next_right_foot_heelstrike]), [0 0 0]+0.01, 'v', 'linewidth', 3);
                    plot(time_forceplate([last_right_foot_pushoff  this_right_foot_pushoff next_right_foot_pushoff]), [0 0 0]+0.01, '^', 'linewidth', 3);
                end                
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
            if index_right == 1 || length(right_touchdown_indices_forceplate) < index_right + 1 || removal_flags(i_trigger) == 1
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
                last_right_foot_pushoff = min(right_pushoff_indices_forceplate(right_pushoff_indices_forceplate >= last_right_foot_heelstrike));
                this_right_foot_pushoff = min(right_pushoff_indices_forceplate(right_pushoff_indices_forceplate >= this_right_foot_heelstrike));
                next_right_foot_pushoff = min(right_pushoff_indices_forceplate(right_pushoff_indices_forceplate >= next_right_foot_heelstrike));
                
                last_left_foot_heelstrike = min(left_touchdown_indices_forceplate(left_touchdown_indices_forceplate >= last_right_foot_heelstrike));
                this_left_foot_heelstrike = min(left_touchdown_indices_forceplate(left_touchdown_indices_forceplate >= this_right_foot_heelstrike));
                next_left_foot_heelstrike = min(left_touchdown_indices_forceplate(left_touchdown_indices_forceplate >= next_right_foot_heelstrike));
                last_left_foot_pushoff = max(left_pushoff_indices_forceplate(left_pushoff_indices_forceplate <= last_right_foot_pushoff));
                this_left_foot_pushoff = max(left_pushoff_indices_forceplate(left_pushoff_indices_forceplate <= this_right_foot_pushoff));
                next_left_foot_pushoff = max(left_pushoff_indices_forceplate(left_pushoff_indices_forceplate <= next_right_foot_pushoff));
                
                
                % notify if events are not sorted properly
                if ~issorted ...
                      ( ...
                        [ ...
                          last_right_foot_heelstrike ...
                          last_left_foot_pushoff ...
                          last_left_foot_heelstrike ...
                          last_right_foot_pushoff ...
                          this_right_foot_heelstrike ...
                          this_left_foot_pushoff ...
                          this_left_foot_heelstrike ...
                          this_right_foot_pushoff ...
                          next_right_foot_heelstrike ...
                          next_left_foot_pushoff ...
                          next_left_foot_heelstrike ...
                          next_right_foot_pushoff ...
                        ] ...
                      );
                  disp(['Trial ' num2str(i_trial) ': Problem with order of events, please check trigger at ' num2str(time_labview(trigger_indices_labview(i_trigger)))]);
                    
                end
                
                if visualize
                    plot(time_forceplate([last_left_foot_heelstrike this_left_foot_heelstrike next_left_foot_heelstrike]), [0 0 0]-0.01, 'v', 'linewidth', 3);
                    plot(time_forceplate([last_left_foot_pushoff  this_left_foot_pushoff next_left_foot_pushoff]), [0 0 0]-0.01, '^', 'linewidth', 3);
                    plot(time_forceplate([last_right_foot_heelstrike this_right_foot_heelstrike next_right_foot_heelstrike]), [0 0 0]+0.01, 'v', 'linewidth', 3);
                    plot(time_forceplate([last_right_foot_pushoff  this_right_foot_pushoff next_right_foot_pushoff]), [0 0 0]+0.01, '^', 'linewidth', 3);
                end            
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
              ) || removal_flags(i_trigger) == 1
            % not all events are present
            removal_flags(i_trigger) = 1;
            last_left_foot_heelstrike = NaN;
            this_left_foot_heelstrike = NaN;
            next_left_foot_heelstrike = NaN;
            last_right_foot_heelstrike = NaN;
            this_right_foot_heelstrike = NaN;
            next_right_foot_heelstrike = NaN;
            
        else % all events are present
            % mark relevant event delimiters depending on wait time and triggering foot
            if strcmp(trigger_foot, 'right')
                % triggered by right foot heelstrike
                if strcmp(stimulus_type, 'none')
                    condition_stance_foot_list{i_trigger, 1} = 'RIGHT';
                    stretch_start_indices_forceplate(i_trigger, 1) = this_left_foot_pushoff;
                    stretch_end_indices_forceplate(i_trigger, 1) = this_left_foot_heelstrike;
                else
                    % this step
                    condition_index_list{i_trigger, 1} = 'ONE';
                    condition_stance_foot_list{i_trigger, 1} = 'STANCE_RIGHT';
                    % next step
                    condition_index_list{i_trigger, 2} = 'TWO';
                    condition_stance_foot_list{i_trigger, 2} = 'STANCE_LEFT';
                    % use two previous steps as control
                    condition_index_list{i_trigger, 3} = 'CONTROL';
                    condition_stance_foot_list{i_trigger, 3} = 'STANCE_RIGHT';
                    condition_index_list{i_trigger, 4} = 'CONTROL';
                    condition_stance_foot_list{i_trigger, 4} = 'STANCE_LEFT';
                    
                    if strcmp(stretches_to_analyze, 'single stance')
                        % this step
                        stretch_start_indices_forceplate(i_trigger, 1) = this_left_foot_pushoff;
                        stretch_end_indices_forceplate(i_trigger, 1) = this_left_foot_heelstrike;
                        % next step
                        stretch_start_indices_forceplate(i_trigger, 2) = this_right_foot_pushoff;
                        stretch_end_indices_forceplate(i_trigger, 2) = next_right_foot_heelstrike;
                        % use two previous steps as control
                        stretch_start_indices_forceplate(i_trigger, 3) = last_left_foot_pushoff;
                        stretch_end_indices_forceplate(i_trigger, 3) = last_left_foot_heelstrike;
                        stretch_start_indices_forceplate(i_trigger, 4) = last_right_foot_pushoff;
                        stretch_end_indices_forceplate(i_trigger, 4) = this_right_foot_heelstrike;
                        
                    elseif strcmp(stretches_to_analyze, 'full stride')
                        % this step
                        stretch_start_indices_forceplate(i_trigger, 1) = this_right_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 1) = this_left_foot_heelstrike;
                        % next step
                        stretch_start_indices_forceplate(i_trigger, 2) = this_left_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 2) = next_right_foot_heelstrike;
                        % use two previous steps as control
                        stretch_start_indices_forceplate(i_trigger, 3) = last_right_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 3) = last_left_foot_heelstrike;
                        stretch_start_indices_forceplate(i_trigger, 4) = last_left_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 4) = this_right_foot_heelstrike;
                    end
                    
                    % visualize
                    if visualize
                        plot(time_forceplate([stretch_start_indices_forceplate(i_trigger, 1) stretch_end_indices_forceplate(i_trigger, 1)]), [0 0]-0.01, 'r', 'linewidth', 3);
                        plot(time_forceplate([stretch_start_indices_forceplate(i_trigger, 2) stretch_end_indices_forceplate(i_trigger, 2)]), [0 0]+0.01, 'g', 'linewidth', 3);
                        plot(time_forceplate([stretch_start_indices_forceplate(i_trigger, 3) stretch_end_indices_forceplate(i_trigger, 3)]), [0 0]-0.01, 'b', 'linewidth', 3);
                        plot(time_forceplate([stretch_start_indices_forceplate(i_trigger, 4) stretch_end_indices_forceplate(i_trigger, 4)]), [0 0]+0.01, 'm', 'linewidth', 3);
                    end

                end
            elseif strcmp(trigger_foot, 'left')
                % triggered by left foot heelstrike
                if strcmp(stimulus_type, 'none')
                    condition_stance_foot_list{i_trigger, 1} = 'LEFT';
                    stretch_start_indices_forceplate(i_trigger, 1) = this_right_foot_pushoff;
                    stretch_end_indices_forceplate(i_trigger, 1) = this_right_foot_heelstrike;
                else
                    % this step
                    condition_index_list{i_trigger, 1} = 'ONE';
                    condition_stance_foot_list{i_trigger, 1} = 'STANCE_LEFT';
                    % next step
                    condition_index_list{i_trigger, 2} = 'TWO';
                    condition_stance_foot_list{i_trigger, 2} = 'STANCE_RIGHT';
                    % use two previous steps as control
                    condition_index_list{i_trigger, 3} = 'CONTROL';
                    condition_stance_foot_list{i_trigger, 3} = 'STANCE_LEFT';
                    condition_index_list{i_trigger, 4} = 'CONTROL';
                    condition_stance_foot_list{i_trigger, 4} = 'STANCE_RIGHT';
                    if strcmp(stretches_to_analyze, 'single stance')
                        % this step
                        stretch_start_indices_forceplate(i_trigger, 1) = this_right_foot_pushoff;
                        stretch_end_indices_forceplate(i_trigger, 1) = this_right_foot_heelstrike;
                        % next step
                        stretch_start_indices_forceplate(i_trigger, 2) = this_left_foot_pushoff;
                        stretch_end_indices_forceplate(i_trigger, 2) = next_left_foot_heelstrike;
                        % use two previous steps as control
                        stretch_start_indices_forceplate(i_trigger, 3) = last_right_foot_pushoff;
                        stretch_end_indices_forceplate(i_trigger, 3) = last_right_foot_heelstrike;
                        stretch_start_indices_forceplate(i_trigger, 4) = last_left_foot_pushoff;
                        stretch_end_indices_forceplate(i_trigger, 4) = this_left_foot_heelstrike;
                        
                    elseif strcmp(stretches_to_analyze, 'full stride')
                        % this step
                        stretch_start_indices_forceplate(i_trigger, 1) = this_left_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 1) = this_right_foot_heelstrike;
                        % next step
                        stretch_start_indices_forceplate(i_trigger, 2) = this_right_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 2) = next_left_foot_heelstrike;
                        % use two previous steps as control
                        stretch_start_indices_forceplate(i_trigger, 3) = last_left_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 3) = last_right_foot_heelstrike;
                        stretch_start_indices_forceplate(i_trigger, 4) = last_right_foot_heelstrike;
                        stretch_end_indices_forceplate(i_trigger, 4) = this_left_foot_heelstrike;
                    end
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
%                     % this step
%                     stretch_start_indices_forceplate(i_trigger, 1) = this_right_foot_pushoff;
%                     stretch_end_indices_forceplate(i_trigger, 1) = this_right_foot_heelstrike;
%                     
%                     % next step
%                     stretch_start_indices_forceplate(i_trigger, 2) = this_left_foot_pushoff;
%                     stretch_end_indices_forceplate(i_trigger, 2) = next_left_foot_heelstrike;
%                     
%                     % use two previous steps as control
%                     stretch_start_indices_forceplate(i_trigger, 3) = last_right_foot_pushoff;
%                     stretch_end_indices_forceplate(i_trigger, 3) = last_right_foot_heelstrike;
%                     
%                     stretch_start_indices_forceplate(i_trigger, 4) = last_left_foot_pushoff;
%                     stretch_end_indices_forceplate(i_trigger, 4) = this_left_foot_heelstrike;
                    
                    % visualize
                    if visualize
                        plot(time_forceplate([stretch_start_indices_forceplate(i_trigger, 1) stretch_end_indices_forceplate(i_trigger, 1)]), [0 0]+0.01, 'r', 'linewidth', 3);
                        plot(time_forceplate([stretch_start_indices_forceplate(i_trigger, 2) stretch_end_indices_forceplate(i_trigger, 2)]), [0 0]-0.01, 'g', 'linewidth', 3);
                        plot(time_forceplate([stretch_start_indices_forceplate(i_trigger, 3) stretch_end_indices_forceplate(i_trigger, 3)]), [0 0]+0.01, 'b', 'linewidth', 3);                    
                        plot(time_forceplate([stretch_start_indices_forceplate(i_trigger, 4) stretch_end_indices_forceplate(i_trigger, 4)]), [0 0]-0.01, 'm', 'linewidth', 3);
                    end                
%                     if strcmp(condition_delay_list{i_trigger, 1}, '0ms') || strcmp(condition_delay_list{i_trigger, 1}, '150ms')
%                         % this step is of interest
%                         condition_stance_foot_list{i_trigger, 1} = 'LEFT';
%                         condition_stance_foot_list{i_trigger, 2} = 'LEFT';
%                         stretch_start_indices_forceplate(i_trigger, 1) = this_left_foot_heelstrike;
%                         stretch_end_indices_forceplate(i_trigger, 1) = this_right_foot_heelstrike;
%                         stretch_start_indices_forceplate(i_trigger, 2) = last_left_foot_heelstrike;
%                         stretch_end_indices_forceplate(i_trigger, 2) = last_right_foot_heelstrike;
%                     elseif strcmp(condition_delay_list{i_trigger, 1}, '450ms')
%                         % next step is of interest
%                         condition_stance_foot_list{i_trigger, 1} = 'RIGHT';
%                         condition_stance_foot_list{i_trigger, 2} = 'RIGHT';
%                         stretch_start_indices_forceplate(i_trigger, 1) = this_right_foot_heelstrike;
%                         stretch_end_indices_forceplate(i_trigger, 1) = next_left_foot_heelstrike;
%                         stretch_start_indices_forceplate(i_trigger, 2) = last_right_foot_heelstrike;
%                         stretch_end_indices_forceplate(i_trigger, 2) = this_left_foot_heelstrike;
%                     else
%                         disp(['Trial ' num2str(i_trial) ': something went wrong at time ' num2str(time_forceplate(trigger_indices_labview(i_trigger))) ' - delay not defined']);
%                     end
                end
            else
                removal_flags(i_trigger) = 1;
                condition_stance_foot_list{i_trigger, 1} = 'UNCLEAR';
                condition_stance_foot_list{i_trigger, 2} = 'UNCLEAR';
                condition_stance_foot_list{i_trigger, 3} = 'UNCLEAR';
                condition_stance_foot_list{i_trigger, 4} = 'UNCLEAR';
            end
        end
    end

    % remove flagged triggers
    unflagged_indices = ~removal_flags;
    if ~strcmp(stimulus_type, 'none')
        % if we have stimulation, deal with labview data (...? what exactly is going on here?)
        stim_start_indices_labview = stim_start_indices_labview(unflagged_indices, :);
    end
    stretch_start_indices_forceplate = stretch_start_indices_forceplate(unflagged_indices, :);
    stretch_end_indices_forceplate = stretch_end_indices_forceplate(unflagged_indices, :);
    condition_stance_foot_list = condition_stance_foot_list(unflagged_indices, :);
    condition_perturbation_list = condition_perturbation_list(unflagged_indices, :);
    condition_delay_list = condition_delay_list(unflagged_indices, :);
    condition_index_list = condition_index_list(unflagged_indices, :);
    closest_heelstrike_distance_times = closest_heelstrike_distance_times(unflagged_indices, :);

    %% form stretches
    %
    % Now take the stretches of data between the previously identified events
    % Result is lists of start/end indices and conditions for each stretch
    % 
    % 
    
    
    % form stretches
    if ~strcmp(stimulus_type, 'none')
        stim_start_indices_labview = [stim_start_indices_labview(:, 1); stim_start_indices_labview(:, 2); stim_start_indices_labview(:, 3); stim_start_indices_labview(:, 4)];
        stretch_start_indices_forceplate = [stretch_start_indices_forceplate(:, 1); stretch_start_indices_forceplate(:, 2); stretch_start_indices_forceplate(:, 3); stretch_start_indices_forceplate(:, 4)];
        stretch_end_indices_forceplate = [stretch_end_indices_forceplate(:, 1); stretch_end_indices_forceplate(:, 2); stretch_end_indices_forceplate(:, 3); stretch_end_indices_forceplate(:, 4)];
        condition_stance_foot_list = [condition_stance_foot_list(:, 1); condition_stance_foot_list(:, 2); condition_stance_foot_list(:, 3); condition_stance_foot_list(:, 4)];
        condition_perturbation_list = [condition_perturbation_list(:, 1); condition_perturbation_list(:, 2); condition_perturbation_list(:, 3); condition_perturbation_list(:, 4)];
        condition_delay_list = [condition_delay_list(:, 1); condition_delay_list(:, 2); condition_delay_list(:, 3); condition_delay_list(:, 4)];
        condition_index_list = [condition_index_list(:, 1); condition_index_list(:, 2); condition_index_list(:, 3); condition_index_list(:, 4)];
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
        condition_perturbation_list = condition_perturbation_list(unflagged_indices, :);
        condition_delay_list = condition_delay_list(unflagged_indices, :);
        condition_index_list = condition_index_list(unflagged_indices, :);
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
        
        start_index_forceplate = stretch_start_indices_forceplate(i_stretch);
        end_index_forceplate = stretch_end_indices_forceplate(i_stretch);
        stretch_start_time = time_forceplate(start_index_forceplate);
        stretch_end_time = time_forceplate(end_index_forceplate);
        if remove_crossover_steps
            % determine if this step should be removed due to crossover of the contralateral foot onto the ipsilateral forceplate
            if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')
%                 stance_foot_cop_stretch = right_copx_trajectory(start_index_forceplate : end_index_forceplate);
%                 touchdown_index_forceplate = stretch_start_index_forceplate + find(stance_foot_cop_stretch~=0, 1, 'first') - 1;
                stance_foot_fz_step = right_fz_trajectory(start_index_forceplate : end_index_forceplate);
                swing_foot_fz_step = left_fz_trajectory(start_index_forceplate : end_index_forceplate);
%                 stance_foot_fz_prelude = right_fz_trajectory(stretch_start_index_forceplate-number_of_prelude_time_steps : touchdown_index_forceplate-1);
            elseif strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')
%                 stance_foot_cop_stretch = left_copx_trajectory(start_index_forceplate : end_index_forceplate);
%                 touchdown_index_forceplate = stretch_start_index_forceplate + find(stance_foot_cop_stretch~=0, 1, 'first') - 1;
                stance_foot_fz_step = left_fz_trajectory(start_index_forceplate : end_index_forceplate);
                swing_foot_fz_step = right_fz_trajectory(start_index_forceplate : end_index_forceplate);
%                 stance_foot_fz_prelude = left_fz_trajectory(touchdown_index_forceplate-number_of_prelude_time_steps : touchdown_index_forceplate-1);
            end
%             touchdown_time_forceplate = time_forceplate(touchdown_index_forceplate);
            time_extracted_forceplate = time_forceplate(start_index_forceplate : end_index_forceplate);
            
            % check if this stretch is usable
            number_of_indices_per_step(i_stretch) = length(time_extracted_forceplate);
            number_of_swing_foot_fz_zeros_stretch = sum(abs(swing_foot_fz_step) < swing_foot_fz_zero_threshold);
            number_of_swing_foot_fz_nonzeros_stretch = sum(abs(swing_foot_fz_step) >= swing_foot_fz_zero_threshold);
            
%             number_of_swing_foot_fz_zeros_prelude = sum(-stance_foot_fz_prelude < swing_foot_fz_zero_threshold);
            
            if number_of_swing_foot_fz_nonzeros_stretch > swing_foot_nonzero_stretch_length_threshold_indices
                removal_flags(i_stretch) = 1;
                disp(['excluding stretch index ' num2str(i_stretch) ' - cross over at time ' num2str(stretch_start_time)]);
            end
            
            
            
            
            
            
%             if number_of_swing_foot_fz_zeros_prelude < swing_foot_zero_stretch_length_threshold_indices
%                 removal_flags(i_stretch) = 1;
%                 disp(['excluding stretch index ' num2str(i_stretch) ' - cross over at time ' num2str(time_forceplate(stretch_start_index_forceplate))]);
%             end
%             if isempty(touchdown_index_forceplate)
%                 touchdown_index_forceplate = stretch_start_index_forceplate;
%                 touchdown_time_forceplate = time_forceplate(touchdown_index_forceplate);
%                 removal_flags(i_stretch) = 1;
%             end
        end
        
        
%         plot(stance_foot_fz_step, 'displayname', num2str(i_stretch));
%         plot(swing_foot_fz_step, 'displayname', num2str(i_stretch));
%         plot(stance_foot_cop_stretch, 'displayname', num2str(i_stretch));
        

        % time
        [~, start_index_mocap] = min(abs(time_mocap - stretch_start_time));
        [~, end_index_mocap] = min(abs(time_mocap - stretch_end_time));
        [~, start_index_labview] = min(abs(time_labview - stretch_start_time));
        [~, end_index_labview] = min(abs(time_labview - stretch_end_time));
        if emg_present
            [~, start_index_emg] = min(abs(time_emg - stretch_start_time));
            [~, end_index_emg] = min(abs(time_emg - stretch_end_time));
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
                stim_start_time = time_labview(stim_start_indices_labview(i_stretch));
                stim_start_time_relative_to_stretch(i_stretch) = stim_start_time - stretch_start_time;
            end
        end
        

        % extract force plate data
        left_cop_x_extracted_stretch = left_copx_trajectory(start_index_forceplate : end_index_forceplate);
        right_cop_x_extracted_stretch = right_copx_trajectory(start_index_forceplate : end_index_forceplate);

        if remove_crossover_steps && strcmp(stretches_to_analyze, 'single stance')
            % check for data correctness - one force plate data stretch should have no zeros
            if any(left_cop_x_extracted_stretch==0) & any(right_cop_x_extracted_stretch==0)
                removal_flags(i_stretch) = 1;
            end
        end
    end

    % remove flagged stretches
    unflagged_indices = ~removal_flags;
    condition_stance_foot_list = condition_stance_foot_list(unflagged_indices);
    condition_perturbation_list = condition_perturbation_list(unflagged_indices);
    condition_delay_list = condition_delay_list(unflagged_indices);
    condition_index_list = condition_index_list(unflagged_indices);
    start_indices_mocap = start_indices_mocap(unflagged_indices);
    end_indices_mocap = end_indices_mocap(unflagged_indices);
    start_indices_forceplate = start_indices_forceplate(unflagged_indices);
    end_indices_forceplate = end_indices_forceplate(unflagged_indices);
    start_indices_labview = start_indices_labview(unflagged_indices);
    end_indices_labview = end_indices_labview(unflagged_indices);
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
        'condition_perturbation_list', ...
        'condition_delay_list', ...
        'condition_index_list', ...
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

























