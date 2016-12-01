function analyzeStimulusResponse(varargin)
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    
    save_data                           = 1;
    
    % setup analysis pattern - this should be done by editing a file
    process_data_balance                = 1;
    process_data_armswing               = 1;
    process_data_forceplate             = 1;
    process_data_emg                    = 1;
    process_data_angles                 = 0;
    process_data_torques                = 0;
    
    wait_times = [0 0.150 0.450];
    wait_time_labels = {'0ms', '150ms', '450ms'};

    condition_labels = {'stance foot', 'perturbation', 'delay', 'index'};
    
    %% define conditions
    % order of condition is STANCE_FOOT, PERTURBATION, DELAY, INDEX
    % STANCE_FOOT: RIGHT, LEFT
    % PERTURBATION: POSITIVE, NEGATIVE
    % DELAY: 0, 150, 450
    % INDEX: ONE, TWO
    %
    % this should be loaded from a file


    % for phase-dependent GVS
    conditions_control = ...
      {
        'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'walking'; ...
        'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'CONTROL', 'walking'; ...
      };

    conditions_to_analyze = ...
      {
        'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'ONE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_RIGHT', '150ms', 'ONE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_RIGHT', '450ms', 'ONE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'ONE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_LEFT', '150ms', 'ONE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_LEFT', '450ms', 'ONE', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'TWO', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_RIGHT', '150ms', 'TWO', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_RIGHT', '450ms', 'TWO', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'TWO', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_LEFT', '150ms', 'TWO', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_LEFT', '450ms', 'TWO', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'THREE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_RIGHT', '150ms', 'THREE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_RIGHT', '450ms', 'THREE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'THREE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_LEFT', '150ms', 'THREE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_LEFT', '450ms', 'THREE', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'FOUR', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_RIGHT', '150ms', 'FOUR', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_RIGHT', '450ms', 'FOUR', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'FOUR', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_LEFT', '150ms', 'FOUR', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_LEFT', '450ms', 'FOUR', 'walking'; ...
      };

%     % for vision
%     conditions_control = ...
%       {
%         'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'walking'; ...
%         'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'CONTROL', 'walking'; ...
%       };
% 
%     conditions_to_analyze = ...
%       {
%         'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'ONE', 'walking'; ...
%         'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'ONE', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'ONE', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'ONE', 'walking'; ...
%         'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'TWO', 'walking'; ...
%         'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'TWO', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'TWO', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'TWO', 'walking'; ...
%         'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'THREE', 'walking'; ...
%         'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'THREE', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'THREE', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'THREE', 'walking'; ...
%         'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'FOUR', 'walking'; ...
%         'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'FOUR', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'FOUR', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'FOUR', 'walking'; ...
%       };
    
    number_of_conditions_control = size(conditions_control, 1);
    number_of_conditions_to_analyze = size(conditions_to_analyze, 1);
    second_to_minute = 1/60;

    number_of_time_steps_normalized = 100;
    
    % determine which control condition applies to which perturbation condition
    applicable_control_condition_indices = zeros(number_of_conditions_to_analyze, 1);
    for i_condition = 1 : number_of_conditions_to_analyze
        if strcmp(conditions_to_analyze(i_condition, 1), 'STANCE_LEFT')
            applicable_control_condition_indices(i_condition) = 1;
        elseif strcmp(conditions_to_analyze(i_condition, 1), 'STANCE_RIGHT')
            applicable_control_condition_indices(i_condition) = 2;
        end
    end

    %% extract data
    %
    % This section loads the trial data from disk and extracts the stretches of interest.
    % Stretch start and end points are defined depending upon experimental conditions.
    % Stretch data is time-normalized and stored in arrays with size (number_of_time_steps_normalized x number_of_stretches)
    % named <variable name>_normalized_all
    %

    load('subjectInfo.mat', 'date', 'subject_id');
    
    % initialize containers
%     stretch_length_indices_forceplate = [];
    condition_stance_foot_list_all = {};
    condition_perturbation_list_all = {};
    condition_delay_list_all = {};
    condition_index_list_all = {};
    condition_experimental_list_all = {};
    origin_trial_list_all = [];
    origin_start_time_list_all = [];
    origin_end_time_list_all = [];
    step_times_all = [];
    pushoff_times_all = [];
    cadence_all = [];
    stim_start_time_relative_to_stretch_all = [];    
    
    if process_data_balance
        step_width_all = [];
        step_length_all = [];
        step_speed_all = [];
        
        c7_x_pos_normalized_all = [];
        c7_y_pos_normalized_all = [];
        c7_z_pos_normalized_all = [];
        lpsi_x_pos_normalized_all = [];
        lpsi_y_pos_normalized_all = [];
        lpsi_z_pos_normalized_all = [];
        rpsi_x_pos_normalized_all = [];
        rpsi_y_pos_normalized_all = [];
        rpsi_z_pos_normalized_all = [];
        lheel_x_pos_normalized_all = [];
        lheel_y_pos_normalized_all = [];
        lheel_z_pos_normalized_all = [];
        rheel_x_pos_normalized_all = [];
        rheel_y_pos_normalized_all = [];
        rheel_z_pos_normalized_all = [];
        
        lheel_x_pos_stancefoot_normalized_all = [];
        rheel_x_pos_stancefoot_normalized_all = [];
        lheel_x_pos_mpsis_normalized_all = [];
        rheel_x_pos_mpsis_normalized_all = [];
        
        foot_placement_world_all = [];
        foot_placement_stancefoot_all = [];
        foot_placement_mpsis_all = [];
        
        % angles
        lleg_angle_ml_normalized_all = [];
        rleg_angle_ml_normalized_all = [];
        trunk_angle_ap_normalized_all = [];
        trunk_angle_ml_normalized_all = [];
        
    end
    if process_data_armswing
        lelb_x_pos_normalized_all = [];
        lelb_y_pos_normalized_all = [];
        lelb_z_pos_normalized_all = [];
        lwra_x_pos_normalized_all = [];
        lwra_y_pos_normalized_all = [];
        lwra_z_pos_normalized_all = [];
        lwrb_x_pos_normalized_all = [];
        lwrb_y_pos_normalized_all = [];
        lwrb_z_pos_normalized_all = [];
        relb_x_pos_normalized_all = [];
        relb_y_pos_normalized_all = [];
        relb_z_pos_normalized_all = [];
        rwra_x_pos_normalized_all = [];
        rwra_y_pos_normalized_all = [];
        rwra_z_pos_normalized_all = [];
        rwrb_x_pos_normalized_all = [];
        rwrb_y_pos_normalized_all = [];
        rwrb_z_pos_normalized_all = [];
        
        % angles
        linclination_normalized_all = [];
        rinclination_normalized_all = [];
    end
    if process_data_forceplate
        cop_x_normalized_all = [];
        cop_x_stancefoot_normalized_all = [];
        cop_x_mpsis_normalized_all = [];
        f_x_normalized_all = [];
        f_z_normalized_all = [];
        m_y_normalized_all = [];
    end
    if process_data_emg
        lglutmed_normalized_all = [];
        ltibiant_normalized_all = [];
        lgastroc_normalized_all = [];
        lperolng_normalized_all = [];
        rglutmed_normalized_all = [];
        rtibiant_normalized_all = [];
        rgastroc_normalized_all = [];
        rperolng_normalized_all = [];
    end
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            % load data
            condition = condition_list{i_condition};
            load(['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'relevantDataStretches')]);
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
            number_of_stretches_trial = length(condition_stance_foot_list);

            condition_stance_foot_list_trial = condition_stance_foot_list;
            condition_perturbation_list_trial = condition_perturbation_list;
            condition_delay_list_trial = condition_delay_list;
            condition_index_list_trial = condition_index_list;
            condition_experimental_list_trial = condition_experimental_list;
            origin_trial_list_trial = zeros(number_of_stretches_trial, 1);
            origin_start_time_list_trial = zeros(number_of_stretches_trial, 1);
            origin_end_time_list_trial = zeros(number_of_stretches_trial, 1);
            step_times_trial = zeros(1, number_of_stretches_trial);
            pushoff_times_trial = zeros(1, number_of_stretches_trial);
            cadence_trial = zeros(1, number_of_stretches_trial);
            
            if process_data_balance

                % define markers and indices
                c7_marker = find(strcmp(marker_headers, 'C7'));
                lpsi_marker = find(strcmp(marker_headers, 'LPSI'));
                rpsi_marker = find(strcmp(marker_headers, 'RPSI'));
                lheel_marker = find(strcmp(marker_headers, 'LHEE'));
                rheel_marker = find(strcmp(marker_headers, 'RHEE'));

                c7_marker_indices = reshape([(c7_marker - 1) * 3 + 1; (c7_marker - 1) * 3 + 2; (c7_marker - 1) * 3 + 3], 1, length(c7_marker)*3);
                lpsi_marker_indices = reshape([(lpsi_marker - 1) * 3 + 1; (lpsi_marker - 1) * 3 + 2; (lpsi_marker - 1) * 3 + 3], 1, length(lpsi_marker)*3);
                rpsi_marker_indices = reshape([(rpsi_marker - 1) * 3 + 1; (rpsi_marker - 1) * 3 + 2; (rpsi_marker - 1) * 3 + 3], 1, length(rpsi_marker)*3);
                lheel_marker_indices = reshape([(lheel_marker - 1) * 3 + 1; (lheel_marker - 1) * 3 + 2; (lheel_marker - 1) * 3 + 3], 1, length(lheel_marker)*3);
                rheel_marker_indices = reshape([(rheel_marker - 1) * 3 + 1; (rheel_marker - 1) * 3 + 2; (rheel_marker - 1) * 3 + 3], 1, length(rheel_marker)*3);

                % rename relevant trajectories
                c7_x_pos_trajectory = marker_trajectories(:, c7_marker_indices(1));
                c7_y_pos_trajectory = marker_trajectories(:, c7_marker_indices(2));
                c7_z_pos_trajectory = marker_trajectories(:, c7_marker_indices(3));
                lpsi_x_pos_trajectory = marker_trajectories(:, lpsi_marker_indices(1));
                lpsi_y_pos_trajectory = marker_trajectories(:, lpsi_marker_indices(2));
                lpsi_z_pos_trajectory = marker_trajectories(:, lpsi_marker_indices(3));
                rpsi_x_pos_trajectory = marker_trajectories(:, rpsi_marker_indices(1));
                rpsi_y_pos_trajectory = marker_trajectories(:, rpsi_marker_indices(2));
                rpsi_z_pos_trajectory = marker_trajectories(:, rpsi_marker_indices(3));
                lheel_x_pos_trajectory = marker_trajectories(:, lheel_marker_indices(1));
                lheel_y_pos_trajectory = marker_trajectories(:, lheel_marker_indices(2));
                lheel_z_pos_trajectory = marker_trajectories(:, lheel_marker_indices(3));
                rheel_x_pos_trajectory = marker_trajectories(:, rheel_marker_indices(1));
                rheel_y_pos_trajectory = marker_trajectories(:, rheel_marker_indices(2));
                rheel_z_pos_trajectory = marker_trajectories(:, rheel_marker_indices(3));
                
                % initialize containers
                step_width_trial = zeros(1, number_of_stretches_trial);
                step_length_trial = zeros(1, number_of_stretches_trial);
                step_speed_trial = zeros(1, number_of_stretches_trial);
                c7_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                c7_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                c7_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lpsi_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lpsi_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lpsi_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rpsi_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rpsi_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rpsi_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lheel_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lheel_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lheel_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rheel_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rheel_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rheel_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lheel_x_pos_stancefoot_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rheel_x_pos_stancefoot_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lheel_x_pos_mpsis_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rheel_x_pos_mpsis_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lleg_angle_ml_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rleg_angle_ml_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                trunk_angle_ap_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                trunk_angle_ml_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                foot_placement_world_trial = zeros(1, number_of_stretches_trial);
                foot_placement_stancefoot_trial = zeros(1, number_of_stretches_trial);
                foot_placement_mpsis_trial = zeros(1, number_of_stretches_trial);
                
            end            
            if process_data_armswing
                
                % define markers and indices
                lelb_marker = find(strcmp(marker_headers, 'LELB'));
                lwra_marker = find(strcmp(marker_headers, 'LWRA'));
                lwrb_marker = find(strcmp(marker_headers, 'LWRB'));
                relb_marker = find(strcmp(marker_headers, 'RELB'));
                rwra_marker = find(strcmp(marker_headers, 'RWRA'));
                rwrb_marker = find(strcmp(marker_headers, 'RWRB'));
                
                lelb_marker_indices = reshape([(lelb_marker - 1) * 3 + 1; (lelb_marker - 1) * 3 + 2; (lelb_marker - 1) * 3 + 3], 1, length(lelb_marker)*3);
                lwra_marker_indices = reshape([(lwra_marker - 1) * 3 + 1; (lwra_marker - 1) * 3 + 2; (lwra_marker - 1) * 3 + 3], 1, length(lwra_marker)*3);
                lwrb_marker_indices = reshape([(lwrb_marker - 1) * 3 + 1; (lwrb_marker - 1) * 3 + 2; (lwrb_marker - 1) * 3 + 3], 1, length(lwrb_marker)*3);
                relb_marker_indices = reshape([(relb_marker - 1) * 3 + 1; (relb_marker - 1) * 3 + 2; (relb_marker - 1) * 3 + 3], 1, length(relb_marker)*3);
                rwra_marker_indices = reshape([(rwra_marker - 1) * 3 + 1; (rwra_marker - 1) * 3 + 2; (rwra_marker - 1) * 3 + 3], 1, length(rwra_marker)*3);
                rwrb_marker_indices = reshape([(rwrb_marker - 1) * 3 + 1; (rwrb_marker - 1) * 3 + 2; (rwrb_marker - 1) * 3 + 3], 1, length(rwrb_marker)*3);
                
                % rename relevant trajectories
                lelb_x_pos_trajectory = marker_trajectories(:, lelb_marker_indices(1));
                lelb_y_pos_trajectory = marker_trajectories(:, lelb_marker_indices(2));
                lelb_z_pos_trajectory = marker_trajectories(:, lelb_marker_indices(3));
                lwra_x_pos_trajectory = marker_trajectories(:, lwra_marker_indices(1));
                lwra_y_pos_trajectory = marker_trajectories(:, lwra_marker_indices(2));
                lwra_z_pos_trajectory = marker_trajectories(:, lwra_marker_indices(3));
                lwrb_x_pos_trajectory = marker_trajectories(:, lwrb_marker_indices(1));
                lwrb_y_pos_trajectory = marker_trajectories(:, lwrb_marker_indices(2));
                lwrb_z_pos_trajectory = marker_trajectories(:, lwrb_marker_indices(3));
                relb_x_pos_trajectory = marker_trajectories(:, relb_marker_indices(1));
                relb_y_pos_trajectory = marker_trajectories(:, relb_marker_indices(2));
                relb_z_pos_trajectory = marker_trajectories(:, relb_marker_indices(3));
                rwra_x_pos_trajectory = marker_trajectories(:, rwra_marker_indices(1));
                rwra_y_pos_trajectory = marker_trajectories(:, rwra_marker_indices(2));
                rwra_z_pos_trajectory = marker_trajectories(:, rwra_marker_indices(3));
                rwrb_x_pos_trajectory = marker_trajectories(:, rwrb_marker_indices(1));
                rwrb_y_pos_trajectory = marker_trajectories(:, rwrb_marker_indices(2));
                rwrb_z_pos_trajectory = marker_trajectories(:, rwrb_marker_indices(3));

                % initialize containers
                lelb_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lelb_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lelb_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lwra_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lwra_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lwra_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lwrb_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lwrb_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lwrb_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                relb_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                relb_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                relb_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rwra_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rwra_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rwra_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rwrb_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rwrb_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rwrb_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                
                linclination_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rinclination_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            
            end
            if process_data_forceplate
                % XXX this should be replaced with a system where the subject info file stores the available data for each trial
                forceplate_file_name = ['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'forceplateTrajectories') '.mat'];
                if exist(forceplate_file_name)
                    load(forceplate_file_name);
                    % initialize containers
                    force_plate_data_available = 1;
                else
                    force_plate_data_available = 0;
                    disp(['Warning: forceplate data not available for trial ' num2str(i_trial)]);
                end
                cop_x_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                cop_x_stancefoot_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                cop_x_mpsis_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                f_x_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                f_z_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                m_y_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            end
            if process_data_emg
                emg_file_name = ['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'emgTrajectories') '.mat'];
                load(emg_file_name);
                
                lglutmed_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LGLUTMED'));
                ltibiant_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LTIBIANT'));
                lgastroc_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LGASTROC'));
                lperolng_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LPEROLNG'));
                rglutmed_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RGLUTMED'));
                rtibiant_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RTIBIANT'));
                rgastroc_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RGASTROC'));
                rperolng_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RPEROLNG'));
                
                lglutmed_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                ltibiant_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lgastroc_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lperolng_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rglutmed_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rtibiant_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rgastroc_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rperolng_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);

            end
            
            % go through each stretch to extract data, time-normalize it and store it in lists
            for i_stretch = 1 : number_of_stretches_trial
                % time
                [~, start_index_mocap] = min(abs(time_mocap - stretch_start_times(i_stretch)));
                [~, pushoff_index_mocap] = min(abs(time_mocap - stretch_pushoff_times(i_stretch)));
                [~, end_index_mocap] = min(abs(time_mocap - stretch_end_times(i_stretch)));
                step_times_trial(i_stretch) = time_mocap(end_index_mocap) - time_mocap(start_index_mocap);
                pushoff_times_trial(i_stretch) = time_mocap(pushoff_index_mocap) - time_mocap(start_index_mocap);
                cadence_trial(i_stretch) = (step_times_trial(i_stretch) * second_to_minute)^(-1);
            
                % balance data
                if process_data_balance
                    % define times
                    time_extracted_mocap = time_mocap(start_index_mocap : end_index_mocap);
                    
                    % extract
                    c7_x_pos_extracted_stretch = c7_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    c7_y_pos_extracted_stretch = c7_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    c7_z_pos_extracted_stretch = c7_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    lpsi_x_pos_extracted_stretch = lpsi_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    lpsi_y_pos_extracted_stretch = lpsi_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    lpsi_z_pos_extracted_stretch = lpsi_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    rpsi_x_pos_extracted_stretch = rpsi_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    rpsi_y_pos_extracted_stretch = rpsi_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    rpsi_z_pos_extracted_stretch = rpsi_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    lheel_x_pos_extracted_stretch = lheel_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    lheel_y_pos_extracted_stretch = lheel_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    lheel_z_pos_extracted_stretch = lheel_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    rheel_x_pos_extracted_stretch = rheel_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    rheel_y_pos_extracted_stretch = rheel_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    rheel_z_pos_extracted_stretch = rheel_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    
                    % calculate vectors
                    mpsi_x_pos_extracted_stretch = (lpsi_x_pos_extracted_stretch + rpsi_x_pos_extracted_stretch) * 0.5;
                    mpsi_y_pos_extracted_stretch = (lpsi_y_pos_extracted_stretch + rpsi_y_pos_extracted_stretch) * 0.5;
                    mpsi_z_pos_extracted_stretch = (lpsi_z_pos_extracted_stretch + rpsi_z_pos_extracted_stretch) * 0.5;
                    lleg_vector_x = mpsi_x_pos_extracted_stretch - lheel_x_pos_extracted_stretch;
                    lleg_vector_y = mpsi_y_pos_extracted_stretch - lheel_y_pos_extracted_stretch;
                    lleg_vector_z = mpsi_z_pos_extracted_stretch - lheel_z_pos_extracted_stretch;
                    rleg_vector_x = mpsi_x_pos_extracted_stretch - rheel_x_pos_extracted_stretch;
                    rleg_vector_y = mpsi_y_pos_extracted_stretch - rheel_y_pos_extracted_stretch;
                    rleg_vector_z = mpsi_z_pos_extracted_stretch - rheel_z_pos_extracted_stretch;
                    trunk_vector_x = c7_x_pos_extracted_stretch - mpsi_x_pos_extracted_stretch;
                    trunk_vector_y = c7_y_pos_extracted_stretch - mpsi_y_pos_extracted_stretch;
                    trunk_vector_z = c7_z_pos_extracted_stretch - mpsi_z_pos_extracted_stretch;

                    % calculate angles
                    trunk_angle_ap_extracted_stretch = rad2deg(atan2(trunk_vector_y, trunk_vector_z));
                    trunk_angle_ml_extracted_stretch = rad2deg(atan2(trunk_vector_x, trunk_vector_z));
                    lleg_angle_ml_extracted_stretch = rad2deg(atan2(lleg_vector_x, lleg_vector_z));
                    rleg_angle_ml_extracted_stretch = rad2deg(atan2(rleg_vector_x, rleg_vector_z));
                    
                    
                    % define stance foot heel as spatial point of reference
                    if strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_RIGHT')
                        stance_foot_heel_x_initial = rheel_x_pos_extracted_stretch(1);
                        stance_foot_heel_y_initial = rheel_y_pos_extracted_stretch(1);
                    elseif strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_LEFT')
                        stance_foot_heel_x_initial = lheel_x_pos_extracted_stretch(1);
                        stance_foot_heel_y_initial = lheel_y_pos_extracted_stretch(1);
                    else
                        error('stance condition should be either "STANCE_LEFT or STANCE_RIGHT"');
                    end
                  
                    % calculate step width and length
                    if strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_RIGHT')
                        step_width_trial(i_stretch) = abs(lheel_x_pos_extracted_stretch(end) - rheel_x_pos_extracted_stretch(1));
                        step_length_trial(i_stretch) = lheel_y_pos_extracted_stretch(end) - rheel_y_pos_extracted_stretch(1);
                        foot_placement_world_trial(i_stretch) = lheel_x_pos_extracted_stretch(end);
                        foot_placement_stancefoot_trial(i_stretch) = lheel_x_pos_extracted_stretch(end) - rheel_x_pos_extracted_stretch(1);
                        foot_placement_mpsis_trial(i_stretch) = lheel_x_pos_extracted_stretch(end) - mpsi_x_pos_extracted_stretch(end);
                    elseif strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_LEFT')
                        step_width_trial(i_stretch) = abs(rheel_x_pos_extracted_stretch(end) - lheel_x_pos_extracted_stretch(1));
                        step_length_trial(i_stretch) = rheel_y_pos_extracted_stretch(end) - lheel_y_pos_extracted_stretch(1);
                        foot_placement_world_trial(i_stretch) = rheel_x_pos_extracted_stretch(end);
                        foot_placement_stancefoot_trial(i_stretch) = rheel_x_pos_extracted_stretch(end) - lheel_x_pos_extracted_stretch(1);
                        foot_placement_mpsis_trial(i_stretch) = rheel_x_pos_extracted_stretch(end) - mpsi_x_pos_extracted_stretch(end);
                    end

                    % normalize mocap data in time
                    time_normalized_mocap = linspace(time_extracted_mocap(1), time_extracted_mocap(end), number_of_time_steps_normalized);
                    lpsi_x_pos_normalized_stretch = spline(time_extracted_mocap, lpsi_x_pos_extracted_stretch, time_normalized_mocap);
                    lheel_x_pos_normalized_stretch = spline(time_extracted_mocap, lheel_x_pos_extracted_stretch, time_normalized_mocap);
                    lheel_y_pos_normalized_stretch = spline(time_extracted_mocap, lheel_y_pos_extracted_stretch, time_normalized_mocap);
                    rpsi_x_pos_normalized_stretch = spline(time_extracted_mocap, rpsi_x_pos_extracted_stretch, time_normalized_mocap);
                    rheel_x_pos_normalized_stretch = spline(time_extracted_mocap, rheel_x_pos_extracted_stretch, time_normalized_mocap);
                    rheel_y_pos_normalized_stretch = spline(time_extracted_mocap, rheel_y_pos_extracted_stretch, time_normalized_mocap);
                    mpsi_x_pos_normalized_stretch = spline(time_extracted_mocap, mpsi_x_pos_extracted_stretch, time_normalized_mocap);
                    mpsi_y_pos_normalized_stretch = spline(time_extracted_mocap, mpsi_y_pos_extracted_stretch, time_normalized_mocap);

                    trunk_angle_ap_normalized_stretch = spline(time_extracted_mocap, trunk_angle_ap_extracted_stretch, time_normalized_mocap);
                    trunk_angle_ml_normalized_stretch = spline(time_extracted_mocap, trunk_angle_ml_extracted_stretch, time_normalized_mocap);
                    lleg_angle_ml_normalized_stretch = spline(time_extracted_mocap, lleg_angle_ml_extracted_stretch, time_normalized_mocap);
                    rleg_angle_ml_normalized_stretch = spline(time_extracted_mocap, rleg_angle_ml_extracted_stretch, time_normalized_mocap);
                    
                    % store
                    lpsi_x_pos_normalized_trial(:, i_stretch) = lpsi_x_pos_normalized_stretch;
                    rpsi_x_pos_normalized_trial(:, i_stretch) = rpsi_x_pos_normalized_stretch;
                    lheel_x_pos_normalized_trial(:, i_stretch) = lheel_x_pos_normalized_stretch;
                    rheel_x_pos_normalized_trial(:, i_stretch) = rheel_x_pos_normalized_stretch;
                    lheel_y_pos_normalized_trial(:, i_stretch) = lheel_y_pos_normalized_stretch;
                    rheel_y_pos_normalized_trial(:, i_stretch) = rheel_y_pos_normalized_stretch;
                    trunk_angle_ap_normalized_trial(:, i_stretch) = trunk_angle_ap_normalized_stretch;
                    trunk_angle_ml_normalized_trial(:, i_stretch) = trunk_angle_ml_normalized_stretch;
                    lleg_angle_ml_normalized_trial(:, i_stretch) = lleg_angle_ml_normalized_stretch;
                    rleg_angle_ml_normalized_trial(:, i_stretch) = rleg_angle_ml_normalized_stretch;
                    
                    lheel_x_pos_stancefoot_normalized_trial(:, i_stretch) = lheel_x_pos_normalized_stretch - stance_foot_heel_x_initial;
                    rheel_x_pos_stancefoot_normalized_trial(:, i_stretch) = rheel_x_pos_normalized_stretch - stance_foot_heel_x_initial;
                    lheel_x_pos_mpsis_normalized_trial(:, i_stretch) = lheel_x_pos_normalized_stretch - mpsi_x_pos_normalized_stretch;
                    rheel_x_pos_mpsis_normalized_trial(:, i_stretch) = rheel_x_pos_normalized_stretch - mpsi_y_pos_normalized_stretch;

                    % time and origin
                    origin_trial_list_trial(i_stretch) = i_trial;
                    origin_start_time_list_trial(i_stretch) = time_mocap(start_index_mocap);
                    origin_end_time_list_trial(i_stretch) = time_mocap(end_index_mocap);
                    
                    % calculate step time and speed
                    step_speed_trial(i_stretch) = step_length_trial(i_stretch) / step_times_trial(i_stretch);
                end
                
                % armswing data
                if process_data_armswing
                    % define times
                    [~, start_index_mocap] = min(abs(time_mocap - stretch_start_times(i_stretch)));
                    [~, end_index_mocap] = min(abs(time_mocap - stretch_end_times(i_stretch)));
                    time_extracted_mocap = time_mocap(start_index_mocap : end_index_mocap);
                    
                    % extract
                    lelb_x_pos_extracted_stretch = lelb_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    lelb_y_pos_extracted_stretch = lelb_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    lelb_z_pos_extracted_stretch = lelb_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    lwra_x_pos_extracted_stretch = lwra_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    lwra_y_pos_extracted_stretch = lwra_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    lwra_z_pos_extracted_stretch = lwra_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    lwrb_x_pos_extracted_stretch = lwrb_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    lwrb_y_pos_extracted_stretch = lwrb_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    lwrb_z_pos_extracted_stretch = lwrb_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    relb_x_pos_extracted_stretch = relb_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    relb_y_pos_extracted_stretch = relb_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    relb_z_pos_extracted_stretch = relb_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    rwra_x_pos_extracted_stretch = rwra_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    rwra_y_pos_extracted_stretch = rwra_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    rwra_z_pos_extracted_stretch = rwra_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    rwrb_x_pos_extracted_stretch = rwrb_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    rwrb_y_pos_extracted_stretch = rwrb_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    rwrb_z_pos_extracted_stretch = rwrb_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    
                    % calculate vectors
                    lwrm_x_pos = (lwra_x_pos_extracted_stretch + lwrb_x_pos_extracted_stretch) * 0.5;
                    lwrm_y_pos = (lwra_y_pos_extracted_stretch + lwrb_y_pos_extracted_stretch) * 0.5;
                    lwrm_z_pos = (lwra_z_pos_extracted_stretch + lwrb_z_pos_extracted_stretch) * 0.5;
                    larm_vector_x = lelb_x_pos_extracted_stretch - lwrm_x_pos;
                    larm_vector_y = lelb_y_pos_extracted_stretch - lwrm_y_pos;
                    larm_vector_z = lelb_z_pos_extracted_stretch - lwrm_z_pos;
                    rwrm_x_pos = (rwra_x_pos_extracted_stretch + rwrb_x_pos_extracted_stretch) * 0.5;
                    rwrm_y_pos = (rwra_y_pos_extracted_stretch + rwrb_y_pos_extracted_stretch) * 0.5;
                    rwrm_z_pos = (rwra_z_pos_extracted_stretch + rwrb_z_pos_extracted_stretch) * 0.5;
                    rarm_vector_x = relb_x_pos_extracted_stretch - rwrm_x_pos;
                    rarm_vector_y = relb_y_pos_extracted_stretch - rwrm_y_pos;
                    rarm_vector_z = relb_z_pos_extracted_stretch - rwrm_z_pos;
                    
                    % calculate angles
                    vertical = [0; 0; 1];
                    linclination_extracted_stretch = zeros(length(time_extracted_mocap), 1);
                    rinclination_extracted_stretch = zeros(length(time_extracted_mocap), 1);
                    for i_time = 1 : length(time_extracted_mocap)
                        left_arm_vector = normVector([larm_vector_x(i_time); larm_vector_y(i_time); larm_vector_z(i_time)]);
                        left_arm_vector_projected_length = norm(left_arm_vector(1:2));
                        linclination_extracted_stretch(i_time) = rad2deg(atan2(left_arm_vector_projected_length, left_arm_vector(3)));
                        
                        right_arm_vector = normVector([rarm_vector_x(i_time); rarm_vector_y(i_time); rarm_vector_z(i_time)]);
                        right_arm_vector_projected_length = norm(right_arm_vector(1:2));
                        rinclination_extracted_stretch(i_time) = rad2deg(atan2(right_arm_vector_projected_length, right_arm_vector(3)));
                        
                    end
                    
                    % define stance foot heel as spatial point of reference
                    if strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_RIGHT')
                        stance_foot_heel_x_initial = rheel_x_pos_extracted_stretch(1);
                        stance_foot_heel_y_initial = rheel_y_pos_extracted_stretch(1);
                    elseif strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_LEFT')
                        stance_foot_heel_x_initial = lheel_x_pos_extracted_stretch(1);
                        stance_foot_heel_y_initial = lheel_y_pos_extracted_stretch(1);
                    else
                        error('stance condition should be either "STANCE_LEFT or STANCE_RIGHT"');
                    end
                  
                    % normalize data in time
                    time_normalized_mocap = linspace(time_extracted_mocap(1), time_extracted_mocap(end), number_of_time_steps_normalized);
                    linclination_normalized_stretch = spline(time_extracted_mocap, linclination_extracted_stretch, time_normalized_mocap);
                    rinclination_normalized_stretch = spline(time_extracted_mocap, rinclination_extracted_stretch, time_normalized_mocap);
                    
                    % store
                    linclination_normalized_trial(:, i_stretch) = linclination_normalized_stretch;
                    rinclination_normalized_trial(:, i_stretch) = rinclination_normalized_stretch;
                end

                % forceplate data
                if process_data_forceplate
                    if force_plate_data_available
                        % define times
                        [~, start_index_forceplate] = min(abs(time_forceplate - stretch_start_times(i_stretch)));
                        [~, end_index_forceplate] = min(abs(time_forceplate - stretch_end_times(i_stretch)));
                        time_extracted_forceplate = time_forceplate(start_index_forceplate : end_index_forceplate);

                        % extract
% XXX                        
% left_forceplate_cop_world = left_forceplate_cop_Acw;
% right_forceplate_cop_world = right_forceplate_cop_Acw;
% left_forceplate_wrench_world = left_forceplate_wrench_Acw;
% right_forceplate_wrench_world = right_forceplate_wrench_Acw;
% total_forceplate_wrench_world = left_forceplate_wrench_Acw + right_forceplate_wrench_Acw;
% 
% copx_trajectory = - total_forceplate_wrench_world(:, 5) ./ total_forceplate_wrench_world(:, 3);
% copy_trajectory = total_forceplate_wrench_world(:, 4) ./ total_forceplate_wrench_world(:, 3);
% total_forceplate_cop_world = [copx_trajectory copy_trajectory];
% 
% % re-zero CoP for low loads
% fz_threshold = 60;
% total_forceplate_low_load_indicator = (-total_forceplate_wrench_world(:, 3) < fz_threshold);
% total_forceplate_cop_world(total_forceplate_low_load_indicator, :) = 0;

                        
                        cop_x_extracted_stretch = total_forceplate_cop_world(start_index_forceplate : end_index_forceplate, 1);
                        f_x_extracted_stretch = total_forceplate_wrench_world(start_index_forceplate : end_index_forceplate, 1);
                        f_z_extracted_stretch = total_forceplate_wrench_world(start_index_forceplate : end_index_forceplate, 3);
                        m_y_extracted_stretch = total_forceplate_wrench_world(start_index_forceplate : end_index_forceplate, 5);

                        % normalize
                        time_normalized_forceplate = linspace(time_extracted_forceplate(1), time_extracted_forceplate(end), number_of_time_steps_normalized);
                        cop_x_normalized_stretch = spline(time_extracted_forceplate, cop_x_extracted_stretch, time_normalized_forceplate);
                        f_x_normalized_stretch = spline(time_extracted_forceplate, f_x_extracted_stretch, time_normalized_forceplate);
                        f_z_normalized_stretch = spline(time_extracted_forceplate, f_z_extracted_stretch, time_normalized_forceplate);
                        m_y_normalized_stretch = spline(time_extracted_forceplate, m_y_extracted_stretch, time_normalized_forceplate);

                        % use stance foot heel as spatial reference and store
                        cop_x_normalized_trial(:, i_stretch) = cop_x_normalized_stretch;
                        cop_x_stancefoot_normalized_trial(:, i_stretch) = cop_x_normalized_stretch - stance_foot_heel_x_initial;
                        cop_x_mpsis_normalized_trial(:, i_stretch) = cop_x_normalized_stretch - mpsi_x_pos_normalized_stretch;
                        f_x_normalized_trial(:, i_stretch) = f_x_normalized_stretch;
                        f_z_normalized_trial(:, i_stretch) = f_z_normalized_stretch;
                        m_y_normalized_trial(:, i_stretch) = m_y_normalized_stretch;
                    end
                end
                
                % emg data
                if process_data_emg
                    % define times
                    [~, start_index_emg] = min(abs(time_emg - stretch_start_times(i_stretch)));
                    [~, end_index_emg] = min(abs(time_emg - stretch_end_times(i_stretch)));
                    time_extracted_emg = time_emg(start_index_emg : end_index_emg);
                    
                    % extract
                    lglutmed_extracted_stretch = lglutmed_trajectory(start_index_emg : end_index_emg, 1);
                    ltibiant_extracted_stretch = ltibiant_trajectory(start_index_emg : end_index_emg, 1);
                    lgastroc_extracted_stretch = lgastroc_trajectory(start_index_emg : end_index_emg, 1);
                    lperolng_extracted_stretch = lperolng_trajectory(start_index_emg : end_index_emg, 1);
                    rglutmed_extracted_stretch = rglutmed_trajectory(start_index_emg : end_index_emg, 1);
                    rtibiant_extracted_stretch = rtibiant_trajectory(start_index_emg : end_index_emg, 1);
                    rgastroc_extracted_stretch = rgastroc_trajectory(start_index_emg : end_index_emg, 1);
                    rperolng_extracted_stretch = rperolng_trajectory(start_index_emg : end_index_emg, 1);
                    
                    % normalize
                    time_normalized_emg = linspace(time_extracted_emg(1), time_extracted_emg(end), number_of_time_steps_normalized);
                    lglutmed_normalized_stretch = spline(time_extracted_emg, lglutmed_extracted_stretch, time_normalized_emg);
                    ltibiant_normalized_stretch = spline(time_extracted_emg, ltibiant_extracted_stretch, time_normalized_emg);
                    lgastroc_normalized_stretch = spline(time_extracted_emg, lgastroc_extracted_stretch, time_normalized_emg);
                    lperolng_normalized_stretch = spline(time_extracted_emg, lperolng_extracted_stretch, time_normalized_emg);
                    rglutmed_normalized_stretch = spline(time_extracted_emg, rglutmed_extracted_stretch, time_normalized_emg);
                    rtibiant_normalized_stretch = spline(time_extracted_emg, rtibiant_extracted_stretch, time_normalized_emg);
                    rgastroc_normalized_stretch = spline(time_extracted_emg, rgastroc_extracted_stretch, time_normalized_emg);
                    rperolng_normalized_stretch = spline(time_extracted_emg, rperolng_extracted_stretch, time_normalized_emg);
                    
                    lglutmed_normalized_trial(:, i_stretch) = lglutmed_normalized_stretch;
                    ltibiant_normalized_trial(:, i_stretch) = ltibiant_normalized_stretch;
                    lgastroc_normalized_trial(:, i_stretch) = lgastroc_normalized_stretch;
                    lperolng_normalized_trial(:, i_stretch) = lperolng_normalized_stretch;
                    rglutmed_normalized_trial(:, i_stretch) = rglutmed_normalized_stretch;
                    rtibiant_normalized_trial(:, i_stretch) = rtibiant_normalized_stretch;
                    rgastroc_normalized_trial(:, i_stretch) = rgastroc_normalized_stretch;
                    rperolng_normalized_trial(:, i_stretch) = rperolng_normalized_stretch;
                    
                end
            end
            
            % append trial containers to total containers
            condition_stance_foot_list_all = [condition_stance_foot_list_all; condition_stance_foot_list_trial];
            condition_perturbation_list_all = [condition_perturbation_list_all; condition_perturbation_list_trial];
            condition_delay_list_all = [condition_delay_list_all; condition_delay_list_trial];
            condition_index_list_all = [condition_index_list_all; condition_index_list_trial];
            condition_experimental_list_all = [condition_experimental_list_all; condition_experimental_list_trial];
            origin_trial_list_all = [origin_trial_list_all; origin_trial_list_trial];
            origin_start_time_list_all = [origin_start_time_list_all; origin_start_time_list_trial];
            origin_end_time_list_all = [origin_end_time_list_all; origin_end_time_list_trial];
            step_times_all = [step_times_all step_times_trial];
            pushoff_times_all = [pushoff_times_all pushoff_times_trial];
            cadence_all = [cadence_all cadence_trial];

            if process_data_balance
                step_width_all = [step_width_all step_width_trial];
                step_length_all = [step_length_all step_length_trial];
                step_speed_all = [step_speed_all step_speed_trial];
                
                lpsi_x_pos_normalized_all = [lpsi_x_pos_normalized_all lpsi_x_pos_normalized_trial];
                rpsi_x_pos_normalized_all = [rpsi_x_pos_normalized_all rpsi_x_pos_normalized_trial];
                lheel_x_pos_normalized_all = [lheel_x_pos_normalized_all lheel_x_pos_normalized_trial];
                rheel_x_pos_normalized_all = [rheel_x_pos_normalized_all rheel_x_pos_normalized_trial];
                lheel_y_pos_normalized_all = [lheel_y_pos_normalized_all lheel_y_pos_normalized_trial];
                rheel_y_pos_normalized_all = [rheel_y_pos_normalized_all rheel_y_pos_normalized_trial];
                
                lheel_x_pos_stancefoot_normalized_all = [lheel_x_pos_stancefoot_normalized_all lheel_x_pos_stancefoot_normalized_trial];
                rheel_x_pos_stancefoot_normalized_all = [rheel_x_pos_stancefoot_normalized_all rheel_x_pos_stancefoot_normalized_trial];
                lheel_x_pos_mpsis_normalized_all = [lheel_x_pos_mpsis_normalized_all lheel_x_pos_mpsis_normalized_trial];
                rheel_x_pos_mpsis_normalized_all = [rheel_x_pos_mpsis_normalized_all rheel_x_pos_mpsis_normalized_trial];

                trunk_angle_ap_normalized_all = [trunk_angle_ap_normalized_all trunk_angle_ap_normalized_trial];
                trunk_angle_ml_normalized_all = [trunk_angle_ml_normalized_all trunk_angle_ml_normalized_trial];
                lleg_angle_ml_normalized_all = [lleg_angle_ml_normalized_all lleg_angle_ml_normalized_trial];
                rleg_angle_ml_normalized_all = [rleg_angle_ml_normalized_all rleg_angle_ml_normalized_trial];
                
                foot_placement_world_all = [foot_placement_world_all foot_placement_world_trial];
                foot_placement_stancefoot_all = [foot_placement_stancefoot_all foot_placement_stancefoot_trial];
                foot_placement_mpsis_all = [foot_placement_mpsis_all foot_placement_mpsis_trial];
            end

            if process_data_armswing
                linclination_normalized_all = [linclination_normalized_all linclination_normalized_trial];
                rinclination_normalized_all = [rinclination_normalized_all rinclination_normalized_trial];
            end
                
            if process_data_forceplate
                cop_x_normalized_all = [cop_x_normalized_all cop_x_normalized_trial];
                cop_x_stancefoot_normalized_all = [cop_x_stancefoot_normalized_all cop_x_stancefoot_normalized_trial];
                cop_x_mpsis_normalized_all = [cop_x_mpsis_normalized_all cop_x_mpsis_normalized_trial];
                f_x_normalized_all = [f_x_normalized_all f_x_normalized_trial];
                f_z_normalized_all = [f_z_normalized_all f_z_normalized_trial];
                m_y_normalized_all = [m_y_normalized_all m_y_normalized_trial];
            end

            if process_data_emg
                lglutmed_normalized_all = [lglutmed_normalized_all lglutmed_normalized_trial];
                ltibiant_normalized_all = [ltibiant_normalized_all ltibiant_normalized_trial];
                lgastroc_normalized_all = [lgastroc_normalized_all lgastroc_normalized_trial];
                lperolng_normalized_all = [lperolng_normalized_all lperolng_normalized_trial];
                rglutmed_normalized_all = [rglutmed_normalized_all rglutmed_normalized_trial];
                rtibiant_normalized_all = [rtibiant_normalized_all rtibiant_normalized_trial];
                rgastroc_normalized_all = [rgastroc_normalized_all rgastroc_normalized_trial];
                rperolng_normalized_all = [rperolng_normalized_all rperolng_normalized_trial];
            end
            
            disp(['Condition ' condition ', Trial ' num2str(i_trial) ' done extracted ' num2str(length(step_times_trial)) ' stretches']);
        end
    end
    number_of_stretches = length(step_times_all);
    
    
    %% extract conditions
    %
    % Compile information about which stretch belongs to which condition
    %
    % allocate indicators
    conditions_control_indicators = false(number_of_stretches, number_of_conditions_control);
    conditions_to_analyze_indicators = false(number_of_stretches, number_of_conditions_to_analyze);

    % extract control condition indicators
    for i_condition = 1 : number_of_conditions_control
        stance_foot_indicator = strcmp(condition_stance_foot_list_all, conditions_control(i_condition, 1));
        perturbation_indicator = strcmp(condition_perturbation_list_all, conditions_control(i_condition, 2));
        delay_indicator = strcmp(condition_delay_list_all, conditions_control(i_condition, 3));
        index_indicator = strcmp(condition_index_list_all, conditions_control(i_condition, 4));

        this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator;
        conditions_control_indicators(:, i_condition) = this_condition_indicator;
    end

    % extract indicators for conditions to analyze
    for i_condition = 1 : number_of_conditions_to_analyze
        stance_foot_indicator = strcmp(condition_stance_foot_list_all, conditions_to_analyze(i_condition, 1));
        perturbation_indicator = strcmp(condition_perturbation_list_all, conditions_to_analyze(i_condition, 2));
        delay_indicator = strcmp(condition_delay_list_all, conditions_to_analyze(i_condition, 3));
        index_indicator = strcmp(condition_index_list_all, conditions_to_analyze(i_condition, 4));

        this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator;
        conditions_to_analyze_indicators(:, i_condition) = this_condition_indicator;
    end

    % give feedback about number of trials per condition
    trials_per_condition_control = sum(conditions_control_indicators)';
    conditions_control_with_number = conditions_control;
    for i_condition = 1 : number_of_conditions_control
        conditions_control_with_number{i_condition, 5} = num2str(trials_per_condition_control(i_condition));
    end
    conditions_control_with_labels = [condition_labels 'number of stretches'; conditions_control_with_number];

    trials_per_condition_to_analyze = sum(conditions_to_analyze_indicators)';
    conditions_to_analyze_with_number = conditions_to_analyze;
    for i_condition = 1 : number_of_conditions_to_analyze
        conditions_to_analyze_with_number{i_condition, 5} = num2str(trials_per_condition_to_analyze(i_condition));
    end
    conditions_to_analyze_with_labels = [condition_labels 'number of stretches'; conditions_to_analyze_with_number];

    disp('Control conditions:')
    disp(conditions_control_with_labels);

    disp('Conditions to analyze:')
    disp(conditions_to_analyze_with_labels);

    disp(['Number of control stretches: ' num2str(sum(trials_per_condition_control))]);
    disp(['Number of stimulus stretches: ' num2str(sum(trials_per_condition_to_analyze))]);
    disp(['Number of unassigned stretches: ' num2str(number_of_stretches - sum(trials_per_condition_control) - sum(trials_per_condition_to_analyze))]);

    
    %% calculate responses
    if process_data_balance
        % calculate control means
        step_width_control_means = zeros(1, number_of_conditions_control);
        foot_placement_world_control_means = zeros(1, number_of_conditions_control);
        foot_placement_stancefoot_control_means = zeros(1, number_of_conditions_control);
        foot_placement_mpsis_control_means = zeros(1, number_of_conditions_control);
        lheel_x_pos_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rheel_x_pos_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lheel_y_pos_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rheel_y_pos_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lheel_x_pos_stancefoot_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rheel_x_pos_stancefoot_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lheel_x_pos_mpsis_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rheel_x_pos_mpsis_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        trunk_angle_ml_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lleg_angle_ml_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rleg_angle_ml_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            step_width_control_means(:, i_condition) = mean(step_width_all(:, condition_indicator), 2);
            foot_placement_world_control_means(:, i_condition) = mean(foot_placement_world_all(:, condition_indicator), 2);
            foot_placement_stancefoot_control_means(:, i_condition) = mean(foot_placement_stancefoot_all(:, condition_indicator), 2);
            foot_placement_mpsis_control_means(:, i_condition) = mean(foot_placement_mpsis_all(:, condition_indicator), 2);
            lheel_x_pos_control_means(:, i_condition) = mean(lheel_x_pos_normalized_all(:, condition_indicator), 2);
            rheel_x_pos_control_means(:, i_condition) = mean(rheel_x_pos_normalized_all(:, condition_indicator), 2);
            lheel_y_pos_control_means(:, i_condition) = mean(lheel_y_pos_normalized_all(:, condition_indicator), 2);
            rheel_y_pos_control_means(:, i_condition) = mean(rheel_y_pos_normalized_all(:, condition_indicator), 2);
            lheel_x_pos_stancefoot_control_means(:, i_condition) = mean(lheel_x_pos_stancefoot_normalized_all(:, condition_indicator), 2);
            rheel_x_pos_stancefoot_control_means(:, i_condition) = mean(rheel_x_pos_stancefoot_normalized_all(:, condition_indicator), 2);
            lheel_x_pos_mpsis_control_means(:, i_condition) = mean(lheel_x_pos_mpsis_normalized_all(:, condition_indicator), 2);
            rheel_x_pos_mpsis_control_means(:, i_condition) = mean(rheel_x_pos_mpsis_normalized_all(:, condition_indicator), 2);
            trunk_angle_ml_control_means(:, i_condition) = mean(trunk_angle_ml_normalized_all(:, condition_indicator), 2);
            lleg_angle_ml_control_means(:, i_condition) = mean(lleg_angle_ml_normalized_all(:, condition_indicator), 2);
            rleg_angle_ml_control_means(:, i_condition) = mean(rleg_angle_ml_normalized_all(:, condition_indicator), 2);
        end
        
        % calculate stimulus responses
        step_width_response = zeros(1, length(step_width_all));
        foot_placement_world_response = zeros(1, length(step_width_all));
        foot_placement_stancefoot_response = zeros(1, length(step_width_all));
        foot_placement_mpsis_response = zeros(1, length(step_width_all));
        lheel_x_pos_response = zeros(size(lheel_x_pos_normalized_all));
        rheel_x_pos_response = zeros(size(rheel_x_pos_normalized_all));
        lheel_y_pos_response = zeros(size(lheel_y_pos_normalized_all));
        rheel_y_pos_response = zeros(size(rheel_y_pos_normalized_all));
        lheel_x_pos_stancefoot_response = zeros(size(lheel_x_pos_normalized_all));
        rheel_x_pos_stancefoot_response = zeros(size(rheel_x_pos_normalized_all));
        lheel_x_pos_mpsis_response = zeros(size(lheel_x_pos_normalized_all));
        rheel_x_pos_mpsis_response = zeros(size(rheel_x_pos_normalized_all));
        trunk_angle_ml_response = zeros(size(trunk_angle_ml_normalized_all));
        lleg_angle_ml_response = zeros(size(lleg_angle_ml_normalized_all));
        rleg_angle_ml_response = zeros(size(rleg_angle_ml_normalized_all));
        
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            step_width_response(condition_indicator) = step_width_all(condition_indicator) - step_width_control_means(:, i_condition);
            foot_placement_world_response(condition_indicator) = foot_placement_world_all(condition_indicator) - foot_placement_world_control_means(:, i_condition);
            foot_placement_stancefoot_response(condition_indicator) = foot_placement_stancefoot_all(condition_indicator) - foot_placement_stancefoot_control_means(:, i_condition);
            foot_placement_mpsis_response(condition_indicator) = foot_placement_mpsis_all(condition_indicator) - foot_placement_mpsis_control_means(:, i_condition);
            lheel_x_pos_response(:, condition_indicator) = lheel_x_pos_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_control_means(:, i_condition), 1, sum(condition_indicator));
            rheel_x_pos_response(:, condition_indicator) = rheel_x_pos_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_control_means(:, i_condition), 1, sum(condition_indicator));
            lheel_y_pos_response(:, condition_indicator) = lheel_y_pos_normalized_all(:, condition_indicator) - repmat(lheel_y_pos_control_means(:, i_condition), 1, sum(condition_indicator));
            rheel_y_pos_response(:, condition_indicator) = rheel_y_pos_normalized_all(:, condition_indicator) - repmat(rheel_y_pos_control_means(:, i_condition), 1, sum(condition_indicator));
            lheel_x_pos_stancefoot_response(:, condition_indicator) = lheel_x_pos_stancefoot_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_stancefoot_control_means(:, i_condition), 1, sum(condition_indicator));
            rheel_x_pos_stancefoot_response(:, condition_indicator) = rheel_x_pos_stancefoot_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_stancefoot_control_means(:, i_condition), 1, sum(condition_indicator));
            lheel_x_pos_mpsis_response(:, condition_indicator) = lheel_x_pos_mpsis_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_mpsis_control_means(:, i_condition), 1, sum(condition_indicator));
            rheel_x_pos_mpsis_response(:, condition_indicator) = rheel_x_pos_mpsis_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_mpsis_control_means(:, i_condition), 1, sum(condition_indicator));
            trunk_angle_ml_response(:, condition_indicator) = trunk_angle_ml_normalized_all(:, condition_indicator) - repmat(trunk_angle_ml_control_means(:, i_condition), 1, sum(condition_indicator));
            lleg_angle_ml_response(:, condition_indicator) = lleg_angle_ml_normalized_all(:, condition_indicator) - repmat(lleg_angle_ml_control_means(:, i_condition), 1, sum(condition_indicator));
            rleg_angle_ml_response(:, condition_indicator) = rleg_angle_ml_normalized_all(:, condition_indicator) - repmat(rleg_angle_ml_control_means(:, i_condition), 1, sum(condition_indicator));
        end
        for i_condition = 1 : number_of_conditions_to_analyze
            condition_indicator = conditions_to_analyze_indicators(:, i_condition);
            step_width_response(condition_indicator) = step_width_all(condition_indicator) - step_width_control_means(:, applicable_control_condition_indices(i_condition));
            foot_placement_world_response(condition_indicator) = foot_placement_world_all(condition_indicator) - foot_placement_world_control_means(:, applicable_control_condition_indices(i_condition));
            foot_placement_stancefoot_response(condition_indicator) = foot_placement_stancefoot_all(condition_indicator) - foot_placement_stancefoot_control_means(:, applicable_control_condition_indices(i_condition));
            foot_placement_mpsis_response(condition_indicator) = foot_placement_mpsis_all(condition_indicator) - foot_placement_mpsis_control_means(:, applicable_control_condition_indices(i_condition));
            lheel_x_pos_response(:, condition_indicator) = lheel_x_pos_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rheel_x_pos_response(:, condition_indicator) = rheel_x_pos_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lheel_y_pos_response(:, condition_indicator) = lheel_y_pos_normalized_all(:, condition_indicator) - repmat(lheel_y_pos_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rheel_y_pos_response(:, condition_indicator) = rheel_y_pos_normalized_all(:, condition_indicator) - repmat(rheel_y_pos_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lheel_x_pos_stancefoot_response(:, condition_indicator) = lheel_x_pos_stancefoot_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_stancefoot_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rheel_x_pos_stancefoot_response(:, condition_indicator) = rheel_x_pos_stancefoot_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_stancefoot_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lheel_x_pos_mpsis_response(:, condition_indicator) = lheel_x_pos_mpsis_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_mpsis_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rheel_x_pos_mpsis_response(:, condition_indicator) = rheel_x_pos_mpsis_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_mpsis_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            trunk_angle_ml_response(:, condition_indicator) = trunk_angle_ml_normalized_all(:, condition_indicator) - repmat(trunk_angle_ml_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lleg_angle_ml_response(:, condition_indicator) = lleg_angle_ml_normalized_all(:, condition_indicator) - repmat(lleg_angle_ml_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rleg_angle_ml_response(:, condition_indicator) = rleg_angle_ml_normalized_all(:, condition_indicator) - repmat(rleg_angle_ml_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
        end
        
    end
    
    if process_data_forceplate
        % calculate control means
        cop_x_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        cop_x_stancefoot_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        cop_x_mpsis_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        f_x_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        f_z_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        m_y_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        % removed these for now, as we're currently not using them. They're easy to add in later
%         lcop_x_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         rcop_x_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         fxl_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         fzl_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         myl_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         fxr_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         fzr_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         myr_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            cop_x_control_means(:, i_condition) = mean(cop_x_normalized_all(:, condition_indicator), 2);
            cop_x_stancefoot_control_means(:, i_condition) = mean(cop_x_stancefoot_normalized_all(:, condition_indicator), 2);
            cop_x_mpsis_control_means(:, i_condition) = mean(cop_x_mpsis_normalized_all(:, condition_indicator), 2);
            f_x_control_means(:, i_condition) = mean(f_x_normalized_all(:, condition_indicator), 2);
            f_z_control_means(:, i_condition) = mean(f_z_normalized_all(:, condition_indicator), 2);
            m_y_control_means(:, i_condition) = mean(m_y_normalized_all(:, condition_indicator), 2);
%             lcop_x_control_means(:, i_condition) = mean(lcop_x_normalized_all(:, condition_indicator), 2);
%             rcop_x_control_means(:, i_condition) = mean(rcop_x_normalized_all(:, condition_indicator), 2);
%             fxl_control_means(:, i_condition) = mean(fxl_normalized_all(:, condition_indicator), 2);
%             fzl_control_means(:, i_condition) = mean(fzl_normalized_all(:, condition_indicator), 2);
%             myl_control_means(:, i_condition) = mean(myl_normalized_all(:, condition_indicator), 2);
%             fxr_control_means(:, i_condition) = mean(fxr_normalized_all(:, condition_indicator), 2);
%             fzr_control_means(:, i_condition) = mean(fzr_normalized_all(:, condition_indicator), 2);
%             myr_control_means(:, i_condition) = mean(myr_normalized_all(:, condition_indicator), 2);
        end
        
        % calculate stimulus responses
        cop_x_response = zeros(size(cop_x_normalized_all));
        cop_x_stancefoot_response = zeros(size(cop_x_normalized_all));
        cop_x_mpsis_response = zeros(size(cop_x_normalized_all));
        f_x_response = zeros(size(f_x_normalized_all));
        f_z_response = zeros(size(f_z_normalized_all));
        m_y_response = zeros(size(m_y_normalized_all));
%         lcop_x_response = zeros(size(lcop_x_normalized_all));
%         rcop_x_response = zeros(size(lcop_x_normalized_all));
%         fxl_response = zeros(size(fxl_normalized_all));
%         fzl_response = zeros(size(fzl_normalized_all));
%         myl_response = zeros(size(myl_normalized_all));
%         fxr_response = zeros(size(fxl_normalized_all));
%         fzr_response = zeros(size(fzl_normalized_all));
%         myr_response = zeros(size(myl_normalized_all));
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            cop_x_response(:, condition_indicator) = cop_x_normalized_all(:, condition_indicator) - repmat(cop_x_control_means(:, i_condition), 1, sum(condition_indicator));
            cop_x_stancefoot_response(:, condition_indicator) = cop_x_stancefoot_normalized_all(:, condition_indicator) - repmat(cop_x_stancefoot_control_means(:, i_condition), 1, sum(condition_indicator));
            cop_x_mpsis_response(:, condition_indicator) = cop_x_mpsis_normalized_all(:, condition_indicator) - repmat(cop_x_mpsis_control_means(:, i_condition), 1, sum(condition_indicator));
            f_x_response(:, condition_indicator) = f_x_normalized_all(:, condition_indicator) - repmat(f_x_control_means(:, i_condition), 1, sum(condition_indicator));
            f_z_response(:, condition_indicator) = f_z_normalized_all(:, condition_indicator) - repmat(f_z_control_means(:, i_condition), 1, sum(condition_indicator));
            m_y_response(:, condition_indicator) = m_y_normalized_all(:, condition_indicator) - repmat(m_y_control_means(:, i_condition), 1, sum(condition_indicator));
%             fxl_response(:, condition_indicator) = fxl_normalized_all(:, condition_indicator) - repmat(fxl_control_means(:, i_condition), 1, sum(condition_indicator));
%             fzl_response(:, condition_indicator) = fzl_normalized_all(:, condition_indicator) - repmat(fzl_control_means(:, i_condition), 1, sum(condition_indicator));
%             myl_response(:, condition_indicator) = myl_normalized_all(:, condition_indicator) - repmat(myl_control_means(:, i_condition), 1, sum(condition_indicator));
%             lcop_x_response(:, condition_indicator) = lcop_x_normalized_all(:, condition_indicator) - repmat(lcop_x_control_means(:, i_condition), 1, sum(condition_indicator));
%             fxr_response(:, condition_indicator) = fxr_normalized_all(:, condition_indicator) - repmat(fxr_control_means(:, i_condition), 1, sum(condition_indicator));
%             fzr_response(:, condition_indicator) = fzr_normalized_all(:, condition_indicator) - repmat(fzr_control_means(:, i_condition), 1, sum(condition_indicator));
%             myr_response(:, condition_indicator) = myr_normalized_all(:, condition_indicator) - repmat(myr_control_means(:, i_condition), 1, sum(condition_indicator));
%             rcop_x_response(:, condition_indicator) = rcop_x_normalized_all(:, condition_indicator) - repmat(rcop_x_control_means(:, i_condition), 1, sum(condition_indicator));
        end
        for i_condition = 1 : number_of_conditions_to_analyze
            condition_indicator = conditions_to_analyze_indicators(:, i_condition);
            cop_x_response(:, condition_indicator) = cop_x_normalized_all(:, condition_indicator) - repmat(cop_x_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            cop_x_stancefoot_response(:, condition_indicator) = cop_x_stancefoot_normalized_all(:, condition_indicator) - repmat(cop_x_stancefoot_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            cop_x_mpsis_response(:, condition_indicator) = cop_x_mpsis_normalized_all(:, condition_indicator) - repmat(cop_x_mpsis_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            f_x_response(:, condition_indicator) = f_x_normalized_all(:, condition_indicator) - repmat(f_x_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            f_z_response(:, condition_indicator) = f_z_normalized_all(:, condition_indicator) - repmat(f_z_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            m_y_response(:, condition_indicator) = m_y_normalized_all(:, condition_indicator) - repmat(m_y_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             fxl_response(:, condition_indicator) = fxl_normalized_all(:, condition_indicator) - repmat(fxl_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             fzl_response(:, condition_indicator) = fzl_normalized_all(:, condition_indicator) - repmat(fzl_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             myl_response(:, condition_indicator) = myl_normalized_all(:, condition_indicator) - repmat(myl_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             lcop_x_response(:, condition_indicator) = lcop_x_normalized_all(:, condition_indicator) - repmat(lcop_x_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             fxr_response(:, condition_indicator) = fxr_normalized_all(:, condition_indicator) - repmat(fxr_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             fzr_response(:, condition_indicator) = fzr_normalized_all(:, condition_indicator) - repmat(fzr_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             myr_response(:, condition_indicator) = myr_normalized_all(:, condition_indicator) - repmat(myr_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             rcop_x_response(:, condition_indicator) = rcop_x_normalized_all(:, condition_indicator) - repmat(rcop_x_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
        end
        
    end
    
    if process_data_emg
        % calculate control means
        lglutmed_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        ltibiant_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lgastroc_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lperolng_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rglutmed_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rtibiant_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rgastroc_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rperolng_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            lglutmed_control_means(:, i_condition) = mean(lglutmed_normalized_all(:, condition_indicator), 2);
            ltibiant_control_means(:, i_condition) = mean(ltibiant_normalized_all(:, condition_indicator), 2);
            lgastroc_control_means(:, i_condition) = mean(lgastroc_normalized_all(:, condition_indicator), 2);
            lperolng_control_means(:, i_condition) = mean(lperolng_normalized_all(:, condition_indicator), 2);
            rglutmed_control_means(:, i_condition) = mean(rglutmed_normalized_all(:, condition_indicator), 2);
            rtibiant_control_means(:, i_condition) = mean(rtibiant_normalized_all(:, condition_indicator), 2);
            rgastroc_control_means(:, i_condition) = mean(rgastroc_normalized_all(:, condition_indicator), 2);
            rperolng_control_means(:, i_condition) = mean(rperolng_normalized_all(:, condition_indicator), 2);
        end
        
        % calculate stimulus responses
        lglutmed_response = zeros(size(lglutmed_normalized_all));
        ltibiant_response = zeros(size(ltibiant_normalized_all));
        lgastroc_response = zeros(size(lgastroc_normalized_all));
        lperolng_response = zeros(size(lperolng_normalized_all));
        rglutmed_response = zeros(size(rglutmed_normalized_all));
        rtibiant_response = zeros(size(rtibiant_normalized_all));
        rgastroc_response = zeros(size(rgastroc_normalized_all));
        rperolng_response = zeros(size(rperolng_normalized_all));
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            lglutmed_response(:, condition_indicator) = lglutmed_normalized_all(:, condition_indicator) - repmat(lglutmed_control_means(:, i_condition), 1, sum(condition_indicator));
            ltibiant_response(:, condition_indicator) = ltibiant_normalized_all(:, condition_indicator) - repmat(ltibiant_control_means(:, i_condition), 1, sum(condition_indicator));
            lgastroc_response(:, condition_indicator) = lgastroc_normalized_all(:, condition_indicator) - repmat(lgastroc_control_means(:, i_condition), 1, sum(condition_indicator));
            lperolng_response(:, condition_indicator) = lperolng_normalized_all(:, condition_indicator) - repmat(lperolng_control_means(:, i_condition), 1, sum(condition_indicator));
            rglutmed_response(:, condition_indicator) = rglutmed_normalized_all(:, condition_indicator) - repmat(rglutmed_control_means(:, i_condition), 1, sum(condition_indicator));
            rtibiant_response(:, condition_indicator) = rtibiant_normalized_all(:, condition_indicator) - repmat(rtibiant_control_means(:, i_condition), 1, sum(condition_indicator));
            rgastroc_response(:, condition_indicator) = rgastroc_normalized_all(:, condition_indicator) - repmat(rgastroc_control_means(:, i_condition), 1, sum(condition_indicator));
            rperolng_response(:, condition_indicator) = rperolng_normalized_all(:, condition_indicator) - repmat(rperolng_control_means(:, i_condition), 1, sum(condition_indicator));
        end
        for i_condition = 1 : number_of_conditions_to_analyze
            condition_indicator = conditions_to_analyze_indicators(:, i_condition);
            lglutmed_response(:, condition_indicator) = lglutmed_normalized_all(:, condition_indicator) - repmat(lglutmed_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            ltibiant_response(:, condition_indicator) = ltibiant_normalized_all(:, condition_indicator) - repmat(ltibiant_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lgastroc_response(:, condition_indicator) = lgastroc_normalized_all(:, condition_indicator) - repmat(lgastroc_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lperolng_response(:, condition_indicator) = lperolng_normalized_all(:, condition_indicator) - repmat(lperolng_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rglutmed_response(:, condition_indicator) = rglutmed_normalized_all(:, condition_indicator) - repmat(rglutmed_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rtibiant_response(:, condition_indicator) = rtibiant_normalized_all(:, condition_indicator) - repmat(rtibiant_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rgastroc_response(:, condition_indicator) = rgastroc_normalized_all(:, condition_indicator) - repmat(rgastroc_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rperolng_response(:, condition_indicator) = rperolng_normalized_all(:, condition_indicator) - repmat(rperolng_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
        end        
    end
    
    
    %% save data
    if save_data

        save ...
          ( ...
            ['analysis' filesep makeFileName(date, subject_id, 'resultsConditions')], ...
            'number_of_time_steps_normalized', ...
            'origin_trial_list_all', ...
            'origin_start_time_list_all', ...
            'origin_end_time_list_all', ...
            'conditions_control', ...
            'conditions_to_analyze', ...
            'condition_stance_foot_list_all', ...
            'condition_perturbation_list_all', ...
            'condition_delay_list_all', ...
            'condition_index_list_all', ...
            'condition_experimental_list_all', ...
            'conditions_control_indicators', ...
            'conditions_to_analyze_indicators', ...
            'step_times_all', ...
            'pushoff_times_all', ...
            'cadence_all' ...
          );

        if process_data_balance
            save ...
              ( ...
                ['analysis' filesep makeFileName(date, subject_id, 'resultsBalance')], ...
                'step_width_all', ...
                'step_length_all', ...
                'step_speed_all', ...
                'foot_placement_world_all', ...
                'foot_placement_stancefoot_all', ...
                'foot_placement_mpsis_all', ...
                'lheel_x_pos_normalized_all', ...
                'rheel_x_pos_normalized_all', ...
                'lheel_y_pos_normalized_all', ...
                'rheel_y_pos_normalized_all', ...
                'lheel_x_pos_stancefoot_normalized_trial', ...
                'rheel_x_pos_stancefoot_normalized_trial', ...
                'lheel_x_pos_mpsis_normalized_trial', ...
                'rheel_x_pos_mpsis_normalized_trial', ...
                'trunk_angle_ap_normalized_all', ...
                'trunk_angle_ml_normalized_all', ...
                'lleg_angle_ml_normalized_all', ...
                'rleg_angle_ml_normalized_all', ...
                'step_width_response', ...
                'foot_placement_world_response', ...
                'foot_placement_stancefoot_response', ...
                'foot_placement_mpsis_response', ...
                'lheel_x_pos_response', ...
                'rheel_x_pos_response', ...
                'lheel_y_pos_response', ...
                'rheel_y_pos_response', ...
                'lheel_x_pos_stancefoot_response', ...
                'rheel_x_pos_stancefoot_response', ...
                'lheel_x_pos_mpsis_response', ...
                'rheel_x_pos_mpsis_response', ...
                'trunk_angle_ml_response', ...
                'lleg_angle_ml_response', ...
                'rleg_angle_ml_response' ...
              );
        end

        if process_data_armswing
            save ...
              ( ...
                ['analysis' filesep makeFileName(date, subject_id, 'resultsArmswing')], ...
                'linclination_normalized_all', ...
                'rinclination_normalized_all' ...
              );
        end

        if process_data_forceplate
            save ...
              ( ...
                ['analysis' filesep makeFileName(date, subject_id, 'resultsForceplate')], ...
                'cop_x_normalized_all', ...
                'cop_x_response', ...
                'cop_x_stancefoot_normalized_all', ...
                'cop_x_stancefoot_response', ...
                'cop_x_mpsis_normalized_all', ...
                'cop_x_mpsis_response', ...
                'f_x_normalized_all', ...
                'f_x_response', ...
                'f_z_normalized_all', ...
                'f_z_response', ...
                'm_y_normalized_all', ...
                'm_y_response' ...
              );
        end
        
        if process_data_emg
            save ...
              ( ...
                ['analysis' filesep makeFileName(date, subject_id, 'resultsEmg')], ...
                'lglutmed_normalized_all', ...
                'ltibiant_normalized_all', ...
                'lgastroc_normalized_all', ...
                'lperolng_normalized_all', ...
                'rglutmed_normalized_all', ...
                'rtibiant_normalized_all', ...
                'rgastroc_normalized_all', ...
                'rperolng_normalized_all' ...
              );
        end        
        
    end    
    
    
    return
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % after here comes old code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    for i_trial = trials_to_process
        % load data
        load(makeFileName(date, subject_id, 'walking', i_trial, 'stepEvents'));
        load(makeFileName(date, subject_id, 'walking', i_trial, 'relevantDataStretches'));
        number_of_stretches_trial = length(condition_stance_foot_list);

        condition_stance_foot_list_trial = condition_stance_foot_list;
        condition_perturbation_list_trial = condition_perturbation_list;
        condition_delay_list_trial = condition_delay_list;
        condition_index_list_trial = condition_index_list;
        origin_trial_list_trial = zeros(number_of_stretches_trial, 1);
        origin_start_time_list_trial = zeros(number_of_stretches_trial, 1);
        origin_end_time_list_trial = zeros(number_of_stretches_trial, 1);
        step_times_trial = zeros(number_of_stretches_trial, 1);
        stim_start_time_relative_to_stretch_trial = stim_start_time_relative_to_stretch;

        if process_data_marker
            load(makeFileName(date, subject_id, 'walking', i_trial, 'markerTrajectories'));
            filter_order_low = 4;
            cutoff_frequency_low = 20; % in Hz
            [b_lowpass, a_lowpass] = butter(filter_order_low, cutoff_frequency_low/(sampling_rate_mocap/2), 'low');
            
            % define markers and indices
            c7_marker = find(strcmp(marker_headers, 'C7'));
            lpsi_marker = find(strcmp(marker_headers, 'LPSI'));
            rpsi_marker = find(strcmp(marker_headers, 'RPSI'));
            lheel_marker = find(strcmp(marker_headers, 'LHEE'));
            rheel_marker = find(strcmp(marker_headers, 'RHEE'));

            c7_marker_indices = reshape([(c7_marker - 1) * 3 + 1; (c7_marker - 1) * 3 + 2; (c7_marker - 1) * 3 + 3], 1, length(c7_marker)*3);
            lpsi_marker_indices = reshape([(lpsi_marker - 1) * 3 + 1; (lpsi_marker - 1) * 3 + 2; (lpsi_marker - 1) * 3 + 3], 1, length(lpsi_marker)*3);
            rpsi_marker_indices = reshape([(rpsi_marker - 1) * 3 + 1; (rpsi_marker - 1) * 3 + 2; (rpsi_marker - 1) * 3 + 3], 1, length(rpsi_marker)*3);
            lheel_marker_indices = reshape([(lheel_marker - 1) * 3 + 1; (lheel_marker - 1) * 3 + 2; (lheel_marker - 1) * 3 + 3], 1, length(lheel_marker)*3);
            rheel_marker_indices = reshape([(rheel_marker - 1) * 3 + 1; (rheel_marker - 1) * 3 + 2; (rheel_marker - 1) * 3 + 3], 1, length(rheel_marker)*3);

            % rename relevant trajectories
            c7_x_pos_trajectory = marker_trajectories(:, c7_marker_indices(1));
            c7_y_pos_trajectory = marker_trajectories(:, c7_marker_indices(2));
            c7_z_pos_trajectory = marker_trajectories(:, c7_marker_indices(3));
            lpsi_x_pos_trajectory = marker_trajectories(:, lpsi_marker_indices(1));
            lpsi_y_pos_trajectory = marker_trajectories(:, lpsi_marker_indices(2));
            lpsi_z_pos_trajectory = marker_trajectories(:, lpsi_marker_indices(3));
            rpsi_x_pos_trajectory = marker_trajectories(:, rpsi_marker_indices(1));
            rpsi_y_pos_trajectory = marker_trajectories(:, rpsi_marker_indices(2));
            rpsi_z_pos_trajectory = marker_trajectories(:, rpsi_marker_indices(3));
            lheel_x_pos_trajectory = marker_trajectories(:, lheel_marker_indices(1));
            lheel_y_pos_trajectory = marker_trajectories(:, lheel_marker_indices(2));
            lheel_z_pos_trajectory = marker_trajectories(:, lheel_marker_indices(3));
            rheel_x_pos_trajectory = marker_trajectories(:, rheel_marker_indices(1));
            rheel_y_pos_trajectory = marker_trajectories(:, rheel_marker_indices(2));
            rheel_z_pos_trajectory = marker_trajectories(:, rheel_marker_indices(3));
%             lpsi_x_vel_trajectory = deriveByTime(filtfilt(b_lowpass, a_lowpass, spline(time_mocap, lpsi_x_pos_trajectory, time_mocap)), sampling_rate_mocap^(-1));
%             rpsi_x_vel_trajectory = deriveByTime(filtfilt(b_lowpass, a_lowpass, spline(time_mocap, rpsi_x_pos_trajectory, time_mocap)), sampling_rate_mocap^(-1));
%             lpsi_x_acc_trajectory = deriveByTime(filtfilt(b_lowpass, a_lowpass, lpsi_x_vel_trajectory), sampling_rate_mocap^(-1));
%             rpsi_x_acc_trajectory = deriveByTime(filtfilt(b_lowpass, a_lowpass, rpsi_x_vel_trajectory), sampling_rate_mocap^(-1));
            
            % initialize containers
            c7_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            c7_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            c7_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            lpsi_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            lpsi_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            lpsi_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            rpsi_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            rpsi_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            rpsi_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            lheel_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            lheel_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            lheel_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            rheel_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            rheel_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            rheel_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
%             lpsi_x_vel_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
%             rpsi_x_vel_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
%             pelvis_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
%             pelvis_x_vel_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
%             pelvis_x_acc_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            
            lleg_angle_ml_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            rleg_angle_ml_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            trunk_angle_ml_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        end
        if process_data_forceplate
            load(makeFileName(date, subject_id, 'walking', i_trial, 'forceplateTrajectories'));

            % initialize containers
            cop_x_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            lcop_x_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            rcop_x_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            fxl_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            fzl_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            myl_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            fxr_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            fzr_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            myr_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
        end
        if process_data_emg
            load(makeFileName(date, subject_id, 'walking', i_trial, 'emgTrajectories'));
            
            % rename relevant trajectories
            lglutmed_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LGLUTMED'));
            ltibiant_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LTIBIANT'));
            lgastroc_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LGASTROC'));
            lperolng_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LPEROLNG'));
            rglutmed_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RGLUTMED'));
            rtibiant_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RTIBIANT'));
            rgastroc_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RGASTROC'));
            rperolng_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RPEROLNG'));
            
            % initialize containers
            lglutmed_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            ltibiant_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            lgastroc_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            lperolng_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            rglutmed_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            rtibiant_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            rgastroc_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            rperolng_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            
        end
        if process_data_angles
            load(makeFileName(date, subject_id, 'walking', i_trial, 'angleTrajectories'));
            load(makeFileName(date, subject_id, 'walking', i_trial, 'kinematicTrajectories'));
            
            % initialize containers
            number_of_joints = plant.numberOfJoints;
            joint_angles_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial, number_of_joints);
            joint_velocities_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial, number_of_joints);
            joint_accelerations_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial, number_of_joints);
        end    
        if process_data_torques
            load(makeFileName(date, subject_id, 'walking', i_trial, 'inverseDynamics_hingeConstraints'));
        end    


        for i_stretch = 1 : number_of_stretches_trial

            % mocap data
            if process_data_marker
                % extract
                start_index_mocap = start_indices_mocap(i_stretch);
                end_index_mocap = end_indices_mocap(i_stretch);
                time_extracted_mocap = time_mocap(start_index_mocap : end_index_mocap);
                c7_x_pos_extracted_stretch = c7_x_pos_trajectory(start_index_mocap : end_index_mocap);
                c7_y_pos_extracted_stretch = c7_y_pos_trajectory(start_index_mocap : end_index_mocap);
                c7_z_pos_extracted_stretch = c7_z_pos_trajectory(start_index_mocap : end_index_mocap);
                lpsi_x_pos_extracted_stretch = lpsi_x_pos_trajectory(start_index_mocap : end_index_mocap);
                lpsi_y_pos_extracted_stretch = lpsi_y_pos_trajectory(start_index_mocap : end_index_mocap);
                lpsi_z_pos_extracted_stretch = lpsi_z_pos_trajectory(start_index_mocap : end_index_mocap);
                rpsi_x_pos_extracted_stretch = rpsi_x_pos_trajectory(start_index_mocap : end_index_mocap);
                rpsi_y_pos_extracted_stretch = rpsi_y_pos_trajectory(start_index_mocap : end_index_mocap);
                rpsi_z_pos_extracted_stretch = rpsi_z_pos_trajectory(start_index_mocap : end_index_mocap);
                lheel_x_pos_extracted_stretch = lheel_x_pos_trajectory(start_index_mocap : end_index_mocap);
                lheel_y_pos_extracted_stretch = lheel_y_pos_trajectory(start_index_mocap : end_index_mocap);
                lheel_z_pos_extracted_stretch = lheel_z_pos_trajectory(start_index_mocap : end_index_mocap);
                rheel_x_pos_extracted_stretch = rheel_x_pos_trajectory(start_index_mocap : end_index_mocap);
                rheel_y_pos_extracted_stretch = rheel_y_pos_trajectory(start_index_mocap : end_index_mocap);
                rheel_z_pos_extracted_stretch = rheel_z_pos_trajectory(start_index_mocap : end_index_mocap);
                
%                 lpsi_x_vel_extracted_stretch = lpsi_x_vel_trajectory(start_index_mocap : end_index_mocap);
%                 rpsi_x_vel_extracted_stretch = rpsi_x_vel_trajectory(start_index_mocap : end_index_mocap);
%                 lpsi_x_acc_extracted_stretch = lpsi_x_acc_trajectory(start_index_mocap : end_index_mocap);
%                 rpsi_x_acc_extracted_stretch = rpsi_x_acc_trajectory(start_index_mocap : end_index_mocap);
                
                
                % calculate vectors
                mpsi_x_pos = (lpsi_x_pos_extracted_stretch + rpsi_x_pos_extracted_stretch) * 0.5;
                mpsi_y_pos = (lpsi_y_pos_extracted_stretch + rpsi_y_pos_extracted_stretch) * 0.5;
                mpsi_z_pos = (lpsi_z_pos_extracted_stretch + rpsi_z_pos_extracted_stretch) * 0.5;
                lleg_vector_x = mpsi_x_pos - lheel_x_pos_extracted_stretch;
                lleg_vector_y = mpsi_y_pos - lheel_y_pos_extracted_stretch;
                lleg_vector_z = mpsi_z_pos - lheel_z_pos_extracted_stretch;
                rleg_vector_x = mpsi_x_pos - rheel_x_pos_extracted_stretch;
                rleg_vector_y = mpsi_y_pos - rheel_y_pos_extracted_stretch;
                rleg_vector_z = mpsi_z_pos - rheel_z_pos_extracted_stretch;
                trunk_vector_x = c7_x_pos_extracted_stretch - mpsi_x_pos;
                trunk_vector_y = c7_y_pos_extracted_stretch - mpsi_y_pos;
                trunk_vector_z = c7_z_pos_extracted_stretch - mpsi_z_pos;
                
                % calculate angles
                trunk_angle_ml_extracted_stretch = rad2deg(atan2(trunk_vector_x, trunk_vector_z));
                lleg_angle_ml_extracted_stretch = rad2deg(atan2(lleg_vector_x, lleg_vector_z));
                rleg_angle_ml_extracted_stretch = rad2deg(atan2(rleg_vector_x, rleg_vector_z));
                
                % define stance foot heel as spatial point of reference
                if strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_RIGHT')
                    stance_foot_heel_x_initial = rheel_x_pos_extracted_stretch(1);
                    stance_foot_heel_y_initial = rheel_y_pos_extracted_stretch(1);
                elseif strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_LEFT')
                    stance_foot_heel_x_initial = lheel_x_pos_extracted_stretch(1);
                    stance_foot_heel_y_initial = lheel_y_pos_extracted_stretch(1);
                end

                % normalize mocap data in time
                time_normalized_mocap = linspace(time_extracted_mocap(1), time_extracted_mocap(end), number_of_time_steps_normalized);
                lpsi_x_pos_normalized_stretch = spline(time_extracted_mocap, lpsi_x_pos_extracted_stretch, time_normalized_mocap);
                lheel_x_pos_normalized_stretch = spline(time_extracted_mocap, lheel_x_pos_extracted_stretch, time_normalized_mocap);
                lheel_y_pos_normalized_stretch = spline(time_extracted_mocap, lheel_y_pos_extracted_stretch, time_normalized_mocap);
                rpsi_x_pos_normalized_stretch = spline(time_extracted_mocap, rpsi_x_pos_extracted_stretch, time_normalized_mocap);
                rheel_x_pos_normalized_stretch = spline(time_extracted_mocap, rheel_x_pos_extracted_stretch, time_normalized_mocap);
                rheel_y_pos_normalized_stretch = spline(time_extracted_mocap, rheel_y_pos_extracted_stretch, time_normalized_mocap);
                
                trunk_angle_ml_normalized_stretch = spline(time_extracted_mocap, trunk_angle_ml_extracted_stretch, time_normalized_mocap);
                lleg_angle_ml_normalized_stretch = spline(time_extracted_mocap, lleg_angle_ml_extracted_stretch, time_normalized_mocap);
                rleg_angle_ml_normalized_stretch = spline(time_extracted_mocap, rleg_angle_ml_extracted_stretch, time_normalized_mocap);
                
%                 lpsi_x_vel_normalized_stretch = spline(time_extracted_mocap, lpsi_x_vel_extracted_stretch, time_normalized_mocap);
%                 rpsi_x_vel_normalized_stretch = spline(time_extracted_mocap, rpsi_x_vel_extracted_stretch, time_normalized_mocap);
%                 lpsi_x_acc_normalized_stretch = spline(time_extracted_mocap, lpsi_x_acc_extracted_stretch, time_normalized_mocap);
%                 rpsi_x_acc_normalized_stretch = spline(time_extracted_mocap, rpsi_x_acc_extracted_stretch, time_normalized_mocap);

                if use_stance_foot_as_reference
                    reference_x = stance_foot_heel_x_initial;
                    reference_y = stance_foot_heel_y_initial;
                else
                    reference_x = 0;
                    reference_y = 0;
                end

                % use stance foot heel as reference and store
                lpsi_x_pos_normalized_trial(:, i_stretch) = lpsi_x_pos_normalized_stretch - reference_x;
                rpsi_x_pos_normalized_trial(:, i_stretch) = rpsi_x_pos_normalized_stretch - reference_x;
                lheel_x_pos_normalized_trial(:, i_stretch) = lheel_x_pos_normalized_stretch - reference_x;
                rheel_x_pos_normalized_trial(:, i_stretch) = rheel_x_pos_normalized_stretch - reference_x;
                lheel_y_pos_normalized_trial(:, i_stretch) = lheel_y_pos_normalized_stretch - reference_y;
                rheel_y_pos_normalized_trial(:, i_stretch) = rheel_y_pos_normalized_stretch - reference_y;
                trunk_angle_ml_normalized_trial(:, i_stretch) = trunk_angle_ml_normalized_stretch;
                lleg_angle_ml_normalized_trial(:, i_stretch) = lleg_angle_ml_normalized_stretch;
                rleg_angle_ml_normalized_trial(:, i_stretch) = rleg_angle_ml_normalized_stretch;
                
%                 lpsi_x_vel_normalized_trial(:, i_stretch) = lpsi_x_vel_normalized_stretch;
%                 rpsi_x_vel_normalized_trial(:, i_stretch) = rpsi_x_vel_normalized_stretch;
%                 pelvis_x_pos_normalized_trial(:, i_stretch) = mean([lpsi_x_pos_normalized_stretch; rpsi_x_pos_normalized_stretch]) - stance_foot_heel_x_initial;
%                 pelvis_x_vel_normalized_trial(:, i_stretch) = mean([lpsi_x_vel_normalized_stretch; rpsi_x_vel_normalized_stretch]);
%                 pelvis_x_acc_normalized_trial(:, i_stretch) = mean([lpsi_x_acc_normalized_stretch; rpsi_x_acc_normalized_stretch]);

                % time and origin
                origin_trial_list_trial(i_stretch) = i_trial;
                origin_start_time_list_trial(i_stretch) = time_mocap(start_index_mocap);
                origin_end_time_list_trial(i_stretch) = time_mocap(end_index_mocap);
                step_times_trial(i_stretch) = time_mocap(end_index_mocap) - time_mocap(start_index_mocap);
            end

            % forceplate data
            if process_data_forceplate
                % define times
                start_index_forceplate = start_indices_forceplate(i_stretch);
                end_index_forceplate = end_indices_forceplate(i_stretch);
                
                % extract
                time_extracted_forceplate = time_forceplate(start_index_forceplate : end_index_forceplate);
                cop_x_extracted_stretch = total_forceplate_cop_world(start_index_forceplate : end_index_forceplate, 1);
                lcop_x_extracted_stretch = left_forceplate_cop_world(start_index_forceplate : end_index_forceplate, 1);
                rcop_x_extracted_stretch = right_forceplate_cop_world(start_index_forceplate : end_index_forceplate, 1);
                fxl_extracted_stretch = left_forceplate_wrench_world(start_index_forceplate : end_index_forceplate, 1);
                fzl_extracted_stretch = left_forceplate_wrench_world(start_index_forceplate : end_index_forceplate, 3);
                myl_extracted_stretch = left_forceplate_wrench_world(start_index_forceplate : end_index_forceplate, 5);
                fxr_extracted_stretch = right_forceplate_wrench_world(start_index_forceplate : end_index_forceplate, 1);
                fzr_extracted_stretch = right_forceplate_wrench_world(start_index_forceplate : end_index_forceplate, 3);
                myr_extracted_stretch = right_forceplate_wrench_world(start_index_forceplate : end_index_forceplate, 5);
                
                % normalize
                time_normalized_forceplate = linspace(time_extracted_forceplate(1), time_extracted_forceplate(end), number_of_time_steps_normalized);
                cop_x_normalized_stretch = spline(time_extracted_forceplate, cop_x_extracted_stretch, time_normalized_forceplate);
                lcop_x_normalized_stretch = spline(time_extracted_forceplate, lcop_x_extracted_stretch, time_normalized_forceplate);
                rcop_x_normalized_stretch = spline(time_extracted_forceplate, rcop_x_extracted_stretch, time_normalized_forceplate);
                fxl_normalized_stretch = spline(time_extracted_forceplate, fxl_extracted_stretch, time_normalized_forceplate);
                fzl_normalized_stretch = spline(time_extracted_forceplate, fzl_extracted_stretch, time_normalized_forceplate);
                myl_normalized_stretch = spline(time_extracted_forceplate, myl_extracted_stretch, time_normalized_forceplate);
                fxr_normalized_stretch = spline(time_extracted_forceplate, fxr_extracted_stretch, time_normalized_forceplate);
                fzr_normalized_stretch = spline(time_extracted_forceplate, fzr_extracted_stretch, time_normalized_forceplate);
                myr_normalized_stretch = spline(time_extracted_forceplate, myr_extracted_stretch, time_normalized_forceplate);

                % use stance foot heel as reference and store
                cop_x_normalized_trial(:, i_stretch) = cop_x_normalized_stretch - stance_foot_heel_x_initial;
                lcop_x_normalized_trial(:, i_stretch) = lcop_x_normalized_stretch - stance_foot_heel_x_initial;
                rcop_x_normalized_trial(:, i_stretch) = rcop_x_normalized_stretch - stance_foot_heel_x_initial;
                
                fxl_normalized_trial(:, i_stretch) = fxl_normalized_stretch;
                fzl_normalized_trial(:, i_stretch) = fzl_normalized_stretch;
                myl_normalized_trial(:, i_stretch) = myl_normalized_stretch;
                fxr_normalized_trial(:, i_stretch) = fxr_normalized_stretch;
                fzr_normalized_trial(:, i_stretch) = fzr_normalized_stretch;
                myr_normalized_trial(:, i_stretch) = myr_normalized_stretch;
            end    
            
            % emg data
            if process_data_emg
                % define times
                start_index_emg = start_indices_emg(i_stretch);
                end_index_emg = end_indices_emg(i_stretch);
                
                % extract
                time_extracted_emg = time_emg(start_index_emg : end_index_emg);
                lglutmed_extracted_stretch = lglutmed_trajectory(start_index_emg : end_index_emg, 1);
                ltibiant_extracted_stretch = ltibiant_trajectory(start_index_emg : end_index_emg, 1);
                lgastroc_extracted_stretch = lgastroc_trajectory(start_index_emg : end_index_emg, 1);
                lperolng_extracted_stretch = lperolng_trajectory(start_index_emg : end_index_emg, 1);
                rglutmed_extracted_stretch = rglutmed_trajectory(start_index_emg : end_index_emg, 1);
                rtibiant_extracted_stretch = rtibiant_trajectory(start_index_emg : end_index_emg, 1);
                rgastroc_extracted_stretch = rgastroc_trajectory(start_index_emg : end_index_emg, 1);
                rperolng_extracted_stretch = rperolng_trajectory(start_index_emg : end_index_emg, 1);
                
                % normalize
                time_normalized_emg = linspace(time_extracted_emg(1), time_extracted_emg(end), number_of_time_steps_normalized);
                lglutmed_normalized_stretch = spline(time_extracted_emg, lglutmed_extracted_stretch, time_normalized_emg);
                ltibiant_normalized_stretch = spline(time_extracted_emg, ltibiant_extracted_stretch, time_normalized_emg);
                lgastroc_normalized_stretch = spline(time_extracted_emg, lgastroc_extracted_stretch, time_normalized_emg);
                lperolng_normalized_stretch = spline(time_extracted_emg, lperolng_extracted_stretch, time_normalized_emg);
                rglutmed_normalized_stretch = spline(time_extracted_emg, rglutmed_extracted_stretch, time_normalized_emg);
                rtibiant_normalized_stretch = spline(time_extracted_emg, rtibiant_extracted_stretch, time_normalized_emg);
                rgastroc_normalized_stretch = spline(time_extracted_emg, rgastroc_extracted_stretch, time_normalized_emg);
                rperolng_normalized_stretch = spline(time_extracted_emg, rperolng_extracted_stretch, time_normalized_emg);
                
                % store
                lglutmed_normalized_trial(:, i_stretch) = lglutmed_normalized_stretch;
                ltibiant_normalized_trial(:, i_stretch) = ltibiant_normalized_stretch;
                lgastroc_normalized_trial(:, i_stretch) = lgastroc_normalized_stretch;
                lperolng_normalized_trial(:, i_stretch) = lperolng_normalized_stretch;
                rglutmed_normalized_trial(:, i_stretch) = rglutmed_normalized_stretch;
                rtibiant_normalized_trial(:, i_stretch) = rtibiant_normalized_stretch;
                rgastroc_normalized_trial(:, i_stretch) = rgastroc_normalized_stretch;
                rperolng_normalized_trial(:, i_stretch) = rperolng_normalized_stretch;
            end

            % angle data
            if process_data_angles
                % extract
                start_index_mocap = start_indices_mocap(i_stretch);
                end_index_mocap = end_indices_mocap(i_stretch);
                time_extracted_mocap = time_mocap(start_index_mocap : end_index_mocap);
                joint_angles_extracted_stretch = joint_angle_trajectories_belt(start_index_mocap : end_index_mocap, :);
                joint_velocities_extracted_stretch = joint_velocity_trajectories_belt(start_index_mocap : end_index_mocap, :);
                joint_accelerations_extracted_stretch = joint_acceleration_trajectories_belt(start_index_mocap : end_index_mocap, :);
                
                % normalize mocap data in time
                time_normalized_mocap = linspace(time_extracted_mocap(1), time_extracted_mocap(end), number_of_time_steps_normalized);
                joint_angles_normalized_stretch = zeros(length(time_normalized_mocap), number_of_joints);
                joint_velocities_normalized_stretch = zeros(length(time_normalized_mocap), number_of_joints);
                joint_accelerations_normalized_stretch = zeros(length(time_normalized_mocap), number_of_joints);
                for i_joint = 1 : number_of_joints
                    joint_angles_normalized_stretch(:, i_joint) = spline(time_extracted_mocap, joint_angles_extracted_stretch(:, i_joint), time_normalized_mocap);
                    joint_velocities_normalized_stretch(:, i_joint) = spline(time_extracted_mocap, joint_velocities_extracted_stretch(:, i_joint), time_normalized_mocap);
                    joint_accelerations_normalized_stretch(:, i_joint) = spline(time_extracted_mocap, joint_accelerations_extracted_stretch(:, i_joint), time_normalized_mocap);
                end
                
                % store
                joint_angles_normalized_trial(:, i_stretch, :) = joint_angles_normalized_stretch;
                joint_velocities_normalized_trial(:, i_stretch, :) = joint_velocities_normalized_stretch;
                joint_accelerations_normalized_trial(:, i_stretch, :) = joint_accelerations_normalized_stretch;
            end



        end

        % append trial containers to total containers
        condition_stance_foot_list_all = [condition_stance_foot_list_all; condition_stance_foot_list_trial];
        condition_perturbation_list_all = [condition_perturbation_list_all; condition_perturbation_list_trial];
        condition_delay_list_all = [condition_delay_list_all; condition_delay_list_trial];
        condition_index_list_all = [condition_index_list_all; condition_index_list_trial];
        origin_trial_list_all = [origin_trial_list_all; origin_trial_list_trial];
        origin_start_time_list_all = [origin_start_time_list_all; origin_start_time_list_trial];
        origin_end_time_list_all = [origin_end_time_list_all; origin_end_time_list_trial];
        step_times_all = [step_times_all; step_times_trial];
        stim_start_time_relative_to_stretch_all = [stim_start_time_relative_to_stretch_all; stim_start_time_relative_to_stretch_trial];
        if process_data_marker
            lpsi_x_pos_normalized_all = [lpsi_x_pos_normalized_all lpsi_x_pos_normalized_trial];
            rpsi_x_pos_normalized_all = [rpsi_x_pos_normalized_all rpsi_x_pos_normalized_trial];
            lheel_x_pos_normalized_all = [lheel_x_pos_normalized_all lheel_x_pos_normalized_trial];
            rheel_x_pos_normalized_all = [rheel_x_pos_normalized_all rheel_x_pos_normalized_trial];
            lheel_y_pos_normalized_all = [lheel_y_pos_normalized_all lheel_y_pos_normalized_trial];
            rheel_y_pos_normalized_all = [rheel_y_pos_normalized_all rheel_y_pos_normalized_trial];
            
            trunk_angle_ml_normalized_all = [trunk_angle_ml_normalized_all trunk_angle_ml_normalized_trial];
            lleg_angle_ml_normalized_all = [lleg_angle_ml_normalized_all lleg_angle_ml_normalized_trial];
            rleg_angle_ml_normalized_all = [rleg_angle_ml_normalized_all rleg_angle_ml_normalized_trial];
            
%             lpsi_x_vel_normalized_all = [lpsi_x_vel_normalized_all lpsi_x_vel_normalized_trial];
%             rpsi_x_vel_normalized_all = [rpsi_x_vel_normalized_all rpsi_x_vel_normalized_trial];
%             pelvis_x_pos_normalized_all = [pelvis_x_pos_normalized_all pelvis_x_pos_normalized_trial];
%             pelvis_x_vel_normalized_all = [pelvis_x_vel_normalized_all pelvis_x_vel_normalized_trial];
%             pelvis_x_acc_normalized_all = [pelvis_x_acc_normalized_all pelvis_x_acc_normalized_trial];
            
            
        end
        if process_data_forceplate
            cop_x_normalized_all = [cop_x_normalized_all cop_x_normalized_trial];
            lcop_x_normalized_all = [lcop_x_normalized_all lcop_x_normalized_trial];
            rcop_x_normalized_all = [rcop_x_normalized_all rcop_x_normalized_trial];
            fxl_normalized_all = [fxl_normalized_all fxl_normalized_trial];
            fzl_normalized_all = [fzl_normalized_all fzl_normalized_trial];
            myl_normalized_all = [myl_normalized_all myl_normalized_trial];
            fxr_normalized_all = [fxr_normalized_all fxr_normalized_trial];
            fzr_normalized_all = [fzr_normalized_all fzr_normalized_trial];
            myr_normalized_all = [myr_normalized_all myr_normalized_trial];
        end
        if process_data_emg
            lglutmed_normalized_all = [lglutmed_normalized_all lglutmed_normalized_trial];
            ltibiant_normalized_all = [ltibiant_normalized_all ltibiant_normalized_trial];
            lgastroc_normalized_all = [lgastroc_normalized_all lgastroc_normalized_trial];
            lperolng_normalized_all = [lperolng_normalized_all lperolng_normalized_trial];
            rglutmed_normalized_all = [rglutmed_normalized_all rglutmed_normalized_trial];
            rtibiant_normalized_all = [rtibiant_normalized_all rtibiant_normalized_trial];
            rgastroc_normalized_all = [rgastroc_normalized_all rgastroc_normalized_trial];
            rperolng_normalized_all = [rperolng_normalized_all rperolng_normalized_trial];
        end
        if process_data_angles
            joint_angles_normalized_all = [joint_angles_normalized_all joint_angles_normalized_trial];
            joint_velocities_normalized_all = [joint_velocities_normalized_all joint_velocities_normalized_trial];
            joint_accelerations_normalized_all = [joint_accelerations_normalized_all joint_accelerations_normalized_trial];
        end

        disp(['Trial ' num2str(i_trial) ' extracted']);
    end
    
    step_time_mean = mean(step_times_all);
    time_normalized = linspace(0, step_time_mean, number_of_time_steps_normalized);
    number_of_stretches = length(step_times_all);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% old stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% calculate responses
if calculate_responses
    
    if process_data_marker
        % calculate control means
        lheel_x_pos_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rheel_x_pos_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lheel_y_pos_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rheel_y_pos_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        trunk_angle_ml_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lleg_angle_ml_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rleg_angle_ml_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        
%         pelvis_x_pos_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         pelvis_x_vel_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         pelvis_x_acc_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            lheel_x_pos_control_means(:, i_condition) = mean(lheel_x_pos_normalized_all(:, condition_indicator), 2);
            rheel_x_pos_control_means(:, i_condition) = mean(rheel_x_pos_normalized_all(:, condition_indicator), 2);
            lheel_y_pos_control_means(:, i_condition) = mean(lheel_y_pos_normalized_all(:, condition_indicator), 2);
            rheel_y_pos_control_means(:, i_condition) = mean(rheel_y_pos_normalized_all(:, condition_indicator), 2);
            trunk_angle_ml_control_means(:, i_condition) = mean(trunk_angle_ml_normalized_all(:, condition_indicator), 2);
            lleg_angle_ml_control_means(:, i_condition) = mean(lleg_angle_ml_normalized_all(:, condition_indicator), 2);
            rleg_angle_ml_control_means(:, i_condition) = mean(rleg_angle_ml_normalized_all(:, condition_indicator), 2);
%             pelvis_x_pos_control_means(:, i_condition) = mean(pelvis_x_pos_normalized_all(:, condition_indicator), 2);
%             pelvis_x_vel_control_means(:, i_condition) = mean(pelvis_x_vel_normalized_all(:, condition_indicator), 2);
%             pelvis_x_acc_control_means(:, i_condition) = mean(pelvis_x_acc_normalized_all(:, condition_indicator), 2);
        end
        
        % calculate stimulus responses
        lheel_x_pos_response = zeros(size(lheel_x_pos_normalized_all));
        rheel_x_pos_response = zeros(size(rheel_x_pos_normalized_all));
        lheel_y_pos_response = zeros(size(lheel_y_pos_normalized_all));
        rheel_y_pos_response = zeros(size(rheel_y_pos_normalized_all));
        trunk_angle_ml_response = zeros(size(trunk_angle_ml_normalized_all));
        lleg_angle_ml_response = zeros(size(lleg_angle_ml_normalized_all));
        rleg_angle_ml_response = zeros(size(rleg_angle_ml_normalized_all));
        
%         pelvis_x_pos_response = zeros(size(pelvis_x_pos_normalized_all));
%         pelvis_x_vel_response = zeros(size(pelvis_x_vel_normalized_all));
%         pelvis_x_acc_response = zeros(size(pelvis_x_acc_normalized_all));
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            lheel_x_pos_response(:, condition_indicator) = lheel_x_pos_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_control_means(:, i_condition), 1, sum(condition_indicator));
            rheel_x_pos_response(:, condition_indicator) = rheel_x_pos_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_control_means(:, i_condition), 1, sum(condition_indicator));
            lheel_y_pos_response(:, condition_indicator) = lheel_y_pos_normalized_all(:, condition_indicator) - repmat(lheel_y_pos_control_means(:, i_condition), 1, sum(condition_indicator));
            rheel_y_pos_response(:, condition_indicator) = rheel_y_pos_normalized_all(:, condition_indicator) - repmat(rheel_y_pos_control_means(:, i_condition), 1, sum(condition_indicator));
            
            trunk_angle_ml_response(:, condition_indicator) = trunk_angle_ml_normalized_all(:, condition_indicator) - repmat(trunk_angle_ml_control_means(:, i_condition), 1, sum(condition_indicator));
            lleg_angle_ml_response(:, condition_indicator) = lleg_angle_ml_normalized_all(:, condition_indicator) - repmat(lleg_angle_ml_control_means(:, i_condition), 1, sum(condition_indicator));
            rleg_angle_ml_response(:, condition_indicator) = rleg_angle_ml_normalized_all(:, condition_indicator) - repmat(rleg_angle_ml_control_means(:, i_condition), 1, sum(condition_indicator));
%             pelvis_x_pos_response(:, condition_indicator) = pelvis_x_pos_normalized_all(:, condition_indicator) - repmat(pelvis_x_pos_control_means(:, i_condition), 1, sum(condition_indicator));
%             pelvis_x_vel_response(:, condition_indicator) = pelvis_x_vel_normalized_all(:, condition_indicator) - repmat(pelvis_x_vel_control_means(:, i_condition), 1, sum(condition_indicator));
%             pelvis_x_acc_response(:, condition_indicator) = pelvis_x_acc_normalized_all(:, condition_indicator) - repmat(pelvis_x_acc_control_means(:, i_condition), 1, sum(condition_indicator));
        end
        for i_condition = 1 : number_of_conditions_to_analyze
            condition_indicator = conditions_to_analyze_indicators(:, i_condition);
            lheel_x_pos_response(:, condition_indicator) = lheel_x_pos_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rheel_x_pos_response(:, condition_indicator) = rheel_x_pos_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lheel_y_pos_response(:, condition_indicator) = lheel_y_pos_normalized_all(:, condition_indicator) - repmat(lheel_y_pos_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rheel_y_pos_response(:, condition_indicator) = rheel_y_pos_normalized_all(:, condition_indicator) - repmat(rheel_y_pos_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            trunk_angle_ml_response(:, condition_indicator) = trunk_angle_ml_normalized_all(:, condition_indicator) - repmat(trunk_angle_ml_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lleg_angle_ml_response(:, condition_indicator) = lleg_angle_ml_normalized_all(:, condition_indicator) - repmat(lleg_angle_ml_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rleg_angle_ml_response(:, condition_indicator) = rleg_angle_ml_normalized_all(:, condition_indicator) - repmat(rleg_angle_ml_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             pelvis_x_pos_response(:, condition_indicator) = pelvis_x_pos_normalized_all(:, condition_indicator) - repmat(pelvis_x_pos_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             pelvis_x_vel_response(:, condition_indicator) = pelvis_x_vel_normalized_all(:, condition_indicator) - repmat(pelvis_x_vel_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             pelvis_x_acc_response(:, condition_indicator) = pelvis_x_acc_normalized_all(:, condition_indicator) - repmat(pelvis_x_acc_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
        end
        
% the following is still relevant, but I'm not quite sure how to automate it 
%{
        % estimate stimulus response
        step_response_all = zeros(number_of_stretches, 1);
        stim_response_all = zeros(number_of_stretches, 1);
        [ ...
          Jacobian_stanceR, ...
          correlation_c_stanceR, ...
          correlation_p_stanceR, ...
          step_response_all(conditions_stanceR_stimPos_0ms), ...
          stim_response_all(conditions_stanceR_stimPos_0ms) ...
        ] = ...
        calculateStepResponse ...
        ( ...
          lheel_x_pos_normalized_all(:, conditions_stanceR_stimPos_0ms), ...
          lheel_x_pos_normalized_all(:, conditions_stanceR_stimNo), ...
          pelvis_x_pos_normalized_all(:, conditions_stanceR_stimPos_0ms), ...
          pelvis_x_pos_normalized_all(:, conditions_stanceR_stimNo), ...
          pelvis_x_vel_normalized_all(:, conditions_stanceR_stimPos_0ms), ...
          pelvis_x_vel_normalized_all(:, conditions_stanceR_stimNo), ...
          round(number_of_time_steps_normalized/2), ...
          number_of_time_steps_normalized ...
        );
        [ ...
          ~, ...
          ~, ...
          ~, ...
          step_response_all(conditions_stanceR_stimNeg_0ms), ...
          stim_response_all(conditions_stanceR_stimNeg_0ms) ...
        ] = ...
        calculateStepResponse ...
        ( ...
          lheel_x_pos_normalized_all(:, conditions_stanceR_stimNeg_0ms), ...
          lheel_x_pos_normalized_all(:, conditions_stanceR_stimNo), ...
          pelvis_x_pos_normalized_all(:, conditions_stanceR_stimNeg_0ms), ...
          pelvis_x_pos_normalized_all(:, conditions_stanceR_stimNo), ...
          pelvis_x_vel_normalized_all(:, conditions_stanceR_stimNeg_0ms), ...
          pelvis_x_vel_normalized_all(:, conditions_stanceR_stimNo), ...
          round(number_of_time_steps_normalized/2), ...
          number_of_time_steps_normalized ...
        );        
    
        [ ...
          ~, ...
          ~, ...
          ~, ...
          step_response_all(conditions_stanceR_stimPos_150ms), ...
          stim_response_all(conditions_stanceR_stimPos_150ms) ...
        ] = ...
        calculateStepResponse ...
        ( ...
          lheel_x_pos_normalized_all(:, conditions_stanceR_stimPos_150ms), ...
          lheel_x_pos_normalized_all(:, conditions_stanceR_stimNo), ...
          pelvis_x_pos_normalized_all(:, conditions_stanceR_stimPos_150ms), ...
          pelvis_x_pos_normalized_all(:, conditions_stanceR_stimNo), ...
          pelvis_x_vel_normalized_all(:, conditions_stanceR_stimPos_150ms), ...
          pelvis_x_vel_normalized_all(:, conditions_stanceR_stimNo), ...
          round(number_of_time_steps_normalized/2), ...
          number_of_time_steps_normalized ...
        );
        [ ...
          ~, ...
          ~, ...
          ~, ...
          step_response_all(conditions_stanceR_stimNeg_150ms), ...
          stim_response_all(conditions_stanceR_stimNeg_150ms) ...
        ] = ...
        calculateStepResponse ...
        ( ...
          lheel_x_pos_normalized_all(:, conditions_stanceR_stimNeg_150ms), ...
          lheel_x_pos_normalized_all(:, conditions_stanceR_stimNo), ...
          pelvis_x_pos_normalized_all(:, conditions_stanceR_stimNeg_150ms), ...
          pelvis_x_pos_normalized_all(:, conditions_stanceR_stimNo), ...
          pelvis_x_vel_normalized_all(:, conditions_stanceR_stimNeg_150ms), ...
          pelvis_x_vel_normalized_all(:, conditions_stanceR_stimNo), ...
          round(number_of_time_steps_normalized/2), ...
          number_of_time_steps_normalized ...
        );
    
        [ ...
          Jacobian_stanceL, ...
          correlation_c_stanceL, ...
          correlation_p_stanceL, ...
          step_response_all(conditions_stanceL_stimPos_450ms), ...
          stim_response_all(conditions_stanceL_stimPos_450ms) ...
        ] = ...
        calculateStepResponse ...
        ( ...
          rheel_x_pos_normalized_all(:, conditions_stanceL_stimPos_450ms), ...
          rheel_x_pos_normalized_all(:, conditions_stanceL_stimNo), ...
          pelvis_x_pos_normalized_all(:, conditions_stanceL_stimPos_450ms), ...
          pelvis_x_pos_normalized_all(:, conditions_stanceL_stimNo), ...
          pelvis_x_vel_normalized_all(:, conditions_stanceL_stimPos_450ms), ...
          pelvis_x_vel_normalized_all(:, conditions_stanceL_stimNo), ...
          round(number_of_time_steps_normalized/2), ...
          number_of_time_steps_normalized ...
        );
        [ ...
          ~, ...
          ~, ...
          ~, ...
          step_response_all(conditions_stanceL_stimNeg_450ms), ...
          stim_response_all(conditions_stanceL_stimNeg_450ms) ...
        ] = ...
        calculateStepResponse ...
        ( ...
          rheel_x_pos_normalized_all(:, conditions_stanceL_stimNeg_450ms), ...
          rheel_x_pos_normalized_all(:, conditions_stanceL_stimNo), ...
          pelvis_x_pos_normalized_all(:, conditions_stanceL_stimNeg_450ms), ...
          pelvis_x_pos_normalized_all(:, conditions_stanceL_stimNo), ...
          pelvis_x_vel_normalized_all(:, conditions_stanceL_stimNeg_450ms), ...
          pelvis_x_vel_normalized_all(:, conditions_stanceL_stimNo), ...
          round(number_of_time_steps_normalized/2), ...
          number_of_time_steps_normalized ...
        );        
%}
    end
    
    if process_data_forceplate
        % calculate control means
        cop_x_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lcop_x_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rcop_x_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        fxl_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        fzl_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        myl_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        fxr_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        fzr_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        myr_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            cop_x_control_means(:, i_condition) = mean(cop_x_normalized_all(:, condition_indicator), 2);
            lcop_x_control_means(:, i_condition) = mean(lcop_x_normalized_all(:, condition_indicator), 2);
            rcop_x_control_means(:, i_condition) = mean(rcop_x_normalized_all(:, condition_indicator), 2);
            fxl_control_means(:, i_condition) = mean(fxl_normalized_all(:, condition_indicator), 2);
            fzl_control_means(:, i_condition) = mean(fzl_normalized_all(:, condition_indicator), 2);
            myl_control_means(:, i_condition) = mean(myl_normalized_all(:, condition_indicator), 2);
            fxr_control_means(:, i_condition) = mean(fxr_normalized_all(:, condition_indicator), 2);
            fzr_control_means(:, i_condition) = mean(fzr_normalized_all(:, condition_indicator), 2);
            myr_control_means(:, i_condition) = mean(myr_normalized_all(:, condition_indicator), 2);
        end
        
        % calculate stimulus responses
        cop_x_response = zeros(size(cop_x_normalized_all));
        lcop_x_response = zeros(size(lcop_x_normalized_all));
        rcop_x_response = zeros(size(lcop_x_normalized_all));
        fxl_response = zeros(size(fxl_normalized_all));
        fzl_response = zeros(size(fzl_normalized_all));
        myl_response = zeros(size(myl_normalized_all));
        fxr_response = zeros(size(fxl_normalized_all));
        fzr_response = zeros(size(fzl_normalized_all));
        myr_response = zeros(size(myl_normalized_all));
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            cop_x_response(:, condition_indicator) = cop_x_normalized_all(:, condition_indicator) - repmat(cop_x_control_means(:, i_condition), 1, sum(condition_indicator));
            fxl_response(:, condition_indicator) = fxl_normalized_all(:, condition_indicator) - repmat(fxl_control_means(:, i_condition), 1, sum(condition_indicator));
            fzl_response(:, condition_indicator) = fzl_normalized_all(:, condition_indicator) - repmat(fzl_control_means(:, i_condition), 1, sum(condition_indicator));
            myl_response(:, condition_indicator) = myl_normalized_all(:, condition_indicator) - repmat(myl_control_means(:, i_condition), 1, sum(condition_indicator));
            lcop_x_response(:, condition_indicator) = lcop_x_normalized_all(:, condition_indicator) - repmat(lcop_x_control_means(:, i_condition), 1, sum(condition_indicator));
            fxr_response(:, condition_indicator) = fxr_normalized_all(:, condition_indicator) - repmat(fxr_control_means(:, i_condition), 1, sum(condition_indicator));
            fzr_response(:, condition_indicator) = fzr_normalized_all(:, condition_indicator) - repmat(fzr_control_means(:, i_condition), 1, sum(condition_indicator));
            myr_response(:, condition_indicator) = myr_normalized_all(:, condition_indicator) - repmat(myr_control_means(:, i_condition), 1, sum(condition_indicator));
            rcop_x_response(:, condition_indicator) = rcop_x_normalized_all(:, condition_indicator) - repmat(rcop_x_control_means(:, i_condition), 1, sum(condition_indicator));
        end
        for i_condition = 1 : number_of_conditions_to_analyze
            condition_indicator = conditions_to_analyze_indicators(:, i_condition);
            cop_x_response(:, condition_indicator) = cop_x_normalized_all(:, condition_indicator) - repmat(cop_x_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            fxl_response(:, condition_indicator) = fxl_normalized_all(:, condition_indicator) - repmat(fxl_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            fzl_response(:, condition_indicator) = fzl_normalized_all(:, condition_indicator) - repmat(fzl_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            myl_response(:, condition_indicator) = myl_normalized_all(:, condition_indicator) - repmat(myl_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lcop_x_response(:, condition_indicator) = lcop_x_normalized_all(:, condition_indicator) - repmat(lcop_x_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            fxr_response(:, condition_indicator) = fxr_normalized_all(:, condition_indicator) - repmat(fxr_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            fzr_response(:, condition_indicator) = fzr_normalized_all(:, condition_indicator) - repmat(fzr_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            myr_response(:, condition_indicator) = myr_normalized_all(:, condition_indicator) - repmat(myr_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rcop_x_response(:, condition_indicator) = rcop_x_normalized_all(:, condition_indicator) - repmat(rcop_x_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
        end
        
%         % control means
%         fxl_mean_stanceL_control = mean(fxl_normalized_all(:, conditions_stanceL_stimNo), 2);
%         fzl_mean_stanceL_control = mean(fzl_normalized_all(:, conditions_stanceL_stimNo), 2);
%         myl_mean_stanceL_control = mean(myl_normalized_all(:, conditions_stanceL_stimNo), 2);
%         lcop_x_mean_stanceL_control = mean(lcop_x_normalized_all(:, conditions_stanceL_stimNo), 2);
%         fxr_mean_stanceR_control = mean(fxr_normalized_all(:, conditions_stanceR_stimNo), 2);
%         fzr_mean_stanceR_control = mean(fzr_normalized_all(:, conditions_stanceR_stimNo), 2);
%         myr_mean_stanceR_control = mean(myr_normalized_all(:, conditions_stanceR_stimNo), 2);
%         rcop_x_mean_stanceR_control = mean(rcop_x_normalized_all(:, conditions_stanceR_stimNo), 2);
%         
%         % response
%         fxl_stanceL_response = fxl_normalized_all - repmat(fxl_mean_stanceL_control, 1, number_of_stretches);
%         fzl_stanceL_response = fzl_normalized_all - repmat(fzl_mean_stanceL_control, 1, number_of_stretches);
%         myl_stanceL_response = myl_normalized_all - repmat(myl_mean_stanceL_control, 1, number_of_stretches);
%         lcop_x_stanceL_response = lcop_x_normalized_all - repmat(lcop_x_mean_stanceL_control, 1, number_of_stretches);
%         fxr_stanceR_response = fxr_normalized_all - repmat(fxr_mean_stanceR_control, 1, number_of_stretches);
%         fzr_stanceR_response = fzr_normalized_all - repmat(fzr_mean_stanceR_control, 1, number_of_stretches);
%         myr_stanceR_response = myr_normalized_all - repmat(myr_mean_stanceR_control, 1, number_of_stretches);
%         rcop_x_stanceR_response = rcop_x_normalized_all - repmat(rcop_x_mean_stanceR_control, 1, number_of_stretches);
    end
    
    if process_data_emg
        % calculate control means
        lglutmed_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        ltibiant_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lgastroc_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lperolng_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rglutmed_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rtibiant_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rgastroc_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rperolng_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            lglutmed_control_means(:, i_condition) = mean(lglutmed_normalized_all(:, condition_indicator), 2);
            ltibiant_control_means(:, i_condition) = mean(ltibiant_normalized_all(:, condition_indicator), 2);
            lgastroc_control_means(:, i_condition) = mean(lgastroc_normalized_all(:, condition_indicator), 2);
            lperolng_control_means(:, i_condition) = mean(lperolng_normalized_all(:, condition_indicator), 2);
            rglutmed_control_means(:, i_condition) = mean(rglutmed_normalized_all(:, condition_indicator), 2);
            rtibiant_control_means(:, i_condition) = mean(rtibiant_normalized_all(:, condition_indicator), 2);
            rgastroc_control_means(:, i_condition) = mean(rgastroc_normalized_all(:, condition_indicator), 2);
            rperolng_control_means(:, i_condition) = mean(rperolng_normalized_all(:, condition_indicator), 2);
        end
        
        % calculate stimulus responses
        lglutmed_response = zeros(size(lglutmed_normalized_all));
        ltibiant_response = zeros(size(ltibiant_normalized_all));
        lgastroc_response = zeros(size(lgastroc_normalized_all));
        lperolng_response = zeros(size(lperolng_normalized_all));
        rglutmed_response = zeros(size(rglutmed_normalized_all));
        rtibiant_response = zeros(size(rtibiant_normalized_all));
        rgastroc_response = zeros(size(rgastroc_normalized_all));
        rperolng_response = zeros(size(rperolng_normalized_all));
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            lglutmed_response(:, condition_indicator) = lglutmed_normalized_all(:, condition_indicator) - repmat(lglutmed_control_means(:, i_condition), 1, sum(condition_indicator));
            ltibiant_response(:, condition_indicator) = ltibiant_normalized_all(:, condition_indicator) - repmat(ltibiant_control_means(:, i_condition), 1, sum(condition_indicator));
            lgastroc_response(:, condition_indicator) = lgastroc_normalized_all(:, condition_indicator) - repmat(lgastroc_control_means(:, i_condition), 1, sum(condition_indicator));
            lperolng_response(:, condition_indicator) = lperolng_normalized_all(:, condition_indicator) - repmat(lperolng_control_means(:, i_condition), 1, sum(condition_indicator));
            rglutmed_response(:, condition_indicator) = rglutmed_normalized_all(:, condition_indicator) - repmat(rglutmed_control_means(:, i_condition), 1, sum(condition_indicator));
            rtibiant_response(:, condition_indicator) = rtibiant_normalized_all(:, condition_indicator) - repmat(rtibiant_control_means(:, i_condition), 1, sum(condition_indicator));
            rgastroc_response(:, condition_indicator) = rgastroc_normalized_all(:, condition_indicator) - repmat(rgastroc_control_means(:, i_condition), 1, sum(condition_indicator));
            rperolng_response(:, condition_indicator) = rperolng_normalized_all(:, condition_indicator) - repmat(rperolng_control_means(:, i_condition), 1, sum(condition_indicator));
        end
        for i_condition = 1 : number_of_conditions_to_analyze
            condition_indicator = conditions_to_analyze_indicators(:, i_condition);
            lglutmed_response(:, condition_indicator) = lglutmed_normalized_all(:, condition_indicator) - repmat(lglutmed_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            ltibiant_response(:, condition_indicator) = ltibiant_normalized_all(:, condition_indicator) - repmat(ltibiant_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lgastroc_response(:, condition_indicator) = lgastroc_normalized_all(:, condition_indicator) - repmat(lgastroc_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lperolng_response(:, condition_indicator) = lperolng_normalized_all(:, condition_indicator) - repmat(lperolng_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rglutmed_response(:, condition_indicator) = rglutmed_normalized_all(:, condition_indicator) - repmat(rglutmed_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rtibiant_response(:, condition_indicator) = rtibiant_normalized_all(:, condition_indicator) - repmat(rtibiant_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rgastroc_response(:, condition_indicator) = rgastroc_normalized_all(:, condition_indicator) - repmat(rgastroc_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rperolng_response(:, condition_indicator) = rperolng_normalized_all(:, condition_indicator) - repmat(rperolng_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
        end        
    end
        
    if process_data_angles
        % control means
        joint_angles_mean_stanceL_control = squeeze(mean(joint_angles_normalized_all(:, conditions_stanceL_stimNo, :), 2));
        joint_angles_mean_stanceR_control = squeeze(mean(joint_angles_normalized_all(:, conditions_stanceR_stimNo, :), 2));
        joint_velocities_mean_stanceL_control = squeeze(mean(joint_velocities_normalized_all(:, conditions_stanceL_stimNo, :), 2));
        joint_velocities_mean_stanceR_control = squeeze(mean(joint_velocities_normalized_all(:, conditions_stanceR_stimNo, :), 2));
        joint_accelerations_mean_stanceL_control = squeeze(mean(joint_accelerations_normalized_all(:, conditions_stanceL_stimNo, :), 2));
        joint_accelerations_mean_stanceR_control = squeeze(mean(joint_accelerations_normalized_all(:, conditions_stanceR_stimNo, :), 2));
        
        % response
        joint_angles_stanceL_response = joint_angles_normalized_all - permute(repmat(joint_angles_mean_stanceL_control, 1, 1, number_of_stretches), [1 3 2]);
        joint_angles_stanceR_response = joint_angles_normalized_all - permute(repmat(joint_angles_mean_stanceR_control, 1, 1, number_of_stretches), [1 3 2]);
        joint_velocities_stanceL_response = joint_velocities_normalized_all - permute(repmat(joint_velocities_mean_stanceL_control, 1, 1, number_of_stretches), [1 3 2]);
        joint_velocities_stanceR_response = joint_velocities_normalized_all - permute(repmat(joint_velocities_mean_stanceR_control, 1, 1, number_of_stretches), [1 3 2]);
        joint_accelerations_stanceL_response = joint_accelerations_normalized_all - permute(repmat(joint_accelerations_mean_stanceL_control, 1, 1, number_of_stretches), [1 3 2]);
        joint_accelerations_stanceR_response = joint_accelerations_normalized_all - permute(repmat(joint_accelerations_mean_stanceR_control, 1, 1, number_of_stretches), [1 3 2]);
        
%         i_joint = 11;
%         plot(squeeze(joint_angles_normalized_all(:, conditions_stanceL_control, i_joint)))
%         plot(squeeze(joint_angles_normalized_all(:, conditions_stanceR_control, i_joint)))
%         plot(squeeze(joint_velocities_normalized_all(:, conditions_stanceL_control, i_joint)))
%         plot(squeeze(joint_velocities_normalized_all(:, conditions_stanceR_control, i_joint)))
%         plot(squeeze(joint_accelerations_normalized_all(:, conditions_stanceL_control, i_joint)))
%         plot(squeeze(joint_accelerations_normalized_all(:, conditions_stanceR_control, i_joint)))
  
%         i_joint = 11;
%         plot(mean(squeeze(joint_angles_stanceL_response(:, conditions_stanceL_control, i_joint)), 2))
%         plot(mean(squeeze(joint_angles_right_response(:, conditions_stanceR_control, i_joint)), 2))
%         plot(mean(squeeze(joint_velocities_stanceL_response(:, conditions_stanceL_control, i_joint)), 2))
%         plot(mean(squeeze(joint_velocities_right_response(:, conditions_stanceR_control, i_joint)), 2))
%         plot(mean(squeeze(joint_accelerations_stanceL_response(:, conditions_stanceL_control, i_joint)), 2))
%         plot(mean(squeeze(joint_accelerations_right_response(:, conditions_stanceR_control, i_joint)), 2))

    end
end

%% save data
if save_data
    % save conditions 
    save ...
      ( ...
        makeFileName(date, subject_id, 'resultsConditions'), ...
        'time_normalized', ...
        'origin_trial_list_all', ...
        'origin_start_time_list_all', ...
        'origin_end_time_list_all', ...
        'condition_stance_foot_list_all', ...
        'condition_perturbation_list_all', ...
        'condition_delay_list_all', ...
        'condition_index_list_all', ...
        'conditions_control_indicators', ...
        'conditions_to_analyze_indicators' ...
      );
  
    if process_data_balance
        save ...
          ( ...
            makeFileName(date, subject_id, 'resultsMarker'), ...
            'lheel_x_pos_normalized_all', ...
            'rheel_x_pos_normalized_all', ...
            'lheel_y_pos_normalized_all', ...
            'rheel_y_pos_normalized_all', ...
            'trunk_angle_ml_normalized_all', ...
            'lleg_angle_ml_normalized_all', ...
            'rleg_angle_ml_normalized_all', ...
            'lheel_x_pos_response', ...
            'rheel_x_pos_response', ...
            'lheel_y_pos_response', ...
            'rheel_y_pos_response', ...
            'trunk_angle_ml_response', ...
            'lleg_angle_ml_response', ...
            'rleg_angle_ml_response' ...
          );
%             'pelvis_x_pos_normalized_all', ...
%             'pelvis_x_vel_normalized_all', ...
%             'pelvis_x_acc_normalized_all', ...
%             'pelvis_x_pos_response', ...
%             'pelvis_x_vel_response', ...
%             'pelvis_x_acc_response' ...
%             'lheel_x_pos_control_means', ...
%             'rheel_x_pos_control_means', ...
%             'lheel_y_pos_control_means', ...
%             'rheel_y_pos_control_means', ...
%             'pelvis_x_pos_control_means', ...
%             'pelvis_x_vel_control_means', ...
%             'pelvis_x_acc_control_means', ...
    end  
    
    if process_data_forceplate
        save ...
          ( ...
            makeFileName(date, subject_id, 'resultsForceplate'), ...
            'cop_x_normalized_all', ...
            'fxl_normalized_all', ...
            'fzl_normalized_all', ...
            'myl_normalized_all', ...
            'lcop_x_normalized_all', ...
            'fxr_normalized_all', ...
            'fzr_normalized_all', ...
            'myr_normalized_all', ...
            'rcop_x_normalized_all', ...
            'cop_x_response', ...
            'fxl_response', ...
            'fzl_response', ...
            'myl_response', ...
            'lcop_x_response', ...
            'fxr_response', ...
            'fzr_response', ...
            'myr_response', ...
            'rcop_x_response' ...
          );
    end
    
    if process_data_emg
        save ...
          ( ...
            makeFileName(date, subject_id, 'resultsEmg'), ...
            'lglutmed_normalized_all', ...
            'ltibiant_normalized_all', ...
            'lgastroc_normalized_all', ...
            'lperolng_normalized_all', ...
            'rglutmed_normalized_all', ...
            'rtibiant_normalized_all', ...
            'rgastroc_normalized_all', ...
            'rperolng_normalized_all', ...
            'lglutmed_response', ...
            'ltibiant_response', ...
            'lgastroc_response', ...
            'lperolng_response', ...
            'rglutmed_response', ...
            'rtibiant_response', ...
            'rgastroc_response', ...
            'rperolng_response' ...
          )
    end    
end
    


%% calculate_strategy_directions
if calculate_strategy_directions
    % define references
    plant.jointAngles = zeros(plant.numberOfJoints, 1);
    plant.jointVelocities = zeros(plant.numberOfJoints, 1);
    plant.updateInternals;
    left_ankle_scs_to_world_rotation_reference = plant.endEffectorTransformations{3}(1:3, 1:3);
    right_ankle_scs_to_world_rotation_reference = plant.endEffectorTransformations{6}(1:3, 1:3);
    
    % left
    left_ankle_strategy = zeros(number_of_time_steps_normalized, number_of_joints);
    left_hip_strategy = zeros(number_of_time_steps_normalized, number_of_joints);
    left_ankle_acceleration_strategy = zeros(number_of_time_steps_normalized, number_of_joints);
    left_hip_acceleration_strategy = zeros(number_of_time_steps_normalized, number_of_joints);

    phi_left_trajectory_left = zeros(number_of_time_steps_normalized, 1);
    rho_left_trajectory_left = zeros(number_of_time_steps_normalized, 1);
    phi_right_trajectory_left = zeros(number_of_time_steps_normalized, 1);
    rho_right_trajectory_left = zeros(number_of_time_steps_normalized, 1);
    phi_left_trajectory_right = zeros(number_of_time_steps_normalized, 1);
    rho_left_trajectory_right = zeros(number_of_time_steps_normalized, 1);
    phi_right_trajectory_right = zeros(number_of_time_steps_normalized, 1);
    rho_right_trajectory_right = zeros(number_of_time_steps_normalized, 1);
    
    A_left_trajectory = cell(number_of_time_steps_normalized, 1);
    M_left_trajectory = cell(number_of_time_steps_normalized, 1);
    P_left_trajectory = cell(number_of_time_steps_normalized, 1);
    J_c_left_trajectory = cell(number_of_time_steps_normalized, 1);
    
    A_right_trajectory = cell(number_of_time_steps_normalized, 1);
    
    disp('Calculating ankle and hip strategy directions for left stance foot');
    tic
    for i_time = 1 : number_of_time_steps_normalized

        % update plant
        plant.jointAngles = joint_angles_mean_left_control(i_time, :)';
        plant.updateInternals;

        M = plant.inertiaMatrix;
        C = plant.coriolisMatrix;
        N = plant.gravitationalTorqueMatrix;
        
        left_ankle_scs_transformation_current = plant.endEffectorTransformations{3};
        left_ankle_scs_to_world_rotation_current = left_ankle_scs_transformation_current(1:3, 1:3);
        left_ankle_scs_rotation_reference_to_current = left_ankle_scs_to_world_rotation_reference^(-1) * left_ankle_scs_to_world_rotation_current;
        left_euler_angles = eulerAnglesFromRotationMatrixYZX(left_ankle_scs_rotation_reference_to_current);
        gamma_left = left_euler_angles(1);
        phi_left = left_euler_angles(2);
        rho_left = left_euler_angles(3);

        % right ankle euler angles
        right_ankle_scs_transformation_current = plant.endEffectorTransformations{6};
        right_ankle_scs_to_world_rotation_current = right_ankle_scs_transformation_current(1:3, 1:3);
        right_ankle_scs_rotation_reference_to_current = right_ankle_scs_to_world_rotation_reference^(-1) * right_ankle_scs_to_world_rotation_current;
        right_euler_angles = eulerAnglesFromRotationMatrixYZX(right_ankle_scs_rotation_reference_to_current);
        gamma_right = right_euler_angles(1);
        phi_right = right_euler_angles(2);
        rho_right = right_euler_angles(3);
        
        phi_left_trajectory_left(i_time) = phi_left;
        rho_left_trajectory_left(i_time) = rho_left;
        phi_right_trajectory_left(i_time) = phi_right;
        rho_right_trajectory_left(i_time) = rho_right;

        % left foot is stance foot
        if phi_left < 0
            left_foot_constraint_number = 1;
        elseif phi_left >= 0
            left_foot_constraint_number = 2;
        else
            error('fun times, some number is neither smaller nor larger than 0.');
        end
        right_foot_constraint_number = 0;

        [A, A_dot] = ...
            createConstraintMatrix_hingeConstraints ...
              ( ...
                plant, ...
                left_foot_constraint_number, ...
                right_foot_constraint_number ...
              );

        P = eye(plant.numberOfJoints) - A' * (A * M^(-1) * A')^(-1) * A * M^(-1);
        P_tilde = [P; eye(6, plant.numberOfJoints)]; % P_tilde relates torques to accelerations, accounting for the constraints
        
        J_c = eye(2, 3) * plant.calculateCenterOfMassJacobian;
        c_two_dot_des = [1; 0];
        
        A_left_trajectory{i_time} = A;
        M_left_trajectory{i_time} = M;
        P_left_trajectory{i_time} = P;
        J_c_left_trajectory{i_time} = J_c;
        
        ankle_strategy = pinv ...
          ( ...
            [J_c * M^(-1) * P; ...                              % desired CoM acceleration
            [zeros(4, 6) eye(4, 32)] * M^(-1) * P; ...          % desired left leg accelerations
            [zeros(6, 12) eye(6, 26)] * M^(-1) * P; ...         % desired right leg accelerations
            [zeros(20, 18) eye(20, 20)] * M^(-1) * P; ...       % desired upper body accelerations
            eye(6, number_of_joints)] ...                       % torques in free dofs
          ) * ...
          [ ...
            c_two_dot_des; ...                                  % desired CoM acceleration
            zeros(4, 1); ...                                    % desired left leg accelerations
            zeros(6, 1); ...                                    % desired right leg accelerations
            zeros(20, 1); ...                                   % desired upper body accelerations
            zeros(6, 1) ...                                     % torques in free dofs
          ];
        ankle_strategy_lambda = (A*M^(-1)*A')^(-1) * (A*M^(-1)*(ankle_strategy));
        ankle_strategy_constraint_torques = A'*ankle_strategy_lambda;
        ankle_strategy_accelerations = M^(-1)*(ankle_strategy - ankle_strategy_constraint_torques);
        ankle_strategy_accelerations_normed = normVector(ankle_strategy_accelerations);
        c_two_dot_resulting_ankle_strategy = J_c * ankle_strategy_accelerations;

        hip_strategy = pinv ...
          ( ...
            [J_c * M^(-1) * P; ...                              % desired CoM acceleration
            ankle_strategy'; ...                                 % desired ankle strategy torque
            [zeros(6, 12) eye(6, 26)] * M^(-1) * P; ...         % desired right leg accelerations
            [zeros(20, 18) eye(20, 20)] * M^(-1) * P; ...       % desired upper body accelerations
            eye(6, number_of_joints)] ...                       % torques in free dofs
          ) * ...
          [ ...
            c_two_dot_des; ...                                  % desired CoM acceleration
            0; ...                                              % desired ankle strategy torque
            zeros(6, 1); ...                                    % desired right leg accelerations
            zeros(20, 1); ...                                   % desired upper body accelerations
            zeros(6, 1) ...                                     % torques in free dofs
          ];
        hip_strategy_lambda = (A*M^(-1)*A')^(-1) * (A*M^(-1)*(hip_strategy));
        hip_strategy_constraint_torques = A'*hip_strategy_lambda;
        hip_strategy_accelerations = M^(-1)*(hip_strategy - hip_strategy_constraint_torques);
        hip_strategy_accelerations_normed = normVector(hip_strategy_accelerations);
        c_two_dot_resulting_hip_strategy = J_c * hip_strategy_accelerations;
        
        left_ankle_strategy(i_time, :) = ankle_strategy;
        left_hip_strategy(i_time, :) = hip_strategy;
        
        % acceleration level
        free_joints = 1 : 6;
        left_ankle_joints = 11 : 12;
        left_leg_joints = 7 : 12;
        A_all_except_left_ankle = eye(number_of_joints);
        A_all_except_left_ankle([free_joints left_ankle_joints], :) = [];
        A_all_except_left_leg = eye(number_of_joints);
        A_all_except_left_leg([free_joints left_leg_joints], :) = [];
        
        ankle_acceleration = pinv ...
          ( ...
            [J_c; ...                                       % desired CoM acceleration
            A; ...                                          % acceleration against constraints
            A_all_except_left_ankle] ...                    % acceleration in joints other than ankle
          ) * ...
          [ ...
            c_two_dot_des; ...                              % desired CoM acceleration                   
            zeros(size(A, 1), 1); ...                       % acceleration against constraints
            zeros(size(A_all_except_left_ankle, 1), 1) ...  % acceleration in joints other than ankle
          ];

        hip_acceleration = pinv ...
          ( ...
            [J_c; ...                                       % desired CoM acceleration
            ankle_acceleration'; ...                        % desired ankle acceleration
            A; ...                                          % acceleration against constraints
            A_all_except_left_leg] ...                      % acceleration in joints other than ankle
          ) * ...
          [ ...
            c_two_dot_des; ...                              % desired CoM acceleration    
            0; ...                                          % desired ankle acceleration
            zeros(size(A, 1), 1); ...                       % acceleration against constraints
            zeros(size(A_all_except_left_leg, 1), 1) ...    % acceleration in joints other than ankle
          ];
        
        left_ankle_acceleration_strategy(i_time, :) = ankle_acceleration;
        left_hip_acceleration_strategy(i_time, :) = hip_acceleration;
    end
    disp('done');
    toc
    
    disp('Calculating ankle and hip strategy directions for right stance foot');
    tic
    % right
    right_ankle_strategy = zeros(number_of_time_steps_normalized, number_of_joints);
    right_hip_strategy = zeros(number_of_time_steps_normalized, number_of_joints);
    right_ankle_acceleration_strategy = zeros(number_of_time_steps_normalized, number_of_joints);
    right_hip_acceleration_strategy = zeros(number_of_time_steps_normalized, number_of_joints);
    for i_time = 1 : number_of_time_steps_normalized

        % update plant
        plant.jointAngles = joint_angles_mean_right_control(i_time, :)';
        plant.updateInternals;

        M = plant.inertiaMatrix;
        C = plant.coriolisMatrix;
        N = plant.gravitationalTorqueMatrix;
        
        left_ankle_scs_transformation_current = plant.endEffectorTransformations{3};
        left_ankle_scs_to_world_rotation_current = left_ankle_scs_transformation_current(1:3, 1:3);
        left_ankle_scs_rotation_reference_to_current = left_ankle_scs_to_world_rotation_reference^(-1) * left_ankle_scs_to_world_rotation_current;
        left_euler_angles = eulerAnglesFromRotationMatrixYZX(left_ankle_scs_rotation_reference_to_current);
        gamma_left = left_euler_angles(1);
        phi_left = left_euler_angles(2);
        rho_left = left_euler_angles(3);

        % right ankle euler angles
        right_ankle_scs_transformation_current = plant.endEffectorTransformations{6};
        right_ankle_scs_to_world_rotation_current = right_ankle_scs_transformation_current(1:3, 1:3);
        right_ankle_scs_rotation_reference_to_current = right_ankle_scs_to_world_rotation_reference^(-1) * right_ankle_scs_to_world_rotation_current;
        right_euler_angles = eulerAnglesFromRotationMatrixYZX(right_ankle_scs_rotation_reference_to_current);
        gamma_right = right_euler_angles(1);
        phi_right = right_euler_angles(2);
        rho_right = right_euler_angles(3);
        
        phi_left_trajectory_right(i_time) = phi_left;
        rho_left_trajectory_right(i_time) = rho_left;
        phi_right_trajectory_right(i_time) = phi_right;
        rho_right_trajectory_right(i_time) = rho_right;

        % right foot is stance foot
        if phi_right < 0
            right_foot_constraint_number = 1;
        elseif phi_right >= 0
            right_foot_constraint_number = 2;
        else
            error('fun times, some number is neither smaller nor larger than 0.');
        end
        left_foot_constraint_number = 0;

        [A, A_dot] = ...
            createConstraintMatrix_hingeConstraints ...
              ( ...
                plant, ...
                left_foot_constraint_number, ...
                right_foot_constraint_number ...
              );

        P = eye(plant.numberOfJoints) - A' * (A * M^(-1) * A')^(-1) * A * M^(-1);
        P_tilde = [P; eye(6, plant.numberOfJoints)]; % P_tilde relates torques to accelerations, accounting for the constraints
        
        J_c = eye(2, 3) * plant.calculateCenterOfMassJacobian;
        c_two_dot_des = [1; 0];
        
        ankle_strategy = pinv ...
          ( ...
            [J_c * M^(-1) * P; ...                              % desired CoM acceleration
            [zeros(6, 6) eye(6, 32)] * M^(-1) * P; ...          % desired left leg accelerations
            [zeros(4, 12) eye(4, 26)] * M^(-1) * P; ...         % desired right leg accelerations
            [zeros(20, 18) eye(20, 20)] * M^(-1) * P; ...       % desired upper body accelerations
            eye(6, number_of_joints)] ...                       % torques in free dofs
          ) * ...
          [ ...
            c_two_dot_des; ...                                  % desired CoM acceleration
            zeros(6, 1); ...                                    % desired left leg accelerations
            zeros(4, 1); ...                                    % desired right leg accelerations
            zeros(20, 1); ...                                   % desired upper body accelerations
            zeros(6, 1) ...                                     % torques in free dofs
          ];
        ankle_strategy_lambda = (A*M^(-1)*A')^(-1) * (A*M^(-1)*(ankle_strategy));
        ankle_strategy_constraint_torques = A'*ankle_strategy_lambda;
        ankle_strategy_accelerations = M^(-1)*(ankle_strategy - ankle_strategy_constraint_torques);
        ankle_strategy_accelerations_normed = normVector(ankle_strategy_accelerations);
        c_two_dot_resulting_ankle_strategy = J_c * ankle_strategy_accelerations;

        hip_strategy = pinv ...
          ( ...
            [J_c * M^(-1) * P; ...                              % desired CoM acceleration
            ankle_strategy'; ...                                 % desired ankle strategy torque
            [zeros(6, 12) eye(6, 26)] * M^(-1) * P; ...         % desired right leg accelerations
            [zeros(20, 18) eye(20, 20)] * M^(-1) * P; ...       % desired upper body accelerations
            eye(6, number_of_joints)] ...                       % torques in free dofs
          ) * ...
          [ ...
            c_two_dot_des; ...                                  % desired CoM acceleration
            0; ...                                              % desired ankle strategy torque
            zeros(6, 1); ...                                    % desired right leg accelerations
            zeros(20, 1); ...                                   % desired upper body accelerations
            zeros(6, 1) ...                                     % torques in free dofs
          ];
        hip_strategy_lambda = (A*M^(-1)*A')^(-1) * (A*M^(-1)*(hip_strategy));
        hip_strategy_constraint_torques = A'*hip_strategy_lambda;
        hip_strategy_accelerations = M^(-1)*(hip_strategy - hip_strategy_constraint_torques);
        hip_strategy_accelerations_normed = normVector(hip_strategy_accelerations);
        c_two_dot_resulting_hip_strategy = J_c * hip_strategy_accelerations;
        
        right_ankle_strategy(i_time, :) = ankle_strategy;
        right_hip_strategy(i_time, :) = hip_strategy;
        
        % acceleration level
        free_joints = 1 : 6;
        right_ankle_joints = 17 : 18;
        right_leg_joints = 13 : 18;
        A_all_except_right_ankle = eye(number_of_joints);
        A_all_except_right_ankle([free_joints right_ankle_joints], :) = [];
        A_all_except_right_leg = eye(number_of_joints);
        A_all_except_right_leg([free_joints right_leg_joints], :) = [];
        
        ankle_acceleration = pinv ...
          ( ...
            [J_c; ...                                       % desired CoM acceleration
            A; ...                                          % acceleration against constraints
            A_all_except_right_ankle] ...                    % acceleration in joints other than ankle
          ) * ...
          [ ...
            c_two_dot_des; ...                              % desired CoM acceleration                   
            zeros(size(A, 1), 1); ...                       % acceleration against constraints
            zeros(size(A_all_except_right_ankle, 1), 1) ...  % acceleration in joints other than ankle
          ];

        hip_acceleration = pinv ...
          ( ...
            [J_c; ...                                       % desired CoM acceleration
            ankle_acceleration'; ...                        % desired ankle acceleration
            A; ...                                          % acceleration against constraints
            A_all_except_right_leg] ...                      % acceleration in joints other than ankle
          ) * ...
          [ ...
            c_two_dot_des; ...                              % desired CoM acceleration    
            0; ...                                          % desired ankle acceleration
            zeros(size(A, 1), 1); ...                       % acceleration against constraints
            zeros(size(A_all_except_right_leg, 1), 1) ...    % acceleration in joints other than ankle
          ];
        
        right_ankle_acceleration_strategy(i_time, :) = ankle_acceleration;
        right_hip_acceleration_strategy(i_time, :) = hip_acceleration;
    end
    disp('done');
    toc
    
end

%% calculate strategy responses
if calculate_strategy_responses
    ankle_acceleration_strategy_left_response = zeros(number_of_time_steps_normalized, number_of_stretches);
    hip_acceleration_strategy_left_response = zeros(number_of_time_steps_normalized, number_of_stretches);
    ankle_acceleration_strategy_right_response = zeros(number_of_time_steps_normalized, number_of_stretches);
    hip_acceleration_strategy_right_response = zeros(number_of_time_steps_normalized, number_of_stretches);
    
    for i_time = 1 : number_of_time_steps_normalized
        left_ankle_acceleration_strategy_normed = normVector(left_ankle_acceleration_strategy(i_time, :));
        left_hip_acceleration_strategy_normed = normVector(left_hip_acceleration_strategy(i_time, :));
        right_ankle_acceleration_strategy_normed = normVector(right_ankle_acceleration_strategy(i_time, :));
        right_hip_acceleration_strategy_normed = normVector(right_hip_acceleration_strategy(i_time, :));
        
        for i_stretch = 1 : number_of_stretches
            % left
            delta_joint_acceleration_left = squeeze(joint_accelerations_left_response(i_time, i_stretch, :));
            ankle_acceleration_strategy_left_response(i_time, i_stretch) = dot(delta_joint_acceleration_left, left_ankle_acceleration_strategy_normed);
            hip_acceleration_strategy_left_response(i_time, i_stretch) = dot(delta_joint_acceleration_left, left_hip_acceleration_strategy_normed);
            % right
            delta_joint_acceleration_right = squeeze(joint_accelerations_right_response(i_time, i_stretch, :));
            ankle_acceleration_strategy_right_response(i_time, i_stretch) = dot(delta_joint_acceleration_right, right_ankle_acceleration_strategy_normed);
            hip_acceleration_strategy_right_response(i_time, i_stretch) = dot(delta_joint_acceleration_right, right_hip_acceleration_strategy_normed);
            
            
            
            
        end
    end
    
    
end

%% calculate stats
% if calculate_stats
%     % this is not needed anymore, I believe
%     if process_data_marker
%         % absolute data
%         lheel_x_pos_civ_stanceR_control = tinv(0.975, sum(conditions_stanceR_stimNo)-1) * std(lheel_x_pos_normalized_all(:, conditions_stanceR_stimNo), 1, 2)/sqrt(sum(conditions_stanceR_stimNo));
%         lheel_x_pos_mean_stanceR_stimPos_0ms = mean(lheel_x_pos_normalized_all(:, conditions_stanceR_stimPos_0ms), 2);
%         lheel_x_pos_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(lheel_x_pos_normalized_all(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
%         lheel_x_pos_mean_stanceR_stimNeg_0ms = mean(lheel_x_pos_normalized_all(:, conditions_stanceR_stimNeg_0ms), 2);
%         lheel_x_pos_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(lheel_x_pos_normalized_all(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
%         lheel_x_pos_mean_stanceR_stimPos_150ms = mean(lheel_x_pos_normalized_all(:, conditions_stanceR_stimPos_150ms), 2);
%         lheel_x_pos_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(lheel_x_pos_normalized_all(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
%         lheel_x_pos_mean_stanceR_stimNeg_150ms = mean(lheel_x_pos_normalized_all(:, conditions_stanceR_stimNeg_150ms), 2);
%         lheel_x_pos_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(lheel_x_pos_normalized_all(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
%         lheel_x_pos_mean_stanceR_stimPos_450ms = mean(lheel_x_pos_normalized_all(:, conditions_stanceR_stimPos_450ms), 2);
%         lheel_x_pos_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(lheel_x_pos_normalized_all(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
%         lheel_x_pos_mean_stanceR_stimNeg_450ms = mean(lheel_x_pos_normalized_all(:, conditions_stanceR_stimNeg_450ms), 2);
%         lheel_x_pos_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(lheel_x_pos_normalized_all(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));
% 
%         rheel_x_pos_civ_stanceL_control = tinv(0.975, sum(conditions_stanceL_stimNo)-1) * std(rheel_x_pos_normalized_all(:, conditions_stanceL_stimNo), 1, 2)/sqrt(sum(conditions_stanceL_stimNo));
%         rheel_x_pos_mean_stanceL_stimPos_0ms = mean(rheel_x_pos_normalized_all(:, conditions_stanceL_stimPos_0ms), 2);
%         rheel_x_pos_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(rheel_x_pos_normalized_all(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
%         rheel_x_pos_mean_stanceL_stimNeg_0ms = mean(rheel_x_pos_normalized_all(:, conditions_stanceL_stimNeg_0ms), 2);
%         rheel_x_pos_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(rheel_x_pos_normalized_all(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
%         rheel_x_pos_mean_stanceL_stimPos_150ms = mean(rheel_x_pos_normalized_all(:, conditions_stanceL_stimPos_150ms), 2);
%         rheel_x_pos_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(rheel_x_pos_normalized_all(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
%         rheel_x_pos_mean_stanceL_stimNeg_150ms = mean(rheel_x_pos_normalized_all(:, conditions_stanceL_stimNeg_150ms), 2);
%         rheel_x_pos_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(rheel_x_pos_normalized_all(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
%         rheel_x_pos_mean_stanceL_stimPos_450ms = mean(rheel_x_pos_normalized_all(:, conditions_stanceL_stimPos_450ms), 2);
%         rheel_x_pos_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(rheel_x_pos_normalized_all(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
%         rheel_x_pos_mean_stanceL_stimNeg_450ms = mean(rheel_x_pos_normalized_all(:, conditions_stanceL_stimNeg_450ms), 2);
%         rheel_x_pos_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(rheel_x_pos_normalized_all(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));
%         
%         % responses
%         lheel_x_pos_response_mean_stanceR_stimPos_0ms = mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2);
%         lheel_x_pos_response_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
%         lheel_x_pos_response_mean_stanceR_stimNeg_0ms = mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2);
%         lheel_x_pos_response_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
%         lheel_x_pos_response_mean_stanceR_stimPos_150ms = mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2);
%         lheel_x_pos_response_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
%         lheel_x_pos_response_mean_stanceR_stimNeg_150ms = mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2);
%         lheel_x_pos_response_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
%         lheel_x_pos_response_mean_stanceR_stimPos_450ms = mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2);
%         lheel_x_pos_response_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
%         lheel_x_pos_response_mean_stanceR_stimNeg_450ms = mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2);
%         lheel_x_pos_response_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));
% 
%         rheel_x_pos_response_mean_stanceL_stimPos_0ms = mean(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2);
%         rheel_x_pos_response_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
%         rheel_x_pos_response_mean_stanceL_stimNeg_0ms = mean(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2);
%         rheel_x_pos_response_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
%         rheel_x_pos_response_mean_stanceL_stimPos_150ms = mean(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2);
%         rheel_x_pos_response_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
%         rheel_x_pos_response_mean_stanceL_stimNeg_150ms = mean(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2);
%         rheel_x_pos_response_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
%         rheel_x_pos_response_mean_stanceL_stimPos_450ms = mean(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2);
%         rheel_x_pos_response_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
%         rheel_x_pos_response_mean_stanceL_stimNeg_450ms = mean(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2);
%         rheel_x_pos_response_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));
% 
%         lheel_y_pos_response_mean_stanceR_stimPos_0ms = mean(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2);
%         lheel_y_pos_response_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
%         lheel_y_pos_response_mean_stanceR_stimNeg_0ms = mean(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2);
%         lheel_y_pos_response_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
%         lheel_y_pos_response_mean_stanceR_stimPos_150ms = mean(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2);
%         lheel_y_pos_response_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
%         lheel_y_pos_response_mean_stanceR_stimNeg_150ms = mean(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2);
%         lheel_y_pos_response_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
%         lheel_y_pos_response_mean_stanceR_stimPos_450ms = mean(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2);
%         lheel_y_pos_response_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
%         lheel_y_pos_response_mean_stanceR_stimNeg_450ms = mean(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2);
%         lheel_y_pos_response_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));
% 
%         rheel_y_pos_response_mean_stanceL_stimPos_0ms = mean(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2);
%         rheel_y_pos_response_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
%         rheel_y_pos_response_mean_stanceL_stimNeg_0ms = mean(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2);
%         rheel_y_pos_response_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
%         rheel_y_pos_response_mean_stanceL_stimPos_150ms = mean(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2);
%         rheel_y_pos_response_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
%         rheel_y_pos_response_mean_stanceL_stimNeg_150ms = mean(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2);
%         rheel_y_pos_response_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
%         rheel_y_pos_response_mean_stanceL_stimPos_450ms = mean(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2);
%         rheel_y_pos_response_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
%         rheel_y_pos_response_mean_stanceL_stimNeg_450ms = mean(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2);
%         rheel_y_pos_response_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));
% 
%     end
% 
%     if process_data_forceplate
%         % absolute data
%         fxl_civ_stanceL_control = tinv(0.975, sum(conditions_stanceL_stimNo)-1) * std(fxl_normalized_all(:, conditions_stanceL_stimNo), 1, 2)/sqrt(sum(conditions_stanceL_stimNo));
%         fxl_mean_stanceL_stimPos_0ms = mean(fxl_normalized_all(:, conditions_stanceL_stimPos_0ms), 2);
%         fxl_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(fxl_normalized_all(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
%         fxl_mean_stanceL_stimNeg_0ms = mean(fxl_normalized_all(:, conditions_stanceL_stimNeg_0ms), 2);
%         fxl_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(fxl_normalized_all(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
%         fxl_mean_stanceL_stimPos_150ms = mean(fxl_normalized_all(:, conditions_stanceL_stimPos_150ms), 2);
%         fxl_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(fxl_normalized_all(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
%         fxl_mean_stanceL_stimNeg_150ms = mean(fxl_normalized_all(:, conditions_stanceL_stimNeg_150ms), 2);
%         fxl_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(fxl_normalized_all(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
%         fxl_mean_stanceL_stimPos_450ms = mean(fxl_normalized_all(:, conditions_stanceL_stimPos_450ms), 2);
%         fxl_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(fxl_normalized_all(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
%         fxl_mean_stanceL_stimNeg_450ms = mean(fxl_normalized_all(:, conditions_stanceL_stimNeg_450ms), 2);
%         fxl_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(fxl_normalized_all(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));
%         
%         fzl_civ_stanceL_control = tinv(0.975, sum(conditions_stanceL_stimNo)-1) * std(fzl_normalized_all(:, conditions_stanceL_stimNo), 1, 2)/sqrt(sum(conditions_stanceL_stimNo));
%         fzl_mean_stanceL_stimPos_0ms = mean(fzl_normalized_all(:, conditions_stanceL_stimPos_0ms), 2);
%         fzl_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(fzl_normalized_all(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
%         fzl_mean_stanceL_stimNeg_0ms = mean(fzl_normalized_all(:, conditions_stanceL_stimNeg_0ms), 2);
%         fzl_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(fzl_normalized_all(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
%         fzl_mean_stanceL_stimPos_150ms = mean(fzl_normalized_all(:, conditions_stanceL_stimPos_150ms), 2);
%         fzl_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(fzl_normalized_all(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
%         fzl_mean_stanceL_stimNeg_150ms = mean(fzl_normalized_all(:, conditions_stanceL_stimNeg_150ms), 2);
%         fzl_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(fzl_normalized_all(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
%         fzl_mean_stanceL_stimPos_450ms = mean(fzl_normalized_all(:, conditions_stanceL_stimPos_450ms), 2);
%         fzl_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(fzl_normalized_all(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
%         fzl_mean_stanceL_stimNeg_450ms = mean(fzl_normalized_all(:, conditions_stanceL_stimNeg_450ms), 2);
%         fzl_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(fzl_normalized_all(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));
% 
%         myl_civ_stanceL_control = tinv(0.975, sum(conditions_stanceL_stimNo)-1) * std(myl_normalized_all(:, conditions_stanceL_stimNo), 1, 2)/sqrt(sum(conditions_stanceL_stimNo));
%         myl_mean_stanceL_stimPos_0ms = mean(myl_normalized_all(:, conditions_stanceL_stimPos_0ms), 2);
%         myl_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(myl_normalized_all(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
%         myl_mean_stanceL_stimNeg_0ms = mean(myl_normalized_all(:, conditions_stanceL_stimNeg_0ms), 2);
%         myl_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(myl_normalized_all(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
%         myl_mean_stanceL_stimPos_150ms = mean(myl_normalized_all(:, conditions_stanceL_stimPos_150ms), 2);
%         myl_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(myl_normalized_all(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
%         myl_mean_stanceL_stimNeg_150ms = mean(myl_normalized_all(:, conditions_stanceL_stimNeg_150ms), 2);
%         myl_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(myl_normalized_all(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
%         myl_mean_stanceL_stimPos_450ms = mean(myl_normalized_all(:, conditions_stanceL_stimPos_450ms), 2);
%         myl_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(myl_normalized_all(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
%         myl_mean_stanceL_stimNeg_450ms = mean(myl_normalized_all(:, conditions_stanceL_stimNeg_450ms), 2);
%         myl_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(myl_normalized_all(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));
%         
%         lcop_x_civ_stanceL_control = tinv(0.975, sum(conditions_stanceL_stimNo)-1) * std(lcop_x_normalized_all(:, conditions_stanceL_stimNo), 1, 2)/sqrt(sum(conditions_stanceL_stimNo));
%         lcop_x_mean_stanceL_stimPos_0ms = mean(lcop_x_normalized_all(:, conditions_stanceL_stimPos_0ms), 2);
%         lcop_x_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(lcop_x_normalized_all(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
%         lcop_x_mean_stanceL_stimNeg_0ms = mean(lcop_x_normalized_all(:, conditions_stanceL_stimNeg_0ms), 2);
%         lcop_x_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(lcop_x_normalized_all(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
%         lcop_x_mean_stanceL_stimPos_150ms = mean(lcop_x_normalized_all(:, conditions_stanceL_stimPos_150ms), 2);
%         lcop_x_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(lcop_x_normalized_all(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
%         lcop_x_mean_stanceL_stimNeg_150ms = mean(lcop_x_normalized_all(:, conditions_stanceL_stimNeg_150ms), 2);
%         lcop_x_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(lcop_x_normalized_all(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
%         lcop_x_mean_stanceL_stimPos_450ms = mean(lcop_x_normalized_all(:, conditions_stanceL_stimPos_450ms), 2);
%         lcop_x_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(lcop_x_normalized_all(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
%         lcop_x_mean_stanceL_stimNeg_450ms = mean(lcop_x_normalized_all(:, conditions_stanceL_stimNeg_450ms), 2);
%         lcop_x_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(lcop_x_normalized_all(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));
% 
%         fxr_civ_stanceR_control = tinv(0.975, sum(conditions_stanceR_stimNo)-1) * std(fxr_normalized_all(:, conditions_stanceR_stimNo), 1, 2)/sqrt(sum(conditions_stanceR_stimNo));
%         fxr_mean_stanceR_stimPos_0ms = mean(fxr_normalized_all(:, conditions_stanceR_stimPos_0ms), 2);
%         fxr_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(fxr_normalized_all(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
%         fxr_mean_stanceR_stimNeg_0ms = mean(fxr_normalized_all(:, conditions_stanceR_stimNeg_0ms), 2);
%         fxr_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(fxr_normalized_all(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
%         fxr_mean_stanceR_stimPos_150ms = mean(fxr_normalized_all(:, conditions_stanceR_stimPos_150ms), 2);
%         fxr_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(fxr_normalized_all(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
%         fxr_mean_stanceR_stimNeg_150ms = mean(fxr_normalized_all(:, conditions_stanceR_stimNeg_150ms), 2);
%         fxr_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(fxr_normalized_all(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
%         fxr_mean_stanceR_stimPos_450ms = mean(fxr_normalized_all(:, conditions_stanceR_stimPos_450ms), 2);
%         fxr_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(fxr_normalized_all(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
%         fxr_mean_stanceR_stimNeg_450ms = mean(fxr_normalized_all(:, conditions_stanceR_stimNeg_450ms), 2);
%         fxr_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(fxr_normalized_all(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));
% 
%         fzr_civ_stanceR_control = tinv(0.975, sum(conditions_stanceR_stimNo)-1) * std(fzr_normalized_all(:, conditions_stanceR_stimNo), 1, 2)/sqrt(sum(conditions_stanceR_stimNo));
%         fzr_mean_stanceR_stimPos_0ms = mean(fzr_normalized_all(:, conditions_stanceR_stimPos_0ms), 2);
%         fzr_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(fzr_normalized_all(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
%         fzr_mean_stanceR_stimNeg_0ms = mean(fzr_normalized_all(:, conditions_stanceR_stimNeg_0ms), 2);
%         fzr_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(fzr_normalized_all(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
%         fzr_mean_stanceR_stimPos_150ms = mean(fzr_normalized_all(:, conditions_stanceR_stimPos_150ms), 2);
%         fzr_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(fzr_normalized_all(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
%         fzr_mean_stanceR_stimNeg_150ms = mean(fzr_normalized_all(:, conditions_stanceR_stimNeg_150ms), 2);
%         fzr_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(fzr_normalized_all(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
%         fzr_mean_stanceR_stimPos_450ms = mean(fzr_normalized_all(:, conditions_stanceR_stimPos_450ms), 2);
%         fzr_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(fzr_normalized_all(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
%         fzr_mean_stanceR_stimNeg_450ms = mean(fzr_normalized_all(:, conditions_stanceR_stimNeg_450ms), 2);
%         fzr_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(fzr_normalized_all(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));
% 
%         myr_civ_stanceR_control = tinv(0.975, sum(conditions_stanceR_stimNo)-1) * std(myr_normalized_all(:, conditions_stanceR_stimNo), 1, 2)/sqrt(sum(conditions_stanceR_stimNo));
%         myr_mean_stanceR_stimPos_0ms = mean(myr_normalized_all(:, conditions_stanceR_stimPos_0ms), 2);
%         myr_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(myr_normalized_all(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
%         myr_mean_stanceR_stimNeg_0ms = mean(myr_normalized_all(:, conditions_stanceR_stimNeg_0ms), 2);
%         myr_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(myr_normalized_all(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
%         myr_mean_stanceR_stimPos_150ms = mean(myr_normalized_all(:, conditions_stanceR_stimPos_150ms), 2);
%         myr_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(myr_normalized_all(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
%         myr_mean_stanceR_stimNeg_150ms = mean(myr_normalized_all(:, conditions_stanceR_stimNeg_150ms), 2);
%         myr_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(myr_normalized_all(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
%         myr_mean_stanceR_stimPos_450ms = mean(myr_normalized_all(:, conditions_stanceR_stimPos_450ms), 2);
%         myr_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(myr_normalized_all(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
%         myr_mean_stanceR_stimNeg_450ms = mean(myr_normalized_all(:, conditions_stanceR_stimNeg_450ms), 2);
%         myr_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(myr_normalized_all(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));
% 
%         rcop_x_civ_stanceR_control = tinv(0.975, sum(conditions_stanceR_stimNo)-1) * std(rcop_x_normalized_all(:, conditions_stanceR_stimNo), 1, 2)/sqrt(sum(conditions_stanceR_stimNo));
%         rcop_x_mean_stanceR_stimPos_0ms = mean(rcop_x_normalized_all(:, conditions_stanceR_stimPos_0ms), 2);
%         rcop_x_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(rcop_x_normalized_all(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
%         rcop_x_mean_stanceR_stimNeg_0ms = mean(rcop_x_normalized_all(:, conditions_stanceR_stimNeg_0ms), 2);
%         rcop_x_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(rcop_x_normalized_all(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
%         rcop_x_mean_stanceR_stimPos_150ms = mean(rcop_x_normalized_all(:, conditions_stanceR_stimPos_150ms), 2);
%         rcop_x_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(rcop_x_normalized_all(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
%         rcop_x_mean_stanceR_stimNeg_150ms = mean(rcop_x_normalized_all(:, conditions_stanceR_stimNeg_150ms), 2);
%         rcop_x_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(rcop_x_normalized_all(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
%         rcop_x_mean_stanceR_stimPos_450ms = mean(rcop_x_normalized_all(:, conditions_stanceR_stimPos_450ms), 2);
%         rcop_x_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(rcop_x_normalized_all(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
%         rcop_x_mean_stanceR_stimNeg_450ms = mean(rcop_x_normalized_all(:, conditions_stanceR_stimNeg_450ms), 2);
%         rcop_x_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(rcop_x_normalized_all(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));
% 
%         % responses
%         fxl_response_mean_stanceL_stimPos_0ms = mean(fxl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2);
%         fxl_response_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(fxl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
%         fxl_response_mean_stanceL_stimNeg_0ms = mean(fxl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2);
%         fxl_response_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(fxl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
%         fxl_response_mean_stanceL_stimPos_150ms = mean(fxl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2);
%         fxl_response_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(fxl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
%         fxl_response_mean_stanceL_stimNeg_150ms = mean(fxl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2);
%         fxl_response_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(fxl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
%         fxl_response_mean_stanceL_stimPos_450ms = mean(fxl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2);
%         fxl_response_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(fxl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
%         fxl_response_mean_stanceL_stimNeg_450ms = mean(fxl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2);
%         fxl_response_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(fxl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));
% 
%         fzl_response_mean_stanceL_stimPos_0ms = mean(fzl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2);
%         fzl_response_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(fzl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
%         fzl_response_mean_stanceL_stimNeg_0ms = mean(fzl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2);
%         fzl_response_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(fzl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
%         fzl_response_mean_stanceL_stimPos_150ms = mean(fzl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2);
%         fzl_response_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(fzl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
%         fzl_response_mean_stanceL_stimNeg_150ms = mean(fzl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2);
%         fzl_response_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(fzl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
%         fzl_response_mean_stanceL_stimPos_450ms = mean(fzl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2);
%         fzl_response_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(fzl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
%         fzl_response_mean_stanceL_stimNeg_450ms = mean(fzl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2);
%         fzl_response_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(fzl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));
% 
%         myl_response_mean_stanceL_stimPos_0ms = mean(myl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2);
%         myl_response_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(myl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
%         myl_response_mean_stanceL_stimNeg_0ms = mean(myl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2);
%         myl_response_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(myl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
%         myl_response_mean_stanceL_stimPos_150ms = mean(myl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2);
%         myl_response_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(myl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
%         myl_response_mean_stanceL_stimNeg_150ms = mean(myl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2);
%         myl_response_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(myl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
%         myl_response_mean_stanceL_stimPos_450ms = mean(myl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2);
%         myl_response_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(myl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
%         myl_response_mean_stanceL_stimNeg_450ms = mean(myl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2);
%         myl_response_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(myl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));
% 
%         lcop_x_response_mean_stanceL_stimPos_0ms = mean(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2);
%         lcop_x_response_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
%         lcop_x_response_mean_stanceL_stimNeg_0ms = mean(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2);
%         lcop_x_response_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
%         lcop_x_response_mean_stanceL_stimPos_150ms = mean(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2);
%         lcop_x_response_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
%         lcop_x_response_mean_stanceL_stimNeg_150ms = mean(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2);
%         lcop_x_response_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
%         lcop_x_response_mean_stanceL_stimPos_450ms = mean(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2);
%         lcop_x_response_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
%         lcop_x_response_mean_stanceL_stimNeg_450ms = mean(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2);
%         lcop_x_response_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));
% 
%         fxr_response_mean_stanceR_stimPos_0ms = mean(fxr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2);
%         fxr_response_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(fxr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
%         fxr_response_mean_stanceR_stimNeg_0ms = mean(fxr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2);
%         fxr_response_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(fxr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
%         fxr_response_mean_stanceR_stimPos_150ms = mean(fxr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2);
%         fxr_response_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(fxr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
%         fxr_response_mean_stanceR_stimNeg_150ms = mean(fxr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2);
%         fxr_response_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(fxr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
%         fxr_response_mean_stanceR_stimPos_450ms = mean(fxr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2);
%         fxr_response_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(fxr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
%         fxr_response_mean_stanceR_stimNeg_450ms = mean(fxr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2);
%         fxr_response_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(fxr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));
% 
%         fzr_response_mean_stanceR_stimPos_0ms = mean(fzr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2);
%         fzr_response_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(fzr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
%         fzr_response_mean_stanceR_stimNeg_0ms = mean(fzr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2);
%         fzr_response_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(fzr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
%         fzr_response_mean_stanceR_stimPos_150ms = mean(fzr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2);
%         fzr_response_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(fzr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
%         fzr_response_mean_stanceR_stimNeg_150ms = mean(fzr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2);
%         fzr_response_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(fzr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
%         fzr_response_mean_stanceR_stimPos_450ms = mean(fzr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2);
%         fzr_response_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(fzr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
%         fzr_response_mean_stanceR_stimNeg_450ms = mean(fzr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2);
%         fzr_response_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(fzr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));
% 
%         myr_response_mean_stanceR_stimPos_0ms = mean(myr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2);
%         myr_response_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(myr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
%         myr_response_mean_stanceR_stimNeg_0ms = mean(myr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2);
%         myr_response_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(myr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
%         myr_response_mean_stanceR_stimPos_150ms = mean(myr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2);
%         myr_response_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(myr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
%         myr_response_mean_stanceR_stimNeg_150ms = mean(myr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2);
%         myr_response_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(myr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
%         myr_response_mean_stanceR_stimPos_450ms = mean(myr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2);
%         myr_response_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(myr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
%         myr_response_mean_stanceR_stimNeg_450ms = mean(myr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2);
%         myr_response_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(myr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));
%     end
%     
%     if process_data_angles
%         % absolute data
%         joint_angles_civ_stanceL_control = tinv(0.975, sum(conditions_stanceL_stimNo)-1) * squeeze(std(joint_angles_normalized_all(:, conditions_stanceL_stimNo, :), 1, 2))/sqrt(sum(conditions_stanceL_stimNo));
%         joint_angles_mean_stanceL_stimPos_0ms = squeeze(mean(joint_angles_normalized_all(:, conditions_stanceL_stimPos_0ms, :), 2));
%         joint_angles_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * squeeze(std(joint_angles_normalized_all(:, conditions_stanceL_stimPos_0ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimPos_0ms));
%         joint_angles_mean_stanceL_stimNeg_0ms = squeeze(mean(joint_angles_normalized_all(:, conditions_stanceL_stimNeg_0ms, :), 2));
%         joint_angles_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * squeeze(std(joint_angles_normalized_all(:, conditions_stanceL_stimNeg_0ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimNeg_0ms));
%         joint_angles_mean_stanceL_stimPos_150ms = squeeze(mean(joint_angles_normalized_all(:, conditions_stanceL_stimPos_150ms, :), 2));
%         joint_angles_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * squeeze(std(joint_angles_normalized_all(:, conditions_stanceL_stimPos_150ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimPos_150ms));
%         joint_angles_mean_stanceL_stimNeg_150ms = squeeze(mean(joint_angles_normalized_all(:, conditions_stanceL_stimNeg_150ms, :), 2));
%         joint_angles_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * squeeze(std(joint_angles_normalized_all(:, conditions_stanceL_stimNeg_150ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimNeg_150ms));
%         joint_angles_mean_stanceL_stimPos_450ms = squeeze(mean(joint_angles_normalized_all(:, conditions_stanceL_stimPos_450ms, :), 2));
%         joint_angles_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * squeeze(std(joint_angles_normalized_all(:, conditions_stanceL_stimPos_450ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimPos_450ms));
%         joint_angles_mean_stanceL_stimNeg_450ms = squeeze(mean(joint_angles_normalized_all(:, conditions_stanceL_stimNeg_450ms, :), 2));
%         joint_angles_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * squeeze(std(joint_angles_normalized_all(:, conditions_stanceL_stimNeg_450ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimNeg_450ms));
% 
%         joint_angles_civ_stanceR_control = tinv(0.975, sum(conditions_stanceR_stimNo)-1) * squeeze(std(joint_angles_normalized_all(:, conditions_stanceR_stimNo, :), 1, 2))/sqrt(sum(conditions_stanceR_stimNo));
%         joint_angles_mean_stanceR_stimPos_0ms = squeeze(mean(joint_angles_normalized_all(:, conditions_stanceR_stimPos_0ms, :), 2));
%         joint_angles_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * squeeze(std(joint_angles_normalized_all(:, conditions_stanceR_stimPos_0ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimPos_0ms));
%         joint_angles_mean_stanceR_stimNeg_0ms = squeeze(mean(joint_angles_normalized_all(:, conditions_stanceR_stimNeg_0ms, :), 2));
%         joint_angles_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * squeeze(std(joint_angles_normalized_all(:, conditions_stanceR_stimNeg_0ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimNeg_0ms));
%         joint_angles_mean_stanceR_stimPos_150ms = squeeze(mean(joint_angles_normalized_all(:, conditions_stanceR_stimPos_150ms, :), 2));
%         joint_angles_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * squeeze(std(joint_angles_normalized_all(:, conditions_stanceR_stimPos_150ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimPos_150ms));
%         joint_angles_mean_stanceR_stimNeg_150ms = squeeze(mean(joint_angles_normalized_all(:, conditions_stanceR_stimNeg_150ms, :), 2));
%         joint_angles_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * squeeze(std(joint_angles_normalized_all(:, conditions_stanceR_stimNeg_150ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimNeg_150ms));
%         joint_angles_mean_stanceR_stimPos_450ms = squeeze(mean(joint_angles_normalized_all(:, conditions_stanceR_stimPos_450ms, :), 2));
%         joint_angles_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * squeeze(std(joint_angles_normalized_all(:, conditions_stanceR_stimPos_450ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimPos_450ms));
%         joint_angles_mean_stanceR_stimNeg_450ms = squeeze(mean(joint_angles_normalized_all(:, conditions_stanceR_stimNeg_450ms, :), 2));
%         joint_angles_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * squeeze(std(joint_angles_normalized_all(:, conditions_stanceR_stimNeg_450ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimNeg_450ms));
% 
%         % responses
%         joint_angles_response_mean_stanceL_stimPos_0ms = squeeze(mean(joint_angles_stanceL_response(:, conditions_stanceL_stimPos_0ms, :), 2));
%         joint_angles_response_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * squeeze(std(joint_angles_stanceL_response(:, conditions_stanceL_stimPos_0ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimPos_0ms));
%         joint_angles_response_mean_stanceL_stimNeg_0ms = squeeze(mean(joint_angles_stanceL_response(:, conditions_stanceL_stimNeg_0ms, :), 2));
%         joint_angles_response_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * squeeze(std(joint_angles_stanceL_response(:, conditions_stanceL_stimNeg_0ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimNeg_0ms));
%         joint_angles_response_mean_stanceL_stimPos_150ms = squeeze(mean(joint_angles_stanceL_response(:, conditions_stanceL_stimPos_150ms, :), 2));
%         joint_angles_response_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * squeeze(std(joint_angles_stanceL_response(:, conditions_stanceL_stimPos_150ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimPos_150ms));
%         joint_angles_response_mean_stanceL_stimNeg_150ms = squeeze(mean(joint_angles_stanceL_response(:, conditions_stanceL_stimNeg_150ms, :), 2));
%         joint_angles_response_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * squeeze(std(joint_angles_stanceL_response(:, conditions_stanceL_stimNeg_150ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimNeg_150ms));
%         joint_angles_response_mean_stanceL_stimPos_450ms = squeeze(mean(joint_angles_stanceL_response(:, conditions_stanceL_stimPos_450ms, :), 2));
%         joint_angles_response_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * squeeze(std(joint_angles_stanceL_response(:, conditions_stanceL_stimPos_450ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimPos_450ms));
%         joint_angles_response_mean_stanceL_stimNeg_450ms = squeeze(mean(joint_angles_stanceL_response(:, conditions_stanceL_stimNeg_450ms, :), 2));
%         joint_angles_response_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * squeeze(std(joint_angles_stanceL_response(:, conditions_stanceL_stimNeg_450ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimNeg_450ms));
% 
%         joint_angles_response_mean_stanceR_stimPos_0ms = squeeze(mean(joint_angles_stanceR_response(:, conditions_stanceR_stimPos_0ms, :), 2));
%         joint_angles_response_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * squeeze(std(joint_angles_stanceR_response(:, conditions_stanceR_stimPos_0ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimPos_0ms));
%         joint_angles_response_mean_stanceR_stimNeg_0ms = squeeze(mean(joint_angles_stanceR_response(:, conditions_stanceR_stimNeg_0ms, :), 2));
%         joint_angles_response_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * squeeze(std(joint_angles_stanceR_response(:, conditions_stanceR_stimNeg_0ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimNeg_0ms));
%         joint_angles_response_mean_stanceR_stimPos_150ms = squeeze(mean(joint_angles_stanceR_response(:, conditions_stanceR_stimPos_150ms, :), 2));
%         joint_angles_response_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * squeeze(std(joint_angles_stanceR_response(:, conditions_stanceR_stimPos_150ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimPos_150ms));
%         joint_angles_response_mean_stanceR_stimNeg_150ms = squeeze(mean(joint_angles_stanceR_response(:, conditions_stanceR_stimNeg_150ms, :), 2));
%         joint_angles_response_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * squeeze(std(joint_angles_stanceR_response(:, conditions_stanceR_stimNeg_150ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimNeg_150ms));
%         joint_angles_response_mean_stanceR_stimPos_450ms = squeeze(mean(joint_angles_stanceR_response(:, conditions_stanceR_stimPos_450ms, :), 2));
%         joint_angles_response_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * squeeze(std(joint_angles_stanceR_response(:, conditions_stanceR_stimPos_450ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimPos_450ms));
%         joint_angles_response_mean_stanceR_stimNeg_450ms = squeeze(mean(joint_angles_stanceR_response(:, conditions_stanceR_stimNeg_450ms, :), 2));
%         joint_angles_response_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * squeeze(std(joint_angles_stanceR_response(:, conditions_stanceR_stimNeg_450ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimNeg_450ms));
% 
%      
%     end
%     
%     % stim start times
%     stim_start_times_mean_stanceL_0ms = mean(stim_start_time_relative_to_stretch_all(conditions_stanceL_0ms));
%     stim_start_times_mean_stanceR_0ms = mean(stim_start_time_relative_to_stretch_all(conditions_stanceR_0ms));
%     stim_start_times_mean_stanceL_150ms = mean(stim_start_time_relative_to_stretch_all(conditions_stanceL_150ms));
%     stim_start_times_mean_stanceR_150ms = mean(stim_start_time_relative_to_stretch_all(conditions_stanceR_150ms));
%     stim_start_times_mean_stanceL_450ms = mean(stim_start_time_relative_to_stretch_all(conditions_stanceL_450ms));
%     stim_start_times_mean_stanceR_450ms = mean(stim_start_time_relative_to_stretch_all(conditions_stanceR_450ms));
%     
% end


%% save data old
% if save_data_old
%     save ...
%       ( ...
%         makeFileName(date, subject_id, 'resultsConditions'), ...
%         'time_normalized', ...
%         'origin_trial_list_all', ...
%         'origin_start_time_list_all', ...
%         'origin_end_time_list_all', ...
%         'condition_stance_foot_list_all', ...
%         'condition_perturbation_list_all', ...
%         'condition_delay_list_all', ...
%         'conditions_stanceL_stimNo', ...
%         'conditions_stanceR_stimNo', ...
%         'conditions_stanceL_stimPos_0ms', ...
%         'conditions_stanceL_stimNeg_0ms', ...
%         'conditions_stanceR_stimPos_0ms', ...
%         'conditions_stanceR_stimNeg_0ms', ...
%         'conditions_stanceL_stimPos_150ms', ...
%         'conditions_stanceL_stimNeg_150ms', ...
%         'conditions_stanceR_stimPos_150ms', ...
%         'conditions_stanceR_stimNeg_150ms', ...
%         'conditions_stanceL_stimPos_450ms', ...
%         'conditions_stanceL_stimNeg_450ms', ...
%         'conditions_stanceR_stimPos_450ms', ...
%         'conditions_stanceR_stimNeg_450ms', ...
%         'conditions_stanceL_0ms', ...
%         'conditions_stanceR_0ms', ...
%         'conditions_stanceL_150ms', ...
%         'conditions_stanceR_150ms', ...
%         'conditions_stanceL_450ms', ...
%         'conditions_stanceR_450ms' ...
%       );
%     
%     if process_data_marker
%         save ...
%           ( ...
%             makeFileName(date, subject_id, 'resultsMarker'), ...
%             'step_response_all', ...
%             'stim_response_all', ...
%             'lpsi_x_pos_normalized_all', ...
%             'rpsi_x_pos_normalized_all', ...
%             'lheel_x_pos_normalized_all', ...
%             'rheel_x_pos_normalized_all', ...
%             'lheel_y_pos_normalized_all', ...
%             'rheel_y_pos_normalized_all', ...
%             'lheel_x_pos_stanceR_response', ...
%             'rheel_x_pos_stanceL_response', ...
%             'lheel_y_pos_stanceR_response', ...
%             'rheel_y_pos_stanceL_response' ...
%           );
%     end
%             
%     if process_data_forceplate
%         save ...
%           ( ...
%             makeFileName(date, subject_id, 'resultsForceplate'), ...
%             'fxl_normalized_all', ...
%             'fzl_normalized_all', ...
%             'myl_normalized_all', ...
%             'lcop_x_normalized_all', ...
%             'fxr_normalized_all', ...
%             'fzr_normalized_all', ...
%             'myr_normalized_all', ...
%             'rcop_x_normalized_all', ...
%             'lcop_x_stanceL_response', ...
%             'rcop_x_stanceR_response' ...
%           );
%     end
%     
%     if process_data_emg
%         save ...
%           ( ...
%             makeFileName(date, subject_id, 'resultsEmg'), ...
%             'lglutmed_normalized_all', ...
%             'ltibiant_normalized_all', ...
%             'lperolng_normalized_all', ...
%             'rglutmed_normalized_all', ...
%             'rtibiant_normalized_all', ...
%             'rperolng_normalized_all', ...
%             'lglutmed_stanceL_response', ...
%             'ltibiant_stanceL_response', ...
%             'lperolng_stanceL_response', ...
%             'rglutmed_stanceL_response', ...
%             'rtibiant_stanceL_response', ...
%             'rperolng_stanceL_response', ...
%             'lglutmed_stanceR_response', ...
%             'ltibiant_stanceR_response', ...
%             'lperolng_stanceR_response', ...
%             'rglutmed_stanceR_response', ...
%             'rtibiant_stanceR_response', ...
%             'rperolng_stanceR_response' ...
%           )
%     end
%     
% end
    

return    

%% angle-related stuff
if do_joint_angle_absolute_plots_left
    joints_to_plot = [8 12 14 18];
%     joints_to_plot = [20 23];
    for i_joint = joints_to_plot
        figure; axes; hold on; title([plant.jointLabels{i_joint} ', LEFT 0ms']); set(gca, 'Fontsize', 12)
        plot(time_normalized, squeeze(joint_angles_normalized_all(:, conditions_stanceL_stimNo, i_joint)), 'color', color_left_control);
        plot(time_normalized, squeeze(joint_angles_normalized_all(:, conditions_stanceL_stimPos_0ms, i_joint)), 'color', color_left_positive);
        plot(time_normalized, squeeze(joint_angles_normalized_all(:, conditions_stanceL_stimNeg_0ms, i_joint)), 'color', color_left_negative);
        plot(time_normalized, mean(squeeze(joint_angles_normalized_all(:, conditions_stanceL_stimNo, i_joint)), 2), 'color', color_left_control, 'linewidth', 5);
        plot(time_normalized, mean(squeeze(joint_angles_normalized_all(:, conditions_stanceL_stimPos_0ms, i_joint)), 2), 'color', color_left_positive, 'linewidth', 5);
        plot(time_normalized, mean(squeeze(joint_angles_normalized_all(:, conditions_stanceL_stimNeg_0ms, i_joint)), 2), 'color', color_left_negative, 'linewidth', 5);
    end    
end

if do_joint_angle_absolute_plots_right
    joints_to_plot = [8 12 14 18];
%     joints_to_plot = [20 23];
    for i_joint = joints_to_plot
        figure; axes; hold on; title([plant.jointLabels{i_joint} ', RIGHT 0ms']); set(gca, 'Fontsize', 12)
        plot(time_normalized, squeeze(joint_angles_normalized_all(:, conditions_stanceR_stimNo, i_joint)), 'color', color_right_control);
        plot(time_normalized, squeeze(joint_angles_normalized_all(:, conditions_stanceR_stimPos_0ms, i_joint)), 'color', color_right_positive);
        plot(time_normalized, squeeze(joint_angles_normalized_all(:, conditions_stanceR_stimNeg_0ms, i_joint)), 'color', color_right_negative);
        plot(time_normalized, mean(squeeze(joint_angles_normalized_all(:, conditions_stanceR_stimNo, i_joint)), 2), 'color', color_right_control, 'linewidth', 5);
        plot(time_normalized, mean(squeeze(joint_angles_normalized_all(:, conditions_stanceR_stimPos_0ms, i_joint)), 2), 'color', color_right_positive, 'linewidth', 5);
        plot(time_normalized, mean(squeeze(joint_angles_normalized_all(:, conditions_stanceR_stimNeg_0ms, i_joint)), 2), 'color', color_right_negative, 'linewidth', 5);
    end    
end


if do_acceleration_strategy_response_plots
%     figure; axes; hold on; title('left ankle strategy, 0ms'); set(gca, 'Fontsize', 12)
%     plot(time_normalized, ankle_acceleration_strategy_left_response(:, conditions_left_control), 'color', color_left_control);
%     plot(time_normalized, ankle_acceleration_strategy_left_response(:, conditions_left_positive_0ms), 'color', color_left_positive);
%     plot(time_normalized, ankle_acceleration_strategy_left_response(:, conditions_left_negative_0ms), 'color', color_left_negative);
%     
%     figure; axes; hold on; title('left hip strategy, 0ms'); set(gca, 'Fontsize', 12)
%     plot(time_normalized, hip_acceleration_strategy_left_response(:, conditions_left_control), 'color', color_left_control);
%     plot(time_normalized, hip_acceleration_strategy_left_response(:, conditions_left_positive_0ms), 'color', color_left_positive);
%     plot(time_normalized, hip_acceleration_strategy_left_response(:, conditions_left_negative_0ms), 'color', color_left_negative);
    
    figure; axes; hold on; title('right ankle strategy, 0ms'); set(gca, 'Fontsize', 12)
    plot(time_normalized, ankle_acceleration_strategy_right_response(:, conditions_stanceR_stimNo), 'color', color_right_control);
    plot(time_normalized, ankle_acceleration_strategy_right_response(:, conditions_stanceR_stimPos_0ms), 'color', color_right_positive);
    plot(time_normalized, ankle_acceleration_strategy_right_response(:, conditions_stanceR_stimNeg_0ms), 'color', color_right_negative);
    plot(time_normalized, mean(ankle_acceleration_strategy_right_response(:, conditions_stanceR_stimNo), 2), 'color', color_right_control, 'linewidth', 5);
    plot(time_normalized, mean(ankle_acceleration_strategy_right_response(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_right_positive, 'linewidth', 5);
    plot(time_normalized, mean(ankle_acceleration_strategy_right_response(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_right_negative, 'linewidth', 5);
    
    figure; axes; hold on; title('right hip strategy, 0ms'); set(gca, 'Fontsize', 12)
    plot(time_normalized, hip_acceleration_strategy_right_response(:, conditions_stanceR_stimNo), 'color', color_right_control);
    plot(time_normalized, hip_acceleration_strategy_right_response(:, conditions_stanceR_stimPos_0ms), 'color', color_right_positive);
    plot(time_normalized, hip_acceleration_strategy_right_response(:, conditions_stanceR_stimNeg_0ms), 'color', color_right_negative);
    plot(time_normalized, mean(hip_acceleration_strategy_right_response(:, conditions_stanceR_stimNo), 2), 'color', color_right_control, 'linewidth', 5);
    plot(time_normalized, mean(hip_acceleration_strategy_right_response(:, conditions_stanceR_stimPos_0ms), 2), 'color', color_right_positive, 'linewidth', 5);
    plot(time_normalized, mean(hip_acceleration_strategy_right_response(:, conditions_stanceR_stimNeg_0ms), 2), 'color', color_right_negative, 'linewidth', 5);
end














