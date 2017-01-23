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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% compare the kinematic tree against the kinematic chain

function analyzeGaitParameters(varargin)
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    
    save_data                           = 1;
    use_stance_foot_as_reference        = 1;

    
    
    % setup analysis pattern - this should be done by editing a file
    process_data_balance                = 1;
    process_data_armswing               = 1;
    process_data_forceplate             = 1;
    process_data_emg                    = 0;
    process_data_angles                 = 0;
    process_data_torques                = 0;

    wait_times = [0 0.150 0.450];
    wait_time_labels = {'0ms', '150ms', '450ms'};

    % ArmSense
    condition_labels = {'stance foot', 'experimental'};
    conditions_to_analyze = ...
      {
        'STANCE_LEFT', 'baselineOG'; ...
        'STANCE_RIGHT', 'baselineOG'; ...
        'STANCE_LEFT', 'baselineTM'; ...
        'STANCE_RIGHT', 'baselineTM'; ...
        'STANCE_LEFT', 'feedback'; ...
        'STANCE_RIGHT', 'feedback'; ...
        'STANCE_LEFT', 'postTM'; ...
        'STANCE_RIGHT', 'postTM'; ...
        'STANCE_LEFT', 'postOG'; ...
        'STANCE_RIGHT', 'postOG'; ...
      };
    number_of_conditions_to_analyze = size(conditions_to_analyze, 1);

    % subconcussion pilot
    condition_labels = {'stance foot', 'experimental'};
    conditions_to_analyze = ...
      {
        'STANCE_LEFT', 'MP'; ...
        'STANCE_RIGHT', 'MP'; ...
        'STANCE_LEFT', 'None'; ...
        'STANCE_RIGHT', 'None'; ...
      };
    number_of_conditions_to_analyze = size(conditions_to_analyze, 1);
    
%     % Obstacle
%     condition_labels = {'stance foot', 'experimental'};
%     conditions_to_analyze = ...
%       {
%         'STANCE_LEFT', 'NEAR'; ...
%         'STANCE_RIGHT', 'NEAR'; ...
%         'STANCE_BOTH', 'NEAR'; ...
%         'STANCE_LEFT', 'FAR'; ...
%         'STANCE_RIGHT', 'FAR'; ...
%         'STANCE_BOTH', 'FAR'; ...
%         'STANCE_LEFT', 'NO'; ...
%         'STANCE_RIGHT', 'NO'; ...
%         'STANCE_BOTH', 'NO'; ...
%       };
%     number_of_conditions_to_analyze = size(conditions_to_analyze, 1);
    
    number_of_time_steps_normalized = 100;
    perSecond_to_perMinute = 60;
    second_to_minute = 1/60;

    %% extract data
    load('subjectInfo.mat', 'date', 'subject_id');
    
    % initialize containers
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
        lheel_x_pos_com_normalized_all = [];
        rheel_x_pos_com_normalized_all = [];
        
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
        lcop_x_normalized_all = [];
        rcop_x_normalized_all = [];
        fxl_normalized_all = [];
        fzl_normalized_all = [];
        myl_normalized_all = [];
        fxr_normalized_all = [];
        fzr_normalized_all = [];
        myr_normalized_all = [];
    end
    if process_data_emg
        lglutmed_emg_normalized_all = [];
        ltibiant_emg_normalized_all = [];
        lperolng_emg_normalized_all = [];
        rglutmed_emg_normalized_all = [];
        rtibiant_emg_normalized_all = [];
        rperolng_emg_normalized_all = [];
    end
    if process_data_angles
        joint_angles_normalized_all = [];
        joint_velocities_normalized_all = [];
        joint_accelerations_normalized_all = [];
    end

    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            % load data
            condition = condition_list{i_condition};
            load(['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'relevantDataStretches')]);
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'kinematicTrajectories')]);
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

            % calculate dependent variables and initialize containers
            if process_data_balance
                % define markers and indices
                c7_marker = find(strcmp(marker_headers, 'C7'));
                lpsi_marker = find(strcmp(marker_headers, 'LPSI'));
                rpsi_marker = find(strcmp(marker_headers, 'RPSI'));
                lheel_marker = find(strcmp(marker_headers, 'LHEE'));
                rheel_marker = find(strcmp(marker_headers, 'RHEE'));
                com_body = find(strcmp(com_labels, 'BODYCOM'));

                c7_marker_indices = reshape([(c7_marker - 1) * 3 + 1; (c7_marker - 1) * 3 + 2; (c7_marker - 1) * 3 + 3], 1, length(c7_marker)*3);
                lpsi_marker_indices = reshape([(lpsi_marker - 1) * 3 + 1; (lpsi_marker - 1) * 3 + 2; (lpsi_marker - 1) * 3 + 3], 1, length(lpsi_marker)*3);
                rpsi_marker_indices = reshape([(rpsi_marker - 1) * 3 + 1; (rpsi_marker - 1) * 3 + 2; (rpsi_marker - 1) * 3 + 3], 1, length(rpsi_marker)*3);
                lheel_marker_indices = reshape([(lheel_marker - 1) * 3 + 1; (lheel_marker - 1) * 3 + 2; (lheel_marker - 1) * 3 + 3], 1, length(lheel_marker)*3);
                rheel_marker_indices = reshape([(rheel_marker - 1) * 3 + 1; (rheel_marker - 1) * 3 + 2; (rheel_marker - 1) * 3 + 3], 1, length(rheel_marker)*3);
                com_body_indices = reshape([(com_body - 1) * 3 + 1; (com_body - 1) * 3 + 2; (com_body - 1) * 3 + 3], 1, length(com_body) *3);

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
                com_x_pos_trajectory = com_trajectories(:, com_body_indices(1));
                com_y_pos_trajectory = com_trajectories(:, com_body_indices(2));
                com_z_pos_trajectory = com_trajectories(:, com_body_indices(3));
                
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
                lheel_x_pos_com_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rheel_x_pos_com_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lleg_angle_ml_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rleg_angle_ml_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                trunk_angle_ap_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                trunk_angle_ml_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                
                
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
                forceplate_file_name = ['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'forceplateTrajectories')];
                if exist(forceplate_file_name)
                    load(forceplate_file_name);
                    % initialize containers
                    force_plate_data_available = 1;
                else
                    force_plate_data_available = 0;
                end
                cop_x_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            end
            if process_data_emg
                load(makeFileName(date, subject_id, condition, i_trial, 'emgTrajectories'));
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
                    com_x_pos_extracted_stretch = com_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    com_y_pos_extracted_stretch = com_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    
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
                    elseif strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_BOTH')
                        stance_foot_heel_x_initial = 0;
                        stance_foot_heel_y_initial = 0;
                    else
                        error('stance condition should be either "STANCE_LEFT", "STANCE_RIGHT" or "STANCE_BOTH"');
                    end
                  
                    % calculate step width and length
                    if strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_RIGHT')
                        step_width_trial(i_stretch) = abs(lheel_x_pos_extracted_stretch(end) - rheel_x_pos_extracted_stretch(1));
                        step_length_trial(i_stretch) = lheel_y_pos_extracted_stretch(end) - rheel_y_pos_extracted_stretch(1);
                    elseif strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_LEFT')
                        step_width_trial(i_stretch) = abs(rheel_x_pos_extracted_stretch(end) - lheel_x_pos_extracted_stretch(1));
                        step_length_trial(i_stretch) = rheel_y_pos_extracted_stretch(end) - lheel_y_pos_extracted_stretch(1);
                    elseif strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_LEFT')
                        step_width_trial(i_stretch) = 0;
                        step_length_trial(i_stretch) = 0;
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
                    com_x_pos_normalized_stretch = spline(time_extracted_mocap, com_x_pos_extracted_stretch, time_normalized_mocap);
                    com_y_pos_normalized_stretch = spline(time_extracted_mocap, com_y_pos_extracted_stretch, time_normalized_mocap);
                    
                    trunk_angle_ap_normalized_stretch = spline(time_extracted_mocap, trunk_angle_ap_extracted_stretch, time_normalized_mocap);
                    trunk_angle_ml_normalized_stretch = spline(time_extracted_mocap, trunk_angle_ml_extracted_stretch, time_normalized_mocap);
                    lleg_angle_ml_normalized_stretch = spline(time_extracted_mocap, lleg_angle_ml_extracted_stretch, time_normalized_mocap);
                    rleg_angle_ml_normalized_stretch = spline(time_extracted_mocap, rleg_angle_ml_extracted_stretch, time_normalized_mocap);
                    
                    % store
                    com_x_pos_normalized_trial(:, i_stretch) = com_x_pos_normalized_stretch;
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
                    rheel_x_pos_mpsis_normalized_trial(:, i_stretch) = rheel_x_pos_normalized_stretch - mpsi_x_pos_normalized_stretch;
                    lheel_x_pos_com_normalized_trial(:, i_stretch) = lheel_x_pos_normalized_stretch - com_x_pos_normalized_stretch;
                    rheel_x_pos_com_normalized_trial(:, i_stretch) = rheel_x_pos_normalized_stretch - com_x_pos_normalized_stretch;

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
                    elseif strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_BOTH')
                        stance_foot_heel_x_initial = 0;
                        stance_foot_heel_y_initial = 0;
                    else
                        error('stance condition should be either "STANCE_LEFT", "STANCE_RIGHT" or "STANCE_BOTH"');
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
                        cop_x_extracted_stretch = total_forceplate_cop_Acw(start_index_forceplate : end_index_forceplate, 1);

                        % normalize
                        time_normalized_forceplate = linspace(time_extracted_forceplate(1), time_extracted_forceplate(end), number_of_time_steps_normalized);
                        cop_x_normalized_stretch = spline(time_extracted_forceplate, cop_x_extracted_stretch, time_normalized_forceplate);

                        % use stance foot heel as reference and store
                        cop_x_normalized_trial(:, i_stretch) = cop_x_normalized_stretch - stance_foot_heel_x_initial;
                    end
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
                lheel_x_pos_com_normalized_all = [lheel_x_pos_com_normalized_all lheel_x_pos_com_normalized_trial];
                rheel_x_pos_com_normalized_all = [rheel_x_pos_com_normalized_all rheel_x_pos_com_normalized_trial];
                
                trunk_angle_ap_normalized_all = [trunk_angle_ap_normalized_all trunk_angle_ap_normalized_trial];
                trunk_angle_ml_normalized_all = [trunk_angle_ml_normalized_all trunk_angle_ml_normalized_trial];
                lleg_angle_ml_normalized_all = [lleg_angle_ml_normalized_all lleg_angle_ml_normalized_trial];
                rleg_angle_ml_normalized_all = [rleg_angle_ml_normalized_all rleg_angle_ml_normalized_trial];
                
                
            end

            if process_data_armswing
                linclination_normalized_all = [linclination_normalized_all linclination_normalized_trial];
                rinclination_normalized_all = [rinclination_normalized_all rinclination_normalized_trial];
            end
                
            if process_data_forceplate
                cop_x_normalized_all = [cop_x_normalized_all cop_x_normalized_trial];
            end


            disp(['Condition ' condition ', Trial ' num2str(i_trial) ' done extracted ' num2str(length(step_times_trial)) ' stretches']);
        end
    end
    step_time_mean = mean(step_times_all);
    time_normalized = linspace(0, step_time_mean, number_of_time_steps_normalized);
    number_of_stretches = length(step_times_all);

    %% extract conditions
    % allocate indicators
    conditions_to_analyze_indicators = false(number_of_stretches, number_of_conditions_to_analyze);

    % extract indicators for conditions to analyze
    for i_condition = 1 : number_of_conditions_to_analyze
        stance_foot_indicator = strcmp(condition_stance_foot_list_all, conditions_to_analyze(i_condition, 1));
        experimental_indicator = strcmp(condition_experimental_list_all, conditions_to_analyze(i_condition, 2));
        
        this_condition_indicator = stance_foot_indicator & experimental_indicator;
        conditions_to_analyze_indicators(:, i_condition) = this_condition_indicator;
    end
    
    % give feedback about number of trials per condition
    trials_per_condition_to_analyze = sum(conditions_to_analyze_indicators)';
    conditions_to_analyze_with_number = conditions_to_analyze;
    for i_condition = 1 : number_of_conditions_to_analyze
        conditions_to_analyze_with_number{i_condition, length(condition_labels)+1} = num2str(trials_per_condition_to_analyze(i_condition));
    end
    conditions_to_analyze_with_labels = [condition_labels 'number of stretches'; conditions_to_analyze_with_number];
    
    disp('Conditions to analyze:')
    disp(conditions_to_analyze_with_labels);
    
    disp(['Number of unassigned stretches: ' num2str(number_of_stretches - sum(trials_per_condition_to_analyze))]);    

    %% calculate stats - old code
    if false

        if process_data_balance
            % absolute data
            lheel_x_pos_civ_stanceR = tinv(0.975, sum(conditions_stanceR)-1) * std(lheel_x_pos_normalized_all(:, conditions_stanceR), 1, 2)/sqrt(sum(conditions_stanceR));
            rheel_x_pos_civ_stanceL = tinv(0.975, sum(conditions_stanceL)-1) * std(rheel_x_pos_normalized_all(:, conditions_stanceL), 1, 2)/sqrt(sum(conditions_stanceL));

            step_width_left_mean = mean(step_width_all(conditions_stanceL));
            step_width_left_std = std(step_width_all(conditions_stanceL));
            step_width_left_civ = tinv(0.975, sum(conditions_stanceL)-1) * std(step_width_all(conditions_stanceL), 1, 2)/sqrt(sum(conditions_stanceL));
            step_width_right_mean = mean(step_width_all(conditions_stanceR));
            step_width_right_std = std(step_width_all(conditions_stanceR));
            step_width_right_civ = tinv(0.975, sum(conditions_stanceR)-1) * std(step_width_all(conditions_stanceR), 1, 2)/sqrt(sum(conditions_stanceR));

        end

        if process_data_forceplate
            % absolute data
            lcop_x_civ_stanceL = tinv(0.975, sum(conditions_stanceL)-1) * std(lcop_x_normalized_all(:, conditions_stanceL), 1, 2)/sqrt(sum(conditions_stanceL));
            rcop_x_civ_stanceR = tinv(0.975, sum(conditions_stanceR)-1) * std(rcop_x_normalized_all(:, conditions_stanceR), 1, 2)/sqrt(sum(conditions_stanceR));
            lcop_x_std_stanceL = std(lcop_x_normalized_all(:, conditions_stanceL), 1, 2);
            rcop_x_std_stanceR = std(rcop_x_normalized_all(:, conditions_stanceR), 1, 2);
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
            'condition_stance_foot_list_all', ...
            'condition_perturbation_list_all', ...
            'condition_delay_list_all', ...
            'condition_index_list_all', ...
            'condition_experimental_list_all', ...
            'conditions_to_analyze_indicators', ... % shouldn't be needed, just leave here for current version of plotGaitParameters to still work
            'condition_labels', ...
            'conditions_to_analyze', ...
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
                'lheel_x_pos_normalized_all', ...
                'rheel_x_pos_normalized_all', ...
                'lheel_y_pos_normalized_all', ...
                'rheel_y_pos_normalized_all', ...
                'lheel_x_pos_stancefoot_normalized_all', ...
                'rheel_x_pos_stancefoot_normalized_all', ...
                'lheel_x_pos_mpsis_normalized_all', ...
                'rheel_x_pos_mpsis_normalized_all', ...
                'lheel_x_pos_com_normalized_all',...
                'rheel_x_pos_com_normalized_all',...
                'trunk_angle_ap_normalized_all', ...
                'trunk_angle_ml_normalized_all', ...
                'lleg_angle_ml_normalized_all', ...
                'rleg_angle_ml_normalized_all', ...
                'step_speed_all' ...
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
                'cop_x_normalized_all' ...
              );
        end
    end
    
end




