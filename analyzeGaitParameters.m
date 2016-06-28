
extract_data                        = 1;
extract_conditions                  = 1;
calculate_responses                 = 1;
calculate_stats                     = 1;

save_data                           = 1;

process_data_marker                 = 1;
process_data_forceplate             = 1;
process_data_emg                    = 0;
process_data_angles                 = 0;
process_data_torques                = 0;

do_joint_angle_absolute_plots_left  = 0;
do_joint_angle_absolute_plots_right = 0;
do_joint_angle_response_plots       = 0;

do_acceleration_strategy_response_plots     = 0;

wait_times = [0 0.150 0.450];
wait_time_labels = {'0ms', '150ms', '450ms'};
load subjectInfo.mat;
% load(makeFileName(date, subject_id, 'model'));

trials_to_process = 5 : 16;
% trials_to_process = 9;
% trials_to_process = [1:11 13:14];


number_of_time_steps_normalized = 150;

experimental_conditions = {'ON', 'OFF', 'POST'};
experimental_condition_by_trial = ...
  { ...
    'PRE', ...
    'PRE', ...
    'PRE', ...
    'ON', ...
    'OFF', ...
    'OFF', ...
    'ON', ...
    'OFF', ...
    'OFF', ...
    'OFF', ...
    'ON', ...
    'ON', ...
    'ON', ...
    'POST', ...
    'POST', ...
    'POST' ...
  };


% experimental_condition_by_trial = ...
%   { ...
%     'NARROW', ...
%     'WIDE', ...
%     'NARROW', ...
%     'NARROW', ...
%     'NARROW', ...
%     'WIDE', ...
%     'WIDE', ...
%     'WIDE', ...
%     'WIDE', ...
%     'NARROW', ...
%     'WIDE', ...
%     'NARROW', ...
%     'NARROW', ...
%     'WIDE' ...
%   };

%% extract data
if extract_data
    % initialize containers
    stretch_length_indices_forceplate = [];
    condition_stance_foot_list_total = {};
    condition_experimental_list_total = {};
    origin_trial_list_total = [];
    origin_start_time_list_total = [];
    origin_end_time_list_total = [];
    step_times_total = [];
    if process_data_marker
        step_width_total = [];
        step_length_total = [];
        lasi_x_pos_normalized_total = [];
        rasi_x_pos_normalized_total = [];
        lpsi_x_pos_normalized_total = [];
        rpsi_x_pos_normalized_total = [];
        lasi_x_vel_normalized_total = [];
        rasi_x_vel_normalized_total = [];
        lpsi_x_vel_normalized_total = [];
        rpsi_x_vel_normalized_total = [];
        pelvis_x_pos_normalized_total = [];
        pelvis_x_vel_normalized_total = [];
        lheel_x_pos_normalized_total = [];
        rheel_x_pos_normalized_total = [];
        lheel_y_pos_normalized_total = [];
        rheel_y_pos_normalized_total = [];
    end
    if process_data_forceplate
        lcop_x_normalized_total = [];
        rcop_x_normalized_total = [];
        fxl_normalized_total = [];
        fzl_normalized_total = [];
        myl_normalized_total = [];
        fxr_normalized_total = [];
        fzr_normalized_total = [];
        myr_normalized_total = [];
    end
    if process_data_emg
        lglutmed_emg_normalized_total = [];
        ltibiant_emg_normalized_total = [];
        lperolng_emg_normalized_total = [];
        rglutmed_emg_normalized_total = [];
        rtibiant_emg_normalized_total = [];
        rperolng_emg_normalized_total = [];
    end
    if process_data_angles
        joint_angles_normalized_total = [];
        joint_velocities_normalized_total = [];
        joint_accelerations_normalized_total = [];
    end
    
    for i_trial = trials_to_process
        % load data
        load(makeFileName(date, subject_id, 'walking', i_trial, 'stepEvents'));
        load(makeFileName(date, subject_id, 'walking', i_trial, 'relevantDataStretches'));
        number_of_stretches_trial = length(condition_stance_foot_list);

        condition_stance_foot_list_trial = condition_stance_foot_list;
        condition_experimental_list_trial = cell(size(condition_stance_foot_list));
        condition_experimental_list_trial(1 : end) = experimental_condition_by_trial(i_trial);
        origin_trial_list_trial = zeros(number_of_stretches_trial, 1);
        origin_start_time_list_trial = zeros(number_of_stretches_trial, 1);
        origin_end_time_list_trial = zeros(number_of_stretches_trial, 1);
        step_times_trial = zeros(number_of_stretches_trial, 1);
        stim_start_time_relative_to_stretch_trial = stim_start_time_relative_to_stretch;

        if process_data_marker
            load(makeFileName(date, subject_id, 'walking', i_trial, 'markerTrajectories'));
            
            % define markers and indices
            lasi_marker = find(strcmp(marker_headers, 'LASI'));
            rasi_marker = find(strcmp(marker_headers, 'RASI'));
            lpsi_marker = find(strcmp(marker_headers, 'LPSI'));
            rpsi_marker = find(strcmp(marker_headers, 'RPSI'));
            lheel_marker = find(strcmp(marker_headers, 'LHEE'));
            rheel_marker = find(strcmp(marker_headers, 'RHEE'));

            lasi_marker_indices = reshape([(lasi_marker - 1) * 3 + 1; (lasi_marker - 1) * 3 + 2; (lasi_marker - 1) * 3 + 3], 1, length(lasi_marker)*3);
            rasi_marker_indices = reshape([(rasi_marker - 1) * 3 + 1; (rasi_marker - 1) * 3 + 2; (rasi_marker - 1) * 3 + 3], 1, length(rasi_marker)*3);
            lpsi_marker_indices = reshape([(lpsi_marker - 1) * 3 + 1; (lpsi_marker - 1) * 3 + 2; (lpsi_marker - 1) * 3 + 3], 1, length(lpsi_marker)*3);
            rpsi_marker_indices = reshape([(rpsi_marker - 1) * 3 + 1; (rpsi_marker - 1) * 3 + 2; (rpsi_marker - 1) * 3 + 3], 1, length(rpsi_marker)*3);
            lheel_marker_indices = reshape([(lheel_marker - 1) * 3 + 1; (lheel_marker - 1) * 3 + 2; (lheel_marker - 1) * 3 + 3], 1, length(lheel_marker)*3);
            rheel_marker_indices = reshape([(rheel_marker - 1) * 3 + 1; (rheel_marker - 1) * 3 + 2; (rheel_marker - 1) * 3 + 3], 1, length(rheel_marker)*3);


            % rename relevant trajectories
            lasi_x_pos_trajectory = marker_trajectories(:, lasi_marker_indices(1));
            rasi_x_pos_trajectory = marker_trajectories(:, rasi_marker_indices(1));
            lpsi_x_pos_trajectory = marker_trajectories(:, lpsi_marker_indices(1));
            rpsi_x_pos_trajectory = marker_trajectories(:, rpsi_marker_indices(1));
            lasi_x_vel_trajectory = deriveByTime(lasi_x_pos_trajectory, sampling_rate_mocap^(-1));
            rasi_x_vel_trajectory = deriveByTime(rasi_x_pos_trajectory, sampling_rate_mocap^(-1));
            lpsi_x_vel_trajectory = deriveByTime(lpsi_x_pos_trajectory, sampling_rate_mocap^(-1));
            rpsi_x_vel_trajectory = deriveByTime(rpsi_x_pos_trajectory, sampling_rate_mocap^(-1));
            lheel_x_pos_trajectory = marker_trajectories(:, lheel_marker_indices(1));
            rheel_x_pos_trajectory = marker_trajectories(:, rheel_marker_indices(1));
            lheel_y_pos_trajectory = marker_trajectories(:, lheel_marker_indices(2));
            rheel_y_pos_trajectory = marker_trajectories(:, rheel_marker_indices(2));

            % initialize containers
            step_width_trial = zeros(1, number_of_stretches_trial);
            step_length_trial = zeros(1, number_of_stretches_trial);
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
            lheel_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            rheel_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            lheel_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            rheel_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);

        end
        if process_data_forceplate
            load(makeFileName(date, subject_id, 'walking', i_trial, 'forceplateTrajectories'));

            % initialize containers
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
                lasi_x_pos_extracted_stretch = lasi_x_pos_trajectory(start_index_mocap : end_index_mocap);
                rasi_x_pos_extracted_stretch = rasi_x_pos_trajectory(start_index_mocap : end_index_mocap);
                lpsi_x_pos_extracted_stretch = lpsi_x_pos_trajectory(start_index_mocap : end_index_mocap);
                rpsi_x_pos_extracted_stretch = rpsi_x_pos_trajectory(start_index_mocap : end_index_mocap);
                lasi_x_vel_extracted_stretch = lasi_x_vel_trajectory(start_index_mocap : end_index_mocap);
                rasi_x_vel_extracted_stretch = rasi_x_vel_trajectory(start_index_mocap : end_index_mocap);
                lpsi_x_vel_extracted_stretch = lpsi_x_vel_trajectory(start_index_mocap : end_index_mocap);
                rpsi_x_vel_extracted_stretch = rpsi_x_vel_trajectory(start_index_mocap : end_index_mocap);
                lheel_x_pos_extracted_stretch = lheel_x_pos_trajectory(start_index_mocap : end_index_mocap);
                rheel_x_pos_extracted_stretch = rheel_x_pos_trajectory(start_index_mocap : end_index_mocap);
                lheel_y_pos_extracted_stretch = lheel_y_pos_trajectory(start_index_mocap : end_index_mocap);
                rheel_y_pos_extracted_stretch = rheel_y_pos_trajectory(start_index_mocap : end_index_mocap);

                % define stance foot heel as spatial point of reference
                if strcmp(condition_stance_foot_list_trial{i_stretch}, 'RIGHT')
                    stance_foot_heel_x_initial = rheel_x_pos_extracted_stretch(1);
                    stance_foot_heel_y_initial = rheel_y_pos_extracted_stretch(1);
                elseif strcmp(condition_stance_foot_list_trial{i_stretch}, 'LEFT')
                    stance_foot_heel_x_initial = lheel_x_pos_extracted_stretch(1);
                    stance_foot_heel_y_initial = lheel_y_pos_extracted_stretch(1);
                end

                % calculate step width and length
                if strcmp(condition_stance_foot_list_trial{i_stretch}, 'RIGHT')
                    step_width_trial(i_stretch) = lheel_x_pos_extracted_stretch(end) - rheel_x_pos_extracted_stretch(end);
                    step_length_trial(i_stretch) = lheel_y_pos_extracted_stretch(end) - rheel_y_pos_extracted_stretch(end);
                elseif strcmp(condition_stance_foot_list_trial{i_stretch}, 'LEFT')
                    step_width_trial(i_stretch) = rheel_x_pos_extracted_stretch(end) - lheel_x_pos_extracted_stretch(end);
                    step_length_trial(i_stretch) = rheel_y_pos_extracted_stretch(end) - lheel_y_pos_extracted_stretch(end);
                end

                % normalize mocap data in time
                time_normalized_mocap = linspace(time_extracted_mocap(1), time_extracted_mocap(end), number_of_time_steps_normalized);
                lasi_x_pos_normalized_stretch = spline(time_extracted_mocap, lasi_x_pos_extracted_stretch, time_normalized_mocap);
                rasi_x_pos_normalized_stretch = spline(time_extracted_mocap, rasi_x_pos_extracted_stretch, time_normalized_mocap);
                lpsi_x_pos_normalized_stretch = spline(time_extracted_mocap, lpsi_x_pos_extracted_stretch, time_normalized_mocap);
                rpsi_x_pos_normalized_stretch = spline(time_extracted_mocap, rpsi_x_pos_extracted_stretch, time_normalized_mocap);
                lasi_x_vel_normalized_stretch = spline(time_extracted_mocap, lasi_x_vel_extracted_stretch, time_normalized_mocap);
                rasi_x_vel_normalized_stretch = spline(time_extracted_mocap, rasi_x_vel_extracted_stretch, time_normalized_mocap);
                lpsi_x_vel_normalized_stretch = spline(time_extracted_mocap, lpsi_x_vel_extracted_stretch, time_normalized_mocap);
                rpsi_x_vel_normalized_stretch = spline(time_extracted_mocap, rpsi_x_vel_extracted_stretch, time_normalized_mocap);
                lheel_x_pos_normalized_stretch = spline(time_extracted_mocap, lheel_x_pos_extracted_stretch, time_normalized_mocap);
                rheel_x_pos_normalized_stretch = spline(time_extracted_mocap, rheel_x_pos_extracted_stretch, time_normalized_mocap);
                lheel_y_pos_normalized_stretch = spline(time_extracted_mocap, lheel_y_pos_extracted_stretch, time_normalized_mocap);
                rheel_y_pos_normalized_stretch = spline(time_extracted_mocap, rheel_y_pos_extracted_stretch, time_normalized_mocap);

                % use stance foot heel as reference and store
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
                lheel_x_pos_normalized_trial(:, i_stretch) = lheel_x_pos_normalized_stretch - stance_foot_heel_x_initial;
                rheel_x_pos_normalized_trial(:, i_stretch) = rheel_x_pos_normalized_stretch - stance_foot_heel_x_initial;
                lheel_y_pos_normalized_trial(:, i_stretch) = lheel_y_pos_normalized_stretch - stance_foot_heel_y_initial;
                rheel_y_pos_normalized_trial(:, i_stretch) = rheel_y_pos_normalized_stretch - stance_foot_heel_y_initial;
                
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
                lcop_x_extracted_stretch = left_forceplate_cop_Acw(start_index_forceplate : end_index_forceplate, 1);
                fxl_extracted_stretch = left_forceplate_wrench_Acw(start_index_forceplate : end_index_forceplate, 1);
                fzl_extracted_stretch = left_forceplate_wrench_Acw(start_index_forceplate : end_index_forceplate, 3);
                myl_extracted_stretch = left_forceplate_wrench_Acw(start_index_forceplate : end_index_forceplate, 5);
                fxr_extracted_stretch = right_forceplate_wrench_Acw(start_index_forceplate : end_index_forceplate, 1);
                fzr_extracted_stretch = right_forceplate_wrench_Acw(start_index_forceplate : end_index_forceplate, 3);
                myr_extracted_stretch = right_forceplate_wrench_Acw(start_index_forceplate : end_index_forceplate, 5);
                rcop_x_extracted_stretch = right_forceplate_cop_Acw(start_index_forceplate : end_index_forceplate, 1);
                time_extracted_forceplate = time_forceplate(start_index_forceplate : end_index_forceplate);
                
                % normalize
                time_normalized_forceplate = linspace(time_extracted_forceplate(1), time_extracted_forceplate(end), number_of_time_steps_normalized);
                fxl_normalized_stretch = spline(time_extracted_forceplate, fxl_extracted_stretch, time_normalized_forceplate);
                fzl_normalized_stretch = spline(time_extracted_forceplate, fzl_extracted_stretch, time_normalized_forceplate);
                myl_normalized_stretch = spline(time_extracted_forceplate, myl_extracted_stretch, time_normalized_forceplate);
                lcop_x_normalized_stretch = spline(time_extracted_forceplate, lcop_x_extracted_stretch, time_normalized_forceplate);
                fxr_normalized_stretch = spline(time_extracted_forceplate, fxr_extracted_stretch, time_normalized_forceplate);
                fzr_normalized_stretch = spline(time_extracted_forceplate, fzr_extracted_stretch, time_normalized_forceplate);
                myr_normalized_stretch = spline(time_extracted_forceplate, myr_extracted_stretch, time_normalized_forceplate);
                rcop_x_normalized_stretch = spline(time_extracted_forceplate, rcop_x_extracted_stretch, time_normalized_forceplate);

                % use stance foot heel as reference and store
                lcop_x_normalized_trial(:, i_stretch) = lcop_x_normalized_stretch - stance_foot_heel_x_initial;
                rcop_x_normalized_trial(:, i_stretch) = rcop_x_normalized_stretch - stance_foot_heel_x_initial;
%                 lcop_x_normalized_trial(:, i_stretch) = lcop_x_normalized_stretch;
%                 rcop_x_normalized_trial(:, i_stretch) = rcop_x_normalized_stretch;
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
        condition_stance_foot_list_total = [condition_stance_foot_list_total; condition_stance_foot_list_trial];
        condition_experimental_list_total = [condition_experimental_list_total; condition_experimental_list_trial];
        origin_trial_list_total = [origin_trial_list_total; origin_trial_list_trial];
        origin_start_time_list_total = [origin_start_time_list_total; origin_start_time_list_trial];
        origin_end_time_list_total = [origin_end_time_list_total; origin_end_time_list_trial];
        step_times_total = [step_times_total; step_times_trial];

        if process_data_marker
            step_width_total = [step_width_total step_width_trial];
            step_length_total = [step_length_total step_length_trial];
            lasi_x_pos_normalized_total = [lasi_x_pos_normalized_total lasi_x_pos_normalized_trial];
            rasi_x_pos_normalized_total = [rasi_x_pos_normalized_total rasi_x_pos_normalized_trial];
            lpsi_x_pos_normalized_total = [lpsi_x_pos_normalized_total lpsi_x_pos_normalized_trial];
            rpsi_x_pos_normalized_total = [rpsi_x_pos_normalized_total rpsi_x_pos_normalized_trial];
            lasi_x_vel_normalized_total = [lasi_x_vel_normalized_total lasi_x_vel_normalized_trial];
            rasi_x_vel_normalized_total = [rasi_x_vel_normalized_total rasi_x_vel_normalized_trial];
            lpsi_x_vel_normalized_total = [lpsi_x_vel_normalized_total lpsi_x_vel_normalized_trial];
            rpsi_x_vel_normalized_total = [rpsi_x_vel_normalized_total rpsi_x_vel_normalized_trial];
            pelvis_x_pos_normalized_total = [pelvis_x_pos_normalized_total pelvis_x_pos_normalized_trial];
            pelvis_x_vel_normalized_total = [pelvis_x_vel_normalized_total pelvis_x_vel_normalized_trial];
            lheel_x_pos_normalized_total = [lheel_x_pos_normalized_total lheel_x_pos_normalized_trial];
            rheel_x_pos_normalized_total = [rheel_x_pos_normalized_total rheel_x_pos_normalized_trial];
            lheel_y_pos_normalized_total = [lheel_y_pos_normalized_total lheel_y_pos_normalized_trial];
            rheel_y_pos_normalized_total = [rheel_y_pos_normalized_total rheel_y_pos_normalized_trial];
        end

        if process_data_forceplate
            fxl_normalized_total = [fxl_normalized_total fxl_normalized_trial];
            fzl_normalized_total = [fzl_normalized_total fxl_normalized_trial];
            myl_normalized_total = [myl_normalized_total fxl_normalized_trial];
            lcop_x_normalized_total = [lcop_x_normalized_total lcop_x_normalized_trial];
            fxr_normalized_total = [fxr_normalized_total fxr_normalized_trial];
            fzr_normalized_total = [fzr_normalized_total fxr_normalized_trial];
            myr_normalized_total = [myr_normalized_total fxr_normalized_trial];
            rcop_x_normalized_total = [rcop_x_normalized_total rcop_x_normalized_trial];
        end

        if process_data_angles
            joint_angles_normalized_total = [joint_angles_normalized_total joint_angles_normalized_trial];
            joint_velocities_normalized_total = [joint_velocities_normalized_total joint_velocities_normalized_trial];
            joint_accelerations_normalized_total = [joint_accelerations_normalized_total joint_accelerations_normalized_trial];
        end

        disp(['Trial ' num2str(i_trial) ' extracted']);
    end
    
    step_time_mean = mean(step_times_total);
    time_normalized = linspace(0, step_time_mean, number_of_time_steps_normalized);
end

%% extract conditions
if extract_conditions
    % extract conditions
    conditions_stanceL = strcmp(condition_stance_foot_list_total, 'LEFT');
    conditions_stanceR = strcmp(condition_stance_foot_list_total, 'RIGHT');
    

    conditions_all = conditions_stanceL | conditions_stanceR;
    
    disp(['Number of data points LEFT:       ' num2str(length(find(conditions_stanceL)))]);
    disp(['Number of data points RIGHT:      ' num2str(length(find(conditions_stanceR)))]);
    
    for i_condition = 1 : length(experimental_conditions)
        disp(['Number of data points ' experimental_conditions{i_condition} ': ' num2str(length(find(strcmp(condition_experimental_list_total, experimental_conditions{i_condition}))))]);
    end
    
end

%% calculate responses
if calculate_responses
    % calculate control means
    number_of_stretches = length(step_times_total);
    
    if process_data_marker
        % control means
        lheel_x_pos_mean_stanceL = mean(lheel_x_pos_normalized_total(:, conditions_stanceL), 2);
        lheel_x_pos_mean_stanceR = mean(lheel_x_pos_normalized_total(:, conditions_stanceR), 2);
        rheel_x_pos_mean_stanceL = mean(rheel_x_pos_normalized_total(:, conditions_stanceL), 2);
        rheel_x_pos_mean_stanceR = mean(rheel_x_pos_normalized_total(:, conditions_stanceR), 2);
        lheel_y_pos_mean_stanceL = mean(lheel_y_pos_normalized_total(:, conditions_stanceL), 2);
        lheel_y_pos_mean_stanceR = mean(lheel_y_pos_normalized_total(:, conditions_stanceR), 2);
        rheel_y_pos_mean_stanceL = mean(rheel_y_pos_normalized_total(:, conditions_stanceL), 2);
        rheel_y_pos_mean_stanceR = mean(rheel_y_pos_normalized_total(:, conditions_stanceR), 2);
        
        % responses
        lheel_x_pos_stanceR_response = lheel_x_pos_normalized_total - repmat(lheel_x_pos_mean_stanceR, 1, number_of_stretches);
        rheel_x_pos_stanceL_response = rheel_x_pos_normalized_total - repmat(rheel_x_pos_mean_stanceL, 1, number_of_stretches);
        lheel_y_pos_stanceR_response = lheel_y_pos_normalized_total - repmat(lheel_y_pos_mean_stanceR, 1, number_of_stretches);
        rheel_y_pos_stanceL_response = rheel_y_pos_normalized_total - repmat(rheel_y_pos_mean_stanceL, 1, number_of_stretches);
    end
    
    if process_data_forceplate
        % control means
        lcop_x_mean_stanceL = mean(lcop_x_normalized_total(:, conditions_stanceL), 2);
        rcop_x_mean_stanceR = mean(rcop_x_normalized_total(:, conditions_stanceR), 2);
        
        % response
        lcop_x_stanceL_response = lcop_x_normalized_total - repmat(lcop_x_mean_stanceL, 1, number_of_stretches);
        rcop_x_stanceR_response = rcop_x_normalized_total - repmat(rcop_x_mean_stanceR, 1, number_of_stretches);
    end
    
    if process_data_angles
        % control means
        joint_angles_mean_stanceL = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceL, :), 2));
        joint_angles_mean_stanceR = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceR, :), 2));
        joint_velocities_mean_stanceL = squeeze(mean(joint_velocities_normalized_total(:, conditions_stanceL, :), 2));
        joint_velocities_mean_stanceR = squeeze(mean(joint_velocities_normalized_total(:, conditions_stanceR, :), 2));
        joint_accelerations_mean_stanceL = squeeze(mean(joint_accelerations_normalized_total(:, conditions_stanceL, :), 2));
        joint_accelerations_mean_stanceR = squeeze(mean(joint_accelerations_normalized_total(:, conditions_stanceR, :), 2));
        
        % response
        joint_angles_stanceL_response = joint_angles_normalized_total - permute(repmat(joint_angles_mean_stanceL, 1, 1, number_of_stretches), [1 3 2]);
        joint_angles_stanceR_response = joint_angles_normalized_total - permute(repmat(joint_angles_mean_stanceR, 1, 1, number_of_stretches), [1 3 2]);
        joint_velocities_stanceL_response = joint_velocities_normalized_total - permute(repmat(joint_velocities_mean_stanceL, 1, 1, number_of_stretches), [1 3 2]);
        joint_velocities_stanceR_response = joint_velocities_normalized_total - permute(repmat(joint_velocities_mean_stanceR, 1, 1, number_of_stretches), [1 3 2]);
        joint_accelerations_stanceL_response = joint_accelerations_normalized_total - permute(repmat(joint_accelerations_mean_stanceL, 1, 1, number_of_stretches), [1 3 2]);
        joint_accelerations_stanceR_response = joint_accelerations_normalized_total - permute(repmat(joint_accelerations_mean_stanceR, 1, 1, number_of_stretches), [1 3 2]);
    end
end

%% calculate stats
if calculate_stats
    
    if process_data_marker
        % absolute data
        lheel_x_pos_civ_stanceR = tinv(0.975, sum(conditions_stanceR)-1) * std(lheel_x_pos_normalized_total(:, conditions_stanceR), 1, 2)/sqrt(sum(conditions_stanceR));
        rheel_x_pos_civ_stanceL = tinv(0.975, sum(conditions_stanceL)-1) * std(rheel_x_pos_normalized_total(:, conditions_stanceL), 1, 2)/sqrt(sum(conditions_stanceL));

        step_width_left_mean = mean(step_width_total(conditions_stanceL));
        step_width_left_std = std(step_width_total(conditions_stanceL));
        step_width_left_civ = tinv(0.975, sum(conditions_stanceL)-1) * std(step_width_total(conditions_stanceL), 1, 2)/sqrt(sum(conditions_stanceL));
        step_width_right_mean = mean(step_width_total(conditions_stanceR));
        step_width_right_std = std(step_width_total(conditions_stanceR));
        step_width_right_civ = tinv(0.975, sum(conditions_stanceR)-1) * std(step_width_total(conditions_stanceR), 1, 2)/sqrt(sum(conditions_stanceR));
        
    end

    if process_data_forceplate
        % absolute data
        lcop_x_civ_stanceL = tinv(0.975, sum(conditions_stanceL)-1) * std(lcop_x_normalized_total(:, conditions_stanceL), 1, 2)/sqrt(sum(conditions_stanceL));
        rcop_x_civ_stanceR = tinv(0.975, sum(conditions_stanceR)-1) * std(rcop_x_normalized_total(:, conditions_stanceR), 1, 2)/sqrt(sum(conditions_stanceR));
        lcop_x_std_stanceL = std(lcop_x_normalized_total(:, conditions_stanceL), 1, 2);
        rcop_x_std_stanceR = std(rcop_x_normalized_total(:, conditions_stanceR), 1, 2);
    end
    
    
end

%% save data
if save_data

    save ...
      ( ...
        makeFileName(date, subject_id, 'gaitParametersConditions'), ...
        'time_normalized', ...
        'origin_trial_list_total', ...
        'origin_start_time_list_total', ...
        'origin_end_time_list_total', ...
        'experimental_conditions', ...
        'condition_stance_foot_list_total', ...
        'condition_experimental_list_total', ...
        'conditions_stanceL', ...
        'conditions_stanceR' ...
      );
    
    if process_data_marker
        save ...
          ( ...
            makeFileName(date, subject_id, 'gaitParametersMarker'), ...
            'step_width_total', ...
            'step_length_total', ...
            'lasi_x_pos_normalized_total', ...
            'rasi_x_pos_normalized_total', ...
            'lpsi_x_pos_normalized_total', ...
            'rpsi_x_pos_normalized_total', ...
            'lasi_x_vel_normalized_total', ...
            'rasi_x_vel_normalized_total', ...
            'lpsi_x_vel_normalized_total', ...
            'rpsi_x_vel_normalized_total', ...
            'pelvis_x_pos_normalized_total', ...
            'pelvis_x_vel_normalized_total', ...
            'lheel_x_pos_normalized_total', ...
            'rheel_x_pos_normalized_total', ...
            'lheel_y_pos_normalized_total', ...
            'rheel_y_pos_normalized_total' ...
          );
    end
    
    if process_data_forceplate
        save ...
          ( ...
            makeFileName(date, subject_id, 'gaitParametersForceplate'), ...
            'lcop_x_normalized_total', ...
            'rcop_x_normalized_total', ...
            'rcop_x_mean_stanceR', ...
            'rcop_x_civ_stanceR', ...
            'rcop_x_std_stanceR', ...
            'lcop_x_mean_stanceL', ...
            'lcop_x_civ_stanceL', ...
            'lcop_x_std_stanceL' ...
          );
    end
end
    

return    




