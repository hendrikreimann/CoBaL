extract_data                        = 1;
extract_conditions                  = 1;
calculate_responses                 = 1;
calculate_strategy_directions       = 0;
calculate_strategy_responses        = 0;
calculate_stats                     = 0;

save_data                           = 1;

process_data_marker                 = 1;
process_data_forceplate             = 1;
process_data_emg                    = 0;
process_data_angles                 = 0;
process_data_torques                = 0;

do_acceleration_strategy_response_plots     = 0;

wait_times = [0 0.150 0.450];
wait_time_labels = {'0ms', '150ms', '450ms'};
load subjectInfo.mat;
% load(makeFileName(date, subject_id, 'model'));

% trials_to_process = 1 : 20;
trials_to_process = 3 : 43;
% trials_to_process = 3;
% trials_to_process = [2 4:21];


number_of_time_steps_normalized = 100;


%% extract data
if extract_data
    % initialize containers
    stretch_length_indices_forceplate = [];
    condition_stance_foot_list_total = {};
    condition_polarity_list_total = {};
    condition_delay_list_total = {};
    origin_trial_list_total = [];
    origin_start_time_list_total = [];
    origin_end_time_list_total = [];
    step_times_total = [];
    stim_start_time_relative_to_stretch_total = [];    
    if process_data_marker
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
        condition_polarity_list_trial = condition_polarity_list;
        condition_delay_list_trial = condition_delay_list;
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
                
                fxl_normalized_trial(:, i_stretch) = fxl_normalized_stretch;
                fzl_normalized_trial(:, i_stretch) = fzl_normalized_stretch;
                myl_normalized_trial(:, i_stretch) = myl_normalized_stretch;
                fxr_normalized_trial(:, i_stretch) = fxr_normalized_stretch;
                fzr_normalized_trial(:, i_stretch) = fzr_normalized_stretch;
                myr_normalized_trial(:, i_stretch) = myr_normalized_stretch;
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
        condition_polarity_list_total = [condition_polarity_list_total; condition_polarity_list_trial];
        condition_delay_list_total = [condition_delay_list_total; condition_delay_list_trial];
        origin_trial_list_total = [origin_trial_list_total; origin_trial_list_trial];
        origin_start_time_list_total = [origin_start_time_list_total; origin_start_time_list_trial];
        origin_end_time_list_total = [origin_end_time_list_total; origin_end_time_list_trial];
        step_times_total = [step_times_total; step_times_trial];
        stim_start_time_relative_to_stretch_total = [stim_start_time_relative_to_stretch_total; stim_start_time_relative_to_stretch_trial];

        if process_data_marker
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
            fzl_normalized_total = [fzl_normalized_total fzl_normalized_trial];
            myl_normalized_total = [myl_normalized_total myl_normalized_trial];
            lcop_x_normalized_total = [lcop_x_normalized_total lcop_x_normalized_trial];
            fxr_normalized_total = [fxr_normalized_total fxr_normalized_trial];
            fzr_normalized_total = [fzr_normalized_total fzr_normalized_trial];
            myr_normalized_total = [myr_normalized_total myr_normalized_trial];
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
    conditions_stanceL_stimNo = strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_polarity_list_total, 'CONTROL');
    conditions_stanceR_stimNo = strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_polarity_list_total, 'CONTROL');
    conditions_stanceL_stimPos_0ms = strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_polarity_list_total, 'POSITIVE') & strcmp(condition_delay_list_total, '0ms');
    conditions_stanceL_stimNeg_0ms = strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_polarity_list_total, 'NEGATIVE') & strcmp(condition_delay_list_total, '0ms');
    conditions_stanceR_stimPos_0ms = strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_polarity_list_total, 'POSITIVE') & strcmp(condition_delay_list_total, '0ms');
    conditions_stanceR_stimNeg_0ms = strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_polarity_list_total, 'NEGATIVE') & strcmp(condition_delay_list_total, '0ms');
    conditions_stanceL_stimPos_150ms = strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_polarity_list_total, 'POSITIVE') & strcmp(condition_delay_list_total, '150ms');
    conditions_stanceL_stimNeg_150ms = strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_polarity_list_total, 'NEGATIVE') & strcmp(condition_delay_list_total, '150ms');
    conditions_stanceR_stimPos_150ms = strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_polarity_list_total, 'POSITIVE') & strcmp(condition_delay_list_total, '150ms');
    conditions_stanceR_stimNeg_150ms = strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_polarity_list_total, 'NEGATIVE') & strcmp(condition_delay_list_total, '150ms');
    conditions_stanceL_stimPos_450ms = strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_polarity_list_total, 'POSITIVE') & strcmp(condition_delay_list_total, '450ms');
    conditions_stanceL_stimNeg_450ms = strcmp(condition_stance_foot_list_total, 'LEFT') & strcmp(condition_polarity_list_total, 'NEGATIVE') & strcmp(condition_delay_list_total, '450ms');
    conditions_stanceR_stimPos_450ms = strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_polarity_list_total, 'POSITIVE') & strcmp(condition_delay_list_total, '450ms');
    conditions_stanceR_stimNeg_450ms = strcmp(condition_stance_foot_list_total, 'RIGHT') & strcmp(condition_polarity_list_total, 'NEGATIVE') & strcmp(condition_delay_list_total, '450ms');

    conditions_stanceL_0ms = conditions_stanceL_stimPos_0ms | conditions_stanceL_stimNeg_0ms;
    conditions_stanceR_0ms = conditions_stanceR_stimPos_0ms | conditions_stanceR_stimNeg_0ms;
    conditions_stanceL_150ms = conditions_stanceL_stimPos_150ms | conditions_stanceL_stimNeg_150ms;
    conditions_stanceR_150ms = conditions_stanceR_stimPos_150ms | conditions_stanceR_stimNeg_150ms;
    conditions_stanceL_450ms = conditions_stanceL_stimPos_450ms | conditions_stanceL_stimNeg_450ms;
    conditions_stanceR_450ms = conditions_stanceR_stimPos_450ms | conditions_stanceR_stimNeg_450ms;
    
    conditions_all = conditions_stanceL_stimNo | conditions_stanceR_stimNo | conditions_stanceL_stimPos_0ms | conditions_stanceL_stimNeg_0ms | conditions_stanceR_stimPos_0ms | conditions_stanceR_stimNeg_0ms;
    
    disp(['Number of stretches LEFT 0ms:       ' num2str(length(find(conditions_stanceL_0ms)))]);
    disp(['Number of stretches RIGHT 0ms:      ' num2str(length(find(conditions_stanceR_0ms)))]);
    disp(['Number of stretches LEFT 150ms:     ' num2str(length(find(conditions_stanceL_150ms)))]);
    disp(['Number of stretches RIGHT 150ms:    ' num2str(length(find(conditions_stanceR_150ms)))]);
    disp(['Number of stretches LEFT 450ms:     ' num2str(length(find(conditions_stanceL_450ms)))]);
    disp(['Number of stretches RIGHT 450ms:    ' num2str(length(find(conditions_stanceR_450ms)))]);
    
end

%% calculate responses
if calculate_responses
    % calculate control means
    number_of_stretches = length(step_times_total);
    
    if process_data_marker
        % control means
        lheel_x_pos_mean_stanceL_control = mean(lheel_x_pos_normalized_total(:, conditions_stanceL_stimNo), 2);
        lheel_x_pos_mean_stanceR_control = mean(lheel_x_pos_normalized_total(:, conditions_stanceR_stimNo), 2);
        rheel_x_pos_mean_stanceL_control = mean(rheel_x_pos_normalized_total(:, conditions_stanceL_stimNo), 2);
        rheel_x_pos_mean_stanceR_control = mean(rheel_x_pos_normalized_total(:, conditions_stanceR_stimNo), 2);
        lheel_y_pos_mean_stanceL_control = mean(lheel_y_pos_normalized_total(:, conditions_stanceL_stimNo), 2);
        lheel_y_pos_mean_stanceR_control = mean(lheel_y_pos_normalized_total(:, conditions_stanceR_stimNo), 2);
        rheel_y_pos_mean_stanceL_control = mean(rheel_y_pos_normalized_total(:, conditions_stanceL_stimNo), 2);
        rheel_y_pos_mean_stanceR_control = mean(rheel_y_pos_normalized_total(:, conditions_stanceR_stimNo), 2);
        
        % responses
        lheel_x_pos_stanceR_response = lheel_x_pos_normalized_total - repmat(lheel_x_pos_mean_stanceR_control, 1, number_of_stretches);
        rheel_x_pos_stanceL_response = rheel_x_pos_normalized_total - repmat(rheel_x_pos_mean_stanceL_control, 1, number_of_stretches);
        lheel_y_pos_stanceR_response = lheel_y_pos_normalized_total - repmat(lheel_y_pos_mean_stanceR_control, 1, number_of_stretches);
        rheel_y_pos_stanceL_response = rheel_y_pos_normalized_total - repmat(rheel_y_pos_mean_stanceL_control, 1, number_of_stretches);
    end
    
    if process_data_forceplate
        % control means
        fxl_mean_stanceL_control = mean(fxl_normalized_total(:, conditions_stanceL_stimNo), 2);
        fzl_mean_stanceL_control = mean(fzl_normalized_total(:, conditions_stanceL_stimNo), 2);
        myl_mean_stanceL_control = mean(myl_normalized_total(:, conditions_stanceL_stimNo), 2);
        lcop_x_mean_stanceL_control = mean(lcop_x_normalized_total(:, conditions_stanceL_stimNo), 2);
        fxr_mean_stanceR_control = mean(fxr_normalized_total(:, conditions_stanceR_stimNo), 2);
        fzr_mean_stanceR_control = mean(fzr_normalized_total(:, conditions_stanceR_stimNo), 2);
        myr_mean_stanceR_control = mean(myr_normalized_total(:, conditions_stanceR_stimNo), 2);
        rcop_x_mean_stanceR_control = mean(rcop_x_normalized_total(:, conditions_stanceR_stimNo), 2);
        
        % response
        fxl_stanceL_response = fxl_normalized_total - repmat(fxl_mean_stanceL_control, 1, number_of_stretches);
        fzl_stanceL_response = fzl_normalized_total - repmat(fzl_mean_stanceL_control, 1, number_of_stretches);
        myl_stanceL_response = myl_normalized_total - repmat(myl_mean_stanceL_control, 1, number_of_stretches);
        lcop_x_stanceL_response = lcop_x_normalized_total - repmat(lcop_x_mean_stanceL_control, 1, number_of_stretches);
        fxr_stanceR_response = fxr_normalized_total - repmat(fxr_mean_stanceR_control, 1, number_of_stretches);
        fzr_stanceR_response = fzr_normalized_total - repmat(fzr_mean_stanceR_control, 1, number_of_stretches);
        myr_stanceR_response = myr_normalized_total - repmat(myr_mean_stanceR_control, 1, number_of_stretches);
        rcop_x_stanceR_response = rcop_x_normalized_total - repmat(rcop_x_mean_stanceR_control, 1, number_of_stretches);
    end
    
    if process_data_angles
        % control means
        joint_angles_mean_stanceL_control = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceL_stimNo, :), 2));
        joint_angles_mean_stanceR_control = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceR_stimNo, :), 2));
        joint_velocities_mean_stanceL_control = squeeze(mean(joint_velocities_normalized_total(:, conditions_stanceL_stimNo, :), 2));
        joint_velocities_mean_stanceR_control = squeeze(mean(joint_velocities_normalized_total(:, conditions_stanceR_stimNo, :), 2));
        joint_accelerations_mean_stanceL_control = squeeze(mean(joint_accelerations_normalized_total(:, conditions_stanceL_stimNo, :), 2));
        joint_accelerations_mean_stanceR_control = squeeze(mean(joint_accelerations_normalized_total(:, conditions_stanceR_stimNo, :), 2));
        
        % response
        joint_angles_stanceL_response = joint_angles_normalized_total - permute(repmat(joint_angles_mean_stanceL_control, 1, 1, number_of_stretches), [1 3 2]);
        joint_angles_stanceR_response = joint_angles_normalized_total - permute(repmat(joint_angles_mean_stanceR_control, 1, 1, number_of_stretches), [1 3 2]);
        joint_velocities_stanceL_response = joint_velocities_normalized_total - permute(repmat(joint_velocities_mean_stanceL_control, 1, 1, number_of_stretches), [1 3 2]);
        joint_velocities_stanceR_response = joint_velocities_normalized_total - permute(repmat(joint_velocities_mean_stanceR_control, 1, 1, number_of_stretches), [1 3 2]);
        joint_accelerations_stanceL_response = joint_accelerations_normalized_total - permute(repmat(joint_accelerations_mean_stanceL_control, 1, 1, number_of_stretches), [1 3 2]);
        joint_accelerations_stanceR_response = joint_accelerations_normalized_total - permute(repmat(joint_accelerations_mean_stanceR_control, 1, 1, number_of_stretches), [1 3 2]);
        
%         i_joint = 11;
%         plot(squeeze(joint_angles_normalized_total(:, conditions_stanceL_control, i_joint)))
%         plot(squeeze(joint_angles_normalized_total(:, conditions_stanceR_control, i_joint)))
%         plot(squeeze(joint_velocities_normalized_total(:, conditions_stanceL_control, i_joint)))
%         plot(squeeze(joint_velocities_normalized_total(:, conditions_stanceR_control, i_joint)))
%         plot(squeeze(joint_accelerations_normalized_total(:, conditions_stanceL_control, i_joint)))
%         plot(squeeze(joint_accelerations_normalized_total(:, conditions_stanceR_control, i_joint)))
  
%         i_joint = 11;
%         plot(mean(squeeze(joint_angles_stanceL_response(:, conditions_stanceL_control, i_joint)), 2))
%         plot(mean(squeeze(joint_angles_right_response(:, conditions_stanceR_control, i_joint)), 2))
%         plot(mean(squeeze(joint_velocities_stanceL_response(:, conditions_stanceL_control, i_joint)), 2))
%         plot(mean(squeeze(joint_velocities_right_response(:, conditions_stanceR_control, i_joint)), 2))
%         plot(mean(squeeze(joint_accelerations_stanceL_response(:, conditions_stanceL_control, i_joint)), 2))
%         plot(mean(squeeze(joint_accelerations_right_response(:, conditions_stanceR_control, i_joint)), 2))

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
if calculate_stats
    
    if process_data_marker
        % absolute data
        lheel_x_pos_civ_stanceR_control = tinv(0.975, sum(conditions_stanceR_stimNo)-1) * std(lheel_x_pos_normalized_total(:, conditions_stanceR_stimNo), 1, 2)/sqrt(sum(conditions_stanceR_stimNo));
        lheel_x_pos_mean_stanceR_stimPos_0ms = mean(lheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_0ms), 2);
        lheel_x_pos_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(lheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
        lheel_x_pos_mean_stanceR_stimNeg_0ms = mean(lheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2);
        lheel_x_pos_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(lheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
        lheel_x_pos_mean_stanceR_stimPos_150ms = mean(lheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_150ms), 2);
        lheel_x_pos_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(lheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
        lheel_x_pos_mean_stanceR_stimNeg_150ms = mean(lheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2);
        lheel_x_pos_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(lheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
        lheel_x_pos_mean_stanceR_stimPos_450ms = mean(lheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_450ms), 2);
        lheel_x_pos_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(lheel_x_pos_normalized_total(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
        lheel_x_pos_mean_stanceR_stimNeg_450ms = mean(lheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2);
        lheel_x_pos_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(lheel_x_pos_normalized_total(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));

        rheel_x_pos_civ_stanceL_control = tinv(0.975, sum(conditions_stanceL_stimNo)-1) * std(rheel_x_pos_normalized_total(:, conditions_stanceL_stimNo), 1, 2)/sqrt(sum(conditions_stanceL_stimNo));
        rheel_x_pos_mean_stanceL_stimPos_0ms = mean(rheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_0ms), 2);
        rheel_x_pos_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(rheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
        rheel_x_pos_mean_stanceL_stimNeg_0ms = mean(rheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2);
        rheel_x_pos_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(rheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
        rheel_x_pos_mean_stanceL_stimPos_150ms = mean(rheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_150ms), 2);
        rheel_x_pos_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(rheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
        rheel_x_pos_mean_stanceL_stimNeg_150ms = mean(rheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2);
        rheel_x_pos_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(rheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
        rheel_x_pos_mean_stanceL_stimPos_450ms = mean(rheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_450ms), 2);
        rheel_x_pos_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(rheel_x_pos_normalized_total(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
        rheel_x_pos_mean_stanceL_stimNeg_450ms = mean(rheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2);
        rheel_x_pos_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(rheel_x_pos_normalized_total(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));
        
        % responses
        lheel_x_pos_response_mean_stanceR_stimPos_0ms = mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2);
        lheel_x_pos_response_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
        lheel_x_pos_response_mean_stanceR_stimNeg_0ms = mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2);
        lheel_x_pos_response_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
        lheel_x_pos_response_mean_stanceR_stimPos_150ms = mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2);
        lheel_x_pos_response_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
        lheel_x_pos_response_mean_stanceR_stimNeg_150ms = mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2);
        lheel_x_pos_response_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
        lheel_x_pos_response_mean_stanceR_stimPos_450ms = mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2);
        lheel_x_pos_response_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
        lheel_x_pos_response_mean_stanceR_stimNeg_450ms = mean(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2);
        lheel_x_pos_response_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(lheel_x_pos_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));

        rheel_x_pos_response_mean_stanceL_stimPos_0ms = mean(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2);
        rheel_x_pos_response_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
        rheel_x_pos_response_mean_stanceL_stimNeg_0ms = mean(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2);
        rheel_x_pos_response_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
        rheel_x_pos_response_mean_stanceL_stimPos_150ms = mean(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2);
        rheel_x_pos_response_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
        rheel_x_pos_response_mean_stanceL_stimNeg_150ms = mean(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2);
        rheel_x_pos_response_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
        rheel_x_pos_response_mean_stanceL_stimPos_450ms = mean(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2);
        rheel_x_pos_response_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
        rheel_x_pos_response_mean_stanceL_stimNeg_450ms = mean(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2);
        rheel_x_pos_response_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(rheel_x_pos_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));

        lheel_y_pos_response_mean_stanceR_stimPos_0ms = mean(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2);
        lheel_y_pos_response_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
        lheel_y_pos_response_mean_stanceR_stimNeg_0ms = mean(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2);
        lheel_y_pos_response_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
        lheel_y_pos_response_mean_stanceR_stimPos_150ms = mean(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2);
        lheel_y_pos_response_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
        lheel_y_pos_response_mean_stanceR_stimNeg_150ms = mean(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2);
        lheel_y_pos_response_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
        lheel_y_pos_response_mean_stanceR_stimPos_450ms = mean(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2);
        lheel_y_pos_response_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
        lheel_y_pos_response_mean_stanceR_stimNeg_450ms = mean(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2);
        lheel_y_pos_response_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(lheel_y_pos_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));

        rheel_y_pos_response_mean_stanceL_stimPos_0ms = mean(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2);
        rheel_y_pos_response_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
        rheel_y_pos_response_mean_stanceL_stimNeg_0ms = mean(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2);
        rheel_y_pos_response_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
        rheel_y_pos_response_mean_stanceL_stimPos_150ms = mean(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2);
        rheel_y_pos_response_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
        rheel_y_pos_response_mean_stanceL_stimNeg_150ms = mean(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2);
        rheel_y_pos_response_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
        rheel_y_pos_response_mean_stanceL_stimPos_450ms = mean(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2);
        rheel_y_pos_response_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
        rheel_y_pos_response_mean_stanceL_stimNeg_450ms = mean(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2);
        rheel_y_pos_response_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(rheel_y_pos_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));

    end

    if process_data_forceplate
        % absolute data
        fxl_civ_stanceL_control = tinv(0.975, sum(conditions_stanceL_stimNo)-1) * std(fxl_normalized_total(:, conditions_stanceL_stimNo), 1, 2)/sqrt(sum(conditions_stanceL_stimNo));
        fxl_mean_stanceL_stimPos_0ms = mean(fxl_normalized_total(:, conditions_stanceL_stimPos_0ms), 2);
        fxl_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(fxl_normalized_total(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
        fxl_mean_stanceL_stimNeg_0ms = mean(fxl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2);
        fxl_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(fxl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
        fxl_mean_stanceL_stimPos_150ms = mean(fxl_normalized_total(:, conditions_stanceL_stimPos_150ms), 2);
        fxl_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(fxl_normalized_total(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
        fxl_mean_stanceL_stimNeg_150ms = mean(fxl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2);
        fxl_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(fxl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
        fxl_mean_stanceL_stimPos_450ms = mean(fxl_normalized_total(:, conditions_stanceL_stimPos_450ms), 2);
        fxl_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(fxl_normalized_total(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
        fxl_mean_stanceL_stimNeg_450ms = mean(fxl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2);
        fxl_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(fxl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));
        
        fzl_civ_stanceL_control = tinv(0.975, sum(conditions_stanceL_stimNo)-1) * std(fzl_normalized_total(:, conditions_stanceL_stimNo), 1, 2)/sqrt(sum(conditions_stanceL_stimNo));
        fzl_mean_stanceL_stimPos_0ms = mean(fzl_normalized_total(:, conditions_stanceL_stimPos_0ms), 2);
        fzl_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(fzl_normalized_total(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
        fzl_mean_stanceL_stimNeg_0ms = mean(fzl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2);
        fzl_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(fzl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
        fzl_mean_stanceL_stimPos_150ms = mean(fzl_normalized_total(:, conditions_stanceL_stimPos_150ms), 2);
        fzl_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(fzl_normalized_total(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
        fzl_mean_stanceL_stimNeg_150ms = mean(fzl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2);
        fzl_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(fzl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
        fzl_mean_stanceL_stimPos_450ms = mean(fzl_normalized_total(:, conditions_stanceL_stimPos_450ms), 2);
        fzl_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(fzl_normalized_total(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
        fzl_mean_stanceL_stimNeg_450ms = mean(fzl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2);
        fzl_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(fzl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));

        myl_civ_stanceL_control = tinv(0.975, sum(conditions_stanceL_stimNo)-1) * std(myl_normalized_total(:, conditions_stanceL_stimNo), 1, 2)/sqrt(sum(conditions_stanceL_stimNo));
        myl_mean_stanceL_stimPos_0ms = mean(myl_normalized_total(:, conditions_stanceL_stimPos_0ms), 2);
        myl_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(myl_normalized_total(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
        myl_mean_stanceL_stimNeg_0ms = mean(myl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2);
        myl_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(myl_normalized_total(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
        myl_mean_stanceL_stimPos_150ms = mean(myl_normalized_total(:, conditions_stanceL_stimPos_150ms), 2);
        myl_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(myl_normalized_total(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
        myl_mean_stanceL_stimNeg_150ms = mean(myl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2);
        myl_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(myl_normalized_total(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
        myl_mean_stanceL_stimPos_450ms = mean(myl_normalized_total(:, conditions_stanceL_stimPos_450ms), 2);
        myl_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(myl_normalized_total(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
        myl_mean_stanceL_stimNeg_450ms = mean(myl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2);
        myl_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(myl_normalized_total(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));
        
        lcop_x_civ_stanceL_control = tinv(0.975, sum(conditions_stanceL_stimNo)-1) * std(lcop_x_normalized_total(:, conditions_stanceL_stimNo), 1, 2)/sqrt(sum(conditions_stanceL_stimNo));
        lcop_x_mean_stanceL_stimPos_0ms = mean(lcop_x_normalized_total(:, conditions_stanceL_stimPos_0ms), 2);
        lcop_x_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(lcop_x_normalized_total(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
        lcop_x_mean_stanceL_stimNeg_0ms = mean(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_0ms), 2);
        lcop_x_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
        lcop_x_mean_stanceL_stimPos_150ms = mean(lcop_x_normalized_total(:, conditions_stanceL_stimPos_150ms), 2);
        lcop_x_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(lcop_x_normalized_total(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
        lcop_x_mean_stanceL_stimNeg_150ms = mean(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_150ms), 2);
        lcop_x_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
        lcop_x_mean_stanceL_stimPos_450ms = mean(lcop_x_normalized_total(:, conditions_stanceL_stimPos_450ms), 2);
        lcop_x_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(lcop_x_normalized_total(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
        lcop_x_mean_stanceL_stimNeg_450ms = mean(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_450ms), 2);
        lcop_x_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(lcop_x_normalized_total(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));

        fxr_civ_stanceR_control = tinv(0.975, sum(conditions_stanceR_stimNo)-1) * std(fxr_normalized_total(:, conditions_stanceR_stimNo), 1, 2)/sqrt(sum(conditions_stanceR_stimNo));
        fxr_mean_stanceR_stimPos_0ms = mean(fxr_normalized_total(:, conditions_stanceR_stimPos_0ms), 2);
        fxr_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(fxr_normalized_total(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
        fxr_mean_stanceR_stimNeg_0ms = mean(fxr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2);
        fxr_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(fxr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
        fxr_mean_stanceR_stimPos_150ms = mean(fxr_normalized_total(:, conditions_stanceR_stimPos_150ms), 2);
        fxr_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(fxr_normalized_total(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
        fxr_mean_stanceR_stimNeg_150ms = mean(fxr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2);
        fxr_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(fxr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
        fxr_mean_stanceR_stimPos_450ms = mean(fxr_normalized_total(:, conditions_stanceR_stimPos_450ms), 2);
        fxr_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(fxr_normalized_total(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
        fxr_mean_stanceR_stimNeg_450ms = mean(fxr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2);
        fxr_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(fxr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));

        fzr_civ_stanceR_control = tinv(0.975, sum(conditions_stanceR_stimNo)-1) * std(fzr_normalized_total(:, conditions_stanceR_stimNo), 1, 2)/sqrt(sum(conditions_stanceR_stimNo));
        fzr_mean_stanceR_stimPos_0ms = mean(fzr_normalized_total(:, conditions_stanceR_stimPos_0ms), 2);
        fzr_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(fzr_normalized_total(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
        fzr_mean_stanceR_stimNeg_0ms = mean(fzr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2);
        fzr_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(fzr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
        fzr_mean_stanceR_stimPos_150ms = mean(fzr_normalized_total(:, conditions_stanceR_stimPos_150ms), 2);
        fzr_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(fzr_normalized_total(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
        fzr_mean_stanceR_stimNeg_150ms = mean(fzr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2);
        fzr_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(fzr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
        fzr_mean_stanceR_stimPos_450ms = mean(fzr_normalized_total(:, conditions_stanceR_stimPos_450ms), 2);
        fzr_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(fzr_normalized_total(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
        fzr_mean_stanceR_stimNeg_450ms = mean(fzr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2);
        fzr_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(fzr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));

        myr_civ_stanceR_control = tinv(0.975, sum(conditions_stanceR_stimNo)-1) * std(myr_normalized_total(:, conditions_stanceR_stimNo), 1, 2)/sqrt(sum(conditions_stanceR_stimNo));
        myr_mean_stanceR_stimPos_0ms = mean(myr_normalized_total(:, conditions_stanceR_stimPos_0ms), 2);
        myr_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(myr_normalized_total(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
        myr_mean_stanceR_stimNeg_0ms = mean(myr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2);
        myr_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(myr_normalized_total(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
        myr_mean_stanceR_stimPos_150ms = mean(myr_normalized_total(:, conditions_stanceR_stimPos_150ms), 2);
        myr_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(myr_normalized_total(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
        myr_mean_stanceR_stimNeg_150ms = mean(myr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2);
        myr_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(myr_normalized_total(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
        myr_mean_stanceR_stimPos_450ms = mean(myr_normalized_total(:, conditions_stanceR_stimPos_450ms), 2);
        myr_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(myr_normalized_total(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
        myr_mean_stanceR_stimNeg_450ms = mean(myr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2);
        myr_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(myr_normalized_total(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));

        rcop_x_civ_stanceR_control = tinv(0.975, sum(conditions_stanceR_stimNo)-1) * std(rcop_x_normalized_total(:, conditions_stanceR_stimNo), 1, 2)/sqrt(sum(conditions_stanceR_stimNo));
        rcop_x_mean_stanceR_stimPos_0ms = mean(rcop_x_normalized_total(:, conditions_stanceR_stimPos_0ms), 2);
        rcop_x_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(rcop_x_normalized_total(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
        rcop_x_mean_stanceR_stimNeg_0ms = mean(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_0ms), 2);
        rcop_x_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
        rcop_x_mean_stanceR_stimPos_150ms = mean(rcop_x_normalized_total(:, conditions_stanceR_stimPos_150ms), 2);
        rcop_x_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(rcop_x_normalized_total(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
        rcop_x_mean_stanceR_stimNeg_150ms = mean(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_150ms), 2);
        rcop_x_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
        rcop_x_mean_stanceR_stimPos_450ms = mean(rcop_x_normalized_total(:, conditions_stanceR_stimPos_450ms), 2);
        rcop_x_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(rcop_x_normalized_total(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
        rcop_x_mean_stanceR_stimNeg_450ms = mean(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_450ms), 2);
        rcop_x_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(rcop_x_normalized_total(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));

        % responses
        fxl_response_mean_stanceL_stimPos_0ms = mean(fxl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2);
        fxl_response_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(fxl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
        fxl_response_mean_stanceL_stimNeg_0ms = mean(fxl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2);
        fxl_response_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(fxl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
        fxl_response_mean_stanceL_stimPos_150ms = mean(fxl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2);
        fxl_response_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(fxl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
        fxl_response_mean_stanceL_stimNeg_150ms = mean(fxl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2);
        fxl_response_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(fxl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
        fxl_response_mean_stanceL_stimPos_450ms = mean(fxl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2);
        fxl_response_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(fxl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
        fxl_response_mean_stanceL_stimNeg_450ms = mean(fxl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2);
        fxl_response_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(fxl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));

        fzl_response_mean_stanceL_stimPos_0ms = mean(fzl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2);
        fzl_response_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(fzl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
        fzl_response_mean_stanceL_stimNeg_0ms = mean(fzl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2);
        fzl_response_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(fzl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
        fzl_response_mean_stanceL_stimPos_150ms = mean(fzl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2);
        fzl_response_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(fzl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
        fzl_response_mean_stanceL_stimNeg_150ms = mean(fzl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2);
        fzl_response_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(fzl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
        fzl_response_mean_stanceL_stimPos_450ms = mean(fzl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2);
        fzl_response_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(fzl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
        fzl_response_mean_stanceL_stimNeg_450ms = mean(fzl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2);
        fzl_response_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(fzl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));

        myl_response_mean_stanceL_stimPos_0ms = mean(myl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2);
        myl_response_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(myl_stanceL_response(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
        myl_response_mean_stanceL_stimNeg_0ms = mean(myl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2);
        myl_response_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(myl_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
        myl_response_mean_stanceL_stimPos_150ms = mean(myl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2);
        myl_response_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(myl_stanceL_response(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
        myl_response_mean_stanceL_stimNeg_150ms = mean(myl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2);
        myl_response_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(myl_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
        myl_response_mean_stanceL_stimPos_450ms = mean(myl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2);
        myl_response_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(myl_stanceL_response(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
        myl_response_mean_stanceL_stimNeg_450ms = mean(myl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2);
        myl_response_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(myl_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));

        lcop_x_response_mean_stanceL_stimPos_0ms = mean(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_0ms), 2);
        lcop_x_response_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * std(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_0ms));
        lcop_x_response_mean_stanceL_stimNeg_0ms = mean(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 2);
        lcop_x_response_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * std(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_0ms));
        lcop_x_response_mean_stanceL_stimPos_150ms = mean(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_150ms), 2);
        lcop_x_response_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * std(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_150ms));
        lcop_x_response_mean_stanceL_stimNeg_150ms = mean(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 2);
        lcop_x_response_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * std(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_150ms));
        lcop_x_response_mean_stanceL_stimPos_450ms = mean(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_450ms), 2);
        lcop_x_response_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * std(lcop_x_stanceL_response(:, conditions_stanceL_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimPos_450ms));
        lcop_x_response_mean_stanceL_stimNeg_450ms = mean(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 2);
        lcop_x_response_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * std(lcop_x_stanceL_response(:, conditions_stanceL_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceL_stimNeg_450ms));

        fxr_response_mean_stanceR_stimPos_0ms = mean(fxr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2);
        fxr_response_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(fxr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
        fxr_response_mean_stanceR_stimNeg_0ms = mean(fxr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2);
        fxr_response_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(fxr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
        fxr_response_mean_stanceR_stimPos_150ms = mean(fxr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2);
        fxr_response_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(fxr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
        fxr_response_mean_stanceR_stimNeg_150ms = mean(fxr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2);
        fxr_response_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(fxr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
        fxr_response_mean_stanceR_stimPos_450ms = mean(fxr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2);
        fxr_response_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(fxr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
        fxr_response_mean_stanceR_stimNeg_450ms = mean(fxr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2);
        fxr_response_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(fxr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));

        fzr_response_mean_stanceR_stimPos_0ms = mean(fzr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2);
        fzr_response_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(fzr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
        fzr_response_mean_stanceR_stimNeg_0ms = mean(fzr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2);
        fzr_response_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(fzr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
        fzr_response_mean_stanceR_stimPos_150ms = mean(fzr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2);
        fzr_response_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(fzr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
        fzr_response_mean_stanceR_stimNeg_150ms = mean(fzr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2);
        fzr_response_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(fzr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
        fzr_response_mean_stanceR_stimPos_450ms = mean(fzr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2);
        fzr_response_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(fzr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
        fzr_response_mean_stanceR_stimNeg_450ms = mean(fzr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2);
        fzr_response_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(fzr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));

        myr_response_mean_stanceR_stimPos_0ms = mean(myr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 2);
        myr_response_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * std(myr_stanceR_response(:, conditions_stanceR_stimPos_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_0ms));
        myr_response_mean_stanceR_stimNeg_0ms = mean(myr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 2);
        myr_response_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * std(myr_stanceR_response(:, conditions_stanceR_stimNeg_0ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_0ms));
        myr_response_mean_stanceR_stimPos_150ms = mean(myr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 2);
        myr_response_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * std(myr_stanceR_response(:, conditions_stanceR_stimPos_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_150ms));
        myr_response_mean_stanceR_stimNeg_150ms = mean(myr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 2);
        myr_response_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * std(myr_stanceR_response(:, conditions_stanceR_stimNeg_150ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_150ms));
        myr_response_mean_stanceR_stimPos_450ms = mean(myr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 2);
        myr_response_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * std(myr_stanceR_response(:, conditions_stanceR_stimPos_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimPos_450ms));
        myr_response_mean_stanceR_stimNeg_450ms = mean(myr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 2);
        myr_response_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * std(myr_stanceR_response(:, conditions_stanceR_stimNeg_450ms), 1, 2)/sqrt(sum(conditions_stanceR_stimNeg_450ms));
    end
    
    if process_data_angles
        % absolute data
        joint_angles_civ_stanceL_control = tinv(0.975, sum(conditions_stanceL_stimNo)-1) * squeeze(std(joint_angles_normalized_total(:, conditions_stanceL_stimNo, :), 1, 2))/sqrt(sum(conditions_stanceL_stimNo));
        joint_angles_mean_stanceL_stimPos_0ms = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceL_stimPos_0ms, :), 2));
        joint_angles_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * squeeze(std(joint_angles_normalized_total(:, conditions_stanceL_stimPos_0ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimPos_0ms));
        joint_angles_mean_stanceL_stimNeg_0ms = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceL_stimNeg_0ms, :), 2));
        joint_angles_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * squeeze(std(joint_angles_normalized_total(:, conditions_stanceL_stimNeg_0ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimNeg_0ms));
        joint_angles_mean_stanceL_stimPos_150ms = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceL_stimPos_150ms, :), 2));
        joint_angles_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * squeeze(std(joint_angles_normalized_total(:, conditions_stanceL_stimPos_150ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimPos_150ms));
        joint_angles_mean_stanceL_stimNeg_150ms = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceL_stimNeg_150ms, :), 2));
        joint_angles_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * squeeze(std(joint_angles_normalized_total(:, conditions_stanceL_stimNeg_150ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimNeg_150ms));
        joint_angles_mean_stanceL_stimPos_450ms = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceL_stimPos_450ms, :), 2));
        joint_angles_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * squeeze(std(joint_angles_normalized_total(:, conditions_stanceL_stimPos_450ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimPos_450ms));
        joint_angles_mean_stanceL_stimNeg_450ms = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceL_stimNeg_450ms, :), 2));
        joint_angles_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * squeeze(std(joint_angles_normalized_total(:, conditions_stanceL_stimNeg_450ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimNeg_450ms));

        joint_angles_civ_stanceR_control = tinv(0.975, sum(conditions_stanceR_stimNo)-1) * squeeze(std(joint_angles_normalized_total(:, conditions_stanceR_stimNo, :), 1, 2))/sqrt(sum(conditions_stanceR_stimNo));
        joint_angles_mean_stanceR_stimPos_0ms = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceR_stimPos_0ms, :), 2));
        joint_angles_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * squeeze(std(joint_angles_normalized_total(:, conditions_stanceR_stimPos_0ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimPos_0ms));
        joint_angles_mean_stanceR_stimNeg_0ms = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceR_stimNeg_0ms, :), 2));
        joint_angles_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * squeeze(std(joint_angles_normalized_total(:, conditions_stanceR_stimNeg_0ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimNeg_0ms));
        joint_angles_mean_stanceR_stimPos_150ms = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceR_stimPos_150ms, :), 2));
        joint_angles_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * squeeze(std(joint_angles_normalized_total(:, conditions_stanceR_stimPos_150ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimPos_150ms));
        joint_angles_mean_stanceR_stimNeg_150ms = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceR_stimNeg_150ms, :), 2));
        joint_angles_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * squeeze(std(joint_angles_normalized_total(:, conditions_stanceR_stimNeg_150ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimNeg_150ms));
        joint_angles_mean_stanceR_stimPos_450ms = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceR_stimPos_450ms, :), 2));
        joint_angles_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * squeeze(std(joint_angles_normalized_total(:, conditions_stanceR_stimPos_450ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimPos_450ms));
        joint_angles_mean_stanceR_stimNeg_450ms = squeeze(mean(joint_angles_normalized_total(:, conditions_stanceR_stimNeg_450ms, :), 2));
        joint_angles_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * squeeze(std(joint_angles_normalized_total(:, conditions_stanceR_stimNeg_450ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimNeg_450ms));

        % responses
        joint_angles_response_mean_stanceL_stimPos_0ms = squeeze(mean(joint_angles_stanceL_response(:, conditions_stanceL_stimPos_0ms, :), 2));
        joint_angles_response_civ_stanceL_stimPos_0ms = tinv(0.975, sum(conditions_stanceL_stimPos_0ms)-1) * squeeze(std(joint_angles_stanceL_response(:, conditions_stanceL_stimPos_0ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimPos_0ms));
        joint_angles_response_mean_stanceL_stimNeg_0ms = squeeze(mean(joint_angles_stanceL_response(:, conditions_stanceL_stimNeg_0ms, :), 2));
        joint_angles_response_civ_stanceL_stimNeg_0ms = tinv(0.975, sum(conditions_stanceL_stimNeg_0ms)-1) * squeeze(std(joint_angles_stanceL_response(:, conditions_stanceL_stimNeg_0ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimNeg_0ms));
        joint_angles_response_mean_stanceL_stimPos_150ms = squeeze(mean(joint_angles_stanceL_response(:, conditions_stanceL_stimPos_150ms, :), 2));
        joint_angles_response_civ_stanceL_stimPos_150ms = tinv(0.975, sum(conditions_stanceL_stimPos_150ms)-1) * squeeze(std(joint_angles_stanceL_response(:, conditions_stanceL_stimPos_150ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimPos_150ms));
        joint_angles_response_mean_stanceL_stimNeg_150ms = squeeze(mean(joint_angles_stanceL_response(:, conditions_stanceL_stimNeg_150ms, :), 2));
        joint_angles_response_civ_stanceL_stimNeg_150ms = tinv(0.975, sum(conditions_stanceL_stimNeg_150ms)-1) * squeeze(std(joint_angles_stanceL_response(:, conditions_stanceL_stimNeg_150ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimNeg_150ms));
        joint_angles_response_mean_stanceL_stimPos_450ms = squeeze(mean(joint_angles_stanceL_response(:, conditions_stanceL_stimPos_450ms, :), 2));
        joint_angles_response_civ_stanceL_stimPos_450ms = tinv(0.975, sum(conditions_stanceL_stimPos_450ms)-1) * squeeze(std(joint_angles_stanceL_response(:, conditions_stanceL_stimPos_450ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimPos_450ms));
        joint_angles_response_mean_stanceL_stimNeg_450ms = squeeze(mean(joint_angles_stanceL_response(:, conditions_stanceL_stimNeg_450ms, :), 2));
        joint_angles_response_civ_stanceL_stimNeg_450ms = tinv(0.975, sum(conditions_stanceL_stimNeg_450ms)-1) * squeeze(std(joint_angles_stanceL_response(:, conditions_stanceL_stimNeg_450ms, :), 1, 2))/sqrt(sum(conditions_stanceL_stimNeg_450ms));

        joint_angles_response_mean_stanceR_stimPos_0ms = squeeze(mean(joint_angles_stanceR_response(:, conditions_stanceR_stimPos_0ms, :), 2));
        joint_angles_response_civ_stanceR_stimPos_0ms = tinv(0.975, sum(conditions_stanceR_stimPos_0ms)-1) * squeeze(std(joint_angles_stanceR_response(:, conditions_stanceR_stimPos_0ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimPos_0ms));
        joint_angles_response_mean_stanceR_stimNeg_0ms = squeeze(mean(joint_angles_stanceR_response(:, conditions_stanceR_stimNeg_0ms, :), 2));
        joint_angles_response_civ_stanceR_stimNeg_0ms = tinv(0.975, sum(conditions_stanceR_stimNeg_0ms)-1) * squeeze(std(joint_angles_stanceR_response(:, conditions_stanceR_stimNeg_0ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimNeg_0ms));
        joint_angles_response_mean_stanceR_stimPos_150ms = squeeze(mean(joint_angles_stanceR_response(:, conditions_stanceR_stimPos_150ms, :), 2));
        joint_angles_response_civ_stanceR_stimPos_150ms = tinv(0.975, sum(conditions_stanceR_stimPos_150ms)-1) * squeeze(std(joint_angles_stanceR_response(:, conditions_stanceR_stimPos_150ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimPos_150ms));
        joint_angles_response_mean_stanceR_stimNeg_150ms = squeeze(mean(joint_angles_stanceR_response(:, conditions_stanceR_stimNeg_150ms, :), 2));
        joint_angles_response_civ_stanceR_stimNeg_150ms = tinv(0.975, sum(conditions_stanceR_stimNeg_150ms)-1) * squeeze(std(joint_angles_stanceR_response(:, conditions_stanceR_stimNeg_150ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimNeg_150ms));
        joint_angles_response_mean_stanceR_stimPos_450ms = squeeze(mean(joint_angles_stanceR_response(:, conditions_stanceR_stimPos_450ms, :), 2));
        joint_angles_response_civ_stanceR_stimPos_450ms = tinv(0.975, sum(conditions_stanceR_stimPos_450ms)-1) * squeeze(std(joint_angles_stanceR_response(:, conditions_stanceR_stimPos_450ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimPos_450ms));
        joint_angles_response_mean_stanceR_stimNeg_450ms = squeeze(mean(joint_angles_stanceR_response(:, conditions_stanceR_stimNeg_450ms, :), 2));
        joint_angles_response_civ_stanceR_stimNeg_450ms = tinv(0.975, sum(conditions_stanceR_stimNeg_450ms)-1) * squeeze(std(joint_angles_stanceR_response(:, conditions_stanceR_stimNeg_450ms, :), 1, 2))/sqrt(sum(conditions_stanceR_stimNeg_450ms));

     
    end
    
    % stim start times
    stim_start_times_mean_stanceL_0ms = mean(stim_start_time_relative_to_stretch_total(conditions_stanceL_0ms));
    stim_start_times_mean_stanceR_0ms = mean(stim_start_time_relative_to_stretch_total(conditions_stanceR_0ms));
    stim_start_times_mean_stanceL_150ms = mean(stim_start_time_relative_to_stretch_total(conditions_stanceL_150ms));
    stim_start_times_mean_stanceR_150ms = mean(stim_start_time_relative_to_stretch_total(conditions_stanceR_150ms));
    stim_start_times_mean_stanceL_450ms = mean(stim_start_time_relative_to_stretch_total(conditions_stanceL_450ms));
    stim_start_times_mean_stanceR_450ms = mean(stim_start_time_relative_to_stretch_total(conditions_stanceR_450ms));
    
end

%% save data
if save_data

    save ...
      ( ...
        makeFileName(date, subject_id, 'resultsConditions'), ...
        'time_normalized', ...
        'origin_trial_list_total', ...
        'origin_start_time_list_total', ...
        'origin_end_time_list_total', ...
        'conditions_stanceL_stimNo', ...
        'conditions_stanceR_stimNo', ...
        'conditions_stanceL_stimPos_0ms', ...
        'conditions_stanceL_stimNeg_0ms', ...
        'conditions_stanceR_stimPos_0ms', ...
        'conditions_stanceR_stimNeg_0ms', ...
        'conditions_stanceL_stimPos_150ms', ...
        'conditions_stanceL_stimNeg_150ms', ...
        'conditions_stanceR_stimPos_150ms', ...
        'conditions_stanceR_stimNeg_150ms', ...
        'conditions_stanceL_stimPos_450ms', ...
        'conditions_stanceL_stimNeg_450ms', ...
        'conditions_stanceR_stimPos_450ms', ...
        'conditions_stanceR_stimNeg_450ms', ...
        'conditions_stanceL_0ms', ...
        'conditions_stanceR_0ms', ...
        'conditions_stanceL_150ms', ...
        'conditions_stanceR_150ms', ...
        'conditions_stanceL_450ms', ...
        'conditions_stanceR_450ms' ...
      );
    
    if process_data_marker
        save ...
          ( ...
            makeFileName(date, subject_id, 'resultsMarker'), ...
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
            'rheel_y_pos_normalized_total', ...
            'lheel_x_pos_stanceR_response', ...
            'rheel_x_pos_stanceL_response', ...
            'lheel_y_pos_stanceR_response', ...
            'rheel_y_pos_stanceL_response' ...
          );
    end
%             'rheel_x_pos_mean_stanceL_control', ...
%             'lheel_x_pos_mean_stanceR_control', ...
%             'rheel_x_pos_civ_stanceL_control', ...
%             'lheel_x_pos_civ_stanceR_control', ... 
%             'lheel_x_pos_mean_stanceR_stimPos_0ms', ...
%             'lheel_x_pos_mean_stanceR_stimNeg_0ms', ...
%             'lheel_x_pos_mean_stanceR_stimPos_150ms', ...
%             'lheel_x_pos_mean_stanceR_stimNeg_150ms', ...
%             'lheel_x_pos_mean_stanceR_stimPos_450ms', ...
%             'lheel_x_pos_mean_stanceR_stimNeg_450ms', ...
%             'lheel_x_pos_civ_stanceR_stimPos_0ms', ...
%             'lheel_x_pos_civ_stanceR_stimNeg_0ms', ...
%             'lheel_x_pos_civ_stanceR_stimPos_150ms', ...
%             'lheel_x_pos_civ_stanceR_stimNeg_150ms', ...
%             'lheel_x_pos_civ_stanceR_stimPos_450ms', ...
%             'lheel_x_pos_civ_stanceR_stimNeg_450ms', ...
%             'rheel_x_pos_mean_stanceL_stimPos_0ms', ...
%             'rheel_x_pos_mean_stanceL_stimNeg_0ms', ...
%             'rheel_x_pos_mean_stanceL_stimPos_150ms', ...
%             'rheel_x_pos_mean_stanceL_stimNeg_150ms', ...
%             'rheel_x_pos_mean_stanceL_stimPos_450ms', ...
%             'rheel_x_pos_mean_stanceL_stimNeg_450ms', ...
%             'rheel_x_pos_civ_stanceL_stimPos_0ms', ...
%             'rheel_x_pos_civ_stanceL_stimNeg_0ms', ...
%             'rheel_x_pos_civ_stanceL_stimPos_150ms', ...
%             'rheel_x_pos_civ_stanceL_stimNeg_150ms', ...
%             'rheel_x_pos_civ_stanceL_stimPos_450ms', ...
%             'rheel_x_pos_civ_stanceL_stimNeg_450ms', ...
%             'lheel_x_pos_response_mean_stanceR_stimPos_0ms', ...
%             'lheel_x_pos_response_mean_stanceR_stimNeg_0ms', ...
%             'lheel_x_pos_response_mean_stanceR_stimPos_150ms', ...
%             'lheel_x_pos_response_mean_stanceR_stimNeg_150ms', ...
%             'lheel_x_pos_response_mean_stanceR_stimPos_450ms', ...
%             'lheel_x_pos_response_mean_stanceR_stimNeg_450ms', ...
%             'lheel_x_pos_response_civ_stanceR_stimPos_0ms', ...
%             'lheel_x_pos_response_civ_stanceR_stimNeg_0ms', ...
%             'lheel_x_pos_response_civ_stanceR_stimPos_150ms', ...
%             'lheel_x_pos_response_civ_stanceR_stimNeg_150ms', ...
%             'lheel_x_pos_response_civ_stanceR_stimPos_450ms', ...
%             'lheel_x_pos_response_civ_stanceR_stimNeg_450ms', ...
%             'rheel_x_pos_response_mean_stanceL_stimPos_0ms', ...
%             'rheel_x_pos_response_mean_stanceL_stimNeg_0ms', ...
%             'rheel_x_pos_response_mean_stanceL_stimPos_150ms', ...
%             'rheel_x_pos_response_mean_stanceL_stimNeg_150ms', ...
%             'rheel_x_pos_response_mean_stanceL_stimPos_450ms', ...
%             'rheel_x_pos_response_mean_stanceL_stimNeg_450ms', ...
%             'rheel_x_pos_response_civ_stanceL_stimPos_0ms', ...
%             'rheel_x_pos_response_civ_stanceL_stimNeg_0ms', ...
%             'rheel_x_pos_response_civ_stanceL_stimPos_150ms', ...
%             'rheel_x_pos_response_civ_stanceL_stimNeg_150ms', ...
%             'rheel_x_pos_response_civ_stanceL_stimPos_450ms', ...
%             'rheel_x_pos_response_civ_stanceL_stimNeg_450ms' ...
            
    if process_data_forceplate
        save ...
          ( ...
            makeFileName(date, subject_id, 'resultsForceplate'), ...
            'fxl_normalized_total', 'fzl_normalized_total', 'myl_normalized_total', 'lcop_x_normalized_total', ...
            'fxr_normalized_total', 'fzr_normalized_total', 'myr_normalized_total', 'rcop_x_normalized_total', ...
            'lcop_x_stanceL_response', 'rcop_x_stanceR_response' ...
          );
%             'rcop_x_mean_stanceR_control', 'lcop_x_mean_stanceL_control', ...
%             'rcop_x_civ_stanceR_control', 'lcop_x_civ_stanceL_control', ...
%             'rcop_x_mean_stanceR_stimPos_0ms', 'rcop_x_mean_stanceR_stimNeg_0ms', 'rcop_x_mean_stanceR_stimPos_150ms', 'rcop_x_mean_stanceR_stimNeg_150ms', 'lcop_x_mean_stanceL_stimPos_450ms', 'lcop_x_mean_stanceL_stimNeg_450ms', ...
%             'rcop_x_civ_stanceR_stimPos_0ms', 'rcop_x_civ_stanceR_stimNeg_0ms', 'rcop_x_civ_stanceR_stimPos_150ms', 'rcop_x_civ_stanceR_stimNeg_150ms', 'lcop_x_civ_stanceL_stimPos_450ms', 'lcop_x_civ_stanceL_stimNeg_450ms', ...
%             'rcop_x_response_mean_stanceR_stimPos_0ms', 'rcop_x_response_mean_stanceR_stimNeg_0ms', 'rcop_x_response_mean_stanceR_stimPos_150ms', 'rcop_x_response_mean_stanceR_stimNeg_150ms', 'lcop_x_response_mean_stanceL_stimPos_450ms', 'lcop_x_response_mean_stanceL_stimNeg_450ms', ...
%             'rcop_x_response_civ_stanceR_stimPos_0ms', 'rcop_x_response_civ_stanceR_stimNeg_0ms', 'rcop_x_response_civ_stanceR_stimPos_150ms', 'rcop_x_response_civ_stanceR_stimNeg_150ms', 'lcop_x_response_civ_stanceL_stimPos_450ms', 'lcop_x_response_civ_stanceL_stimNeg_450ms' ...
    end
end
    

return    


%% angle-related stuff
if do_joint_angle_absolute_plots_left
    joints_to_plot = [8 12 14 18];
%     joints_to_plot = [20 23];
    for i_joint = joints_to_plot
        figure; axes; hold on; title([plant.jointLabels{i_joint} ', LEFT 0ms']); set(gca, 'Fontsize', 12)
        plot(time_normalized, squeeze(joint_angles_normalized_total(:, conditions_stanceL_stimNo, i_joint)), 'color', color_left_control);
        plot(time_normalized, squeeze(joint_angles_normalized_total(:, conditions_stanceL_stimPos_0ms, i_joint)), 'color', color_left_positive);
        plot(time_normalized, squeeze(joint_angles_normalized_total(:, conditions_stanceL_stimNeg_0ms, i_joint)), 'color', color_left_negative);
        plot(time_normalized, mean(squeeze(joint_angles_normalized_total(:, conditions_stanceL_stimNo, i_joint)), 2), 'color', color_left_control, 'linewidth', 5);
        plot(time_normalized, mean(squeeze(joint_angles_normalized_total(:, conditions_stanceL_stimPos_0ms, i_joint)), 2), 'color', color_left_positive, 'linewidth', 5);
        plot(time_normalized, mean(squeeze(joint_angles_normalized_total(:, conditions_stanceL_stimNeg_0ms, i_joint)), 2), 'color', color_left_negative, 'linewidth', 5);
    end    
end

if do_joint_angle_absolute_plots_right
    joints_to_plot = [8 12 14 18];
%     joints_to_plot = [20 23];
    for i_joint = joints_to_plot
        figure; axes; hold on; title([plant.jointLabels{i_joint} ', RIGHT 0ms']); set(gca, 'Fontsize', 12)
        plot(time_normalized, squeeze(joint_angles_normalized_total(:, conditions_stanceR_stimNo, i_joint)), 'color', color_right_control);
        plot(time_normalized, squeeze(joint_angles_normalized_total(:, conditions_stanceR_stimPos_0ms, i_joint)), 'color', color_right_positive);
        plot(time_normalized, squeeze(joint_angles_normalized_total(:, conditions_stanceR_stimNeg_0ms, i_joint)), 'color', color_right_negative);
        plot(time_normalized, mean(squeeze(joint_angles_normalized_total(:, conditions_stanceR_stimNo, i_joint)), 2), 'color', color_right_control, 'linewidth', 5);
        plot(time_normalized, mean(squeeze(joint_angles_normalized_total(:, conditions_stanceR_stimPos_0ms, i_joint)), 2), 'color', color_right_positive, 'linewidth', 5);
        plot(time_normalized, mean(squeeze(joint_angles_normalized_total(:, conditions_stanceR_stimNeg_0ms, i_joint)), 2), 'color', color_right_negative, 'linewidth', 5);
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














