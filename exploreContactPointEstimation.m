
extract_data                            = 0;
estimate_contact_point_trajectories     = 0;
fit_rollover_shape                      = 1;


trial_to_process = 1;



if extract_data
    % load data
    load subjectInfo.mat;
    load(makeFileName(date, subject_id, 'walking', trial_to_process, 'markerTrajectories'));
    load(makeFileName(date, subject_id, 'walking', trial_to_process, 'labviewTrajectories'));
    load(makeFileName(date, subject_id, 'walking', trial_to_process, 'forceplateTrajectories'));
    load(makeFileName(date, subject_id, 'walking', trial_to_process, 'stepEvents'));

    time_steps_to_consider = 1 : 1000;
    time_steps_to_consider = 1 : length(time_mocap);
    number_of_time_steps = length(time_steps_to_consider);
    
    % calculate belt coordinates
    belt_speed_trajectory = mean([belt_speed_left_trajectory belt_speed_right_trajectory], 2);
    delta_t = diff(time_labview);
    belt_position_trajectory_labview = zeros(size(belt_speed_trajectory));
    for i_time = 2 : length(belt_speed_trajectory)
        belt_position_trajectory_labview(i_time) = belt_position_trajectory_labview(i_time-1) + delta_t(i_time-1) * belt_speed_trajectory(i_time-1);
    end
    belt_position_trajectory_mocap = spline(time_labview, belt_position_trajectory_labview, time_mocap);

    % extract data
    left_ankle_marker = find(strcmp(marker_headers, 'LANK'));
    left_ankle_marker_indices = reshape([(left_ankle_marker - 1) * 3 + 1; (left_ankle_marker - 1) * 3 + 2; (left_ankle_marker - 1) * 3 + 3], 1, length(left_ankle_marker)*3);
    left_ankle_marker_yz_pos_trajectory = marker_trajectories(time_steps_to_consider, left_ankle_marker_indices(2:3));
    left_ankle_marker_yz_pos_trajectory(:, 1) = left_ankle_marker_yz_pos_trajectory(:, 1) + belt_position_trajectory_mocap(time_steps_to_consider)';

    left_toes_marker = find(strcmp(marker_headers, 'LTOE'));
    left_toes_marker_indices = reshape([(left_toes_marker - 1) * 3 + 1; (left_toes_marker - 1) * 3 + 2; (left_toes_marker - 1) * 3 + 3], 1, length(left_toes_marker)*3);
    left_toes_marker_yz_pos_trajectory = marker_trajectories(time_steps_to_consider, left_toes_marker_indices(2:3));
    left_toes_marker_yz_pos_trajectory(:, 1) = left_toes_marker_yz_pos_trajectory(:, 1) + belt_position_trajectory_mocap(time_steps_to_consider)';

    % calculate derivatives
    filter_order = 4;
    cutoff_frequency = 20; % cutoff frequency, in Hz
    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate_mocap/2));	% set filter parameters for butterworth filter: 2=order of filter;
    left_ankle_marker_yz_vel_trajectory = deriveByTime(filtfilt(b, a, left_ankle_marker_yz_pos_trajectory), 1/sampling_rate_mocap);
    left_ankle_marker_yz_acc_trajectory = deriveByTime(filtfilt(b, a, left_ankle_marker_yz_vel_trajectory), 1/sampling_rate_mocap);
    left_toes_marker_yz_vel_trajectory = deriveByTime(filtfilt(b, a, left_toes_marker_yz_pos_trajectory), 1/sampling_rate_mocap);
    left_toes_marker_yz_acc_trajectory = deriveByTime(filtfilt(b, a, left_toes_marker_yz_vel_trajectory), 1/sampling_rate_mocap);
    
    % set irrelevant data points to NaN
    irrelevant_data_points = ~left_contact_indicators_mocap(time_steps_to_consider);
    left_ankle_marker_yz_pos_trajectory(irrelevant_data_points, :) = NaN;
    left_ankle_marker_yz_vel_trajectory(irrelevant_data_points, :) = NaN;
    left_ankle_marker_yz_acc_trajectory(irrelevant_data_points, :) = NaN;
    left_toes_marker_yz_pos_trajectory(irrelevant_data_points, :) = NaN;
    left_toes_marker_yz_vel_trajectory(irrelevant_data_points, :) = NaN;
    left_toes_marker_yz_acc_trajectory(irrelevant_data_points, :) = NaN;
    
    % spline smooth
%     figure; axes; hold on; axis equal; title('ankle trajectory, spline smoothed, world frame')
    i_time = 1;
    while i_time < number_of_time_steps
        % find next stretch start
        stretch_start = i_time-1 + find(~irrelevant_data_points(i_time:end), 1, 'first');
        stretch_end = stretch_start-1 + find(irrelevant_data_points(stretch_start:end), 1, 'first') - 1;
        
        if ~isempty(stretch_start) && ~isempty(stretch_end)
            % extract and fit data
            curvefit = fit(left_ankle_marker_yz_pos_trajectory(stretch_start:stretch_end, 1), left_ankle_marker_yz_pos_trajectory(stretch_start:stretch_end, 2), 'smoothingspline', 'SmoothingParam', 0.9999999);
%             plot(fitobject,left_ankle_marker_yz_pos_trajectory(stretch_start:stretch_end, 1), left_ankle_marker_yz_pos_trajectory(stretch_start:stretch_end, 2)); 
%             hold on
            
            % evaluate and differentiate fit
            fit_y_data = left_ankle_marker_yz_pos_trajectory(stretch_start:stretch_end, 1);
            fit_z_data = curvefit(fit_y_data);
            fit_dz_by_dy = differentiate(curvefit, fit_y_data);
            dy = 0.1;
            fit_dy_data = ones(size(fit_dz_by_dy)) * dy;
            fit_dz_data = fit_dz_by_dy * dy;
            
            % save fit in place of original data
            left_ankle_marker_yz_pos_trajectory(stretch_start:stretch_end, 1) = fit_y_data;
            left_ankle_marker_yz_pos_trajectory(stretch_start:stretch_end, 2) = fit_z_data;
            left_ankle_marker_yz_vel_trajectory(stretch_start:stretch_end, 1) = fit_dy_data;
            left_ankle_marker_yz_vel_trajectory(stretch_start:stretch_end, 2) = fit_dz_data;
            
            
%             example_index = round(length(stretch_start:stretch_end)/2) + 100;
%             example_y = fit_y_data(example_index);
%             example_z = fit_z_data(example_index);
%             example_dz_by_dy = fit_dz_by_dy(example_index);
            
%             plot(fit_y_data, fit_z_data, '-');
%             hold on
%             plot(example_y, example_z, 'bo');
%             plot(example_y + [-1 1]*dy, example_z + [-1 1]*dy*example_dz_by_dy, 'b-');
            
%             set(ankle_pos_plot, 'xdata', left_ankle_point_position(1), 'ydata', left_ankle_point_position(2));
%             set(ankle_vel_plot, 'xdata', left_ankle_point_position(1) + left_ankle_point_velocity_normed(1)*[-1 1]*vel_scaler, 'ydata', left_ankle_point_position(2) + left_ankle_point_velocity_normed(2)*[-1 1]*vel_scaler);
            
            
            
            i_time = stretch_end + 1;
        else
            i_time = number_of_time_steps;
        end
    end

    
    
end



%% estimate contact point
if estimate_contact_point_trajectories
    left_contact_point_world_trajectory = zeros(number_of_time_steps, 2) * NaN;
    left_contact_point_foot_trajectory = zeros(number_of_time_steps, 2) * NaN;
    left_toes_position_foot_trajectory = zeros(number_of_time_steps, 2) * NaN;
    left_theta_trajectory = zeros(number_of_time_steps, 1) * NaN;
    
    
    % prep visualization
    figure; axes; hold on; axis equal; title(['ankle trajectory, world frame - ' subject_id])
    set(gca, 'xlim', [1.55 1.95], 'ylim', [-0.01 0.3]); 
    ankle_trajectory_plot = plot(left_ankle_marker_yz_pos_trajectory(:, 1), left_ankle_marker_yz_pos_trajectory(:, 2));
    ankle_pos_plot = plot(0, 0, 'x', 'markersize', 12);
    ankle_vel_plot = plot(0, 0, '-', 'markersize', 12);
    ankle_vel_normal_plot = plot(0, 0, '-', 'markersize', 12);
    contact_plot = plot(0, 0, 'o', 'markersize', 12);
    vel_scaler = 10;
    
    figure; axes; hold on; axis equal; title(['contact trajectory, foot frame - ' subject_id])
    contact_trajectory_foot_plot = plot(0, 0, '-', 'linewidth', 2);
    
    
    for i_time = 1 : number_of_time_steps
        left_ankle_point_position_world = left_ankle_marker_yz_pos_trajectory(i_time, :)';
        left_toes_position_world = left_toes_marker_yz_pos_trajectory(i_time, :)';
        left_ankle_point_velocity_normed = normVector(left_ankle_marker_yz_vel_trajectory(i_time, :)');
        left_ankle_point_velocity_normal = [left_ankle_point_velocity_normed(2); -left_ankle_point_velocity_normed(1)];
        contact_point_world = lineIntersection(left_ankle_point_position_world, [0; 0], left_ankle_point_velocity_normal, [1; 0]);
        left_contact_point_world_trajectory(i_time, :) = contact_point_world;
        
        left_ankle_to_contact_vector = contact_point_world - left_ankle_point_position_world;
        left_ankle_to_toes_vector = left_toes_position_world - left_ankle_point_position_world;
        left_theta_trajectory(i_time, :) = -atan2(left_ankle_to_toes_vector(2), left_ankle_to_toes_vector(1));
        R_foot_to_world = [cos(left_theta_trajectory(i_time, :)) sin(left_theta_trajectory(i_time, :)); -sin(left_theta_trajectory(i_time, :)), cos(left_theta_trajectory(i_time, :))];
        T_foot_to_world = [R_foot_to_world left_ankle_point_position_world; 0 0 1];
        
        if ~isnan(left_ankle_marker_yz_pos_trajectory(i_time, 1))
            left_contact_point_foot_trajectory(i_time, :) = eye(2, 3) * T_foot_to_world^(-1) * [contact_point_world; 1];
            left_toes_position_foot_trajectory(i_time, :) = eye(2, 3) * T_foot_to_world^(-1) * [left_toes_position_world; 1];
        end
        set(ankle_pos_plot, 'xdata', left_ankle_point_position_world(1), 'ydata', left_ankle_point_position_world(2));
        set(ankle_vel_plot, 'xdata', left_ankle_point_position_world(1) + left_ankle_point_velocity_normed(1)*[-1 1]*vel_scaler, 'ydata', left_ankle_point_position_world(2) + left_ankle_point_velocity_normed(2)*[-1 1]*vel_scaler);
        set(ankle_vel_normal_plot, 'xdata', left_ankle_point_position_world(1) + left_ankle_point_velocity_normal(1)*[-1 1]*vel_scaler, 'ydata', left_ankle_point_position_world(2) + left_ankle_point_velocity_normal(2)*[-1 1]*vel_scaler);
        set(contact_plot, 'xdata', contact_point_world(1), 'ydata', contact_point_world(2));
        
        set(contact_trajectory_foot_plot, 'xdata', left_contact_point_foot_trajectory(:, 1), 'ydata', left_contact_point_foot_trajectory(:, 2));
        
%         drawnow;
    end
end
   
%% fit contact surface with smoothing splines
if fit_rollover_shape
    left_toes_position_foot = nanmean(left_toes_position_foot_trajectory);
    rollover_shape = left_contact_point_foot_trajectory;
    nan_indices = isnan(left_contact_point_foot_trajectory(:, 1));
    rollover_shape(nan_indices, :) = [];
    
    rollover_shape_x = rollover_shape(:, 1);
    rollover_shape_y = rollover_shape(:, 2);
    
%     curvefit = fit(rollover_shape(:, 1), rollover_shape(:, 2), 'smoothingspline', 'SmoothingParam', 0.999999);
    fit_y_data = left_contact_point_foot_trajectory(:, 1);
    fit_z_data = curvefit(fit_y_data);
    
    figure; axes; hold on; axis equal
    plot(left_contact_point_foot_trajectory(:, 1), left_contact_point_foot_trajectory(:, 2));
    plot(fit_y_data, fit_z_data, '-', 'linewidth', 3);
    
    plot(0, 0, 'o', 'markersize', 10, 'color', 'k', 'markerfacecolor', 'k')
    plot(left_toes_position_foot(1), left_toes_position_foot(2), 'o', 'markersize', 10, 'color', 'b', 'markerfacecolor', 'b')
    
    legend('contact point trajectory', 'contact point spline smoothed', 'ankle', 'toes')
end


















