
% load data
trial_number = 2;
load subjectInfo.mat;
load(makeFileName(date, subject_id, 'model'));
load(makeFileName(date, subject_id, 'walking', trial_number, 'markerTrajectories'));
load(makeFileName(date, subject_id, 'walking', trial_number, 'angleTrajectories'));
load(makeFileName(date, subject_id, 'walking', trial_number, 'kinematicTrajectories'));
load(makeFileName(date, subject_id, 'walking', trial_number, 'stepEvents'));



number_of_time_steps = size(T_left_ankle_to_world_trajectory, 1);

calculate_trajectories              = 1;
find_phi_constraint                 = 1;
find_rho_constraint_by_phi          = 0;
save_results                        = 1;

number_of_points_of_interest = 200;
number_of_phi_values = 200;
number_of_joints = plant.numberOfJoints;

plant.jointAngles = zeros(number_of_joints, 1);
plant.updateKinematics;
left_ankle_scs_to_world_rotation_reference = plant.endEffectorTransformations{3}(1:3, 1:3);
right_ankle_scs_to_world_rotation_reference = plant.endEffectorTransformations{6}(1:3, 1:3);

    
%% calculate trajectories
if calculate_trajectories

    % calculate z-trajectories for multiple points along the sphere z-axis
    phi_left_trajectory = zeros(number_of_time_steps, 1);
    rho_left_trajectory = zeros(number_of_time_steps, 1);
    gamma_left_trajectory = zeros(number_of_time_steps, 1);
    phi_right_trajectory = zeros(number_of_time_steps, 1);
    rho_right_trajectory = zeros(number_of_time_steps, 1);
    gamma_right_trajectory = zeros(number_of_time_steps, 1);

    points_of_interest_y_local = linspace(-0.02, 0.08, number_of_points_of_interest);
    points_of_interest_sphere = [zeros(1, number_of_points_of_interest); points_of_interest_y_local; zeros(1, number_of_points_of_interest); ones(1, number_of_points_of_interest)]; % coordinates of the points of interest in sphere coordinates
    left_points_of_interest_world_trajectories_x = zeros(number_of_time_steps, size(points_of_interest_sphere, 2));
    left_points_of_interest_world_trajectories_y = zeros(number_of_time_steps, size(points_of_interest_sphere, 2));
    left_points_of_interest_world_trajectories_z = zeros(number_of_time_steps, size(points_of_interest_sphere, 2));
    right_points_of_interest_world_trajectories_x = zeros(number_of_time_steps, size(points_of_interest_sphere, 2));
    right_points_of_interest_world_trajectories_y = zeros(number_of_time_steps, size(points_of_interest_sphere, 2));
    right_points_of_interest_world_trajectories_z = zeros(number_of_time_steps, size(points_of_interest_sphere, 2));
    for i_time = 1 : number_of_time_steps
        % left ankle euler angles
        left_ankle_scs_transformation_current = reshape(T_left_ankle_to_world_trajectory(i_time, :), 4, 4);
        left_ankle_scs_to_world_rotation_current = left_ankle_scs_transformation_current(1:3, 1:3);
        left_ankle_scs_rotation_reference_to_current = left_ankle_scs_to_world_rotation_reference^(-1) * left_ankle_scs_to_world_rotation_current;
        left_euler_angles = eulerAnglesFromRotationMatrixYZX(left_ankle_scs_rotation_reference_to_current);
        gamma_left_trajectory(i_time) = left_euler_angles(1);
        phi_left_trajectory(i_time) = left_euler_angles(2);
        rho_left_trajectory(i_time) = left_euler_angles(3);

        % left ankle points of interest
        left_points_of_interest_world_current = left_ankle_scs_transformation_current * points_of_interest_sphere;
        left_points_of_interest_world_trajectories_x(i_time, :) = left_points_of_interest_world_current(1, :);
        left_points_of_interest_world_trajectories_y(i_time, :) = left_points_of_interest_world_current(2, :);
        left_points_of_interest_world_trajectories_z(i_time, :) = left_points_of_interest_world_current(3, :);

        % right ankle euler angles
        right_ankle_scs_transformation_current = reshape(T_right_ankle_to_world_trajectory(i_time, :), 4, 4);
        right_ankle_scs_to_world_rotation_current = right_ankle_scs_transformation_current(1:3, 1:3);
        right_ankle_scs_rotation_reference_to_current = right_ankle_scs_to_world_rotation_reference^(-1) * right_ankle_scs_to_world_rotation_current;
        right_euler_angles = eulerAnglesFromRotationMatrixYZX(right_ankle_scs_rotation_reference_to_current);
        gamma_right_trajectory(i_time) = right_euler_angles(1);
        phi_right_trajectory(i_time) = right_euler_angles(2);
        rho_right_trajectory(i_time) = right_euler_angles(3);

        % right ankle points of interest
        right_points_of_interest_world_current = right_ankle_scs_transformation_current * points_of_interest_sphere;
        right_points_of_interest_world_trajectories_x(i_time, :) = right_points_of_interest_world_current(1, :);
        right_points_of_interest_world_trajectories_y(i_time, :) = right_points_of_interest_world_current(2, :);
        right_points_of_interest_world_trajectories_z(i_time, :) = right_points_of_interest_world_current(3, :);
    end

    % get relevant data points
    relevant_data_points_left = left_contact_indicators_mocap;
    irrelevant_data_points_left = ~relevant_data_points_left;
    relevant_data_points_right = right_contact_indicators_mocap;
    irrelevant_data_points_right = ~relevant_data_points_right;

    phi_left_relevant_trajectory = phi_left_trajectory;   phi_left_relevant_trajectory(irrelevant_data_points_left) = NaN;
    rho_left_relevant_trajectory = rho_left_trajectory;   rho_left_relevant_trajectory(irrelevant_data_points_left) = NaN;
    phi_right_relevant_trajectory = phi_right_trajectory;   phi_right_relevant_trajectory(irrelevant_data_points_right) = NaN;
    rho_right_relevant_trajectory = rho_right_trajectory;   rho_right_relevant_trajectory(irrelevant_data_points_right) = NaN;

    phi_dot_left_trajectory = deriveByTime(phi_left_trajectory, 1/sampling_rate_mocap);
    rho_dot_left_trajectory = deriveByTime(rho_left_trajectory, 1/sampling_rate_mocap);
    gamma_dot_left_trajectory = deriveByTime(gamma_left_trajectory, 1/sampling_rate_mocap);
    phi_dot_left_relevant_trajectory = deriveByTime(phi_left_relevant_trajectory, 1/sampling_rate_mocap);
    rho_dot_left_relevant_trajectory = deriveByTime(rho_left_relevant_trajectory, 1/sampling_rate_mocap);
    phi_dot_right_trajectory = deriveByTime(phi_right_trajectory, 1/sampling_rate_mocap);
    rho_dot_right_trajectory = deriveByTime(rho_right_trajectory, 1/sampling_rate_mocap);
    gamma_dot_right_trajectory = deriveByTime(gamma_right_trajectory, 1/sampling_rate_mocap);
    phi_dot_right_relevant_trajectory = deriveByTime(phi_right_relevant_trajectory, 1/sampling_rate_mocap);
    rho_dot_right_relevant_trajectory = deriveByTime(rho_right_relevant_trajectory, 1/sampling_rate_mocap);

    body_velocity_left_relevant_trajectory = V_body_left_ankle;   body_velocity_left_relevant_trajectory(irrelevant_data_points_left, :) = NaN;
    V_body_left_ankle_relevant_trajectory = V_body_left_ankle; V_body_left_ankle_relevant_trajectory(irrelevant_data_points_left, :) = NaN;
    body_velocity_right_relevant_trajectory = V_body_right_ankle;   body_velocity_right_relevant_trajectory(irrelevant_data_points_right, :) = NaN;
    V_body_right_ankle_relevant_trajectory = V_body_right_ankle; V_body_right_ankle_relevant_trajectory(irrelevant_data_points_right, :) = NaN;

    body_velocity_left_relevant_trajectory_norms = sum(body_velocity_left_relevant_trajectory.^2, 2).^(0.5);
    body_velocity_left_relevant_trajectory_normed = body_velocity_left_relevant_trajectory ./ repmat(body_velocity_left_relevant_trajectory_norms, 1, 6);

%     % plot phi and rho
%     figure; hold on;
%     plot(time_mocap, phi_left_trajectory);
%     plot(time_mocap, rho_left_trajectory);
%     plot(time_mocap, phi_right_trajectory);
%     plot(time_mocap, rho_right_trajectory);
%     legend('\phi left', '\rho left', '\phi right', '\rho right')

    % plot left relevant
    figure; hold on;
    plot(time_mocap, phi_left_trajectory);
    plot(time_mocap, rho_left_trajectory);
    plot(time_mocap, phi_left_relevant_trajectory, 'linewidth', 3);
    plot(time_mocap, rho_left_relevant_trajectory, 'linewidth', 3);
    legend('\phi left', '\rho left', '\phi left relevant', '\rho left relevant')

    % plot right relevant
    figure; hold on;
    plot(time_mocap, phi_right_trajectory);
    plot(time_mocap, rho_right_trajectory);
    plot(time_mocap, phi_right_relevant_trajectory, 'linewidth', 3);
    plot(time_mocap, rho_right_relevant_trajectory, 'linewidth', 3);
    legend('\phi right', '\rho right', '\phi right relevant', '\rho right relevant')

%     % plot body velocity
%     figure; V_body_left_ankle_axes = axes; hold on; title('left ankle body velocity')
%     plot(V_body_left_ankle(:, 1:6));
%     legend('V_1', 'V_2', 'V_3', 'V_4', 'V_5', 'V_6');
% 
%     figure; V_body_right_ankle_axes = axes; hold on; title('right ankle body velocity')
%     plot(V_body_right_ankle(:, 1:6));
%     legend('V_1', 'V_2', 'V_3', 'V_4', 'V_5', 'V_6');

%     linkaxes([V_body_left_ankle_axes V_body_right_ankle_axes], 'x')

end
    
    
%% find_phi_constraint
if find_phi_constraint
    plot_relevant_body_velocity_data_left    = 1;
    plot_relevant_body_velocity_data_right   = 1;
    
    % find data points for each constraint
    phi_threshold = 0.05;
    left_heelstrike_indices = phi_left_relevant_trajectory > phi_threshold;
    left_pushoff_indices = phi_left_relevant_trajectory < -phi_threshold;
    right_heelstrike_indices = phi_right_relevant_trajectory > phi_threshold;
    right_pushoff_indices = phi_right_relevant_trajectory < -phi_threshold;
    
    body_velocity_left_heelstrike_trajectory = V_body_left_ankle;   body_velocity_left_heelstrike_trajectory(~left_heelstrike_indices, :) = NaN;
    body_velocity_left_pushoff_trajectory = V_body_left_ankle;   body_velocity_left_pushoff_trajectory(~left_pushoff_indices, :) = NaN;
    body_velocity_right_heelstrike_trajectory = V_body_right_ankle;   body_velocity_right_heelstrike_trajectory(~right_heelstrike_indices, :) = NaN;
    body_velocity_right_pushoff_trajectory = V_body_right_ankle;   body_velocity_right_pushoff_trajectory(~right_pushoff_indices, :) = NaN;
    
    % fit constraints
    [V_body_left_1_fit_heelstrike, V_body_left_1_gof_heelstrike] = fit([phi_left_trajectory(left_heelstrike_indices), rho_left_trajectory(left_heelstrike_indices)], V_body_left_ankle(left_heelstrike_indices, 1), 'poly55');
    [V_body_left_2_fit_heelstrike, V_body_left_2_gof_heelstrike] = fit([phi_left_trajectory(left_heelstrike_indices), rho_left_trajectory(left_heelstrike_indices)], V_body_left_ankle(left_heelstrike_indices, 2), 'poly55');
    [V_body_left_3_fit_heelstrike, V_body_left_3_gof_heelstrike] = fit([phi_left_trajectory(left_heelstrike_indices), rho_left_trajectory(left_heelstrike_indices)], V_body_left_ankle(left_heelstrike_indices, 3), 'poly55');
    [V_body_left_4_fit_heelstrike, V_body_left_4_gof_heelstrike] = fit([phi_left_trajectory(left_heelstrike_indices), rho_left_trajectory(left_heelstrike_indices)], V_body_left_ankle(left_heelstrike_indices, 4), 'poly55');
    [V_body_left_5_fit_heelstrike, V_body_left_5_gof_heelstrike] = fit([phi_left_trajectory(left_heelstrike_indices), rho_left_trajectory(left_heelstrike_indices)], V_body_left_ankle(left_heelstrike_indices, 5), 'poly55');
    [V_body_left_6_fit_heelstrike, V_body_left_6_gof_heelstrike] = fit([phi_left_trajectory(left_heelstrike_indices), rho_left_trajectory(left_heelstrike_indices)], V_body_left_ankle(left_heelstrike_indices, 6), 'poly55');
    [V_body_left_1_fit_pushoff, V_body_left_1_gof_pushoff] = fit([phi_left_trajectory(left_pushoff_indices), rho_left_trajectory(left_pushoff_indices)], V_body_left_ankle(left_pushoff_indices, 1), 'poly55');
    [V_body_left_2_fit_pushoff, V_body_left_2_gof_pushoff] = fit([phi_left_trajectory(left_pushoff_indices), rho_left_trajectory(left_pushoff_indices)], V_body_left_ankle(left_pushoff_indices, 2), 'poly55');
    [V_body_left_3_fit_pushoff, V_body_left_3_gof_pushoff] = fit([phi_left_trajectory(left_pushoff_indices), rho_left_trajectory(left_pushoff_indices)], V_body_left_ankle(left_pushoff_indices, 3), 'poly55');
    [V_body_left_4_fit_pushoff, V_body_left_4_gof_pushoff] = fit([phi_left_trajectory(left_pushoff_indices), rho_left_trajectory(left_pushoff_indices)], V_body_left_ankle(left_pushoff_indices, 4), 'poly55');
    [V_body_left_5_fit_pushoff, V_body_left_5_gof_pushoff] = fit([phi_left_trajectory(left_pushoff_indices), rho_left_trajectory(left_pushoff_indices)], V_body_left_ankle(left_pushoff_indices, 5), 'poly55');
    [V_body_left_6_fit_pushoff, V_body_left_6_gof_pushoff] = fit([phi_left_trajectory(left_pushoff_indices), rho_left_trajectory(left_pushoff_indices)], V_body_left_ankle(left_pushoff_indices, 6), 'poly55');
    
    [V_body_right_1_fit_heelstrike, V_body_right_1_gof_heelstrike] = fit([phi_right_trajectory(right_heelstrike_indices), rho_right_trajectory(right_heelstrike_indices)], V_body_right_ankle(right_heelstrike_indices, 1), 'poly55');
    [V_body_right_2_fit_heelstrike, V_body_right_2_gof_heelstrike] = fit([phi_right_trajectory(right_heelstrike_indices), rho_right_trajectory(right_heelstrike_indices)], V_body_right_ankle(right_heelstrike_indices, 2), 'poly55');
    [V_body_right_3_fit_heelstrike, V_body_right_3_gof_heelstrike] = fit([phi_right_trajectory(right_heelstrike_indices), rho_right_trajectory(right_heelstrike_indices)], V_body_right_ankle(right_heelstrike_indices, 3), 'poly55');
    [V_body_right_4_fit_heelstrike, V_body_right_4_gof_heelstrike] = fit([phi_right_trajectory(right_heelstrike_indices), rho_right_trajectory(right_heelstrike_indices)], V_body_right_ankle(right_heelstrike_indices, 4), 'poly55');
    [V_body_right_5_fit_heelstrike, V_body_right_5_gof_heelstrike] = fit([phi_right_trajectory(right_heelstrike_indices), rho_right_trajectory(right_heelstrike_indices)], V_body_right_ankle(right_heelstrike_indices, 5), 'poly55');
    [V_body_right_6_fit_heelstrike, V_body_right_6_gof_heelstrike] = fit([phi_right_trajectory(right_heelstrike_indices), rho_right_trajectory(right_heelstrike_indices)], V_body_right_ankle(right_heelstrike_indices, 6), 'poly55');
    [V_body_right_1_fit_pushoff, V_body_right_1_gof_pushoff] = fit([phi_right_trajectory(right_pushoff_indices), rho_right_trajectory(right_pushoff_indices)], V_body_right_ankle(right_pushoff_indices, 1), 'poly55');
    [V_body_right_2_fit_pushoff, V_body_right_2_gof_pushoff] = fit([phi_right_trajectory(right_pushoff_indices), rho_right_trajectory(right_pushoff_indices)], V_body_right_ankle(right_pushoff_indices, 2), 'poly55');
    [V_body_right_3_fit_pushoff, V_body_right_3_gof_pushoff] = fit([phi_right_trajectory(right_pushoff_indices), rho_right_trajectory(right_pushoff_indices)], V_body_right_ankle(right_pushoff_indices, 3), 'poly55');
    [V_body_right_4_fit_pushoff, V_body_right_4_gof_pushoff] = fit([phi_right_trajectory(right_pushoff_indices), rho_right_trajectory(right_pushoff_indices)], V_body_right_ankle(right_pushoff_indices, 4), 'poly55');
    [V_body_right_5_fit_pushoff, V_body_right_5_gof_pushoff] = fit([phi_right_trajectory(right_pushoff_indices), rho_right_trajectory(right_pushoff_indices)], V_body_right_ankle(right_pushoff_indices, 5), 'poly55');
    [V_body_right_6_fit_pushoff, V_body_right_6_gof_pushoff] = fit([phi_right_trajectory(right_pushoff_indices), rho_right_trajectory(right_pushoff_indices)], V_body_right_ankle(right_pushoff_indices, 6), 'poly55');
    
    if plot_relevant_body_velocity_data_left
        % visualize - normalized body velocity vs. phi and rho
        figure; axes; hold on; title('v_1');
        plot3(phi_left_trajectory, rho_left_trajectory, V_body_left_ankle(:, 1));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_relevant_trajectory(:, 1));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_heelstrike_trajectory(:, 1));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_pushoff_trajectory(:, 1));
        plot3(phi_left_trajectory(left_heelstrike_indices), rho_left_trajectory(left_heelstrike_indices), feval(V_body_left_1_fit_heelstrike, [phi_left_trajectory(left_heelstrike_indices) rho_left_trajectory(left_heelstrike_indices)]), 'o');
        plot3(phi_left_trajectory(left_pushoff_indices), rho_left_trajectory(left_pushoff_indices), feval(V_body_left_1_fit_pushoff, [phi_left_trajectory(left_pushoff_indices) rho_left_trajectory(left_pushoff_indices)]), 'o');
        xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
%         legend('all', 'contact', 'heelstrike', 'pushoff', 'heelstrike fit', 'pushoff fit')
        uicontrol('style', 'text', 'position', [10 25 600 30], 'string', ['HS: R-square = ' num2str(V_body_left_1_gof_heelstrike.rsquare) ', RMS error = ' num2str(V_body_left_1_gof_heelstrike.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')
        uicontrol('style', 'text', 'position', [10 5 600 30], 'string', ['PO: R-square = ' num2str(V_body_left_1_gof_pushoff.rsquare) ', RMS error = ' num2str(V_body_left_1_gof_pushoff.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')

        figure; axes; hold on; title('v_2');
        plot3(phi_left_trajectory, rho_left_trajectory, V_body_left_ankle(:, 2));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_relevant_trajectory(:, 2));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_heelstrike_trajectory(:, 2));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_pushoff_trajectory(:, 2));
        plot3(phi_left_trajectory(left_heelstrike_indices), rho_left_trajectory(left_heelstrike_indices), feval(V_body_left_2_fit_heelstrike, [phi_left_trajectory(left_heelstrike_indices) rho_left_trajectory(left_heelstrike_indices)]), 'o');
        plot3(phi_left_trajectory(left_pushoff_indices), rho_left_trajectory(left_pushoff_indices), feval(V_body_left_2_fit_pushoff, [phi_left_trajectory(left_pushoff_indices) rho_left_trajectory(left_pushoff_indices)]), 'o');
        xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
%         legend('all', 'contact', 'heelstrike', 'pushoff', 'heelstrike fit', 'pushoff fit')
        uicontrol('style', 'text', 'position', [10 25 600 30], 'string', ['HS: R-square = ' num2str(V_body_left_2_gof_heelstrike.rsquare) ', RMS error = ' num2str(V_body_left_2_gof_heelstrike.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')
        uicontrol('style', 'text', 'position', [10 5 600 30], 'string', ['PO: R-square = ' num2str(V_body_left_2_gof_pushoff.rsquare) ', RMS error = ' num2str(V_body_left_2_gof_pushoff.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')

        figure; axes; hold on; title('v_3');
        plot3(phi_left_trajectory, rho_left_trajectory, V_body_left_ankle(:, 3));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_relevant_trajectory(:, 3));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_heelstrike_trajectory(:, 3));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_pushoff_trajectory(:, 3));
        plot3(phi_left_trajectory(left_heelstrike_indices), rho_left_trajectory(left_heelstrike_indices), feval(V_body_left_3_fit_heelstrike, [phi_left_trajectory(left_heelstrike_indices) rho_left_trajectory(left_heelstrike_indices)]), 'o');
        plot3(phi_left_trajectory(left_pushoff_indices), rho_left_trajectory(left_pushoff_indices), feval(V_body_left_3_fit_pushoff, [phi_left_trajectory(left_pushoff_indices) rho_left_trajectory(left_pushoff_indices)]), 'o');
        xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
%         legend('all', 'contact', 'heelstrike', 'pushoff', 'heelstrike fit', 'pushoff fit')
        uicontrol('style', 'text', 'position', [10 25 600 30], 'string', ['HS: R-square = ' num2str(V_body_left_3_gof_heelstrike.rsquare) ', RMS error = ' num2str(V_body_left_3_gof_heelstrike.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')
        uicontrol('style', 'text', 'position', [10 5 600 30], 'string', ['PO: R-square = ' num2str(V_body_left_3_gof_pushoff.rsquare) ', RMS error = ' num2str(V_body_left_3_gof_pushoff.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')

        figure; axes; hold on; title('v_4');
        plot3(phi_left_trajectory, rho_left_trajectory, V_body_left_ankle(:, 4));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_relevant_trajectory(:, 4));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_heelstrike_trajectory(:, 4));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_pushoff_trajectory(:, 4));
        plot3(phi_left_trajectory(left_heelstrike_indices), rho_left_trajectory(left_heelstrike_indices), feval(V_body_left_4_fit_heelstrike, [phi_left_trajectory(left_heelstrike_indices) rho_left_trajectory(left_heelstrike_indices)]), 'o');
        plot3(phi_left_trajectory(left_pushoff_indices), rho_left_trajectory(left_pushoff_indices), feval(V_body_left_4_fit_pushoff, [phi_left_trajectory(left_pushoff_indices) rho_left_trajectory(left_pushoff_indices)]), 'o');
        xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
%         legend('all', 'contact', 'heelstrike', 'pushoff', 'heelstrike fit', 'pushoff fit')
        uicontrol('style', 'text', 'position', [10 25 600 30], 'string', ['HS: R-square = ' num2str(V_body_left_4_gof_heelstrike.rsquare) ', RMS error = ' num2str(V_body_left_4_gof_heelstrike.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')
        uicontrol('style', 'text', 'position', [10 5 600 30], 'string', ['PO: R-square = ' num2str(V_body_left_4_gof_pushoff.rsquare) ', RMS error = ' num2str(V_body_left_4_gof_pushoff.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')

        figure; axes; hold on; title('v_5');
        plot3(phi_left_trajectory, rho_left_trajectory, V_body_left_ankle(:, 5));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_relevant_trajectory(:, 5));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_heelstrike_trajectory(:, 5));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_pushoff_trajectory(:, 5));
        plot3(phi_left_trajectory(left_heelstrike_indices), rho_left_trajectory(left_heelstrike_indices), feval(V_body_left_5_fit_heelstrike, [phi_left_trajectory(left_heelstrike_indices) rho_left_trajectory(left_heelstrike_indices)]), 'o');
        plot3(phi_left_trajectory(left_pushoff_indices), rho_left_trajectory(left_pushoff_indices), feval(V_body_left_5_fit_pushoff, [phi_left_trajectory(left_pushoff_indices) rho_left_trajectory(left_pushoff_indices)]), 'o');
        xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
%         legend('all', 'contact', 'heelstrike', 'pushoff', 'heelstrike fit', 'pushoff fit')
        uicontrol('style', 'text', 'position', [10 25 600 30], 'string', ['HS: R-square = ' num2str(V_body_left_5_gof_heelstrike.rsquare) ', RMS error = ' num2str(V_body_left_5_gof_heelstrike.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')
        uicontrol('style', 'text', 'position', [10 5 600 30], 'string', ['PO: R-square = ' num2str(V_body_left_5_gof_pushoff.rsquare) ', RMS error = ' num2str(V_body_left_5_gof_pushoff.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')

        figure; axes; hold on; title('v_6');
        plot3(phi_left_trajectory, rho_left_trajectory, V_body_left_ankle(:, 6));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_relevant_trajectory(:, 6));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_heelstrike_trajectory(:, 6));
        plot3(phi_left_trajectory, rho_left_trajectory, body_velocity_left_pushoff_trajectory(:, 6));
        plot3(phi_left_trajectory(left_heelstrike_indices), rho_left_trajectory(left_heelstrike_indices), feval(V_body_left_6_fit_heelstrike, [phi_left_trajectory(left_heelstrike_indices) rho_left_trajectory(left_heelstrike_indices)]), 'o');
        plot3(phi_left_trajectory(left_pushoff_indices), rho_left_trajectory(left_pushoff_indices), feval(V_body_left_6_fit_pushoff, [phi_left_trajectory(left_pushoff_indices) rho_left_trajectory(left_pushoff_indices)]), 'o');
        xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
%         legend('all', 'contact', 'heelstrike', 'pushoff', 'heelstrike fit', 'pushoff fit')
        uicontrol('style', 'text', 'position', [10 25 600 30], 'string', ['HS: R-square = ' num2str(V_body_left_6_gof_heelstrike.rsquare) ', RMS error = ' num2str(V_body_left_6_gof_heelstrike.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')
        uicontrol('style', 'text', 'position', [10 5 600 30], 'string', ['PO: R-square = ' num2str(V_body_left_6_gof_pushoff.rsquare) ', RMS error = ' num2str(V_body_left_6_gof_pushoff.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')
    end
    
    if plot_relevant_body_velocity_data_right
        % visualize - normalized body velocity vs. phi and rho
        figure; axes; hold on; title('v_1');
        plot3(phi_right_trajectory, rho_right_trajectory, V_body_right_ankle(:, 1));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_relevant_trajectory(:, 1));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_heelstrike_trajectory(:, 1));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_pushoff_trajectory(:, 1));
        plot3(phi_right_trajectory(right_heelstrike_indices), rho_right_trajectory(right_heelstrike_indices), feval(V_body_right_1_fit_heelstrike, [phi_right_trajectory(right_heelstrike_indices) rho_right_trajectory(right_heelstrike_indices)]), 'o');
        plot3(phi_right_trajectory(right_pushoff_indices), rho_right_trajectory(right_pushoff_indices), feval(V_body_right_1_fit_pushoff, [phi_right_trajectory(right_pushoff_indices) rho_right_trajectory(right_pushoff_indices)]), 'o');
        xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
%         legend('all', 'contact', 'heelstrike', 'pushoff', 'heelstrike fit', 'pushoff fit')
        uicontrol('style', 'text', 'position', [10 25 600 30], 'string', ['HS: R-square = ' num2str(V_body_right_1_gof_heelstrike.rsquare) ', RMS error = ' num2str(V_body_right_1_gof_heelstrike.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')
        uicontrol('style', 'text', 'position', [10 5 600 30], 'string', ['PO: R-square = ' num2str(V_body_right_1_gof_pushoff.rsquare) ', RMS error = ' num2str(V_body_right_1_gof_pushoff.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')

        figure; axes; hold on; title('v_2');
        plot3(phi_right_trajectory, rho_right_trajectory, V_body_right_ankle(:, 2));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_relevant_trajectory(:, 2));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_heelstrike_trajectory(:, 2));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_pushoff_trajectory(:, 2));
        plot3(phi_right_trajectory(right_heelstrike_indices), rho_right_trajectory(right_heelstrike_indices), feval(V_body_right_2_fit_heelstrike, [phi_right_trajectory(right_heelstrike_indices) rho_right_trajectory(right_heelstrike_indices)]), 'o');
        plot3(phi_right_trajectory(right_pushoff_indices), rho_right_trajectory(right_pushoff_indices), feval(V_body_right_2_fit_pushoff, [phi_right_trajectory(right_pushoff_indices) rho_right_trajectory(right_pushoff_indices)]), 'o');
        xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
%         legend('all', 'contact', 'heelstrike', 'pushoff', 'heelstrike fit', 'pushoff fit')
        uicontrol('style', 'text', 'position', [10 25 600 30], 'string', ['HS: R-square = ' num2str(V_body_right_2_gof_heelstrike.rsquare) ', RMS error = ' num2str(V_body_right_2_gof_heelstrike.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')
        uicontrol('style', 'text', 'position', [10 5 600 30], 'string', ['PO: R-square = ' num2str(V_body_right_2_gof_pushoff.rsquare) ', RMS error = ' num2str(V_body_right_2_gof_pushoff.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')

        figure; axes; hold on; title('v_3');
        plot3(phi_right_trajectory, rho_right_trajectory, V_body_right_ankle(:, 3));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_relevant_trajectory(:, 3));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_heelstrike_trajectory(:, 3));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_pushoff_trajectory(:, 3));
        plot3(phi_right_trajectory(right_heelstrike_indices), rho_right_trajectory(right_heelstrike_indices), feval(V_body_right_3_fit_heelstrike, [phi_right_trajectory(right_heelstrike_indices) rho_right_trajectory(right_heelstrike_indices)]), 'o');
        plot3(phi_right_trajectory(right_pushoff_indices), rho_right_trajectory(right_pushoff_indices), feval(V_body_right_3_fit_pushoff, [phi_right_trajectory(right_pushoff_indices) rho_right_trajectory(right_pushoff_indices)]), 'o');
        xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
%         legend('all', 'contact', 'heelstrike', 'pushoff', 'heelstrike fit', 'pushoff fit')
        uicontrol('style', 'text', 'position', [10 25 600 30], 'string', ['HS: R-square = ' num2str(V_body_right_3_gof_heelstrike.rsquare) ', RMS error = ' num2str(V_body_right_3_gof_heelstrike.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')
        uicontrol('style', 'text', 'position', [10 5 600 30], 'string', ['PO: R-square = ' num2str(V_body_right_3_gof_pushoff.rsquare) ', RMS error = ' num2str(V_body_right_3_gof_pushoff.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')

        figure; axes; hold on; title('v_4');
        plot3(phi_right_trajectory, rho_right_trajectory, V_body_right_ankle(:, 4));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_relevant_trajectory(:, 4));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_heelstrike_trajectory(:, 4));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_pushoff_trajectory(:, 4));
        plot3(phi_right_trajectory(right_heelstrike_indices), rho_right_trajectory(right_heelstrike_indices), feval(V_body_right_4_fit_heelstrike, [phi_right_trajectory(right_heelstrike_indices) rho_right_trajectory(right_heelstrike_indices)]), 'o');
        plot3(phi_right_trajectory(right_pushoff_indices), rho_right_trajectory(right_pushoff_indices), feval(V_body_right_4_fit_pushoff, [phi_right_trajectory(right_pushoff_indices) rho_right_trajectory(right_pushoff_indices)]), 'o');
        xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
%         legend('all', 'contact', 'heelstrike', 'pushoff', 'heelstrike fit', 'pushoff fit')
        uicontrol('style', 'text', 'position', [10 25 600 30], 'string', ['HS: R-square = ' num2str(V_body_right_4_gof_heelstrike.rsquare) ', RMS error = ' num2str(V_body_right_4_gof_heelstrike.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')
        uicontrol('style', 'text', 'position', [10 5 600 30], 'string', ['PO: R-square = ' num2str(V_body_right_4_gof_pushoff.rsquare) ', RMS error = ' num2str(V_body_right_4_gof_pushoff.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')

        figure; axes; hold on; title('v_5');
        plot3(phi_right_trajectory, rho_right_trajectory, V_body_right_ankle(:, 5));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_relevant_trajectory(:, 5));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_heelstrike_trajectory(:, 5));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_pushoff_trajectory(:, 5));
        plot3(phi_right_trajectory(right_heelstrike_indices), rho_right_trajectory(right_heelstrike_indices), feval(V_body_right_5_fit_heelstrike, [phi_right_trajectory(right_heelstrike_indices) rho_right_trajectory(right_heelstrike_indices)]), 'o');
        plot3(phi_right_trajectory(right_pushoff_indices), rho_right_trajectory(right_pushoff_indices), feval(V_body_right_5_fit_pushoff, [phi_right_trajectory(right_pushoff_indices) rho_right_trajectory(right_pushoff_indices)]), 'o');
        xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
%         legend('all', 'contact', 'heelstrike', 'pushoff', 'heelstrike fit', 'pushoff fit')
        uicontrol('style', 'text', 'position', [10 25 600 30], 'string', ['HS: R-square = ' num2str(V_body_right_5_gof_heelstrike.rsquare) ', RMS error = ' num2str(V_body_right_5_gof_heelstrike.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')
        uicontrol('style', 'text', 'position', [10 5 600 30], 'string', ['PO: R-square = ' num2str(V_body_right_5_gof_pushoff.rsquare) ', RMS error = ' num2str(V_body_right_5_gof_pushoff.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')

        figure; axes; hold on; title('v_6');
        plot3(phi_right_trajectory, rho_right_trajectory, V_body_right_ankle(:, 6));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_relevant_trajectory(:, 6));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_heelstrike_trajectory(:, 6));
        plot3(phi_right_trajectory, rho_right_trajectory, body_velocity_right_pushoff_trajectory(:, 6));
        plot3(phi_right_trajectory(right_heelstrike_indices), rho_right_trajectory(right_heelstrike_indices), feval(V_body_right_6_fit_heelstrike, [phi_right_trajectory(right_heelstrike_indices) rho_right_trajectory(right_heelstrike_indices)]), 'o');
        plot3(phi_right_trajectory(right_pushoff_indices), rho_right_trajectory(right_pushoff_indices), feval(V_body_right_6_fit_pushoff, [phi_right_trajectory(right_pushoff_indices) rho_right_trajectory(right_pushoff_indices)]), 'o');
        xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
%         legend('all', 'contact', 'heelstrike', 'pushoff', 'heelstrike fit', 'pushoff fit')
        uicontrol('style', 'text', 'position', [10 25 600 30], 'string', ['HS: R-square = ' num2str(V_body_right_6_gof_heelstrike.rsquare) ', RMS error = ' num2str(V_body_right_6_gof_heelstrike.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')
        uicontrol('style', 'text', 'position', [10 5 600 30], 'string', ['PO: R-square = ' num2str(V_body_right_6_gof_pushoff.rsquare) ', RMS error = ' num2str(V_body_right_6_gof_pushoff.rmse)], 'fontsize', 14, 'HorizontalAlignment', 'left')
    end
end
    
if false % was: find_phi_constraint
%     V_phi_body_trajectory = zeros(number_of_time_steps, 6);
%     V_phi_body_direction_trajectory = zeros(number_of_time_steps, 6);
%     V_phi_body_normed_trajectory = zeros(number_of_time_steps, 6);
%     body_velocity_normed_trajectory = zeros(number_of_time_steps, 6);
    
    time_steps_to_fit = 1 : 3000;
%     for i_time = 1 : number_of_time_steps
%     for i_time = time_steps_to_fit
%         phi_dot = phi_dot_left_relevant_trajectory(i_time);
%         
%         % calculate z-offset of rho rotation center
%         z_offset = feval(rho_cor_fit, phi_left_trajectory(i_time));
%         rho_rotation_center_local = [0; 0; z_offset; 1];
%         
%         % calculate V_rho_body
%         plant.jointAngles = angle_trajectories(i_time, :)';
%         plant.updateKinematics;
%         J_body = plant.bodyJacobians{3};
%         ankle_transformation = reshape(T_left_ankle_to_world_trajectory(i_time, :), 4, 4);
%         rho_rotation_center_world = ankle_transformation * rho_rotation_center_local;
%         
%         p_virtual_contact = [rho_rotation_center_world(1:2); 0; 1];
%         J_body_virtual_contact = plant.calculateArbitraryFrameBodyJacobian([eye(4, 3) p_virtual_contact], 12);
%         A_contact = J_body_virtual_contact([1:3 6], :);
%         
%         P_A_contact = A_contact' * (A_contact * A_contact')^(-1) * A_contact; % projection onto the space spanned by the columns of A
%         P_A_contact_orth = eye(number_of_joints) - P_A_contact; % projection onto the null space of A
% 
%         % ACHTUNG: this is hard-coded for the left leg now, change it later (left ankle in-eversion is 12th joint)
%         joint_velocity_rho_change_unconstrained = [zeros(11, 1); 1; zeros(26, 1)];
%         joint_velocity_rho_change_constrained = (P_A_contact_orth * joint_velocity_rho_change_unconstrained);
%         body_velocity_rho_change_constrained = J_body * joint_velocity_rho_change_constrained;
% %         V_rho_body = body_velocity_rho_change_constrained * 1 / body_velocity_rho_change_constrained(5) * rho_dot;
%         body_velocity_rho = V_body_left_ankle(i_time, 4);
%         V_rho_body = body_velocity_rho_change_constrained * 1 / body_velocity_rho_change_constrained(4) * body_velocity_rho;
%         
%         % calculate V_phi_body
%         V_total_body = V_body_left_ankle(i_time, :)';
%         V_phi_body = V_total_body - V_rho_body;
%         V_phi_body_direction = V_phi_body * phi_dot^(-1);
%         
%         V_phi_body_trajectory(i_time, :) = V_phi_body;
%         V_phi_body_direction_trajectory(i_time, :) = V_phi_body_direction;
%         V_phi_body_normed_trajectory(i_time, :) = normVector(V_phi_body_direction);
%         
%         body_velocity_normed_trajectory(i_time, :) = normVector(V_body_left_ankle(i_time, :));
%         
%         if i_time/10  == floor(i_time / 10)
%             disp(num2str(i_time))
%         end
%     end
    
    % normalize direction
%     V_phi_body_normed_trajectory(V_phi_body_normed_trajectory(:, 1)<0, :) = - V_phi_body_normed_trajectory(V_phi_body_normed_trajectory(:, 1)<0, :);
%     body_velocity_normed_trajectory(body_velocity_normed_trajectory(:, 1)<0, :) = - body_velocity_normed_trajectory(body_velocity_normed_trajectory(:, 1)<0, :);

%     relevant_data_points = abs(phi_dot_trajectory_relevant) > 1.0;
%     irrelevant_data_points = ~relevant_data_points;
%     phi_dot_trajectory_double_relevant = phi_dot_trajectory_relevant; phi_dot_trajectory_double_relevant(irrelevant_data_points) = NaN;
%     body_velocity_trajectory_relevant = V_body_left_ankle; body_velocity_trajectory_relevant(irrelevant_data_points_left, :) = NaN;
%     body_velocity_normed_trajectory_relevant = body_velocity_normed_trajectory; body_velocity_normed_trajectory_relevant(irrelevant_data_points_left, :) = NaN;
%     V_phi_body_trajectory_relevant = V_phi_body_trajectory; V_phi_body_trajectory_relevant(irrelevant_data_points_left, :) = NaN;
%     V_phi_body_direction_trajectory_relevant = V_phi_body_direction_trajectory; V_phi_body_direction_trajectory_relevant(irrelevant_data_points_left, :) = NaN;
%     V_phi_body_normed_trajectory_relevant = V_phi_body_normed_trajectory; V_phi_body_normed_trajectory_relevant(irrelevant_data_points_left, :) = NaN;

    %%
%     figure; axes; hold on; title('v_1');
%     plot(V_phi_body_trajectory(time_steps_to_fit, 1));
%     plot(V_phi_body_trajectory_relevant(time_steps_to_fit, 1));
%     plot(V_phi_body_direction_trajectory_relevant(time_steps_to_fit, 1), 'linewidth', 2);
%     plot(phi_dot_left_trajectory(time_steps_to_fit, 1));
% %     plot(V_phi_body_normed_trajectory_relevant(time_steps_to_fit, 1));
%     legend('body', 'normed');
% %     legend('body', 'direction', 'normed');

%     figure; axes; hold on; title('v_2');
%     plot(V_phi_body_trajectory(time_steps_to_fit, 2));
%     plot(V_phi_body_trajectory_relevant(time_steps_to_fit, 2));
%     plot(V_phi_body_direction_trajectory_relevant(time_steps_to_fit, 2), 'linewidth', 2);
%     plot(V_phi_body_normed_trajectory_relevant(time_steps_to_fit, 2));
%     
%     figure; axes; hold on; title('v_3');
%     plot(V_phi_body_trajectory(time_steps_to_fit, 3));
%     plot(V_phi_body_trajectory_relevant(time_steps_to_fit, 3));
%     plot(V_phi_body_direction_trajectory_relevant(time_steps_to_fit, 3), 'linewidth', 2);
%     plot(V_phi_body_normed_trajectory_relevant(time_steps_to_fit, 3));
%     
%     figure; axes; hold on; title('v_4');
%     plot(V_phi_body_trajectory(time_steps_to_fit, 4));
%     plot(V_phi_body_trajectory_relevant(time_steps_to_fit, 4));
%     plot(V_phi_body_direction_trajectory_relevant(time_steps_to_fit, 4), 'linewidth', 2);
%     plot(V_phi_body_normed_trajectory_relevant(time_steps_to_fit, 4));
%     
%     figure; axes; hold on; title('v_5');
%     plot(V_phi_body_trajectory(time_steps_to_fit, 5));
%     plot(V_phi_body_trajectory_relevant(time_steps_to_fit, 5));
%     plot(V_phi_body_direction_trajectory_relevant(time_steps_to_fit, 5), 'linewidth', 2);
%     plot(V_phi_body_normed_trajectory_relevant(time_steps_to_fit, 5));
%     
%     figure; axes; hold on; title('v_6');
%     plot(V_phi_body_trajectory(time_steps_to_fit, 6));
%     plot(V_phi_body_trajectory_relevant(time_steps_to_fit, 6));
%     plot(V_phi_body_direction_trajectory_relevant(time_steps_to_fit, 6), 'linewidth', 2);
%     plot(V_phi_body_normed_trajectory_relevant(time_steps_to_fit, 6));
%     return
    %%
    
    
%     V_phi_body_direction_norm_trajectory = sum(V_phi_body_direction_trajectory.^2, 2).^(0.5);
%     relevant_data_points = V_phi_body_direction_norm_trajectory > 0.1;
    
%     figure; axes; hold on
%     plot(phi_dot_trajectory);
%     plot(phi_dot_trajectory_relevant);
%     plot(phi_dot_trajectory_double_relevant);
%     
%     figure; axes; hold on
%     plot(rho_dot_trajectory);
%     plot(rho_dot_trajectory_relevant);
%     
%     plot(V_phi_body_normed_trajectory(:, 1));
%     plot(V_phi_body_normed_trajectory_relevant(:, 1));
    
    

    return
    
    
%     figure; axes; hold on; title('v_1');
%     plot3(phi_trajectory, rho_trajectory, V_phi_body_trajectory(:, 1));
%     plot3(phi_trajectory, rho_trajectory, V_phi_body_trajectory_relevant(:, 1), 'linewidth', 2);
% %     plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 1));
%     plot3(phi_trajectory, rho_trajectory, V_phi_body_normed_trajectory_relevant(:, 1));
%     xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
%     
%     figure; axes; hold on; title('v_2');
% %     plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 2), 'linewidth', 2);
%     plot3(phi_trajectory, rho_trajectory, V_phi_body_normed_trajectory(:, 2));
%     plot3(phi_trajectory, rho_trajectory, body_velocity_normed_trajectory(:, 2));
%     xlabel('\phi'); ylabel('\rho'); zlabel('v_2')
%     
%     figure; axes; hold on; title('v_3');
% %     plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 3), 'linewidth', 2);
%     plot3(phi_trajectory, rho_trajectory, V_phi_body_normed_trajectory(:, 3));
%     plot3(phi_trajectory, rho_trajectory, body_velocity_normed_trajectory(:, 3));
%     xlabel('\phi'); ylabel('\rho'); zlabel('v_3')
%     
%     figure; axes; hold on; title('v_4');
% %     plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 4), 'linewidth', 2);
%     plot3(phi_trajectory, rho_trajectory, V_phi_body_normed_trajectory(:, 4));
%     plot3(phi_trajectory, rho_trajectory, body_velocity_normed_trajectory(:, 4));
%     xlabel('\phi'); ylabel('\rho'); zlabel('v_4')
%     
%     figure; axes; hold on; title('v_5');
% %     plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 5), 'linewidth', 2);
%     plot3(phi_trajectory, rho_trajectory, V_phi_body_normed_trajectory(:, 5));
%     plot3(phi_trajectory, rho_trajectory, body_velocity_normed_trajectory(:, 5));
%     xlabel('\phi'); ylabel('\rho'); zlabel('v_5')
%     
%     figure; axes; hold on; title('v_6');
% %     plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 6), 'linewidth', 2);
%     plot3(phi_trajectory, rho_trajectory, V_phi_body_normed_trajectory(:, 6));
%     plot3(phi_trajectory, rho_trajectory, body_velocity_normed_trajectory(:, 6));
%     xlabel('\phi'); ylabel('\rho'); zlabel('v_6')
    
    
    % approximate constraints
%     data_points_for_fit = 1 : 1 : length(phi_trajectory);
    data_points_for_fit = time_steps_to_fit;
    data_points_for_fit(isnan(phi_trajectory(data_points_for_fit))) = [];
    data_points_for_fit(isnan(rho_trajectory(data_points_for_fit))) = [];
    data_points_for_fit(isnan(V_phi_body_normed_trajectory(data_points_for_fit))) = [];
    V_phi_body_1_fit = fit([phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 1), 'poly55');
    V_phi_body_2_fit = fit([phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 2), 'poly55');
    V_phi_body_3_fit = fit([phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 3), 'poly55');
    V_phi_body_4_fit = fit([phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 4), 'poly55');
    V_phi_body_5_fit = fit([phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 5), 'poly55');
    V_phi_body_6_fit = fit([phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 6), 'poly55');

%     figure; plot(V_phi_body_1_fit, [phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 1));
%     figure; plot(V_phi_body_2_fit, [phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 2));
%     figure; plot(V_phi_body_3_fit, [phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 3));
%     figure; plot(V_phi_body_4_fit, [phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 4));
%     figure; plot(V_phi_body_5_fit, [phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 5));
%     figure; plot(V_phi_body_6_fit, [phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 6));
    
end        


%% save results
if save_results
%     constraints_file_name = makeFileName(date, subject_id, 'walking', i_trial, 'constraints');
    constraints_file_name = makeFileName(date, subject_id, 'constraints');
    save ...
      ( ...
        constraints_file_name, ...
        'V_body_left_1_fit_heelstrike', ...
        'V_body_left_2_fit_heelstrike', ...
        'V_body_left_3_fit_heelstrike', ...
        'V_body_left_4_fit_heelstrike', ...
        'V_body_left_5_fit_heelstrike', ...
        'V_body_left_6_fit_heelstrike', ...
        'V_body_left_1_fit_pushoff', ...
        'V_body_left_2_fit_pushoff', ...
        'V_body_left_3_fit_pushoff', ...
        'V_body_left_4_fit_pushoff', ...
        'V_body_left_5_fit_pushoff', ...
        'V_body_left_6_fit_pushoff', ...
        'V_body_right_1_fit_heelstrike', ...
        'V_body_right_2_fit_heelstrike', ...
        'V_body_right_3_fit_heelstrike', ...
        'V_body_right_4_fit_heelstrike', ...
        'V_body_right_5_fit_heelstrike', ...
        'V_body_right_6_fit_heelstrike', ...
        'V_body_right_1_fit_pushoff', ...
        'V_body_right_2_fit_pushoff', ...
        'V_body_right_3_fit_pushoff', ...
        'V_body_right_4_fit_pushoff', ...
        'V_body_right_5_fit_pushoff', ...
        'V_body_right_6_fit_pushoff' ...
      );
    
end


    %% find rho constraint by phi
    if find_rho_constraint_by_phi
        
        % fit surfaces to phi-rho-z data
        data_points_for_fit_left = find(left_contact_indicators_mocap);
        phi_left_fitpoints = phi_left_trajectory(data_points_for_fit_left);
        rho_left_fitpoints = rho_left_trajectory(data_points_for_fit_left);
        poi_z_left_fits = cell(1, number_of_points_of_interest);
        for i_poi = 1 : number_of_points_of_interest
            p_z = left_points_of_interest_world_trajectories_z(data_points_for_fit_left, i_poi);
            poi_z_left_fits{i_poi} = fit([phi_left_fitpoints, rho_left_fitpoints], p_z, 'poly55');
%             figure; plot(poi_z_left_fits{i_poi}, [phi_left_fitpoints, rho_left_fitpoints], left_points_of_interest_world_trajectories_z(data_points_for_fit_left, i_poi));
%             xlabel('\phi'); ylabel('\rho'); zlabel('z'); 
        end
        
        % for a grid of each phi values, find the point of interest that moves least under changes of rho
        phi_left_values = linspace(min(phi_left_fitpoints), max(phi_left_fitpoints), number_of_phi_values);
        left_z_values = cell(number_of_phi_values, number_of_points_of_interest);
        poi_z_left_fits_by_rho_rms = zeros(number_of_phi_values, number_of_points_of_interest);
        
%         figure; axes; hold on
%         plot(rho_trajectory);
        phi_window_radius = 0.025;
        for i_phi = 1 : length(phi_left_values)
            % for each point of interest, compare z-values to straight line
            rho_windows = zeros(number_of_phi_values, 2);
            phi_window = phi_left_values(i_phi) + [-1 1] * phi_window_radius;
            phi_window_indices = find(phi_window(1) <= phi_left_trajectory & phi_left_trajectory <= phi_window(2));
%             plot(phi_window_indices, rho_trajectory(phi_window_indices), '-')                
            rho_windows(i_phi, :) = [min(rho_left_trajectory(phi_window_indices)) max(rho_left_trajectory(phi_window_indices))];
%             plot([i_phi i_phi], rho_windows(i_phi, :));
            rho_values = linspace(rho_windows(i_phi, 1), rho_windows(i_phi, 2), 20);

    %         figure; hold on
            for i_poi = 1 : number_of_points_of_interest
                left_z_values{i_phi, i_poi} = feval(poi_z_left_fits{i_poi}, [ones(length(rho_values), 1)*phi_left_values(i_phi), rho_values']);
                z_values_rms = sum((left_z_values{i_phi, i_poi} - mean(left_z_values{i_phi, i_poi})).^2).^(0.5);
                poi_z_left_fits_by_rho_rms(i_phi, i_poi) = z_values_rms;
            end



            % plot z-curves
%             figure; axes; hold on; title(['\phi = ' num2str(phi_values(i_phi))]);
%             for i_poi = 1 : number_of_points_of_interest
%                 plot(rho_values, z_values{i_phi, i_poi});
%                 plot(rho_values, mean(z_values{i_phi, i_poi})*ones(size(rho_values)));
%             end
        end

%         figure; plot(poi_z_fits_by_rho_rms', '-x'); title('different colors are different values of \phi')
%         xlabel('point of interest'); ylabel('error')


%         [~, poi_z_fits_by_rho_derivative_error_minimum_indices] = min(poi_z_fits_by_rho_derivative_error, [], 2);
        [poi_z_fits_by_rho_rms_minima, poi_z_fits_by_rho_rms_minimum_indices] = min(poi_z_left_fits_by_rho_rms, [], 2);

        options = fitoptions('smoothingspline');
        options.SmoothingParam = 0.9999;
        rho_cor_fit = fit(phi_left_values', points_of_interest_y_local(poi_z_fits_by_rho_rms_minimum_indices)', 'smoothingspline', options);
        
        % plot minimum by phi value
        figure; plot(rho_cor_fit, phi_left_values, points_of_interest_y_local(poi_z_fits_by_rho_rms_minimum_indices));
        xlabel('\phi'); ylabel('z'); 

        figure; plot(phi_left_values, poi_z_fits_by_rho_rms_minima);
        xlabel('\phi'); ylabel('rho rms minimum'); 

%%


%         % visualize
%         phi_values = -1 : 0.1 : 0;
%         rho_values = -0.3 : 0.1 : 0.3;
%         [phi_grid, rho_grid] = meshgrid(phi_values, rho_values);
%         poi_z_grids = cell(1, number_of_points_of_interest);
%         for i_poi = 1 : number_of_points_of_interest
%             poi_z_grids{i_poi} = feval(poi_z_fits{i_poi}, phi_grid, rho_grid);
%         end
%         pois_to_plot = 1 : number_of_points_of_interest;
%         figure; axes; hold on;
%         xlabel('phi'); ylabel('rho'); 
%         for i_poi = pois_to_plot
%             surf(phi_grid, rho_grid, poi_z_grids{i_poi});
%             plot3(phi_trajectory, rho_trajectory, points_of_interest_world_trajectories_z(:, i_poi));
%         end
% 
%         % fix a phi value
%         phi = -0.5;
% 
%         % evaluate changes of z_i by rho for this phi value
%         rho_values = (-0.3 : 0.01 : 0.3)';
%         phi_values = phi * ones(size(rho_values));
%         z_values = cell(1, number_of_points_of_interest);
%         for i_poi = 1 : number_of_points_of_interest
%             z_values{i_poi} = feval(poi_z_fits{i_poi}, [phi_values, rho_values]);
%         end
% 
%         figure; axes; hold on; axis equal
%         for i_poi = 1 : number_of_points_of_interest
%             plot(rho_values, z_values{i_poi});
%         end



        % fit surfaces to phi-rho-z data
        data_points_for_fit_right = find(right_contact_indicators_mocap);
        phi_right_fitpoints = phi_right_trajectory(data_points_for_fit_right);
        rho_right_fitpoints = rho_right_trajectory(data_points_for_fit_right);
        poi_z_right_fits = cell(1, number_of_points_of_interest);
        for i_poi = 1 : number_of_points_of_interest
            p_z = right_points_of_interest_world_trajectories_z(data_points_for_fit_right, i_poi);
            poi_z_right_fits{i_poi} = fit([phi_right_fitpoints, rho_right_fitpoints], p_z, 'poly55');
%             figure; plot(poi_z_right_fits{i_poi}, [phi_right_fitpoints, rho_right_fitpoints], right_points_of_interest_world_trajectories_z(data_points_for_fit_right, i_poi));
%             xlabel('\phi'); ylabel('\rho'); zlabel('z'); 
        end
        
        % for a grid of each phi values, find the point of interest that moves least under changes of rho
        phi_right_values = linspace(min(phi_right_fitpoints), max(phi_right_fitpoints), number_of_phi_values);
        left_z_values = cell(number_of_phi_values, number_of_points_of_interest);
        poi_z_left_fits_by_rho_rms = zeros(number_of_phi_values, number_of_points_of_interest);
        
%         figure; axes; hold on
%         plot(rho_trajectory);
        phi_window_radius = 0.025;
        for i_phi = 1 : length(phi_left_values)
            % for each point of interest, compare z-values to straight line
            rho_windows = zeros(number_of_phi_values, 2);
            phi_window = phi_right_values(i_phi) + [-1 1] * phi_window_radius;
            phi_window_indices = find(phi_window(1) <= phi_right_trajectory & phi_right_trajectory <= phi_window(2));
%             plot(phi_window_indices, rho_trajectory(phi_window_indices), '-')                
            rho_windows(i_phi, :) = [min(rho_right_trajectory(phi_window_indices)) max(rho_right_trajectory(phi_window_indices))];
%             plot([i_phi i_phi], rho_windows(i_phi, :));
            rho_values = linspace(rho_windows(i_phi, 1), rho_windows(i_phi, 2), 20);

    %         figure; hold on
            for i_poi = 1 : number_of_points_of_interest
                right_z_values{i_phi, i_poi} = feval(poi_z_right_fits{i_poi}, [ones(length(rho_values), 1)*phi_right_values(i_phi), rho_values']);
                z_values_rms = sum((right_z_values{i_phi, i_poi} - mean(right_z_values{i_phi, i_poi})).^2).^(0.5);
                poi_z_right_fits_by_rho_rms(i_phi, i_poi) = z_values_rms;
            end



            % plot z-curves
            figure; axes; hold on; title(['\phi = ' num2str(phi_values(i_phi))]);
            for i_poi = 1 : number_of_points_of_interest
                plot(rho_values, z_values{i_phi, i_poi});
                plot(rho_values, mean(z_values{i_phi, i_poi})*ones(size(rho_values)));
            end
        end

%         figure; plot(poi_z_fits_by_rho_rms', '-x'); title('different colors are different values of \phi')
%         xlabel('point of interest'); ylabel('error')


%         [~, poi_z_fits_by_rho_derivative_error_minimum_indices] = min(poi_z_fits_by_rho_derivative_error, [], 2);
        [poi_z_fits_by_rho_rms_minima, poi_z_fits_by_rho_rms_minimum_indices] = min(poi_z_right_fits_by_rho_rms, [], 2);

        options = fitoptions('smoothingspline');
        options.SmoothingParam = 0.9999;
        rho_cor_fit = fit(phi_right_values', points_of_interest_y_local(poi_z_fits_by_rho_rms_minimum_indices)', 'smoothingspline', options);
        
        % plot minimum by phi value
        figure; plot(rho_cor_fit, phi_right_values, points_of_interest_y_local(poi_z_fits_by_rho_rms_minimum_indices));
        xlabel('\phi'); ylabel('z'); 

        figure; plot(phi_right_values, poi_z_fits_by_rho_rms_minima);
        xlabel('\phi'); ylabel('rho rms minimum'); 
    end    
    
    

    
    































return
%% old
    
    ankle_z_trajectory = joint_frame_trajectory(:, 15);
    
    %% calculate elevation angle
    ankle_x_azimuth_trajectory = zeros(number_of_time_steps, 1);
    ankle_x_elevation_trajectory = zeros(number_of_time_steps, 1);
    ankle_z_azimuth_trajectory = zeros(number_of_time_steps, 1);
    ankle_z_elevation_trajectory = zeros(number_of_time_steps, 1);

    for i_time = 1:number_of_time_steps
        T_ankle_to_world_current = reshape(joint_frame_trajectory(i_time, :), 4, 4);
        
        % calculate spherical coordinates of the ankle frame x- and z-axes
        ankle_frame_x_axis_current = T_ankle_to_world_current(1:3, 1);
        [ankle_x_azimuth_trajectory(i_time), ankle_x_elevation_trajectory(i_time)] = cart2sph(ankle_frame_x_axis_current(1), ankle_frame_x_axis_current(2), ankle_frame_x_axis_current(3));
        ankle_frame_z_axis_current = T_ankle_to_world_current(1:3, 3);
        [ankle_z_azimuth_trajectory(i_time), ankle_z_elevation_trajectory(i_time)] = cart2sph(ankle_frame_z_axis_current(1), ankle_frame_z_axis_current(2), ankle_frame_z_axis_current(3));
    end
    
    phi = ankle_x_elevation_trajectory;
    phi_dot = deriveByTime(phi, sampling_rate^(-1));
    rho = ankle_z_elevation_trajectory;
    rho_dot = deriveByTime(rho, sampling_rate^(-1));
    % body velocity is defined as g^(-1) * \dot g
    % I am interested in change by the elevation angle phi
    % so I use the chain rule dg/dt = dg/dphi * dphi/dt
    % then multiply by (dphi/dt)^(-1)
    % see note from 4.4.2016
    d_g_by_d_phi_vee_trajectory = zeros(number_of_time_steps, 6);
    d_g_by_d_rho_vee_trajectory = zeros(number_of_time_steps, 6);
    for i_time = 1:number_of_time_steps
        T_ankle_to_world_current = reshape(joint_frame_trajectory(i_time, :), 4, 4); % this is g
        body_velocity = body_velocity_trajectory(i_time, :);
        body_velocity_wedge = wedgeTwist(body_velocity);

        d_g_by_d_phi = T_ankle_to_world_current * body_velocity_wedge * phi_dot(i_time)^(-1);
        d_g_by_d_phi_vee = veeTwist(d_g_by_d_phi);
        d_g_by_d_phi_vee_trajectory(i_time, :) = d_g_by_d_phi_vee;
        
        d_g_by_d_rho = T_ankle_to_world_current * body_velocity_wedge * rho_dot(i_time)^(-1);
        d_g_by_d_rho_vee = veeTwist(d_g_by_d_rho);
        d_g_by_d_rho_vee_trajectory(i_time, :) = d_g_by_d_rho_vee;
    end    
    
    % psi is old and probably useless
    psi = body_velocity_trajectory .* repmat(phi_dot.^(-1), 1, 6);
    velocity_threshold = 0.5;
    velocity_threshold = 1.5;

%     relevant_data_points = abs(phi_dot)>velocity_threshold & contact_indicator_trajectory;
    relevant_data_points = contact_indicator_trajectory;
    irrelevant_data_points = ~relevant_data_points;
    phi_relevant = phi;         phi_relevant(irrelevant_data_points) = NaN;
    phi_dot_relevant = phi_dot; phi_dot_relevant(irrelevant_data_points) = NaN;
    rho_relevant = rho;         rho_relevant(irrelevant_data_points) = NaN;
    psi_relevant = psi;         psi_relevant(irrelevant_data_points, :) = NaN;
    body_velocity_trajectory_relevant = body_velocity_trajectory;       body_velocity_trajectory_relevant(irrelevant_data_points, :) = NaN;
    d_g_by_d_phi_vee_relevant = d_g_by_d_phi_vee_trajectory;            d_g_by_d_phi_vee_relevant(irrelevant_data_points, :) = NaN;
    d_g_by_d_rho_vee_relevant = d_g_by_d_rho_vee_trajectory;            d_g_by_d_rho_vee_relevant(irrelevant_data_points, :) = NaN;
    
    nan_data_points = any(isnan(psi_relevant), 2);
    phi_nanless = phi(~nan_data_points);
    psi_nanless = psi(~nan_data_points, :);


    
    % visualize to check - time
    figure; axes; hold on
    plot(phi);
    plot(rho);
    contact_indicator_trajectory_nanned = contact_indicator_trajectory; contact_indicator_trajectory_nanned(contact_indicator_trajectory_nanned==0) = NaN;
    plot(contact_indicator_trajectory_nanned-1, 'linewidth', 5);
    plot(ankle_z_trajectory);
    legend('phi', 'rho', 'contact', 'ankle_z')
    
    % visualize - ankle_z vs. phi and rho
    figure; axes; hold on; title('v_1');
    plot3(phi, rho, ankle_z_trajectory(:, 1));
    plot3(phi_relevant, rho_relevant, ankle_z_trajectory(:, 1));
    xlabel('\phi'); ylabel('\rho'); zlabel('ankle_z')
    
    % visualize - body velocity vs. phi and rho
    figure; axes; hold on; title('v_1');
    plot3(phi, rho, body_velocity_trajectory(:, 1));
    plot3(phi_relevant, rho_relevant, body_velocity_trajectory_relevant(:, 1));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_1')
    figure; axes; hold on; title('v_2');
    plot3(phi, rho, body_velocity_trajectory(:, 2));
    plot3(phi_relevant, rho_relevant, body_velocity_trajectory_relevant(:, 2));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_2')
    figure; axes; hold on; title('v_3');
    plot3(phi, rho, body_velocity_trajectory(:, 3));
    plot3(phi_relevant, rho_relevant, body_velocity_trajectory_relevant(:, 3));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_3')
%     figure; axes; hold on; title('v_4');
%     plot3(phi_relevant, rho_relevant, body_velocity_trajectory_relevant(:, 4));
%     figure; axes; hold on; title('v_5');
%     plot3(phi_relevant, rho_relevant, body_velocity_trajectory_relevant(:, 5));
%     figure; axes; hold on; title('v_6');
%     plot3(phi_relevant, rho_relevant, body_velocity_trajectory_relevant(:, 6));
    
%     % visualize - psi vs. phi
%     figure; axes; hold on; title('\psi_1');
%     plot(phi_relevant, psi_relevant(:, 1));
%     figure; axes; hold on; title('\psi_2');
%     plot(phi_relevant, psi_relevant(:, 2));
%     figure; axes; hold on; title('\psi_3');
%     plot(phi_relevant, psi_relevant(:, 3));
%     figure; axes; hold on; title('\psi_4');
%     plot(phi_relevant, psi_relevant(:, 4));
%     figure; axes; hold on; title('\psi_5');
%     plot(phi_relevant, psi_relevant(:, 5));
%     figure; axes; hold on; title('\psi_6');
%     plot(phi_relevant, psi_relevant(:, 6));

%     % visualize - body velocity vs. phi
%     figure; axes; hold on; title('v_1');
%     plot(phi_relevant, d_g_by_d_phi_vee_relevant(:, 1));
%     figure; axes; hold on; title('v_2');
%     plot(phi_relevant, d_g_by_d_phi_vee_relevant(:, 2));
%     figure; axes; hold on; title('v_3');
%     plot(phi_relevant, d_g_by_d_phi_vee_relevant(:, 3));
%     figure; axes; hold on; title('v_4');
%     plot(phi_relevant, d_g_by_d_phi_vee_relevant(:, 4));
%     figure; axes; hold on; title('v_5');
%     plot(phi_relevant, d_g_by_d_phi_vee_relevant(:, 5));
%     figure; axes; hold on; title('v_6');
%     plot(phi_relevant, d_g_by_d_phi_vee_relevant(:, 6));
%     
%     % visualize - body velocity vs. rho
%     figure; axes; hold on; title('v_1');
%     plot(rho_relevant, d_g_by_d_rho_vee_relevant(:, 1));
%     figure; axes; hold on; title('v_2');
%     plot(rho_relevant, d_g_by_d_rho_vee_relevant(:, 2));
%     figure; axes; hold on; title('v_3');
%     plot(rho_relevant, d_g_by_d_rho_vee_relevant(:, 3));
%     figure; axes; hold on; title('v_4');
%     plot(rho_relevant, d_g_by_d_rho_vee_relevant(:, 4));
%     figure; axes; hold on; title('v_5');
%     plot(rho_relevant, d_g_by_d_rho_vee_relevant(:, 5));
%     figure; axes; hold on; title('v_6');
%     plot(rho_relevant, d_g_by_d_rho_vee_relevant(:, 6));
%     distFig
    
%     % visualize - body velocity vs. phi and rho
%     figure; axes; hold on; title('v_1');
%     plot3(phi_relevant, rho_relevant, d_g_by_d_phi_vee_relevant(:, 1));
%     figure; axes; hold on; title('v_2');
%     plot3(phi_relevant, rho_relevant, d_g_by_d_phi_vee_relevant(:, 2));
%     figure; axes; hold on; title('v_3');
%     plot3(phi_relevant, rho_relevant, d_g_by_d_phi_vee_relevant(:, 3));
%     figure; axes; hold on; title('v_4');
%     plot3(phi_relevant, rho_relevant, d_g_by_d_phi_vee_relevant(:, 4));
%     figure; axes; hold on; title('v_5');
%     plot3(phi_relevant, rho_relevant, d_g_by_d_phi_vee_relevant(:, 5));
%     figure; axes; hold on; title('v_6');
%     plot3(phi_relevant, rho_relevant, d_g_by_d_phi_vee_relevant(:, 6));



    phi_trajectory = phi_relevant;
    rho_trajectory = rho_relevant;
    body_velocity_trajectory = body_velocity_trajectory_relevant;
    return
    
    % polynomial fits
    fit_order = 4;
    phi_values = linspace(min(phi) - 0.05*(max(phi) - min(phi)), max(phi) +  0.05*(max(phi) - min(phi)), 200);
%     phi_values = linspace(-0.15, 0.7, 200);
    polyfit_translation_x = polyfit(phi_nanless, psi_nanless(:, 1), fit_order);
    polyfit_translation_x_values = polyval(polyfit_translation_x, phi_values);
    polyfit_translation_y = polyfit(phi_nanless, psi_nanless(:, 2), fit_order);
    polyfit_translation_y_values = polyval(polyfit_translation_y, phi_values);
    polyfit_translation_z = polyfit(phi_nanless, psi_nanless(:, 3), fit_order);
    polyfit_translation_z_values = polyval(polyfit_translation_z, phi_values);
    polyfit_rotation_x = polyfit(phi_nanless, psi_nanless(:, 4), fit_order);
    polyfit_rotation_x_values = polyval(polyfit_rotation_x, phi_values);
    polyfit_rotation_y = polyfit(phi_nanless, psi_nanless(:, 5), fit_order);
    polyfit_rotation_y_values = polyval(polyfit_rotation_y, phi_values);
    polyfit_rotation_z = polyfit(phi_nanless, psi_nanless(:, 6), fit_order);
    polyfit_rotation_z_values = polyval(polyfit_rotation_z, phi_values);
    
    
    
    
    
    figure; axes; hold on; title('\phi');
    time = (1 : number_of_time_steps) * sampling_rate^(-1);
    xlabel('time'); ylabel('\phi');
    plot(time, phi, 'linewidth', 2);
    
    figure; axes; hold on; title('d\phi/dt');
    xlabel('time'); ylabel('d\phi/dt');
    plot(time, phi_dot, 'linewidth', 1);
    plot(time, phi_dot_relevant, 'linewidth', 2);
    
    
    
    % plot the stuff
    figure; axes; hold on; title('body velocity, translational');
    plot(ankle_x_elevation_trajectory, body_velocity_trajectory(:, 1:3));
    xlabel('\phi'); ylabel('body velocity');
    legend('x', 'y', 'z')
    
    figure; axes; hold on; title('body velocity, rotational');
    plot(ankle_x_elevation_trajectory, body_velocity_trajectory(:, 4:6));
    xlabel('\phi'); ylabel('body velocity');
    legend('x', 'y', 'z')
    
    % plot the stuff
    figure; axes; hold on; title('body velocity times (d\phi/dt)^{-1}, translational x');
    plot(phi_relevant, psi_relevant(:, 1), 'linewidth', 1);
    plot(phi_values, polyfit_translation_x_values, 'linewidth', 2);
    ylim([-0.25, 0.25]);
    figure; axes; hold on; title('body velocity times (d\phi/dt)^{-1}, translational y');
    plot(phi_relevant, psi_relevant(:, 2), 'linewidth', 1);
    plot(phi_values, polyfit_translation_y_values, 'linewidth', 2);
    ylim([-0.25, 0.25]);
    figure; axes; hold on; title('body velocity times (d\phi/dt)^{-1}, translational z');
    plot(phi_relevant, psi_relevant(:, 3), 'linewidth', 1);
    plot(phi_values, polyfit_translation_z_values, 'linewidth', 2);
    ylim([-0.25, 0.25]);
    
    figure; axes; hold on; title('body velocity times (d\phi/dt)^{-1}, rotational x');
    plot(phi_relevant, psi_relevant(:, 4), 'linewidth', 1);
    plot(phi_values, polyfit_rotation_x_values, 'linewidth', 2);
    ylim([-1.5, 1.5]);
    figure; axes; hold on; title('body velocity times (d\phi/dt)^{-1}, rotational y');
    plot(phi_relevant, psi_relevant(:, 5), 'linewidth', 1);
    plot(phi_values, polyfit_rotation_y_values, 'linewidth', 2);
    ylim([-1.5, 1.5]);
    figure; axes; hold on; title('body velocity times (d\phi/dt)^{-1}, rotational z');
    plot(phi_relevant, psi_relevant(:, 6), 'linewidth', 1);
    plot(phi_values, polyfit_rotation_z_values, 'linewidth', 2);
    ylim([-1.5, 1.5]);
    
    
    
    
    
    
    
    polyfit_azimuth = [];
    polyfit_elevation = [];
    return
    
    
    
    
    
    

    joint_center_velocity_trajectory_world = centdiff(joint_center_trajectory_world, sampling_rate^(-1));

    % transform velocity into ankle frame
    joint_center_velocity_trajectory_ankle = zeros(size(joint_center_velocity_trajectory_world));
    for i_time = 1:number_of_time_steps
        T_ankle_to_world_current = reshape(joint_frame_trajectory(i_time, :), 4, 4);
        if any(isnan(T_ankle_to_world_current))
            joint_center_velocity_trajectory_ankle(i_time, :) = NaN;
        else
            joint_center_velocity_world_current = [joint_center_velocity_trajectory_world(i_time, :)'; 0];
            joint_center_velocity_ankle_current = T_ankle_to_world_current^(-1) * joint_center_velocity_world_current;
            joint_center_velocity_trajectory_ankle(i_time, :) = joint_center_velocity_ankle_current(1:3);
        end
    end

    % exclude the velocity vectors with too small magnitude
%     speed_threshold = 0.03;
    speed_threshold = 0.06;

    % replace the time steps with too small speed with NaNs
    joint_center_speed_trajectory = sum(joint_center_velocity_trajectory_ankle.^2, 2).^0.5;
    joint_center_velocity_trajectory_fast = joint_center_velocity_trajectory_ankle;
    joint_center_velocity_trajectory_fast(joint_center_speed_trajectory < speed_threshold, :) = NaN;
    
    % norm and flip those that go in the wrong direction
    joint_center_velocity_trajectory_fast_nanless = joint_center_velocity_trajectory_fast(~any(isnan(joint_center_velocity_trajectory_fast), 2), :);
    [V, ~] = eig(joint_center_velocity_trajectory_fast_nanless'*joint_center_velocity_trajectory_fast_nanless);
    main_direction = V(:, end);
    joint_center_velocity_trajectory_fast_normed = joint_center_velocity_trajectory_fast .* repmat(joint_center_speed_trajectory, 1, 3).^(-1);
    joint_center_velocity_trajectory_fast_normed_flipped = joint_center_velocity_trajectory_fast_normed;
    
    joint_center_velocity_dot_main_direction = joint_center_velocity_trajectory_fast_normed_flipped * main_direction;
    joint_center_velocity_trajectory_fast_normed_flipped(joint_center_velocity_dot_main_direction<0, :) = -joint_center_velocity_trajectory_fast_normed_flipped(joint_center_velocity_dot_main_direction<0, :);

%     joint_center_velocity_trajectory_fast_normed_flipped(joint_center_velocity_trajectory_fast_normed(:, 1)<0, :) = -joint_center_velocity_trajectory_fast_normed_flipped(joint_center_velocity_trajectory_fast_normed(:, 1)<0, :);
%     joint_center_velocity_trajectory_fast_normed_flipped(ankle_elevation_velocity_trajectory<0, :) = -joint_center_velocity_trajectory_fast_normed_flipped(ankle_elevation_velocity_trajectory<0, :);

    % transform into azimuth and elevation
    [joint_center_velocity_azimuth_trajectory, joint_center_velocity_elevation_trajectory] = cart2sph(joint_center_velocity_trajectory_fast_normed_flipped(:, 1), joint_center_velocity_trajectory_fast_normed_flipped(:, 2), joint_center_velocity_trajectory_fast_normed_flipped(:, 3));

    % fit dependency of velocity direction on ankle elevation angle in spherical coordinates
    x = ankle_elevation_trajectory(~isnan(joint_center_velocity_azimuth_trajectory));
    y = joint_center_velocity_azimuth_trajectory(~isnan(joint_center_velocity_azimuth_trajectory));
    polyfit_azimuth = polyfit(ankle_elevation_trajectory(~isnan(joint_center_velocity_azimuth_trajectory)), joint_center_velocity_azimuth_trajectory(~isnan(joint_center_velocity_azimuth_trajectory)), 5);
    polyfit_elevation = polyfit(ankle_elevation_trajectory(~isnan(joint_center_velocity_elevation_trajectory)), joint_center_velocity_elevation_trajectory(~isnan(joint_center_velocity_elevation_trajectory)), 5);
    ankle_elevation_values = linspace(min(ankle_elevation_trajectory), max(ankle_elevation_trajectory), 100);
    polyfit_azimuth_values = polyval(polyfit_azimuth, ankle_elevation_values);
    polyfit_elevation_values = polyval(polyfit_elevation, ankle_elevation_values);






%     figure; axes; hold on; axis equal; title('velocity, world');
%     xlabel('x'); ylabel('y'), zlabel('z');
%     xdata = joint_center_velocity_trajectory_world(:, 1)';
%     ydata = joint_center_velocity_trajectory_world(:, 2)';
%     zdata = joint_center_velocity_trajectory_world(:, 3)';
%     colordata = ankle_elevation_trajectory';
%     surface([xdata; xdata],[ydata; ydata],[zdata; zdata],[colordata; colordata], 'facecol', 'no', 'edgecol', 'interp', 'linew', 2);
% 
%     figure; axes; hold on; axis equal; title('velocity, ankle');
%     xlabel('x'); ylabel('y'), zlabel('z');
%     xdata = joint_center_velocity_trajectory_ankle(:, 1)';
%     ydata = joint_center_velocity_trajectory_ankle(:, 2)';
%     zdata = joint_center_velocity_trajectory_ankle(:, 3)';
%     colordata = ankle_elevation_trajectory';
%     surface([xdata; xdata],[ydata; ydata],[zdata; zdata],[colordata; colordata], 'facecol', 'no', 'edgecol', 'interp', 'linew', 2);




    figure; axes; hold on; axis equal; title('velocity, ankle');
    xlabel('x'); ylabel('y'), zlabel('z');
    xdata = joint_center_velocity_trajectory_fast(:, 1)';
    ydata = joint_center_velocity_trajectory_fast(:, 2)';
    zdata = joint_center_velocity_trajectory_fast(:, 3)';
    colordata = ankle_elevation_trajectory';
    surface([xdata; xdata],[ydata; ydata],[zdata; zdata],[colordata; colordata], 'facecol', 'no', 'edgecol', 'interp', 'linew', 2);
    ellipsoid(0,0,0,speed_threshold,speed_threshold,speed_threshold,20)
    alpha(0.3)

    figure; axes; hold on; axis equal; title('velocity, ankle');
    xlabel('x'); ylabel('y'), zlabel('z');
    xdata = joint_center_velocity_trajectory_fast_normed_flipped(:, 1)';
    ydata = joint_center_velocity_trajectory_fast_normed_flipped(:, 2)';
    zdata = joint_center_velocity_trajectory_fast_normed_flipped(:, 3)';
    colordata = ankle_elevation_trajectory';
    surface([xdata; xdata],[ydata; ydata],[zdata; zdata],[colordata; colordata], 'facecol', 'no', 'edgecol', 'interp', 'linew', 2);
    % ellipsoid(0,0,0,1,1,1,20)
    % alpha(0.3)


    figure; axes; hold on; xlabel('phi'); ylabel('velocity, azimuth')
    plot(ankle_elevation_trajectory, joint_center_velocity_azimuth_trajectory);
    plot(ankle_elevation_values, polyfit_azimuth_values);
    figure; axes; hold on; xlabel('phi'); ylabel('velocity, elevation')
    plot(ankle_elevation_trajectory, joint_center_velocity_elevation_trajectory);
    plot(ankle_elevation_values, polyfit_elevation_values);    
    
    
    
    
    
    
    
% end
