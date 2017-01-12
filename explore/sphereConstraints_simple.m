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

% ellipsoidConstraintExploration

% define plant
r_x = 1;
r_y = 1.0;
r_z = 1.0;

generate_trajectories               = 0;
approximate_constraints             = 0;
visualize_constraints               = 0;
find_rho_constraint_by_phi          = 0;
find_phi_constraint                 = 1;
calculate_allowed_body_velocities   = 1;
show_stick_figure                   = 0;

% movement_mode = 'position';
movement_mode = 'velocity';
% movement_mode = 'axis';

number_of_sets = 100;

plant = GeneralKinematicTree ...
( ...
  {[0; 0; r_z]; [0; 0; r_z]; [0; 0; r_z]; [0; 0; r_z]; [0; 0; r_z]; [0; 0; r_z]}, ...
  {[1; 0; 0], [0; 1; 0], [0; 0; 1], [0; 0; 1], [1; 0; 0], [0; 1; 0]}, ...
  [2 2 2 1 1 1], ...
  [1 1 1 1 1 1], ...
  {[eye(3) [0; 0; r_z]; [0 0 0 1]]}, ...
  {[0; 0; r_z]; [0; 0; r_z]; [0; 0; r_z]; [0; 0; r_z]; [0; 0; r_z]; [0; 0; r_z]}, ...
  {eye(3), eye(3), eye(3), eye(3), eye(3), eye(3)}, ...
  {zeros(6, 6), zeros(6, 6), zeros(6, 6), zeros(6, 6), zeros(6, 6), eye(6)} ...
);
plant.updateInternals();
number_of_joints = plant.numberOfJoints;





%% generate trajectories
if generate_trajectories
    total_time = 1;
    time_step = 0.01;
    time = (time_step : time_step : total_time)';
    number_of_time_steps_set = length(time);

    starting_joint_angles = zeros(1, number_of_joints);
    % starting_joint_angles(4) = pi/4;
    % starting_joint_angles = rand(number_of_joints, 1);

    amplitudes = [0.2 -0.1 -0.15 -0.15 0.2 -0.1] * 3;
    phase_offsets = [0.1 -0.3 0.15 0.35 0.4 -0.1];
    frequencies = [.5 1 1.5 1.3 1.8 1.1] * 0.3;
    starting_angle_set = zeros(1, 6);

    % meandering
%     amplitudes = [1 -0.7 0.5 2 1 0] * 1;
%     phase_offsets = [-0.3 0.1 0.15 0.35 0.4 -0.1];
%     frequencies = [0.02 .3 1.5 1.3 1.8 1.1] * 1;
%     starting_angle_set = zeros(5, 6);
%     starting_angle_set(2, 4) = 0.3;
%     starting_angle_set(3, 4) = 0.6;
%     starting_angle_set(4, 4) = 0.9;
%     starting_angle_set(5, 4) = 1.2;
%     number_of_sets = size(starting_angle_set, 1);
    
    % foot-like
    amplitudes = [0 0 0 0 1 5] * 1;
    phase_offsets = [0 0 0 0 0 0.5];
    frequencies = [0 0 0 0 0.25 0.25] * 1;
    starting_angle_set = zeros(number_of_sets, 6);
    starting_angle_set(:, 5) = 0.5;
    inversion = linspace(-1, 1, number_of_sets);
    starting_angle_set(:, 4) = rand(number_of_sets, 1)*2-1;
    sinusoids = (-cos((repmat(time, 1, number_of_joints)+repmat(phase_offsets, number_of_time_steps_set, 1))*2*pi.*repmat(frequencies, number_of_time_steps_set, 1))) .* repmat(amplitudes, number_of_time_steps_set, 1);
    velocity_factor = rand(1, number_of_sets) + 0.5;

    if strcmp(movement_mode, 'position')
        joint_angle_trajectories = repmat(starting_joint_angles, number_of_time_steps_set, 1) + sinusoids;
    elseif strcmp(movement_mode, 'velocity')
        joint_angle_trajectories = zeros(number_of_time_steps_set * number_of_sets, number_of_joints);
        transformation_trajectories = zeros(number_of_time_steps_set * number_of_sets, 16);
        body_velocity_trajectory = zeros(number_of_time_steps_set * number_of_sets, 6);
        for i_set = 1 : number_of_sets
            amplitudes = [0 0 0 0 3 5*inversion(i_set)] * 1;
            sinusoids = (-cos((repmat(time, 1, number_of_joints)+repmat(phase_offsets, number_of_time_steps_set, 1))*2*pi.*repmat(frequencies, number_of_time_steps_set, 1))) .* repmat(amplitudes, number_of_time_steps_set, 1);
            starting_joint_angles = starting_angle_set(i_set, :)';

            joint_angle_trajectories_set = zeros(number_of_time_steps_set, number_of_joints);
            joint_angle_trajectories_set(1, :) = starting_joint_angles;
            transformation_trajectories_set = zeros(number_of_time_steps_set, 16);
            body_velocity_trajectory_set = zeros(number_of_time_steps_set, 6);
            
            joint_velocity_trajectories_unconstrained = sinusoids * velocity_factor(i_set);
            joint_velocity_trajectories_constrained = zeros(number_of_time_steps_set, number_of_joints);
            
            plant.jointAngles = starting_joint_angles;
            plant.updateConfiguration;
            transformation_trajectories_set(1, :) = reshape(plant.endEffectorTransformations{1}, 1, 16);
            for i_time = 2 : number_of_time_steps_set
                p_contact = [plant.jointAngles(1:2); 0; 1];
                J_contact = plant.calculateArbitraryPointJacobian(p_contact, plant.numberOfJoints, 'world');
                A = [J_contact; [0 0 0 1 0 0];];
                P_A = A' * (A * A')^(-1) * A; % projection onto the space spanned by A
                P_A_orth = eye(number_of_joints) - P_A;


                joint_velocity_unconstrained = joint_velocity_trajectories_unconstrained(i_time, :);
                joint_velocity_constrained = (P_A_orth * joint_velocity_unconstrained')';

            %         joint_velocity_constrained = joint_velocity_unconstrained;

                % euler step
                joint_angle_trajectories_set(i_time, :) = joint_angle_trajectories_set(i_time-1, :) + time_step * joint_velocity_constrained;

                plant.jointAngles = joint_angle_trajectories_set(i_time, :)';
                plant.updateKinematics;

                transformation_trajectories_set(i_time, :) = reshape(plant.endEffectorTransformations{1}, 1, 16);
                body_velocity_trajectory_set(i_time, :) = plant.bodyJacobians{1} * joint_velocity_constrained';
            end
            joint_angle_trajectories((i_set-1)*number_of_time_steps_set+1 : i_set*number_of_time_steps_set, :) = joint_angle_trajectories_set;
            transformation_trajectories((i_set-1)*number_of_time_steps_set+1 : i_set*number_of_time_steps_set, :) = transformation_trajectories_set;
            body_velocity_trajectory((i_set-1)*number_of_time_steps_set+1 : i_set*number_of_time_steps_set, :) = body_velocity_trajectory_set;
        end
    elseif strcmp(movement_mode, 'axis')
        joint_angle_trajectories = zeros(number_of_time_steps_set * number_of_sets, number_of_joints);
        transformation_trajectories = zeros(number_of_time_steps_set * number_of_sets, 16);
        tic
        for i_set = 1 : number_of_sets
            starting_joint_angles = starting_angle_set(i_set, :);

            joint_angle_trajectories_set = zeros(number_of_time_steps_set, number_of_joints);
            joint_angle_trajectories_set(1, :) = starting_joint_angles;
            axis_of_rotation_trajectory = [sinusoids(:, 1:2) zeros(number_of_time_steps_set, 1)] ./ repmat(sum([sinusoids(:, 1:2) zeros(number_of_time_steps_set, 1)].^2, 2).^(0.5), 1, 3);
            angular_velocity_trajectory = sinusoids(:, 4);

        %     axis_of_rotation_trajectory = repmat([sqrt(2)/2 sqrt(2)/2 0], number_of_time_steps, 1);
            angular_velocity_trajectory = ones(number_of_time_steps_set, 1) * 0.1;

            joint_velocity_trajectories = zeros(number_of_time_steps_set, number_of_joints);
            transformation_trajectories_set = zeros(number_of_time_steps_set, 16);

            plant.jointAngles = starting_joint_angles;
            plant.updateConfiguration;
            transformation_trajectories_set(1, :) = reshape(plant.endEffectorTransformations{1}, 1, 16);
            for i_time = 2 : number_of_time_steps_set
                omega = axis_of_rotation_trajectory(i_time, :)';
                R_current = plant.endEffectorTransformations{1}(1:3, 1:3);
                R_euler_step = expAxis(omega, angular_velocity_trajectory(i_time) * time_step);
                R_new = R_euler_step * R_current;

                theta_4to6_current = rotm2eul(R_current, 'ZYX');
                theta_4to6_new = rotm2eul(R_new, 'ZYX');

                % constraints
                center_movement_direction = cross(omega, [0; 0; 1]);
                center_velocity = center_movement_direction * angular_velocity_trajectory(i_time);
                theta_1to3_current = plant.jointAngles(1:3)';
                theta_1to3_new = theta_1to3_current + time_step * center_velocity;

                joint_angle_trajectories_set(i_time, :) = [theta_1to3_new; theta_4to6_new'];
                plant.jointAngles = joint_angle_trajectories_set(i_time, :);
                plant.updateKinematics;

                transformation_trajectories_set(i_time, :) = reshape(plant.endEffectorTransformations{1}, 1, 16);
            end
            joint_angle_trajectories((i_set-1)*number_of_time_steps_set+1 : i_set*number_of_time_steps_set, :) = joint_angle_trajectories_set;
            transformation_trajectories((i_set-1)*number_of_time_steps_set+1 : i_set*number_of_time_steps_set, :) = transformation_trajectories_set;

            disp(['finished set ' num2str(i_set)]);
        end
        toc
    end
    plot(joint_angle_trajectories)
    legend('show')
end

%% approximate constraints
if approximate_constraints
    % calculate trajectories for points of interest and elevation angles
    number_of_time_steps = size(transformation_trajectories, 1);
    sphere_x_azimuth_trajectory = zeros(number_of_time_steps, 1);
    phi_trajectory = zeros(number_of_time_steps, 1);
    sphere_y_azimuth_trajectory = zeros(number_of_time_steps, 1);
    rho_trajectory = zeros(number_of_time_steps, 1);
    gamma_trajectory = zeros(number_of_time_steps, 1);
    
    points_of_interest_sphere = [eye(3); ones(1, 3)]; % coordinates of the points of interest in sphere coordinates
    points_of_interest_world_trajectories = zeros(number_of_time_steps, size(points_of_interest_sphere, 2), 3);
    points_of_interest_world_trajectories_x = zeros(number_of_time_steps, size(points_of_interest_sphere, 2));
    points_of_interest_world_trajectories_y = zeros(number_of_time_steps, size(points_of_interest_sphere, 2));
    points_of_interest_world_trajectories_z = zeros(number_of_time_steps, size(points_of_interest_sphere, 2));
    for i_time = 1 : number_of_time_steps
        % get transformation
        sphere_transformation = reshape(transformation_trajectories(i_time, :), 4, 4);
        sphere_rotation = sphere_transformation(1:3, 1:3);
%         euler_angles = rotm2eul(sphere_rotation, 'ZYX');
        euler_angles = eulerAnglesFromRotationMatrixZXY(sphere_rotation);
        gamma_trajectory(i_time) = euler_angles(1);
        phi_trajectory(i_time) = euler_angles(2);
        rho_trajectory(i_time) = euler_angles(3);
        
        % calculate points of interest
        points_of_interest_world_current = sphere_transformation * points_of_interest_sphere;
        points_of_interest_world_trajectories(i_time, :, :) = points_of_interest_world_current(1:3, :);
        
        points_of_interest_world_trajectories_x(i_time, :) = points_of_interest_world_current(1, 1:3);
        points_of_interest_world_trajectories_y(i_time, :) = points_of_interest_world_current(2, 1:3);
        points_of_interest_world_trajectories_z(i_time, :) = points_of_interest_world_current(3, 1:3);
        
    end
    
    phi_dot_trajectory = deriveByTime(phi_trajectory, time_step);
    rho_dot_trajectory = deriveByTime(rho_trajectory, time_step);
    
    d_g_by_d_phi_vee_trajectory = zeros(number_of_time_steps, 6);
    d_g_by_d_rho_vee_trajectory = zeros(number_of_time_steps, 6);
    for i_time = 1:number_of_time_steps
        sphere_transformation = reshape(transformation_trajectories(i_time, :), 4, 4);
        body_velocity = body_velocity_trajectory(i_time, :);
        body_velocity_wedge = wedgeTwist(body_velocity);

        d_g_by_d_phi = sphere_transformation * body_velocity_wedge * phi_dot_trajectory(i_time)^(-1);
        d_g_by_d_phi_vee = veeTwist(d_g_by_d_phi);
        d_g_by_d_phi_vee_trajectory(i_time, :) = d_g_by_d_phi_vee;
        
        d_g_by_d_rho = sphere_transformation * body_velocity_wedge * rho_dot_trajectory(i_time)^(-1);
        d_g_by_d_rho_vee = veeTwist(d_g_by_d_rho);
        d_g_by_d_rho_vee_trajectory(i_time, :) = d_g_by_d_rho_vee;
    end    
    
    % extract relevant data points
    relevant_data_points = abs(phi_trajectory) > 0.01 & abs(phi_trajectory - 1) > 0.01;
    irrelevant_data_points = ~relevant_data_points;
    phi_relevant_trajectory = phi_trajectory;   phi_relevant_trajectory(irrelevant_data_points) = NaN;
    phi_dot_trajectory_relevant = phi_dot_trajectory;      phi_dot_trajectory_relevant(irrelevant_data_points) = NaN;
    rho_relevant_trajectory = rho_trajectory;              rho_relevant_trajectory(irrelevant_data_points) = NaN;
    body_velocity_trajectory_relevant = body_velocity_trajectory;       body_velocity_trajectory_relevant(irrelevant_data_points, :) = NaN;
    d_g_by_d_phi_vee_relevant = d_g_by_d_phi_vee_trajectory;            d_g_by_d_phi_vee_relevant(irrelevant_data_points, :) = NaN;
    d_g_by_d_rho_vee_relevant = d_g_by_d_rho_vee_trajectory;            d_g_by_d_rho_vee_relevant(irrelevant_data_points, :) = NaN;
    
    
    
    
    
    
    
    
    data_points_for_fit = 1 : 1 : number_of_time_steps;
    phi_fitpoints = phi_trajectory(data_points_for_fit);
    rho_fitpoints = rho_trajectory(data_points_for_fit);
    p_1_x = points_of_interest_world_trajectories_x(data_points_for_fit, 1);
    p_1_y = points_of_interest_world_trajectories_y(data_points_for_fit, 1);
    p_1_z = points_of_interest_world_trajectories_z(data_points_for_fit, 1);
    p_2_x = points_of_interest_world_trajectories_x(data_points_for_fit, 2);
    p_2_y = points_of_interest_world_trajectories_y(data_points_for_fit, 2);
    p_2_z = points_of_interest_world_trajectories_z(data_points_for_fit, 2);
    p_3_x = points_of_interest_world_trajectories_x(data_points_for_fit, 3);
    p_3_y = points_of_interest_world_trajectories_y(data_points_for_fit, 3);
    p_3_z = points_of_interest_world_trajectories_z(data_points_for_fit, 3);

    p_1_z_fit = fit([phi_fitpoints, rho_fitpoints], p_1_z, 'poly55');
    p_2_z_fit = fit([phi_fitpoints, rho_fitpoints], p_2_z, 'poly55');
    p_3_z_fit = fit([phi_fitpoints, rho_fitpoints], p_3_z, 'poly55');
    
%     figure; plot(p_1_z_fit, [phi_fitpoints, rho_fitpoints], p_1_z);
%     figure; plot(p_2_z_fit, [phi_fitpoints, rho_fitpoints], p_2_z);
%     figure; plot(p_3_z_fit, [phi_fitpoints, rho_fitpoints], p_3_z);
    
%     figure; axes; hold on; title('parameter trajectories')
%     plot(phi_trajectory);
%     plot(rho_trajectory);
%     plot(gamma_trajectory);
%     legend('phi', 'rho', 'gamma')
%     
%     figure; axes; hold on; title('joint angle trajectories')
%     plot(joint_angle_trajectories(:, 4:6));
%     legend('\theta_4', '\theta_5', '\theta_6')
    
%     figure; axes; hold on; title('parameter trajectories')
%     plot(phi_trajectory, 'linewidth', 2);
%     plot(rho_trajectory, 'linewidth', 2);
%     plot(gamma_trajectory, 'linewidth', 2);
%     plot(joint_angle_trajectories(:, 4:6));
%     legend('\phi', '\rho', '\gamma', '\theta_4', '\theta_5', '\theta_6')
    
%     figure; axes; hold on; title('points of interest: x-paths')
%     plot3(phi_trajectory, rho_trajectory, points_of_interest_world_trajectories_x);
%     xlabel('\phi'); ylabel('\rho'); zlabel('x-coordinate')
%     
%     figure; axes; hold on; title('points of interest: y-paths')
%     plot3(phi_trajectory, rho_trajectory, points_of_interest_world_trajectories_y);
%     xlabel('\phi'); ylabel('\rho'); zlabel('x-coordinate')
%     
%     figure; axes; hold on; title('points of interest: z-paths')
%     plot3(phi_trajectory, rho_trajectory, points_of_interest_world_trajectories_z);
%     xlabel('\phi'); ylabel('\rho'); zlabel('x-coordinate')
    
%     figure; axes; hold on
%     plot(phi_trajectory)
%     plot(phi_relevant_trajectory)




%     % visualize - body velocity vs. phi and rho
%     figure; axes; hold on; title('v_1');
%     plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 1));
%     plot3(phi_relevant_trajectory, rho_relevant_trajectory, body_velocity_trajectory_relevant(:, 1));
%     xlabel('\phi'); ylabel('\rho'); zlabel('v_1')
%     figure; axes; hold on; title('v_2');
%     plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 2));
%     plot3(phi_relevant_trajectory, rho_relevant_trajectory, body_velocity_trajectory_relevant(:, 2));
%     xlabel('\phi'); ylabel('\rho'); zlabel('v_2')
%     figure; axes; hold on; title('v_3');
%     plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 3));
%     plot3(phi_relevant_trajectory, rho_relevant_trajectory, body_velocity_trajectory_relevant(:, 3));
%     xlabel('\phi'); ylabel('\rho'); zlabel('v_3')    
%     figure; axes; hold on; title('v_4');
%     plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 4));
%     plot3(phi_relevant_trajectory, rho_relevant_trajectory, body_velocity_trajectory_relevant(:, 4));
%     figure; axes; hold on; title('v_5');
%     plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 5));
%     plot3(phi_relevant_trajectory, rho_relevant_trajectory, body_velocity_trajectory_relevant(:, 5));
%     figure; axes; hold on; title('v_6');
%     plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 6));
%     plot3(phi_relevant_trajectory, rho_relevant_trajectory, body_velocity_trajectory_relevant(:, 6));
    
    % visualize - normalized body velocity vs. phi and rho
%     body_velocity_trajectory_norm = sum(body_velocity_trajectory.^2, 2).^(0.5);
%     body_velocity_trajectory_normed = body_velocity_trajectory .* repmat(body_velocity_trajectory_norm.^(-1), 1, 6);
%     body_velocity_trajectory_normed_relevant = body_velocity_trajectory_normed;       body_velocity_trajectory_normed_relevant(irrelevant_data_points, :) = NaN;
%     figure; axes; hold on; title('v_1');
%     plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory_normed(:, 1));
%     plot3(phi_relevant_trajectory, rho_relevant_trajectory, body_velocity_trajectory_normed_relevant(:, 1));
%     xlabel('\phi'); ylabel('\rho'); zlabel('v_1')
%     figure; axes; hold on; title('v_2');
%     plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory_normed(:, 2));
%     plot3(phi_relevant_trajectory, rho_relevant_trajectory, body_velocity_trajectory_normed_relevant(:, 2));
%     xlabel('\phi'); ylabel('\rho'); zlabel('v_2')
%     figure; axes; hold on; title('v_3');
%     plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory_normed(:, 3));
%     plot3(phi_relevant_trajectory, rho_relevant_trajectory, body_velocity_trajectory_normed_relevant(:, 3));
%     xlabel('\phi'); ylabel('\rho'); zlabel('v_3')    
%     figure; axes; hold on; title('v_4');
%     plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory_normed(:, 4));
%     plot3(phi_relevant_trajectory, rho_relevant_trajectory, body_velocity_trajectory_normed_relevant(:, 4));
%     figure; axes; hold on; title('v_5');
%     plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory_normed(:, 5));
%     plot3(phi_relevant_trajectory, rho_relevant_trajectory, body_velocity_trajectory_normed_relevant(:, 5));
%     figure; axes; hold on; title('v_6');
%     plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory_normed(:, 6));
%     plot3(phi_relevant_trajectory, rho_relevant_trajectory, body_velocity_trajectory_normed_relevant(:, 6));
end

%% visualize_constraints
if visualize_constraints
    % set configuration and update
    time_index = 5000;
    theta = joint_angle_trajectories(time_index, :)';
    plant.jointAngles = theta;
    plant.updateConfiguration;
    
    % calculate phi, rho and z-values
    sphere_transformation = plant.endEffectorTransformations{1};

    sphere_frame_x_axis_current = sphere_transformation(1:3, 1);
    [~, phi] = cart2sph(sphere_frame_x_axis_current(1), sphere_frame_x_axis_current(2), sphere_frame_x_axis_current(3));
    sphere_frame_y_axis_current = sphere_transformation(1:3, 2);
    [~, rho] = cart2sph(sphere_frame_y_axis_current(1), sphere_frame_y_axis_current(2), sphere_frame_y_axis_current(3));
    
    cos_gamma = dot(sphere_frame_x_axis_current, [1; 0; 0]);
    gamma = acos(cos_gamma);


    points_of_interest_world_current = sphere_transformation * points_of_interest_sphere;
    points_of_interest_world = points_of_interest_world_current(1:3, :);

    points_of_interest_world_x = points_of_interest_world_current(1, 1:3);
    points_of_interest_world_y = points_of_interest_world_current(2, 1:3);
    points_of_interest_world_z = points_of_interest_world_current(3, 1:3);
    
    % calculate derivatives
    [d_z_1_by_d_phi, d_z_1_by_d_rho] = differentiate(p_1_z_fit, phi, rho);
    [d_z_2_by_d_phi, d_z_2_by_d_rho] = differentiate(p_2_z_fit, phi, rho);
    [d_z_3_by_d_phi, d_z_3_by_d_rho] = differentiate(p_3_z_fit, phi, rho);
    
    % create constraint matrices for z_i
    p_1_local = [1; 0; 0; 1];
    p_2_local = [0; 1; 0; 1];
    p_3_local = [0; 0; 1; 1];
    J_p_1 = plant.calculateArbitraryPointJacobian(p_1_local, number_of_joints, 'local');
    J_p_2 = plant.calculateArbitraryPointJacobian(p_2_local, number_of_joints, 'local');
    J_p_3 = plant.calculateArbitraryPointJacobian(p_3_local, number_of_joints, 'local');
    d_z_1_by_d_x = J_p_1(3, 1); d_z_1_by_d_y = J_p_1(3, 2); d_z_1_by_d_z = J_p_1(3, 3);
    d_z_2_by_d_x = J_p_2(3, 1); d_z_2_by_d_y = J_p_2(3, 2); d_z_2_by_d_z = J_p_2(3, 3);
    d_z_3_by_d_x = J_p_3(3, 1); d_z_3_by_d_y = J_p_3(3, 2); d_z_3_by_d_z = J_p_3(3, 3);
    
    
    
    
    
    
    
    
    % visualize
    phi_values = -1 : 0.1 : 1;
    rho_values = -1 : 0.1 : 1;
    [phi_grid, rho_grid] = meshgrid(phi_values, rho_values);
    p_1_z_grid = feval(p_1_z_fit, phi_grid, rho_grid);
    p_2_z_grid = feval(p_2_z_fit, phi_grid, rho_grid);
    p_3_z_grid = feval(p_3_z_fit, phi_grid, rho_grid);
    
    figure; axes; hold on
    surf(phi_grid, rho_grid, p_1_z_grid);
    plot3(phi_trajectory, rho_trajectory, points_of_interest_world_trajectories_z);
    plot3(phi, rho, points_of_interest_world_z(1), 'o', 'markersize', 10, 'linewidth', 3);
    plot3(phi + [-1 1], rho + [0 0], points_of_interest_world_z(1) + [-1 1] * d_z_1_by_d_phi, '-', 'linewidth', 3);
    plot3(phi + [0 0], rho + [-1 1], points_of_interest_world_z(1) + [-1 1] * d_z_1_by_d_rho, '-', 'linewidth', 3);

%     figure; axes; hold on
%     surf(phi_grid, rho_grid, p_2_z_grid);
%     plot3(phi_trajectory, rho_trajectory, points_of_interest_world_trajectories_z);
%     plot3(phi, rho, points_of_interest_world_z(2), 'o', 'markersize', 10, 'linewidth', 3);
%     plot3(phi + [-1 1], rho + [0 0], points_of_interest_world_z(2) + [-1 1] * d_z_2_by_d_phi, '-', 'linewidth', 3);
%     plot3(phi + [0 0], rho + [-1 1], points_of_interest_world_z(2) + [-1 1] * d_z_2_by_d_rho, '-', 'linewidth', 3);
%     
%     figure; axes; hold on
%     surf(phi_grid, rho_grid, p_3_z_grid);
%     plot3(phi_trajectory, rho_trajectory, points_of_interest_world_trajectories_z);
%     plot3(phi, rho, points_of_interest_world_z(3), 'o', 'markersize', 10, 'linewidth', 3);
%     plot3(phi + [-1 1], rho + [0 0], points_of_interest_world_z(3) + [-1 1] * d_z_3_by_d_phi, '-', 'linewidth', 3);
%     plot3(phi + [0 0], rho + [-1 1], points_of_interest_world_z(3) + [-1 1] * d_z_3_by_d_rho, '-', 'linewidth', 3);
    
    
end

%% find rho constraint by phi
if find_rho_constraint_by_phi
    % calculate z-trajectories for multiple points along the sphere z-axis
    number_of_time_steps = size(transformation_trajectories, 1);
    sphere_x_azimuth_trajectory = zeros(number_of_time_steps, 1);
    phi_trajectory = zeros(number_of_time_steps, 1);
    sphere_y_azimuth_trajectory = zeros(number_of_time_steps, 1);
    rho_trajectory = zeros(number_of_time_steps, 1);
    gamma_trajectory = zeros(number_of_time_steps, 1);
    
    points_of_interest_z_local = -1 : 0.25 : 1;
    number_of_points_of_interest = length(points_of_interest_z_local);
    points_of_interest_sphere = [zeros(2, number_of_points_of_interest); points_of_interest_z_local; ones(1, number_of_points_of_interest)]; % coordinates of the points of interest in sphere coordinates
    points_of_interest_world_trajectories_x = zeros(number_of_time_steps, size(points_of_interest_sphere, 2));
    points_of_interest_world_trajectories_y = zeros(number_of_time_steps, size(points_of_interest_sphere, 2));
    points_of_interest_world_trajectories_z = zeros(number_of_time_steps, size(points_of_interest_sphere, 2));
    for i_time = 1 : number_of_time_steps
        sphere_transformation = reshape(transformation_trajectories(i_time, :), 4, 4);
        sphere_rotation = sphere_transformation(1:3, 1:3);
%         euler_angles = rotm2eul(sphere_rotation, 'ZYX');
        euler_angles = eulerAnglesFromRotationMatrixZXY(sphere_rotation);
        gamma_trajectory(i_time) = euler_angles(1);
        phi_trajectory(i_time) = euler_angles(2);
        rho_trajectory(i_time) = euler_angles(3);
        
        points_of_interest_world_current = sphere_transformation * points_of_interest_sphere;
        
        points_of_interest_world_trajectories_x(i_time, :) = points_of_interest_world_current(1, :);
        points_of_interest_world_trajectories_y(i_time, :) = points_of_interest_world_current(2, :);
        points_of_interest_world_trajectories_z(i_time, :) = points_of_interest_world_current(3, :);
    end
    
    % fit surfaces to phi-rho-z data
    data_points_for_fit = 1 : 1 : number_of_time_steps;
    phi_fitpoints = phi_trajectory(data_points_for_fit);
    rho_fitpoints = rho_trajectory(data_points_for_fit);
    poi_z_fits = cell(1, number_of_points_of_interest);
    for i_poi = 1 : number_of_points_of_interest
        p_z = points_of_interest_world_trajectories_z(data_points_for_fit, i_poi);
        poi_z_fits{i_poi} = fit([phi_fitpoints, rho_fitpoints], p_z, 'poly55');
%         figure; plot(poi_z_fits{i_poi}, [phi_fitpoints, rho_fitpoints], points_of_interest_world_trajectories_z(data_points_for_fit, i_poi));
%         xlabel('\phi'); ylabel('\rho'); zlabel('z'); 
    end
    % 
    phi_values = -0.3 : 0.1 : 0.3;
    number_of_phi_values = length(phi_values);
    rho_values = -0.3 : 0.01 : 0.3;
    
    z_values = cell(number_of_phi_values, number_of_points_of_interest);
    poi_z_fits_by_rho = cell(number_of_phi_values, number_of_points_of_interest);
    poi_z_fits_by_rho_error = zeros(number_of_phi_values, number_of_points_of_interest);
    poi_z_fits_by_rho_derivative = zeros(number_of_phi_values, number_of_points_of_interest);
    poi_z_fits_by_rho_derivative_error = zeros(number_of_phi_values, number_of_points_of_interest);
    for i_phi = 1 : length(phi_values)
        % for each point of interest, calculate derivative at 0 and use as fit
        
        
%         figure; hold on
        for i_poi = 1 : number_of_points_of_interest
            z_values{i_phi, i_poi} = feval(poi_z_fits{i_poi}, [ones(length(rho_values), 1)*phi_values(i_phi), rho_values']);
            poi_z_fits_by_rho{i_phi, i_poi} = fit(rho_values', z_values{i_phi, i_poi}, 'poly1');
            error = feval(poi_z_fits_by_rho{i_phi, i_poi}, rho_values) - z_values{i_phi, i_poi};
            poi_z_fits_by_rho_error(i_phi, i_poi) = sum(error.^2)^(0.5);
            
            
            [d_z_by_d_phi, poi_z_fits_by_rho_derivative(i_phi, i_poi)] = differentiate(poi_z_fits{i_poi}, phi_values(i_phi), 0);
            poi_z_fits_by_rho_derivative(i_phi, i_poi) = 0; % for some reason the differentiate function returns NaN, but it should be 0
            first_order_fit_values = feval(poi_z_fits{i_poi}, [phi_values(i_phi), 0]) + rho_values * poi_z_fits_by_rho_derivative(i_phi, i_poi);
            error = first_order_fit_values - z_values{i_phi, i_poi}';
            poi_z_fits_by_rho_derivative_error(i_phi, i_poi) = sum(error.^2)^(0.5);
            
%             plot(rho_values, z_values{i_phi, i_poi});
%             plot(rho_values, first_order_fit_values);
            
%             figure; plot(poi_z_fits{i_poi}, [phi_fitpoints, rho_fitpoints], points_of_interest_world_trajectories_z(data_points_for_fit, i_poi));
%             xlabel('\phi'); ylabel('\rho'); zlabel('z'); 
        end
        
        
        
%         % plot z-curves
%         figure; axes; hold on;
%         for i_poi = 1 : number_of_points_of_interest
%             plot(rho_values, z_values{i_phi, i_poi});
%             plot(rho_values, feval(poi_z_fits_by_rho{i_phi, i_poi}, rho_values));
%         end
    end
    
    figure; plot(poi_z_fits_by_rho_derivative_error', '-x'); title('different colors are different values of \phi')
    xlabel('point of interest'); ylabel('error')
    
    
    [~, poi_z_fits_by_rho_derivative_error_minimum_indices] = min(poi_z_fits_by_rho_derivative_error, [], 2);
    center_of_rho_rotation_z_estimate = mean(points_of_interest_sphere(3, poi_z_fits_by_rho_derivative_error_minimum_indices));
    rho_rotation_center_local = [0; 0; center_of_rho_rotation_z_estimate; 1];
    
    
    
    
    
    
    
    
    
    
    
    % visualize
    phi_values = -1 : 0.1 : 0;
    rho_values = -0.3 : 0.1 : 0.3;
    [phi_grid, rho_grid] = meshgrid(phi_values, rho_values);
    poi_z_grids = cell(1, number_of_points_of_interest);
    for i_poi = 1 : number_of_points_of_interest
        poi_z_grids{i_poi} = feval(poi_z_fits{i_poi}, phi_grid, rho_grid);
    end
    pois_to_plot = 1 : number_of_points_of_interest;
    figure; axes; hold on;
    xlabel('phi'); ylabel('rho'); 
    for i_poi = pois_to_plot
        surf(phi_grid, rho_grid, poi_z_grids{i_poi});
        plot3(phi_trajectory, rho_trajectory, points_of_interest_world_trajectories_z(:, i_poi));
    end
    
    % fix a phi value
    phi = -0.5;
    
    % evaluate changes of z_i by rho for this phi value
    rho_values = (-0.3 : 0.01 : 0.3)';
    phi_values = phi * ones(size(rho_values));
    z_values = cell(1, number_of_points_of_interest);
    for i_poi = 1 : number_of_points_of_interest
        z_values{i_poi} = feval(poi_z_fits{i_poi}, [phi_values, rho_values]);
    end
    
    figure; axes; hold on; axis equal
    for i_poi = 1 : number_of_points_of_interest
        plot(rho_values, z_values{i_poi});
    end
end

%% find_phi_constraint
if find_phi_constraint
    rho_rotation_center_local = [0; 0; 0; 1]; % eventually this might depend upon phi, but for now, we can assume it's fixed
    
    number_of_time_steps = size(transformation_trajectories, 1);
    sphere_x_azimuth_trajectory = zeros(number_of_time_steps, 1);
    phi_trajectory = zeros(number_of_time_steps, 1);
    sphere_y_azimuth_trajectory = zeros(number_of_time_steps, 1);
    rho_trajectory = zeros(number_of_time_steps, 1);
    gamma_trajectory = zeros(number_of_time_steps, 1);
    
    for i_time = 1 : number_of_time_steps
        sphere_transformation = reshape(transformation_trajectories(i_time, :), 4, 4);
        sphere_rotation = sphere_transformation(1:3, 1:3);
%         euler_angles = rotm2eul(sphere_rotation, 'ZYX');
        euler_angles = eulerAnglesFromRotationMatrixZXY(sphere_rotation);
        gamma_trajectory(i_time) = euler_angles(1);
        phi_trajectory(i_time) = euler_angles(2);
        rho_trajectory(i_time) = euler_angles(3);
    end    
    
    irrelevant_data_points = number_of_time_steps_set : number_of_time_steps_set : number_of_time_steps;
    phi_relevant_trajectory = phi_trajectory;   phi_relevant_trajectory(irrelevant_data_points) = NaN;
    phi_dot_trajectory_relevant = deriveByTime(phi_relevant_trajectory, time_step);
    rho_relevant_trajectory = rho_trajectory;   rho_relevant_trajectory(irrelevant_data_points) = NaN;
    rho_dot_trajectory_relevant = deriveByTime(rho_relevant_trajectory, time_step);
    
    
    phi_dot_trajectory = deriveByTime(phi_trajectory, time_step);
    rho_dot_trajectory = deriveByTime(rho_trajectory, time_step);
    
    
%     figure; axes; hold on
%     plot(phi_dot_trajectory_relevant, 'linewidth', 2);
%     plot(rho_dot_trajectory_relevant, 'linewidth', 2);
%     plot(body_velocity_trajectory);
%     legend('\phi', '\rho', 'v_1', 'v_2', 'v_3', 'v_4', 'v_5', 'v_6')
    
    
    V_phi_body_trajectory = zeros(number_of_time_steps, 6);
    V_phi_body_direction_trajectory = zeros(number_of_time_steps, 6);
    V_phi_body_normed_trajectory = zeros(number_of_time_steps, 6);
    for i_time = 1 : number_of_time_steps
        phi_dot = phi_dot_trajectory_relevant(i_time);
        rho_dot = rho_dot_trajectory_relevant(i_time);
        
        % calculate V_rho_body
        plant.jointAngles = joint_angle_trajectories(i_time, :)';
        plant.updateKinematics;
        J_body = plant.bodyJacobians{1};
        sphere_transformation = reshape(transformation_trajectories(i_time, :), 4, 4);
        rho_rotation_center_world = sphere_transformation * rho_rotation_center_local;
        p_virtual_contact = [plant.jointAngles(1:2); 0; 1];
        J_virtual_contact = plant.calculateArbitraryPointJacobian(p_virtual_contact, plant.numberOfJoints, 'world');
        
        A = [J_virtual_contact; [0 0 0 1 0 0];];
        A = J_virtual_contact;
        P_A = A' * (A * A')^(-1) * A; % projection onto the space spanned by the columns of A
        P_A_orth = eye(number_of_joints) - P_A; % projection onto the null space of A

        joint_velocity_rho_change_unconstrained = [0; 0; 0; 0; 0; 1];
        joint_velocity_rho_change_constrained = (P_A_orth * joint_velocity_rho_change_unconstrained);
        body_velocity_rho_change_constrained = J_body * joint_velocity_rho_change_constrained;
%         V_rho_body = body_velocity_rho_change_constrained * 1 / body_velocity_rho_change_constrained(5) * rho_dot;
        V_rho_body = body_velocity_rho_change_constrained * 1 / body_velocity_rho_change_constrained(5) * body_velocity_trajectory(i_time, 5);
        
        % calculate V_phi_body
        V_total_body = body_velocity_trajectory(i_time, :)';
        V_phi_body = V_total_body - V_rho_body;
        V_phi_body_direction = V_phi_body * phi_dot^(-1);
        
        V_phi_body_trajectory(i_time, :) = V_phi_body;
        V_phi_body_direction_trajectory(i_time, :) = V_phi_body_direction;
        V_phi_body_normed_trajectory(i_time, :) = normVector(V_phi_body_direction);
    end    
    
    relevant_data_points = abs(phi_dot_trajectory_relevant) > 1.0;
    irrelevant_data_points = ~relevant_data_points;
    phi_dot_trajectory_double_relevant = phi_dot_trajectory_relevant; phi_dot_trajectory_double_relevant(irrelevant_data_points) = NaN;
    V_phi_body_direction_trajectory_relevant = V_phi_body_direction_trajectory; V_phi_body_direction_trajectory_relevant(irrelevant_data_points, :) = NaN;
    V_phi_body_normed_trajectory_relevant = V_phi_body_normed_trajectory; V_phi_body_normed_trajectory_relevant(irrelevant_data_points, :) = NaN;
    
%     figure; axes; hold on
%     plot(phi_dot_trajectory);
%     plot(phi_dot_trajectory_relevant);
%     plot(phi_dot_trajectory_double_relevant);
%     
%     figure; axes; hold on
%     plot(rho_dot_trajectory);
%     plot(rho_dot_trajectory_relevant);
    
%     plot(V_phi_body_normed_trajectory(:, 1));
%     plot(V_phi_body_normed_trajectory_relevant(:, 1));
    
    
    % visualize - normalized body velocity vs. phi and rho
    figure; axes; hold on; title('v_1');
%     plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 1));
%     plot3(phi_trajectory, rho_trajectory, V_phi_body_trajectory(:, 1));
%     plot3(phi_trajectory, rho_trajectory, V_phi_body_normed_trajectory(:, 1));
    plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 1));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_1')
    
    figure; axes; hold on; title('v_2');
    plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 2));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_2')
    
    figure; axes; hold on; title('v_3');
    plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 3));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_3')
    
    figure; axes; hold on; title('v_4');
    plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 4));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_4')
    
    figure; axes; hold on; title('v_5');
    plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 5));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_5')
    
    figure; axes; hold on; title('v_6');
    plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 6));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_6')
    
end

%% calculate_allowed_body_velocities
if calculate_allowed_body_velocities

    phi_values = -0.8 : 0.1 : 0.5;
    rho_values = -0.5 : 0.05 : 0.5;
    number_of_phi_values = length(phi_values);
    number_of_rho_values = length(rho_values);
    body_velocity_phi_change_grid = zeros(number_of_phi_values, number_of_rho_values, 6);
    body_velocity_rho_change_grid = zeros(number_of_phi_values, number_of_rho_values, 6);
    for i_phi = 1 : number_of_phi_values
        for i_rho = 1 : number_of_rho_values
            plant.jointAngles = [0; 0; 0; 0; phi_values(i_phi); rho_values(i_rho)];
            plant.updateKinematics;
            
            p_contact = [plant.jointAngles(1:2); 0; 1];
            J_contact = plant.calculateArbitraryPointJacobian(p_contact, plant.numberOfJoints, 'world');
            A = [J_contact; [0 0 0 1 0 0];];
            P_A = A' * (A * A')^(-1) * A; % projection onto the space spanned by A
            P_A_orth = eye(number_of_joints) - P_A;

            joint_velocity_phi_change_unconstrained = [0; 0; 0; 0; 1; 0];
            joint_velocity_phi_change_constrained = (P_A_orth * joint_velocity_phi_change_unconstrained);
            joint_velocity_rho_change_unconstrained = [0; 0; 0; 0; 0; 1];
            joint_velocity_rho_change_constrained = (P_A_orth * joint_velocity_rho_change_unconstrained);
            
            J_body = plant.bodyJacobians{1};
            body_velocity_phi_change = plant.bodyJacobians{1} * joint_velocity_phi_change_constrained;
            body_velocity_phi_change_grid(i_phi, i_rho, :) = normVector(body_velocity_phi_change);
            body_velocity_rho_change = plant.bodyJacobians{1} * joint_velocity_rho_change_constrained;
            body_velocity_rho_change_grid(i_phi, i_rho, :) = normVector(body_velocity_rho_change);
            
        end
    end
    [phi_grid, rho_grid] = meshgrid(phi_values, rho_values);
    
    figure; axes; hold on; title('V_1 - allowed')
    surf(phi_grid, rho_grid, squeeze(body_velocity_phi_change_grid(:, :, 1))');
    surf(phi_grid, rho_grid, squeeze(body_velocity_rho_change_grid(:, :, 1))');
    plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 1));
    plot3(phi_trajectory, rho_trajectory, V_phi_body_normed_trajectory_relevant(:, 1));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_1')
    
    figure; axes; hold on; title('V_2 - allowed')
    surf(phi_grid, rho_grid, squeeze(body_velocity_phi_change_grid(:, :, 2))');
    surf(phi_grid, rho_grid, squeeze(body_velocity_rho_change_grid(:, :, 2))');
    plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 2));
    plot3(phi_trajectory, rho_trajectory, V_phi_body_normed_trajectory_relevant(:, 2));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_1')
    
    figure; axes; hold on; title('V_3 - allowed')
    surf(phi_grid, rho_grid, squeeze(body_velocity_phi_change_grid(:, :, 3))');
    surf(phi_grid, rho_grid, squeeze(body_velocity_rho_change_grid(:, :, 3))');
    plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 3));
    plot3(phi_trajectory, rho_trajectory, V_phi_body_normed_trajectory_relevant(:, 3));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_1')
    
    figure; axes; hold on; title('V_4 - allowed')
    surf(phi_grid, rho_grid, squeeze(body_velocity_phi_change_grid(:, :, 4))');
    surf(phi_grid, rho_grid, squeeze(body_velocity_rho_change_grid(:, :, 4))');
    plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 4));
    plot3(phi_trajectory, rho_trajectory, V_phi_body_normed_trajectory_relevant(:, 4));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_1')
    
    figure; axes; hold on; title('V_5 - allowed')
    surf(phi_grid, rho_grid, squeeze(body_velocity_phi_change_grid(:, :, 5))');
    surf(phi_grid, rho_grid, squeeze(body_velocity_rho_change_grid(:, :, 5))');
    plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 5));
    plot3(phi_trajectory, rho_trajectory, V_phi_body_normed_trajectory_relevant(:, 5));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_1')
    
    figure; axes; hold on; title('V_6 - allowed')
    surf(phi_grid, rho_grid, squeeze(body_velocity_phi_change_grid(:, :, 6))');
    surf(phi_grid, rho_grid, squeeze(body_velocity_rho_change_grid(:, :, 6))');
    plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 6));
    plot3(phi_trajectory, rho_trajectory, V_phi_body_normed_trajectory_relevant(:, 6));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_1')
    
end

%% show stick figure
if show_stick_figure
    scene_limits = 3*[-1 1; -1 1; -1 1];
    stick_figure = showMarkerComparisonStickFigure(plant, joint_angle_trajectories, [], scene_limits);

    resolution = 20;
    [x,y,z] = ...
      ellipsoid ...
      ( ...
        0, 0, 0, ...
        r_x, r_y, r_z, ...
        resolution ...
      );
    shape_data = ones(4, resolution+1, resolution+1);
    shape_data(1, :, :) = x;
    shape_data(2, :, :) = y;
    shape_data(3, :, :) = z;

    stick_figure.linkMassShapeData{6} = shape_data;
    stick_figure.linkMassShapeSurfs(6) = surf ...
      ( stick_figure.sceneAxes, ...
        x, y, z ...
      );
    stick_figure.update;
end






