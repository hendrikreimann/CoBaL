

create_plant                        = 1;
generate_trajectories               = 1;
find_rho_constraint_by_phi          = 1;
find_phi_constraint                 = 1;
check_constraints                   = 1;
show_stick_figure                   = 0;

number_of_sets = 25;

%% create plant
if create_plant
    joint_positions = ...
      { ...
        [0; 0; 1], [0; 0; 1], [0; 0; 1], [0; 0; 1], [0; 0; 1], [0; 0; 1] ... % pelvis free body
        [0.1; 0; 0.9], [0.1; 0; 0.9], [0.1; 0; 0.9], ... % hip
        [0.1; 0; 0.5], ...  % knee
        [0.1; 0; 0.1], [0.1; 0; 0.1], ... % ankle
        [0; 0; 1.2], [0; 0; 1.2], [0; 0; 1.2] ... % lumbar
      };
    joint_axes = ...
      { ...
        [1; 0; 0], [0; 1; 0], [0; 0; 1], [1; 0; 0], [0; 1; 0], [0; 0; 1] ... % pelvis free body
        [1; 0; 0], [0; -1; 0], [0; 0; 1], ... % hip
        [-1; 0; 0], ...  % knee
        [1; 0; 0], [0; 1; 0], ... % ankle
        [1; 0; 0], [0; 1; 0], [0; 0; 1] ... % lumbar
      };
    joint_types = ...
      [ ...
        2 2 2 1 1 1 ... % virtual dofs
        1 1 1 1 1 1 ... % leg
        1 1 1 ...       % lumbar
      ];
    link_com_positions = ...
      { ...
        [0; 0; 1], [0; 0; 1], [0; 0; 1], [0; 0; 1], [0; 0; 1], [0; 0; 1] ... % pelvis free body
        [0.1; 0; 0.7], [0.1; 0; 0.7], [0.1; 0; 0.7], ... % hip
        [0.1; 0; 0.3], ...  % knee
        [0.1; 0; 0.1], [0.1; 0; 0.1], ... % ankle
        [0; 0; 1.5], [0; 0; 1.5], [0; 0; 1.5] ... % lumbar
      };
    link_orientations = ...
      { ...
        eye(3); eye(3); eye(3); eye(3); eye(3); eye(3); ... % pelvis
        eye(3); eye(3); eye(3); eye(3); eye(3); eye(3); ... % leg
        eye(3); eye(3); eye(3);  ... % torso
      };
    generalized_inertia_matrix_pelvis           = [10*eye(3) zeros(3); zeros(3) 0.01*eye(3)];
    generalized_inertia_matrix_thigh            = [5*eye(3) zeros(3); zeros(3) 0.01*eye(3)];
    generalized_inertia_matrix_shank            = [3*eye(3) zeros(3); zeros(3) 0.01*eye(3)];
    generalized_inertia_matrix_foot             = [1*eye(3) zeros(3); zeros(3) 0.01*eye(3)];
    generalized_inertia_matrix_torso            = [50*eye(3) zeros(3); zeros(3) 0.01*eye(3)];
    generalized_link_inertia_matrices = ...
      { ...
        zeros(6, 6), zeros(6, 6), zeros(6, 6), zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_pelvis, ...
        zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_thigh, generalized_inertia_matrix_shank, zeros(6, 6), generalized_inertia_matrix_foot, ...
        zeros(6, 6), zeros(6, 6), generalized_inertia_matrix_torso ...
      };
    end_effector_transformations = ...
      { ...
        [eye(3), [0.1; 0; 0]; 0 0 0 1], ... % ankle
        [eye(3), [0; 0; 1.8]; 0 0 0 1], ... % head
      };
    branch_matrix = ...
      [ ...
        1 1 1 1 1 1    1 1 1 1 1 1   0 0 0; ... % leg until heel
        1 1 1 1 1 1    0 0 0 0 0 0   1 1 1; ... % torso until head
      ];
  
    plant = GeneralKinematicTree ...
    ( ...
      joint_positions, ...
      joint_axes, ...
      joint_types, ...
      branch_matrix, ...
      end_effector_transformations, ...
      link_com_positions, ...
      link_orientations, ...
      generalized_link_inertia_matrices ...
    );
    plant.updateInternals();
    number_of_joints = plant.numberOfJoints;
    
    plant.jointLabels = ...
      { ...
        'pelvis, x-translation', 'pelvis, y-translation', 'pelvis, z-translation', 'pelvis, x-rotation', 'pelvis, y-rotation', 'pelvis, z-rotation', ...
        'hip flexion/extension', 'hip ab/adduction', 'hip external/internal rotation', 'knee flexion/extension', 'ankle dorsi/plantarflexion', 'ankle inversion/eversion', ...
        'lumbar joint - backward/forward bending', 'lumbar joint - sideways bending (right/left)', 'lumbar joint - internal rotation (left/right)', ...
      };
end

%% generate trajectories
if generate_trajectories
    total_time = 1;
    time_step = 0.01;
    time = (time_step : time_step : total_time)';
    number_of_time_steps_set = length(time);
    inversion = linspace(-1, 1, number_of_sets);
    velocity_factor = rand(1, number_of_sets)*0.4 + 0.8;
    foot_sphere_radius = 0.1;

    joint_angle_trajectories = zeros(number_of_time_steps_set * number_of_sets, number_of_joints);
    transformation_trajectories = zeros(number_of_time_steps_set * number_of_sets, 16);
    body_velocity_trajectory = zeros(number_of_time_steps_set * number_of_sets, 6);
    for i_set = 1 : number_of_sets
        % set start
        plant.jointAngles = zeros(number_of_joints, 1);
        plant.jointAngles(1) = (rand-0.5)*0.5; % randomize initial dorsiflexion
        plant.jointAngles(2) = (rand-0.5)*0.5; % randomize initial dorsiflexion
        plant.jointAngles(6) = (rand-0.5)*2.5; % randomize initial dorsiflexion
        plant.jointAngles(7) = 0.15 + (rand-0.5)*0.1; % randomize initial dorsiflexion
        plant.jointAngles(10) = 0.3 + (rand-0.5)*0.2; % randomize initial dorsiflexion
        plant.jointAngles(11) = 0.8 + (rand-0.5)*0.1; % randomize initial dorsiflexion
        plant.jointAngles(12) = 0 + (rand-0.5)*0.3;   % randomize initial inversion
        plant.updateKinematics;
        
        % correct foot to ground
        current_ankle_height = plant.jointTransformations{12}(3, 4);
        plant.jointAngles(3) = foot_sphere_radius - current_ankle_height;
        plant.updateKinematics;
        new_ankle_height = plant.jointTransformations{12}(3, 4);
        
        amplitudes = [zeros(1, 10) -1.5 2*inversion(i_set) zeros(1, 3)] * 1;
        phase_offsets = [zeros(1, 10) 0 0.5 zeros(1, 3)];
        frequencies = [zeros(1, 10) 0.5 0.5 zeros(1, 3)] * 1;
        offsets = [zeros(1, 10) 1.2 0 zeros(1, 3)] * 1;
        sinusoids = ...
          ( ...
            -cos( ...
                  (repmat(time, 1, number_of_joints) + repmat(phase_offsets, number_of_time_steps_set, 1)) ...
                  *2*pi.*repmat(frequencies, number_of_time_steps_set, 1) ...
                ) ...
            + repmat(offsets, number_of_time_steps_set, 1) ...
          ) ...
            .* repmat(amplitudes, number_of_time_steps_set, 1);
%         plot(sinusoids(:, 11:12));

        joint_angle_trajectories_set = zeros(number_of_time_steps_set, number_of_joints);
        joint_angle_trajectories_set(1, :) = plant.jointAngles;
        transformation_trajectories_set = zeros(number_of_time_steps_set, 16);
        body_velocity_trajectory_set = zeros(number_of_time_steps_set, 6);

        joint_velocity_trajectories_unconstrained = sinusoids * velocity_factor(i_set);
        joint_velocity_trajectories_constrained = zeros(number_of_time_steps_set, number_of_joints);

        transformation_trajectories_set(1, :) = reshape(plant.endEffectorTransformations{1}, 1, 16);
        for i_time = 2 : number_of_time_steps_set
                p_ankle = plant.jointTransformations{12}(:, 4);
                p_contact = p_ankle - [0; 0; foot_sphere_radius; 1];
                J_body_virtual_contact = plant.calculateArbitraryFrameBodyJacobian([eye(4, 3) p_contact], 12);
                A = J_body_virtual_contact([1:3 6], :);
                
                P_A = A' * (A * A')^(-1) * A; % projection onto the space spanned by A
                P_A_orth = eye(number_of_joints) - P_A;

                joint_velocity_unconstrained = joint_velocity_trajectories_unconstrained(i_time, :);
                joint_velocity_constrained = (P_A_orth * joint_velocity_unconstrained')';

                joint_velocity_constrained = joint_velocity_constrained * joint_velocity_unconstrained(11)/joint_velocity_constrained(11);

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
    
    
    
    
    
    
    
    
    
    
    

end

%% find rho constraint by phi
if find_rho_constraint_by_phi
    % calculate z-trajectories for multiple points along the sphere z-axis
    number_of_time_steps = size(transformation_trajectories, 1);
    phi_trajectory = zeros(number_of_time_steps, 1);
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
    body_velocity_relevant_trajectory = body_velocity_trajectory;   body_velocity_relevant_trajectory(irrelevant_data_points, :) = NaN;
    
    phi_dot_trajectory = deriveByTime(phi_trajectory, time_step);
    rho_dot_trajectory = deriveByTime(rho_trajectory, time_step);
    
    
%     figure; axes; hold on
%     plot(phi_dot_trajectory_relevant, 'linewidth', 2);
%     plot(rho_dot_trajectory_relevant, 'linewidth', 2);
%     plot(body_velocity_trajectory);
%     legend('d\phi/dt', 'd\rho/dt', 'v_1', 'v_2', 'v_3', 'v_4', 'v_5', 'v_6')
    
    
    V_phi_body_trajectory = zeros(number_of_time_steps, 6);
    V_phi_body_direction_trajectory = zeros(number_of_time_steps, 6);
    V_phi_body_normed_trajectory = zeros(number_of_time_steps, 6);
    body_velocity_normed_trajectory = zeros(number_of_time_steps, 6);
    for i_time = 1 : number_of_time_steps
        phi_dot = phi_dot_trajectory_relevant(i_time);
        rho_dot = rho_dot_trajectory_relevant(i_time);
        
        % calculate V_rho_body
        plant.jointAngles = joint_angle_trajectories(i_time, :)';
        plant.updateKinematics;
        J_body = plant.bodyJacobians{1};
        sphere_transformation = reshape(transformation_trajectories(i_time, :), 4, 4);
        rho_rotation_center_world = sphere_transformation * rho_rotation_center_local;
        
        p_virtual_contact = [rho_rotation_center_world(1:2); 0; 1];
        J_body_virtual_contact = plant.calculateArbitraryFrameBodyJacobian([eye(4, 3) p_virtual_contact], 12);
        A = J_body_virtual_contact([1:3 6], :);
        
        P_A = A' * (A * A')^(-1) * A; % projection onto the space spanned by the columns of A
        P_A_orth = eye(number_of_joints) - P_A; % projection onto the null space of A

        joint_velocity_rho_change_unconstrained = [zeros(11, 1); 1; zeros(3, 1)];
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
        
        body_velocity_normed_trajectory(i_time, :) = normVector(body_velocity_trajectory(i_time, :));
    end    
    body_velocity_normed_trajectory = -body_velocity_normed_trajectory;

    relevant_data_points = abs(phi_dot_trajectory_relevant) > 1.0;
    irrelevant_data_points = ~relevant_data_points;
    phi_dot_trajectory_double_relevant = phi_dot_trajectory_relevant; phi_dot_trajectory_double_relevant(irrelevant_data_points) = NaN;
    V_phi_body_direction_trajectory_relevant = V_phi_body_direction_trajectory; V_phi_body_direction_trajectory_relevant(irrelevant_data_points, :) = NaN;
    V_phi_body_normed_trajectory_relevant = V_phi_body_normed_trajectory; V_phi_body_normed_trajectory_relevant(irrelevant_data_points, :) = NaN;

%     figure; axes; hold on; title('v_1');
% %     plot(V_phi_body_trajectory(:, 1));
%     plot(V_phi_body_direction_trajectory(:, 1), 'linewidth', 2);
%     plot(V_phi_body_normed_trajectory(:, 1));
%     legend('body', 'direction', 'normed');
%     
%     figure; axes; hold on; title('v_2');
% %     plot(V_phi_body_trajectory(:, 2));
%     plot(V_phi_body_direction_trajectory(:, 2), 'linewidth', 2);
%     plot(V_phi_body_normed_trajectory(:, 2));
%     
%     figure; axes; hold on; title('v_3');
% %     plot(V_phi_body_trajectory(:, 3));
%     plot(V_phi_body_direction_trajectory(:, 3), 'linewidth', 2);
%     plot(V_phi_body_normed_trajectory(:, 3));
%     
%     figure; axes; hold on; title('v_4');
% %     plot(V_phi_body_trajectory(:, 4));
%     plot(V_phi_body_direction_trajectory(:, 4), 'linewidth', 2);
%     plot(V_phi_body_normed_trajectory(:, 4));
%     
%     figure; axes; hold on; title('v_5');
% %     plot(V_phi_body_trajectory(:, 5));
%     plot(V_phi_body_direction_trajectory(:, 5), 'linewidth', 2);
%     plot(V_phi_body_normed_trajectory(:, 5));
%     
%     figure; axes; hold on; title('v_6');
% %     plot(V_phi_body_trajectory(:, 6));
%     plot(V_phi_body_direction_trajectory(:, 6), 'linewidth', 2);
%     plot(V_phi_body_normed_trajectory(:, 6));
    
    
    
    
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
    
    
%     % visualize - normalized body velocity vs. phi and rho
%     figure; axes; hold on; title('v_1');
% %     plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 1));
% %     plot3(phi_relevant_trajectory, rho_relevant_trajectory, body_velocity_relevant_trajectory(:, 1));
% %     plot3(phi_trajectory, rho_trajectory, V_phi_body_trajectory(:, 1));
% %     plot3(phi_trajectory, rho_trajectory, V_phi_body_direction_trajectory_relevant(:, 1), 'linewidth', 2);
%     plot3(phi_trajectory, rho_trajectory, V_phi_body_normed_trajectory(:, 1));
%     plot3(phi_trajectory, rho_trajectory, body_velocity_normed_trajectory(:, 1));
%     xlabel('\phi'); ylabel('\rho'); zlabel('v_1')
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
    data_points_for_fit = 1 : 1 : length(phi_trajectory);
    data_points_for_fit(isnan(phi_trajectory(data_points_for_fit))) = [];
    data_points_for_fit(isnan(rho_trajectory(data_points_for_fit))) = [];
    data_points_for_fit(isnan(V_phi_body_normed_trajectory(data_points_for_fit))) = [];
    V_phi_body_1_fit = fit([phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 1), 'poly55');
    V_phi_body_2_fit = fit([phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 2), 'poly55');
    V_phi_body_3_fit = fit([phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 3), 'poly55');
    V_phi_body_4_fit = fit([phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 4), 'poly55');
    V_phi_body_5_fit = fit([phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 5), 'poly55');
    V_phi_body_6_fit = fit([phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 6), 'poly55');

    figure; plot(V_phi_body_1_fit, [phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 1));
    figure; plot(V_phi_body_2_fit, [phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 2));
    figure; plot(V_phi_body_3_fit, [phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 3));
    figure; plot(V_phi_body_4_fit, [phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 4));
    figure; plot(V_phi_body_5_fit, [phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 5));
    figure; plot(V_phi_body_6_fit, [phi_trajectory(data_points_for_fit), rho_trajectory(data_points_for_fit)], V_phi_body_normed_trajectory(data_points_for_fit, 6));
    
end

%% check constraints
if check_constraints
    % set configuration
    phi_check = 0.4;
    rho_check = 0.1;
    
    plant.jointAngles = zeros(number_of_joints, 1);
    plant.jointAngles(11) = phi_check;
    plant.jointAngles(12) = rho_check;
    plant.updateKinematics;
    
    % calculate actual constraint
    p_ankle = plant.jointTransformations{12}(:, 4);
    p_contact = p_ankle - [0; 0; foot_sphere_radius; 1];
    J_body_virtual_contact = plant.calculateArbitraryFrameBodyJacobian([eye(4, 3) p_contact], 12);
    A_actual = J_body_virtual_contact([1:3 6], :);

    % calculate allowed body velocities
    rho_rotation_center_local = [0; 0; 0; 1]; % eventually this might depend upon phi, but for now, we can assume it's fixed
    J_body = plant.bodyJacobians{1};
    sphere_transformation = reshape(plant.jointTransformations{12}, 4, 4);
    rho_rotation_center_world = sphere_transformation * rho_rotation_center_local;
    p_virtual_contact = [rho_rotation_center_world(1:2); 0; 1];
    J_body_virtual_contact = plant.calculateArbitraryFrameBodyJacobian([eye(4, 3) p_virtual_contact], 12);
    A_rho = J_body_virtual_contact([1:3 6], :);
    P_A_rho = A_rho' * (A_rho * A_rho')^(-1) * A_rho; % projection onto the space spanned by the columns of A
    P_A_rho_orth = eye(number_of_joints) - P_A_rho; % projection onto the null space of A

    joint_velocity_rho_change_unconstrained = [zeros(11, 1); 1; zeros(3, 1)];
    joint_velocity_rho_change_constrained = (P_A_rho_orth * joint_velocity_rho_change_unconstrained);
    body_velocity_rho_change_constrained = J_body * joint_velocity_rho_change_constrained;
    V_rho_body = body_velocity_rho_change_constrained * 1 / body_velocity_rho_change_constrained(5);
    
    
    V_phi_body = zeros(6, 1);
    V_phi_body(1) = feval(V_phi_body_1_fit,[phi_check, rho_check]);
    V_phi_body(2) = feval(V_phi_body_2_fit,[phi_check, rho_check]);
    V_phi_body(3) = feval(V_phi_body_3_fit,[phi_check, rho_check]);
    V_phi_body(4) = feval(V_phi_body_4_fit,[phi_check, rho_check]);
    V_phi_body(5) = feval(V_phi_body_5_fit,[phi_check, rho_check]);
    V_phi_body(6) = feval(V_phi_body_6_fit,[phi_check, rho_check]);
    
    % calculate allowed joint velocities
        % is this restrained enough? Maybe I introduce an unwanted element here?
    W_phi = pinv(J_body) * V_phi_body;
    W_rho = pinv(J_body) * V_rho_body;
    [U_body, S_body, V_body] = svd(J_body);
    B = V_body(:, 7 : 15); % this is the null space of the body Jacobian, i.e. the space of joint changes that do not move the foot
    C = [B W_phi W_rho]; % this is the space of unconstrained joint changes
    [U_constraint, S_constraint, V_constraint] = svd(C');
    A_constraint = V_constraint(:, 12:15)';
    
    
    % compare this with the solution from createConstraintMatrix_bodyVelocityConstraints_24.m
    V_body_allowed = [V_phi_body V_rho_body];
    [~, ~, V] = svd(V_body_allowed');
    C = V(:, 3:6)'; % orthogonal complement of the allowed body velocity directions
    A_check = C * J_body;
    
    
    
    % A_actual, A_constraint and A_check should all span the same subspace
end



%% show stick figure
if show_stick_figure
    scene_limits = 1*[-1 1; -1 1; -0.1 1.9];
    
    % show stick figure of geometric model
%     stick_figure = KinematicTreeController(plant, scene_limits, 'ellipsoid');
    
    
%     joint_angle_trajectories = zeros(2, number_of_joints);
    stick_figure = showMarkerComparisonStickFigure(plant, joint_angle_trajectories, [], scene_limits);

    resolution = 20;
    [x,y,z] = ...
      ellipsoid ...
      ( ...
        0, 0, 0, ...
        foot_sphere_radius, foot_sphere_radius, foot_sphere_radius, ...
        resolution ...
      );
    shape_data = ones(4, resolution+1, resolution+1);
    shape_data(1, :, :) = x;
    shape_data(2, :, :) = y;
    shape_data(3, :, :) = z;

    stick_figure.linkMassShapeData{12} = shape_data;
    stick_figure.linkMassShapeSurfs(12) = surf ...
      ( stick_figure.sceneAxes, ...
        x, y, z ...
      );
    stick_figure.update;
end

