% function ...
%   [ ...
%     polyfit_translation_x, ...
%     polyfit_translation_y, ...
%     polyfit_translation_z, ...
%     polyfit_rotation_x, ...
%     polyfit_rotation_y, ...
%     polyfit_rotation_z ...
%   ] ...
% function [phi_trajectory, rho_trajectory, ankle_z_trajectory] ...
%   = estimateBodyVelocityConstraint ...
%   ( ...
%     joint_frame_trajectory, ...
%     body_velocity_trajectory, ...
%     contact_indicator_trajectory, ...
%     sampling_rate ...
%   )

% load data
trial_number = 2;
load subjectInfo.mat;
load(makeFileName(date, subject_id, 'model'));
load(makeFileName(date, subject_id, 'walking', trial_number, 'markerTrajectories'));
load(makeFileName(date, subject_id, 'walking', trial_number, 'angleTrajectories'));
load(makeFileName(date, subject_id, 'walking', trial_number, 'kinematicTrajectories'));
load(makeFileName(date, subject_id, 'walking', trial_number, 'stepEvents'));



    number_of_time_steps = size(T_left_ankle_to_world_trajectory, 1);

    calculate_trajectories              = 0;
    find_rho_constraint_by_phi          = 0;
    find_phi_constraint                 = 1;

    number_of_points_of_interest = 200;
    number_of_phi_values = 200;
    number_of_joints = plant.numberOfJoints;
    
    %% calculate trajectories
    if calculate_trajectories
        
        % calculate z-trajectories for multiple points along the sphere z-axis
        phi_trajectory = zeros(number_of_time_steps, 1);
        rho_trajectory = zeros(number_of_time_steps, 1);
        gamma_trajectory = zeros(number_of_time_steps, 1);

        points_of_interest_z_local = linspace(-0.02, 0.08, number_of_points_of_interest);
        points_of_interest_sphere = [zeros(2, number_of_points_of_interest); points_of_interest_z_local; ones(1, number_of_points_of_interest)]; % coordinates of the points of interest in sphere coordinates
        points_of_interest_world_trajectories_x = zeros(number_of_time_steps, size(points_of_interest_sphere, 2));
        points_of_interest_world_trajectories_y = zeros(number_of_time_steps, size(points_of_interest_sphere, 2));
        points_of_interest_world_trajectories_z = zeros(number_of_time_steps, size(points_of_interest_sphere, 2));
        for i_time = 1 : number_of_time_steps
            ankle_scs_transformation_current = reshape(T_left_ankle_to_world_trajectory(i_time, :), 4, 4);
            ankle_scs_to_world_rotation_current = ankle_scs_transformation_current(1:3, 1:3);

            euler_angles = eulerAnglesFromRotationMatrixYZX(ankle_scs_to_world_rotation_current);
            gamma_trajectory(i_time) = euler_angles_abc(1);
            phi_trajectory(i_time) = euler_angles_abc(2);
            rho_trajectory(i_time) = euler_angles_abc(3);

            points_of_interest_world_current = ankle_scs_transformation_current * points_of_interest_sphere;

            points_of_interest_world_trajectories_x(i_time, :) = points_of_interest_world_current(1, :);
            points_of_interest_world_trajectories_y(i_time, :) = points_of_interest_world_current(2, :);
            points_of_interest_world_trajectories_z(i_time, :) = points_of_interest_world_current(3, :);
        end
        
        % get relevant data points
        relevant_data_points = left_contact_indicators_mocap;
        irrelevant_data_points = ~relevant_data_points;

        phi_trajectory_relevant = phi_trajectory; phi_trajectory_relevant(irrelevant_data_points) = NaN;
        rho_trajectory_relevant = rho_trajectory; rho_trajectory_relevant(irrelevant_data_points) = NaN;
        
        
        phi_relevant_trajectory = phi_trajectory;   phi_relevant_trajectory(irrelevant_data_points) = NaN;
        rho_relevant_trajectory = rho_trajectory;   rho_relevant_trajectory(irrelevant_data_points) = NaN;

        phi_dot_trajectory = deriveByTime(phi_trajectory, 1/sampling_rate_mocap);
        rho_dot_trajectory = deriveByTime(rho_trajectory, 1/sampling_rate_mocap);
        gamma_dot_trajectory = deriveByTime(gamma_trajectory, 1/sampling_rate_mocap);
        phi_dot_trajectory_relevant = deriveByTime(phi_relevant_trajectory, 1/sampling_rate_mocap);
        rho_dot_trajectory_relevant = deriveByTime(rho_relevant_trajectory, 1/sampling_rate_mocap);
        
        body_velocity_relevant_trajectory = V_body_left_ankle;   body_velocity_relevant_trajectory(irrelevant_data_points, :) = NaN;
        V_body_left_ankle_trajectory_relevant = V_body_left_ankle; V_body_left_ankle_trajectory_relevant(irrelevant_data_points, :) = NaN;
        
        % compare parameter velocities with body velocity
        figure; axes; hold on
        plot(V_body_left_ankle(:, 4:6));
        plot(gamma_dot_trajectory);
        plot(phi_dot_trajectory);
        plot(rho_dot_trajectory);
        legend('V_4', 'V_5', 'V_6', 'd\gamma / dt', 'd\phi / dt', 'd\rho / dt');
    end
    
    %% find rho constraint by phi
    if find_rho_constraint_by_phi
        
        % fit surfaces to phi-rho-z data
        data_points_for_fit = find(left_contact_indicators_mocap);
        phi_fitpoints = phi_trajectory(data_points_for_fit);
        rho_fitpoints = rho_trajectory(data_points_for_fit);
        poi_z_fits = cell(1, number_of_points_of_interest);
        for i_poi = 1 : number_of_points_of_interest
            p_z = points_of_interest_world_trajectories_z(data_points_for_fit, i_poi);
            poi_z_fits{i_poi} = fit([phi_fitpoints, rho_fitpoints], p_z, 'poly55');
%             figure; plot(poi_z_fits{i_poi}, [phi_fitpoints, rho_fitpoints], points_of_interest_world_trajectories_z(data_points_for_fit, i_poi));
%             xlabel('\phi'); ylabel('\rho'); zlabel('z'); 
        end
        
        % for a grid of each phi values, find the point of interest that moves least under changes of rho
        phi_values = linspace(min(phi_fitpoints), max(phi_fitpoints), number_of_phi_values);
        
        z_values = cell(number_of_phi_values, number_of_points_of_interest);
        poi_z_fits_by_rho = cell(number_of_phi_values, number_of_points_of_interest);
        poi_z_fits_by_rho_error = zeros(number_of_phi_values, number_of_points_of_interest);
        poi_z_fits_by_rho_derivative = zeros(number_of_phi_values, number_of_points_of_interest);
        poi_z_fits_by_rho_derivative_error = zeros(number_of_phi_values, number_of_points_of_interest);
        poi_z_fits_by_rho_rms = zeros(number_of_phi_values, number_of_points_of_interest);
        
%         figure; axes; hold on
%         plot(rho_trajectory);
        phi_window_radius = 0.025;
        for i_phi = 1 : length(phi_values)
            % for each point of interest, compare z-values to straight line
            rho_windows = zeros(number_of_phi_values, 2);
            phi_window = phi_values(i_phi) + [-1 1] * phi_window_radius;
            phi_window_indices = find(phi_window(1) <= phi_trajectory & phi_trajectory <= phi_window(2));
%             plot(phi_window_indices, rho_trajectory(phi_window_indices), '-')                
            rho_windows(i_phi, :) = [min(rho_trajectory(phi_window_indices)) max(rho_trajectory(phi_window_indices))];
%             plot([i_phi i_phi], rho_windows(i_phi, :));
            rho_values = linspace(rho_windows(i_phi, 1), rho_windows(i_phi, 2), 20);

    %         figure; hold on
            for i_poi = 1 : number_of_points_of_interest
                z_values{i_phi, i_poi} = feval(poi_z_fits{i_poi}, [ones(length(rho_values), 1)*phi_values(i_phi), rho_values']);
                z_values_rms = sum((z_values{i_phi, i_poi} - mean(z_values{i_phi, i_poi})).^2).^(0.5);
                poi_z_fits_by_rho_rms(i_phi, i_poi) = z_values_rms;
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
        [poi_z_fits_by_rho_rms_minima, poi_z_fits_by_rho_rms_minimum_indices] = min(poi_z_fits_by_rho_rms, [], 2);

        options = fitoptions('smoothingspline');
        options.SmoothingParam = 0.9999;
        rho_cor_fit = fit(phi_values', points_of_interest_z_local(poi_z_fits_by_rho_rms_minimum_indices)', 'smoothingspline', options);
        
        % plot minimum by phi value
        figure; plot(rho_cor_fit, phi_values, points_of_interest_z_local(poi_z_fits_by_rho_rms_minimum_indices));
        xlabel('\phi'); ylabel('z'); 

        figure; plot(phi_values, poi_z_fits_by_rho_rms_minima);
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
    end    
    
    
%% find_phi_constraint
if find_phi_constraint
    
    % calculate parameter trajectories (ACHTUNG: do I need to repeat this? This was already done in the previous section
%     phi_trajectory = zeros(number_of_time_steps, 1);
%     rho_trajectory = zeros(number_of_time_steps, 1);
%     gamma_trajectory = zeros(number_of_time_steps, 1);
%     for i_time = 1 : number_of_time_steps
%         ankle_scs_transformation_current = reshape(joint_frame_trajectory(i_time, :), 4, 4);
%         ankle_scs_to_world_rotation_current = ankle_scs_transformation_current(1:3, 1:3);
% 
%         ankle_scs_current_to_reference_rotation = ankle_scs_to_world_rotation_reference^(-1) * ankle_scs_to_world_rotation_current;
%         ankle_abc_current_to_reference_rotation = scs_to_abc_rotation * ankle_scs_current_to_reference_rotation * scs_to_abc_rotation^(-1);
% 
%         euler_angles_scs = eulerAnglesFromRotationMatrixZXY(ankle_scs_current_to_reference_rotation);
%         euler_angles_abc = eulerAnglesFromRotationMatrixZXY(ankle_abc_current_to_reference_rotation);
%         gamma_trajectory(i_time) = euler_angles_abc(1);
%         phi_trajectory(i_time) = euler_angles_abc(2);
%         rho_trajectory(i_time) = euler_angles_abc(3);
%     end    
    
    
%     figure; axes; hold on
%     plot(phi_trajectory)
%     plot(phi_trajectory_relevant)
    
    
    
    V_phi_body_trajectory = zeros(number_of_time_steps, 6);
    V_phi_body_direction_trajectory = zeros(number_of_time_steps, 6);
    V_phi_body_normed_trajectory = zeros(number_of_time_steps, 6);
    body_velocity_normed_trajectory = zeros(number_of_time_steps, 6);
    
    time_steps_to_fit = 1 : 3000;
%     for i_time = 1 : number_of_time_steps
    for i_time = time_steps_to_fit
        phi_dot = phi_dot_trajectory_relevant(i_time);
        
        % calculate z-offset of rho rotation center
        z_offset = feval(rho_cor_fit, phi_trajectory(i_time));
        rho_rotation_center_local = [0; 0; z_offset; 1];
        
        % calculate V_rho_body
        plant.jointAngles = angle_trajectories(i_time, :)';
        plant.updateKinematics;
        J_body = plant.bodyJacobians{3};
        ankle_transformation = reshape(T_left_ankle_to_world_trajectory(i_time, :), 4, 4);
        rho_rotation_center_world = ankle_transformation * rho_rotation_center_local;
        
        p_virtual_contact = [rho_rotation_center_world(1:2); 0; 1];
        J_body_virtual_contact = plant.calculateArbitraryFrameBodyJacobian([eye(4, 3) p_virtual_contact], 12);
        A_contact = J_body_virtual_contact([1:3 6], :);
        
        P_A_contact = A_contact' * (A_contact * A_contact')^(-1) * A_contact; % projection onto the space spanned by the columns of A
        P_A_contact_orth = eye(number_of_joints) - P_A_contact; % projection onto the null space of A

        % ACHTUNG: this is hard-coded for the left leg now, change it later (left ankle in-eversion is 12th joint)
        joint_velocity_rho_change_unconstrained = [zeros(11, 1); 1; zeros(26, 1)];
        joint_velocity_rho_change_constrained = (P_A_contact_orth * joint_velocity_rho_change_unconstrained);
        body_velocity_rho_change_constrained = J_body * joint_velocity_rho_change_constrained;
%         V_rho_body = body_velocity_rho_change_constrained * 1 / body_velocity_rho_change_constrained(5) * rho_dot;
        body_velocity_rho = V_body_left_ankle(i_time, 4);
        V_rho_body = body_velocity_rho_change_constrained * 1 / body_velocity_rho_change_constrained(4) * body_velocity_rho;
        
        % calculate V_phi_body
        V_total_body = V_body_left_ankle(i_time, :)';
        V_phi_body = V_total_body - V_rho_body;
        V_phi_body_direction = V_phi_body * phi_dot^(-1);
        
        V_phi_body_trajectory(i_time, :) = V_phi_body;
        V_phi_body_direction_trajectory(i_time, :) = V_phi_body_direction;
        V_phi_body_normed_trajectory(i_time, :) = normVector(V_phi_body_direction);
        
        body_velocity_normed_trajectory(i_time, :) = normVector(V_body_left_ankle(i_time, :));
        
        if i_time/10  == floor(i_time / 10)
            disp(num2str(i_time))
        end
    end
    
    % normalize direction
%     V_phi_body_normed_trajectory(V_phi_body_normed_trajectory(:, 1)<0, :) = - V_phi_body_normed_trajectory(V_phi_body_normed_trajectory(:, 1)<0, :);
%     body_velocity_normed_trajectory(body_velocity_normed_trajectory(:, 1)<0, :) = - body_velocity_normed_trajectory(body_velocity_normed_trajectory(:, 1)<0, :);

%     relevant_data_points = abs(phi_dot_trajectory_relevant) > 1.0;
%     irrelevant_data_points = ~relevant_data_points;
%     phi_dot_trajectory_double_relevant = phi_dot_trajectory_relevant; phi_dot_trajectory_double_relevant(irrelevant_data_points) = NaN;
    body_velocity_trajectory_relevant = V_body_left_ankle; body_velocity_trajectory_relevant(irrelevant_data_points, :) = NaN;
    body_velocity_normed_trajectory_relevant = body_velocity_normed_trajectory; body_velocity_normed_trajectory_relevant(irrelevant_data_points, :) = NaN;
    V_phi_body_trajectory_relevant = V_phi_body_trajectory; V_phi_body_trajectory_relevant(irrelevant_data_points, :) = NaN;
    V_phi_body_direction_trajectory_relevant = V_phi_body_direction_trajectory; V_phi_body_direction_trajectory_relevant(irrelevant_data_points, :) = NaN;
    V_phi_body_normed_trajectory_relevant = V_phi_body_normed_trajectory; V_phi_body_normed_trajectory_relevant(irrelevant_data_points, :) = NaN;

    %%
    figure; axes; hold on; title('v_1');
    plot(V_phi_body_trajectory(time_steps_to_fit, 1));
    plot(V_phi_body_trajectory_relevant(time_steps_to_fit, 1));
    plot(V_phi_body_direction_trajectory_relevant(time_steps_to_fit, 1), 'linewidth', 2);
    plot(phi_dot_trajectory(time_steps_to_fit, 1));
%     plot(V_phi_body_normed_trajectory_relevant(time_steps_to_fit, 1));
    legend('body', 'normed');
%     legend('body', 'direction', 'normed');
    return
    figure; axes; hold on; title('v_2');
    plot(V_phi_body_trajectory(time_steps_to_fit, 2));
    plot(V_phi_body_trajectory_relevant(time_steps_to_fit, 2));
    plot(V_phi_body_direction_trajectory_relevant(time_steps_to_fit, 2), 'linewidth', 2);
    plot(V_phi_body_normed_trajectory_relevant(time_steps_to_fit, 2));
    
    figure; axes; hold on; title('v_3');
    plot(V_phi_body_trajectory(time_steps_to_fit, 3));
    plot(V_phi_body_trajectory_relevant(time_steps_to_fit, 3));
    plot(V_phi_body_direction_trajectory_relevant(time_steps_to_fit, 3), 'linewidth', 2);
    plot(V_phi_body_normed_trajectory_relevant(time_steps_to_fit, 3));
    
    figure; axes; hold on; title('v_4');
    plot(V_phi_body_trajectory(time_steps_to_fit, 4));
    plot(V_phi_body_trajectory_relevant(time_steps_to_fit, 4));
    plot(V_phi_body_direction_trajectory_relevant(time_steps_to_fit, 4), 'linewidth', 2);
    plot(V_phi_body_normed_trajectory_relevant(time_steps_to_fit, 4));
    
    figure; axes; hold on; title('v_5');
    plot(V_phi_body_trajectory(time_steps_to_fit, 5));
    plot(V_phi_body_trajectory_relevant(time_steps_to_fit, 5));
    plot(V_phi_body_direction_trajectory_relevant(time_steps_to_fit, 5), 'linewidth', 2);
    plot(V_phi_body_normed_trajectory_relevant(time_steps_to_fit, 5));
    
    figure; axes; hold on; title('v_6');
    plot(V_phi_body_trajectory(time_steps_to_fit, 6));
    plot(V_phi_body_trajectory_relevant(time_steps_to_fit, 6));
    plot(V_phi_body_direction_trajectory_relevant(time_steps_to_fit, 6), 'linewidth', 2);
    plot(V_phi_body_normed_trajectory_relevant(time_steps_to_fit, 6));
    return
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
    
    
    % visualize - normalized body velocity vs. phi and rho
    figure; axes; hold on; title('v_1');
    plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 1));
    plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory_relevant(:, 1));
    plot3(phi_trajectory, rho_trajectory, body_velocity_normed_trajectory_relevant(:, 1));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);

    figure; axes; hold on; title('v_2');
    plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 2));
    plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory_relevant(:, 2));
    plot3(phi_trajectory, rho_trajectory, body_velocity_normed_trajectory_relevant(:, 2));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
    
    figure; axes; hold on; title('v_3');
    plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 3));
    plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory_relevant(:, 3));
    plot3(phi_trajectory, rho_trajectory, body_velocity_normed_trajectory_relevant(:, 3));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
    
    figure; axes; hold on; title('v_4');
    plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 4));
    plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory_relevant(:, 4));
    plot3(phi_trajectory, rho_trajectory, body_velocity_normed_trajectory_relevant(:, 4));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
    
    figure; axes; hold on; title('v_5');
    plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 5));
    plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory_relevant(:, 5));
    plot3(phi_trajectory, rho_trajectory, body_velocity_normed_trajectory_relevant(:, 5));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
    
    figure; axes; hold on; title('v_6');
    plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 6));
    plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory_relevant(:, 6));
    plot3(phi_trajectory, rho_trajectory, body_velocity_normed_trajectory_relevant(:, 6));
    xlabel('\phi'); ylabel('\rho'); zlabel('v_1'); view([0 -1 0]);
    
    
    
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
