% function ...
%   [ ...
%     polyfit_translation_x, ...
%     polyfit_translation_y, ...
%     polyfit_translation_z, ...
%     polyfit_rotation_x, ...
%     polyfit_rotation_y, ...
%     polyfit_rotation_z ...
%   ] ...
function [phi_trajectory, rho_trajectory, ankle_z_trajectory] ...
  = estimateBodyVelocityConstraint ...
  ( ...
    joint_frame_trajectory, ...
    body_velocity_trajectory, ...
    contact_indicator_trajectory, ...
    sampling_rate ...
  )
    number_of_time_steps = size(body_velocity_trajectory, 1);

    
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
    
    
    
    
    
    
    
end
