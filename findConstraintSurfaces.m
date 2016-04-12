
use_parallel    = 1;

find_body_constraints               = 0;


% load data
trial_number = 2;
load subjectInfo.mat;
load(makeFileName(date, subject_id, 'model'));
load(makeFileName(date, subject_id, 'walking', trial_number, 'markerTrajectories'));
load(makeFileName(date, subject_id, 'walking', trial_number, 'angleTrajectories'));
load(makeFileName(date, subject_id, 'walking', trial_number, 'kinematicTrajectories'));
load(makeFileName(date, subject_id, 'walking', trial_number, 'stepEvents'));







[phi_trajectory, rho_trajectory, ankle_z_trajectory] ...
= estimateBodyVelocityConstraint(T_left_ankle_to_world_trajectory, V_body_left_ankle, left_contact_indicators_mocap, sampling_rate_mocap);



% set irrelevant data points to NaN
phi_threshold = 10.1;
relevant_data_points = phi_trajectory < phi_threshold;
irrelevant_data_points = ~relevant_data_points;
phi_relevant = phi_trajectory;         phi_relevant(irrelevant_data_points) = NaN;
rho_relevant = rho_trajectory;         rho_relevant(irrelevant_data_points) = NaN;
ankle_z_relevant = ankle_z_trajectory;     ankle_z_relevant(irrelevant_data_points, :) = NaN;

nan_data_points = any(isnan(phi_relevant), 2);
phi_nanless = phi_relevant(~nan_data_points);
rho_nanless = rho_relevant(~nan_data_points);
ankle_z_nanless = ankle_z_relevant(~nan_data_points);

return
% visualize - body velocity vs. phi
figure; axes; hold on; title('v_1');
plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 1));
plot3(phi_relevant, rho_relevant, body_velocity_relevant(:, 1));

figure; axes; hold on; title('v_2');
plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 2));
plot3(phi_relevant, rho_relevant, body_velocity_relevant(:, 2));

figure; axes; hold on; title('v_3');
plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 3));
plot3(phi_relevant, rho_relevant, body_velocity_relevant(:, 3));

figure; axes; hold on; title('v_4');
plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 4));
plot3(phi_relevant, rho_relevant, body_velocity_relevant(:, 4));

figure; axes; hold on; title('v_5');
plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 5));
plot3(phi_relevant, rho_relevant, body_velocity_relevant(:, 5));

figure; axes; hold on; title('v_6');
plot3(phi_trajectory, rho_trajectory, body_velocity_trajectory(:, 6));
plot3(phi_relevant, rho_relevant, body_velocity_relevant(:, 6));


return
    
    

% find left heel constraint
[ ...
  polyfit_translation_x, ...
  polyfit_translation_y, ...
  polyfit_translation_z, ...
  polyfit_rotation_x, ...
  polyfit_rotation_y, ...
  polyfit_rotation_z ...
] ...
= estimateBodyVelocityConstraint(T_left_ankle_to_world_trajectory, V_body_left_ankle, left_contact_indicators_mocap, sampling_rate_mocap);
polyfit_body_velocity_left_heel_roll = ...
[ ...
  polyfit_translation_x; ...
  polyfit_translation_y; ...
  polyfit_translation_z; ...
  polyfit_rotation_x; ...
  polyfit_rotation_y; ...
  polyfit_rotation_z ...
];













if find_body_constraints
    
    
    
    
    
    
    
    % find right heel constraint
    [ ...
      polyfit_translation_x, ...
      polyfit_translation_y, ...
      polyfit_translation_z, ...
      polyfit_rotation_x, ...
      polyfit_rotation_y, ...
      polyfit_rotation_z ...
    ] ...
    = estimateBodyVelocityConstraint(T_ankle_to_world_trajectory, V_body_left_ankle);

    polyfit_body_velocity_right_heel_roll = ...
    [ ...
      polyfit_translation_x; ...
      polyfit_translation_y; ...
      polyfit_translation_z; ...
      polyfit_rotation_x; ...
      polyfit_rotation_y; ...
      polyfit_rotation_z ...
    ];

    % find right toes constraint
    file_name = createFileName(date, subject_id, study_label, num2str(right_toes_reference_file_index), 'rightAnkleFrameTrajectory');
    load([data_root directorySeparator file_name]);
    file_name = createFileName(date, subject_id, study_label, num2str(right_toes_reference_file_index), 'rightAnkleFrameBodyVelocity');
    load([data_root directorySeparator file_name]);
    [ ...
      polyfit_translation_x, ...
      polyfit_translation_y, ...
      polyfit_translation_z, ...
      polyfit_rotation_x, ...
      polyfit_rotation_y, ...
      polyfit_rotation_z ...
    ] ...
    = estimateBodyVelocityConstraint(T_ankle_to_world_trajectory, V_body_left_ankle);
    polyfit_body_velocity_right_toes_roll = ...
    [ ...
      polyfit_translation_x; ...
      polyfit_translation_y; ...
      polyfit_translation_z; ...
      polyfit_rotation_x; ...
      polyfit_rotation_y; ...
      polyfit_rotation_z ...
    ];

    % find left heel constraint
    file_name = createFileName(date, subject_id, study_label, num2str(left_heel_reference_file_index), 'leftAnkleFrameTrajectory');
    load([data_root directorySeparator file_name]);
    file_name = createFileName(date, subject_id, study_label, num2str(left_heel_reference_file_index), 'leftAnkleFrameBodyVelocity');
    load([data_root directorySeparator file_name]);
    [ ...
      polyfit_translation_x, ...
      polyfit_translation_y, ...
      polyfit_translation_z, ...
      polyfit_rotation_x, ...
      polyfit_rotation_y, ...
      polyfit_rotation_z ...
    ] ...
    = estimateBodyVelocityConstraint(T_ankle_to_world_trajectory, V_body_left_ankle);
    polyfit_body_velocity_left_heel_roll = ...
    [ ...
      polyfit_translation_x; ...
      polyfit_translation_y; ...
      polyfit_translation_z; ...
      polyfit_rotation_x; ...
      polyfit_rotation_y; ...
      polyfit_rotation_z ...
    ];

    % find left toes constraint
    file_name = createFileName(date, subject_id, study_label, num2str(left_toes_reference_file_index), 'leftAnkleFrameTrajectory');
    load([data_root directorySeparator file_name]);
    file_name = createFileName(date, subject_id, study_label, num2str(left_toes_reference_file_index), 'leftAnkleFrameBodyVelocity');
    load([data_root directorySeparator file_name]);
    [ ...
      polyfit_translation_x, ...
      polyfit_translation_y, ...
      polyfit_translation_z, ...
      polyfit_rotation_x, ...
      polyfit_rotation_y, ...
      polyfit_rotation_z ...
    ] ...
    = estimateBodyVelocityConstraint(T_ankle_to_world_trajectory, V_body_left_ankle);
    polyfit_body_velocity_left_toes_roll = ...
    [ ...
      polyfit_translation_x; ...
      polyfit_translation_y; ...
      polyfit_translation_z; ...
      polyfit_rotation_x; ...
      polyfit_rotation_y; ...
      polyfit_rotation_z ...
    ];


    distFig('rows', 4)
    
    file_name = createFileName(date, subject_id, study_label, num2str(1), 'bodyVelocityConstraints');
    save( ...
          [data_root directorySeparator file_name], ...
          'polyfit_body_velocity_right_heel_roll', ...
          'polyfit_body_velocity_right_toes_roll', ...
          'polyfit_body_velocity_left_heel_roll', ...
          'polyfit_body_velocity_left_toes_roll' ...
        );    
end






return

%% here comes old stuff, kept for historical reasons as it might be useful at some point
% identify the plane of movement
circle_reconstr = cell(1, 4);
normals = zeros(3, 3);
for i_marker = 2 : 4
    [center, normal, radius, V] = circleFit(right_foot_markers_trajectory(:, (i_marker-1)*3+1 : (i_marker-1)*3+3));

    rotation = [V normal];
    transformation = [rotation center(1:3); 0 0 0 1];
    circle_reconstr_plane = [cos(0:0.01:2*pi)*radius; sin(0:0.01:2*pi)*radius; zeros(size(0:0.01:2*pi)); ones(size(0:0.01:2*pi))];
    circle_reconstr{i_marker} = transformation * circle_reconstr_plane;

    if dot(normal, [1; 0; 0]) < 0;
        normal = -normal;
    end
    normals(:, i_marker-1) = normal;
end
plane_normal = normVector(mean(normals, 2));    

% calculate cluster to world reference transformation
number_of_markers = length(right_foot_markers_reference) / 3;
right_foot_cluster_reference_world = reshape(right_foot_markers_reference, 3, number_of_markers);
p_right_foot_cluster_reference_to_world = mean(right_foot_cluster_reference_world, 2);
R_right_foot_cluster_reference_to_world = eye(3, 3);
T_right_foot_cluster_reference_to_world = [R_right_foot_cluster_reference_to_world p_right_foot_cluster_reference_to_world; 0 0 0 1];

foot_azimuth_trajectory = zeros(number_of_time_frames, 1);
foot_elevation_trajectory = zeros(number_of_time_frames, 1);
foot_radius_trajectory = zeros(number_of_time_frames, 1);

ankle_azimuth_trajectory = zeros(number_of_time_frames, 1);
ankle_elevation_trajectory = zeros(number_of_time_frames, 1);
ankle_radius_trajectory = zeros(number_of_time_frames, 1);

for i_time = 1:number_of_time_frames
    % calculate cluster to world transformation
    right_foot_cluster_current_world = reshape(right_foot_markers_trajectory(i_time, :), 3, number_of_markers);
    p_right_foot_cluster_current_to_world = mean(right_foot_cluster_current_world, 2);
    R_right_foot_cluster_current_to_world = rotationMatrixFromMarkers(right_foot_cluster_reference_world, right_foot_cluster_current_world);
    T_right_foot_cluster_current_to_world = [R_right_foot_cluster_current_to_world p_right_foot_cluster_current_to_world; 0 0 0 1];

    % calculate current position of joint center
    joint_center_current_world = T_right_foot_cluster_current_to_world * [right_ankle_cor_rightfoot; 1];
    joint_center_trajectory_world(i_time, :) = joint_center_current_world(1:3);

    % calculate angle of the cluster frame
    marker_one_current_position = right_foot_cluster_current_world(:, 1);
    marker_two_current_position = right_foot_cluster_current_world(:, 4);
    main_axis_current = marker_two_current_position - marker_one_current_position;
    [foot_azimuth_trajectory(i_time), foot_elevation_trajectory(i_time), foot_radius_trajectory(i_time)] = cart2sph(main_axis_current(1), main_axis_current(2), main_axis_current(3));
    
    % calculate angle of the ankle frame
    ankle_transformation_current = reshape(T_ankle_to_world_trajectory(i_time, :), 4, 4);
    ankle_frame_x_axis_current = ankle_transformation_current(1:3, 1);
    [ankle_azimuth_trajectory(i_time), ankle_elevation_trajectory(i_time), ankle_radius_trajectory(i_time)] = cart2sph(ankle_frame_x_axis_current(1), ankle_frame_x_axis_current(2), ankle_frame_x_axis_current(3));
end
joint_center_velocity_trajectory_world = centdiff(joint_center_trajectory_world, sampling_rate^(-1));

% transform velocity into ankle frame
joint_center_velocity_trajectory_ankle = zeros(size(joint_center_velocity_trajectory_world));
for i_time = 1:number_of_time_frames
    ankle_transformation_current = reshape(T_ankle_to_world_trajectory(i_time, :), 4, 4);
    if any(isnan(ankle_transformation_current))
        joint_center_velocity_trajectory_ankle(i_time, :) = NaN;
    else
        joint_center_velocity_world_current = [joint_center_velocity_trajectory_world(i_time, :)'; 0];
        joint_center_velocity_ankle_current = ankle_transformation_current^(-1) * joint_center_velocity_world_current;
        joint_center_velocity_trajectory_ankle(i_time, :) = joint_center_velocity_ankle_current(1:3);
    end
end

% exclude the velocity vectors with too small magnitude
speed_threshold = 0.03;
speed_threshold = 0.05;
% joint_center_speed_trajectory = sum(joint_center_velocity_trajectory_ankle.^2, 2).^0.5;
% joint_center_velocity_trajectory_fast = joint_center_velocity_trajectory_ankle(joint_center_speed_trajectory > speed_threshold, :);
% joint_center_velocity_trajectory_fast_normed = joint_center_velocity_trajectory_fast .* repmat(joint_center_speed_trajectory(joint_center_speed_trajectory > speed_threshold), 1, 3).^(-1);
% joint_center_velocity_trajectory_fast_normed_flipped = joint_center_velocity_trajectory_fast_normed;
% joint_center_velocity_trajectory_fast_normed_flipped(joint_center_velocity_trajectory_fast_normed(:, 1)<0, :) = -joint_center_velocity_trajectory_fast_normed_flipped(joint_center_velocity_trajectory_fast_normed(:, 1)<0, :);
% ankle_elevation_trajectory_fast = ankle_elevation_trajectory(joint_center_speed_trajectory > speed_threshold);
% number_of_time_frames_fast = length(ankle_elevation_trajectory_fast);

% replace the time steps with too small speed with NaNs
joint_center_speed_trajectory = sum(joint_center_velocity_trajectory_ankle.^2, 2).^0.5;
joint_center_velocity_trajectory_fast = joint_center_velocity_trajectory_ankle;
joint_center_velocity_trajectory_fast(joint_center_speed_trajectory < speed_threshold, :) = NaN;
joint_center_velocity_trajectory_fast_normed = joint_center_velocity_trajectory_fast .* repmat(joint_center_speed_trajectory, 1, 3).^(-1);
joint_center_velocity_trajectory_fast_normed_flipped = joint_center_velocity_trajectory_fast_normed;
joint_center_velocity_trajectory_fast_normed_flipped(joint_center_velocity_trajectory_fast_normed(:, 1)<0, :) = -joint_center_velocity_trajectory_fast_normed_flipped(joint_center_velocity_trajectory_fast_normed(:, 1)<0, :);

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







% project ankle CoR onto the plane of movement and make polynomial fit
[~, ~, V] = svd(plane_normal');
V_plane = V(:, 2:3);
joint_center_trajectory_transformed = (V * joint_center_trajectory_world')';
polyfit_y = polyfit(foot_elevation_trajectory, joint_center_trajectory_transformed(:, 2), 5);
polyfit_z = polyfit(foot_elevation_trajectory, joint_center_trajectory_transformed(:, 3), 5);
polyfit_y_derivative = polyder(polyfit_y);
polyfit_z_derivative = polyder(polyfit_z);

% evaluate fits
elevation_values = linspace(min(foot_elevation_trajectory), max(foot_elevation_trajectory), 100);
polyfit_y_values = polyval(polyfit_y, elevation_values);
polyfit_y_derivatives = polyval(polyfit_y_derivative, elevation_values);
polyfit_z_values = polyval(polyfit_z, elevation_values);
polyfit_z_derivatives = polyval(polyfit_z_derivative, elevation_values);

% plot
% % 3d paths
% figure; axes; hold on; axis equal;
% xlabel('x'); ylabel('y'), zlabel('z');
% plot3(joint_center_trajectory_world(:, 1), joint_center_trajectory_world(:, 2), joint_center_trajectory_world(:, 3));
% for i_marker = 1 : number_of_markers
%     plot3(right_foot_markers_trajectory(:, (i_marker-1)*3+1), right_foot_markers_trajectory(:, (i_marker-1)*3+2), right_foot_markers_trajectory(:, (i_marker-1)*3+3));
% end
% for i_marker = 2 : number_of_markers
%     plot3(circle_reconstr{i_marker}(1, :), circle_reconstr{i_marker}(2, :), circle_reconstr{i_marker}(3, :))
% end

% 3d paths
% figure; axes; hold on; axis equal; title('ankle CoR path')
% xlabel('x'); ylabel('y'), zlabel('z');
% plot3(joint_center_trajectory_transformed(:, 1), joint_center_trajectory_transformed(:, 2), joint_center_trajectory_transformed(:, 3));

% 3d velocities
% figure; axes; hold on; axis equal; title('velocity, world');
% xlabel('x'); ylabel('y'), zlabel('z');
% plot3(joint_center_velocity_trajectory_world(:, 1), joint_center_velocity_trajectory_world(:, 2), joint_center_velocity_trajectory_world(:, 3), 'x');

% figure; axes; hold on; axis equal; title('velocity, ankle');
% xlabel('x'); ylabel('y'), zlabel('z');
% huedata = (ankle_elevation_trajectory - min(ankle_elevation_trajectory)) * 1 / (max(ankle_elevation_trajectory) - min(ankle_elevation_trajectory));
% for i_time = 1 : number_of_time_frames
%     if ~isnan(joint_center_velocity_trajectory_ankle(i_time, 1))
%         color = hsv2rgb(huedata(i_time), 1, 1);
%         plot3(joint_center_velocity_trajectory_ankle(i_time, 1), joint_center_velocity_trajectory_ankle(i_time, 2), joint_center_velocity_trajectory_ankle(i_time, 3), 'x', 'color', color);
%     end    
% end
% % plot3(joint_center_velocity_trajectory_ankle(:, 1), joint_center_velocity_trajectory_ankle(:, 2), joint_center_velocity_trajectory_ankle(:, 3), 'x');
% % xdata = joint_center_velocity_trajectory(:, 1)';
% % ydata = joint_center_velocity_trajectory(:, 2)';
% % zdata = joint_center_velocity_trajectory(:, 3)';
% % colordata = ankle_elevation_trajectory';
% % surface([xdata; xdata],[ydata; ydata],[zdata; zdata],[colordata; colordata], 'facecol', 'no', 'edgecol', 'interp', 'linew', 2);
% ellipsoid(0,0,0,speed_threshold,speed_threshold,speed_threshold,20)
% alpha(0.3)

% figure; axes; hold on; axis equal; title('velocity, ankle');
% xlabel('x'); ylabel('y'), zlabel('z');
% % xdata = joint_center_velocity_trajectory_fast(:, 1)';
% % ydata = joint_center_velocity_trajectory_fast(:, 2)';
% % zdata = joint_center_velocity_trajectory_fast(:, 3)';
% % colordata = ankle_elevation_trajectory_fast';
% % surface([xdata; xdata],[ydata; ydata],[zdata; zdata],[colordata; colordata], 'facecol', 'no', 'edgecol', 'interp', 'linew', 2);
% huedata = (ankle_elevation_trajectory_fast - min(ankle_elevation_trajectory_fast)) * 1 / (max(ankle_elevation_trajectory_fast) - min(ankle_elevation_trajectory_fast));
% for i_time = 1 : number_of_time_frames_fast
%     if ~isnan(joint_center_velocity_trajectory_fast(i_time, 1))
%         color = hsv2rgb(huedata(i_time), 1, 1);
%         plot3(joint_center_velocity_trajectory_fast(i_time, 1), joint_center_velocity_trajectory_fast(i_time, 2), joint_center_velocity_trajectory_fast(i_time, 3), 'x', 'color', color);
%     end    
% end
% ellipsoid(0,0,0,speed_threshold,speed_threshold,speed_threshold,20)
% alpha(0.3)

% figure; axes; hold on; axis equal; title('velocity, ankle');
% xlabel('x'); ylabel('y'), zlabel('z');
% xdata = joint_center_velocity_trajectory_fast(:, 1)';
% ydata = joint_center_velocity_trajectory_fast(:, 2)';
% zdata = joint_center_velocity_trajectory_fast(:, 3)';
% colordata = ankle_elevation_trajectory';
% surface([xdata; xdata],[ydata; ydata],[zdata; zdata],[colordata; colordata], 'facecol', 'no', 'edgecol', 'interp', 'linew', 2);
% ellipsoid(0,0,0,speed_threshold,speed_threshold,speed_threshold,20)
% alpha(0.3)
% 
% figure; axes; hold on; axis equal; title('velocity, ankle');
% xlabel('x'); ylabel('y'), zlabel('z');
% xdata = joint_center_velocity_trajectory_fast_normed(:, 1)';
% ydata = joint_center_velocity_trajectory_fast_normed(:, 2)';
% zdata = joint_center_velocity_trajectory_fast_normed(:, 3)';
% colordata = ankle_elevation_trajectory';
% surface([xdata; xdata],[ydata; ydata],[zdata; zdata],[colordata; colordata], 'facecol', 'no', 'edgecol', 'interp', 'linew', 2);
% ellipsoid(0,0,0,1,1,1,20)
% alpha(0.3)

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




% [joint_center_trajectory, parameter_angle_trajectory, joint_center_trajectory_transformed] = estimateConstraintSurface ...
%   ( ...
%     right_ankle_cor_rightfoot, ...
%     right_foot_markers_reference, ...
%     right_foot_markers_trajectory, ...
%     [1, 4] ...
%   );

% figure; axes; hold on; xlabel('elevation, foot'); ylabel('ankle CoR')
% plot(foot_elevation_trajectory, joint_center_trajectory_world(:, 1));
% figure; axes; hold on; xlabel('elevation, foot'); ylabel('ankle CoR')
% plot(foot_elevation_trajectory, joint_center_trajectory_world(:, 2));
% figure; axes; hold on; xlabel('elevation, foot'); ylabel('ankle CoR')
% plot(foot_elevation_trajectory, joint_center_trajectory_world(:, 3));
% 
% figure; axes; hold on; xlabel('elevation, ankle'); ylabel('ankle CoR')
% plot(ankle_elevation_trajectory, joint_center_trajectory_world(:, 1));
% figure; axes; hold on; xlabel('elevation, ankle'); ylabel('ankle CoR')
% plot(ankle_elevation_trajectory, joint_center_trajectory_world(:, 2));
% figure; axes; hold on; xlabel('elevation, ankle'); ylabel('ankle CoR')
% plot(ankle_elevation_trajectory, joint_center_trajectory_world(:, 3));

% figure; axes; hold on; xlabel('elevation, foot'); ylabel('ankle CoR velocity, world')
% plot(foot_elevation_trajectory, joint_center_velocity_trajectory_world(:, 1));
% plot(foot_elevation_trajectory, joint_center_velocity_trajectory_world(:, 2));
% plot(foot_elevation_trajectory, joint_center_velocity_trajectory_world(:, 3));
% 
% figure; axes; hold on; xlabel('elevation, foot'); ylabel('ankle CoR velocity, ankle')
% plot(foot_elevation_trajectory, joint_center_velocity_trajectory_ankle(:, 1));
% plot(foot_elevation_trajectory, joint_center_velocity_trajectory_ankle(:, 2));
% plot(foot_elevation_trajectory, joint_center_velocity_trajectory_ankle(:, 3));






% % find right toes constraint surface
% right_foot_markers_trajectory = right_toes_reference(:, right_foot_markers_indices);
% [joint_center_trajectory, parameter_angle_trajectory] = estimateConstraintSurface ...
%   ( ...
%     right_ankle_cor_rightfoot, ...
%     right_foot_markers_reference, ...
%     right_foot_markers_trajectory, ...
%     [1, 4] ...
%   );
trajectory_x = joint_center_trajectory_transformed(:, 1);
trajectory_y = joint_center_trajectory_transformed(:, 2);
trajectory_z = joint_center_trajectory_transformed(:, 3);

