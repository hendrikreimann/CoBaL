function optimizedJointAngles = optimizeJointAngles ...
( ...
  plant, ...
  pelvisChain, ...
  leftLegChain, ...
  rightLegChain, ...
  trunkChain, ...
  markerTrajectories ...
)

%
% ---------------------------------------------------------------------
%
% optimizeJointAngles.m
%
% optimizes the joint angles by minimizing the squared distance between measured and reconstructed marker positions
% 
% Hendrik Reimann, 2016
% Department of Kinesiology
% Temple University
% hendrikreimann@gmail.com
%
% ---------------------------------------------------------------------
%
%
% ---------------------------------------------------------------------

number_of_time_steps = size(markerTrajectories, 1);
number_of_joints = plant.numberOfJoints;
number_of_markers = size(markerTrajectories, 2) / 3;

current_marker_positions_reconstr = [];
optimizedJointAngles = zeros(number_of_time_steps, number_of_joints);
virtual_joints = 1 : 6;
left_leg_joints = 7 : 12;
right_leg_joints = 13 : 18;
trunk_and_neck_joints = 19 : 24;
marker_map = plant.markerExportMap;
virtual_joint_markers = find(marker_map(1, :)==6);
left_thigh_markers = find(marker_map(1, :)==9);
left_shank_markers = find(marker_map(1, :)==10);
left_foot_markers = find(marker_map(1, :)==12);
left_leg_markers = [left_thigh_markers left_shank_markers left_foot_markers];
right_thigh_markers = find(marker_map(1, :)==15);
right_shank_markers = find(marker_map(1, :)==16);
right_foot_markers = find(marker_map(1, :)==18);
right_leg_markers = [right_thigh_markers right_shank_markers right_foot_markers];
trunk_markers = find(marker_map(1, :)==21);

head_markers = find(marker_map(1, :)==24);
trunk_and_head_markers = [trunk_markers head_markers];

marker_map_pelvis = pelvisChain.markerExportMap;
marker_map_right_leg = rightLegChain.markerExportMap;
marker_map_left_leg = leftLegChain.markerExportMap;
marker_map_trunk = trunkChain.markerExportMap;
virtual_joint_markers_modular = find(marker_map_pelvis(1, :)==6);
left_leg_markers_modular = [find(marker_map_left_leg(1, :)==3) find(marker_map_left_leg(1, :)==4) find(marker_map_left_leg(1, :)==6)];
right_leg_markers_modular = [find(marker_map_right_leg(1, :)==3) find(marker_map_right_leg(1, :)==4) find(marker_map_right_leg(1, :)==6)];
trunk_and_head_markers_modular = [find(marker_map_trunk(1, :)==3) find(marker_map_trunk(1, :)==6)];

trunk_and_head_markers(6 : 19) = []; % Achtung! hard-coded to remove the arm markers
trunk_and_head_markers_modular(6 : 19) = [];

options = optimset ...
    ( ...
        'GradObj', 'off', ...
        'Display','off', ...
        'LargeScale', 'off', ...
        'DerivativeCheck', 'on', ...
        'UseParallel', 'always' ...
    );

for i_time = 1 : number_of_time_steps
    theta_0 = zeros(plant.numberOfJoints, 1);
    current_marker_positions_measured = markerTrajectories(i_time, :);
    if any(isnan(current_marker_positions_measured)) % check whether NaNs are present
        theta_opt = zeros(size(theta_0)) * NaN;
    else % if not, optimize
        
        % use modular plant
        theta_virtual = fminunc(@objfun_virtual_modular, theta_0(virtual_joints), options);
        
        pelvisChain.jointAngles = theta_virtual;
        pelvisChain.updateConfiguration();
        pelvis_to_world_poe = pelvisChain.productsOfExponentials{6};
        theta_left_leg = fminunc(@objfun_left_leg_modular, theta_0(left_leg_joints), options);
        theta_right_leg = fminunc(@objfun_right_leg_modular, theta_0(right_leg_joints), options);
        theta_trunk = fminunc(@objfun_trunk_modular, theta_0(trunk_and_neck_joints), options);
        theta_opt = [theta_virtual; theta_left_leg; theta_right_leg; theta_trunk];
        
    end
    optimizedJointAngles(i_time, :) = theta_opt;
end


function f = objfun_virtual_modular(theta_virtual)
    pelvisChain.jointAngles = theta_virtual;
    pelvisChain.updateKinematics();
    current_marker_positions_reconstr = pelvisChain.exportMarkerPositions();

    % calculate reconstruction error
    marker_reconstruction_error = zeros(1, number_of_markers);
    for i_marker = 1 : length(virtual_joint_markers_modular)
        marker_indices = 3*(virtual_joint_markers(i_marker)-1)+1 : 3*(virtual_joint_markers(i_marker)-1)+3;
        marker_indices_modular = 3*(virtual_joint_markers_modular(i_marker)-1)+1 : 3*(virtual_joint_markers_modular(i_marker)-1)+3;
        
        marker_position_measured = current_marker_positions_measured(marker_indices);
        marker_position_reconstr = current_marker_positions_reconstr(marker_indices_modular);
        
        error_vector = marker_position_measured - marker_position_reconstr;
        marker_reconstruction_error(i_marker) = norm(error_vector);
    end
    
    % calculate sum of squares
    f = sum(marker_reconstruction_error.^2);
end

function f = objfun_right_leg_modular(theta_right_leg)
    rightLegChain.jointAngles = theta_right_leg;
    rightLegChain.updateKinematics();
    current_marker_positions_reconstr = rightLegChain.exportMarkerPositions();

    % calculate reconstruction error
    marker_reconstruction_error = zeros(1, number_of_markers);
    for i_marker = 1 : length(right_leg_markers_modular)
        marker_indices = 3*(right_leg_markers(i_marker)-1)+1 : 3*(right_leg_markers(i_marker)-1)+3;
        marker_indices_modular = 3*(right_leg_markers_modular(i_marker)-1)+1 : 3*(right_leg_markers_modular(i_marker)-1)+3;
        
        marker_position_measured = current_marker_positions_measured(marker_indices);
        marker_position_reconstr_pelvis = current_marker_positions_reconstr(marker_indices_modular);
        
        marker_position_reconstr_world_homogeneous = pelvis_to_world_poe * [marker_position_reconstr_pelvis'; 1];
        marker_position_reconstr_world = marker_position_reconstr_world_homogeneous(1:3)';
        
        error_vector = marker_position_measured - marker_position_reconstr_world;
        marker_reconstruction_error(i_marker) = norm(error_vector);
    end
    
    % calculate sum of squares
    f = sum(marker_reconstruction_error.^2);
end

function f = objfun_left_leg_modular(theta_left_leg)
    leftLegChain.jointAngles = theta_left_leg;
    leftLegChain.updateKinematics();
    current_marker_positions_reconstr = leftLegChain.exportMarkerPositions();

    % calculate reconstruction error
    marker_reconstruction_error = zeros(1, number_of_markers);
    for i_marker = 1 : length(left_leg_markers_modular)
        marker_indices = 3*(left_leg_markers(i_marker)-1)+1 : 3*(left_leg_markers(i_marker)-1)+3;
        marker_indices_modular = 3*(left_leg_markers_modular(i_marker)-1)+1 : 3*(left_leg_markers_modular(i_marker)-1)+3;
        
        marker_position_measured = current_marker_positions_measured(marker_indices);
        marker_position_reconstr_pelvis = current_marker_positions_reconstr(marker_indices_modular);
        
        marker_position_reconstr_world_homogeneous = pelvis_to_world_poe * [marker_position_reconstr_pelvis'; 1];
        marker_position_reconstr_world = marker_position_reconstr_world_homogeneous(1:3)';
        
        error_vector = marker_position_measured - marker_position_reconstr_world;
        marker_reconstruction_error(i_marker) = norm(error_vector);
    end
    
    % calculate sum of squares
    f = sum(marker_reconstruction_error.^2);
end

function f = objfun_trunk_modular(theta_trunk)
    trunkChain.jointAngles = theta_trunk;
    trunkChain.updateKinematics();
    current_marker_positions_reconstr = trunkChain.exportMarkerPositions();

    % calculate reconstruction error
    marker_reconstruction_error = zeros(1, number_of_markers);
    for i_marker = 1 : length(trunk_and_head_markers_modular)
        marker_indices = 3*(trunk_and_head_markers(i_marker)-1)+1 : 3*(trunk_and_head_markers(i_marker)-1)+3;
        marker_indices_modular = 3*(trunk_and_head_markers_modular(i_marker)-1)+1 : 3*(trunk_and_head_markers_modular(i_marker)-1)+3;
        
        marker_position_measured = current_marker_positions_measured(marker_indices);
        marker_position_reconstr_pelvis = current_marker_positions_reconstr(marker_indices_modular);
        
        marker_position_reconstr_world_homogeneous = pelvis_to_world_poe * [marker_position_reconstr_pelvis'; 1];
        marker_position_reconstr_world = marker_position_reconstr_world_homogeneous(1:3)';
        
        error_vector = marker_position_measured - marker_position_reconstr_world;
        marker_reconstruction_error(i_marker) = norm(error_vector);
    end
    
    % calculate sum of squares
    f = sum(marker_reconstruction_error.^2);
end

end
