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

% optimizes the joint angles by minimizing the squared distance between measured and reconstructed marker positions

function optimizedJointAngles = optimizeJointAngles ...
  ( ...
    complete_tree, ...
    pelvis_chain, ...
    left_leg_chain, ...
    right_leg_chain, ...
    torso_chain, ...
    left_arm_chain, ...
    right_arm_chain, ...
    markerTrajectories, ...
    joint_angle_trajectories, ...
    weight_matrix ...
  )

    number_of_time_steps = size(markerTrajectories, 1);
    number_of_joints = complete_tree.numberOfJoints;
    number_of_markers = size(markerTrajectories, 2) / 3;

    current_marker_positions_reconstr = [];
    optimizedJointAngles = zeros(number_of_time_steps, number_of_joints);

    pelvis_joints = complete_tree.getJointGroup('pelvis');
    left_leg_joints = complete_tree.getJointGroup('left leg');
    right_leg_joints = complete_tree.getJointGroup('right leg');
    torso_joints = complete_tree.getJointGroup('torso');
    left_arm_joints = complete_tree.getJointGroup('left arm');
    right_arm_joints = complete_tree.getJointGroup('right arm');
    left_hip_joints = complete_tree.getJointGroup('left hip');
    left_knee_joints = complete_tree.getJointGroup('left knee');
    left_ankle_joints = complete_tree.getJointGroup('left ankle');
    right_hip_joints = complete_tree.getJointGroup('right hip');
    right_knee_joints = complete_tree.getJointGroup('right knee');
    right_ankle_joints = complete_tree.getJointGroup('right ankle');
    lumbar_joints = complete_tree.getJointGroup('lumbar');
    cervix_joints = complete_tree.getJointGroup('cervix');
    left_shoulder_joints = complete_tree.getJointGroup('left shoulder');
    left_elbow_joints = complete_tree.getJointGroup('left elbow');
    left_wrist_joints = complete_tree.getJointGroup('left wrist');
    right_shoulder_joints = complete_tree.getJointGroup('right shoulder');
    right_elbow_joints = complete_tree.getJointGroup('right elbow');
    right_wrist_joints = complete_tree.getJointGroup('right wrist');

    marker_map = complete_tree.markerExportMap;
    pelvis_markers = find(marker_map(1, :)==pelvis_joints(end));
    left_thigh_markers = find(marker_map(1, :)==left_hip_joints(end));
    left_shank_markers = find(marker_map(1, :)==left_knee_joints(end));
    left_foot_markers = find(marker_map(1, :)==left_ankle_joints(end));
    left_leg_markers = [left_thigh_markers left_shank_markers left_foot_markers];
    right_thigh_markers = find(marker_map(1, :)==right_hip_joints(end));
    right_shank_markers = find(marker_map(1, :)==right_knee_joints(end));
    right_foot_markers = find(marker_map(1, :)==right_ankle_joints(end));
    right_leg_markers = [right_thigh_markers right_shank_markers right_foot_markers];
    trunk_markers = find(marker_map(1, :)==lumbar_joints(end));
    left_upperarm_markers = find(marker_map(1, :)==left_shoulder_joints(end));
    left_lowerarm_markers = find(marker_map(1, :)==left_elbow_joints(end));
    left_hand_markers = find(marker_map(1, :)==left_wrist_joints(end));
    left_arm_markers = [left_upperarm_markers left_lowerarm_markers left_hand_markers];
    right_upperarm_markers = find(marker_map(1, :)==right_shoulder_joints(end));
    right_lowerarm_markers = find(marker_map(1, :)==right_elbow_joints(end));
    right_hand_markers = find(marker_map(1, :)==right_wrist_joints(end));
    right_arm_markers = [right_upperarm_markers right_lowerarm_markers right_hand_markers];
    head_markers = find(marker_map(1, :)==cervix_joints(end));
    trunk_and_head_markers = [trunk_markers head_markers];

    last_left_hip_joint_index_in_limb_plant = left_hip_joints(end) - (left_leg_joints(1) - 1);
    last_left_knee_joint_index_in_limb_plant = left_knee_joints(end) - (left_leg_joints(1) - 1);
    last_left_ankle_joint_index_in_limb_plant = left_ankle_joints(end) - (left_leg_joints(1) - 1);
    last_right_hip_joint_index_in_limb_plant = right_hip_joints(end) - (right_leg_joints(1) - 1);
    last_right_knee_joint_index_in_limb_plant = right_knee_joints(end) - (right_leg_joints(1) - 1);
    last_right_ankle_joint_index_in_limb_plant = right_ankle_joints(end) - (right_leg_joints(1) - 1);
    last_lumbar_joint_index_in_limb_plant = lumbar_joints(end) - (torso_joints(1) - 1);
    last_cervix_joint_index_in_limb_plant = cervix_joints(end) - (torso_joints(1) - 1);
    last_left_shoulder_joint_index_in_limb_plant = left_shoulder_joints(end) - (left_arm_joints(1) - 1);
    last_left_elbow_joint_index_in_limb_plant = left_elbow_joints(end) - (left_arm_joints(1) - 1);
    last_left_wrist_joint_index_in_limb_plant = left_wrist_joints(end) - (left_arm_joints(1) - 1);
    last_right_shoulder_joint_index_in_limb_plant = right_shoulder_joints(end) - (right_arm_joints(1) - 1);
    last_right_elbow_joint_index_in_limb_plant = right_elbow_joints(end) - (right_arm_joints(1) - 1);
    last_right_wrist_joint_index_in_limb_plant = right_wrist_joints(end) - (right_arm_joints(1) - 1);
    
    marker_map_pelvis = pelvis_chain.markerExportMap;
    marker_map_left_leg = left_leg_chain.markerExportMap;
    marker_map_right_leg = right_leg_chain.markerExportMap;
    marker_map_trunk_and_head = torso_chain.markerExportMap;
    marker_map_left_arm = left_arm_chain.markerExportMap;
    marker_map_right_arm = right_arm_chain.markerExportMap;
    pelvis_markers_modular = find(marker_map_pelvis(1, :)==pelvis_joints(end));
    left_leg_markers_modular = [find(marker_map_left_leg(1, :)==last_left_hip_joint_index_in_limb_plant) find(marker_map_left_leg(1, :)==last_left_knee_joint_index_in_limb_plant) find(marker_map_left_leg(1, :)==last_left_ankle_joint_index_in_limb_plant)];
    right_leg_markers_modular = [find(marker_map_right_leg(1, :)==last_right_hip_joint_index_in_limb_plant) find(marker_map_right_leg(1, :)==last_right_knee_joint_index_in_limb_plant) find(marker_map_right_leg(1, :)==last_right_ankle_joint_index_in_limb_plant)];
    trunk_and_head_markers_modular = [find(marker_map_trunk_and_head(1, :)==last_lumbar_joint_index_in_limb_plant) find(marker_map_trunk_and_head(1, :)==last_cervix_joint_index_in_limb_plant)];
    left_arm_markers_modular = [find(marker_map_left_arm(1, :)==last_left_shoulder_joint_index_in_limb_plant) find(marker_map_left_arm(1, :)==last_left_elbow_joint_index_in_limb_plant) find(marker_map_left_arm(1, :)==last_left_wrist_joint_index_in_limb_plant)];
    right_arm_markers_modular = [find(marker_map_right_arm(1, :)==last_right_shoulder_joint_index_in_limb_plant) find(marker_map_right_arm(1, :)==last_right_elbow_joint_index_in_limb_plant) find(marker_map_right_arm(1, :)==last_right_wrist_joint_index_in_limb_plant)];

    % trunk_and_head_markers(6 : 19) = []; % Achtung! hard-coded to remove the arm markers
    % trunk_and_head_markers_modular(6 : 19) = [];

    options = optimset ...
        ( ...
            'GradObj', 'off', ...
            'Display','off', ...
            'LargeScale', 'off', ...
            'DerivativeCheck', 'on', ...
            'UseParallel', 'always' ...
        );

    weight_matrix_by_indices = reshape(repmat(weight_matrix, 3, 1), 1, length(weight_matrix)*3);

    for i_time = 1 : number_of_time_steps
        theta_0 = joint_angle_trajectories(i_time, :)';
%         theta_0 = zeros(complete_tree.numberOfJoints, 1);
        current_marker_positions_measured = markerTrajectories(i_time, :);
        current_marker_positions_relevant = current_marker_positions_measured(weight_matrix_by_indices~=0);

        if any(isnan(current_marker_positions_relevant)) % check whether NaNs are present
            theta_opt = zeros(size(theta_0)) * NaN;
        else % if not, optimize

            % use modular plant
            theta_virtual = fminunc(@objfun_virtual_modular, theta_0(pelvis_joints), options);

            pelvis_chain.jointAngles = theta_virtual;
            pelvis_chain.updateConfiguration();
            pelvis_to_world_poe = pelvis_chain.productsOfExponentials{6};
            theta_left_leg = fminunc(@objfun_left_leg_modular, theta_0(left_leg_joints), options);
            theta_right_leg = fminunc(@objfun_right_leg_modular, theta_0(right_leg_joints), options);
            theta_trunk_and_head = fminunc(@objfun_trunk_modular, theta_0(torso_joints), options);

            trunk_to_world_poe = pelvis_chain.productsOfExponentials{6} * torso_chain.productsOfExponentials{3};
            theta_left_arm = fminunc(@objfun_left_arm_modular, theta_0(left_arm_joints), options);
            theta_right_arm = fminunc(@objfun_right_arm_modular, theta_0(right_arm_joints), options);

            theta_opt = [theta_virtual' theta_left_leg' theta_right_leg' theta_trunk_and_head' theta_left_arm' theta_right_arm'];

        end
        optimizedJointAngles(i_time, :) = theta_opt;
    end

    function f = objfun_virtual_modular(theta_virtual)
        pelvis_chain.jointAngles = theta_virtual;
        pelvis_chain.updateKinematics();
        current_marker_positions_reconstr = pelvis_chain.exportMarkerPositions();

        % calculate reconstruction error
        marker_reconstruction_error = zeros(1, number_of_markers);
        for i_marker = 1 : length(pelvis_markers_modular)
            marker_indices = 3*(pelvis_markers(i_marker)-1)+1 : 3*(pelvis_markers(i_marker)-1)+3;
            marker_indices_modular = 3*(pelvis_markers_modular(i_marker)-1)+1 : 3*(pelvis_markers_modular(i_marker)-1)+3;

            marker_position_measured = current_marker_positions_measured(marker_indices);
            marker_position_reconstr = current_marker_positions_reconstr(marker_indices_modular);

            error_vector = marker_position_measured - marker_position_reconstr;
            marker_reconstruction_error(i_marker) = norm(error_vector * weight_matrix(pelvis_markers(i_marker)));
        end

        % calculate sum of squares
        f = sum(marker_reconstruction_error.^2);
    end

    function f = objfun_right_leg_modular(theta_right_leg)
        right_leg_chain.jointAngles = theta_right_leg;
        right_leg_chain.updateKinematics();
        current_marker_positions_reconstr = right_leg_chain.exportMarkerPositions();

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
            marker_reconstruction_error(i_marker) = norm(error_vector * weight_matrix(right_leg_markers(i_marker)));
        end

        % calculate sum of squares
        f = sum(marker_reconstruction_error.^2);
    end

    function f = objfun_left_leg_modular(theta_left_leg)
        left_leg_chain.jointAngles = theta_left_leg;
        left_leg_chain.updateKinematics();
        current_marker_positions_reconstr = left_leg_chain.exportMarkerPositions();

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
            marker_reconstruction_error(i_marker) = norm(error_vector * weight_matrix(left_leg_markers(i_marker)));
        end

        % calculate sum of squares
        f = sum(marker_reconstruction_error.^2);
    end

    function f = objfun_trunk_modular(theta_trunk)
        torso_chain.jointAngles = theta_trunk;
        torso_chain.updateKinematics();
        current_marker_positions_reconstr = torso_chain.exportMarkerPositions();

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
            marker_reconstruction_error(i_marker) = norm(error_vector * weight_matrix(trunk_and_head_markers(i_marker)));
        end

        % calculate sum of squares
        f = sum(marker_reconstruction_error.^2);
    end

    function f = objfun_right_arm_modular(theta_right_arm)
        right_arm_chain.jointAngles = theta_right_arm;
        right_arm_chain.updateKinematics();
        current_marker_positions_reconstr = right_arm_chain.exportMarkerPositions();

        % calculate reconstruction error
        marker_reconstruction_error = zeros(1, number_of_markers);
        for i_marker = 1 : length(right_arm_markers_modular)
            marker_indices = 3*(right_arm_markers(i_marker)-1)+1 : 3*(right_arm_markers(i_marker)-1)+3;
            marker_indices_modular = 3*(right_arm_markers_modular(i_marker)-1)+1 : 3*(right_arm_markers_modular(i_marker)-1)+3;

            marker_position_measured = current_marker_positions_measured(marker_indices);
            marker_position_reconstr_pelvis = current_marker_positions_reconstr(marker_indices_modular);

            marker_position_reconstr_world_homogeneous = trunk_to_world_poe * [marker_position_reconstr_pelvis'; 1];
            marker_position_reconstr_world = marker_position_reconstr_world_homogeneous(1:3)';

            error_vector = marker_position_measured - marker_position_reconstr_world;
            marker_reconstruction_error(i_marker) = norm(error_vector * weight_matrix(right_arm_markers(i_marker)));
        end

        % calculate sum of squares
        f = sum(marker_reconstruction_error.^2);
    end

    function f = objfun_left_arm_modular(theta_left_arm)
        left_arm_chain.jointAngles = theta_left_arm;
        left_arm_chain.updateKinematics();
        current_marker_positions_reconstr = left_arm_chain.exportMarkerPositions();

        % calculate reconstruction error
        marker_reconstruction_error = zeros(1, number_of_markers);
        for i_marker = 1 : length(left_arm_markers_modular)
            marker_indices = 3*(left_arm_markers(i_marker)-1)+1 : 3*(left_arm_markers(i_marker)-1)+3;
            marker_indices_modular = 3*(left_arm_markers_modular(i_marker)-1)+1 : 3*(left_arm_markers_modular(i_marker)-1)+3;

            marker_position_measured = current_marker_positions_measured(marker_indices);
            marker_position_reconstr_pelvis = current_marker_positions_reconstr(marker_indices_modular);

            marker_position_reconstr_world_homogeneous = trunk_to_world_poe * [marker_position_reconstr_pelvis'; 1];
            marker_position_reconstr_world = marker_position_reconstr_world_homogeneous(1:3)';

            error_vector = marker_position_measured - marker_position_reconstr_world;
            marker_reconstruction_error(i_marker) = norm(error_vector * weight_matrix(left_arm_markers(i_marker)));
        end

        % calculate sum of squares
        f = sum(marker_reconstruction_error.^2);
    end

end
