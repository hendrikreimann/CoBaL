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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% this function optimizes the kinematic variables to better fit the marker data

% input: 
% subjectInfo.mat
% subjectModel.mat
% markerTrajectories
% kinematicTrajectories
%
% output:
% file kinematicTrajectories.mat, containing
% - joint_center_trajectories
% - com_trajectories
% - com_labels
% - joint_angle_trajectories


function optimizeKinematicTrajectories(varargin)
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    
    % load settings
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    
    load('subjectInfo.mat', 'date', 'subject_id');
    load('subjectModel.mat');
    
    number_of_joints = kinematic_tree.numberOfJoints; %#ok<NODEF>
    
    %% determine weights
    weight_matrix = ones(1, length(kinematic_tree.markerLabels)); % TODO: make this a setting
    marker_weight_table = study_settings.get('marker_weights');
    for i_weight = 1 : size(marker_weight_table, 1)
        this_weight_marker_label = marker_weight_table{i_weight, 1};
        this_weight = str2double(marker_weight_table{i_weight, 2});
        
        this_weight_marker_index = find(strcmp(kinematic_tree.markerLabels, this_weight_marker_label), 1, 'first');
        weight_matrix(this_weight_marker_index) = this_weight;
        
    end
    
    %% create limb plants

    % extract reference positions
    kinematic_tree.jointAngles = zeros(kinematic_tree.numberOfJoints, 1);
    kinematic_tree.updateConfiguration();
    joint_positions = cell(1, number_of_joints);
    joint_axes = cell(1, number_of_joints);
    link_com_positions = cell(1, number_of_joints);
    link_orientations = cell(1, number_of_joints);
    generalized_link_inertia_matrices = kinematic_tree.generalizedInertiaMatrices;
    for i_joint = 1 : number_of_joints
        joint_positions{i_joint} = kinematic_tree.jointTransformations{i_joint}(1:3, 4);
        if norm(kinematic_tree.referenceJointTwists{i_joint}(4:6)) ~= 0
            joint_axes{i_joint} = kinematic_tree.referenceJointTwists{i_joint}(4:6);
        else
            joint_axes{i_joint} = kinematic_tree.referenceJointTwists{i_joint}(1:3);
        end
        link_com_positions{i_joint} = kinematic_tree.linkTransformations{i_joint}(1:3, 4);
        link_orientations{i_joint} = kinematic_tree.linkTransformations{i_joint}(1:3, 1:3);
    end
    
    % extract relevant data
    pelvis_joints = kinematic_tree.getJointGroup('pelvis');
    left_leg_joints = kinematic_tree.getJointGroup('left leg');
    right_leg_joints = kinematic_tree.getJointGroup('right leg');
    torso_joints = kinematic_tree.getJointGroup('torso');
    left_arm_joints = kinematic_tree.getJointGroup('left arm');
    right_arm_joints = kinematic_tree.getJointGroup('right arm');
    
    left_hip_joints = kinematic_tree.getJointGroup('left hip');
    left_knee_joints = kinematic_tree.getJointGroup('left knee');
    left_ankle_joints = kinematic_tree.getJointGroup('left ankle');
    right_hip_joints = kinematic_tree.getJointGroup('right hip');
    right_knee_joints = kinematic_tree.getJointGroup('right knee');
    right_ankle_joints = kinematic_tree.getJointGroup('right ankle');
    lumbar_joints = kinematic_tree.getJointGroup('lumbar');
    cervix_joints = kinematic_tree.getJointGroup('cervix');
    left_shoulder_joints = kinematic_tree.getJointGroup('left shoulder');
    left_elbow_joints = kinematic_tree.getJointGroup('left elbow');
    left_wrist_joints = kinematic_tree.getJointGroup('left wrist');
    right_shoulder_joints = kinematic_tree.getJointGroup('right shoulder');
    right_elbow_joints = kinematic_tree.getJointGroup('right elbow');
    right_wrist_joints = kinematic_tree.getJointGroup('right wrist');

    left_hip_position = joint_positions{left_leg_joints(1)};
    right_hip_position = joint_positions{right_leg_joints(1)};
    lumbar_joint_position = joint_positions{torso_joints(1)};
    
    left_toes_eef_index = kinematic_tree.getEndEffectorIndex('left toes');
    right_toes_eef_index = kinematic_tree.getEndEffectorIndex('right toes');
    head_eef_index = kinematic_tree.getEndEffectorIndex('head');
    left_hand_eef_index = kinematic_tree.getEndEffectorIndex('left hand');
    right_hand_eef_index = kinematic_tree.getEndEffectorIndex('right hand');

    % create pelvis tree
    pelvis_chain = GeneralKinematicTree ...
    ( ...
      joint_positions(pelvis_joints), ...
      joint_axes(pelvis_joints), ...
      [2 2 2 1 1 1], ...                                                        % types
      ones(3, size(pelvis_joints, 2)), ...                                      % branch matrix
      {left_hip_position, right_hip_position, lumbar_joint_position}, ...       % end-effectors
      link_com_positions(pelvis_joints), ...
      link_orientations(pelvis_joints), ...
      generalized_link_inertia_matrices(pelvis_joints) ...
    );

    % create left leg tree
    left_leg_chain = GeneralKinematicTree ...
    ( ...
      joint_positions(left_leg_joints), ...
      joint_axes(left_leg_joints), ...
      ones(size(left_leg_joints)), ...                                          % types
      ones(size(left_leg_joints)), ...                                          % branch matrix
      kinematic_tree.endEffectorPositions(left_toes_eef_index), ...             % end-effector
      link_com_positions(left_leg_joints), ...
      link_orientations(left_leg_joints), ...
      generalized_link_inertia_matrices(left_leg_joints) ...
    );

    % create right leg tree
    right_leg_chain = GeneralKinematicTree ...
    ( ...
      joint_positions(right_leg_joints), ...
      joint_axes(right_leg_joints), ...
      ones(size(right_leg_joints)), ...                                         % types
      ones(size(right_leg_joints)), ...                                         % branch matrix
      kinematic_tree.endEffectorPositions(right_toes_eef_index), ...            % end-effector
      link_com_positions(right_leg_joints), ...
      link_orientations(right_leg_joints), ...
      generalized_link_inertia_matrices(right_leg_joints) ...
    );

    % create torso tree
    torso_chain = GeneralKinematicTree ...
    ( ...
      joint_positions(torso_joints), ...
      joint_axes(torso_joints), ...
      ones(size(torso_joints)), ...                                             % types
      ones(size(torso_joints)), ...                                             % branch matrix
      kinematic_tree.endEffectorPositions(head_eef_index), ...                  % end-effector
      link_com_positions(torso_joints), ...
      link_orientations(torso_joints), ...
      generalized_link_inertia_matrices(torso_joints) ...
    );

    % create left arm tree
    left_arm_chain = GeneralKinematicTree ...
    ( ...
      joint_positions(left_arm_joints), ...
      joint_axes(left_arm_joints), ...
      ones(size(left_arm_joints)), ...                                          % types
      ones(size(left_arm_joints)), ...                                          % branch matrix
      kinematic_tree.endEffectorPositions(left_hand_eef_index), ...             % end-effector
      link_com_positions(left_arm_joints), ...
      link_orientations(left_arm_joints), ...
      generalized_link_inertia_matrices(left_arm_joints) ...
    );

    % create right arm tree
    right_arm_chain = GeneralKinematicTree ...
    ( ...
      joint_positions(right_arm_joints), ...
      joint_axes(right_arm_joints), ...
      ones(size(torso_joints)), ...                                             % types
      ones(size(torso_joints)), ...                                             % branch matrix
      kinematic_tree.endEffectorPositions(right_hand_eef_index), ...            % end-effector
      link_com_positions(right_arm_joints), ...
      link_orientations(right_arm_joints), ...
      generalized_link_inertia_matrices(right_arm_joints) ...
    );

    % add markers to pelvis tree
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{pelvis_joints(end)}, 2)
        pelvis_chain.addMarker(6, kinematic_tree.markerReferencePositions{pelvis_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(pelvis_joints(end), i_marker));
    end

    % add markers to left leg tree
    last_left_hip_joint_index_in_limb_plant = left_hip_joints(end) - (left_leg_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{left_hip_joints(end)}, 2)
        left_leg_chain.addMarker(last_left_hip_joint_index_in_limb_plant, kinematic_tree.markerReferencePositions{left_hip_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(left_hip_joints(end), i_marker));
    end
    last_left_knee_joint_index_in_limb_plant = left_knee_joints(end) - (left_leg_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{left_knee_joints(end)}, 2)
        left_leg_chain.addMarker(last_left_knee_joint_index_in_limb_plant, kinematic_tree.markerReferencePositions{left_knee_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(left_knee_joints(end), i_marker));
    end
    last_left_ankle_joint_index_in_limb_plant = left_ankle_joints(end) - (left_leg_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{left_ankle_joints(end)}, 2)
        left_leg_chain.addMarker(last_left_ankle_joint_index_in_limb_plant, kinematic_tree.markerReferencePositions{left_ankle_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(left_ankle_joints(end), i_marker));
    end
    
    % add markers to right leg tree
    last_right_hip_joint_index_in_limb_plant = right_hip_joints(end) - (right_leg_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{right_hip_joints(end)}, 2)
        right_leg_chain.addMarker(last_right_hip_joint_index_in_limb_plant, kinematic_tree.markerReferencePositions{right_hip_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(right_hip_joints(end), i_marker));
    end
    last_right_knee_joint_index_in_limb_plant = right_knee_joints(end) - (right_leg_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{right_knee_joints(end)}, 2)
        right_leg_chain.addMarker(last_right_knee_joint_index_in_limb_plant, kinematic_tree.markerReferencePositions{right_knee_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(right_knee_joints(end), i_marker));
    end
    last_right_ankle_joint_index_in_limb_plant = right_ankle_joints(end) - (right_leg_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{right_ankle_joints(end)}, 2)
        right_leg_chain.addMarker(last_right_ankle_joint_index_in_limb_plant, kinematic_tree.markerReferencePositions{right_ankle_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(right_ankle_joints(end), i_marker));
    end
    
    % add markers to torso tree
    last_lumbar_joint_index_in_limb_plant = lumbar_joints(end) - (torso_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{lumbar_joints(end)}, 2)
        torso_chain.addMarker(last_lumbar_joint_index_in_limb_plant, kinematic_tree.markerReferencePositions{lumbar_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(lumbar_joints(end), i_marker));
    end
    last_cervix_joint_index_in_limb_plant = cervix_joints(end) - (torso_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{cervix_joints(end)}, 2)
        torso_chain.addMarker(last_cervix_joint_index_in_limb_plant, kinematic_tree.markerReferencePositions{cervix_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(cervix_joints(end), i_marker));
    end
    
    % add markers to left arm tree
    last_left_shoulder_joint_index_in_limb_plant = left_shoulder_joints(end) - (left_arm_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{left_shoulder_joints(end)}, 2)
        left_arm_chain.addMarker(last_left_shoulder_joint_index_in_limb_plant, kinematic_tree.markerReferencePositions{left_shoulder_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(left_shoulder_joints(end), i_marker));
    end
    last_left_elbow_joint_index_in_limb_plant = left_elbow_joints(end) - (left_arm_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{left_elbow_joints(end)}, 2)
        left_arm_chain.addMarker(last_left_elbow_joint_index_in_limb_plant, kinematic_tree.markerReferencePositions{left_elbow_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(left_elbow_joints(end), i_marker));
    end
    last_left_wrist_joint_index_in_limb_plant = left_wrist_joints(end) - (left_arm_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{left_wrist_joints(end)}, 2)
        left_arm_chain.addMarker(last_left_wrist_joint_index_in_limb_plant, kinematic_tree.markerReferencePositions{left_wrist_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(left_wrist_joints(end), i_marker));
    end
    
    % add markers to right arm tree
    last_right_shoulder_joint_index_in_limb_plant = right_shoulder_joints(end) - (right_arm_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{right_shoulder_joints(end)}, 2)
        right_arm_chain.addMarker(last_right_shoulder_joint_index_in_limb_plant, kinematic_tree.markerReferencePositions{right_shoulder_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(right_shoulder_joints(end), i_marker));
    end
    last_right_elbow_joint_index_in_limb_plant = right_elbow_joints(end) - (right_arm_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{right_elbow_joints(end)}, 2)
        right_arm_chain.addMarker(last_right_elbow_joint_index_in_limb_plant, kinematic_tree.markerReferencePositions{right_elbow_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(right_elbow_joints(end), i_marker));
    end
    last_right_wrist_joint_index_in_limb_plant = right_wrist_joints(end) - (right_arm_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{right_wrist_joints(end)}, 2)
        right_arm_chain.addMarker(last_right_wrist_joint_index_in_limb_plant, kinematic_tree.markerReferencePositions{right_wrist_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(right_wrist_joints(end), i_marker));
    end
    
    
    
    
    
    
%     % visualize to check the limb chains
%     pelvis_chain.updateConfiguration;
%     left_leg_chain.updateConfiguration;
%     right_leg_chain.updateConfiguration;
%     torso_chain.updateConfiguration;
%     left_arm_chain.updateConfiguration;
%     right_arm_chain.updateConfiguration;
%     
%     scene_bound = [-0.5 0.5; 0.5 1.5; 0 2];
%     stick_figure = KinematicTreeController(kinematic_tree, scene_bound, 'none');
%     stick_figure_pelvis = KinematicTreeController(pelvis_chain, scene_bound, 'none');
%     stick_figure_left_leg = KinematicTreeController(left_leg_chain, scene_bound, 'none');
%     stick_figure_right_leg = KinematicTreeController(right_leg_chain, scene_bound, 'none');
%     stick_figure_torso = KinematicTreeController(torso_chain, scene_bound, 'none');
%     stick_figure_left_arm = KinematicTreeController(left_arm_chain, scene_bound, 'none');
%     stick_figure_right_arm = KinematicTreeController(right_arm_chain, scene_bound, 'none');
%     
%     Link = linkprop([stick_figure.sceneAxes stick_figure_pelvis.sceneAxes], {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'}); 
%     setappdata(gcf, 'StoreTheLink', Link);
%     Link = linkprop([stick_figure.sceneAxes stick_figure_left_leg.sceneAxes], {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'}); 
%     setappdata(gcf, 'StoreTheLink', Link);
%     Link = linkprop([stick_figure.sceneAxes stick_figure_right_leg.sceneAxes], {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'}); 
%     setappdata(gcf, 'StoreTheLink', Link);
%     Link = linkprop([stick_figure.sceneAxes stick_figure_torso.sceneAxes], {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'}); 
%     setappdata(gcf, 'StoreTheLink', Link);
%     Link = linkprop([stick_figure.sceneAxes stick_figure_left_arm.sceneAxes], {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'}); 
%     setappdata(gcf, 'StoreTheLink', Link);
%     Link = linkprop([stick_figure.sceneAxes stick_figure_right_arm.sceneAxes], {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'}); 
%     setappdata(gcf, 'StoreTheLink', Link);
    
    
%     %% test the marker positions from the limb chains
%     % set some configuration
%     theta = 0.1*ones(number_of_joints, 1);
%     kinematic_tree.jointAngles = theta;
%     pelvis_chain.jointAngles = theta(pelvis_joints);
%     left_leg_chain.jointAngles = theta(left_leg_joints);
%     right_leg_chain.jointAngles = theta(right_leg_joints);
%     torso_chain.jointAngles = theta(torso_joints);
% 
%     kinematic_tree.updateConfiguration();
%     pelvis_chain.updateConfiguration();
%     right_leg_chain.updateConfiguration();
%     left_leg_chain.updateConfiguration();
%     torso_chain.updateConfiguration();
% 
%     right_thigh_marker_positions_from_right_leg = right_leg_chain.markerPositions{last_right_hip_joint_index_in_limb_plant};
%     right_thigh_marker_positions_from_full = kinematic_tree.markerPositions{right_hip_joints(end)};
%     pelvis_to_world_poe = pelvis_chain.productsOfExponentials{pelvis_joints(end)};
%     right_thigh_marker_positions_from_right_leg_transformed = pelvis_to_world_poe * right_thigh_marker_positions_from_right_leg; % should be equal to right_thigh_marker_positions_from_full
% 
% 
% 
    
    %% optimize
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            % load data
            condition = condition_list{i_condition};
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'kinematicTrajectories')]);
            
            number_of_time_steps = size(marker_trajectories, 1);

            % optimize
            joint_angle_trajectories_calculated = joint_angle_trajectories;
            
            
            joint_angle_trajectories_optimized = zeros(size(joint_angle_trajectories_calculated));
            
            weight_matrix = ones(1, size(marker_trajectories, 2)/3); % TODO: make this a setting
            
            tic
            time_steps_to_optimize = 1 : number_of_time_steps;
%             time_steps_to_optimize = 1 : 2;
            joint_angle_trajectories_optimized(time_steps_to_optimize, :) = ...
            optimizeJointAngles ...
            ( ...
              kinematic_tree, ...
              pelvis_chain, ...
              left_leg_chain, ...
              right_leg_chain, ...
              torso_chain, ...
              left_arm_chain, ...
              right_arm_chain, ...
              marker_trajectories(time_steps_to_optimize, :), ...
              joint_angle_trajectories_calculated(time_steps_to_optimize, :), ...
              weight_matrix ...
            );
            toc
                
            % calculate CoM
            joint_center_trajectories_calculated = joint_center_trajectories;
            joint_center_trajectories_optimized = zeros(number_of_time_steps, length(joint_center_headers)*3);
            com_trajectories_calculated = com_trajectories;
            com_trajectories_optimized = zeros(number_of_time_steps, length(com_labels)*3);
            for i_time = time_steps_to_optimize
                % set kinematic tree configuration
                theta = joint_angle_trajectories_calculated(i_time, :)';
                kinematic_tree.jointAngles = theta;
                kinematic_tree.updateKinematics;
                
                % calculate joint center positions
                for i_center = 1 : length(joint_center_headers)
                    joint_center = kinematic_tree.getPointOfInterestPosition(joint_center_headers{i_center});
                    joint_center_trajectories_optimized(i_time, (i_center-1)*3 + [1 2 3]) = joint_center;
                end
                
                % get segment CoMs
                for i_segment = 1 : length(com_labels) - 1
                    segment_label = com_labels{i_segment}(1 : end-3);
                    joint_index = getSegmentJointIndex(kinematic_tree, segment_label);
                    segment_com = kinematic_tree.linkTransformations{joint_index}(1:3, 4);
                    com_trajectories_optimized(i_time, (i_segment-1)*3 + [1 2 3]) = segment_com;
                end
                com_trajectories_optimized(i_time, end-2 : end) = kinematic_tree.calculateCenterOfMassPosition;
            end
            
            % save
            variables_to_save = struct;
            variables_to_save.joint_angle_trajectories_calculated = joint_angle_trajectories_calculated;
            variables_to_save.joint_angle_trajectories_optimized = joint_angle_trajectories_optimized;
            variables_to_save.joint_center_trajectories_calculated = joint_center_trajectories_calculated;
            variables_to_save.joint_center_trajectories_optimized = joint_center_trajectories_optimized;
            variables_to_save.com_trajectories_calculated = com_trajectories_calculated;
            variables_to_save.com_trajectories_optimized = com_trajectories_optimized;
            
            save_folder = 'processed';
            save_file_name = makeFileName(date, subject_id, condition, i_trial, 'kinematicTrajectories.mat');
            saveDataToFile([save_folder filesep save_file_name], variables_to_save);
            disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' save_folder filesep save_file_name]);

            addAvailableData('joint_center_trajectories_calculated', 'time_mocap', 'sampling_rate_mocap', 'joint_center_labels', save_folder, save_file_name);
            addAvailableData('joint_center_trajectories_optimized', 'time_mocap', 'sampling_rate_mocap', 'joint_center_labels', save_folder, save_file_name);
            addAvailableData('com_trajectories_calculated', 'time_mocap', 'sampling_rate_mocap', 'com_labels', save_folder, save_file_name);
            addAvailableData('com_trajectories_optimized', 'time_mocap', 'sampling_rate_mocap', 'com_labels', save_folder, save_file_name);
            addAvailableData('joint_angle_trajectories_calculated', 'time_mocap', 'sampling_rate_mocap', 'joint_angle_trajectories_calculated', save_folder, save_file_name);
            addAvailableData('joint_angle_trajectories_optimized', 'time_mocap', 'sampling_rate_mocap', 'joint_angle_trajectories_optimized', save_folder, save_file_name);
        end
    end
end