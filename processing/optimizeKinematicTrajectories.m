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
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'use_parallel', false)
    parse(parser, varargin{:})
    use_parallel = parser.Results.use_parallel;
    
    
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
    
    if use_parallel
        % get or open pool of workers
        poolobject = gcp;
        number_of_labs = poolobject.NumWorkers;
    end    
    
    %% determine weights
    weight_matrix = ones(1, length(kinematic_tree.markerLabels));
    marker_weight_table = study_settings.get('marker_weights');
    for i_weight = 1 : size(marker_weight_table, 1)
        this_weight_marker_label = marker_weight_table{i_weight, 1};
        this_weight = str2double(marker_weight_table{i_weight, 2});
        
        this_weight_marker_index = find(strcmp(kinematic_tree.markerLabels, this_weight_marker_label), 1, 'first');
        weight_matrix(this_weight_marker_index) = this_weight;
    end
    
    %% create limb kinematic_trees

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
    
    left_toes_eef_index = kinematic_tree.getEndEffectorIndex('LTOESEEF');
    right_toes_eef_index = kinematic_tree.getEndEffectorIndex('RTOESEEF');
    head_eef_index = kinematic_tree.getEndEffectorIndex('head');
    left_hand_eef_index = kinematic_tree.getEndEffectorIndex('LHANDEEF');
    right_hand_eef_index = kinematic_tree.getEndEffectorIndex('RHANDEEF');

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
    last_left_hip_joint_index_in_limb_kinematic_tree = left_hip_joints(end) - (left_leg_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{left_hip_joints(end)}, 2)
        left_leg_chain.addMarker(last_left_hip_joint_index_in_limb_kinematic_tree, kinematic_tree.markerReferencePositions{left_hip_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(left_hip_joints(end), i_marker));
    end
    last_left_knee_joint_index_in_limb_kinematic_tree = left_knee_joints(end) - (left_leg_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{left_knee_joints(end)}, 2)
        left_leg_chain.addMarker(last_left_knee_joint_index_in_limb_kinematic_tree, kinematic_tree.markerReferencePositions{left_knee_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(left_knee_joints(end), i_marker));
    end
    last_left_ankle_joint_index_in_limb_kinematic_tree = left_ankle_joints(end) - (left_leg_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{left_ankle_joints(end)}, 2)
        left_leg_chain.addMarker(last_left_ankle_joint_index_in_limb_kinematic_tree, kinematic_tree.markerReferencePositions{left_ankle_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(left_ankle_joints(end), i_marker));
    end
    
    % add markers to right leg tree
    last_right_hip_joint_index_in_limb_kinematic_tree = right_hip_joints(end) - (right_leg_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{right_hip_joints(end)}, 2)
        right_leg_chain.addMarker(last_right_hip_joint_index_in_limb_kinematic_tree, kinematic_tree.markerReferencePositions{right_hip_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(right_hip_joints(end), i_marker));
    end
    last_right_knee_joint_index_in_limb_kinematic_tree = right_knee_joints(end) - (right_leg_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{right_knee_joints(end)}, 2)
        right_leg_chain.addMarker(last_right_knee_joint_index_in_limb_kinematic_tree, kinematic_tree.markerReferencePositions{right_knee_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(right_knee_joints(end), i_marker));
    end
    last_right_ankle_joint_index_in_limb_kinematic_tree = right_ankle_joints(end) - (right_leg_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{right_ankle_joints(end)}, 2)
        right_leg_chain.addMarker(last_right_ankle_joint_index_in_limb_kinematic_tree, kinematic_tree.markerReferencePositions{right_ankle_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(right_ankle_joints(end), i_marker));
    end
    
    % add markers to torso tree
    last_lumbar_joint_index_in_limb_kinematic_tree = lumbar_joints(end) - (torso_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{lumbar_joints(end)}, 2)
        torso_chain.addMarker(last_lumbar_joint_index_in_limb_kinematic_tree, kinematic_tree.markerReferencePositions{lumbar_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(lumbar_joints(end), i_marker));
    end
    last_cervix_joint_index_in_limb_kinematic_tree = cervix_joints(end) - (torso_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{cervix_joints(end)}, 2)
        torso_chain.addMarker(last_cervix_joint_index_in_limb_kinematic_tree, kinematic_tree.markerReferencePositions{cervix_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(cervix_joints(end), i_marker));
    end
    
    % add markers to left arm tree
    last_left_shoulder_joint_index_in_limb_kinematic_tree = left_shoulder_joints(end) - (left_arm_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{left_shoulder_joints(end)}, 2)
        left_arm_chain.addMarker(last_left_shoulder_joint_index_in_limb_kinematic_tree, kinematic_tree.markerReferencePositions{left_shoulder_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(left_shoulder_joints(end), i_marker));
    end
    last_left_elbow_joint_index_in_limb_kinematic_tree = left_elbow_joints(end) - (left_arm_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{left_elbow_joints(end)}, 2)
        left_arm_chain.addMarker(last_left_elbow_joint_index_in_limb_kinematic_tree, kinematic_tree.markerReferencePositions{left_elbow_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(left_elbow_joints(end), i_marker));
    end
    last_left_wrist_joint_index_in_limb_kinematic_tree = left_wrist_joints(end) - (left_arm_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{left_wrist_joints(end)}, 2)
        left_arm_chain.addMarker(last_left_wrist_joint_index_in_limb_kinematic_tree, kinematic_tree.markerReferencePositions{left_wrist_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(left_wrist_joints(end), i_marker));
    end
    
    % add markers to right arm tree
    last_right_shoulder_joint_index_in_limb_kinematic_tree = right_shoulder_joints(end) - (right_arm_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{right_shoulder_joints(end)}, 2)
        right_arm_chain.addMarker(last_right_shoulder_joint_index_in_limb_kinematic_tree, kinematic_tree.markerReferencePositions{right_shoulder_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(right_shoulder_joints(end), i_marker));
    end
    last_right_elbow_joint_index_in_limb_kinematic_tree = right_elbow_joints(end) - (right_arm_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{right_elbow_joints(end)}, 2)
        right_arm_chain.addMarker(last_right_elbow_joint_index_in_limb_kinematic_tree, kinematic_tree.markerReferencePositions{right_elbow_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(right_elbow_joints(end), i_marker));
    end
    last_right_wrist_joint_index_in_limb_kinematic_tree = right_wrist_joints(end) - (right_arm_joints(1) - 1);
    for i_marker = 1 : size(kinematic_tree.markerReferencePositions{right_wrist_joints(end)}, 2)
        right_arm_chain.addMarker(last_right_wrist_joint_index_in_limb_kinematic_tree, kinematic_tree.markerReferencePositions{right_wrist_joints(end)}(1:3, i_marker), kinematic_tree.getMarkerVisualizationColor(right_wrist_joints(end), i_marker));
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
%     right_thigh_marker_positions_from_right_leg = right_leg_chain.markerPositions{last_right_hip_joint_index_in_limb_kinematic_tree};
%     right_thigh_marker_positions_from_full = kinematic_tree.markerPositions{right_hip_joints(end)};
%     pelvis_to_world_poe = pelvis_chain.productsOfExponentials{pelvis_joints(end)};
%     right_thigh_marker_positions_from_right_leg_transformed = pelvis_to_world_poe * right_thigh_marker_positions_from_right_leg; % should be equal to right_thigh_marker_positions_from_full
% 
% 
% 
    
    %% process
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            %% load data
            condition = condition_list{i_condition};
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'kinematicTrajectories')]);
            
%             % modify marker trajectories: set LTHIA marker to NaN to check if treatment of weightless markers is working properly
%             LTHIA_number = find(strcmp(marker_labels, 'LTHIA'));
%             LTHIA_indices = (LTHIA_number-1)*3 + (1:3);
%             marker_trajectories(:, LTHIA_indices) = NaN;
%             joint_angle_trajectories(:, 7) = NaN;
            
            number_of_time_steps = size(marker_trajectories, 1);

            % determine time steps to optimize
%             time_steps_to_optimize = 1 : number_of_time_steps;
%             time_steps_to_optimize = 1 : 20;
%             time_steps_to_optimize = 29999 : 30000;
%             time_steps_to_optimize = 30000;
            
            time_steps_to_optimize = determineTimeStepsToProcess(date, subject_id, condition, i_trial, study_settings.get('data_stretch_padding'));
%             time_steps_to_optimize = determineTimeStepsToProcess(date, subject_id, condition, i_trial, 0);

% cut down time steps to optimize for debugging
% time_steps_to_optimize = time_steps_to_optimize(1006);
% time_steps_to_optimize = time_steps_to_optimize(1005 : 1006);
% time_steps_to_optimize = 3026 : 3030;



            number_of_time_steps_to_optimize = length(time_steps_to_optimize);
            
            %% optimize
            joint_angle_trajectories_calculated = joint_angle_trajectories;
            joint_angle_trajectories_optimized = zeros(size(joint_angle_trajectories_calculated)) * NaN;
            
            TODO: check what is used as initial condition, what happens if this contains NaNs? Maybe that's the
            problem with the padding?
            
            tic
%             disp([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Condition ' condition ', Trial ' num2str(i_trial)])
%             fprintf([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Optimizing joint angles... \n'])
            disp([' - Condition ' condition ', Trial ' num2str(i_trial)])
            fprintf([' - Optimizing joint angles... \n'])
            if use_parallel
                joint_angle_trajectories_optimized_pool = zeros(size(joint_angle_trajectories_optimized));
                


                spmd
                    % create a copy of the kinematic_tree for each worker
                    kinematic_tree_pool = kinematic_tree.copy;
                    pelvis_chain_pool = pelvis_chain.copy;
                    left_leg_chain_pool = left_leg_chain.copy;
                    right_leg_chain_pool = right_leg_chain.copy;
                    torso_chain_pool = torso_chain.copy;
                    left_arm_chain_pool = left_arm_chain.copy;
                    right_arm_chain_pool = right_arm_chain.copy;

                    time_steps_to_optimize_lab = time_steps_to_optimize(labindex : numlabs : number_of_time_steps_to_optimize);

                    joint_angle_trajectories_optimized_pool(time_steps_to_optimize_lab, :) = ...
                    optimizeJointAngles ...
                    ( ...
                      kinematic_tree_pool, ...
                      pelvis_chain_pool, ...
                      left_leg_chain_pool, ...
                      right_leg_chain_pool, ...
                      torso_chain_pool, ...
                      left_arm_chain_pool, ...
                      right_arm_chain_pool, ...
                      marker_trajectories(time_steps_to_optimize_lab, :), ...
                      joint_angle_trajectories_calculated(time_steps_to_optimize, :), ...
                      weight_matrix ...
                    );        
                end

                % reassemble
                for i_lab = 1 : number_of_labs
                    joint_angle_trajectories_optimized_lab = joint_angle_trajectories_optimized_pool{i_lab};
                    joint_angle_trajectories_optimized(time_steps_to_optimize(i_lab : number_of_labs : number_of_time_steps_to_optimize), :) = joint_angle_trajectories_optimized_lab(time_steps_to_optimize(i_lab : number_of_labs : number_of_time_steps_to_optimize), :);
                end                  
            end
            if ~use_parallel
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
            end
            joint_angle_trajectories_optimized = normalizeAngle(joint_angle_trajectories_optimized);
%             fprintf([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - finished\n'])
            fprintf([' - finished\n'])
            toc
                
            %% get joint centers and CoM from the kinematic tree
%             fprintf([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Calculating joint centers and CoM... \n'])
            fprintf([' - Calculating joint centers and CoM... \n'])
            joint_center_trajectories_calculated = joint_center_trajectories;
            joint_center_trajectories_optimized = zeros(number_of_time_steps, length(joint_center_headers)*3);
            com_trajectories_calculated = com_trajectories;
            com_trajectories_optimized = zeros(number_of_time_steps, length(com_labels)*3);
            
            tic
            if use_parallel
                % make variables accessible to workers by declaring them
                joint_center_headers_pool = joint_center_headers;
                joint_center_trajectories_optimized_pool = joint_center_trajectories_optimized;
                com_trajectories_optimized_pool = com_trajectories_optimized;
                com_labels_pool = com_labels;
                spmd
                    kinematic_tree_pool = kinematic_tree.copy;
                    for i_time_index = labindex : numlabs : number_of_time_steps_to_optimize
                        i_time = time_steps_to_optimize(i_time_index);
                        if any(any(isnan(joint_center_trajectories_optimized_pool(i_time, :))))
                            joint_center_trajectories_optimized_pool(i_time, :) = NaN;
                            com_trajectories_optimized_pool(i_time, :) = NaN;
                        else
                            % set kinematic tree configuration
                            theta = joint_angle_trajectories_optimized(i_time, :)';
                            kinematic_tree_pool.jointAngles = theta;
                            kinematic_tree_pool.updateKinematics;

                            % calculate joint center positions
                            for i_center = 1 : length(joint_center_headers_pool)
                                this_joint_label = joint_center_headers_pool{i_center};

                                % try as joint
                                joint_indices = kinematic_tree_pool.getJointGroup(this_joint_label);
                                if ~isempty(joint_indices)
                                    joint_position = kinematic_tree_pool.jointTransformations{joint_indices(end)}(1:3, 4);
                                    joint_center_trajectories_optimized_pool(i_time, (i_center-1)*3 + [1 2 3]) = joint_position;
                                end

                                % try as end-effector
                                point_indices = kinematic_tree_pool.getEndEffectorIndex(this_joint_label);
                                if ~isempty(point_indices)
                                    point_position = kinematic_tree_pool.endEffectorPositions{point_indices(end)};
                                    joint_center_trajectories_optimized_pool(i_time, (i_center-1)*3 + [1 2 3]) = point_position;
                                end
                            end

                            % get segment CoMs
                            for i_segment = 1 : length(com_labels_pool) - 1
                                segment_label = com_labels_pool{i_segment}(1 : end-3);
                                joint_index = getSegmentJointIndex(kinematic_tree_pool, segment_label);
                                segment_com = kinematic_tree_pool.linkTransformations{joint_index}(1:3, 4);
                                com_trajectories_optimized_pool(i_time, (i_segment-1)*3 + [1 2 3]) = segment_com;
                            end
                            com_trajectories_optimized_pool(i_time, end-2 : end) = kinematic_tree_pool.calculateCenterOfMassPosition;                        
                        end
                    end
                end
                % reassemble
                for i_lab = 1 : number_of_labs
                    joint_center_trajectories_optimized_lab = joint_center_trajectories_optimized_pool{i_lab};
                    joint_center_trajectories_optimized(time_steps_to_optimize(i_lab : number_of_labs : number_of_time_steps_to_optimize), :) = joint_center_trajectories_optimized_lab(time_steps_to_optimize(i_lab : number_of_labs : number_of_time_steps_to_optimize), :);
                    com_trajectories_optimized_lab = com_trajectories_optimized_pool{i_lab};
                    com_trajectories_optimized(time_steps_to_optimize(i_lab : number_of_labs : number_of_time_steps_to_optimize), :) = com_trajectories_optimized_lab(time_steps_to_optimize(i_lab : number_of_labs : number_of_time_steps_to_optimize), :);
                end
            end
            if ~use_parallel
                for i_time_step = 1 : length(time_steps_to_process)
                    i_time = time_steps_to_process(i_time_step);
                    if any(any(isnan(joint_angle_trajectories_optimized(i_time, :))))
                        joint_center_trajectories_optimized(i_time, :) = NaN;
                        com_trajectories_optimized(i_time, :) = NaN;
                    else
                    
                    
                        % set kinematic tree configuration
                        theta = joint_angle_trajectories_optimized(i_time, :)';
                        kinematic_tree.jointAngles = theta;
                        kinematic_tree.updateKinematics;

                        % calculate joint center positions
                        for i_center = 1 : length(joint_center_headers)
                            this_joint_label = joint_center_headers{i_center};

                            % try as joint
                            joint_indices = kinematic_tree.getJointGroup(this_joint_label);
                            if ~isempty(joint_indices)
                                joint_position = kinematic_tree.jointTransformations{joint_indices(end)}(1:3, 4);
                                joint_center_trajectories_optimized(i_time, (i_center-1)*3 + [1 2 3]) = joint_position;
                            end

                            % try as end-effector
                            point_indices = kinematic_tree.getEndEffectorIndex(this_joint_label);
                            if ~isempty(point_indices)
                                point_position = kinematic_tree.endEffectorPositions{point_indices(end)};
                                joint_center_trajectories_optimized(i_time, (i_center-1)*3 + [1 2 3]) = point_position;
                            end


        %                     joint_center = kinematic_tree.getPointOfInterestPosition(joint_center_headers{i_center});
        %                     joint_center_trajectories_optimized(i_time, (i_center-1)*3 + [1 2 3]) = joint_center;
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
                    
                    % give progress feedback
                    display_step = 10;
                    if (i_time_step / display_step) == floor(i_time_step / display_step)
                        disp([num2str(i_time_step) '(' num2str(length(time_steps_to_process)) ')']);
                    end                        
                    
                end                
            end
            toc
            
            
            
%             figure; axes; hold on
%             plot(com_trajectories_optimized);
            
            
            
            
            
            
            
            

%             fprintf([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - finished\n'])
            fprintf([' - finished\n'])
%             fprintf('finished\n')

            %% save
            variables_to_save = struct;
            variables_to_save.joint_labels = kinematic_tree.jointLabels;
            variables_to_save.joint_angle_trajectories = joint_angle_trajectories_optimized;
            variables_to_save.joint_angle_trajectories_calculated = joint_angle_trajectories_calculated;
            variables_to_save.joint_angle_trajectories_optimized = joint_angle_trajectories_optimized;
            
            variables_to_save.joint_center_trajectories = joint_center_trajectories_optimized;
            variables_to_save.joint_center_trajectories_calculated = joint_center_trajectories_calculated;
            variables_to_save.joint_center_trajectories_optimized = joint_center_trajectories_optimized;
            
            variables_to_save.com_trajectories = com_trajectories_optimized;
            variables_to_save.com_trajectories_calculated = com_trajectories_calculated;
            variables_to_save.com_trajectories_optimized = com_trajectories_optimized;
            
            save_folder = 'processed';
            save_file_name = makeFileName(date, subject_id, condition, i_trial, 'kinematicTrajectories.mat');
            saveDataToFile([save_folder filesep save_file_name], variables_to_save);
            disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' save_folder filesep save_file_name]);

            addAvailableData('joint_center_trajectories_calculated', 'time_mocap', 'sampling_rate_mocap', 'joint_center_labels', save_folder, save_file_name);
            addAvailableData('com_trajectories_calculated', 'time_mocap', 'sampling_rate_mocap', 'com_labels', save_folder, save_file_name);
            addAvailableData('joint_angle_trajectories_calculated', 'time_mocap', 'sampling_rate_mocap', 'joint_labels', save_folder, save_file_name);
            addAvailableData('joint_center_trajectories_optimized', 'time_mocap', 'sampling_rate_mocap', 'joint_center_labels', save_folder, save_file_name);
            addAvailableData('com_trajectories_optimized', 'time_mocap', 'sampling_rate_mocap', 'com_labels', save_folder, save_file_name);
            addAvailableData('joint_angle_trajectories_optimized', 'time_mocap', 'sampling_rate_mocap', 'joint_labels', save_folder, save_file_name);
            
        end
    end
end











