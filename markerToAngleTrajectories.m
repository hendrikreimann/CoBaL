
use_parallel            = 0;

trials_to_process = 6;

data_directory = '/Users/reimajbi/Neuro/walking/data/20160322_EEJ'; % should this be a function argument?

load subjectInfo.mat;
model_file_name = makeFileName(date, subject_id, 'model');
load(model_file_name);

weight_matrix = ones(1, plant.getNumberOfMarkers);
number_of_joints = plant.numberOfJoints;




%% create limb plants

% extract reference positions
plant.jointAngles = zeros(plant.numberOfJoints, 1);
plant.updateConfiguration();
joint_positions = cell(1, number_of_joints);
joint_axes = cell(1, number_of_joints);
link_com_positions = cell(1, number_of_joints);
link_orientations = cell(1, number_of_joints);
generalized_link_inertia_matrices = plant.generalizedInertiaMatrices;
for i_joint = 1 : number_of_joints
    joint_positions{i_joint} = plant.jointTransformations{i_joint}(1:3, 4);
    if norm(plant.referenceJointTwists{i_joint}(4:6)) ~= 0
        joint_axes{i_joint} = plant.referenceJointTwists{i_joint}(4:6);
    else
        joint_axes{i_joint} = plant.referenceJointTwists{i_joint}(1:3);
    end
    link_com_positions{i_joint} = plant.linkTransformations{i_joint}(1:3, 4);
    link_orientations{i_joint} = plant.linkTransformations{i_joint}(1:3, 1:3);
end

% create chains
pelvis_chain = GeneralKinematicTree ...
( ...
  joint_positions(1:6), ...
  joint_axes(1:6), ...
  [2 2 2 1 1 1], ...                                                        % types
  ones(3, 6), ...                                                           % branch matrix
  {joint_positions{7}, joint_positions{13}, joint_positions{19}}, ...       % end-effectors
  link_com_positions(1:6), ...
  link_orientations(1:6), ...
  generalized_link_inertia_matrices(1:6) ...
);
left_leg_chain = GeneralKinematicTree ...
( ...
  joint_positions(7:12), ...
  joint_axes(7:12), ...
  [1 1 1 1 1 1], ...                                                        % types
  ones(1, 6), ...                                                           % branch matrix
  {plant.endEffectorPositions{1}}, ...                                      % end-effector
  link_com_positions(7:12), ...
  link_orientations(7:12), ...
  generalized_link_inertia_matrices(7:12) ...
);
right_leg_chain = GeneralKinematicTree ...
( ...
  joint_positions(13:18), ...
  joint_axes(13:18), ...
  [1 1 1 1 1 1], ...                                                        % types
  ones(1, 6), ...                                                           % branch matrix
  {plant.endEffectorPositions{4}}, ...                                      % end-effector
  link_com_positions(13:18), ...
  link_orientations(13:18), ...
  generalized_link_inertia_matrices(13:18) ...
);
trunk_chain = GeneralKinematicTree ...
( ...
  joint_positions(19:24), ...
  joint_axes(19:24), ...
  [1 1 1 1 1 1], ...                                                        % types
  ones(1, 6), ...                                                           % branch matrix
  {plant.endEffectorPositions{7}}, ...                                      % end-effector
  link_com_positions(19:24), ...
  link_orientations(19:24), ...
  generalized_link_inertia_matrices(19:24) ...
);
left_arm_chain = GeneralKinematicTree ...
( ...
  joint_positions(25:31), ...
  joint_axes(25:31), ...
  [1 1 1 1 1 1 1], ...                                                        % types
  ones(1, 7), ...                                                           % branch matrix
  {plant.endEffectorPositions{7}}, ...                                      % end-effector
  link_com_positions(25:31), ...
  link_orientations(25:31), ...
  generalized_link_inertia_matrices(25:31) ...
);
right_arm_chain = GeneralKinematicTree ...
( ...
  joint_positions(32:38), ...
  joint_axes(32:38), ...
  [1 1 1 1 1 1 1], ...                                                        % types
  ones(1, 7), ...                                                           % branch matrix
  {plant.endEffectorPositions{8}}, ...                                      % end-effector
  link_com_positions(32:38), ...
  link_orientations(32:38), ...
  generalized_link_inertia_matrices(32:38) ...
);

% add markers
for i_marker = 1 : size(plant.markerReferencePositions{6}, 2)
    pelvis_chain.addMarker(6, plant.markerReferencePositions{6}(1:3, i_marker), plant.getMarkerVisualizationColor(6, i_marker));
end
left_leg_index_offset = 6;
for i_joint = [3 4 6]
    for i_marker = 1 : size(plant.markerReferencePositions{i_joint + left_leg_index_offset}, 2)
        left_leg_chain.addMarker(i_joint, plant.markerReferencePositions{i_joint + left_leg_index_offset}(1:3, i_marker), plant.getMarkerVisualizationColor(i_joint + left_leg_index_offset, i_marker));
    end
end
right_leg_index_offset = 12;
for i_joint = [3 4 6]
    for i_marker = 1 : size(plant.markerReferencePositions{i_joint + right_leg_index_offset}, 2)
        right_leg_chain.addMarker(i_joint, plant.markerReferencePositions{i_joint + right_leg_index_offset}(1:3, i_marker), plant.getMarkerVisualizationColor(i_joint + right_leg_index_offset, i_marker));
    end
end
trunk_index_offset = 18;
for i_joint = [3 6]
    for i_marker = 1 : size(plant.markerReferencePositions{i_joint + trunk_index_offset}, 2)
        trunk_chain.addMarker(i_joint, plant.markerReferencePositions{i_joint + trunk_index_offset}(1:3, i_marker), plant.getMarkerVisualizationColor(i_joint + trunk_index_offset, i_marker));
    end
end
left_arm_index_offset = 24;
for i_joint = [3 5 7]
    for i_marker = 1 : size(plant.markerReferencePositions{i_joint + left_arm_index_offset}, 2)
        left_arm_chain.addMarker(i_joint, plant.markerReferencePositions{i_joint + left_arm_index_offset}(1:3, i_marker), plant.getMarkerVisualizationColor(i_joint + left_arm_index_offset, i_marker));
    end
end
right_arm_index_offset = 31;
for i_joint = [3 5 7]
    for i_marker = 1 : size(plant.markerReferencePositions{i_joint + right_arm_index_offset}, 2)
        right_arm_chain.addMarker(i_joint, plant.markerReferencePositions{i_joint + right_arm_index_offset}(1:3, i_marker), plant.getMarkerVisualizationColor(i_joint + right_arm_index_offset, i_marker));
    end
end

%% test the marker positions from the limb chains
% % set some configuration
% theta = 0.1*ones(24, 1);
% plant.jointAngles = theta;
% pelvis_chain.jointAngles = theta(1:6);
% left_leg_chain.jointAngles = theta(7:12);
% right_leg_chain.jointAngles = theta(13:18);
% trunk_chain.jointAngles = theta(19:24);
% 
% plant.updateConfiguration();
% pelvis_chain.updateConfiguration();
% right_leg_chain.updateConfiguration();
% left_leg_chain.updateConfiguration();
% trunk_chain.updateConfiguration();
% 
% right_thigh_marker_positions_from_right_leg = right_leg_chain.markerPositions{3};
% right_thigh_marker_positions_from_full = plant.markerPositions{15};
% pelvis_to_world_poe = pelvis_chain.productsOfExponentials{6};
% right_thigh_marker_positions_from_right_leg_transformed = pelvis_to_world_poe * right_thigh_marker_positions_from_right_leg; % should be equal to right_thigh_marker_positions_from_full









%% optimize
for i_trial = trials_to_process
    % load data
    marker_trajectories_file_name = makeFileName(date, subject_id, 'walking', i_trial, 'markerTrajectories');
    load([data_directory filesep marker_trajectories_file_name]);
    
    % schedule for optimization
    number_of_time_steps = size(marker_trajectories, 1);
    time_steps_to_optimize = 1 : number_of_time_steps;
    time_steps_to_optimize = 1540 : 1550;
    
    % run optimization
    number_of_time_steps_to_optimize = length(time_steps_to_optimize);
    joint_angle_trajectory_optimized = zeros(number_of_time_steps, number_of_joints);
    tic
    if use_parallel
        joint_angle_trajectory_optimized_pool = zeros(size(joint_angle_trajectory_optimized));
        
        % get or open pool of workers
        poolobject = gcp;
        number_of_labs = poolobject.NumWorkers;

        % create a copy of the plant for each worker
        spmd
            plant_pool = plant.copy;
            pelvis_chain_pool = pelvis_chain.copy;
            right_leg_chain_pool = right_leg_chain.copy;
            left_leg_chain_pool = left_leg_chain.copy;
            trunk_chain_pool = trunk_chain.copy;
            right_arm_chain_pool = right_arm_chain.copy;
            left_arm_chain_pool = left_arm_chain.copy;
            
            time_steps_to_optimize_lab = time_steps_to_optimize(labindex : numlabs : number_of_time_steps_to_optimize);
            
            joint_angle_trajectory_optimized_pool(time_steps_to_optimize_lab, :) = ...
            optimizeJointAngles ...
            ( ...
              plant_pool, ...
              pelvis_chain_pool, ...
              left_leg_chain_pool, ...
              right_leg_chain_pool, ...
              trunk_chain_pool, ...
              left_arm_chain_pool, ...
              right_arm_chain_pool, ...
              marker_trajectories(time_steps_to_optimize_lab, :) ...
            );        
        end
        
        % reassemble
        for i_lab = 1 : number_of_labs
            joint_angle_trajectory_optimized_lab = joint_angle_trajectory_optimized_pool{i_lab};
            joint_angle_trajectory_optimized(time_steps_to_optimize(i_lab : number_of_labs : number_of_time_steps_to_optimize), :) = joint_angle_trajectory_optimized_lab(time_steps_to_optimize(i_lab : number_of_labs : number_of_time_steps_to_optimize), :);
        end        

    else
        joint_angle_trajectory_optimized(time_steps_to_optimize, :) = ...
        optimizeJointAngles ...
        ( ...
          plant, ...
          pelvis_chain, ...
          left_leg_chain, ...
          right_leg_chain, ...
          trunk_chain, ...
          left_arm_chain, ...
          right_arm_chain, ...
          marker_trajectories(time_steps_to_optimize, :) ...
        );
    end
    toc
    
    % normalize
    angle_trajectories = normalizeAngle(joint_angle_trajectory_optimized);
    
    % save
    angle_trajectory_file_name = makeFileName(date, subject_id, 'walking', i_trial, 'angleTrajectories');
    save([data_directory filesep angle_trajectory_file_name], 'angle_trajectories');
    
    disp(['Trial ' num2str(i_trial) ' completed']);
end








