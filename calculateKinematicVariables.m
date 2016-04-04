% reconstruct the marker trajectories from the joint angle trajectories

% flags
reconstruct     = 1;
use_parallel    = 1;

% trials_to_process = 1001;
% trials_to_process = 6;
% trials_to_process = 1001 : 1180;
trials_to_process = 4;

load subjectInfo.mat;
model_file_name = makeFileName(date, subject_id, 'model');
load(model_file_name);

if use_parallel
    poolobject = gcp;
    number_of_labs = poolobject.NumWorkers;
end

disp('calculating kinematic variables');
for i_trial = trials_to_process
    tic

    % load data
    marker_trajectories_file_name = makeFileName(date, subject_id, 'walking', i_trial, 'markerTrajectories');
    load(marker_trajectories_file_name);
    angle_trajectories_file_name = makeFileName(date, subject_id, 'walking', i_trial, 'angleTrajectories');
    load(angle_trajectories_file_name);
    
    number_of_time_steps = size(angle_trajectories, 1);
    
    marker_trajectories_reconstructed = zeros(size(marker_trajectories));
    right_heel_trajectory = zeros(number_of_time_steps, 3);
    right_toes_trajectory = zeros(number_of_time_steps, 3);
    left_heel_trajectory = zeros(number_of_time_steps, 3);
    left_toes_trajectory = zeros(number_of_time_steps, 3);
    pelvis_roll_angle = zeros(number_of_time_steps, 1);
    trunk_roll_angle = zeros(number_of_time_steps, 1);
    head_roll_angle = zeros(number_of_time_steps, 1);
    right_foot_roll_angle = zeros(number_of_time_steps, 1);
    left_foot_roll_angle = zeros(number_of_time_steps, 1);
    T_left_ankle_to_world_trajectory = zeros(number_of_time_steps, 16);
    T_right_ankle_to_world_trajectory = zeros(number_of_time_steps, 16);
    com_trajectory = zeros(number_of_time_steps, 3);

   if use_parallel
        % get or open pool of workers
        marker_trajectories_pool = zeros(size(marker_trajectories_reconstructed));
        left_heel_trajectory_pool = zeros(number_of_time_steps, 3);
        left_toes_trajectory_pool = zeros(number_of_time_steps, 3);
        right_heel_trajectory_pool = zeros(number_of_time_steps, 3);
        right_toes_trajectory_pool = zeros(number_of_time_steps, 3);
        pelvis_roll_angle_trajectory_pool = zeros(number_of_time_steps, 1);
        trunk_roll_angle_trajectory_pool = zeros(number_of_time_steps, 1);
        head_roll_angle_trajectory_pool = zeros(number_of_time_steps, 1);
        left_foot_roll_angle_trajectory_pool = zeros(number_of_time_steps, 1);
        right_foot_roll_angle_trajectory_pool = zeros(number_of_time_steps, 1);
        T_left_ankle_to_world_trajectory_pool = zeros(number_of_time_steps, 16);
        T_right_ankle_to_world_trajectory_pool = zeros(number_of_time_steps, 16);
        com_trajectory_pool = zeros(number_of_time_steps, 3);


        % create a copy of the plant for each worker
        spmd
            plant_pool = plant.copy;
            for i_time = labindex : numlabs : number_of_time_steps
                joint_angles = angle_trajectories(i_time, :)';
                if any(isnan(joint_angles))
                    marker_trajectories_pool(i_time, :) = NaN;
                    left_heel_trajectory_pool(i_time, :) = NaN;
                    left_toes_trajectory_pool(i_time, :) = NaN;
                    right_heel_trajectory_pool(i_time, :) = NaN;
                    right_toes_trajectory_pool(i_time, :) = NaN;
                    pelvis_roll_angle_trajectory_pool(i_time, :) = NaN;
                    trunk_roll_angle_trajectory_pool(i_time, :) = NaN;
                    head_roll_angle_trajectory_pool(i_time, :) = NaN;
                    left_foot_roll_angle_trajectory_pool(i_time, :) = NaN;
                    right_foot_roll_angle_trajectory_pool(i_time, :) = NaN;
                    T_left_ankle_to_world_trajectory_pool(i_time, :) = NaN;
                    T_right_ankle_to_world_trajectory_pool(i_time, :) = NaN;
                    com_trajectory_pool(i_time, :) = NaN;
                else
                    plant_pool.jointAngles = joint_angles;
                    plant_pool.updateKinematics();
                    marker_trajectories_pool(i_time, :) = plant_pool.exportMarkerPositions();
                    left_heel_trajectory_pool(i_time, :) = plant_pool.endEffectorPositions{1};
                    left_toes_trajectory_pool(i_time, :) = plant_pool.endEffectorPositions{2};
                    right_heel_trajectory_pool(i_time, :) = plant_pool.endEffectorPositions{4};
                    right_toes_trajectory_pool(i_time, :) = plant_pool.endEffectorPositions{5};
                    pelvis_roll_angle_trajectory_pool(i_time, :) = -atan2(plant_pool.jointTransformations{6}(1, 3), plant_pool.jointTransformations{6}(3, 3));
                    trunk_roll_angle_trajectory_pool(i_time, :) = -atan2(plant_pool.jointTransformations{21}(1, 3), plant_pool.jointTransformations{21}(3, 3));
                    head_roll_angle_trajectory_pool(i_time, :) = -atan2(plant_pool.jointTransformations{24}(1, 3), plant_pool.jointTransformations{24}(3, 3));
                    left_foot_roll_angle_trajectory_pool(i_time, :) = -atan2(plant_pool.jointTransformations{12}(1, 2), plant_pool.jointTransformations{12}(3, 2));
                    right_foot_roll_angle_trajectory_pool(i_time, :) = -atan2(plant_pool.jointTransformations{18}(1, 2), plant_pool.jointTransformations{18}(3, 2));
                    T_left_ankle_to_world_trajectory_pool(i_time, :) = reshape(plant_pool.endEffectorTransformations{3}, 1, 16);
                    T_right_ankle_to_world_trajectory_pool(i_time, :) = reshape(plant_pool.endEffectorTransformations{6}, 1, 16);
                    com_trajectory_pool(i_time, :) = plant_pool.calculateCenterOfMassPosition();
                end

            end
        end
        for i_lab = 1 : number_of_labs
            marker_trajectories_lab = marker_trajectories_pool{i_lab};
            left_heel_trajectory_lab = left_heel_trajectory_pool{i_lab};
            left_toes_trajectory_lab = left_toes_trajectory_pool{i_lab};
            right_heel_trajectory_lab = right_heel_trajectory_pool{i_lab};
            right_toes_trajectory_lab = right_toes_trajectory_pool{i_lab};
            pelvis_roll_angle_trajectory_lab = pelvis_roll_angle_trajectory_pool{i_lab};
            trunk_roll_angle_trajectory_lab = trunk_roll_angle_trajectory_pool{i_lab};
            head_roll_angle_trajectory_lab = head_roll_angle_trajectory_pool{i_lab};
            left_foot_roll_angle_trajectory_lab = left_foot_roll_angle_trajectory_pool{i_lab};
            right_foot_roll_angle_trajectory_lab = right_foot_roll_angle_trajectory_pool{i_lab};
            T_left_ankle_to_world_trajectory_lab = T_left_ankle_to_world_trajectory_pool{i_lab};
            T_right_ankle_to_world_trajectory_lab = T_right_ankle_to_world_trajectory_pool{i_lab};
            com_trajectory_lab = com_trajectory_pool{i_lab};

            marker_trajectories_reconstructed(i_lab : number_of_labs : number_of_time_steps, :) = marker_trajectories_lab(i_lab : number_of_labs : number_of_time_steps, :);
            left_heel_trajectory(i_lab : number_of_labs : number_of_time_steps, :) = left_heel_trajectory_lab(i_lab : number_of_labs : number_of_time_steps, :);
            left_toes_trajectory(i_lab : number_of_labs : number_of_time_steps, :) = left_toes_trajectory_lab(i_lab : number_of_labs : number_of_time_steps, :);
            right_heel_trajectory(i_lab : number_of_labs : number_of_time_steps, :) = right_heel_trajectory_lab(i_lab : number_of_labs : number_of_time_steps, :);
            right_toes_trajectory(i_lab : number_of_labs : number_of_time_steps, :) = right_toes_trajectory_lab(i_lab : number_of_labs : number_of_time_steps, :);
            pelvis_roll_angle(i_lab : number_of_labs : number_of_time_steps, :) = pelvis_roll_angle_trajectory_lab(i_lab : number_of_labs : number_of_time_steps, :);
            trunk_roll_angle(i_lab : number_of_labs : number_of_time_steps, :) = trunk_roll_angle_trajectory_lab(i_lab : number_of_labs : number_of_time_steps, :);
            head_roll_angle(i_lab : number_of_labs : number_of_time_steps, :) = head_roll_angle_trajectory_lab(i_lab : number_of_labs : number_of_time_steps, :);
            left_foot_roll_angle(i_lab : number_of_labs : number_of_time_steps, :) = left_foot_roll_angle_trajectory_lab(i_lab : number_of_labs : number_of_time_steps, :);
            right_foot_roll_angle(i_lab : number_of_labs : number_of_time_steps, :) = right_foot_roll_angle_trajectory_lab(i_lab : number_of_labs : number_of_time_steps, :);
            T_left_ankle_to_world_trajectory(i_lab : number_of_labs : number_of_time_steps, :) = T_left_ankle_to_world_trajectory_lab(i_lab : number_of_labs : number_of_time_steps, :);
            T_right_ankle_to_world_trajectory(i_lab : number_of_labs : number_of_time_steps, :) = T_right_ankle_to_world_trajectory_lab(i_lab : number_of_labs : number_of_time_steps, :);

            com_trajectory(i_lab : number_of_labs : number_of_time_steps, :) = com_trajectory_lab(i_lab : number_of_labs : number_of_time_steps, :);
        end

    else
%         for i_time = 1 : number_of_time_steps
        for i_time = 1540 : 1550
            joint_angles = angle_trajectories(i_time, :)';
            if any(isnan(joint_angles))
                marker_trajectories_reconstructed(i_time, :) = NaN;
                left_heel_trajectory(i_time, :) = NaN;
                left_toes_trajectory(i_time, :) = NaN;
                right_heel_trajectory(i_time, :) = NaN;
                right_toes_trajectory(i_time, :) = NaN;
                pelvis_roll_angle(i_time, :) = NaN;
                trunk_roll_angle(i_time, :) = NaN;
                head_roll_angle(i_time, :) = NaN;
                left_foot_roll_angle(i_time, :) = NaN;
                right_foot_roll_angle(i_time, :) = NaN;
                T_left_ankle_to_world_trajectory(i_time, :) = NaN;
                T_right_ankle_to_world_trajectory(i_time, :) = NaN;
                com_trajectory(i_time, :) = NaN;
            else
                plant.jointAngles = joint_angles;
                plant.updateKinematics();
                marker_trajectories_reconstructed(i_time, :) = plant.exportMarkerPositions();
                left_heel_trajectory(i_time, :) = plant.endEffectorPositions{1};
                left_toes_trajectory(i_time, :) = plant.endEffectorPositions{2};
                right_heel_trajectory(i_time, :) = plant.endEffectorPositions{4};
                right_toes_trajectory(i_time, :) = plant.endEffectorPositions{5};
                pelvis_roll_angle(i_time, :) = atan2(plant.jointTransformations{6}(2, 3), plant.jointTransformations{6}(3, 3));
                trunk_roll_angle(i_time, :) = atan2(plant.jointTransformations{21}(2, 3), plant.jointTransformations{21}(3, 3));
                head_roll_angle(i_time, :) = atan2(plant.jointTransformations{24}(2, 3), plant.jointTransformations{24}(3, 3));
                right_foot_roll_angle(i_time, :) = atan2(plant.jointTransformations{18}(2, 2), plant.jointTransformations{18}(3, 2));
                left_foot_roll_angle(i_time, :) = atan2(plant.jointTransformations{12}(2, 2), plant.jointTransformations{12}(3, 2));
                T_left_ankle_to_world_trajectory(i_time, :) = reshape(plant.endEffectorTransformations{3}, 1, 16);
                T_right_ankle_to_world_trajectory(i_time, :) = reshape(plant.endEffectorTransformations{6}, 1, 16);
                com_trajectory(i_time, :) = plant.calculateCenterOfMassPosition();
            end
            step = 100;
            if ((i_time / step) == floor(i_time / step)) && ~ ((i_time / (5*step)) == floor(i_time / (5*step)))
                fprintf('.');
            elseif ((i_time / (5*step)) == floor(i_time / (5*step))) && ~ ((i_time / (10*step)) == floor(i_time / (10*step)))
                fprintf(',');
            elseif ((i_time / (10*step)) == floor(i_time / (10*step)))
                fprintf('|');
            end
        end
    end        

    % calculate the reconstruction errors
    number_of_markers = plant.getNumberOfMarkers();
    number_of_time_steps = size(marker_trajectories_reconstructed, 1);
    marker_reconstruction_error = zeros(number_of_time_steps, number_of_markers);
    for i_marker = 1 : number_of_markers
        for i_time = 1 : number_of_time_steps
            marker_position_measured = marker_trajectories(i_time, 3*(i_marker-1)+1 : 3*(i_marker-1)+3);
            marker_position_reconstr = marker_trajectories_reconstructed(i_time, 3*(i_marker-1)+1 : 3*(i_marker-1)+3);
            error_vector = marker_position_measured - marker_position_reconstr;
            marker_reconstruction_error(i_time, i_marker) = norm(error_vector);
        end
    end

    %% save data
    kinematic_trajectories_file_name = makeFileName(date, subject_id, 'walking', i_trial, 'kinematicTrajectories');
    save ...
      ( ...
        kinematic_trajectories_file_name, ...
        'marker_trajectories_reconstructed', ...
        'marker_reconstruction_error', ...
        'right_heel_trajectory', ...
        'right_toes_trajectory', ...
        'left_heel_trajectory', ...
        'left_toes_trajectory', ...
        'pelvis_roll_angle', ...
        'trunk_roll_angle', ...
        'head_roll_angle', ...
        'right_foot_roll_angle', ...
        'left_foot_roll_angle', ...
        'T_left_ankle_to_world_trajectory', ...
        'T_right_ankle_to_world_trajectory', ...
        'com_trajectory' ...
      );

    
    disp(['Trial ' num2str(i_trial) ' completed']);
    toc
end
fprintf(' done\n');



