% reconstruct the marker trajectories from the joint angle trajectories

% flags
reconstruct             = 1;
use_parallel            = 1;
use_filtered_data       = 1;
visualize_derivatives   = 0;
process_all_data        = 1;

% trials_to_process = 1001;
% trials_to_process = 6;
% trials_to_process = 1001 : 1180;
trials_to_process = 3;

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

    %% load data
    load(makeFileName(date, subject_id, 'walking', i_trial, 'markerTrajectories'));
    load(makeFileName(date, subject_id, 'walking', i_trial, 'angleTrajectories'));
    load(makeFileName(date, subject_id, 'walking', i_trial, 'labviewTrajectories'));
    load(makeFileName(date, subject_id, 'walking', i_trial, 'relevantDataStretches'));
    
    filter_order = 2;
    cutoff_frequency = 20; % cutoff frequency, in Hz
    [b_filter, a_filter] = butter(filter_order, cutoff_frequency/(sampling_rate_mocap/2));
    
    number_of_time_steps = size(joint_angle_trajectories, 1);
    number_of_joints = plant.numberOfJoints;
    
    %% shift to belt coordinates
    
    % calculate offset between world coordinates and belt coordinates
    belt_speed_trajectory = mean([belt_speed_left_trajectory belt_speed_right_trajectory], 2);
    delta_t = diff(time_labview);
    belt_position_trajectory_forceplate = zeros(size(belt_speed_trajectory));
    for i_time = 2 : length(belt_speed_trajectory)
        belt_position_trajectory_forceplate(i_time) = belt_position_trajectory_forceplate(i_time-1) + delta_t(i_time-1) * belt_speed_trajectory(i_time-1);
    end
    belt_position_trajectory_mocap = spline(time_labview, belt_position_trajectory_forceplate, time_mocap)';
    
    joint_angle_trajectories_belt = joint_angle_trajectories;
    joint_angle_trajectories_belt(:, 2) = joint_angle_trajectories_belt(:, 2) + belt_position_trajectory_mocap;
    
    % set irrelevant data points to NaN
    number_of_time_steps_mocap = length(time_mocap);
    if process_all_data
        data_points_to_process_mocap = 1 : length(time_mocap);
%         data_points_to_process_mocap = 500 : 600;
    end
    all_data_points = 1 : number_of_time_steps;
    irrelevant_data_points = ~ismember(all_data_points, data_points_to_process_mocap);
    joint_angle_trajectories(irrelevant_data_points, :) = NaN;
    joint_angle_trajectories_belt(irrelevant_data_points, :) = NaN;
    
    %% calculate angle derivatives
    joint_angles_filtered = zeros(size(joint_angle_trajectories)) * NaN;
    joint_angles_unfiltered = zeros(size(joint_angle_trajectories)) * NaN;
    joint_velocities_from_filtered = zeros(size(joint_angle_trajectories)) * NaN;
    joint_velocities_from_unfiltered = zeros(size(joint_angle_trajectories)) * NaN;
    joint_accelerations_from_filtered = zeros(size(joint_angle_trajectories)) * NaN;
    joint_accelerations_from_unfiltered = zeros(size(joint_angle_trajectories)) * NaN;
    
    joint_angles_belt_filtered = zeros(size(joint_angle_trajectories_belt)) * NaN;
    joint_angles_belt_unfiltered = zeros(size(joint_angle_trajectories_belt)) * NaN;
    joint_velocities_belt_from_filtered = zeros(size(joint_angle_trajectories_belt)) * NaN;
    joint_velocities_belt_from_unfiltered = zeros(size(joint_angle_trajectories_belt)) * NaN;
    joint_accelerations_belt_from_filtered = zeros(size(joint_angle_trajectories_belt)) * NaN;
    joint_accelerations_belt_from_unfiltered = zeros(size(joint_angle_trajectories_belt)) * NaN;
    
    % find the gaps
    gap_start_indices = [];
    gap_end_indices = [];
    for i_time = 1 : number_of_time_steps
        if isnan(joint_angle_trajectories_belt(i_time, 1))
            % check if this is the start of a gap
            if (i_time == 1) || (~isnan(joint_angle_trajectories_belt(i_time-1, 1)))
                gap_start_indices = [gap_start_indices; i_time]; %#ok<*AGROW>
            end
            % check if this is the end of a gap
            if (i_time == number_of_time_steps) || (~isnan(joint_angle_trajectories_belt(i_time+1, 1)))
                gap_end_indices = [gap_end_indices; i_time];
            end

        end
    end
    gaps = [];
    for i_gap = 1 : length(gap_start_indices)
        gaps = [gaps; [gap_start_indices(i_gap) gap_end_indices(i_gap)]];
    end

    % find the pieces of data without gaps
    data_stretches_without_gaps = [];
    if isempty(gaps)
        data_stretches_without_gaps = [1 number_of_time_steps];
    else
        if  gaps(1, 1) > 1
            data_stretches_without_gaps = [1 gaps(1, 1)-1];
        end
        for i_gap = 1 : size(gaps, 1) - 1
            % add the stretch after this big gaps to the list
            data_stretches_without_gaps = [data_stretches_without_gaps; gaps(i_gap, 2)+1 gaps(i_gap+1, 1)-1];
        end
        if ~isempty(gaps) && gaps(end, 2) < number_of_time_steps
            data_stretches_without_gaps = [data_stretches_without_gaps; gaps(end, 2)+1 number_of_time_steps];
        end
    end    
    
    % filter and calculate time derivatives for each stretch
    for i_stretch = 1 : size(data_stretches_without_gaps, 1)
        joint_angle_data_stretch = joint_angle_trajectories(data_stretches_without_gaps(i_stretch, 1) : data_stretches_without_gaps(i_stretch, 2), :);
        joint_angle_data_belt_stretch = joint_angle_trajectories_belt(data_stretches_without_gaps(i_stretch, 1) : data_stretches_without_gaps(i_stretch, 2), :);
        if length(joint_angle_data_stretch) > 3*filter_order
            stretch_data_points = data_stretches_without_gaps(i_stretch, 1) : data_stretches_without_gaps(i_stretch, 2);
            for i_joint = 1 : number_of_joints
                % filter data from this stretch and calculate derivatives
                joint_angles_stretch_filtered = filtfilt(b_filter, a_filter, joint_angle_data_stretch(:, i_joint));
                joint_angles_stretch_unfiltered = joint_angle_data_stretch(:, i_joint);
                joint_velocities_stretch_from_filtered = centdiff(joint_angles_stretch_filtered, 1/sampling_rate_mocap);
                joint_velocities_stretch_from_unfiltered = centdiff(joint_angles_stretch_unfiltered, 1/sampling_rate_mocap);
                joint_accelerations_stretch_from_filtered = centdiff(joint_velocities_stretch_from_filtered, 1/sampling_rate_mocap);
                joint_accelerations_stretch_from_unfiltered = centdiff(joint_velocities_stretch_from_unfiltered, 1/sampling_rate_mocap);
                
                joint_angles_belt_stretch_filtered = filtfilt(b_filter, a_filter, joint_angle_data_belt_stretch(:, i_joint));
                joint_angles_belt_stretch_unfiltered = joint_angle_data_belt_stretch(:, i_joint);
                joint_velocities_belt_stretch_from_filtered = centdiff(joint_angles_belt_stretch_filtered, 1/sampling_rate_mocap);
                joint_velocities_belt_stretch_from_unfiltered = centdiff(joint_angles_belt_stretch_unfiltered, 1/sampling_rate_mocap);
                joint_accelerations_belt_stretch_from_filtered = centdiff(joint_velocities_belt_stretch_from_filtered, 1/sampling_rate_mocap);
                joint_accelerations_belt_stretch_from_unfiltered = centdiff(joint_velocities_belt_stretch_from_unfiltered, 1/sampling_rate_mocap);
                
                % insert stretch data into trajectories
                joint_angles_filtered(stretch_data_points, i_joint) = joint_angles_stretch_filtered;
                joint_angles_unfiltered(stretch_data_points, i_joint) = joint_angles_stretch_unfiltered;
                joint_velocities_from_filtered(stretch_data_points, i_joint) = joint_velocities_stretch_from_filtered;
                joint_velocities_from_unfiltered(stretch_data_points, i_joint) = joint_velocities_stretch_from_unfiltered;
                joint_accelerations_from_filtered(stretch_data_points, i_joint) = joint_accelerations_stretch_from_filtered;
                joint_accelerations_from_unfiltered(stretch_data_points, i_joint) = joint_accelerations_stretch_from_unfiltered;
                
                joint_angles_belt_filtered(stretch_data_points, i_joint) = joint_angles_belt_stretch_filtered;
                joint_angles_belt_unfiltered(stretch_data_points, i_joint) = joint_angles_belt_stretch_unfiltered;
                joint_velocities_belt_from_filtered(stretch_data_points, i_joint) = joint_velocities_belt_stretch_from_filtered;
                joint_velocities_belt_from_unfiltered(stretch_data_points, i_joint) = joint_velocities_belt_stretch_from_unfiltered;
                joint_accelerations_belt_from_filtered(stretch_data_points, i_joint) = joint_accelerations_belt_stretch_from_filtered;
                joint_accelerations_belt_from_unfiltered(stretch_data_points, i_joint) = joint_accelerations_belt_stretch_from_unfiltered;
            end
        end
    end
    
    % choose data
    if use_filtered_data
        joint_angle_trajectories = joint_angles_filtered;
        joint_velocity_trajectories = joint_velocities_from_filtered;
        joint_acceleration_trajectories = joint_accelerations_from_filtered;
        
        joint_angle_trajectories_belt = joint_angles_belt_filtered;
        joint_velocity_trajectories_belt = joint_velocities_belt_from_filtered;
        joint_acceleration_trajectories_belt = joint_accelerations_belt_from_filtered;
    else
        joint_angle_trajectories = joint_angles_unfiltered;
        joint_velocity_trajectories = joint_velocities_from_unfiltered;
        joint_acceleration_trajectories = joint_accelerations_from_unfiltered;
        
        joint_angle_trajectories_belt = joint_angles_belt_unfiltered;
        joint_velocity_trajectories_belt = joint_velocities_belt_from_unfiltered;
        joint_acceleration_trajectories_belt = joint_accelerations_belt_from_unfiltered;
    end    
    
    
    % visualize
    if visualize_derivatives
        joints_to_visualize = 2;
        angle_derivative_axes = zeros(3, length(joints_to_visualize));
        for i_joint = joints_to_visualize
            figure; angle_derivative_axes(1, joints_to_visualize==i_joint) = axes; hold on; title([plant.jointLabels{i_joint} ' - angle'])
            plot(time_mocap, joint_angles_unfiltered(:, i_joint), 'x-');
            plot(time_mocap, joint_angles_filtered(:, i_joint));
            legend('unfiltered', 'filtered')

            figure; angle_derivative_axes(2, joints_to_visualize==i_joint) = axes; hold on; title([plant.jointLabels{i_joint} ' - velocity'])
            plot(time_mocap, joint_velocities_from_unfiltered(:, i_joint), 'x-');
            plot(time_mocap, joint_velocities_from_filtered(:, i_joint));
            legend('from unfiltered', 'from filtered')

            figure; angle_derivative_axes(3, joints_to_visualize==i_joint) = axes; hold on; title([plant.jointLabels{i_joint} ' - acceleration'])
            plot(time_mocap, joint_accelerations_from_unfiltered(:, i_joint), 'x-');
            plot(time_mocap, joint_accelerations_from_filtered(:, i_joint));
            legend('from unfiltered', 'from filtered')
        end
        linkaxes(angle_derivative_axes, 'x');
        distFig('cols', length(joints_to_visualize))
    end    
    
    %% calculate model-dependent variables - world coordinates
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
    V_body_left_ankle = zeros(number_of_time_steps, 6);
    V_body_right_ankle = zeros(number_of_time_steps, 6);
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
        V_body_left_ankle_pool = zeros(number_of_time_steps, 6);
        V_body_right_ankle_pool = zeros(number_of_time_steps, 6);
        com_trajectory_pool = zeros(number_of_time_steps, 3);


        % create a copy of the plant for each worker
        spmd
            plant_pool = plant.copy;
            for i_time = labindex : numlabs : number_of_time_steps
                joint_angles = joint_angle_trajectories(i_time, :)';
                joint_velocities = joint_velocity_trajectories(i_time, :)';
                if ~ismember(i_time, data_points_to_process_mocap) || any(isnan(joint_angles))
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
                    V_body_left_ankle_pool(i_time, :) = NaN;
                    V_body_right_ankle_pool(i_time, :) = NaN;
                    com_trajectory_pool(i_time, :) = NaN;
                else
                    plant_pool.jointAngles = joint_angles;
                    plant_pool.jointVelocities = joint_velocities;
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
                    V_body_left_ankle_pool(i_time, :) = plant.bodyJacobians{3} * joint_velocities;
                    V_body_right_ankle_pool(i_time, :) = plant.bodyJacobians{6} * joint_velocities;
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
            V_body_left_ankle_lab = V_body_left_ankle_pool{i_lab};
            V_body_right_ankle_lab = V_body_right_ankle_pool{i_lab};
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
            V_body_left_ankle(i_lab : number_of_labs : number_of_time_steps, :) = V_body_left_ankle_lab(i_lab : number_of_labs : number_of_time_steps, :);
            V_body_right_ankle(i_lab : number_of_labs : number_of_time_steps, :) = V_body_right_ankle_lab(i_lab : number_of_labs : number_of_time_steps, :);
            com_trajectory(i_lab : number_of_labs : number_of_time_steps, :) = com_trajectory_lab(i_lab : number_of_labs : number_of_time_steps, :);
        end

    else
        for i_time = 1 : number_of_time_steps
            joint_angles = joint_angle_trajectories(i_time, :)';
            joint_velocities = joint_velocity_trajectories(i_time, :)';
            if ~ismember(i_time, data_points_to_process_mocap) || any(isnan(joint_angles))
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
                V_body_left_ankle(i_time, :) = NaN;
                V_body_right_ankle(i_time, :) = NaN;
                com_trajectory(i_time, :) = NaN;
            else
                plant.jointAngles = joint_angles;
                plant.jointVelocities = joint_velocities;
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
                V_body_left_ankle(i_time, :) = plant.bodyJacobians{3} * joint_velocities;
                V_body_right_ankle(i_time, :) = plant.bodyJacobians{6} * joint_velocities;
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

    %% calculate model-dependent variables - belt coordinates
    marker_trajectories_reconstructed_belt = zeros(size(marker_trajectories));
    right_heel_trajectory_belt = zeros(number_of_time_steps, 3);
    right_toes_trajectory_belt = zeros(number_of_time_steps, 3);
    left_heel_trajectory_belt = zeros(number_of_time_steps, 3);
    left_toes_trajectory_belt = zeros(number_of_time_steps, 3);
    pelvis_roll_angle_belt = zeros(number_of_time_steps, 1);
    trunk_roll_angle_belt = zeros(number_of_time_steps, 1);
    head_roll_angle_belt = zeros(number_of_time_steps, 1);
    right_foot_roll_angle_belt = zeros(number_of_time_steps, 1);
    left_foot_roll_angle_belt = zeros(number_of_time_steps, 1);
    T_left_ankle_to_world_trajectory_belt = zeros(number_of_time_steps, 16);
    T_right_ankle_to_world_trajectory_belt = zeros(number_of_time_steps, 16);
    V_body_left_ankle_belt = zeros(number_of_time_steps, 6);
    V_body_right_ankle_belt = zeros(number_of_time_steps, 6);
    com_trajectory_belt = zeros(number_of_time_steps, 3);
    
   if use_parallel
        % get or open pool of workers
        marker_trajectories_belt_pool = zeros(size(marker_trajectories_reconstructed_belt));
        left_heel_trajectory_belt_pool = zeros(number_of_time_steps, 3);
        left_toes_trajectory_belt_pool = zeros(number_of_time_steps, 3);
        right_heel_trajectory_belt_pool = zeros(number_of_time_steps, 3);
        right_toes_trajectory_belt_pool = zeros(number_of_time_steps, 3);
        pelvis_roll_angle_trajectory_belt_pool = zeros(number_of_time_steps, 1);
        trunk_roll_angle_trajectory_belt_pool = zeros(number_of_time_steps, 1);
        head_roll_angle_trajectory_belt_pool = zeros(number_of_time_steps, 1);
        left_foot_roll_angle_trajectory_belt_pool = zeros(number_of_time_steps, 1);
        right_foot_roll_angle_trajectory_belt_pool = zeros(number_of_time_steps, 1);
        T_left_ankle_to_world_trajectory_belt_pool = zeros(number_of_time_steps, 16);
        T_right_ankle_to_world_trajectory_belt_pool = zeros(number_of_time_steps, 16);
        V_body_left_ankle_belt_pool = zeros(number_of_time_steps, 6);
        V_body_right_ankle_belt_pool = zeros(number_of_time_steps, 6);
        com_trajectory_belt_pool = zeros(number_of_time_steps, 3);


        % create a copy of the plant for each worker
        spmd
            plant_pool = plant.copy;
            for i_time = labindex : numlabs : number_of_time_steps
                joint_angles = joint_angle_trajectories_belt(i_time, :)';
                joint_velocities = joint_velocity_trajectories_belt(i_time, :)';
                if any(isnan(joint_angles))
                    marker_trajectories_belt_pool(i_time, :) = NaN;
                    left_heel_trajectory_belt_pool(i_time, :) = NaN;
                    left_toes_trajectory_belt_pool(i_time, :) = NaN;
                    right_heel_trajectory_belt_pool(i_time, :) = NaN;
                    right_toes_trajectory_belt_pool(i_time, :) = NaN;
                    pelvis_roll_angle_trajectory_belt_pool(i_time, :) = NaN;
                    trunk_roll_angle_trajectory_belt_pool(i_time, :) = NaN;
                    head_roll_angle_trajectory_belt_pool(i_time, :) = NaN;
                    left_foot_roll_angle_trajectory_belt_pool(i_time, :) = NaN;
                    right_foot_roll_angle_trajectory_belt_pool(i_time, :) = NaN;
                    T_left_ankle_to_world_trajectory_belt_pool(i_time, :) = NaN;
                    T_right_ankle_to_world_trajectory_belt_pool(i_time, :) = NaN;
                    V_body_left_ankle_belt_pool(i_time, :) = NaN;
                    V_body_right_ankle_belt_pool(i_time, :) = NaN;
                    com_trajectory_belt_pool(i_time, :) = NaN;
                else
                    plant_pool.jointAngles = joint_angles;
                    plant_pool.jointVelocities = joint_velocities;
                    plant_pool.updateKinematics();
                    marker_trajectories_belt_pool(i_time, :) = plant_pool.exportMarkerPositions();
                    left_heel_trajectory_belt_pool(i_time, :) = plant_pool.endEffectorPositions{1};
                    left_toes_trajectory_belt_pool(i_time, :) = plant_pool.endEffectorPositions{2};
                    right_heel_trajectory_belt_pool(i_time, :) = plant_pool.endEffectorPositions{4};
                    right_toes_trajectory_belt_pool(i_time, :) = plant_pool.endEffectorPositions{5};
                    pelvis_roll_angle_trajectory_belt_pool(i_time, :) = -atan2(plant_pool.jointTransformations{6}(1, 3), plant_pool.jointTransformations{6}(3, 3));
                    trunk_roll_angle_trajectory_belt_pool(i_time, :) = -atan2(plant_pool.jointTransformations{21}(1, 3), plant_pool.jointTransformations{21}(3, 3));
                    head_roll_angle_trajectory_belt_pool(i_time, :) = -atan2(plant_pool.jointTransformations{24}(1, 3), plant_pool.jointTransformations{24}(3, 3));
                    left_foot_roll_angle_trajectory_belt_pool(i_time, :) = -atan2(plant_pool.jointTransformations{12}(1, 2), plant_pool.jointTransformations{12}(3, 2));
                    right_foot_roll_angle_trajectory_belt_pool(i_time, :) = -atan2(plant_pool.jointTransformations{18}(1, 2), plant_pool.jointTransformations{18}(3, 2));
                    T_left_ankle_to_world_trajectory_belt_pool(i_time, :) = reshape(plant_pool.endEffectorTransformations{3}, 1, 16);
                    T_right_ankle_to_world_trajectory_belt_pool(i_time, :) = reshape(plant_pool.endEffectorTransformations{6}, 1, 16);
                    V_body_left_ankle_belt_pool(i_time, :) = plant.bodyJacobians{3} * joint_velocities;
                    V_body_right_ankle_belt_pool(i_time, :) = plant.bodyJacobians{6} * joint_velocities;
                    com_trajectory_belt_pool(i_time, :) = plant_pool.calculateCenterOfMassPosition();
                end

            end
        end
        for i_lab = 1 : number_of_labs
            marker_trajectories_belt_lab = marker_trajectories_belt_pool{i_lab};
            left_heel_trajectory_belt_lab = left_heel_trajectory_belt_pool{i_lab};
            left_toes_trajectory_belt_lab = left_toes_trajectory_belt_pool{i_lab};
            right_heel_trajectory_belt_lab = right_heel_trajectory_belt_pool{i_lab};
            right_toes_trajectory_belt_lab = right_toes_trajectory_belt_pool{i_lab};
            pelvis_roll_angle_trajectory_belt_lab = pelvis_roll_angle_trajectory_belt_pool{i_lab};
            trunk_roll_angle_trajectory_belt_lab = trunk_roll_angle_trajectory_belt_pool{i_lab};
            head_roll_angle_trajectory_belt_lab = head_roll_angle_trajectory_belt_pool{i_lab};
            left_foot_roll_angle_trajectory_belt_lab = left_foot_roll_angle_trajectory_belt_pool{i_lab};
            right_foot_roll_angle_trajectory_belt_lab = right_foot_roll_angle_trajectory_belt_pool{i_lab};
            T_left_ankle_to_world_trajectory_belt_lab = T_left_ankle_to_world_trajectory_belt_pool{i_lab};
            T_right_ankle_to_world_trajectory_belt_lab = T_right_ankle_to_world_trajectory_belt_pool{i_lab};
            V_body_left_ankle_belt_lab = V_body_left_ankle_belt_pool{i_lab};
            V_body_right_ankle_belt_lab = V_body_right_ankle_belt_pool{i_lab};
            com_trajectory_belt_lab = com_trajectory_belt_pool{i_lab};

            marker_trajectories_reconstructed_belt(i_lab : number_of_labs : number_of_time_steps, :) = marker_trajectories_belt_lab(i_lab : number_of_labs : number_of_time_steps, :);
            left_heel_trajectory_belt(i_lab : number_of_labs : number_of_time_steps, :) = left_heel_trajectory_belt_lab(i_lab : number_of_labs : number_of_time_steps, :);
            left_toes_trajectory_belt(i_lab : number_of_labs : number_of_time_steps, :) = left_toes_trajectory_belt_lab(i_lab : number_of_labs : number_of_time_steps, :);
            right_heel_trajectory_belt(i_lab : number_of_labs : number_of_time_steps, :) = right_heel_trajectory_belt_lab(i_lab : number_of_labs : number_of_time_steps, :);
            right_toes_trajectory_belt(i_lab : number_of_labs : number_of_time_steps, :) = right_toes_trajectory_belt_lab(i_lab : number_of_labs : number_of_time_steps, :);
            pelvis_roll_angle_belt(i_lab : number_of_labs : number_of_time_steps, :) = pelvis_roll_angle_trajectory_belt_lab(i_lab : number_of_labs : number_of_time_steps, :);
            trunk_roll_angle_belt(i_lab : number_of_labs : number_of_time_steps, :) = trunk_roll_angle_trajectory_belt_lab(i_lab : number_of_labs : number_of_time_steps, :);
            head_roll_angle_belt(i_lab : number_of_labs : number_of_time_steps, :) = head_roll_angle_trajectory_belt_lab(i_lab : number_of_labs : number_of_time_steps, :);
            left_foot_roll_angle_belt(i_lab : number_of_labs : number_of_time_steps, :) = left_foot_roll_angle_trajectory_belt_lab(i_lab : number_of_labs : number_of_time_steps, :);
            right_foot_roll_angle_belt(i_lab : number_of_labs : number_of_time_steps, :) = right_foot_roll_angle_trajectory_belt_lab(i_lab : number_of_labs : number_of_time_steps, :);
            T_left_ankle_to_world_trajectory_belt(i_lab : number_of_labs : number_of_time_steps, :) = T_left_ankle_to_world_trajectory_belt_lab(i_lab : number_of_labs : number_of_time_steps, :);
            T_right_ankle_to_world_trajectory_belt(i_lab : number_of_labs : number_of_time_steps, :) = T_right_ankle_to_world_trajectory_belt_lab(i_lab : number_of_labs : number_of_time_steps, :);
            V_body_left_ankle_belt(i_lab : number_of_labs : number_of_time_steps, :) = V_body_left_ankle_belt_lab(i_lab : number_of_labs : number_of_time_steps, :);
            V_body_right_ankle_belt(i_lab : number_of_labs : number_of_time_steps, :) = V_body_right_ankle_belt_lab(i_lab : number_of_labs : number_of_time_steps, :);
            com_trajectory_belt(i_lab : number_of_labs : number_of_time_steps, :) = com_trajectory_belt_lab(i_lab : number_of_labs : number_of_time_steps, :);
        end

    else
        for i_time = 1 : number_of_time_steps
            joint_angles = joint_angle_trajectories_belt(i_time, :)';
            joint_velocities = joint_velocity_trajectories(i_time, :)';
            if any(isnan(joint_angles))
                marker_trajectories_reconstructed_belt(i_time, :) = NaN;
                left_heel_trajectory_belt(i_time, :) = NaN;
                left_toes_trajectory_belt(i_time, :) = NaN;
                right_heel_trajectory_belt(i_time, :) = NaN;
                right_toes_trajectory_belt(i_time, :) = NaN;
                pelvis_roll_angle_belt(i_time, :) = NaN;
                trunk_roll_angle_belt(i_time, :) = NaN;
                head_roll_angle_belt(i_time, :) = NaN;
                left_foot_roll_angle_belt(i_time, :) = NaN;
                right_foot_roll_angle_belt(i_time, :) = NaN;
                T_left_ankle_to_world_trajectory_belt(i_time, :) = NaN;
                T_right_ankle_to_world_trajectory_belt(i_time, :) = NaN;
                V_body_left_ankle_belt(i_time, :) = NaN;
                V_body_right_ankle_belt(i_time, :) = NaN;
                com_trajectory_belt(i_time, :) = NaN;
            else
                plant.jointAngles = joint_angles;
                plant.jointVelocities = joint_velocities;
                plant.updateKinematics();
                marker_trajectories_reconstructed_belt(i_time, :) = plant.exportMarkerPositions();
                left_heel_trajectory_belt(i_time, :) = plant.endEffectorPositions{1};
                left_toes_trajectory_belt(i_time, :) = plant.endEffectorPositions{2};
                right_heel_trajectory_belt(i_time, :) = plant.endEffectorPositions{4};
                right_toes_trajectory_belt(i_time, :) = plant.endEffectorPositions{5};
                pelvis_roll_angle_belt(i_time, :) = atan2(plant.jointTransformations{6}(2, 3), plant.jointTransformations{6}(3, 3));
                trunk_roll_angle_belt(i_time, :) = atan2(plant.jointTransformations{21}(2, 3), plant.jointTransformations{21}(3, 3));
                head_roll_angle_belt(i_time, :) = atan2(plant.jointTransformations{24}(2, 3), plant.jointTransformations{24}(3, 3));
                right_foot_roll_angle_belt(i_time, :) = atan2(plant.jointTransformations{18}(2, 2), plant.jointTransformations{18}(3, 2));
                left_foot_roll_angle_belt(i_time, :) = atan2(plant.jointTransformations{12}(2, 2), plant.jointTransformations{12}(3, 2));
                T_left_ankle_to_world_trajectory_belt(i_time, :) = reshape(plant.endEffectorTransformations{3}, 1, 16);
                T_right_ankle_to_world_trajectory_belt(i_time, :) = reshape(plant.endEffectorTransformations{6}, 1, 16);
                V_body_left_ankle_belt(i_time, :) = plant.bodyJacobians{3} * joint_velocities;
                V_body_right_ankle_belt(i_time, :) = plant.bodyJacobians{6} * joint_velocities;
                com_trajectory_belt(i_time, :) = plant.calculateCenterOfMassPosition();
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
    marker_reconstruction_error_belt = zeros(number_of_time_steps, number_of_markers);
    for i_marker = 1 : number_of_markers
        for i_time = 1 : number_of_time_steps
            marker_position_measured = marker_trajectories(i_time, 3*(i_marker-1)+1 : 3*(i_marker-1)+3);
            marker_position_reconstr = marker_trajectories_reconstructed_belt(i_time, 3*(i_marker-1)+1 : 3*(i_marker-1)+3);
            error_vector = marker_position_measured - marker_position_reconstr;
            marker_reconstruction_error_belt(i_time, i_marker) = norm(error_vector);
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
        'V_body_left_ankle', ...
        'V_body_right_ankle', ...
        'com_trajectory', ...
        'joint_angle_trajectories', ...
        'joint_velocity_trajectories', ...
        'joint_acceleration_trajectories', ...
        'belt_position_trajectory_mocap', ...
        'belt_position_trajectory_forceplate', ...
        'marker_trajectories_reconstructed_belt', ...
        'marker_reconstruction_error_belt', ...
        'right_heel_trajectory_belt', ...
        'right_toes_trajectory_belt', ...
        'left_heel_trajectory_belt', ...
        'left_toes_trajectory_belt', ...
        'pelvis_roll_angle_belt', ...
        'trunk_roll_angle_belt', ...
        'head_roll_angle_belt', ...
        'right_foot_roll_angle_belt', ...
        'left_foot_roll_angle_belt', ...
        'T_left_ankle_to_world_trajectory_belt', ...
        'T_right_ankle_to_world_trajectory_belt', ...
        'V_body_left_ankle_belt', ...
        'V_body_right_ankle_belt', ...
        'com_trajectory_belt', ...
        'joint_angle_trajectories_belt', ...
        'joint_velocity_trajectories_belt', ...
        'joint_acceleration_trajectories_belt' ...
      );

    
    disp(['Trial ' num2str(i_trial) ' completed']);
    toc
end
fprintf(' done\n');



