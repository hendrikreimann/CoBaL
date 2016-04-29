% calculate ground reaction wrenches

use_parallel = 1;
use_smoothed = 1;

use_point_constraints           = 0;
use_hinge_constraints           = 1;
use_body_velocity_constraints   = 0;

plot_ground_reaction_wrenches   = 1;
save_results                    = 1;

trials_to_process = 2;
trials_to_exclude = [];

% load data
load subjectInfo.mat;
load(makeFileName(date, subject_id, 'model'));

% data_points = 1 : numberOfDataPoints;
data_points = 1000 : 6000;
% data_points = 1000 : 1050;
% data_points = 1000 : 11000;

if use_parallel
    poolobject = gcp;
    number_of_labs = poolobject.NumWorkers;
end

for i_trial = trials_to_process
    if ~ismember(i_trial, trials_to_exclude)
        tic
        
        load(makeFileName(date, subject_id, 'walking', i_trial, 'markerTrajectories'));
        load(makeFileName(date, subject_id, 'walking', i_trial, 'angleTrajectories'));
        load(makeFileName(date, subject_id, 'walking', i_trial, 'kinematicTrajectories'));
        load(makeFileName(date, subject_id, 'walking', i_trial, 'stepEvents'));
        load(makeFileName(date, subject_id, 'walking', i_trial, 'forcePlateData'));
        number_of_time_steps = size(joint_angle_trajectories_belt, 1);
        
        label = 'inverseDynamics';
        if use_point_constraints
            label = [label '_pointConstraints'];
        end            
        if use_hinge_constraints
            label = [label '_hingeConstraints'];
        end            
        if use_body_velocity_constraints
            label = [label '_bodyVelocityConstraints'];
        end            
        load(makeFileName(date, subject_id, 'walking', i_trial, label));

        if use_smoothed
            constraint_torque_trajectories_left_current = constraint_torque_trajectories_left_smoothed;
            constraint_torque_trajectories_right_current = constraint_torque_trajectories_right_smoothed;
        else
            constraint_torque_trajectories_left_current = constraint_torque_trajectories_left;
            constraint_torque_trajectories_right_current = constraint_torque_trajectories_right;
        end
        
        right_ground_reaction_wrench_trajectory = zeros(numberOfDataPoints, 6);
        left_ground_reaction_wrench_trajectory = zeros(numberOfDataPoints, 6);

        if use_parallel
            right_ground_reaction_wrench_trajectory_pool = zeros(numberOfDataPoints, 6);
            left_ground_reaction_wrench_trajectory_pool = zeros(numberOfDataPoints, 6);
            spmd
                plant_pool = plant.copy;
                for i_time = data_points(1)+labindex-1 : numlabs : data_points(end)
                    if any(isnan(joint_angle_trajectories_belt(i_time, :)))
                        right_ground_reaction_wrench_trajectory_pool(i_time, :) = NaN;
                        left_ground_reaction_wrench_trajectory_pool(i_time, :) = NaN;
                    else
                        % update model
                        plant_pool.jointAngles = joint_angle_trajectories_belt(i_time, :)';
                        plant_pool.jointVelocities = joint_velocity_trajectories_belt(i_time, :)';
                        plant_pool.jointAccelerations = joint_acceleration_trajectories_belt(i_time, :)';
                        plant_pool.updateKinematics;

                        % calculate ground reaction wrenches
                        treadmill_origin_belt = [0; joint_angle_trajectories_belt(i_time, 2); 0]; % this is the physical location of the world coordinate frame in belt coordinates
                        world_to_belt_transformation = [eye(3), treadmill_origin_belt; 0 0 0 1];
                        left_ground_reaction_wrench_trajectory_pool(i_time, :) = calculateInstantaneousGroundReactionWrench(plant_pool, constraint_torque_trajectories_left_current(i_time, :)', world_to_belt_transformation);
                        right_ground_reaction_wrench_trajectory_pool(i_time, :) = calculateInstantaneousGroundReactionWrench(plant_pool, constraint_torque_trajectories_right_current(i_time, :)', world_to_belt_transformation);
                    end
                end
            end

            % reassemble
            for i_lab = 1 : number_of_labs
                left_ground_reaction_wrench_trajectory_lab = left_ground_reaction_wrench_trajectory_pool{i_lab};
                right_ground_reaction_wrench_trajectory_lab = right_ground_reaction_wrench_trajectory_pool{i_lab};

                left_ground_reaction_wrench_trajectory(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
                    = left_ground_reaction_wrench_trajectory_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
                right_ground_reaction_wrench_trajectory(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
                    = right_ground_reaction_wrench_trajectory_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
            end                
        else
            for i_time = data_points
                if any(isnan(joint_angle_trajectories_belt(i_time, :)))
                    right_ground_reaction_wrench_trajectory(i_time, :) = NaN;
                    left_ground_reaction_wrench_trajectory(i_time, :) = NaN;
                else
                    % TODO: transform the plant back to world coordinates from belt coordinates, so the ground reaction wrench will coincide with the force plate readings 
                    % update model
                    plant.jointAngles = joint_angle_trajectories_belt(i_time, :)';
                    plant.jointVelocities = joint_velocity_trajectories_belt(i_time, :)';
                    plant.jointAccelerations = joint_acceleration_trajectories_belt(i_time, :)';
                    plant.updateKinematics;

                    % calculate ground reaction wrenches
                    treadmill_origin_belt = [0; joint_angle_trajectories_belt(i_time, 2); 0]; % this is the physical location of the world coordinate frame in belt coordinates
                    world_to_belt_transformation = [eye(3), treadmill_origin_belt; 0 0 0 1];
                    left_ground_reaction_wrench_trajectory(i_time, :) = calculateInstantaneousGroundReactionWrench(plant, constraint_torque_trajectories_left_current(i_time, :)', world_to_belt_transformation);
                    right_ground_reaction_wrench_trajectory(i_time, :) = calculateInstantaneousGroundReactionWrench(plant, constraint_torque_trajectories_right_current(i_time, :)', world_to_belt_transformation);
                    
                end
            end
        end
        toc
        
        
        %% save data
        if save_results
            label = 'groundReactionWrenches';
            if use_point_constraints
                label = [label '_pointConstraints'];
            end            
            if use_hinge_constraints
                label = [label '_hingeConstraints'];
            end            
            if use_body_velocity_constraints
                label = [label '_bodyVelocityConstraints'];
            end            

            save( ...
                  makeFileName(date, subject_id, 'walking', i_trial, label), ...
                  'left_ground_reaction_wrench_trajectory', ...
                  'right_ground_reaction_wrench_trajectory' ...
                );
        end
        
        
        %% plot ground reaction wrenches
        if plot_ground_reaction_wrenches
            color_x_a = [1 0 0];
            color_y_a = [0 1 0];
            color_z_a = [0 0 1];
            color_x_b = [0.7 0 0.3];
            color_y_b = [0.3 0.7 0];
            color_z_b = [0 0.3 0.7];
            
            figure; left_grf_axes = axes; hold on; title('left ground reaction forces');
            plot(time_force_plate, left_force_plate_wrench_Acw(:, 1), 'color', color_x_a);
            plot(time_force_plate, left_force_plate_wrench_Acw(:, 2), 'color', color_y_a);
            plot(time_force_plate, left_force_plate_wrench_Acw(:, 3), 'color', color_z_a);
            plot(time_mocap, left_ground_reaction_wrench_trajectory(:, 1), 'color', color_x_b);
            plot(time_mocap, left_ground_reaction_wrench_trajectory(:, 2), 'color', color_y_b);
            plot(time_mocap, left_ground_reaction_wrench_trajectory(:, 3), 'color', color_z_b);
            
            figure; left_grm_axes = axes; hold on; title('left ground reaction moments');
            plot(time_force_plate, left_force_plate_wrench_Acw(:, 4), 'color', color_x_a);
            plot(time_force_plate, left_force_plate_wrench_Acw(:, 5), 'color', color_y_a);
            plot(time_force_plate, left_force_plate_wrench_Acw(:, 6), 'color', color_z_a);
            plot(time_mocap, left_ground_reaction_wrench_trajectory(:, 4), 'color', color_x_b);
            plot(time_mocap, left_ground_reaction_wrench_trajectory(:, 5), 'color', color_y_b);
            plot(time_mocap, left_ground_reaction_wrench_trajectory(:, 6), 'color', color_z_b);
            
            figure; right_grf_axes = axes; hold on; title('right ground reaction forces');
            plot(time_force_plate, right_force_plate_wrench_Acw(:, 1), 'color', color_x_a);
            plot(time_force_plate, right_force_plate_wrench_Acw(:, 2), 'color', color_y_a);
            plot(time_force_plate, right_force_plate_wrench_Acw(:, 3), 'color', color_z_a);
            plot(time_mocap, right_ground_reaction_wrench_trajectory(:, 1), 'color', color_x_b);
            plot(time_mocap, right_ground_reaction_wrench_trajectory(:, 2), 'color', color_y_b);
            plot(time_mocap, right_ground_reaction_wrench_trajectory(:, 3), 'color', color_z_b);
            
            figure; right_grm_axes = axes; hold on; title('right ground reaction moments');
            plot(time_force_plate, right_force_plate_wrench_Acw(:, 4), 'color', color_x_a);
            plot(time_force_plate, right_force_plate_wrench_Acw(:, 5), 'color', color_y_a);
            plot(time_force_plate, right_force_plate_wrench_Acw(:, 6), 'color', color_z_a);
            plot(time_mocap, right_ground_reaction_wrench_trajectory(:, 4), 'color', color_x_b);
            plot(time_mocap, right_ground_reaction_wrench_trajectory(:, 5), 'color', color_y_b);
            plot(time_mocap, right_ground_reaction_wrench_trajectory(:, 6), 'color', color_z_b);
            
            linkaxes([left_grf_axes left_grm_axes right_grf_axes right_grm_axes], 'x')
            distFig('rows', 2, 'tight', true);
        end
    end
end



