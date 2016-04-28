% calculate ground reaction wrenches

use_parallel = 1;

use_point_constraints           = 0;
use_hinge_constraints           = 1;
use_body_velocity_constraints   = 0;

plot_ground_reaction_wrenches   = 1;

trials_to_process = 2;
trials_to_exclude = [];

% load data
load subjectInfo.mat;
load(makeFileName(date, subject_id, 'model'));

% data_points = 1 : numberOfDataPoints;
data_points = 1000 : 6000;
data_points = 1000 : 1050;

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
        number_of_time_steps = size(angle_trajectories, 1);
        
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

        right_ground_reaction_wrench_trajectory_origin = zeros(numberOfDataPoints, 6);
        left_ground_reaction_wrench_trajectory_origin = zeros(numberOfDataPoints, 6);

        if use_parallel
            right_ground_reaction_wrench_trajectory_origin_pool = zeros(numberOfDataPoints, 6);
            left_ground_reaction_wrench_trajectory_origin_pool = zeros(numberOfDataPoints, 6);
            spmd
                plant_pool = plant.copy;
                for i_time = data_points(1)+labindex-1 : numlabs : data_points(end)
                    if any(isnan(angle_trajectories(i_time, :)))
                        right_ground_reaction_wrench_trajectory_origin_pool(i_time, :) = NaN;
                        left_ground_reaction_wrench_trajectory_origin_pool(i_time, :) = NaN;
                    else
                        % update model
                        plant_pool.jointAngles = angle_trajectories(i_time, :)';
                        plant_pool.jointVelocities = joint_velocity_trajectories(i_time, :)';
                        plant_pool.jointAccelerations = joint_acceleration_trajectories(i_time, :)';
                        plant_pool.updateKinematics;

                        % calculate ground reaction wrenches
                        left_ground_reaction_wrench_trajectory_origin_pool(i_time, :) = calculateGroundReactionWrench(plant_pool, constraint_torque_trajectories_left(i_time, :)', eye(4));
                        right_ground_reaction_wrench_trajectory_origin_pool(i_time, :) = calculateGroundReactionWrench(plant_pool, constraint_torque_trajectories_right(i_time, :)', eye(4));
                    end
                end
            end

            % reassemble
            for i_lab = 1 : number_of_labs
                left_ground_reaction_wrench_trajectory_origin_lab = left_ground_reaction_wrench_trajectory_origin_pool{i_lab};
                right_ground_reaction_wrench_trajectory_origin_lab = right_ground_reaction_wrench_trajectory_origin_pool{i_lab};

                left_ground_reaction_wrench_trajectory_origin(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
                    = left_ground_reaction_wrench_trajectory_origin_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
                right_ground_reaction_wrench_trajectory_origin(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
                    = right_ground_reaction_wrench_trajectory_origin_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
            end                
        else
            for i_time = data_points
                if any(isnan(angle_trajectories(i_time, :)))
                    right_ground_reaction_wrench_trajectory_origin(i_time, :) = NaN;
                    left_ground_reaction_wrench_trajectory_origin(i_time, :) = NaN;
                else
                    % TODO: transform the plant back to world coordinates from belt coordinates, so the ground reaction wrench will coincide with the force plate readings 
                    % update model
                    plant.jointAngles = angle_trajectories(i_time, :)';
                    plant.jointVelocities = joint_velocity_trajectories(i_time, :)';
                    plant.jointAccelerations = joint_acceleration_trajectories(i_time, :)';
                    plant.updateKinematics;

                    % calculate ground reaction wrenches
                    left_ground_reaction_wrench_trajectory_origin(i_time, :) = calculateGroundReactionWrench(plant, constraint_torque_trajectories_left(i_time, :)', eye(4));
                    right_ground_reaction_wrench_trajectory_origin(i_time, :) = calculateGroundReactionWrench(plant, constraint_torque_trajectories_right(i_time, :)', eye(4));
                end
            end
        end
        toc
        
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
            plot(time_mocap, left_ground_reaction_wrench_trajectory_origin(:, 1), 'color', color_x_b);
            plot(time_mocap, left_ground_reaction_wrench_trajectory_origin(:, 2), 'color', color_y_b);
            plot(time_mocap, left_ground_reaction_wrench_trajectory_origin(:, 3), 'color', color_z_b);
            
            figure; left_grm_axes = axes; hold on; title('left ground reaction moments');
            plot(time_force_plate, left_force_plate_wrench_Acw(:, 4), 'color', color_x_a);
            plot(time_force_plate, left_force_plate_wrench_Acw(:, 5), 'color', color_y_a);
            plot(time_force_plate, left_force_plate_wrench_Acw(:, 6), 'color', color_z_a);
            plot(time_mocap, left_ground_reaction_wrench_trajectory_origin(:, 4), 'color', color_x_b);
            plot(time_mocap, left_ground_reaction_wrench_trajectory_origin(:, 5), 'color', color_y_b);
            plot(time_mocap, left_ground_reaction_wrench_trajectory_origin(:, 6), 'color', color_z_b);
            
            figure; right_grf_axes = axes; hold on; title('right ground reaction forces');
            plot(time_force_plate, right_force_plate_wrench_Acw(:, 1), 'color', color_x_a);
            plot(time_force_plate, right_force_plate_wrench_Acw(:, 2), 'color', color_y_a);
            plot(time_force_plate, right_force_plate_wrench_Acw(:, 3), 'color', color_z_a);
            plot(time_mocap, right_ground_reaction_wrench_trajectory_origin(:, 1), 'color', color_x_b);
            plot(time_mocap, right_ground_reaction_wrench_trajectory_origin(:, 2), 'color', color_y_b);
            plot(time_mocap, right_ground_reaction_wrench_trajectory_origin(:, 3), 'color', color_z_b);
            
            figure; right_grm_axes = axes; hold on; title('right ground reaction moments');
            plot(time_force_plate, right_force_plate_wrench_Acw(:, 4), 'color', color_x_a);
            plot(time_force_plate, right_force_plate_wrench_Acw(:, 5), 'color', color_y_a);
            plot(time_force_plate, right_force_plate_wrench_Acw(:, 6), 'color', color_z_a);
            plot(time_mocap, right_ground_reaction_wrench_trajectory_origin(:, 4), 'color', color_x_b);
            plot(time_mocap, right_ground_reaction_wrench_trajectory_origin(:, 5), 'color', color_y_b);
            plot(time_mocap, right_ground_reaction_wrench_trajectory_origin(:, 6), 'color', color_z_b);
            
            linkaxes([left_grf_axes left_grm_axes right_grf_axes right_grm_axes], 'x')
            distFig('rows', 2, 'tight', true);
        end        
    end
end



