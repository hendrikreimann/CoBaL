% inverse dynamics

determine_constraint_numbers            = 0;
calculate_dynamic_matrices              = 0;
calculate_torques                       = 0;
filter_torques                          = 0;

plot_joint_torques                      = 1;
plot_joint_torques_smoothed             = 0;

process_all_data                        = 1;
use_parallel                            = 0;

load_results                            = 0;
save_results                            = 1;

use_point_constraints                   = 0;
use_hinge_constraints                   = 1;
use_body_velocity_constraints           = 0;
use_rollover_constraints                = 0;


% select which trials to process
% trials_to_process = 1;

% trials_to_process = 1001 : 1180;
% trials_to_exclude = [1046 1064 1066 1069 1089 1092 1109 1114 1132];
trials_to_process = 1;
trials_to_exclude = [];



% load data
load subjectInfo.mat;
load(makeFileName(date, subject_id, 'model'));

% parameters
number_of_joints = plant.numberOfJoints;

data_points = 1 : 30000;
% data_points = 200 : 900;
% data_points = 450 : 480;
% data_points = 1000 : 11000;
data_points = 1000 : 1050;

if use_parallel && (calculate_dynamic_matrices || calculate_torques)
    poolobject = gcp;
    number_of_labs = poolobject.NumWorkers;
end

if use_body_velocity_constraints
    load(makeFileName(date, subject_id, 'constraints'));
    V_body_left_fits_heelstrike = ...
      { ...
        V_body_left_1_fit_heelstrike, ...
        V_body_left_2_fit_heelstrike, ...
        V_body_left_3_fit_heelstrike, ...
        V_body_left_4_fit_heelstrike, ...
        V_body_left_5_fit_heelstrike, ...
        V_body_left_6_fit_heelstrike ...
      };
    V_body_left_fits_pushoff = ...
      { ...
        V_body_left_1_fit_pushoff, ...
        V_body_left_2_fit_pushoff, ...
        V_body_left_3_fit_pushoff, ...
        V_body_left_4_fit_pushoff, ...
        V_body_left_5_fit_pushoff, ...
        V_body_left_6_fit_pushoff ...
      };
    V_body_right_fits_heelstrike = ...
      { ...
        V_body_right_1_fit_heelstrike, ...
        V_body_right_2_fit_heelstrike, ...
        V_body_right_3_fit_heelstrike, ...
        V_body_right_4_fit_heelstrike, ...
        V_body_right_5_fit_heelstrike, ...
        V_body_right_6_fit_heelstrike ...
      };
    V_body_right_fits_pushoff = ...
      { ...
        V_body_right_1_fit_pushoff, ...
        V_body_right_2_fit_pushoff, ...
        V_body_right_3_fit_pushoff, ...
        V_body_right_4_fit_pushoff, ...
        V_body_right_5_fit_pushoff, ...
        V_body_right_6_fit_pushoff ...
      };
  
end
if use_rollover_constraints
    load(makeFileName(date, subject_id, 'rolloverShapes'));
    rollover_shapes = struct;
    rollover_shapes.left_fit_theta_x = left_fit_theta_x;
    rollover_shapes.left_fit_theta_y = left_fit_theta_y;
    rollover_shapes.right_fit_theta_x = right_fit_theta_x;
    rollover_shapes.right_fit_theta_y = right_fit_theta_y;
end

left_ankle_scs_to_world_rotation_reference = plant.endEffectorTransformations{3}(1:3, 1:3);
right_ankle_scs_to_world_rotation_reference = plant.endEffectorTransformations{6}(1:3, 1:3);

for i_trial = trials_to_process
    if ~ismember(i_trial, trials_to_exclude)
	
        display(['Calculating inverse dynamics: trial ' num2str(i_trial)]);
        load(makeFileName(date, subject_id, 'walking', i_trial, 'markerTrajectories'));
        load(makeFileName(date, subject_id, 'walking', i_trial, 'angleTrajectories'));
        load(makeFileName(date, subject_id, 'walking', i_trial, 'kinematicTrajectories'));
        load(makeFileName(date, subject_id, 'walking', i_trial, 'stepEvents'));
        load(makeFileName(date, subject_id, 'walking', i_trial, 'forceplateTrajectories'));
        number_of_time_steps = size(joint_angle_trajectories_belt, 1);
        if process_all_data
            data_points = 2001 : 3000;
            data_points = 1 : 30000;
            data_points = 3100 : 3200;
%             data_points = 3100 : 3100;
        else
            load(makeFileName(date, subject_id, 'walking', i_trial, 'relevantDataStretches'));
        end
        

        if load_results
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
            if use_rollover_constraints
                label = [label '_rolloverConstraints'];
            end            
            if use_full_constraints
                label = [label '_fullConstraints'];
            end            

            load(makeFileName(date, subject_id, 'walking', i_trial, label));
        end
        
        %% determine_constraint_numbers
        if determine_constraint_numbers
            % 0 = swing
            % 1 = heelstrike
            % 2 = pushoff
            gamma_left_trajectory = zeros(number_of_time_steps, 1);
            phi_left_trajectory = zeros(number_of_time_steps, 1);
            rho_left_trajectory = zeros(number_of_time_steps, 1);
            theta_left_trajectory = zeros(number_of_time_steps, 1);
            T_foot_to_world_left_trajectory = zeros(number_of_time_steps, 16);
            left_foot_constraint_number_trajectory = zeros(number_of_time_steps, 1);
            gamma_right_trajectory = zeros(number_of_time_steps, 1);
            phi_right_trajectory = zeros(number_of_time_steps, 1);
            rho_right_trajectory = zeros(number_of_time_steps, 1);
            theta_right_trajectory = zeros(number_of_time_steps, 1);
            T_foot_to_world_right_trajectory = zeros(number_of_time_steps, 16);
            right_foot_constraint_number_trajectory = zeros(number_of_time_steps, 1);
            
            % extract data
            left_ankle_marker = find(strcmp(marker_headers, 'LANK'));
            left_ankle_marker_indices = reshape([(left_ankle_marker - 1) * 3 + 1; (left_ankle_marker - 1) * 3 + 2; (left_ankle_marker - 1) * 3 + 3], 1, length(left_ankle_marker)*3);
            left_ankle_marker_pos_trajectory = marker_trajectories(:, left_ankle_marker_indices);

            left_toes_marker = find(strcmp(marker_headers, 'LTOE'));
            left_toes_marker_indices = reshape([(left_toes_marker - 1) * 3 + 1; (left_toes_marker - 1) * 3 + 2; (left_toes_marker - 1) * 3 + 3], 1, length(left_toes_marker)*3);
            left_toes_marker_pos_trajectory = marker_trajectories(:, left_toes_marker_indices);

            right_ankle_marker = find(strcmp(marker_headers, 'RANK'));
            right_ankle_marker_indices = reshape([(right_ankle_marker - 1) * 3 + 1; (right_ankle_marker - 1) * 3 + 2; (right_ankle_marker - 1) * 3 + 3], 1, length(right_ankle_marker)*3);
            right_ankle_marker_pos_trajectory = marker_trajectories(:, right_ankle_marker_indices);

            right_toes_marker = find(strcmp(marker_headers, 'RTOE'));
            right_toes_marker_indices = reshape([(right_toes_marker - 1) * 3 + 1; (right_toes_marker - 1) * 3 + 2; (right_toes_marker - 1) * 3 + 3], 1, length(right_toes_marker)*3);
            right_toes_marker_pos_trajectory = marker_trajectories(:, right_toes_marker_indices);
            
            for i_time = data_points
                % right foot segment euler angles
                left_ankle_scs_transformation_current = reshape(T_left_ankle_to_world_trajectory(i_time, :), 4, 4);
                left_ankle_scs_to_world_rotation_current = left_ankle_scs_transformation_current(1:3, 1:3);
                left_ankle_scs_rotation_reference_to_current = left_ankle_scs_to_world_rotation_reference^(-1) * left_ankle_scs_to_world_rotation_current;
                left_euler_angles = eulerAnglesFromRotationMatrixYZX(left_ankle_scs_rotation_reference_to_current);
                gamma_left_trajectory(i_time) = left_euler_angles(1);
                phi_left_trajectory(i_time) = left_euler_angles(2);
                rho_left_trajectory(i_time) = left_euler_angles(3);

                % right foot segment euler angles
                right_ankle_scs_transformation_current = reshape(T_right_ankle_to_world_trajectory(i_time, :), 4, 4);
                right_ankle_scs_to_world_rotation_current = right_ankle_scs_transformation_current(1:3, 1:3);
                right_ankle_scs_rotation_reference_to_current = right_ankle_scs_to_world_rotation_reference^(-1) * right_ankle_scs_to_world_rotation_current;
                right_euler_angles = eulerAnglesFromRotationMatrixYZX(right_ankle_scs_rotation_reference_to_current);
                gamma_right_trajectory(i_time) = right_euler_angles(1);
                phi_right_trajectory(i_time) = right_euler_angles(2);
                rho_right_trajectory(i_time) = right_euler_angles(3);
                
                % left foot theta angle and transformation
                left_ankle_point_position_world = left_ankle_marker_pos_trajectory(i_time, :)';
                left_toes_position_world = left_toes_marker_pos_trajectory(i_time, :)';
                left_ankle_to_toes_vector = left_toes_position_world - left_ankle_point_position_world;
                foot_y = normVector(left_ankle_to_toes_vector);
                foot_z = normVector([0; 0; 1] - dot(foot_y, [0; 0; 1]) * foot_y);
                foot_x = cross(foot_y, foot_z);
                R_foot_to_world_left = [foot_x foot_y foot_z];
                T_foot_to_world_left = [R_foot_to_world_left left_ankle_point_position_world; 0 0 0 1];
                T_foot_to_world_left_trajectory(i_time, :) = reshape(T_foot_to_world_left, 1, 16);
                
                theta_left_trajectory(i_time, :) = -atan2(left_ankle_to_toes_vector(3), left_ankle_to_toes_vector(2));
                
                % right foot theta angle and transformation
                right_ankle_point_position_world = right_ankle_marker_pos_trajectory(i_time, :)';
                right_toes_position_world = right_toes_marker_pos_trajectory(i_time, :)';
                right_ankle_to_toes_vector = right_toes_position_world - right_ankle_point_position_world;
                foot_y = normVector(right_ankle_to_toes_vector);
                foot_z = normVector([0; 0; 1] - dot(foot_y, [0; 0; 1]) * foot_y);
                foot_x = cross(foot_y, foot_z);
                R_foot_to_world_right = [foot_x foot_y foot_z];
                T_foot_to_world_right = [R_foot_to_world_right right_ankle_point_position_world; 0 0 0 1];
                theta_right_trajectory(i_time, :) = -atan2(right_ankle_to_toes_vector(3), right_ankle_to_toes_vector(2));
                T_foot_to_world_right_trajectory(i_time, :) = reshape(T_foot_to_world_right, 1, 16);
                
                

                % left foot
                if left_contact_indicators_mocap(i_time)
                    if phi_left_trajectory(i_time) < 0
                        left_foot_constraint_number_trajectory(i_time) = 1;
                    elseif phi_left_trajectory(i_time) >= 0
                        left_foot_constraint_number_trajectory(i_time) = 2;
                    else
                        error('fun times, some number is neither smaller nor larger than 0.');
                    end
                else
                    left_foot_constraint_number_trajectory(i_time) = 0;
                end
                % right foot
                if right_contact_indicators_mocap(i_time)
                    if phi_right_trajectory(i_time) > 0
                        right_foot_constraint_number_trajectory(i_time) = 1;
                    elseif phi_right_trajectory(i_time) <= 0
                        right_foot_constraint_number_trajectory(i_time) = 2;
                    else
                        error('fun times, some number is neither smaller nor larger than 0.');
                    end
                else
                    right_foot_constraint_number_trajectory(i_time) = 0;
                end
            end
        end
        
        %% calculate_dynamic_matrices
        if calculate_dynamic_matrices
            tic
            fprintf('Calculating the dynamic matrices ... ');

            number_of_time_steps = size(joint_angle_trajectories_belt, 1);
            inertia_matrix_trajectory = cell(number_of_time_steps, 1);
            coriolis_matrix_trajectory = cell(number_of_time_steps, 1);
            gravitation_matrix_trajectory = cell(number_of_time_steps, 1);

            constraint_matrix_trajectory = cell(number_of_time_steps, 1);
            constraint_matrix_dot_trajectory = cell(number_of_time_steps, 1);
            number_of_lambdas = zeros(number_of_time_steps, 1);
            number_of_left_foot_constraints = zeros(number_of_time_steps, 1);
            number_of_right_foot_constraints = zeros(number_of_time_steps, 1);
                
                


            if use_parallel
                inertia_matrix_trajectory_pool = cell(size(inertia_matrix_trajectory));
                coriolis_matrix_trajectory_pool = cell(size(coriolis_matrix_trajectory));
                gravitation_matrix_trajectory_pool = cell(size(gravitation_matrix_trajectory));
                constraint_matrix_trajectory_pool = cell(size(constraint_matrix_trajectory));
                constraint_matrix_dot_trajectory_pool = cell(size(constraint_matrix_dot_trajectory));
                number_of_lambdas_pool = zeros(number_of_time_steps, 1);
                spmd
                    plant_pool = plant.copy;
                    for i_time = data_points(1)+labindex-1 : numlabs : data_points(end)
                        if any(isnan(joint_angle_trajectories_belt(i_time, :)))
                            inertia_matrix_trajectory_pool{i_time} = NaN;
                            coriolis_matrix_trajectory_pool{i_time} = NaN;
                            gravitation_matrix_trajectory_pool{i_time} = NaN;
                            number_of_lambdas_pool(i_time) = NaN;
                            constraint_matrix_trajectory_pool{i_time} = NaN;
                            constraint_matrix_dot_trajectory_pool{i_time} = NaN;

                        else
                            % update model
                            plant_pool.jointAngles = joint_angle_trajectories_belt(i_time, :)';
                            plant_pool.jointVelocities = joint_velocity_trajectories_belt(i_time, :)';
                            plant_pool.jointAccelerations = joint_acceleration_trajectories_belt(i_time, :)';
                            plant_pool.updateInternals;

                            % save dynamic matrices
                            inertia_matrix_trajectory_pool{i_time} = plant_pool.inertiaMatrix;
                            coriolis_matrix_trajectory_pool{i_time} = plant_pool.coriolisMatrix;
                            gravitation_matrix_trajectory_pool{i_time} = plant_pool.gravitationalTorqueMatrix;

                            if use_point_constraints
                                [constraint_matrix_trajectory_pool{i_time}, constraint_matrix_dot_trajectory_pool{i_time}] = ...
                                    createConstraintMatrix_pointConstraints ...
                                      ( ...
                                        plant_pool, ...
                                        left_foot_constraint_number_trajectory(i_time), ...
                                        right_foot_constraint_number_trajectory(i_time) ...
                                      );
                            end
                            if use_hinge_constraints
                                [ ...
                                  constraint_matrix_trajectory_pool{i_time}, ...
                                  constraint_matrix_dot_trajectory_pool{i_time}, ...
                                  number_of_left_foot_constraints(i_time), ...
                                  number_of_right_foot_constraints(i_time) ...
                                ] = ...
                                createConstraintMatrix_hingeConstraints ...
                                  ( ...
                                    plant_pool, ...
                                    left_foot_constraint_number_trajectory(i_time), ...
                                    right_foot_constraint_number_trajectory(i_time) ...
                                  );
                            end
                            if use_body_velocity_constraints
                                [constraint_matrix_trajectory_pool{i_time}, constraint_matrix_dot_trajectory_pool{i_time}] = ...
                                    createConstraintMatrix_bodyVelocityConstraints ...
                                      ( ...
                                        plant_pool, ...
                                        left_foot_constraint_number_trajectory(i_time), ...
                                        right_foot_constraint_number_trajectory(i_time), ...
                                        phi_left_trajectory(i_time), ...
                                        rho_left_trajectory(i_time), ...
                                        phi_right_trajectory(i_time), ...
                                        rho_right_trajectory(i_time), ...
                                        V_body_left_fits_heelstrike, ...
                                        V_body_left_fits_pushoff, ...
                                        V_body_right_fits_heelstrike, ...
                                        V_body_right_fits_pushoff...
                                      );
                            end
                           number_of_lambdas_pool(i_time) = size(constraint_matrix_trajectory_pool{i_time}, 1);
                        end
                    end
                end
                % reassemble
                for i_lab = 1 : number_of_labs
                    inertia_matrix_trajectory_lab = inertia_matrix_trajectory_pool{i_lab};
                    inertia_matrix_trajectory(data_points(1)+i_lab-1 : number_of_labs : data_points(end)) = inertia_matrix_trajectory_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end));
                    coriolis_matrix_trajectory_lab = coriolis_matrix_trajectory_pool{i_lab};
                    coriolis_matrix_trajectory(data_points(1)+i_lab-1 : number_of_labs : data_points(end)) = coriolis_matrix_trajectory_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end));
                    gravitation_matrix_trajectory_lab = gravitation_matrix_trajectory_pool{i_lab};
                    gravitation_matrix_trajectory(data_points(1)+i_lab-1 : number_of_labs : data_points(end)) = gravitation_matrix_trajectory_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end));

                    constraint_matrix_trajectory_lab = constraint_matrix_trajectory_pool{i_lab};
                    constraint_matrix_trajectory(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) = constraint_matrix_trajectory_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
                    constraint_matrix_dot_trajectory_lab = constraint_matrix_dot_trajectory_pool{i_lab};
                    constraint_matrix_dot_trajectory(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) = constraint_matrix_dot_trajectory_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);

                    number_of_lambdas_lab = number_of_lambdas_pool{i_lab};
                    number_of_lambdas(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) = number_of_lambdas_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
                end
            else
                for i_time = data_points
                    if any(isnan(joint_angle_trajectories_belt(i_time, :)))
                        inertia_matrix_trajectory{i_time} = NaN;
                        coriolis_matrix_trajectory{i_time} = NaN;
                        gravitation_matrix_trajectory{i_time} = NaN;
                        constraint_matrix_trajectory{i_time} = NaN;
                        constraint_matrix_dot_trajectory{i_time} = NaN;
                    else
                        % update model
                        plant.jointAngles = joint_angle_trajectories_belt(i_time, :)';
                        plant.jointVelocities = joint_velocity_trajectories_belt(i_time, :)';
                        plant.jointAccelerations = joint_acceleration_trajectories_belt(i_time, :)';
                        plant.updateInternals;
                        
                        % save dynamic matrices
                        inertia_matrix_trajectory{i_time} = plant.inertiaMatrix;
                        coriolis_matrix_trajectory{i_time} = plant.coriolisMatrix;
                        gravitation_matrix_trajectory{i_time} = plant.gravitationalTorqueMatrix;

                        if use_point_constraints
                            [constraint_matrix_trajectory{i_time}, constraint_matrix_dot_trajectory{i_time}] = ...
                                createConstraintMatrix_pointConstraints ...
                                  ( ...
                                    plant, ...
                                    left_foot_constraint_number_trajectory(i_time), ...
                                    right_foot_constraint_number_trajectory(i_time) ...
                                  );
                        end
                        if use_hinge_constraints
                            [ ...
                              constraint_matrix_trajectory{i_time}, ...
                              constraint_matrix_dot_trajectory{i_time}, ...
                              number_of_left_foot_constraints(i_time), ...
                              number_of_right_foot_constraints(i_time) ...
                            ] = ...
                            createConstraintMatrix_hingeConstraints ...
                              ( ...
                                plant, ...
                                left_foot_constraint_number_trajectory(i_time), ...
                                right_foot_constraint_number_trajectory(i_time) ...
                              );
                        end
                        if use_body_velocity_constraints
                            [constraint_matrix_trajectory{i_time}, constraint_matrix_dot_trajectory{i_time}] = ...
                                createConstraintMatrix_bodyVelocityConstraints ...
                                  ( ...
                                    plant, ...
                                    left_foot_constraint_number_trajectory(i_time), ...
                                    right_foot_constraint_number_trajectory(i_time), ...
                                    phi_left_trajectory(i_time), ...
                                    rho_left_trajectory(i_time), ...
                                    phi_right_trajectory(i_time), ...
                                    rho_right_trajectory(i_time), ...
                                    V_body_left_fits_heelstrike, ...
                                    V_body_left_fits_pushoff, ...
                                    V_body_right_fits_heelstrike, ...
                                    V_body_right_fits_pushoff ...
                                  );
                        end
                        if use_rollover_constraints
                            [constraint_matrix_trajectory{i_time}, constraint_matrix_dot_trajectory{i_time}] = ...
                                createConstraintMatrix_rolloverConstraints ...
                                  ( ...
                                    plant, ...
                                    rollover_shapes, ...
                                    left_foot_constraint_number_trajectory(i_time), ...
                                    right_foot_constraint_number_trajectory(i_time), ...
                                    theta_left_trajectory(i_time), ...
                                    theta_right_trajectory(i_time), ...
                                    reshape(T_foot_to_world_left_trajectory(i_time, :), 4, 4), ...
                                    reshape(T_foot_to_world_right_trajectory(i_time, :), 4, 4) ...
                                  );
                            
                        end
                        number_of_lambdas(i_time) = size(constraint_matrix_trajectory{i_time}, 1);

                    end
                    % give progress feedback
                    display_step = 1;
                    last_time_step = data_points(end);
                    if (i_time / display_step) == floor(i_time / display_step)
                        disp([num2str(i_time) '(' num2str(last_time_step) ')']);
                    end
                end
            end
            fprintf(' done\n');
            toc
        end

        %% calculate_torques
        if calculate_torques
            tic
            fprintf('Calculating the constraint torques ... ');
            number_of_joints = plant.numberOfJoints;
            number_of_time_steps = size(joint_angle_trajectories_belt, 1);
            virtual_joints = 1:6;
            constraint_torque_trajectories_all = zeros(number_of_time_steps, number_of_joints);
            constraint_torque_trajectories_right = zeros(number_of_time_steps, number_of_joints);
            constraint_torque_trajectories_left = zeros(number_of_time_steps, number_of_joints);
            lambda_trajectories = cell(number_of_time_steps, 1);
            joint_torque_trajectories = zeros(number_of_time_steps, number_of_joints);
            induced_accelerations_applied_trajectories = zeros(number_of_time_steps, number_of_joints);
            induced_accelerations_gravity_trajectories = zeros(number_of_time_steps, number_of_joints);
            induced_accelerations_movement_trajectories = zeros(number_of_time_steps, number_of_joints);
            induced_accelerations_applied_single_trajectories = zeros(number_of_time_steps, number_of_joints, number_of_joints); % (i_time, i, j)-th entry holds acceleration of the i_th joint from torque at the j_th joint
            violating_velocity_trajectories = zeros(number_of_time_steps, number_of_joints);
            violating_acceleration_trajectories = zeros(number_of_time_steps, number_of_joints);
%             explained_acceleration_trajectories = zeros(number_of_time_steps, number_of_joints);
%             leftover_acceleration_trajectories = zeros(number_of_time_steps, number_of_joints);

            if use_parallel
                constraint_torque_trajectories_all_pool = zeros(number_of_time_steps, number_of_joints);
                constraint_torque_trajectories_right_pool = zeros(number_of_time_steps, number_of_joints);
                constraint_torque_trajectories_left_pool = zeros(number_of_time_steps, number_of_joints);
                lambda_trajectories_pool = cell(number_of_time_steps, 1);
                joint_torque_trajectories_pool = zeros(number_of_time_steps, number_of_joints);
                induced_accelerations_applied_trajectories_pool = zeros(number_of_time_steps, number_of_joints);
                induced_accelerations_gravity_trajectories_pool = zeros(number_of_time_steps, number_of_joints);
                induced_accelerations_movement_trajectories_pool = zeros(number_of_time_steps, number_of_joints);
                induced_accelerations_applied_single_trajectories_pool = zeros(number_of_time_steps, number_of_joints, number_of_joints); 
                violating_velocity_trajectories_pool = zeros(number_of_time_steps, number_of_joints);
                violating_acceleration_trajectories_pool = zeros(number_of_time_steps, number_of_joints);
%                 explained_acceleration_trajectories_pool = zeros(number_of_time_steps, number_of_joints);
%                 leftover_acceleration_trajectories_pool = zeros(number_of_time_steps, number_of_joints);
                spmd
                    for i_time = data_points(1)+labindex-1 : numlabs : data_points(end)
                        if any(isnan(joint_angle_trajectories_belt(i_time, :)))
                            constraint_torque_trajectories_all_pool(i_time, :) = NaN;
                            constraint_torque_trajectories_right_pool(i_time, :) = NaN;
                            constraint_torque_trajectories_left_pool(i_time, :) = NaN;
                            lambda_trajectories_pool{i_time} = NaN;
                            joint_torque_trajectories_pool(i_time, :) = NaN;
                            induced_accelerations_applied_trajectories_pool(i_time, :) = NaN;
                            induced_accelerations_gravity_trajectories_pool(i_time, :) = NaN;
                            induced_accelerations_movement_trajectories_pool(i_time, :) = NaN;
                            induced_accelerations_applied_single_trajectories_pool(i_time, :, :) = NaN;
                            violating_velocity_trajectories_pool(i_time, :, :) = NaN;
                            violating_acceleration_trajectories_pool(i_time, :, :) = NaN;
%                             explained_acceleration_trajectories_pool(i_time, :, :) = NaN;
%                             leftover_acceleration_trajectories_pool(i_time, :, :) = NaN;
                        else
                            M = inertia_matrix_trajectory{i_time};
                            C = coriolis_matrix_trajectory{i_time};
                            N = gravitation_matrix_trajectory{i_time};
                            theta_dot = joint_velocity_trajectories_belt(i_time, :)';
                            theta_two_dot = joint_acceleration_trajectories_belt(i_time, :)';

                            % calculate constraint forces
%                             active_constraint = find(constraint_indicator_trajectories(i_time, :));

                            A = constraint_matrix_trajectory{i_time};
                            A_dot = constraint_matrix_dot_trajectory{i_time};

                            k_c = rank(A);
                            k_v = length(virtual_joints);
                            [~, ~, V_c] = svd(A);
                            C_c = V_c(:, k_c+1:end);
                            B_v = [eye(k_v); zeros(number_of_joints - k_v, k_v)];
                            k_w = rank([B_v C_c]);
                            [~, ~, V_w] = svd([B_v C_c]');
                            C_w = V_w(:, k_w+1:end);
                            D = [B_v C_w];
                            if rank(A) > 0
                                P = eye(plant.numberOfJoints) - A' * (A * M^(-1) * A')^(-1) * A * M^(-1);
                                Q = M*theta_two_dot + P*C*theta_dot + P*N - A'*(A*M^(-1)*A')^(-1)*(A * theta_two_dot);
                            else
                                P = eye(plant.numberOfJoints);
                                Q = M*theta_two_dot + C*theta_dot + N;
                            end
                            E = [P; D'];

                            H = A'*(A*M^(-1)*A')^(-1)*A;
                            T = pinv(E) * [Q; zeros(size(D, 2), 1)];
                            if rank(A) == 0
                                lambda = 0;
                            else
                                lambda = (A*M^(-1)*A')^(-1) *  A*M^(-1)*(T - C*theta_dot - N) - (A*M^(-1)*A')^(-1)*A*theta_two_dot;
                                test_left_side = M*theta_two_dot + C*theta_dot + N + A'*lambda;
                                lambda_T = (A*M^(-1)*A')^(-1) *  A*M^(-1) * T;
                                lambda_C = - (A*M^(-1)*A')^(-1) *  A*M^(-1) * C*theta_dot;
                                lambda_N = - (A*M^(-1)*A')^(-1) *  A*M^(-1) * N;
                            end


                            % calculate induced accelerations
                            if rank(A) == 0
                                induced_acceleration_applied = M^(-1)*T;
                                induced_acceleration_gravity = - M^(-1)*N;
                                induced_acceleration_movement = - M^(-1)*C*theta_dot;
                                induced_accelerations_applied_single = zeros(number_of_joints);
                                for i_joint = 1 : number_of_joints
                                    T_tilde_i = T;
                                    T_tilde_i(i_joint) = 0;
                                    T_bar_i = T - T_tilde_i;
                                    accelerations_induced_by_T_bar_i = M^(-1)*T_bar_i;
                                    induced_accelerations_applied_single(:, i_joint) = accelerations_induced_by_T_bar_i;
                                end
                                induced_accelerations_applied_single_check = sum(induced_accelerations_applied_single, 2);
                            else
                                induced_acceleration_applied = M^(-1)*(T - A'*lambda_T);
                                induced_acceleration_gravity = - M^(-1)*(N + A'*lambda_N);
                                induced_acceleration_movement = - M^(-1)*(C*theta_dot + A'*lambda_C);
                                induced_accelerations_applied_single = zeros(number_of_joints);
                                % calculate lambdas corresponding to each joint torque
                                for i_joint = 1 : number_of_joints
                                    T_tilde_i = T;
                                    T_tilde_i(i_joint) = 0;
                                    T_bar_i = T - T_tilde_i;
                                    lambda_T_bar_i = (A*M^(-1)*A')^(-1) *  A*M^(-1) * T_bar_i;
                                    accelerations_induced_by_T_bar_i = M^(-1)*(T_bar_i - A'*lambda_T_bar_i);
                                    induced_accelerations_applied_single(:, i_joint) = accelerations_induced_by_T_bar_i;
                                end
                                induced_accelerations_applied_single_check = sum(induced_accelerations_applied_single, 2);

                            end
                            
                            % calculate in how far the velocities are violated by the constraints
                            B_c = V_c(:, 1:k_c); % B_c spans the orthogonal complement of the null space of A (those velocities that are ONLY violating the constraints)
                            P_b = B_c*B_c';
                            violating_velocity_trajectories_pool(i_time, :, :) = P_b * theta_dot;
                            violating_acceleration_trajectories_pool(i_time, :) = P_b * theta_two_dot;
                            explained_acceleration = M^(-1)*(T - C*theta_dot - N - A'*lambda);
%                             explained_acceleration_trajectories_pool(i_time, :) = explained_acceleration;
%                             leftover_acceleration_trajectories_pool(i_time, :) = theta_two_dot - explained_acceleration;

                            if use_body_velocity_constraints
                                if left_foot_constraint_number_trajectory(i_time) == 0
                                    number_of_left_foot_constraints = 0;
                                elseif left_foot_constraint_number_trajectory(i_time) == 1
                                    number_of_left_foot_constraints = 4;
                                elseif left_foot_constraint_number_trajectory(i_time) == 2
                                    number_of_left_foot_constraints = 4;
                                else
                                    error('Left foot constraint number must be an integer between 0 and 2')
                                end

                                if right_foot_constraint_number_trajectory(i_time) == 0
                                    number_of_right_foot_constraints = 0;
                                elseif right_foot_constraint_number_trajectory(i_time) == 1
                                    number_of_right_foot_constraints = 4;
                                elseif right_foot_constraint_number_trajectory(i_time) == 2
                                    number_of_right_foot_constraints = 4;
                                else
                                    error('Right foot constraint number must be an integer between 0 and 2')
                                end
                            end
                            if use_hinge_constraints
                                if left_foot_constraint_number_trajectory(i_time) == 0
                                    number_of_left_foot_constraints = 0;
                                elseif left_foot_constraint_number_trajectory(i_time) == 1
                                    number_of_left_foot_constraints = 5;
                                elseif left_foot_constraint_number_trajectory(i_time) == 2
                                    number_of_left_foot_constraints = 5;
                                else
                                    error('Left foot constraint number must be an integer between 0 and 2')
                                end

                                if right_foot_constraint_number_trajectory(i_time) == 0
                                    number_of_right_foot_constraints = 0;
                                elseif right_foot_constraint_number_trajectory(i_time) == 1
                                    number_of_right_foot_constraints = 5;
                                elseif right_foot_constraint_number_trajectory(i_time) == 2
                                    number_of_right_foot_constraints = 5;
                                else
                                    error('Right foot constraint number must be an integer between 0 and 2')
                                end
                                
                            end
                            if use_point_constraints
                                if left_foot_constraint_number_trajectory(i_time) == 0
                                    number_of_left_foot_constraints = 0;
                                elseif left_foot_constraint_number_trajectory(i_time) == 1
                                    number_of_left_foot_constraints = 3;
                                elseif left_foot_constraint_number_trajectory(i_time) == 2
                                    number_of_left_foot_constraints = 3;
                                else
                                    error('Left foot constraint number must be an integer between 0 and 2')
                                end

                                if right_foot_constraint_number_trajectory(i_time) == 0
                                    number_of_right_foot_constraints = 0;
                                elseif right_foot_constraint_number_trajectory(i_time) == 1
                                    number_of_right_foot_constraints = 3;
                                elseif right_foot_constraint_number_trajectory(i_time) == 2
                                    number_of_right_foot_constraints = 3;
                                else
                                    error('Right foot constraint number must be an integer between 0 and 2')
                                end
                            end
                            if number_of_right_foot_constraints + number_of_left_foot_constraints == 0
                                lambda_right = 0;
                                lambda_left = 0;
                            else
                                lambda_right = [lambda(1:number_of_right_foot_constraints); zeros(number_of_left_foot_constraints, 1)];
                                lambda_left = [zeros(number_of_right_foot_constraints, 1); lambda(number_of_right_foot_constraints+1 : number_of_right_foot_constraints+number_of_left_foot_constraints);];
                            end



                            % save to lists
                            constraint_torque_trajectories_all_pool(i_time, :) = A' * lambda;
                            constraint_torque_trajectories_right_pool(i_time, :) = A' * lambda_right;
                            constraint_torque_trajectories_left_pool(i_time, :) = A' * lambda_left;
                            lambda_trajectories_pool{i_time} = lambda;
                            joint_torque_trajectories_pool(i_time, :) = T;
                            induced_accelerations_applied_trajectories_pool(i_time, :) = induced_acceleration_applied;
                            induced_accelerations_gravity_trajectories_pool(i_time, :) = induced_acceleration_gravity;
                            induced_accelerations_movement_trajectories_pool(i_time, :) = induced_acceleration_movement;
                            induced_accelerations_applied_single_trajectories_pool(i_time, :, :) = induced_accelerations_applied_single;

                        end
                    end
                end
                % reassemble
                for i_lab = 1 : number_of_labs
                    constraint_torque_trajectories_all_lab = constraint_torque_trajectories_all_pool{i_lab};
                    constraint_torque_trajectories_right_lab = constraint_torque_trajectories_right_pool{i_lab};
                    constraint_torque_trajectories_left_lab = constraint_torque_trajectories_left_pool{i_lab};
                    lambda_trajectories_lab = lambda_trajectories_pool{i_lab};
                    joint_torque_trajectories_lab = joint_torque_trajectories_pool{i_lab};
%                     right_ground_reaction_wrench_trajectory_origin_lab = right_ground_reaction_wrench_trajectory_origin_pool{i_lab};
%                     left_ground_reaction_wrench_trajectory_origin_lab = left_ground_reaction_wrench_trajectory_origin_pool{i_lab};
                    induced_accelerations_applied_trajectories_lab = induced_accelerations_applied_trajectories_pool{i_lab};
                    induced_accelerations_gravity_trajectories_lab = induced_accelerations_gravity_trajectories_pool{i_lab};
                    induced_accelerations_movement_trajectories_lab = induced_accelerations_movement_trajectories_pool{i_lab};
                    induced_accelerations_applied_single_trajectories_lab = induced_accelerations_applied_single_trajectories_pool{i_lab};
                    violating_velocity_trajectories_lab = violating_velocity_trajectories_pool{i_lab};
                    violating_acceleration_trajectories_lab = violating_acceleration_trajectories_pool{i_lab};
%                     explained_acceleration_trajectories_lab = explained_acceleration_trajectories_pool{i_lab};
%                     leftover_acceleration_trajectories_lab = leftover_acceleration_trajectories_pool{i_lab};

                    constraint_torque_trajectories_all(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
                        = constraint_torque_trajectories_all_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
                    constraint_torque_trajectories_right(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
                        = constraint_torque_trajectories_right_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
                    constraint_torque_trajectories_left(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
                        = constraint_torque_trajectories_left_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
                    lambda_trajectories(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
                        = lambda_trajectories_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
                    joint_torque_trajectories(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
                        = joint_torque_trajectories_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
                    induced_accelerations_applied_trajectories(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
                        = induced_accelerations_applied_trajectories_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
                    induced_accelerations_gravity_trajectories(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
                        = induced_accelerations_gravity_trajectories_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
                    induced_accelerations_movement_trajectories(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
                        = induced_accelerations_movement_trajectories_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
                    induced_accelerations_applied_single_trajectories(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :, :) ...
                        = induced_accelerations_applied_single_trajectories_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :, :);
                    violating_velocity_trajectories(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
                        = violating_velocity_trajectories_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
                    violating_acceleration_trajectories(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
                        = violating_acceleration_trajectories_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
%                     explained_acceleration_trajectories(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
%                         = explained_acceleration_trajectories_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
%                     leftover_acceleration_trajectories(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
%                         = leftover_acceleration_trajectories_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
                end
            else
                for i_time = data_points
                    if any(isnan(joint_angle_trajectories_belt(i_time, :)))
                        constraint_torque_trajectories_all(i_time, :) = NaN;
                        constraint_torque_trajectories_right(i_time, :) = NaN;
                        constraint_torque_trajectories_left(i_time, :) = NaN;
                        lambda_trajectories{i_time} = NaN;
                        joint_torque_trajectories(i_time, :) = NaN;
                        right_ground_reaction_wrench_trajectory_origin(i_time, :) = NaN;
                        left_ground_reaction_wrench_trajectory_origin(i_time, :) = NaN;
                        induced_accelerations_applied_trajectories(i_time, :) = NaN;
                        induced_accelerations_gravity_trajectories(i_time, :) = NaN;
                        induced_accelerations_movement_trajectories(i_time, :) = NaN;
                        induced_accelerations_applied_single_trajectories(i_time, :, :) = NaN;
                        violating_velocity_trajectories(i_time, :, :) = NaN;
                        violating_acceleration_trajectories(i_time, :, :) = NaN;
    %                     explained_acceleration_trajectories(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
    %                         = explained_acceleration_trajectories_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
    %                     leftover_acceleration_trajectories(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :) ...
    %                         = leftover_acceleration_trajectories_lab(data_points(1)+i_lab-1 : number_of_labs : data_points(end), :);
                    else
                        % extract relevant variables
                        M = inertia_matrix_trajectory{i_time};
                        C = coriolis_matrix_trajectory{i_time};
                        N = gravitation_matrix_trajectory{i_time};
                        theta_dot = joint_velocity_trajectories_belt(i_time, :)';
                        theta_two_dot = joint_acceleration_trajectories_belt(i_time, :)';

                        % calculate constraint forces
%                         active_constraint = find(constraint_indicator_trajectories(i_time, :));


                        % get constraint matrices
                        A = constraint_matrix_trajectory{i_time};
                        A_dot = constraint_matrix_dot_trajectory{i_time};
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         % begin old stuff
%                         k_c = rank(A);
%                         k_v = length(virtual_joints);
%                         [~, ~, V_c] = svd(A);
%                         C_c = V_c(:, k_c+1:end); % C_c spans the null space of A
%                         B_v = [eye(k_v); zeros(number_of_joints - k_v, k_v)];
%                         k_w = rank([B_v C_c]);
%                         [~, ~, V_w] = svd([B_v C_c]');
%                         C_w = V_w(:, k_w+1:end);
%                         D = [B_v C_w];
%                         if rank(A) > 0
%                             P = eye(plant.numberOfJoints) - A' * (A * M^(-1) * A')^(-1) * A * M^(-1);
%                             Q = M*theta_two_dot + P*C*theta_dot + P*N - A'*(A*M^(-1)*A')^(-1)*(A * theta_two_dot);
%                             Q = M*theta_two_dot + P*C*theta_dot + P*N + A'*(A*M^(-1)*A')^(-1)*(A_dot * theta_dot);
%                         else
%                             P = eye(plant.numberOfJoints);
%                             Q = M*theta_two_dot + C*theta_dot + N;
%                         end
%                         E = [P; D'];
% 
%                         H = A'*(A*M^(-1)*A')^(-1)*A;
%                         T = pinv(E) * [Q; zeros(size(D, 2), 1)];
%                         if rank(A) == 0
%                             lambda = 0;
%                         else
%                             lambda = (A*M^(-1)*A')^(-1) *  A*M^(-1)*(T - C*theta_dot - N) - (A*M^(-1)*A')^(-1)*A*theta_two_dot;
%                             test_left_side = M*theta_two_dot + C*theta_dot + N + A'*lambda;
%                             lambda_T = (A*M^(-1)*A')^(-1) *  A*M^(-1) * T;
%                             lambda_C = - (A*M^(-1)*A')^(-1) *  A*M^(-1) * C*theta_dot;
%                             lambda_N = - (A*M^(-1)*A')^(-1) *  A*M^(-1) * N;
% 
%     %                                 lambda_leftover = -(A*M^(-1)*A')^(-1) * A * theta_two_dot;
%     %                                 test_right_side = lambda_T + lambda_C + lambda_N + lambda_leftover;
% 
%     %                                 T_pure = T;
%     %                                 T = T_pure - A'*lambda_leftover;
%                         end                        
%                         T_old = T;
%                         
%                         % end old stuff
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        
                        
                        
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % begin new stuff
                        P = eye(plant.numberOfJoints) - A' * (A * M^(-1) * A')^(-1) * A * M^(-1);

                        % force plate data
                        [~, time_index_forceplate] = min(abs(time_forceplate - time_mocap(i_time)));
                        combined_forceplate_wrench_data = total_forceplate_wrench_Acw(time_index_forceplate, :)';
                        if number_of_left_foot_constraints(i_time) ~= 0 & number_of_right_foot_constraints(i_time) == 0
                            Q_left = eye(number_of_left_foot_constraints(i_time));
                            R_left = pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 12)');
                            H_1_left = R_left * A' * Q_left * (A*M^(-1)*A')^(-1)*A*M^(-1);
                            b_1_left = combined_forceplate_wrench_data - R_left * A' * Q_left * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*plant.jointVelocities - N) + A_dot*plant.jointVelocities);
                            
                            H_1 = H_1_left;
                            b_1 = b_1_left;
                        elseif number_of_left_foot_constraints(i_time) == 0 & number_of_right_foot_constraints(i_time) ~= 0
                            Q_right = eye(number_of_right_foot_constraints(i_time));
                            R_right = pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 18)');
                            H_1_right = R_right * A' * Q_right * (A*M^(-1)*A')^(-1)*A*M^(-1);
                            b_1_right = combined_forceplate_wrench_data - R_right * A' * Q_right * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*plant.jointVelocities - N) + A_dot*plant.jointVelocities);
                            
                            H_1 = H_1_right;
                            b_1 = b_1_right;
                        elseif number_of_left_foot_constraints(i_time) ~= 0 & number_of_right_foot_constraints(i_time) ~= 0
                            Q_left = [eye(number_of_left_foot_constraints(i_time)) zeros(number_of_right_foot_constraints(i_time), number_of_left_foot_constraints(i_time)); zeros(number_of_left_foot_constraints(i_time), number_of_right_foot_constraints(i_time)) zeros(number_of_right_foot_constraints(i_time))];
                            R_left = pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 12)');
                            H_1_left = R_left * A' * Q_left * (A*M^(-1)*A')^(-1)*A*M^(-1);

                            Q_right = [zeros(number_of_left_foot_constraints(i_time)) zeros(number_of_right_foot_constraints(i_time), number_of_left_foot_constraints(i_time)); zeros(number_of_left_foot_constraints(i_time), number_of_right_foot_constraints(i_time)) eye(number_of_right_foot_constraints(i_time))];
                            R_right = pinv(plant.calculateArbitraryFrameBodyJacobian(eye(4, 4), 18)');
                            H_1_right = R_right * A' * Q_right * (A*M^(-1)*A')^(-1)*A*M^(-1);
                            
                            b_1 = combined_forceplate_wrench_data ...
                                  - R_left * A' * Q_left * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*plant.jointVelocities - N) + A_dot*plant.jointVelocities) ...
                                  - R_right * A' * Q_right * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*plant.jointVelocities - N) + A_dot*plant.jointVelocities);
                        end
                        
                        
                        % acceleration data
                        

                        % this gets the same results as the old version
                        H_2 = P;
                        b_2 = (M*theta_two_dot + P*C*theta_dot + P*N + A'*(A*M^(-1)*A')^(-1)*(A_dot * theta_dot));

                        % according to equation - new version
                        H_2 = M^(-1)*P;
                        b_2 = theta_two_dot + M^(-1) * P*(C*theta_dot + N) + M^(-1) * A'*(A*M^(-1)*A')^(-1)*A_dot*theta_dot;
                        
                        
                        
                        
                        
                        
                        % no virtual torques
                        k_v = length(virtual_joints);
                        B = [eye(k_v) zeros(k_v, number_of_joints - k_v)];
                        
                        % no workless torques
                        k_c = rank(A);
                        [~, ~, V_c] = svd(A);
                        C_c = V_c(:, k_c+1:end); % C_c spans the null space of A
                        k_w = rank([B' C_c]);
                        [~, ~, V_w] = svd([B' C_c]');
                        C_w = V_w(:, k_w+1:end);                        
                        
                        
                        % combine
                        H = ...
                          [ ...
                            H_2; ...
                            B; ...
                            C_w'; ...
                          ];
%                             H_1; ...
                        b = ...
                          [ ...
                            b_2; ...
                            zeros(k_v, 1); ...
                            zeros(number_of_joints - k_w, 1); ...
                          ];
%                             b_1; ...
                        
% status 31.8.2016: I realized that the old code probably had an error, where I would use Q = M*theta_two_dot + P*C*theta_dot + P*N - A'*(A*M^(-1)*A')^(-1)*(A * theta_two_dot);
% with A * theta_two_dot instead of A_dot * theta_dot, and wrong sign. The corrected version does not change the
% calculated ground reaction forces much. 

% Now, add the GRF constraint for a stretch with only one contact

                        T = pinv(H) * b;
                        
%                         T_check = [T_old T];
%                         T_check';
                        
%                         lambda = (A*M^(-1)*A')^(-1) *  A*M^(-1)*(T - C*theta_dot - N) - (A*M^(-1)*A')^(-1)*A*theta_two_dot; % this also has the silly error using A*theta_two_dot instead of A_dot*theta_dot
                        lambda = (A*M^(-1)*A')^(-1) *  A*M^(-1)*(T - C*theta_dot - N) + (A*M^(-1)*A')^(-1)*A_dot*theta_dot; % corrected here
                        
                        % this is equal now, move forward from here
%                         
%                         % checks
%                         check_1 = H_1 * T - b_1;
%                         check_2 = H_2 * T - b_2;
%                         check_3 = B * T;
%                         
%                         
%                         
% 
% %                         k_c = rank(A);
% %                         [~, ~, V_c] = svd(A);
% %                         C_c = V_c(:, k_c+1:end); % C_c spans the null space of A
% %                         k_w = rank([B_v C_c]);
% %                         [~, ~, V_w] = svd([B_v C_c]');
% %                         C_w = V_w(:, k_w+1:end);
% %                         D = [B_v C_w];
% %                         if rank(A) > 0
% %                             P = eye(plant.numberOfJoints) - A' * (A * M^(-1) * A')^(-1) * A * M^(-1);
% %                             Q = M*theta_two_dot + P*C*theta_dot + P*N - A'*(A*M^(-1)*A')^(-1)*(A * theta_two_dot);
% %                         else
% %                             P = eye(plant.numberOfJoints);
% %                             Q = M*theta_two_dot + C*theta_dot + N;
% %                         end
% %                         E = [P; D'];
% 
% %                         H = A'*(A*M^(-1)*A')^(-1)*A;
% %                         T = pinv(E) * [Q; zeros(size(D, 2), 1)];
% %                         if rank(A) == 0
% %                             lambda = 0;
% %                         else
%                             lambda = (A*M^(-1)*A')^(-1) *  A*M^(-1)*(T - C*theta_dot - N) - (A*M^(-1)*A')^(-1)*A*theta_two_dot;
%                             test_left_side = M*theta_two_dot + C*theta_dot + N + A'*lambda;
%                             lambda_T = (A*M^(-1)*A')^(-1) *  A*M^(-1) * T;
%                             lambda_C = - (A*M^(-1)*A')^(-1) *  A*M^(-1) * C*theta_dot;
%                             lambda_N = - (A*M^(-1)*A')^(-1) *  A*M^(-1) * N;

                        % end new stuff
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
    %                                 lambda_leftover = -(A*M^(-1)*A')^(-1) * A * theta_two_dot;
    %                                 test_right_side = lambda_T + lambda_C + lambda_N + lambda_leftover;

    %                                 T_pure = T;
    %                                 T = T_pure - A'*lambda_leftover;
%                         end

                        % calculate induced accelerations
%                         if rank(A) == 0
%                             induced_acceleration_applied = M^(-1)*T;
%                             induced_acceleration_gravity = - M^(-1)*N;
%                             induced_acceleration_movement = - M^(-1)*C*theta_dot;
%                             induced_accelerations_applied_single = zeros(number_of_joints);
%                             for i_joint = 1 : number_of_joints
%                                 T_tilde_i = T;
%                                 T_tilde_i(i_joint) = 0;
%                                 T_bar_i = T - T_tilde_i;
%                                 accelerations_induced_by_T_bar_i = M^(-1)*T_bar_i;
%                                 induced_accelerations_applied_single(:, i_joint) = accelerations_induced_by_T_bar_i;
%                             end
%                             induced_accelerations_applied_single_check = sum(induced_accelerations_applied_single, 2);
%                         else
%                             induced_acceleration_applied = M^(-1)*(T - A'*lambda_T);
%                             induced_acceleration_gravity = - M^(-1)*(N + A'*lambda_N);
%                             induced_acceleration_movement = - M^(-1)*(C*theta_dot + A'*lambda_C);
%                             induced_accelerations_applied_single = zeros(number_of_joints);
%                             % calculate lambdas corresponding to each joint torque
%                             lambda_T_bar_i_collection = zeros(size(lambda, 1), number_of_joints);
%                             T_bar_i_collection = zeros(number_of_joints, number_of_joints);
%                             for i_joint = 1 : number_of_joints
%                                 T_tilde_i = T;
%                                 T_tilde_i(i_joint) = 0;
%                                 T_bar_i = T - T_tilde_i;
%                                 T_bar_i_collection(:, i_joint) = T_bar_i;
%                                 lambda_T_bar_i = (A*M^(-1)*A')^(-1) *  A*M^(-1) * T_bar_i;
%                                 lambda_T_bar_i_collection(:, i_joint) = lambda_T_bar_i;
%                                 accelerations_induced_by_T_bar_i = M^(-1)*(T_bar_i - A'*lambda_T_bar_i); % maybe this sign is wrong? should this be a + instead?
%                                 induced_accelerations_applied_single(:, i_joint) = accelerations_induced_by_T_bar_i;
%                             end
%                             T_bar_i_check = sum(T_bar_i_collection, 2);
%                             lambda_T_bar_i_check = sum(lambda_T_bar_i_collection, 2);
%                             induced_accelerations_applied_single_check = sum(induced_accelerations_applied_single, 2);
% 
%                         end

                        % calculate in how far the velocities are violated by the constraints
%                         B_c = V_c(:, 1:k_c); % B_c spans the orthogonal complement of the null space of A (those velocities that are ONLY violating the constraints)
%                         P_b = B_c*B_c';
%                         violating_velocity_trajectories(i_time, :) = P_b * theta_dot;
%                         violating_acceleration_trajectories(i_time, :) = P_b * theta_two_dot;
%                         explained_acceleration = M^(-1)*(T - C*theta_dot - N - A'*lambda);
%                         explained_acceleration_trajectories(i_time, :) = explained_acceleration;
%                         leftover_acceleration_trajectories(i_time, :) = theta_two_dot - explained_acceleration;



                        % test
%                         P_c = C_c*C_c';
%                         q_rand = randn(number_of_joints, 1);
%                         q_rand_b = P_b * q_rand;
%                         q_rand_c = P_c * q_rand; % this should lie in the null space of A, i.e. A*q_rand_c = 0




%                         if use_body_velocity_constraints
%                             if left_foot_constraint_number_trajectory(i_time) == 0
%                                 number_of_left_foot_constraints = 0;
%                             elseif left_foot_constraint_number_trajectory(i_time) == 1
%                                 number_of_left_foot_constraints = 4;
%                             elseif left_foot_constraint_number_trajectory(i_time) == 2
%                                 number_of_left_foot_constraints = 4;
%                             else
%                                 error('Left foot constraint number must be an integer between 0 and 2')
%                             end
%                             
%                             if right_foot_constraint_number_trajectory(i_time) == 0
%                                 number_of_right_foot_constraints = 0;
%                             elseif right_foot_constraint_number_trajectory(i_time) == 1
%                                 number_of_right_foot_constraints = 4;
%                             elseif right_foot_constraint_number_trajectory(i_time) == 2
%                                 number_of_right_foot_constraints = 4;
%                             else
%                                 error('Right foot constraint number must be an integer between 0 and 2')
%                             end
%                         end
%                         if use_hinge_constraints
%                             if left_foot_constraint_number_trajectory(i_time) == 0
%                                 number_of_left_foot_constraints = 0;
%                             elseif left_foot_constraint_number_trajectory(i_time) == 1
%                                 number_of_left_foot_constraints = 5;
%                             elseif left_foot_constraint_number_trajectory(i_time) == 2
%                                 number_of_left_foot_constraints = 5;
%                             else
%                                 error('Left foot constraint number must be an integer between 0 and 2')
%                             end
%                             
%                             if right_foot_constraint_number_trajectory(i_time) == 0
%                                 number_of_right_foot_constraints = 0;
%                             elseif right_foot_constraint_number_trajectory(i_time) == 1
%                                 number_of_right_foot_constraints = 5;
%                             elseif right_foot_constraint_number_trajectory(i_time) == 2
%                                 number_of_right_foot_constraints = 5;
%                             else
%                                 error('Right foot constraint number must be an integer between 0 and 2')
%                             end
%                         end
%                         if use_point_constraints
%                             if left_foot_constraint_number_trajectory(i_time) == 0
%                                 number_of_left_foot_constraints = 0;
%                             elseif left_foot_constraint_number_trajectory(i_time) == 1
%                                 number_of_left_foot_constraints = 3;
%                             elseif left_foot_constraint_number_trajectory(i_time) == 2
%                                 number_of_left_foot_constraints = 3;
%                             else
%                                 error('Left foot constraint number must be an integer between 0 and 2')
%                             end
%                             
%                             if right_foot_constraint_number_trajectory(i_time) == 0
%                                 number_of_right_foot_constraints = 0;
%                             elseif right_foot_constraint_number_trajectory(i_time) == 1
%                                 number_of_right_foot_constraints = 3;
%                             elseif right_foot_constraint_number_trajectory(i_time) == 2
%                                 number_of_right_foot_constraints = 3;
%                             else
%                                 error('Right foot constraint number must be an integer between 0 and 2')
%                             end
%                         end
                        if number_of_left_foot_constraints(i_time) + number_of_right_foot_constraints(i_time) == 0
                            lambda_right = 0;
                            lambda_left = 0;
                        else
                            lambda_left = [lambda(1:number_of_left_foot_constraints(i_time)); zeros(number_of_right_foot_constraints(i_time), 1)];
                            lambda_right = [zeros(number_of_left_foot_constraints(i_time), 1); lambda(number_of_left_foot_constraints(i_time)+1 : end);];
                        end

                        % save to lists
                        constraint_torque_trajectories_all(i_time, :) = A' * lambda;
                        constraint_torque_trajectories_right(i_time, :) = A' * lambda_right;
                        constraint_torque_trajectories_left(i_time, :) = A' * lambda_left;
                        lambda_trajectories{i_time} = lambda;
                        joint_torque_trajectories(i_time, :) = T;

%                         induced_accelerations_applied_trajectories(i_time, :) = induced_acceleration_applied;
%                         induced_accelerations_gravity_trajectories(i_time, :) = induced_acceleration_gravity;
%                         induced_accelerations_movement_trajectories(i_time, :) = induced_acceleration_movement;
%                         induced_accelerations_applied_single_trajectories(i_time, :, :) = induced_accelerations_applied_single;
                    end
                    display_step = 10;
                    last_time_step = data_points(end);
                    if (i_time / display_step) == floor(i_time / display_step)
                        disp([num2str(i_time) '(' num2str(last_time_step) ')']);
                    end
                end
            end
            fprintf(' done\n');
            toc
        end
        
        %% filter_torques
        if filter_torques
            gap_length_limit = 0.15; % maximum of allowable gap being splined over, in seconds
            time_steps = 1 : number_of_time_steps;
            gap_length_limit_time_steps = gap_length_limit * sampling_rate_mocap;
            gap_template = ones(round(gap_length_limit * sampling_rate_mocap), 1);

            filter_order = 1;
            cutoff_frequency = 3; % cutoff frequency, in Hz
            cutoff_frequency = 10; % cutoff frequency, in Hz
            cutoff_frequency = 20; % cutoff frequency, in Hz
            [b_filter, a_filter] = butter(filter_order, cutoff_frequency/(sampling_rate_mocap/2));	% set filter parameters for butterworth filter

            data_raw = joint_angle_trajectories_belt(:, 1);

            % find the gaps
            gaps = isnan(data_raw);
            gap_start_indices = [];
            gap_end_indices = [];
            for i_time = 1 : number_of_time_steps
                if isnan(data_raw(i_time))
                    % check if this is the start of a gap
                    if (i_time == 1) || (~isnan(data_raw(i_time-1)))
                        gap_start_indices = [gap_start_indices; i_time]; %#ok<*AGROW>
                    end
                    % check if this is the end of a gap
                    if (i_time == number_of_time_steps) || (~isnan(data_raw(i_time+1)))
                        gap_end_indices = [gap_end_indices; i_time];
                    end

                end
            end

            % find the big gaps
            big_gaps = [];
            for i_gap = 1 : length(gap_start_indices)
                % check if this is a long gap
                if ((gap_end_indices(i_gap) - gap_start_indices(i_gap)) > gap_length_limit_time_steps)
                    big_gaps = [big_gaps; [gap_start_indices(i_gap) gap_end_indices(i_gap)]];
                end
            end

            % find the pieces of data without big gaps
            data_stretches_without_big_gaps = [];
            if isempty(big_gaps)
                data_stretches_without_big_gaps = [1 number_of_time_steps];
            else
                if  big_gaps(1, 1) > 1
                    data_stretches_without_big_gaps = [1 big_gaps(1, 1)-1];
                end
                for i_gap = 1 : size(big_gaps, 1) - 1
                    % add the stretch after this big gaps to the list
                    data_stretches_without_big_gaps = [data_stretches_without_big_gaps; big_gaps(i_gap, 2)+1 big_gaps(i_gap+1, 1)-1];
                end
                if ~isempty(big_gaps) && big_gaps(end, 2) < number_of_time_steps
                    data_stretches_without_big_gaps = [data_stretches_without_big_gaps; big_gaps(end, 2)+1 number_of_time_steps];
                end
                % TODO: check if the limit cases are treated correctly

                % make indicator for plotting
                big_gap_indicator = zeros(size(data_raw));
                for i_gap = 1 : size(big_gaps, 1)
                    big_gap_indicator((time_steps >= big_gaps(i_gap, 1)) & (time_steps <= big_gaps(i_gap, 2))) = 1;
                end
                big_gap_indicator(big_gap_indicator == 0) = NaN;
            end        

%             % spline over the gaps
%             underconstrained_data_points = data_points(number_of_lambdas(data_points) < 6);
%             constraint_torques_trajectory_gaps = constraint_torque_trajectories_all;
%             constraint_torques_right_trajectory_gaps = constraint_torque_trajectories_right;
%             constraint_torques_left_trajectory_gaps = constraint_torque_trajectories_left;
%             joint_torques_trajectory_gaps = joint_torque_trajectories;
%             right_ground_reaction_wrench_trajectory_origin_gaps = right_ground_reaction_wrench_trajectory_origin;
%             left_ground_reaction_wrench_trajectory_origin_gaps = left_ground_reaction_wrench_trajectory_origin;
%             induced_accelerations_applied_gaps = induced_accelerations_applied_trajectories;
%             induced_accelerations_gravity_gaps = induced_accelerations_gravity_trajectories;
%             induced_accelerations_movement_gaps = induced_accelerations_movement_trajectories;
%             induced_accelerations_applied_single_gaps = induced_accelerations_applied_single_trajectories;
%             violating_velocity_trajectory_gaps = violating_velocity_trajectories;
%             violating_acceleration_trajectory_gaps = violating_acceleration_trajectories;
%             explained_acceleration_trajectory_gaps = explained_acceleration_trajectories;
%             leftover_acceleration_trajectory_gaps = leftover_acceleration_trajectories;
% 
%             constraint_torques_trajectory_gaps(underconstrained_data_points, :) = NaN;
%             constraint_torques_right_trajectory_gaps(underconstrained_data_points, :) = NaN;
%             constraint_torques_left_trajectory_gaps(underconstrained_data_points, :) = NaN;
%             joint_torques_trajectory_gaps(underconstrained_data_points, :) = NaN;
%             right_ground_reaction_wrench_trajectory_origin_gaps(underconstrained_data_points, :) = NaN;
%             left_ground_reaction_wrench_trajectory_origin_gaps(underconstrained_data_points, :) = NaN;
%             induced_accelerations_applied_gaps(underconstrained_data_points, :) = NaN;
%             induced_accelerations_gravity_gaps(underconstrained_data_points, :) = NaN;
%             induced_accelerations_movement_gaps(underconstrained_data_points, :) = NaN;
%             induced_accelerations_applied_single_gaps(underconstrained_data_points, :, :) = NaN;
%             violating_velocity_trajectory_gaps(underconstrained_data_points, :) = NaN;
%             violating_acceleration_trajectory_gaps(underconstrained_data_points, :) = NaN;
%             explained_acceleration_trajectory_gaps(underconstrained_data_points, :) = NaN;
%             leftover_acceleration_trajectory_gaps(underconstrained_data_points, :) = NaN;
% 
%             constraint_torques_trajectory_splined = zeros(number_of_time_steps, number_of_joints);
%             constraint_torques_right_trajectory_splined = zeros(number_of_time_steps, number_of_joints);
%             constraint_torques_left_trajectory_splined = zeros(number_of_time_steps, number_of_joints);
%             joint_torques_trajectory_splined = zeros(number_of_time_steps, number_of_joints);
%             right_ground_reaction_wrench_trajectory_origin_splined = zeros(number_of_time_steps, 6);
%             left_ground_reaction_wrench_trajectory_origin_splined = zeros(number_of_time_steps, 6);
%             induced_accelerations_applied_splined = zeros(number_of_time_steps, number_of_joints);
%             induced_accelerations_gravity_splined = zeros(number_of_time_steps, number_of_joints);
%             induced_accelerations_movement_splined = zeros(number_of_time_steps, number_of_joints);
%             induced_accelerations_applied_single_splined = zeros(number_of_time_steps, number_of_joints, number_of_joints);
%             violating_velocity_trajectory_splined = zeros(number_of_time_steps, number_of_joints);
%             violating_acceleration_trajectory_splined = zeros(number_of_time_steps, number_of_joints);
%             explained_acceleration_trajectory_splined = zeros(number_of_time_steps, number_of_joints);
%             leftover_acceleration_trajectory_splined = zeros(number_of_time_steps, number_of_joints);
% 
%             for i_joint = 1 : number_of_joints
%                 % spline
%                 [~, constraint_torques_trajectory_splined(:, i_joint)] = evalc('csaps(time, constraint_torques_trajectory_gaps(:, i_joint), 1, time)');
%                 [~, constraint_torques_right_trajectory_splined(:, i_joint)] = evalc('csaps(time, constraint_torques_right_trajectory_gaps(:, i_joint), 1, time)');
%                 [~, constraint_torques_left_trajectory_splined(:, i_joint)] = evalc('csaps(time, constraint_torques_left_trajectory_gaps(:, i_joint), 1, time)');
%                 [~, joint_torques_trajectory_splined(:, i_joint)] = evalc('csaps(time, joint_torques_trajectory_gaps(:, i_joint), 1, time)');
%                 [~, induced_accelerations_applied_splined(:, i_joint)] = evalc('csaps(time, induced_accelerations_applied_gaps(:, i_joint), 1, time)');
%                 [~, induced_accelerations_gravity_splined(:, i_joint)] = evalc('csaps(time, induced_accelerations_gravity_gaps(:, i_joint), 1, time)');
%                 [~, induced_accelerations_movement_splined(:, i_joint)] = evalc('csaps(time, induced_accelerations_movement_gaps(:, i_joint), 1, time)');
%                 for j_joint = 1 : number_of_joints
%                     [~, induced_accelerations_applied_single_splined(:, i_joint, j_joint)] = evalc('csaps(time, induced_accelerations_applied_single_gaps(:, i_joint, j_joint), 1, time)');
%                 end
%                 [~, violating_velocity_trajectory_splined(:, i_joint)] = evalc('csaps(time, violating_velocity_trajectory_gaps(:, i_joint), 1, time)');
%                 [~, violating_acceleration_trajectory_splined(:, i_joint)] = evalc('csaps(time, violating_acceleration_trajectory_gaps(:, i_joint), 1, time)');
%                 [~, explained_acceleration_trajectory_splined(:, i_joint)] = evalc('csaps(time, explained_acceleration_trajectory_gaps(:, i_joint), 1, time)');
%                 [~, leftover_acceleration_trajectory_splined(:, i_joint)] = evalc('csaps(time, leftover_acceleration_trajectory_gaps(:, i_joint), 1, time)');
%             end
%             for i_dof = 1 : 6
%                 [~, right_ground_reaction_wrench_trajectory_origin_splined(:, i_dof)] = evalc('csaps(time, right_ground_reaction_wrench_trajectory_origin_gaps(:, i_dof), 1, time)');
%                 [~, left_ground_reaction_wrench_trajectory_origin_splined(:, i_dof)] = evalc('csaps(time, left_ground_reaction_wrench_trajectory_origin_gaps(:, i_dof), 1, time)');
%             end            
% 
%             % remove parts in big gaps
%             constraint_torques_trajectory_splined(big_gap_indicator==1, :) = NaN;
%             constraint_torques_trajectory_splined(big_gap_indicator==1, :) = NaN;
%             constraint_torques_right_trajectory_splined(big_gap_indicator==1, :) = NaN;
%             constraint_torques_left_trajectory_splined(big_gap_indicator==1, :) = NaN;
%             joint_torques_trajectory_splined(big_gap_indicator==1, :) = NaN;
%             right_ground_reaction_wrench_trajectory_origin_splined(big_gap_indicator==1, :) = NaN;
%             left_ground_reaction_wrench_trajectory_origin_splined(big_gap_indicator==1, :) = NaN;
%             induced_accelerations_applied_splined(big_gap_indicator==1, :) = NaN;
%             induced_accelerations_gravity_splined(big_gap_indicator==1, :) = NaN;
%             induced_accelerations_movement_splined(big_gap_indicator==1, :) = NaN;
%             induced_accelerations_applied_single_splined(big_gap_indicator==1, :, :) = NaN;
%             violating_velocity_trajectory_splined(big_gap_indicator==1, :) = NaN;
%             violating_acceleration_trajectory_splined(big_gap_indicator==1, :) = NaN;
%             explained_acceleration_trajectory_splined(big_gap_indicator==1, :) = NaN;
%             leftover_acceleration_trajectory_splined(big_gap_indicator==1, :) = NaN;


            % filter the stretches of data between the big gaps
            constraint_torque_trajectories_all_smoothed = zeros(number_of_time_steps, number_of_joints);
            constraint_torque_trajectories_left_smoothed = zeros(number_of_time_steps, number_of_joints);
            constraint_torque_trajectories_right_smoothed = zeros(number_of_time_steps, number_of_joints);
            joint_torque_trajectories_smoothed = zeros(number_of_time_steps, number_of_joints);
%             induced_accelerations_applied_smoothed = zeros(number_of_time_steps, number_of_joints);
%             induced_accelerations_gravity_smoothed = zeros(number_of_time_steps, number_of_joints);
%             induced_accelerations_movement_smoothed = zeros(number_of_time_steps, number_of_joints);
%             induced_accelerations_applied_single_smoothed = zeros(number_of_time_steps, number_of_joints, number_of_joints);
%             violating_velocity_trajectory_smoothed = zeros(number_of_time_steps, number_of_joints);
%             violating_acceleration_trajectory_smoothed = zeros(number_of_time_steps, number_of_joints);
%             explained_acceleration_trajectory_smoothed = zeros(number_of_time_steps, number_of_joints);
%             leftover_acceleration_trajectory_smoothed = zeros(number_of_time_steps, number_of_joints);
            for i_stretch = 1 : size(data_stretches_without_big_gaps, 1)
%                 constraint_torque_trajectories_all_smoothed(data_stretches_without_big_gaps(i_stretch, 1) : data_stretches_without_big_gaps(i_stretch, 2), :) = filtfilt(b_filter, a_filter, constraint_torque_trajectories_all(data_stretches_without_big_gaps(i_stretch, 1) : data_stretches_without_big_gaps(i_stretch, 2), :));
                constraint_torque_trajectories_left_smoothed(data_stretches_without_big_gaps(i_stretch, 1) : data_stretches_without_big_gaps(i_stretch, 2), :) = filtfilt(b_filter, a_filter, constraint_torque_trajectories_left(data_stretches_without_big_gaps(i_stretch, 1) : data_stretches_without_big_gaps(i_stretch, 2), :));
                constraint_torque_trajectories_right_smoothed(data_stretches_without_big_gaps(i_stretch, 1) : data_stretches_without_big_gaps(i_stretch, 2), :) = filtfilt(b_filter, a_filter, constraint_torque_trajectories_right(data_stretches_without_big_gaps(i_stretch, 1) : data_stretches_without_big_gaps(i_stretch, 2), :));
                joint_torque_trajectories_smoothed(data_stretches_without_big_gaps(i_stretch, 1) : data_stretches_without_big_gaps(i_stretch, 2), :) = filtfilt(b_filter, a_filter, joint_torque_trajectories(data_stretches_without_big_gaps(i_stretch, 1) : data_stretches_without_big_gaps(i_stretch, 2), :));
            end
            
        end
        
        %% plot_joint_torques
        if plot_joint_torques
            angle_plot_groups = ...
              { ...
                 1 : 6; ...
                 7 : 12; ...
                 13 : 18; ...
                 19 : 24 ...
               };
            number_of_groups = size(angle_plot_groups, 1);
            group_figures = zeros(1, number_of_groups);
            group_axes = zeros(1, number_of_groups);
            
            for i_group = 1 : number_of_groups
                group_figures(i_group) = figure; group_axes(i_group) = axes; hold on
                for i_joint = angle_plot_groups{i_group}
                    plot(time_mocap, joint_torque_trajectories(:, i_joint), 'linewidth', 2, 'displayname', plant.jointLabels{i_joint})
                end
                legend('show', 'location', 'SE')
            end
            
            linkaxes(group_axes, 'x');
            distFig('rows', number_of_groups, 'Position', 'W', 'only', group_figures);
        end

        %% plot_joint_torques
        if plot_joint_torques
            angle_plot_groups = ...
              { ...
                 1 : 6; ...
                 7 : 12; ...
                 13 : 18; ...
                 19 : 24 ...
               };
            number_of_groups = size(angle_plot_groups, 1);
            group_figures_smoothed = zeros(1, number_of_groups);
            group_axes_smoothed = zeros(1, number_of_groups);
            
            for i_group = 1 : number_of_groups
                group_figures_smoothed(i_group) = figure; group_axes_smoothed(i_group) = axes; hold on
                for i_joint = angle_plot_groups{i_group}
                    plot(time_mocap, joint_torque_trajectories_smoothed(:, i_joint), 'linewidth', 2, 'displayname', plant.jointLabels{i_joint})
                end
                legend('show', 'location', 'SE')
            end
            
            linkaxes(group_axes_smoothed, 'x');
            distFig('rows', number_of_groups, 'Position', 'E', 'only', group_figures_smoothed);
        end
        
%         linkaxes([group_axes group_axes_smoothed], 'x');
        
        %% save results
        if save_results
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

%             inverse_dynamics_file_name = createFileName(date, subject_id, study_label, i_trial, label);
            inverse_dynamics_file_name = makeFileName(date, subject_id, 'walking', i_trial, label);
            
            save( ...
                  inverse_dynamics_file_name, ...
                  'inertia_matrix_trajectory', ...
                  'coriolis_matrix_trajectory', ...
                  'constraint_matrix_trajectory', ...
                  'constraint_matrix_dot_trajectory', ...
                  'gravitation_matrix_trajectory', ...
                  'number_of_lambdas', ...
                  'joint_torque_trajectories', ...
                  'constraint_torque_trajectories_all', ...
                  'constraint_torque_trajectories_left', ...
                  'constraint_torque_trajectories_right', ...
                  'left_foot_constraint_number_trajectory', ...
                  'right_foot_constraint_number_trajectory', ...
                  'joint_torque_trajectories_smoothed', ...
                  'constraint_torque_trajectories_all_smoothed', ...
                  'constraint_torque_trajectories_left_smoothed', ...
                  'constraint_torque_trajectories_right_smoothed' ...
                );
%                   'right_force_plate_wrench_origin', ...
%                   'left_force_plate_wrench_origin', ...
%                   'constraint_torques_trajectory_gaps', ...
%                   'constraint_torques_right_trajectory_gaps', ...
%                   'constraint_torques_left_trajectory_gaps', ...
%                   'joint_torques_trajectory_gaps', ...
%                   'induced_accelerations_applied_gaps', ...
%                   'induced_accelerations_gravity_gaps', ...
%                   'induced_accelerations_movement_gaps', ...
%                   'induced_accelerations_applied_single_gaps', ...
%                   'violating_velocity_trajectory_gaps', ...
%                   'violating_acceleration_trajectory_gaps', ...
%                   'explained_acceleration_trajectory_gaps', ...
%                   'leftover_acceleration_trajectory_gaps', ...
%                   'constraint_torques_right_trajectory_splined', ...
%                   'constraint_torques_left_trajectory_splined', ...
%                   'joint_torques_trajectory_splined', ...
%                   'induced_accelerations_applied_splined', ...
%                   'induced_accelerations_gravity_splined', ...
%                   'induced_accelerations_movement_splined', ...
%                   'induced_accelerations_applied_single_splined', ...
%                   'violating_velocity_trajectory_splined', ...
%                   'violating_acceleration_trajectory_splined', ...
%                   'explained_acceleration_trajectory_splined', ...
%                   'leftover_acceleration_trajectory_splined', ...
%                   'constraint_torques_trajectory_smoothed', ...
%                   'constraint_torques_right_trajectory_smoothed', ...
%                   'constraint_torques_left_trajectory_smoothed', ...
%                   'joint_torques_trajectory_smoothed', ...
%                   'induced_accelerations_applied_smoothed', ...
%                   'induced_accelerations_gravity_smoothed', ...
%                   'induced_accelerations_movement_smoothed', ...
%                   'induced_accelerations_applied_single_smoothed', ...
%                   'violating_velocity_trajectory_smoothed', ...
%                   'violating_acceleration_trajectory_smoothed', ...
%                   'explained_acceleration_trajectory_smoothed', ...
%                   'leftover_acceleration_trajectory_smoothed', ...
%                   'right_ground_reaction_wrench_trajectory_origin_smoothed', ...
%                   'left_ground_reaction_wrench_trajectory_origin_smoothed' ...
        end
        disp('done')
    end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        




end







