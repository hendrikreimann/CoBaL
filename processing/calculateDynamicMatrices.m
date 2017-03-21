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

% this function calculates joint and passive torques

% input: 
%
% output:


function calculateDynamicMatrices(varargin)
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'use_parallel', false)
    addParameter(parser, 'constraint', 'point')
    parse(parser, varargin{:})
    use_parallel = parser.Results.use_parallel;
    constraint = parser.Results.constraint;
    
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
    
    %% optimize
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            %% load data
            condition = condition_list{i_condition};
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'kinematicTrajectories')]);
            load(['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'stepEvents')]);
            [left_forceplate_wrench_world_trajectory, time_forceplate] = loadData(date, subject_id, condition, i_trial, 'left_forceplate_wrench_world');
            right_forceplate_wrench_world_trajectory = loadData(date, subject_id, condition, i_trial, 'right_forceplate_wrench_world');
            
            number_of_time_steps = size(marker_trajectories, 1);

            % determine time steps to optimize
            time_steps_to_process = 1 : number_of_time_steps;
            time_steps_to_process = 1001 : 2000;
            
%             time_steps_to_process = determineTimeStepsToOptimize(date, subject_id, condition, i_trial, study_settings.get('data_stretch_padding'));
%             time_steps_to_process = determineTimeStepsToOptimize(date, subject_id, condition, i_trial, 0);

            disp([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Condition ' condition ', Trial ' num2str(i_trial)])
            fprintf([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Calculating dynamic matrices... \n'])
            
            % determine constraint numbers
            left_foot_constraint_number_trajectory = determineConstraintNumbers(time_marker, left_touchdown_times, left_fullstance_times, left_pushoff_times);
            right_foot_constraint_number_trajectory = determineConstraintNumbers(time_marker, right_touchdown_times, right_fullstance_times, right_pushoff_times);
            
            %% calculate belt space transformation
            [belt_speed_left_trajectory, time_belts] = loadData(date, subject_id, condition, i_trial, 'belt_speed_left_trajectory');
            belt_speed_right_trajectory = loadData(date, subject_id, condition, i_trial, 'belt_speed_right_trajectory');
            belt_speed_trajectory_belts = mean([belt_speed_left_trajectory belt_speed_right_trajectory], 2);
            
            belt_speed_trajectory_mocap = spline(time_belts, belt_speed_trajectory_belts, time_mocap)';
            delta_t = 1 / sampling_rate_mocap;
            belt_position_trajectory_mocap = zeros(size(belt_speed_trajectory_mocap));
            for i_time = 2 : length(belt_speed_trajectory_mocap)
                belt_position_trajectory_mocap(i_time) = belt_position_trajectory_mocap(i_time-1) + delta_t * belt_speed_trajectory_mocap(i_time-1);
            end
            
            % transform joint angles and ground reaction wrenches into belt space
            joint_angle_trajectories_belt = joint_angle_trajectories;
            joint_angle_trajectories_belt(:, 2) = joint_angle_trajectories_belt(:, 2) + belt_position_trajectory_mocap';
            
            belt_position_trajectory_forceplate = spline(time_mocap, belt_position_trajectory_mocap, time_forceplate)';
            left_forceplate_wrench_belt_trajectory = zeros(size(left_forceplate_wrench_world_trajectory));
            right_forceplate_wrench_belt_trajectory = zeros(size(right_forceplate_wrench_world_trajectory));
            for i_time = 1 : length(time_forceplate)
                left_forceplate_wrench_world = left_forceplate_wrench_world_trajectory(i_time, :)';
                right_forceplate_wrench_world = right_forceplate_wrench_world_trajectory(i_time, :)';
                
                % define forceplate rotation and translation
                world_to_Awb_rotation = [1 0 0; 0 1 0; 0 0 1];
                world_to_Awb_translation = [0; -belt_position_trajectory_forceplate(i_time); 0];
                world_to_Awb_trafo = [world_to_Awb_rotation world_to_Awb_translation; 0 0 0 1];
                world_to_Awb_adjoint = rigidToAdjointTransformation(world_to_Awb_trafo);

                % transform
                left_forceplate_wrench_belt = (world_to_Awb_adjoint' * left_forceplate_wrench_world);
                right_forceplate_wrench_belt = (world_to_Awb_adjoint' * right_forceplate_wrench_world);
                
                left_forceplate_wrench_belt_trajectory(i_time, :) = left_forceplate_wrench_belt;
                right_forceplate_wrench_belt_trajectory(i_time, :) = right_forceplate_wrench_belt;
            end            
            
%             figure; axes; hold on
%             plot(left_forceplate_wrench_world_trajectory(:, 1:3), 'linewidth', 3);
%             plot(left_forceplate_wrench_belt_trajectory(:, 1:3));
% 
%             figure; axes; hold on
%             plot(left_forceplate_wrench_world_trajectory(:, 4), 'linewidth', 3);
%             plot(left_forceplate_wrench_belt_trajectory(:, 4));
%             
%             figure; axes; hold on
%             plot(left_forceplate_wrench_world_trajectory(:, 5), 'linewidth', 3);
%             plot(left_forceplate_wrench_belt_trajectory(:, 5));
%             
%             figure; axes; hold on
%             plot(left_forceplate_wrench_world_trajectory(:, 6), 'linewidth', 3);
%             plot(left_forceplate_wrench_belt_trajectory(:, 6));
            
            %% derive joint angle trajectories by time
            if study_settings.get('filter_joint_angle_data')
                filter_order = 4;
                cutoff_frequency = study_settings.get('joint_angle_data_cutoff_frequency'); % in Hz
                [b_joint_angle, a_joint_angle] = butter(filter_order, cutoff_frequency/(sampling_rate_mocap/2));
                joint_angle_trajectories_belt = nanfiltfilt(b_joint_angle, a_joint_angle, joint_angle_trajectories_belt);
            end
            joint_velocity_trajectories_belt = deriveByTime(joint_angle_trajectories_belt, time_mocap);
            if study_settings.get('filter_joint_velocity_data')
                filter_order = 4;
                cutoff_frequency = study_settings.get('joint_velocity_data_cutoff_frequency'); % in Hz
                [b_joint_velocity, a_joint_velocity] = butter(filter_order, cutoff_frequency/(sampling_rate_mocap/2));
                joint_velocity_trajectories_belt = nanfiltfilt(b_joint_velocity, a_joint_velocity, joint_velocity_trajectories_belt);
            end
            joint_acceleration_trajectories_belt = deriveByTime(joint_velocity_trajectories_belt, time_mocap);
            
%             figure; axes; hold on
%             plot(joint_velocity_trajectories_belt(:, 1:3))
%             
%             figure; axes; hold on
%             plot(joint_velocity_trajectories_belt(:, 4:6))
%             
%             figure; axes; hold on
%             plot(joint_velocity_trajectories_belt(:, 7:9))
%             
%             figure; axes; hold on
%             plot(joint_velocity_trajectories_belt(:, 10:13))

%             figure; axes; hold on
%             plot(joint_acceleration_trajectories_belt(:, 1:3))
%             
%             figure; axes; hold on
%             plot(joint_acceleration_trajectories_belt(:, 4:6))
%             
%             figure; axes; hold on
%             plot(joint_acceleration_trajectories_belt(:, 7:9))
%             
%             figure; axes; hold on
%             plot(joint_acceleration_trajectories_belt(:, 10:13))

            %% calculate dynamic matrices
            number_of_time_steps = size(joint_angle_trajectories_belt, 1);
            inertia_matrix_trajectory = cell(number_of_time_steps, 1);
            coriolis_matrix_trajectory = cell(number_of_time_steps, 1);
            gravitation_matrix_trajectory = cell(number_of_time_steps, 1);

            constraint_matrix_trajectory = cell(number_of_time_steps, 1);
            constraint_matrix_dot_trajectory = cell(number_of_time_steps, 1);
            number_of_lambdas = zeros(number_of_time_steps, 1);
            number_of_left_foot_constraints = zeros(number_of_time_steps, 1);
            number_of_right_foot_constraints = zeros(number_of_time_steps, 1);
                
            tic
            if use_parallel
                inertia_matrix_trajectory_pool = cell(size(inertia_matrix_trajectory));
                coriolis_matrix_trajectory_pool = cell(size(coriolis_matrix_trajectory));
                gravitation_matrix_trajectory_pool = cell(size(gravitation_matrix_trajectory));
                constraint_matrix_trajectory_pool = cell(size(constraint_matrix_trajectory));
                constraint_matrix_dot_trajectory_pool = cell(size(constraint_matrix_dot_trajectory));
                number_of_lambdas_pool = zeros(number_of_time_steps, 1);
                plant = kinematic_tree;
                spmd
                    plant_pool = plant.copy;
                    for i_time = time_steps_to_process(1)+labindex-1 : numlabs : time_steps_to_process(end)
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

                            if strcmp(constraint, 'point')
                                [constraint_matrix_trajectory_pool{i_time}, constraint_matrix_dot_trajectory_pool{i_time}] = ...
                                    createConstraintMatrix_pointConstraints ...
                                      ( ...
                                        plant_pool, ...
                                        left_foot_constraint_number_trajectory(i_time), ...
                                        right_foot_constraint_number_trajectory(i_time) ...
                                      );
                            end
                            if strcmp(constraint, 'hinge')
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
                            if strcmp(constraint, 'body_velocity')
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
                    inertia_matrix_trajectory(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end)) = inertia_matrix_trajectory_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end));
                    coriolis_matrix_trajectory_lab = coriolis_matrix_trajectory_pool{i_lab};
                    coriolis_matrix_trajectory(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end)) = coriolis_matrix_trajectory_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end));
                    gravitation_matrix_trajectory_lab = gravitation_matrix_trajectory_pool{i_lab};
                    gravitation_matrix_trajectory(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end)) = gravitation_matrix_trajectory_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end));

                    constraint_matrix_trajectory_lab = constraint_matrix_trajectory_pool{i_lab};
                    constraint_matrix_trajectory(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) = constraint_matrix_trajectory_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
                    constraint_matrix_dot_trajectory_lab = constraint_matrix_dot_trajectory_pool{i_lab};
                    constraint_matrix_dot_trajectory(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) = constraint_matrix_dot_trajectory_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);

                    number_of_lambdas_lab = number_of_lambdas_pool{i_lab};
                    number_of_lambdas(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) = number_of_lambdas_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
                end
            else
                for i_time = time_steps_to_process
                    if any(isnan(joint_angle_trajectories_belt(i_time, :)))
                        inertia_matrix_trajectory{i_time} = NaN;
                        coriolis_matrix_trajectory{i_time} = NaN;
                        gravitation_matrix_trajectory{i_time} = NaN;
                        constraint_matrix_trajectory{i_time} = NaN;
                        constraint_matrix_dot_trajectory{i_time} = NaN;
                    else
                        % update model
                        kinematic_tree.jointAngles = joint_angle_trajectories_belt(i_time, :)';
                        kinematic_tree.jointVelocities = joint_velocity_trajectories_belt(i_time, :)';
                        kinematic_tree.jointAccelerations = joint_acceleration_trajectories_belt(i_time, :)';
                        kinematic_tree.updateInternals;
                        
                        % save dynamic matrices
                        inertia_matrix_trajectory{i_time} = kinematic_tree.inertiaMatrix;
                        coriolis_matrix_trajectory{i_time} = kinematic_tree.coriolisMatrix;
                        gravitation_matrix_trajectory{i_time} = kinematic_tree.gravitationalTorqueMatrix;

                        if strcmp(constraint, 'point')
                            [constraint_matrix_trajectory{i_time}, constraint_matrix_dot_trajectory{i_time}] = ...
                                createConstraintMatrix_pointConstraints ...
                                  ( ...
                                    kinematic_tree, ...
                                    left_foot_constraint_number_trajectory(i_time), ...
                                    right_foot_constraint_number_trajectory(i_time) ...
                                  );
                        end
                        if strcmp(constraint, 'hinge')
                            [ ...
                              constraint_matrix_trajectory{i_time}, ...
                              constraint_matrix_dot_trajectory{i_time}, ...
                              number_of_left_foot_constraints(i_time), ...
                              number_of_right_foot_constraints(i_time) ...
                            ] = ...
                            createConstraintMatrix_hingeConstraints ...
                              ( ...
                                kinematic_tree, ...
                                left_foot_constraint_number_trajectory(i_time), ...
                                right_foot_constraint_number_trajectory(i_time) ...
                              );
                        end
                        if strcmp(constraint, 'body_velocity')
                            [constraint_matrix_trajectory{i_time}, constraint_matrix_dot_trajectory{i_time}] = ...
                                createConstraintMatrix_bodyVelocityConstraints ...
                                  ( ...
                                    kinematic_tree, ...
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
                        if strcmp(constraint, 'rollover')
                            [constraint_matrix_trajectory{i_time}, constraint_matrix_dot_trajectory{i_time}] = ...
                                createConstraintMatrix_rolloverConstraints ...
                                  ( ...
                                    kinematic_tree, ...
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
                    last_time_step = time_steps_to_process(end);
                    if (i_time / display_step) == floor(i_time / display_step)
                        disp([num2str(i_time) '(' num2str(last_time_step) ')']);
                    end
                end
            end
            fprintf(' done\n');
            toc
            
            %% save
            variables_to_save = struct;
            variables_to_save.time_mocap = time_mocap;
            variables_to_save.sampling_rate_mocap = sampling_rate_mocap;
            
            variables_to_save.joint_labels = kinematic_tree.jointLabels;
            variables_to_save.joint_angle_trajectories_belt = joint_angle_trajectories_belt;
            variables_to_save.joint_velocity_trajectories_belt = joint_velocity_trajectories_belt;
            variables_to_save.joint_acceleration_trajectories_belt = joint_acceleration_trajectories_belt;
            
            variables_to_save.inertia_matrix_trajectory = inertia_matrix_trajectory;
            variables_to_save.coriolis_matrix_trajectory = coriolis_matrix_trajectory;
            variables_to_save.gravitation_matrix_trajectory = gravitation_matrix_trajectory;
            variables_to_save.constraint_matrix_trajectory = constraint_matrix_trajectory;
            variables_to_save.constraint_matrix_dot_trajectory = constraint_matrix_dot_trajectory;
            variables_to_save.left_foot_constraint_number_trajectory = left_foot_constraint_number_trajectory;
            variables_to_save.right_foot_constraint_number_trajectory = right_foot_constraint_number_trajectory;
            variables_to_save.belt_position_trajectory_mocap = belt_position_trajectory_mocap';
            
            save_folder = 'processed';
            save_file_name = makeFileName(date, subject_id, condition, i_trial, 'dynamicTrajectories.mat');
            saveDataToFile([save_folder filesep save_file_name], variables_to_save);
            disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' save_folder filesep save_file_name]);
            
            addAvailableData('joint_angle_trajectories_belt', 'time_mocap', 'sampling_rate_mocap', 'joint_labels', save_folder, save_file_name);
            addAvailableData('joint_velocity_trajectories_belt', 'time_mocap', 'sampling_rate_mocap', 'joint_labels', save_folder, save_file_name);
            addAvailableData('joint_acceleration_trajectories_belt', 'time_mocap', 'sampling_rate_mocap', 'joint_labels', save_folder, save_file_name);
            addAvailableData('belt_position_trajectory_mocap', 'time_mocap', 'sampling_rate_mocap', '', save_folder, save_file_name);

            addAvailableData('inertia_matrix_trajectory', 'time_mocap', 'sampling_rate_mocap', '', save_folder, save_file_name);
            addAvailableData('coriolis_matrix_trajectory', 'time_mocap', 'sampling_rate_mocap', '', save_folder, save_file_name);
            addAvailableData('gravitation_matrix_trajectory', 'time_mocap', 'sampling_rate_mocap', '', save_folder, save_file_name);
            addAvailableData('constraint_matrix_trajectory', 'time_mocap', 'sampling_rate_mocap', '', save_folder, save_file_name);
            addAvailableData('constraint_matrix_dot_trajectory', 'time_mocap', 'sampling_rate_mocap', '', save_folder, save_file_name);
            
            
            
            
        end
    end
end

function constraint_numbers = determineConstraintNumbers(time, touchdown_times, fullstance_times, pushoff_times)
    constraint_numbers = zeros(size(time));
    for i_time = 1 : length(time)
        time_now = time(i_time);
        last_touchdown = max(touchdown_times(touchdown_times < time_now));
        last_fullstance = max(fullstance_times(fullstance_times < time_now));
        last_pushoff = max(pushoff_times(pushoff_times < time_now));
        next_touchdown = min(touchdown_times(touchdown_times >= time_now));
        next_fullstance = min(fullstance_times(fullstance_times >= time_now));
        next_pushoff = min(pushoff_times(pushoff_times >= time_now));
        
        if isempty(last_touchdown) && isempty(last_fullstance) && isempty(last_pushoff)
            % are we at the beginning?
            constraint_numbers(i_time) = NaN;
        elseif isempty(next_touchdown) && isempty(next_fullstance) && isempty(next_pushoff)
            % last event unknown, can't decide
            constraint_numbers(i_time) = NaN;
        else
            % insert placeholders
            if isempty(last_touchdown)
                last_touchdown = time(1);
            end
            if isempty(last_fullstance)
                last_fullstance = time(1);
            end
            if isempty(last_pushoff)
                last_pushoff = time(1);
            end
            if isempty(next_touchdown)
                next_touchdown = time(end);
            end
            if isempty(next_fullstance)
                next_fullstance = time(end);
            end
            if isempty(next_pushoff)
                next_pushoff = time(end);
            end
            
            % determine last and next events
            [~, last_event_index] = max([last_touchdown, last_fullstance, last_pushoff]);
            [~, next_event_index] = min([next_touchdown, next_fullstance, next_pushoff]);
            
            % determine constraint number
            if last_event_index == 1 && next_event_index == 2
                constraint_numbers(i_time) = 1;
            elseif last_event_index == 2 && next_event_index == 3
                constraint_numbers(i_time) = 2;
            elseif last_event_index == 3 && next_event_index == 1
                constraint_numbers(i_time) = 0;
            else
                constraint_numbers(i_time) = NaN;
            end
        end
        
        
        
    end
end
