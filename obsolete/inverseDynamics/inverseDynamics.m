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


function inverseDynamics(varargin)
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
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'dynamicTrajectories')]);
            [left_forceplate_wrench_world_trajectory, time_forceplate] = loadData(date, subject_id, condition, i_trial, 'left_forceplate_wrench_world');
            right_forceplate_wrench_world_trajectory = loadData(date, subject_id, condition, i_trial, 'right_forceplate_wrench_world');
            
            number_of_time_steps = size(marker_trajectories, 1);

            % determine time steps to optimize
%             time_steps_to_process = 1 : number_of_time_steps;
%             time_steps_to_process = 1001 : 2000;
            time_steps_to_process = determineTimeStepsToProcess(date, subject_id, condition, i_trial, study_settings.get('data_stretch_padding'));

            % remove time steps to process based on available events
            bad_time_steps_left = time_steps_to_process(isnan(left_foot_constraint_number_trajectory(time_steps_to_process)));
            bad_time_steps_right = time_steps_to_process(isnan(right_foot_constraint_number_trajectory(time_steps_to_process)));
            bad_time_steps = [bad_time_steps_left bad_time_steps_right];
            
            % ignore the bad time steps and re-determinate the time steps to process
            if ~isempty(bad_time_steps)
                bad_times = time_mocap(bad_time_steps);
                variables_to_save = struct;
                if exist('ignore_times', 'var')
                    variables_to_save.ignore_times = [ignore_times; bad_times];
                else
                    variables_to_save.ignore_times = bad_times;
                end
                
                
                save_folder = 'analysis';
                save_file_name = makeFileName(date, subject_id, condition, i_trial, 'stepEvents.mat');
                saveDataToFile([save_folder filesep save_file_name], variables_to_save);

                findRelevantDataStretches('condition', {condition}, 'trial', i_trial);
                time_steps_to_process = determineTimeStepsToProcess(date, subject_id, condition, i_trial, study_settings.get('data_stretch_padding'));
            end   
            
%             disp([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Condition ' condition ', Trial ' num2str(i_trial)])
%             fprintf([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Calculating ... \n'])
            disp([' - Condition ' condition ', Trial ' num2str(i_trial)])
            fprintf([' - Calculating ... \n'])
            

            %% calculate_torques
            tic
            fprintf('Calculating the constraint torques ... ');
            number_of_joints = kinematic_tree.numberOfJoints;
            number_of_time_steps = size(joint_angle_trajectories_belt, 1);
            virtual_joints = 1:6;
            constraint_torque_trajectories_all = zeros(number_of_time_steps, number_of_joints);
            constraint_torque_trajectories_right = zeros(number_of_time_steps, number_of_joints);
            constraint_torque_trajectories_left = zeros(number_of_time_steps, number_of_joints);
            lambda_trajectories = cell(number_of_time_steps, 1);
            joint_torque_trajectories = zeros(number_of_time_steps, number_of_joints);
%             induced_accelerations_applied_trajectories = zeros(number_of_time_steps, number_of_joints);
%             induced_accelerations_gravity_trajectories = zeros(number_of_time_steps, number_of_joints);
%             induced_accelerations_movement_trajectories = zeros(number_of_time_steps, number_of_joints);
%             induced_accelerations_applied_single_trajectories = zeros(number_of_time_steps, number_of_joints, number_of_joints); % (i_time, i, j)-th entry holds acceleration of the i_th joint from torque at the j_th joint
%             violating_velocity_trajectories = zeros(number_of_time_steps, number_of_joints);
%             violating_acceleration_trajectories = zeros(number_of_time_steps, number_of_joints);
%             explained_acceleration_trajectories = zeros(number_of_time_steps, number_of_joints);
%             leftover_acceleration_trajectories = zeros(number_of_time_steps, number_of_joints);


            if use_parallel
                constraint_torque_trajectories_all_pool = zeros(number_of_time_steps, number_of_joints);
                constraint_torque_trajectories_right_pool = zeros(number_of_time_steps, number_of_joints);
                constraint_torque_trajectories_left_pool = zeros(number_of_time_steps, number_of_joints);
                lambda_trajectories_pool = cell(number_of_time_steps, 1);
                joint_torque_trajectories_pool = zeros(number_of_time_steps, number_of_joints);
                % make variables accessible to workers by declaring them
                left_foot_constraint_number_trajectory_pool = left_foot_constraint_number_trajectory;
                right_foot_constraint_number_trajectory_pool = right_foot_constraint_number_trajectory;
                kinematic_tree_pool = kinematic_tree.copy;

% %                 induced_accelerations_applied_trajectories_pool = zeros(number_of_time_steps, number_of_joints);
% %                 induced_accelerations_gravity_trajectories_pool = zeros(number_of_time_steps, number_of_joints);
% %                 induced_accelerations_movement_trajectories_pool = zeros(number_of_time_steps, number_of_joints);
% %                 induced_accelerations_applied_single_trajectories_pool = zeros(number_of_time_steps, number_of_joints, number_of_joints); 
% %                 violating_velocity_trajectories_pool = zeros(number_of_time_steps, number_of_joints);
% %                 violating_acceleration_trajectories_pool = zeros(number_of_time_steps, number_of_joints);
% %                 explained_acceleration_trajectories_pool = zeros(number_of_time_steps, number_of_joints);
% %                 leftover_acceleration_trajectories_pool = zeros(number_of_time_steps, number_of_joints);


                spmd
                    for i_time = time_steps_to_process(1)+labindex-1 : numlabs : time_steps_to_process(end)
                        if any(isnan(joint_angle_trajectories_belt(i_time, :))) | any(isnan(constraint_matrix_trajectory{i_time}))
                            constraint_torque_trajectories_all_pool(i_time, :) = NaN;
                            constraint_torque_trajectories_right_pool(i_time, :) = NaN;
                            constraint_torque_trajectories_left_pool(i_time, :) = NaN;
                            lambda_trajectories_pool{i_time} = NaN;
                            joint_torque_trajectories_pool(i_time, :) = NaN;
%                             induced_accelerations_applied_trajectories_pool(i_time, :) = NaN;
%                             induced_accelerations_gravity_trajectories_pool(i_time, :) = NaN;
%                             induced_accelerations_movement_trajectories_pool(i_time, :) = NaN;
%                             induced_accelerations_applied_single_trajectories_pool(i_time, :, :) = NaN;
%                             violating_velocity_trajectories_pool(i_time, :, :) = NaN;
%                             violating_acceleration_trajectories_pool(i_time, :, :) = NaN;
%                             explained_acceleration_trajectories_pool(i_time, :, :) = NaN;
%                             leftover_acceleration_trajectories_pool(i_time, :, :) = NaN;
                        else
                            M = inertia_matrix_trajectory{i_time};
                            C = coriolis_matrix_trajectory{i_time};
                            N = gravitation_matrix_trajectory{i_time};
                            theta_dot = joint_velocity_trajectories_belt(i_time, :)';
                            theta_two_dot = joint_acceleration_trajectories_belt(i_time, :)';
                            
                            % get constraint matrices
                            A = constraint_matrix_trajectory{i_time};
                            A_dot = constraint_matrix_dot_trajectory{i_time};

                            P = eye(kinematic_tree_pool.numberOfJoints) - A' * (A * M^(-1) * A')^(-1) * A * M^(-1);

                            H_2 = P;
                            b_2 = (M*theta_two_dot + P*C*theta_dot + P*N + A'*(A*M^(-1)*A')^(-1)*(A_dot * theta_dot));

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
                            b = ...
                              [ ...
                                b_2; ...
                                zeros(k_v, 1); ...
                                zeros(number_of_joints - k_w, 1); ...
                              ];

                            T = pinv(H) * b;
                            lambda = (A*M^(-1)*A')^(-1) *  A*M^(-1)*(T - C*theta_dot - N) + (A*M^(-1)*A')^(-1)*A_dot*theta_dot; % corrected here
                            
                            if strcmp(constraint, 'point')
                                if left_foot_constraint_number_trajectory_pool(i_time) == 0
                                    number_of_left_foot_constraints = 0;
                                elseif left_foot_constraint_number_trajectory_pool(i_time) == 1
                                    number_of_left_foot_constraints = 3;
                                elseif left_foot_constraint_number_trajectory_pool(i_time) == 2
                                    number_of_left_foot_constraints = 3;
                                else
                                    error('Left foot constraint number must be an integer between 0 and 2')
                                end

                                if right_foot_constraint_number_trajectory_pool(i_time) == 0
                                    number_of_right_foot_constraints = 0;
                                elseif right_foot_constraint_number_trajectory_pool(i_time) == 1
                                    number_of_right_foot_constraints = 3;
                                elseif right_foot_constraint_number_trajectory_pool(i_time) == 2
                                    number_of_right_foot_constraints = 3;
                                else
                                    error('Right foot constraint number must be an integer between 0 and 2')
                                end
                            end
                            if number_of_left_foot_constraints + number_of_right_foot_constraints == 0
                                lambda_right = 0;
                                lambda_left = 0;
                            else
                                lambda_left = [lambda(1:number_of_left_foot_constraints); zeros(number_of_right_foot_constraints, 1)];
                                lambda_right = [zeros(number_of_left_foot_constraints, 1); lambda(number_of_left_foot_constraints+1 : end);];
                            end

                            % save to lists
                            constraint_torque_trajectories_all_pool(i_time, :) = A' * lambda;
                            constraint_torque_trajectories_right_pool(i_time, :) = A' * lambda_right;
                            constraint_torque_trajectories_left_pool(i_time, :) = A' * lambda_left;
                            lambda_trajectories_pool{i_time} = lambda;
                            joint_torque_trajectories_pool(i_time, :) = T;

%                             induced_accelerations_applied_trajectories_pool(i_time, :) = induced_acceleration_applied;
%                             induced_accelerations_gravity_trajectories_pool(i_time, :) = induced_acceleration_gravity;
%                             induced_accelerations_movement_trajectories_pool(i_time, :) = induced_acceleration_movement;
%                             induced_accelerations_applied_single_trajectories_pool(i_time, :, :) = induced_accelerations_applied_single;
% 
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
% %                     right_ground_reaction_wrench_trajectory_origin_lab = right_ground_reaction_wrench_trajectory_origin_pool{i_lab};
% %                     left_ground_reaction_wrench_trajectory_origin_lab = left_ground_reaction_wrench_trajectory_origin_pool{i_lab};
%                     induced_accelerations_applied_trajectories_lab = induced_accelerations_applied_trajectories_pool{i_lab};
%                     induced_accelerations_gravity_trajectories_lab = induced_accelerations_gravity_trajectories_pool{i_lab};
%                     induced_accelerations_movement_trajectories_lab = induced_accelerations_movement_trajectories_pool{i_lab};
%                     induced_accelerations_applied_single_trajectories_lab = induced_accelerations_applied_single_trajectories_pool{i_lab};
%                     violating_velocity_trajectories_lab = violating_velocity_trajectories_pool{i_lab};
%                     violating_acceleration_trajectories_lab = violating_acceleration_trajectories_pool{i_lab};
% %                     explained_acceleration_trajectories_lab = explained_acceleration_trajectories_pool{i_lab};
% %                     leftover_acceleration_trajectories_lab = leftover_acceleration_trajectories_pool{i_lab};
% 
                    constraint_torque_trajectories_all(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
                        = constraint_torque_trajectories_all_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
                    constraint_torque_trajectories_right(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
                        = constraint_torque_trajectories_right_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
                    constraint_torque_trajectories_left(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
                        = constraint_torque_trajectories_left_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
                    lambda_trajectories(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
                        = lambda_trajectories_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
                    joint_torque_trajectories(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
                        = joint_torque_trajectories_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
%                     induced_accelerations_applied_trajectories(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
%                         = induced_accelerations_applied_trajectories_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
%                     induced_accelerations_gravity_trajectories(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
%                         = induced_accelerations_gravity_trajectories_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
%                     induced_accelerations_movement_trajectories(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
%                         = induced_accelerations_movement_trajectories_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
%                     induced_accelerations_applied_single_trajectories(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :, :) ...
%                         = induced_accelerations_applied_single_trajectories_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :, :);
%                     violating_velocity_trajectories(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
%                         = violating_velocity_trajectories_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
%                     violating_acceleration_trajectories(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
%                         = violating_acceleration_trajectories_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
% %                     explained_acceleration_trajectories(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
% %                         = explained_acceleration_trajectories_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
% %                     leftover_acceleration_trajectories(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
% %                         = leftover_acceleration_trajectories_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
                end
            else
                for i_time = time_steps_to_process
                    if any(isnan(joint_angle_trajectories_belt(i_time, :))) | any(isnan(constraint_matrix_trajectory{i_time}))
                        constraint_torque_trajectories_all(i_time, :) = NaN;
                        constraint_torque_trajectories_right(i_time, :) = NaN;
                        constraint_torque_trajectories_left(i_time, :) = NaN;
                        lambda_trajectories{i_time} = NaN;
                        joint_torque_trajectories(i_time, :) = NaN;
                        right_ground_reaction_wrench_trajectory_origin(i_time, :) = NaN;
                        left_ground_reaction_wrench_trajectory_origin(i_time, :) = NaN;
%                         induced_accelerations_applied_trajectories(i_time, :) = NaN;
%                         induced_accelerations_gravity_trajectories(i_time, :) = NaN;
%                         induced_accelerations_movement_trajectories(i_time, :) = NaN;
%                         induced_accelerations_applied_single_trajectories(i_time, :, :) = NaN;
%                         violating_velocity_trajectories(i_time, :, :) = NaN;
%                         violating_acceleration_trajectories(i_time, :, :) = NaN;
    %                     explained_acceleration_trajectories(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
    %                         = explained_acceleration_trajectories_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
    %                     leftover_acceleration_trajectories(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) ...
    %                         = leftover_acceleration_trajectories_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
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
                        

                        % begin new stuff
                        P = eye(kinematic_tree.numberOfJoints) - A' * (A * M^(-1) * A')^(-1) * A * M^(-1);

%                         % force plate data
%                         [~, time_index_forceplate] = min(abs(time_forceplate - time_mocap(i_time)));
%                         combined_forceplate_wrench_data = total_forceplate_wrench_Acw(time_index_forceplate, :)';
%                         if number_of_left_foot_constraints(i_time) ~= 0 & number_of_right_foot_constraints(i_time) == 0
%                             Q_left = eye(number_of_left_foot_constraints(i_time));
%                             R_left = pinv(kinematic_tree.calculateArbitraryFrameBodyJacobian(eye(4, 4), 12)');
%                             H_1_left = R_left * A' * Q_left * (A*M^(-1)*A')^(-1)*A*M^(-1);
%                             b_1_left = combined_forceplate_wrench_data - R_left * A' * Q_left * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*kinematic_tree.jointVelocities - N) + A_dot*kinematic_tree.jointVelocities);
%                             
%                             H_1 = H_1_left;
%                             b_1 = b_1_left;
%                         elseif number_of_left_foot_constraints(i_time) == 0 & number_of_right_foot_constraints(i_time) ~= 0
%                             Q_right = eye(number_of_right_foot_constraints(i_time));
%                             R_right = pinv(kinematic_tree.calculateArbitraryFrameBodyJacobian(eye(4, 4), 18)');
%                             H_1_right = R_right * A' * Q_right * (A*M^(-1)*A')^(-1)*A*M^(-1);
%                             b_1_right = combined_forceplate_wrench_data - R_right * A' * Q_right * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*kinematic_tree.jointVelocities - N) + A_dot*kinematic_tree.jointVelocities);
%                             
%                             H_1 = H_1_right;
%                             b_1 = b_1_right;
%                         elseif number_of_left_foot_constraints(i_time) ~= 0 & number_of_right_foot_constraints(i_time) ~= 0
%                             Q_left = [eye(number_of_left_foot_constraints(i_time)) zeros(number_of_right_foot_constraints(i_time), number_of_left_foot_constraints(i_time)); zeros(number_of_left_foot_constraints(i_time), number_of_right_foot_constraints(i_time)) zeros(number_of_right_foot_constraints(i_time))];
%                             R_left = pinv(kinematic_tree.calculateArbitraryFrameBodyJacobian(eye(4, 4), 12)');
%                             H_1_left = R_left * A' * Q_left * (A*M^(-1)*A')^(-1)*A*M^(-1);
% 
%                             Q_right = [zeros(number_of_left_foot_constraints(i_time)) zeros(number_of_right_foot_constraints(i_time), number_of_left_foot_constraints(i_time)); zeros(number_of_left_foot_constraints(i_time), number_of_right_foot_constraints(i_time)) eye(number_of_right_foot_constraints(i_time))];
%                             R_right = pinv(kinematic_tree.calculateArbitraryFrameBodyJacobian(eye(4, 4), 18)');
%                             H_1_right = R_right * A' * Q_right * (A*M^(-1)*A')^(-1)*A*M^(-1);
%                             
%                             b_1 = combined_forceplate_wrench_data ...
%                                   - R_left * A' * Q_left * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*kinematic_tree.jointVelocities - N) + A_dot*kinematic_tree.jointVelocities) ...
%                                   - R_right * A' * Q_right * (A*M^(-1)*A')^(-1) * (A*M^(-1)*(- C*kinematic_tree.jointVelocities - N) + A_dot*kinematic_tree.jointVelocities);
%                         end
                        
                        
                        % acceleration data
                        

                        % this gets the same results as the old version

                        H_2 = P;
                        b_2 = (M*theta_two_dot + P*C*theta_dot + P*N + A'*(A*M^(-1)*A')^(-1)*(A_dot * theta_dot));

                        % according to equation - new version
%                         H_2 = M^(-1)*P;
%                         b_2 = theta_two_dot + M^(-1) * P*(C*theta_dot + N) + M^(-1) * A'*(A*M^(-1)*A')^(-1)*A_dot*theta_dot;
                        
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

                        T = pinv(H) * b;

%                         % weight matrix
%                         acceleration_weights = [1e-3 1e-2 1e-1 1 1e1 1e2 1e3];
%                         T_weighted = zeros(size(T, 1), length(acceleration_weights));
%                         for i_weight = 1 : length(acceleration_weights)
%                             W_inv = diag([ones(1, size(H_2, 1))*acceleration_weights(i_weight), ones(1, k_v + number_of_joints - k_w)]);
%                             H_weighted_pinv = (H'*W_inv*H)^(-1)*H'*W_inv;
%                             T_weighted(:, i_weight) = H_weighted_pinv * b;
%                         end
                        
                        
% status 31.8.2016: I realized that the old code probably had an error, where I would use Q = M*theta_two_dot + P*C*theta_dot + P*N - A'*(A*M^(-1)*A')^(-1)*(A * theta_two_dot);
% with A * theta_two_dot instead of A_dot * theta_dot, and wrong sign. The corrected version does not change the
% calculated ground reaction forces much. 

% Now, add the GRF constraint for a stretch with only one contact

                        
                        
                        
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
% %                             P = eye(kinematic_tree.numberOfJoints) - A' * (A * M^(-1) * A')^(-1) * A * M^(-1);
% %                             Q = M*theta_two_dot + P*C*theta_dot + P*N - A'*(A*M^(-1)*A')^(-1)*(A * theta_two_dot);
% %                         else
% %                             P = eye(kinematic_tree.numberOfJoints);
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

                        % new version, adapted from exploreConstrainedMovement script
                        b_acc = (M*theta_two_dot + P*C*theta_dot + P*N + A'*(A*M^(-1)*A')^(-1)*(A_dot * theta_dot)); % same as b_2 above)
                        
                        %% solution two - contact constraints and no torques in the free DoFs
                        H = ...
                          [ ...
                            P; ...
                            eye(k_v, number_of_joints); ...
                          ];
                        b = ...
                          [ ...
                            b_acc; ...
                            zeros(k_v, 1); ...
                          ];

                        T_2 = pinv(H) * b;
                        lambda = (A*M^(-1)*A')^(-1) ...
                            * (A*M^(-1)*(T_2 - C*theta_dot - N) + A_dot*theta_dot);
                        kinematic_tree.externalTorques = T_2;
                        kinematic_tree.constraintTorques = A'*lambda;
                        kinematic_tree.calculateAccelerationsFromExternalTorques;
                        accelerations_check = kinematic_tree.jointAccelerations;

                        % check for workless torques
                        k_constraints_full = rank(A);
                        [~, ~, W_constraints] = svd(A);                                             % use singular value decomposition to get the basis
                        B_constraints = W_constraints(:, 1:k_constraints_full);                     % E_rng contains the base vectors of the range space - the torque combinations that do no work
                        C_constraints = W_constraints(:, k_constraints_full+1:end);                 % E_nul contains the base vectors of the null space - the torque combinations that do work
                        P_constraint_full = B_constraints*B_constraints';                           % K_rng projects onto the range space
                        P_constraint_free = C_constraints*C_constraints';                           % K_nul projects onto the null space
                        T_workpure_2 = P_constraint_free * T_2;
                        T_workless_2 = P_constraint_full * T_2;
                        T_check_2 = T_workpure_2 + T_workless_2;

                        % check for workless torques - but only those that also have no virtual torques
                        B = [eye(k_v) zeros(k_v, number_of_joints - k_v)];
                        k_c = rank(A);
                        [~, ~, V_c] = svd(A);
                        C_c = V_c(:, k_c+1:end); % C_c spans the null space of A
                        k_w = rank([B' C_c]);
                        [~, ~, V_w] = svd([B' C_c]');
                        B_w = V_w(:, 1:k_w);                                                        % B_w columns span the space of torque vectors that do work, or have a virtual component, or both
                        C_w = V_w(:, k_w+1:end);                                                    % C_w columns span the space of torque vectors that do no work and have no virtual component                        
                        P_workpure = B_w*B_w';                                                      % P_workpure projects onto the space spanned by B_w columns
                        P_workless = C_w*C_w';                                                      % P_workless projects onto the space spanned by C_w columns
                        T_workpure_2 = P_workpure * T_2;
                        T_workless_2 = P_workless * T_2;
                        T_check_2 = T_workpure_2 + T_workless_2;

                        % check what happens if I remove workless torques
                        lambda = (A*M^(-1)*A')^(-1) ...
                            * (A*M^(-1)*(T_workpure_2 - C*theta_dot - N) + A_dot*theta_dot);
                        kinematic_tree.externalTorques = T_workpure_2;
                        kinematic_tree.constraintTorques = A'*lambda;
                        kinematic_tree.calculateAccelerationsFromExternalTorques;
                        accelerations_check_2 = kinematic_tree.jointAccelerations;

                        %% solution four - contact constraints, no torques in the free DoFs and no workless torques

                        B = eye(k_v, number_of_joints);
                        k_c = rank(A);
                        [~, ~, V_c] = svd(A);
                        C_c = V_c(:, k_c+1:end); % C_c spans the null space of A
                        k_w = rank([B' C_c]);
                        [~, ~, V_w] = svd([B' C_c]');
                        C_w = V_w(:, k_w+1:end);                        
                        % combine
                        H = ...
                          [ ...
                            P; ...
                            B; ...
                            C_w'; ...
                          ];


                        % this does not work, apparently I have either virtual torques or workless torques

                        b = ...
                          [ ...
                            b_acc; ...
                            zeros(k_v, 1); ...
                            zeros(number_of_joints - k_w, 1); ...
                          ];


                        T_4 = pinv(H) * b;
                        lambda = (A*M^(-1)*A')^(-1) ...
                            * (A*M^(-1)*(T_4 - C*theta_dot - N) + A_dot*theta_dot);
                        kinematic_tree.externalTorques = T_4;
                        kinematic_tree.constraintTorques = A'*lambda;
                        kinematic_tree.calculateAccelerationsFromExternalTorques;
                        accelerations_check = kinematic_tree.jointAccelerations;

                        % check for workless torques
                        k_constraints_full = rank(A);
                        [~, ~, W_constraints] = svd(A);                                             % use singular value decomposition to get the basis
                        B_constraints = W_constraints(:, 1:k_constraints_full);                     % E_rng contains the base vectors of the range space - the torque combinations that do no work
                        C_constraints = W_constraints(:, k_constraints_full+1:end);                 % E_nul contains the base vectors of the null space - the torque combinations that do work
                        P_constraint_full = B_constraints*B_constraints';                           % P_constraint_full projects onto the range space
                        P_constraint_free = C_constraints*C_constraints';                           % P_constraint_free projects onto the null space
                        T_workpure_4 = P_constraint_free * T_4;
                        T_workless_4 = P_constraint_full * T_4;
                        T_check_4 = T_workpure_4 + T_workless_4;

                        % check for workless torques - but only those that also have no virtual torques
                        B = [eye(k_v) zeros(k_v, number_of_joints - k_v)];
                        k_c = rank(A);
                        [~, ~, V_c] = svd(A);
                        C_c = V_c(:, k_c+1:end); % C_c spans the null space of A
                        k_w = rank([B' C_c]);
                        [~, ~, V_w] = svd([B' C_c]');
                        B_w = V_w(:, 1:k_w);                                                        % B_w columns span the space of torque vectors that do work, or have a virtual component, or both
                        C_w = V_w(:, k_w+1:end);                                                    % C_w columns span the space of torque vectors that do no work and have no virtual component                        
                        P_workpure = B_w*B_w';                                                      % P_workpure projects onto the space spanned by B_w columns
                        P_workless = C_w*C_w';                                                      % P_workless projects onto the space spanned by C_w columns
                        T_workpure_4 = P_workpure * T_4;
                        T_workless_4 = P_workless * T_4;
                        T_check_4 = T_workpure_4 + T_workless_4;

                        % check what happens if I remove workless torques
                        lambda = (A*M^(-1)*A')^(-1) ...
                            * (A*M^(-1)*(T_workpure_4 - C*theta_dot - N) + A_dot*theta_dot);
                        kinematic_tree.externalTorques = T_workpure_4;
                        kinematic_tree.constraintTorques = A'*lambda;
                        kinematic_tree.calculateAccelerationsFromExternalTorques;
                        accelerations_check_4 = kinematic_tree.jointAccelerations;




                        if strcmp(constraint, 'point')
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
                        if number_of_left_foot_constraints + number_of_right_foot_constraints == 0
                            lambda_right = 0;
                            lambda_left = 0;
                        else
                            lambda_left = [lambda(1:number_of_left_foot_constraints); zeros(number_of_right_foot_constraints, 1)];
                            lambda_right = [zeros(number_of_left_foot_constraints, 1); lambda(number_of_left_foot_constraints+1 : end);];
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
                    display_step = 100;
                    last_time_step = time_steps_to_process(end);
                    if (i_time / display_step) == floor(i_time / display_step)
                        disp([num2str(i_time) '(' num2str(last_time_step) ')']);
                    end
                end
            end
            fprintf(' done\n');
            toc
            
            %% filter joint torques
            if study_settings.get('filter_joint_torques')
                joint_torque_trajectories_unfiltered = joint_torque_trajectories;
                filter_order = 4;
                cutoff_frequency = study_settings.get('joint_torques_cutoff_frequency'); % in Hz
                [b_marker, a_marker] = butter(filter_order, cutoff_frequency/(sampling_rate_mocap/2));
                joint_torque_trajectories = nanfiltfilt(b_marker, a_marker, joint_torque_trajectories);
            end
            
            
            
            %% save
            variables_to_save = struct;
            
            variables_to_save.joint_torque_trajectories = joint_torque_trajectories;
            variables_to_save.constraint_torque_trajectories_all = constraint_torque_trajectories_all;
            variables_to_save.constraint_torque_trajectories_left = constraint_torque_trajectories_left;
            variables_to_save.constraint_torque_trajectories_right = constraint_torque_trajectories_right;
            
            save_folder = 'processed';
            save_file_name = makeFileName(date, subject_id, condition, i_trial, 'dynamicTrajectories.mat');
            saveDataToFile([save_folder filesep save_file_name], variables_to_save);
            disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' save_folder filesep save_file_name]);
            
            addAvailableData('joint_torque_trajectories', 'time_mocap', 'sampling_rate_mocap', 'joint_labels', save_folder, save_file_name);
            addAvailableData('constraint_torque_trajectories_all', 'time_mocap', 'sampling_rate_mocap', '', save_folder, save_file_name);
            addAvailableData('constraint_torque_trajectories_left', 'time_mocap', 'sampling_rate_mocap', '', save_folder, save_file_name);
            addAvailableData('constraint_torque_trajectories_right', 'time_mocap', 'sampling_rate_mocap', '', save_folder, save_file_name);
    
        end
    end
end











