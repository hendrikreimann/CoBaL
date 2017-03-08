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
    
    %% optimize
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            % load data
            condition = condition_list{i_condition};
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'kinematicTrajectories')]);
            load(['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'stepEvents')]);
            
            number_of_time_steps = size(marker_trajectories, 1);

            % determine time steps to optimize
            time_steps_to_process = 1 : number_of_time_steps;
            time_steps_to_process = 301 : 320;
            
%             time_steps_to_optimize = determineTimeStepsToOptimize(date, subject_id, condition, i_trial, study_settings.get('data_stretch_padding'));
%             time_steps_to_optimize = determineTimeStepsToOptimize(date, subject_id, condition, i_trial, 0);

            disp([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Condition ' condition ', Trial ' num2str(i_trial)])
            fprintf([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Calculating ... \n'])
            
            % determine constraint numbers
            left_foot_constraint_number_trajectory = determineConstraintNumbers(time_marker, left_touchdown_times, left_fullstance_times, left_pushoff_times);
            right_foot_constraint_number_trajectory = determineConstraintNumbers(time_marker, right_touchdown_times, right_fullstance_times, right_pushoff_times);

            % 
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










