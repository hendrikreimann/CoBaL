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


function calculateGroundReactionWrenches(varargin)
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
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            %% load data
            condition = condition_list{i_condition};
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'kinematicTrajectories')]);
            load(['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'stepEvents')]);
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'dynamicTrajectories')]);
            
            number_of_time_steps = size(marker_trajectories, 1);

            % determine time steps to optimize
            time_steps_to_process = 1 : number_of_time_steps;
            time_steps_to_process = 1001 : 2000;
%             time_steps_to_process = 1001 : 1010;
            
%             time_steps_to_optimize = determineTimeStepsToOptimize(date, subject_id, condition, i_trial, study_settings.get('data_stretch_padding'));
%             time_steps_to_optimize = determineTimeStepsToOptimize(date, subject_id, condition, i_trial, 0);

            disp([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Condition ' condition ', Trial ' num2str(i_trial)])
            fprintf([datestr(datetime,'yyyy-mm-dd HH:MM:SS') ' - Calculating ground reaction wrenches ... \n'])
            

            %% calculate
            tic
            left_ground_reaction_wrench_trajectory = zeros(number_of_time_steps, 6);
            right_ground_reaction_wrench_trajectory = zeros(number_of_time_steps, 6);
            if use_parallel
                left_ground_reaction_wrench_trajectory_pool = zeros(size(left_ground_reaction_wrench_trajectory));
                right_ground_reaction_wrench_trajectory_pool = zeros(size(left_ground_reaction_wrench_trajectory));
                kinematic_tree_pool = kinematic_tree.copy;
                joint_velocity_trajectories_belt_pool = joint_velocity_trajectories_belt;
                joint_acceleration_trajectories_belt_pool = joint_acceleration_trajectories_belt;
                constraint_torque_trajectories_left_pool = constraint_torque_trajectories_left;
                constraint_torque_trajectories_right_pool = constraint_torque_trajectories_right;
                belt_position_trajectory_mocap_pool = belt_position_trajectory_mocap;
                spmd
                    for i_time = time_steps_to_process
                        if any(isnan(joint_angle_trajectories_belt(i_time, :)))
                            left_ground_reaction_wrench_trajectory_pool(i_time, :) = NaN;
                            right_ground_reaction_wrench_trajectory_pool(i_time, :) = NaN;
                        else
                            % update model
                            kinematic_tree_pool.jointAngles = joint_angle_trajectories_belt(i_time, :)';
                            kinematic_tree_pool.jointVelocities = joint_velocity_trajectories_belt_pool(i_time, :)';
                            kinematic_tree_pool.jointAccelerations = joint_acceleration_trajectories_belt_pool(i_time, :)';
                            kinematic_tree_pool.updateKinematics;

                            % calculate ground reaction wrenches
                            treadmill_origin_belt = [0; belt_position_trajectory_mocap_pool(i_time); 0]; % this is the physical location of the world coordinate frame in belt coordinates
                            world_to_belt_transformation = [eye(3), treadmill_origin_belt; 0 0 0 1];
                            left_ground_reaction_wrench_trajectory_pool(i_time, :) = calculateInstantaneousGroundReactionWrench(kinematic_tree_pool, constraint_torque_trajectories_left_pool(i_time, :)', world_to_belt_transformation);
                            right_ground_reaction_wrench_trajectory_pool(i_time, :) = calculateInstantaneousGroundReactionWrench(kinematic_tree_pool, constraint_torque_trajectories_right_pool(i_time, :)', world_to_belt_transformation);

                        end
                    end
                end
                % reassemble
                for i_lab = 1 : number_of_labs
                    left_ground_reaction_wrench_trajectory_lab = left_ground_reaction_wrench_trajectory_pool{i_lab};
                    left_ground_reaction_wrench_trajectory(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) = left_ground_reaction_wrench_trajectory_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
                    right_ground_reaction_wrench_trajectory_lab = right_ground_reaction_wrench_trajectory_pool{i_lab};
                    right_ground_reaction_wrench_trajectory(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :) = right_ground_reaction_wrench_trajectory_lab(time_steps_to_process(1)+i_lab-1 : number_of_labs : time_steps_to_process(end), :);
                end
            else
                for i_time = time_steps_to_process
                    if any(isnan(joint_angle_trajectories_belt(i_time, :)))
                        left_ground_reaction_wrench_trajectory(i_time, :) = NaN;
                        right_ground_reaction_wrench_trajectory(i_time, :) = NaN;
                    else
                        % update model
                        kinematic_tree.jointAngles = joint_angle_trajectories_belt(i_time, :)';
                        kinematic_tree.jointVelocities = joint_velocity_trajectories_belt(i_time, :)';
                        kinematic_tree.jointAccelerations = joint_acceleration_trajectories_belt(i_time, :)';
                        kinematic_tree.updateKinematics;

                        % calculate ground reaction wrenches
                        treadmill_origin_belt = [0; belt_position_trajectory_mocap(i_time); 0]; % this is the physical location of the world coordinate frame in belt coordinates
                        world_to_belt_transformation = [eye(3), treadmill_origin_belt; 0 0 0 1];
                        left_ground_reaction_wrench_trajectory(i_time, :) = calculateInstantaneousGroundReactionWrench(kinematic_tree, constraint_torque_trajectories_left(i_time, :)', world_to_belt_transformation);
                        right_ground_reaction_wrench_trajectory(i_time, :) = calculateInstantaneousGroundReactionWrench(kinematic_tree, constraint_torque_trajectories_right(i_time, :)', world_to_belt_transformation);

                    end
                    display_step = 10;
                    last_time_step = time_steps_to_process(end);
                    if (i_time / display_step) == floor(i_time / display_step)
                        disp([num2str(i_time) '(' num2str(last_time_step) ')']);
                    end
                end
            end
            total_ground_reaction_wrench_trajectory = left_ground_reaction_wrench_trajectory + right_ground_reaction_wrench_trajectory;
            fprintf(' done\n');
            toc
            
            %% save
            variables_to_save = struct;
            variables_to_save.left_ground_reaction_wrench_trajectory = left_ground_reaction_wrench_trajectory;
            variables_to_save.right_ground_reaction_wrench_trajectory = right_ground_reaction_wrench_trajectory;
            variables_to_save.total_ground_reaction_wrench_trajectory = total_ground_reaction_wrench_trajectory;

            save_folder = 'processed';
            save_file_name = makeFileName(date, subject_id, condition, i_trial, 'dynamicTrajectories.mat');
            saveDataToFile([save_folder filesep save_file_name], variables_to_save);
            disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' save_folder filesep save_file_name]);
            
            addAvailableData('left_ground_reaction_wrench_trajectory', 'time_mocap', 'sampling_rate_mocap', '', save_folder, save_file_name);
            addAvailableData('right_ground_reaction_wrench_trajectory', 'time_mocap', 'sampling_rate_mocap', '', save_folder, save_file_name);
            addAvailableData('total_ground_reaction_wrench_trajectory', 'time_mocap', 'sampling_rate_mocap', '', save_folder, save_file_name);
    
        end
    end
end











