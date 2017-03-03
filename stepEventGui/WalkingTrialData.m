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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% compare the kinematic tree against the kinematic chain

classdef WalkingTrialData < handle
    properties
        % administration
        data_directory = [];
        condition = [];
        trial_number = [];
        recording_time = 0;
        subject_settings;
        date = [];
        subject_id = [];
        sampling_rate_marker = -1;
        sampling_rate_forceplate = -1;
        kinematic_data_available = 0;
        com_data_available = 0;
        joint_angle_data_available = 0;
        forceplate_data_available = 0;
        
        % time
        time_marker = [];
        time_forceplate = [];
        time_labview = [];
        selected_time = [];
        
        % marker kinematics
        marker_positions = [];
        marker_labels = [];
        joint_center_positions = [];
        joint_center_labels = [];
        com_positions = [];
        com_labels = [];
        joint_angles = [];
        
        left_heel_y_pos = [];
        left_heel_y_vel = [];
        left_heel_y_acc = [];

        right_heel_y_pos = [];
        right_heel_y_vel = [];
        right_heel_y_acc = [];

        left_heel_z_pos = [];
        left_heel_z_vel = [];
        left_heel_z_acc = [];

        right_heel_z_pos = [];
        right_heel_z_vel = [];
        right_heel_z_acc = [];

        left_toes_z_pos = [];
        left_toes_z_vel = [];
        left_toes_z_acc = [];

        right_toes_z_pos = [];
        right_toes_z_vel = [];
        right_toes_z_acc = [];
        
        left_fx = [];
        left_fy = [];
        left_fz = [];
        left_mx = [];
        left_my = [];
        left_mz = [];
        
        right_fx = [];
        right_fy = [];
        right_fz = [];
        right_mx = [];
        right_my = [];
        right_mz = [];
        
        
        left_foot_fx = [];
        left_foot_fy = [];
        left_foot_fz = [];
        left_foot_mx = [];
        left_foot_my = [];
        left_foot_mz = [];
        
        right_foot_fx = [];
        right_foot_fy = [];
        right_foot_fz = [];
        right_foot_mx = [];
        right_foot_my = [];
        right_foot_mz = [];
        stimulus_state = [];
        
        left_arm_angle = [];
        right_arm_angle = [];
        left_leg_angle = [];
        right_leg_angle = [];
        
        left_arm_phase = [];
        right_arm_phase = [];
        left_leg_phase = [];
        right_leg_phase = [];
        
        data_labels_mocap = ...
          { ...
            'left_heel_z_pos', 'left_heel_z_vel', 'left_heel_z_acc', ...
            'left_toes_z_pos', 'left_toes_z_vel', 'left_toes_z_acc', ...
            'right_heel_z_pos', 'right_heel_z_vel', 'right_heel_z_acc', ...
            'right_toes_z_pos', 'right_toes_z_vel', 'right_toes_z_acc', ...
            'left_heel_y_pos', 'left_heel_y_vel', 'left_heel_y_acc', ...
            'left_toes_y_pos', 'left_toes_y_vel', 'left_toes_y_acc', ...
            'right_heel_y_pos', 'right_heel_y_vel', 'right_heel_y_acc', ...
            'right_toes_y_pos', 'right_toes_y_vel', 'right_toes_y_acc', ...
            'left_arm_angle', 'right_arm_angle', 'left_leg_angle', 'right_leg_angle', ...
            'left_arm_phase', 'right_arm_phase', 'left_leg_phase', 'right_leg_phase' ...
          };
        data_labels_forceplate = ...
          {
            'left_fx', 'left_fy', 'left_fz', 'left_mx', 'left_my', 'left_mz', ...
            'right_fx', 'right_fy', 'right_fz', 'right_mx', 'right_my', 'right_mz' ...
            'left_foot_fx', 'left_foot_fy', 'left_foot_fz', 'left_foot_mx', 'left_foot_my', 'left_foot_mz', ...
            'right_foot_fx', 'right_foot_fy', 'right_foot_fz', 'right_foot_mx', 'right_foot_my', 'right_foot_mz' ...
          }
        data_labels_labview = ...
          {
            'stimulus_state' ...
          }
    end
    methods
        function this = WalkingTrialData(dataDirectory, condition, trialNumber)
            this.data_directory = dataDirectory;
            this.condition = condition;
            this.trial_number = trialNumber;
            this.loadSubjectInfo();
            this.loadMarkerTrajectories();
            this.loadForceplateTrajectories();
%             this.loadLabviewTrajectories();

            this.selected_time = this.time_marker(1);
        end
        function loadSubjectInfo(this)
            loaded_subject_info = load([this.data_directory filesep 'subjectInfo.mat']);
            this.date = loaded_subject_info.date;
            this.subject_id = loaded_subject_info.subject_id;
            this.subject_settings = loadSettingsFile('subjectSettings.txt');
            
            if exist([this.data_directory filesep 'subjectModel.mat'], 'file')
                loaded_subject_model = load([this.data_directory filesep 'subjectModel.mat']);
                this.joint_center_labels = loaded_subject_model.joint_center_headers;
            end
        end
        function loadMarkerTrajectories(this)
            [this.marker_positions, this.time_marker, this.sampling_rate_marker, this.marker_labels] = loadData(this.date, this.subject_id, this.condition, this.trial_number, 'marker_trajectories');
            
            % load kinematic data
            [joint_center_trajectories, ~, ~, this.joint_center_labels, this.kinematic_data_available] = loadData(this.date, this.subject_id, this.condition, this.trial_number, 'joint_center_trajectories', 'optional');
            this.joint_center_positions = joint_center_trajectories;
            [this.com_positions, ~, ~, this.com_labels, this.com_data_available] = loadData(this.date, this.subject_id, this.condition, this.trial_number, 'com_trajectories', 'optional');
            [this.joint_angles, ~, ~, ~, this.joint_angle_data_available] = loadData(this.date, this.subject_id, this.condition, this.trial_number, 'joint_angle_trajectories', 'optional');
            
            
            % time
            this.recording_time = this.time_marker(end);
            
            % extract data
            left_heel_marker = find(strcmp(this.marker_labels, 'LHEE'));
            left_toes_marker = find(strcmp(this.marker_labels, 'LTOE'));
            right_heel_marker = find(strcmp(this.marker_labels, 'RHEE'));
            right_toes_marker = find(strcmp(this.marker_labels, 'RTOE'));
            left_heel_marker_indices = reshape([(left_heel_marker - 1) * 3 + 1; (left_heel_marker - 1) * 3 + 2; (left_heel_marker - 1) * 3 + 3], 1, length(left_heel_marker)*3);
            left_toes_marker_indices = reshape([(left_toes_marker - 1) * 3 + 1; (left_toes_marker - 1) * 3 + 2; (left_toes_marker - 1) * 3 + 3], 1, length(left_toes_marker)*3);
            right_heel_marker_indices = reshape([(right_heel_marker - 1) * 3 + 1; (right_heel_marker - 1) * 3 + 2; (right_heel_marker - 1) * 3 + 3], 1, length(right_heel_marker)*3);
            right_toes_marker_indices = reshape([(right_toes_marker - 1) * 3 + 1; (right_toes_marker - 1) * 3 + 2; (right_toes_marker - 1) * 3 + 3], 1, length(right_toes_marker)*3);
            left_heel_z_trajectory = this.marker_positions(:, left_heel_marker_indices(3));
            left_toes_z_trajectory = this.marker_positions(:, left_toes_marker_indices(3));
            right_heel_z_trajectory = this.marker_positions(:, right_heel_marker_indices(3));
            right_toes_z_trajectory = this.marker_positions(:, right_toes_marker_indices(3));
            left_heel_y_trajectory = this.marker_positions(:, left_heel_marker_indices(2));
            left_toes_y_trajectory = this.marker_positions(:, left_toes_marker_indices(2));
            right_heel_y_trajectory = this.marker_positions(:, right_heel_marker_indices(2));
            right_toes_y_trajectory = this.marker_positions(:, right_toes_marker_indices(2));

            % calculate derivatives
            filter_order = 2;
            cutoff_frequency = 20; % cutoff frequency, in Hz
            [b, a] = butter(filter_order, cutoff_frequency/(this.sampling_rate_marker/2));	% set filter parameters for butterworth filter: 2=order of filter;
            left_heel_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, left_heel_z_trajectory), 1/this.sampling_rate_marker);
            right_heel_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, right_heel_z_trajectory), 1/this.sampling_rate_marker);
            left_heel_z_acc_trajectory = deriveByTime(nanfiltfilt(b, a, left_heel_z_vel_trajectory), 1/this.sampling_rate_marker);
            right_heel_z_acc_trajectory = deriveByTime(nanfiltfilt(b, a, right_heel_z_vel_trajectory), 1/this.sampling_rate_marker);
            left_toes_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, left_toes_z_trajectory), 1/this.sampling_rate_marker);
            right_toes_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, right_toes_z_trajectory), 1/this.sampling_rate_marker);
            left_toes_z_acc_trajectory = deriveByTime(nanfiltfilt(b, a, left_toes_z_vel_trajectory), 1/this.sampling_rate_marker);
            right_toes_z_acc_trajectory = deriveByTime(nanfiltfilt(b, a, right_toes_z_vel_trajectory), 1/this.sampling_rate_marker);        
            left_heel_y_vel_trajectory = deriveByTime(nanfiltfilt(b, a, left_heel_y_trajectory), 1/this.sampling_rate_marker);
            right_heel_y_vel_trajectory = deriveByTime(nanfiltfilt(b, a, right_heel_y_trajectory), 1/this.sampling_rate_marker);
            left_heel_y_acc_trajectory = deriveByTime(nanfiltfilt(b, a, left_heel_y_vel_trajectory), 1/this.sampling_rate_marker);
            right_heel_y_acc_trajectory = deriveByTime(nanfiltfilt(b, a, right_heel_y_vel_trajectory), 1/this.sampling_rate_marker);
            left_toes_y_vel_trajectory = deriveByTime(nanfiltfilt(b, a, left_toes_y_trajectory), 1/this.sampling_rate_marker);
            right_toes_y_vel_trajectory = deriveByTime(nanfiltfilt(b, a, right_toes_y_trajectory), 1/this.sampling_rate_marker);
            left_toes_y_acc_trajectory = deriveByTime(nanfiltfilt(b, a, left_toes_y_vel_trajectory), 1/this.sampling_rate_marker);
            right_toes_y_acc_trajectory = deriveByTime(nanfiltfilt(b, a, right_toes_y_vel_trajectory), 1/this.sampling_rate_marker);        

            this.left_heel_z_pos = left_heel_z_trajectory;
            this.left_heel_z_vel = left_heel_z_vel_trajectory;
            this.left_heel_z_acc = left_heel_z_acc_trajectory;

            this.right_heel_z_pos = right_heel_z_trajectory;
            this.right_heel_z_vel = right_heel_z_vel_trajectory;
            this.right_heel_z_acc = right_heel_z_acc_trajectory;

            this.left_toes_z_pos = left_toes_z_trajectory;
            this.left_toes_z_vel = left_toes_z_vel_trajectory;
            this.left_toes_z_acc = left_toes_z_acc_trajectory;

            this.right_toes_z_pos = right_toes_z_trajectory;
            this.right_toes_z_vel = right_toes_z_vel_trajectory;
            this.right_toes_z_acc = right_toes_z_acc_trajectory;        
            
            this.left_heel_y_pos = left_heel_y_trajectory;
            this.left_heel_y_vel = left_heel_y_vel_trajectory;
            this.left_heel_y_acc = left_heel_y_acc_trajectory;

            this.right_heel_y_pos = right_heel_y_trajectory;
            this.right_heel_y_vel = right_heel_y_vel_trajectory;
            this.right_heel_y_acc = right_heel_y_acc_trajectory;
            
            % angles
            this.left_arm_angle = loadData(this.date, this.subject_id, this.condition, this.trial_number, 'left_arm_angle', 'optional');
            this.right_arm_angle = loadData(this.date, this.subject_id, this.condition, this.trial_number, 'right_arm_angle', 'optional');
            this.left_leg_angle = loadData(this.date, this.subject_id, this.condition, this.trial_number, 'left_leg_angle', 'optional');
            this.right_leg_angle = loadData(this.date, this.subject_id, this.condition, this.trial_number, 'right_leg_angle', 'optional');
            this.left_arm_phase = loadData(this.date, this.subject_id, this.condition, this.trial_number, 'left_arm_phase', 'optional');
            this.right_arm_phase = loadData(this.date, this.subject_id, this.condition, this.trial_number, 'right_arm_phase', 'optional');
            this.left_leg_phase = loadData(this.date, this.subject_id, this.condition, this.trial_number, 'left_leg_phase', 'optional');
            this.right_leg_phase = loadData(this.date, this.subject_id, this.condition, this.trial_number, 'right_leg_phase', 'optional');
            
            % joint angles
            joint_angle_trajectories = loadData(this.date, this.subject_id, this.condition, this.trial_number, 'joint_angle_trajectories', 'optional');
        end
        function loadForceplateTrajectories(this)
            file_name_forceplate = [this.data_directory filesep 'processed' filesep makeFileName(this.date, this.subject_id, this.condition, this.trial_number, 'forceplateTrajectories.mat')];
            if exist(file_name_forceplate, 'file')
                loaded_forceplate_trajectories = load(file_name_forceplate);
                this.time_forceplate = loaded_forceplate_trajectories.time_forceplate;

                this.left_fx = loaded_forceplate_trajectories.fxl_trajectory;
                this.left_fy = loaded_forceplate_trajectories.fyl_trajectory;
                this.left_fz = loaded_forceplate_trajectories.fzl_trajectory;
                this.left_mx = loaded_forceplate_trajectories.mxl_trajectory;
                this.left_my = loaded_forceplate_trajectories.myl_trajectory;
                this.left_mz = loaded_forceplate_trajectories.mzl_trajectory;
                this.right_fx = loaded_forceplate_trajectories.fxr_trajectory;
                this.right_fy = loaded_forceplate_trajectories.fyr_trajectory;
                this.right_fz = loaded_forceplate_trajectories.fzr_trajectory;
                this.right_mx = loaded_forceplate_trajectories.mxr_trajectory;
                this.right_my = loaded_forceplate_trajectories.myr_trajectory;
                this.right_mz = loaded_forceplate_trajectories.mzr_trajectory;
                
                this.left_foot_fx = loaded_forceplate_trajectories.left_foot_wrench_world(:, 1);
                this.left_foot_fy = loaded_forceplate_trajectories.left_foot_wrench_world(:, 2);
                this.left_foot_fz = loaded_forceplate_trajectories.left_foot_wrench_world(:, 3);
                this.left_foot_mx = loaded_forceplate_trajectories.left_foot_wrench_world(:, 4);
                this.left_foot_my = loaded_forceplate_trajectories.left_foot_wrench_world(:, 5);
                this.left_foot_mz = loaded_forceplate_trajectories.left_foot_wrench_world(:, 6);
                this.right_foot_fx = loaded_forceplate_trajectories.right_foot_wrench_world(:, 1);
                this.right_foot_fy = loaded_forceplate_trajectories.right_foot_wrench_world(:, 2);
                this.right_foot_fz = loaded_forceplate_trajectories.right_foot_wrench_world(:, 3);
                this.right_foot_mx = loaded_forceplate_trajectories.right_foot_wrench_world(:, 4);
                this.right_foot_my = loaded_forceplate_trajectories.right_foot_wrench_world(:, 5);
                this.right_foot_mz = loaded_forceplate_trajectories.right_foot_wrench_world(:, 6);
                
            else
                this.time_forceplate = [];

                this.left_fx = [];
                this.left_fy = [];
                this.left_fz = [];
                this.left_mx = [];
                this.left_my = [];
                this.left_mz = [];
                this.right_fx = [];
                this.right_fy = [];
                this.right_fz = [];
                this.right_mx = [];
                this.right_my = [];
                this.right_mz = [];
                
                this.left_foot_fx = [];
                this.left_foot_fy = [];
                this.left_foot_fz = [];
                this.left_foot_mx = [];
                this.left_foot_my = [];
                this.left_foot_mz = [];
                this.right_foot_fx = [];
                this.right_foot_fy = [];
                this.right_foot_fz = [];
                this.right_foot_mx = [];
                this.right_foot_my = [];
                this.right_foot_mz = [];
            end
            
            
        end
        
        function time = getTime(this, data_label)
            if any(strcmp(data_label, this.data_labels_mocap))
                time = this.time_marker;
            elseif any(strcmp(data_label, this.data_labels_forceplate))
                time = this.time_forceplate;
            elseif any(strcmp(data_label, this.data_labels_labview))
                time = this.time_labview;
            else
                error('provided unknown data label');
            end
        end
        function data = getData(this, data_label)
            eval(['data = this.' data_label ';']);
        end
        
        function stepSelectedTime(this, direction, stepsize)
            if nargin < 3
                stepsize = 1;
            end
            
            % get current time step
            [~, time_index_mocap] = min(abs(this.time_marker - this.selected_time));
            
            if strcmp(direction, 'back')
                new_time_index_mocap = time_index_mocap - stepsize;
            elseif strcmp(direction, 'forward')
                new_time_index_mocap = time_index_mocap + stepsize;
            else
                error('Direction must be either "back" or "forward"');
            end
            
            % enforce limits
            if new_time_index_mocap < 1
                new_time_index_mocap = 1;
            elseif new_time_index_mocap > length(this.time_marker)
                new_time_index_mocap = length(this.time_marker);
            end
            
            % set result
            this.selected_time = this.time_marker(new_time_index_mocap);
        end
    end
end
            
            