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
        date = [];
        subject_id = [];
        sampling_rate_mocap = -1;
        sampling_rate_forceplate = -1;
        
        % time
        time_mocap = [];
        time_forceplate = [];
        time_labview = [];
        selected_time = [];
        
        % marker kinematics
        marker_positions = [];
        marker_headers = [];
        joint_center_positions = [];
        joint_center_headers = [];
        com_positions = [];
        com_headers = [];
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
        
        stimulus_state = [];
        
        data_labels_mocap = ...
          { ...
            'left_heel_z_pos', 'left_heel_z_vel', 'left_heel_z_acc', ...
            'left_toes_z_pos', 'left_toes_z_vel', 'left_toes_z_acc', ...
            'right_heel_z_pos', 'right_heel_z_vel', 'right_heel_z_acc', ...
            'right_toes_z_pos', 'right_toes_z_vel', 'right_toes_z_acc', ...
            'left_heel_y_pos', 'left_heel_y_vel', 'left_heel_y_acc', ...
            'left_toes_y_pos', 'left_toes_y_vel', 'left_toes_y_acc', ...
            'right_heel_y_pos', 'right_heel_y_vel', 'right_heel_y_acc', ...
            'right_toes_y_pos', 'right_toes_y_vel', 'right_toes_y_acc' ...
          };
        data_labels_forceplate = ...
          {
            'left_fx', 'left_fy', 'left_fz', 'left_mx', 'left_my', 'left_mz', ...
            'right_fx', 'right_fy', 'right_fz', 'right_mx', 'right_my', 'right_mz' ...
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

            this.selected_time = this.time_mocap(1);
        end
        function loadSubjectInfo(this)
            loaded_subject_info = load([this.data_directory filesep 'subjectInfo.mat']);
            this.date = loaded_subject_info.date;
            this.subject_id = loaded_subject_info.subject_id;
            
            if exist([this.data_directory filesep 'subjectModel.mat'], 'file')
                loaded_subject_model = load([this.data_directory filesep 'subjectModel.mat']);
                this.joint_center_headers = loaded_subject_model.joint_center_headers;
            end
        end
        function loadMarkerTrajectories(this)
            % load marker data
            loaded_marker_trajectories = load([this.data_directory filesep 'processed' filesep makeFileName(this.date, this.subject_id, this.condition, this.trial_number, 'markerTrajectories')]);
            this.marker_positions = loaded_marker_trajectories.marker_trajectories;
            this.marker_headers = loaded_marker_trajectories.marker_headers;
            
            % load kinematic data
            com_file_name = [this.data_directory filesep 'processed' filesep makeFileName(this.date, this.subject_id, this.condition, this.trial_number, 'kinematicTrajectories.mat')];
            if exist(com_file_name, 'file')
                loaded_trajectories = load(com_file_name);
                this.joint_center_positions = loaded_trajectories.joint_center_trajectories;
                this.com_positions = loaded_trajectories.com_trajectories;
                this.com_headers = loaded_trajectories.com_labels;
                this.joint_angles = loaded_trajectories.joint_angle_trajectories;
            end
            
            % time
            this.sampling_rate_mocap = loaded_marker_trajectories.sampling_rate_mocap;
            this.recording_time = loaded_marker_trajectories.time_mocap(end);
            time_mocap = loaded_marker_trajectories.time_mocap;
           
            % extract data
            left_heel_marker = find(strcmp(loaded_marker_trajectories.marker_headers, 'LHEE'));
            left_toes_marker = find(strcmp(loaded_marker_trajectories.marker_headers, 'LTOE'));
            right_heel_marker = find(strcmp(loaded_marker_trajectories.marker_headers, 'RHEE'));
            right_toes_marker = find(strcmp(loaded_marker_trajectories.marker_headers, 'RTOE'));
            left_heel_marker_indices = reshape([(left_heel_marker - 1) * 3 + 1; (left_heel_marker - 1) * 3 + 2; (left_heel_marker - 1) * 3 + 3], 1, length(left_heel_marker)*3);
            left_toes_marker_indices = reshape([(left_toes_marker - 1) * 3 + 1; (left_toes_marker - 1) * 3 + 2; (left_toes_marker - 1) * 3 + 3], 1, length(left_toes_marker)*3);
            right_heel_marker_indices = reshape([(right_heel_marker - 1) * 3 + 1; (right_heel_marker - 1) * 3 + 2; (right_heel_marker - 1) * 3 + 3], 1, length(right_heel_marker)*3);
            right_toes_marker_indices = reshape([(right_toes_marker - 1) * 3 + 1; (right_toes_marker - 1) * 3 + 2; (right_toes_marker - 1) * 3 + 3], 1, length(right_toes_marker)*3);
            left_heel_z_trajectory = loaded_marker_trajectories.marker_trajectories(:, left_heel_marker_indices(3));
            left_toes_z_trajectory = loaded_marker_trajectories.marker_trajectories(:, left_toes_marker_indices(3));
            right_heel_z_trajectory = loaded_marker_trajectories.marker_trajectories(:, right_heel_marker_indices(3));
            right_toes_z_trajectory = loaded_marker_trajectories.marker_trajectories(:, right_toes_marker_indices(3));
            left_heel_y_trajectory = loaded_marker_trajectories.marker_trajectories(:, left_heel_marker_indices(2));
            left_toes_y_trajectory = loaded_marker_trajectories.marker_trajectories(:, left_toes_marker_indices(2));
            right_heel_y_trajectory = loaded_marker_trajectories.marker_trajectories(:, right_heel_marker_indices(2));
            right_toes_y_trajectory = loaded_marker_trajectories.marker_trajectories(:, right_toes_marker_indices(2));

            % calculate derivatives
            filter_order = 2;
            cutoff_frequency = 20; % cutoff frequency, in Hz
            [b, a] = butter(filter_order, cutoff_frequency/(this.sampling_rate_mocap/2));	% set filter parameters for butterworth filter: 2=order of filter;
            left_heel_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, left_heel_z_trajectory), 1/this.sampling_rate_mocap);
            right_heel_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, right_heel_z_trajectory), 1/this.sampling_rate_mocap);
            left_heel_z_acc_trajectory = deriveByTime(nanfiltfilt(b, a, left_heel_z_vel_trajectory), 1/this.sampling_rate_mocap);
            right_heel_z_acc_trajectory = deriveByTime(nanfiltfilt(b, a, right_heel_z_vel_trajectory), 1/this.sampling_rate_mocap);
            left_toes_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, left_toes_z_trajectory), 1/this.sampling_rate_mocap);
            right_toes_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, right_toes_z_trajectory), 1/this.sampling_rate_mocap);
            left_toes_z_acc_trajectory = deriveByTime(nanfiltfilt(b, a, left_toes_z_vel_trajectory), 1/this.sampling_rate_mocap);
            right_toes_z_acc_trajectory = deriveByTime(nanfiltfilt(b, a, right_toes_z_vel_trajectory), 1/this.sampling_rate_mocap);        
            left_heel_y_vel_trajectory = deriveByTime(nanfiltfilt(b, a, left_heel_y_trajectory), 1/this.sampling_rate_mocap);
            right_heel_y_vel_trajectory = deriveByTime(nanfiltfilt(b, a, right_heel_y_trajectory), 1/this.sampling_rate_mocap);
            left_heel_y_acc_trajectory = deriveByTime(nanfiltfilt(b, a, left_heel_y_vel_trajectory), 1/this.sampling_rate_mocap);
            right_heel_y_acc_trajectory = deriveByTime(nanfiltfilt(b, a, right_heel_y_vel_trajectory), 1/this.sampling_rate_mocap);
            left_toes_y_vel_trajectory = deriveByTime(nanfiltfilt(b, a, left_toes_y_trajectory), 1/this.sampling_rate_mocap);
            right_toes_y_vel_trajectory = deriveByTime(nanfiltfilt(b, a, right_toes_y_trajectory), 1/this.sampling_rate_mocap);
            left_toes_y_acc_trajectory = deriveByTime(nanfiltfilt(b, a, left_toes_y_vel_trajectory), 1/this.sampling_rate_mocap);
            right_toes_y_acc_trajectory = deriveByTime(nanfiltfilt(b, a, right_toes_y_vel_trajectory), 1/this.sampling_rate_mocap);        

            % package
            this.time_mocap = loaded_marker_trajectories.time_mocap;

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
                
            end
            
            
        end
%         function loadLabviewTrajectories(this)
%             loaded_labview_trajectories = load([this.data_directory filesep makeFileName(this.date, this.subject_id, this.condition, this.trial_number, 'labviewTrajectories')]);
%             
%             this.time_labview = loaded_labview_trajectories.time_labview;
%             this.stimulus_state = loaded_labview_trajectories.stimulus_state_trajectory;
%         end
        
        function time = getTime(this, data_label)
            if any(strcmp(data_label, this.data_labels_mocap))
                time = this.time_mocap;
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
            [~, time_index_mocap] = min(abs(this.time_mocap - this.selected_time));
            
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
            elseif new_time_index_mocap > length(this.time_mocap)
                new_time_index_mocap = length(this.time_mocap);
            end
            
            % set result
            this.selected_time = this.time_mocap(new_time_index_mocap);
        end
    end
end
            
            