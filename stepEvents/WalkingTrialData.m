classdef WalkingTrialData < handle
    properties
        % administration
        data_directory = [];
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
        
        % marker kinematics
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
            'right_toes_z_pos', 'right_toes_z_vel', 'right_toes_z_acc' ...
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
        function this = WalkingTrialData(dataDirectory, trialNumber)
            if nargin > 1
                this.data_directory = dataDirectory;
                this.trial_number = trialNumber;
                this.loadSubjectInfo();
                this.loadMarkerTrajectories();
                this.loadForceplateTrajectories();
                this.loadLabviewTrajectories();
            end
        end
        function loadSubjectInfo(this)
            loaded_subject_info = load([this.data_directory filesep 'subjectInfo.mat']);
            this.recording_time = loaded_subject_info.recordingTime;
            this.date = loaded_subject_info.date;
            this.subject_id = loaded_subject_info.subject_id;
        end
        function loadMarkerTrajectories(this)
            loaded_marker_trajectories = load([this.data_directory filesep makeFileName(this.date, this.subject_id, 'walking', this.trial_number, 'markerTrajectories')]);
            
            this.sampling_rate_mocap = loaded_marker_trajectories.sampling_rate_mocap;

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

            % calculate derivatives
            filter_order = 2;
            cutoff_frequency = 20; % cutoff frequency, in Hz
            [b, a] = butter(filter_order, cutoff_frequency/(this.sampling_rate_mocap/2));	% set filter parameters for butterworth filter: 2=order of filter;
            left_heel_z_vel_trajectory = deriveByTime(filtfilt(b, a, left_heel_z_trajectory), 1/this.sampling_rate_mocap);
            right_heel_z_vel_trajectory = deriveByTime(filtfilt(b, a, right_heel_z_trajectory), 1/this.sampling_rate_mocap);
            left_heel_z_acc_trajectory = deriveByTime(filtfilt(b, a, left_heel_z_vel_trajectory), 1/this.sampling_rate_mocap);
            right_heel_z_acc_trajectory = deriveByTime(filtfilt(b, a, right_heel_z_vel_trajectory), 1/this.sampling_rate_mocap);
            left_toes_z_vel_trajectory = deriveByTime(filtfilt(b, a, left_toes_z_trajectory), 1/this.sampling_rate_mocap);
            right_toes_z_vel_trajectory = deriveByTime(filtfilt(b, a, right_toes_z_trajectory), 1/this.sampling_rate_mocap);
            left_toes_z_acc_trajectory = deriveByTime(filtfilt(b, a, left_toes_z_vel_trajectory), 1/this.sampling_rate_mocap);
            right_toes_z_acc_trajectory = deriveByTime(filtfilt(b, a, right_toes_z_vel_trajectory), 1/this.sampling_rate_mocap);        

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
        end
        function loadForceplateTrajectories(this)
            loaded_forceplate_trajectories = load([this.data_directory filesep makeFileName(this.date, this.subject_id, 'walking', this.trial_number, 'forceplateTrajectories')]);
            
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
            
        end
        function loadLabviewTrajectories(this)
            loaded_labview_trajectories = load([this.data_directory filesep makeFileName(this.date, this.subject_id, 'walking', this.trial_number, 'labviewTrajectories')]);
            
            this.time_labview = loaded_labview_trajectories.time_labview;
            this.stimulus_state = loaded_labview_trajectories.stimulus_foot_state;
        end
        
        function time = getTime(this, data_label)
            if any(strcmp(data_label, this.data_labels_mocap))
                time = this.time_mocap;
            elseif any(strcmp(data_label, this.data_labels_forceplate))
                time = this.time_forceplate;
            elseif any(strcmp(data_label, this.data_labels_labview))
                time = this.time_labview;
            end
        end
        function data = getData(this, data_label)
            eval(['data = this.' data_label ';']);
        end
        
    end
end
            
            