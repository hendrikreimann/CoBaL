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

% this function finds the heelstrike and pushoff events

% input: 
% subjectInfo.mat
% subjectModel.mat
% markerTrajectories
% left_touchdown_method / right_touchdown_method in subjectSettings
%
% output:
% file stepEvents.mat, containing
% - left_pushoff_times
% - left_touchdown_times
% - right_pushoff_times
% - right_touchdown_times


function findStepEvents(varargin)
    % parse arguments
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;

    % figure out folders
    if ~exist('analysis', 'dir')
        mkdir('analysis')
    end
    load('subjectInfo.mat', 'date', 'subject_id');

    % load settings
    subject_settings = loadSettingsFile('subjectSettings.txt');
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = loadSettingsFile(study_settings_file);
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            %% prepare
            % load data
            condition = condition_list{i_condition};
            [marker_trajectories, time_marker, sampling_rate_marker, marker_labels] = loadData(date, subject_id, condition, i_trial, 'marker_trajectories');
            [left_foot_wrench_world, time_left_forceplate, ~, ~, left_forceplate_available] = loadData(date, subject_id, condition, i_trial, 'left_foot_wrench_world', 'optional');
            [right_foot_wrench_world, time_right_forceplate, ~, ~, right_forceplate_available] = loadData(date, subject_id, condition, i_trial, 'right_foot_wrench_world', 'optional');
            if left_forceplate_available && right_forceplate_available
                left_fz_trajectory = left_foot_wrench_world(:, 3);
                right_fz_trajectory = right_foot_wrench_world(:, 3);
            end
            
            % extract data
            LHEE_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LHEE');
            LHEE_z_trajectory = LHEE_trajectory(:, 3);
            LTOE_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LTOE');
            LTOE_z_trajectory = LTOE_trajectory(:, 3);
            
            RHEE_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RHEE');
            RHEE_z_trajectory = RHEE_trajectory(:, 3);
            RTOE_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RTOE');
            RTOE_z_trajectory = RTOE_trajectory(:, 3);
            
            LELB_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LELB');
            LWRA_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LWRA');
            LWRB_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LWRB');
            RELB_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RELB');
            RWRA_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RWRA');
            RWRB_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RWRB');
            
            LANK_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LANK');
            LPSI_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LPSI');
            LASI_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'LASI');
            RANK_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RANK');
            RPSI_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RPSI');
            RASI_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, 'RASI');
            
            % calculate derivatives
            filter_order = 2;
            cutoff_frequency = 20; % cutoff frequency, in Hz
            [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate_marker/2));	% set filter parameters for butterworth filter: 2=order of filter;
            LHEE_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, LHEE_z_trajectory), 1/sampling_rate_marker);
            RHEE_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, RHEE_z_trajectory), 1/sampling_rate_marker);
            LHEE_z_acc_trajectory = deriveByTime(nanfiltfilt(b, a, LHEE_z_vel_trajectory), 1/sampling_rate_marker);
            RHEE_z_acc_trajectory = deriveByTime(nanfiltfilt(b, a, RHEE_z_vel_trajectory), 1/sampling_rate_marker);
            LTOE_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, LTOE_z_trajectory), 1/sampling_rate_marker);
            RTOE_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, RTOE_z_trajectory), 1/sampling_rate_marker);
            
            % struct for saving
            variables_to_save = struct;

            %% find events for left foot
            left_touchdown_times = [];
            if any(strcmp(subject_settings.left_touchdown_method, 'heel_position_minima'))
                % find touch down indices as negative peaks of the heel marker z-position
                [~, left_heel_peak_locations] = findpeaks(-LHEE_z_trajectory, 'MinPeakProminence', subject_settings.left_touchdown_peak_prominence_threshold, 'MinPeakDistance', subject_settings.left_touchdown_peak_distance_threshold * sampling_rate_marker);
                left_touchdown_indices_mocap = left_heel_peak_locations';
                left_touchdown_times = [left_touchdown_times; time_marker(left_touchdown_indices_mocap)];
            end
            if any(strcmp(subject_settings.left_touchdown_method, 'toe_position_minima'))
                % find touch down indices as negative peaks of the heel marker z-position
                [~, left_toe_peak_locations] = findpeaks(-LTOE_z_trajectory, 'MinPeakProminence', subject_settings.left_touchdown_peak_prominence_threshold, 'MinPeakDistance', subject_settings.left_touchdown_peak_distance_threshold * sampling_rate_marker);
                left_touchdown_indices_mocap = left_toe_peak_locations';
                left_touchdown_times = [left_touchdown_times; time_marker(left_touchdown_indices_mocap)];
            end
            if any(strcmp(subject_settings.left_touchdown_method, 'toe_velocity_minima'))
                [~, left_toe_peak_locations] = findpeaks(-LTOE_z_vel_trajectory, 'MinPeakProminence', subject_settings.left_touchdown_peak_prominence_threshold, 'MinPeakDistance', subject_settings.left_touchdown_peak_distance_threshold * sampling_rate_marker);
                left_touchdown_indices_mocap = left_toe_peak_locations';
                left_touchdown_times = [left_touchdown_times; time_marker(left_touchdown_indices_mocap)];
            end
            if any(strcmp(subject_settings.left_touchdown_method, 'first_acceleration_peak_after_mid_swing'))
                % left_touchdown_peak_distance_threshold: 0.8
                % left touchdown peak prominence threshold: 5
                
                % find mid-swing as peaks of heel position
                [~, left_heel_midswing_locations] = findpeaks(LHEE_z_trajectory, 'MinPeakDistance', subject_settings.left_touchdown_peak_distance_threshold * sampling_rate_marker);
                
                % find acceleration peaks
                [~, left_heel_acc_peak_locations] = findpeaks(LHEE_z_acc_trajectory, 'MinPeakProminence', subject_settings.left_touchdown_peak_prominence_threshold);

                % identify acceleration peaks as touchdowns
                left_touchdown_indices_mocap = zeros(size(left_heel_midswing_locations));
                for i_step = 1 : length(left_heel_midswing_locations)
                    touchdown_index_index = left_heel_acc_peak_locations(find(left_heel_acc_peak_locations > left_heel_midswing_locations(i_step), 1, 'first'));
                    if ~isempty(touchdown_index_index)
                        left_touchdown_indices_mocap(i_step) = touchdown_index_index;
                    end
                end
                left_touchdown_indices_mocap(left_touchdown_indices_mocap==0) = [];
                left_touchdown_times = [left_touchdown_times; time_marker(left_touchdown_indices_mocap)];
            end

            left_pushoff_times = [];
            if any(strcmp(subject_settings.left_pushoff_method, 'first_velocity_peak_after_touchdown'))
                % for pushoff, find the first significant toes z-velocity peak after each touchdown
                [~, left_toes_vel_peak_locations] = findpeaks(LTOE_z_vel_trajectory, 'MinPeakProminence', subject_settings.left_pushoff_peak_prominence_threshold);
                left_toes_vel_peak_locations = left_toes_vel_peak_locations';
                left_pushoff_indices_mocap = zeros(size(left_touchdown_indices_mocap));
                for i_touchdown = 1 : length(left_touchdown_indices_mocap)
                    pushoff_index_index = find(left_toes_vel_peak_locations > left_touchdown_indices_mocap(i_touchdown), 1, 'first');
                    if ~isempty(pushoff_index_index)
                        left_pushoff_indices_mocap(i_touchdown) = left_toes_vel_peak_locations(pushoff_index_index);
                    end
                end
                left_pushoff_indices_mocap(left_pushoff_indices_mocap==0) = [];
                left_pushoff_times = [left_pushoff_times; time_marker(left_pushoff_indices_mocap)];
            end
            if any(strcmp(subject_settings.left_pushoff_method, 'forceplate_threshold'))
                left_pushoff_diff_forceplate = diff(sign(abs(left_fz_trajectory) - subject_settings.forceplate_load_threshold));
                left_pushoff_times = [left_pushoff_times; time_left_forceplate(left_pushoff_diff_forceplate~=0)];
            end
            
            % add new variables to be saved
            variables_to_save.left_pushoff_times = left_pushoff_times;
            variables_to_save.left_touchdown_times = left_touchdown_times;

            %% find events for right foot
            right_touchdown_times = [];
            if any(strcmp(subject_settings.right_touchdown_method, 'heel_position_minima'))
                % find touch down indices as negative peaks of the heel marker z-position
                [~, right_heel_peak_locations] = findpeaks(-RHEE_z_trajectory, 'MinPeakProminence', subject_settings.right_touchdown_peak_prominence_threshold, 'MinPeakDistance', subject_settings.right_touchdown_peak_distance_threshold);
                right_touchdown_indices_mocap = right_heel_peak_locations';
                right_touchdown_times = [right_touchdown_times; time_marker(right_touchdown_indices_mocap)];
            end
            if any(strcmp(subject_settings.right_touchdown_method, 'toe_position_minima'))
                % find touch down indices as negative peaks of the heel marker z-position
                [~, right_toe_peak_locations] = findpeaks(-RTOE_z_trajectory, 'MinPeakProminence', subject_settings.right_touchdown_peak_prominence_threshold, 'MinPeakDistance', subject_settings.right_touchdown_peak_distance_threshold);
                right_touchdown_indices_mocap = right_toe_peak_locations';
                right_touchdown_times = [right_touchdown_times; time_marker(right_touchdown_indices_mocap)];
            end
            if any(strcmp(subject_settings.right_touchdown_method, 'toe_velocity_minima'))
                [~, right_toe_peak_locations] = findpeaks(-RTOE_z_vel_trajectory, 'MinPeakProminence', subject_settings.right_touchdown_peak_prominence_threshold, 'MinPeakDistance', subject_settings.right_touchdown_peak_distance_threshold * sampling_rate_marker);
                right_touchdown_indices_mocap = right_toe_peak_locations';
                right_touchdown_times = [right_touchdown_times; time_marker(right_touchdown_indices_mocap)];
            end
            if any(strcmp(subject_settings.right_touchdown_method, 'first_acceleration_peak_after_mid_swing'))
                % right_touchdown_peak_distance_threshold: 0.8
                % right touchdown peak prominence threshold: 5
                
                % find mid-swing as peaks of heel position
                [~, right_heel_midswing_locations] = findpeaks(RHEE_z_trajectory, 'MinPeakDistance', subject_settings.right_touchdown_peak_distance_threshold * sampling_rate_marker);
                
                % find acceleration peaks
                [~, right_heel_acc_peak_locations] = findpeaks(RHEE_z_acc_trajectory, 'MinPeakProminence', subject_settings.right_touchdown_peak_prominence_threshold);

                % identify acceleration peaks as touchdowns
                right_touchdown_indices_mocap = zeros(size(right_heel_midswing_locations));
                for i_step = 1 : length(right_heel_midswing_locations)
                    touchdown_index_index = right_heel_acc_peak_locations(find(right_heel_acc_peak_locations > right_heel_midswing_locations(i_step), 1, 'first'));
                    if ~isempty(touchdown_index_index)
                        right_touchdown_indices_mocap(i_step) = touchdown_index_index;
                    end
                end
                right_touchdown_indices_mocap(right_touchdown_indices_mocap==0) = [];
                right_touchdown_times = [right_touchdown_times; time_marker(right_touchdown_indices_mocap)];
            end
  
            right_pushoff_times = []; 
            if any(strcmp(subject_settings.right_pushoff_method, 'first_velocity_peak_after_touchdown'))
                % for pushoff, find the first significant toes z-velocity peak after each touchdown
                [~, right_toes_vel_peak_locations] = findpeaks(RTOE_z_vel_trajectory, 'MinPeakProminence', subject_settings.right_pushoff_peak_prominence_threshold);
                right_toes_vel_peak_locations = right_toes_vel_peak_locations';
                right_pushoff_indices_mocap = zeros(size(right_touchdown_indices_mocap));
                for i_touchdown = 1 : length(right_touchdown_indices_mocap)
                    pushoff_index_index = find(right_toes_vel_peak_locations > right_touchdown_indices_mocap(i_touchdown), 1, 'first');
                    if ~isempty(pushoff_index_index)
                        right_pushoff_indices_mocap(i_touchdown) = right_toes_vel_peak_locations(pushoff_index_index);
                    end
                end
                right_pushoff_indices_mocap(right_pushoff_indices_mocap==0) = [];
                right_pushoff_times = [right_pushoff_times; time_marker(right_pushoff_indices_mocap)];
            end
            if any(strcmp(subject_settings.right_pushoff_method, 'forceplate_threshold'))
                right_pushoff_diff_forceplate = diff(sign(abs(right_fz_trajectory) - subject_settings.forceplate_load_threshold));
                right_pushoff_times = [right_pushoff_times; time_right_forceplate(right_pushoff_diff_forceplate~=0)];
            end

            % remove duplicates
            left_pushoff_times = unique(left_pushoff_times);
            left_touchdown_times = unique(left_touchdown_times);
            right_pushoff_times = unique(right_pushoff_times);
            right_touchdown_times = unique(right_touchdown_times);

            % determine indices
            left_touchdown_indices_mocap = zeros(size(left_touchdown_times));
            for i_index = 1 : length(left_touchdown_times)
                [~, index_mocap] = min(abs(time_marker - left_touchdown_times(i_index)));
                left_touchdown_indices_mocap(i_index) = index_mocap;
            end
            left_pushoff_indices_mocap = zeros(size(left_pushoff_times));
            for i_index = 1 : length(left_pushoff_times)
                [~, index_mocap] = min(abs(time_marker - left_pushoff_times(i_index)));
                left_pushoff_indices_mocap(i_index) = index_mocap;
            end
            right_touchdown_indices_mocap = zeros(size(right_touchdown_times));
            for i_index = 1 : length(right_touchdown_times)
                [~, index_mocap] = min(abs(time_marker - right_touchdown_times(i_index)));
                right_touchdown_indices_mocap(i_index) = index_mocap;
            end
            right_pushoff_indices_mocap = zeros(size(right_pushoff_times));
            for i_index = 1 : length(right_pushoff_times)
                [~, index_mocap] = min(abs(time_marker - right_pushoff_times(i_index)));
                right_pushoff_indices_mocap(i_index) = index_mocap;
            end
            
            % add new variables to be saved
            variables_to_save.right_pushoff_times = right_pushoff_times;
            variables_to_save.right_touchdown_times = right_touchdown_times;

            %% find events for angles
            % TODO: change conditionals to use a WalkingDataCustodian
            if any(strcmp(study_settings.variables_to_analyze(:, 1), 'left_arm_phase')) || any(strcmp(study_settings.variables_to_analyze(:, 1), 'left_arm_right_leg_relative_phase'))
                % calculate vectors
                left_wrist_center_trajectory = (LWRA_trajectory + LWRB_trajectory) * 0.5;
                left_arm_vector_trajectory = LELB_trajectory - left_wrist_center_trajectory;

                % calculate angles
                larm_angle = rad2deg(atan2(-left_arm_vector_trajectory(:, 2), left_arm_vector_trajectory(:, 3)));

                % find negative peaks
                [~, left_arm_swing_onset_indices] = findpeaks(-larm_angle, 'MinPeakProminence', subject_settings.left_armswing_peak_prominence_threshold, 'MinPeakDistance', subject_settings.left_armswing_peak_distance_threshold * sampling_rate_marker);
                left_arm_swing_onset_times = time_marker(left_arm_swing_onset_indices);

                % add new variables to be saved
                variables_to_save.left_arm_swing_onset_times = left_arm_swing_onset_times;
            end
            if any(strcmp(study_settings.variables_to_analyze(:, 1), 'right_arm_phase')) || any(strcmp(study_settings.variables_to_analyze(:, 1), 'right_arm_left_leg_relative_phase'))
                % calculate vectors
                right_wrist_center_trajectory = (RWRA_trajectory + RWRB_trajectory) * 0.5;
                right_arm_vector_trajectory = RELB_trajectory - right_wrist_center_trajectory;
                
                % calculate angles
                rarm_angle = rad2deg(atan2(-right_arm_vector_trajectory(:, 2), right_arm_vector_trajectory(:, 3)));

                % find negative peaks
                [~, right_arm_swing_onset_indices] = findpeaks(-rarm_angle, 'MinPeakProminence', subject_settings.right_armswing_peak_prominence_threshold, 'MinPeakDistance', subject_settings.right_armswing_peak_distance_threshold * sampling_rate_marker);
                right_arm_swing_onset_times = time_marker(right_arm_swing_onset_indices);

                % add new variables to be saved
                variables_to_save.right_arm_swing_onset_times = right_arm_swing_onset_times;
            end
            if any(strcmp(study_settings.variables_to_analyze(:, 1), 'left_leg_phase')) || any(strcmp(study_settings.variables_to_analyze(:, 1), 'left_arm_right_leg_relative_phase'))
                % calculate vectors
                left_pelvis_center_trajectory = (LPSI_trajectory + LASI_trajectory) * 0.5;
                left_leg_vector_trajectory = left_pelvis_center_trajectory - LANK_trajectory;
                
                % calculate angles
                lleg_angle = rad2deg(atan2(-left_leg_vector_trajectory(:, 2), left_leg_vector_trajectory(:, 3)));

                % find negative peaks
                [~, left_leg_swing_onset_indices] = findpeaks(-lleg_angle, 'MinPeakProminence', subject_settings.left_legswing_peak_prominence_threshold, 'MinPeakDistance', subject_settings.left_legswing_peak_distance_threshold * sampling_rate_marker);
                left_leg_swing_onset_times = time_marker(left_leg_swing_onset_indices);

                % add new variables to be saved
                variables_to_save.left_leg_swing_onset_times = left_leg_swing_onset_times;
            end
            if any(strcmp(study_settings.variables_to_analyze(:, 1), 'right_leg_phase')) || any(strcmp(study_settings.variables_to_analyze(:, 1), 'right_arm_left_leg_relative_phase'))
                % calculate vectors
                right_pelvis_center_trajectory = (RPSI_trajectory + RASI_trajectory) * 0.5;
                right_leg_vector_trajectory = right_pelvis_center_trajectory - RANK_trajectory;
                
                % calculate angles
                rleg_angle = rad2deg(atan2(-right_leg_vector_trajectory(:, 2), right_leg_vector_trajectory(:, 3)));

                % find negative peaks
                [~, right_leg_swing_onset_indices] = findpeaks(-rleg_angle, 'MinPeakProminence', subject_settings.right_legswing_peak_prominence_threshold, 'MinPeakDistance', subject_settings.right_legswing_peak_distance_threshold * sampling_rate_marker);
                right_leg_swing_onset_times = time_marker(right_leg_swing_onset_indices);

                % add new variables to be saved
                variables_to_save.right_leg_swing_onset_times = right_leg_swing_onset_times;
            end
            


            
            % normalize
%             larm_angle_normalized = normalizePeriodicVariable(larm_angle, left_arm_peak_locations);
%             rarm_angle_normalized = normalizePeriodicVariable(rarm_angle, right_arm_peak_locations);
%             lleg_angle_normalized = normalizePeriodicVariable(lleg_angle, left_leg_peak_locations);
%             rleg_angle_normalized = normalizePeriodicVariable(rleg_angle, right_leg_peak_locations);
            
            %% visualize
            color_heelstrike = [1 0 0];
            color_pushoff = [0 1 0];
            color_peak = [0.5, 0.2, 1];

            force_scaler = 2e-4;
            vel_scaler = 1;
            acc_scaler = .2;



            if left_forceplate_available && right_forceplate_available
                left_fz_trajectory_marker = spline(time_left_forceplate, left_fz_trajectory, time_marker);
                right_fz_trajectory_marker = spline(time_right_forceplate, right_fz_trajectory, time_marker);
            end
            if visualize
                step_event_figures = zeros(1, 4);
                
                % left position
                step_event_figures(1) = figure; axes_left = axes; hold on; title('left foot marker positions')
                plot(time_marker, LHEE_z_trajectory, 'linewidth', 1, 'displayname', 'left heel vertical');
                plot(time_marker, LTOE_z_trajectory, 'linewidth', 1, 'displayname', 'left toes vertical');
                plot(time_marker(left_touchdown_indices_mocap), LHEE_z_trajectory(left_touchdown_indices_mocap), 'v', 'linewidth', 2, 'color', color_heelstrike, 'displayname', 'left touchdown');
                plot(time_marker(left_pushoff_indices_mocap), LTOE_z_trajectory(left_pushoff_indices_mocap), '^', 'linewidth', 2, 'color', color_pushoff, 'displayname', 'left pushoff');
                plot(time_marker(left_leg_swing_onset_indices), LTOE_z_trajectory(left_leg_swing_onset_indices), '+', 'linewidth', 2, 'color', color_peak, 'displayname', 'left leg angle peaks');
                if left_forceplate_available
                    plot(time_left_forceplate, left_fz_trajectory*force_scaler, 'displayname', 'left forceplate')
                    h = plot(time_marker(left_touchdown_indices_mocap), left_fz_trajectory_marker(left_touchdown_indices_mocap)*force_scaler, 'v', 'linewidth', 2, 'color', color_heelstrike);
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    h = plot(time_marker(left_pushoff_indices_mocap), left_fz_trajectory_marker(left_pushoff_indices_mocap)*force_scaler, '^', 'linewidth', 2, 'color', color_pushoff);
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                legend('toggle');

                % left derivatives
                step_event_figures(2) = figure; axes_left_derivatives = axes; hold on;  title('left foot marker derivatives (scaled for qualitative overview)')
                plot(time_marker, LHEE_z_acc_trajectory*acc_scaler, 'linewidth', 1, 'displayname', 'left heel acceleration vertical');
                plot(time_marker, LTOE_z_vel_trajectory*vel_scaler, 'linewidth', 1, 'displayname', 'left toes velocity vertical');
                plot(time_marker(left_touchdown_indices_mocap), LHEE_z_acc_trajectory(left_touchdown_indices_mocap)*acc_scaler, 'v', 'linewidth', 2, 'color', color_heelstrike, 'displayname', 'left touchdown');
                plot(time_marker(left_pushoff_indices_mocap), LTOE_z_vel_trajectory(left_pushoff_indices_mocap)*vel_scaler, '^', 'linewidth', 2, 'color', color_pushoff, 'displayname', 'left pushoff');
                if left_forceplate_available
                    plot(time_left_forceplate, left_fz_trajectory*force_scaler, 'displayname', 'left forceplate')
                    h = plot(time_marker(left_touchdown_indices_mocap), left_fz_trajectory_marker(left_touchdown_indices_mocap)*force_scaler, 'v', 'linewidth', 2, 'color', color_heelstrike);
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    h = plot(time_marker(left_pushoff_indices_mocap), left_fz_trajectory_marker(left_pushoff_indices_mocap)*force_scaler, '^', 'linewidth', 2, 'color', color_pushoff);
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                legend('toggle');

                % right position
                step_event_figures(3) = figure; axes_right = axes; hold on; title('right foot marker positions')
                plot(time_marker, RHEE_z_trajectory, 'linewidth', 1, 'displayname', 'right heel vertical');
                plot(time_marker, RTOE_z_trajectory, 'linewidth', 1, 'displayname', 'right toes vertical');
                plot(time_marker(right_touchdown_indices_mocap), RHEE_z_trajectory(right_touchdown_indices_mocap), 'v', 'linewidth', 2, 'color', color_heelstrike, 'displayname', 'right touchdown');
                plot(time_marker(right_pushoff_indices_mocap), RTOE_z_trajectory(right_pushoff_indices_mocap), '^', 'linewidth', 2, 'color', color_pushoff, 'displayname', 'right pushoff');
                if right_forceplate_available
                    plot(time_right_forceplate, right_fz_trajectory*force_scaler, 'displayname', 'right forceplate')
                    h = plot(time_marker(right_touchdown_indices_mocap), right_fz_trajectory_marker(right_touchdown_indices_mocap)*force_scaler, 'v', 'linewidth', 2, 'color', color_heelstrike);
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    h = plot(time_marker(right_pushoff_indices_mocap), right_fz_trajectory_marker(right_pushoff_indices_mocap)*force_scaler, '^', 'linewidth', 2, 'color', color_pushoff);
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                legend('toggle');

                % right derivatives
                step_event_figures(4) = figure; axes_right_derivatives = axes; hold on;  title('right foot marker derivatives (scaled for qualitative overview)')
                plot(time_marker, RHEE_z_acc_trajectory*acc_scaler, 'linewidth', 1, 'displayname', 'right heel acceleration vertical');
                plot(time_marker, RTOE_z_vel_trajectory*vel_scaler, 'linewidth', 1, 'displayname', 'right toes velocity vertical');
                plot(time_marker(right_touchdown_indices_mocap), RHEE_z_acc_trajectory(right_touchdown_indices_mocap)*acc_scaler, 'v', 'linewidth', 2, 'color', color_heelstrike, 'displayname', 'right touchdown');
                plot(time_marker(right_pushoff_indices_mocap), RTOE_z_vel_trajectory(right_pushoff_indices_mocap)*vel_scaler, '^', 'linewidth', 2, 'color', color_pushoff, 'displayname', 'right pushoff');
                if right_forceplate_available
                    plot(time_right_forceplate, right_fz_trajectory*force_scaler, 'displayname', 'right forceplate')
                    h = plot(time_marker(right_touchdown_indices_mocap), right_fz_trajectory_marker(right_touchdown_indices_mocap)*force_scaler, 'v', 'linewidth', 2, 'color', color_heelstrike);
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    h = plot(time_marker(right_pushoff_indices_mocap), right_fz_trajectory_marker(right_pushoff_indices_mocap)*force_scaler, '^', 'linewidth', 2, 'color', color_pushoff);
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                legend('toggle');

                linkaxes([axes_left axes_left_derivatives axes_right axes_right_derivatives], 'x')
                
                distFig(step_event_figures, 'rows', 2)
            end

            % change event variables to column vectors if necessary
            if isrow(left_pushoff_times)
                left_pushoff_times = left_pushoff_times'; %#ok<NASGU>
            end
            if isrow(left_touchdown_times)
                left_touchdown_times = left_touchdown_times'; %#ok<NASGU>
            end
            if isrow(right_pushoff_times)
                right_pushoff_times = right_pushoff_times'; %#ok<NASGU>
            end
            if isrow(right_touchdown_times)
                right_touchdown_times = right_touchdown_times'; %#ok<NASGU>
            end

            %% save
            step_events_file_name = ['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'stepEvents')];
            save(step_events_file_name, '-struct', 'variables_to_save');

            disp(['Finding Step Events: condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' step_events_file_name]);
        end
    end
end















