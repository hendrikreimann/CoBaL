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

% this function calculates the swing phase for arms and legs during walking

% input:
% subjectInfo.mat
% subjectModel.mat
% markerTrajectories
%
% output:
% file kinematicTrajectories.mat, containing
% - joint_center_trajectories
% - com_trajectories
% - com_labels
% - joint_angle_trajectories


function calculateSwingPhase(varargin)
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    
    % load settings
    load('subjectInfo.mat', 'date', 'subject_id');
    subject_settings = loadSettingsFile('subjectSettings.txt');
    
%     for i_condition = 1 : length(condition_list)
    for i_condition = 1
        trials_to_process = trial_number_list{i_condition};
        trials_to_process = 1;
        for i_trial = trials_to_process
            % load data
            condition = condition_list{i_condition};
            [marker_trajectories, time_marker, sampling_rate_marker, marker_labels] = loadData(date, subject_id, condition, i_trial, 'marker_trajectories');
            
            number_of_time_steps = size(marker_trajectories, 1);

            % extract necessary data
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
            
            % calculate vectors
            left_wrist_center_trajectory = (LWRA_trajectory + LWRB_trajectory) * 0.5;
            left_arm_vector_trajectory = LELB_trajectory - left_wrist_center_trajectory;
            left_pelvis_center_trajectory = (LPSI_trajectory + LASI_trajectory) * 0.5;
            left_leg_vector_trajectory = left_pelvis_center_trajectory - LANK_trajectory;
            right_wrist_center_trajectory = (RWRA_trajectory + RWRB_trajectory) * 0.5;
            right_arm_vector_trajectory = RELB_trajectory - right_wrist_center_trajectory;
            right_pelvis_center_trajectory = (RPSI_trajectory + RASI_trajectory) * 0.5;
            right_leg_vector_trajectory = right_pelvis_center_trajectory - RANK_trajectory;

            % calculate angles
            larm_angle = rad2deg(atan2(-left_arm_vector_trajectory(:, 2), left_arm_vector_trajectory(:, 3)));
            rarm_angle = rad2deg(atan2(-right_arm_vector_trajectory(:, 2), right_arm_vector_trajectory(:, 3)));
            lleg_angle = rad2deg(atan2(-left_leg_vector_trajectory(:, 2), left_leg_vector_trajectory(:, 3)));
            rleg_angle = rad2deg(atan2(-right_leg_vector_trajectory(:, 2), right_leg_vector_trajectory(:, 3)));
            
            % fake data
%             larm_angle = sin(time_marker*2*pi);
%             [~, left_arm_peak_locations] = findpeaks(-larm_angle, 'MinPeakProminence', 0.5, 'MinPeakDistance', subject_settings.left_armswing_peak_distance_threshold * sampling_rate_marker);
            
            % find negative peaks
            [~, left_arm_peak_locations] = findpeaks(-larm_angle, 'MinPeakProminence', subject_settings.left_armswing_peak_prominence_threshold, 'MinPeakDistance', subject_settings.left_armswing_peak_distance_threshold * sampling_rate_marker);
            [~, right_arm_peak_locations] = findpeaks(-rarm_angle, 'MinPeakProminence', subject_settings.right_armswing_peak_prominence_threshold, 'MinPeakDistance', subject_settings.right_armswing_peak_distance_threshold * sampling_rate_marker);
            [~, left_leg_peak_locations] = findpeaks(-lleg_angle, 'MinPeakProminence', subject_settings.left_legswing_peak_prominence_threshold, 'MinPeakDistance', subject_settings.left_legswing_peak_distance_threshold * sampling_rate_marker);
            [~, right_leg_peak_locations] = findpeaks(-rleg_angle, 'MinPeakProminence', subject_settings.right_legswing_peak_prominence_threshold, 'MinPeakDistance', subject_settings.right_legswing_peak_distance_threshold * sampling_rate_marker);
            
            % calculate time derivatives
%             filter_order = 4;
%             cutoff_frequency = 10; % cutoff frequency, in Hz
%             [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate_marker/2));	% set filter parameters for butterworth filter: 2=order of filter;
%             larm_angle_dot = deriveByTime(nanfiltfilt(b, a, larm_angle), 1/sampling_rate_marker);
%             rarm_angle_dot = deriveByTime(nanfiltfilt(b, a, rarm_angle), 1/sampling_rate_marker);
%             lleg_angle_dot = deriveByTime(nanfiltfilt(b, a, lleg_angle), 1/sampling_rate_marker);
%             rleg_angle_dot = deriveByTime(nanfiltfilt(b, a, rleg_angle), 1/sampling_rate_marker);

            % normalize
            [larm_angle_normalized, larm_angle_dot_normalized] = normalizePeriodicVariable(larm_angle, time_marker, left_arm_peak_locations);
            [rarm_angle_normalized, rarm_angle_dot_normalized] = normalizePeriodicVariable(rarm_angle, time_marker, right_arm_peak_locations);
            [lleg_angle_normalized, lleg_angle_dot_normalized] = normalizePeriodicVariable(lleg_angle, time_marker, left_leg_peak_locations);
            [rleg_angle_normalized, rleg_angle_dot_normalized] = normalizePeriodicVariable(rleg_angle, time_marker, right_leg_peak_locations);
            
            
            % calculate phase
            larm_phase = atan(larm_angle_dot_normalized ./ larm_angle_normalized);
            rarm_phase = atan(rarm_angle_dot_normalized ./ rarm_angle_normalized);
            lleg_phase = atan(lleg_angle_dot_normalized ./ lleg_angle_normalized);
            rleg_phase = atan(rleg_angle_dot_normalized ./ rleg_angle_normalized);
            
            
%             figure; hold on;
%             plot(time_marker, larm_angle_normalized);
%             plot(time_marker, larm_angle_dot_normalized);
            
%             figure; hold on;
%             plot3(larm_angle_normalized, larm_angle_normalized_dot, time_marker);
            
            % visualize
%             if visualize
            if true
                color_peak = [0.7, 0.2, 0.07];
                
                % left arm
                figure; axes; hold on; title('left armswing')
%                 plot(time_marker, larm_angle, 'linewidth', 1, 'displayname', 'left arm angle');
%                 plot(time_marker(left_arm_peak_locations), larm_angle(left_arm_peak_locations), '+', 'linewidth', 2, 'color', color_peak);
                plot(time_marker, larm_angle_normalized, 'linewidth', 1, 'displayname', 'left arm angle');
                plot(time_marker(left_arm_peak_locations), larm_angle_normalized(left_arm_peak_locations), '+', 'linewidth', 2, 'color', color_peak);
                plot(time_marker, larm_angle_dot_normalized, 'linewidth', 1, 'displayname', 'left arm angle dot');
                plot(time_marker, larm_phase, 'linewidth', 1, 'displayname', 'left arm phase');
                legend('show')
                
                % right arm
                figure; axes; hold on; title('right armswing')
%                 plot(time_marker, rarm_angle, 'linewidth', 1, 'displayname', 'right arm angle');
%                 plot(time_marker(right_arm_peak_locations), rarm_angle(right_arm_peak_locations), '+', 'linewidth', 2, 'color', color_peak);
                plot(time_marker, rarm_angle_normalized, 'linewidth', 1, 'displayname', 'right arm angle');
                plot(time_marker(right_arm_peak_locations), rarm_angle_normalized(right_arm_peak_locations), '+', 'linewidth', 2, 'color', color_peak);
                plot(time_marker, rarm_angle_dot_normalized, 'linewidth', 1, 'displayname', 'right arm angle dot');
                plot(time_marker, rarm_phase, 'linewidth', 1, 'displayname', 'right arm phase');
                
                % left leg
                figure; axes; hold on; title('left legswing')
%                 plot(time_marker, lleg_angle, 'linewidth', 1, 'displayname', 'left leg angle');
%                 plot(time_marker(left_leg_peak_locations), lleg_angle(left_leg_peak_locations), '+', 'linewidth', 2, 'color', color_peak);
                plot(time_marker, lleg_angle_normalized, 'linewidth', 1, 'displayname', 'left leg angle');
                plot(time_marker(left_leg_peak_locations), lleg_angle_normalized(left_leg_peak_locations), '+', 'linewidth', 2, 'color', color_peak);
                plot(time_marker, lleg_angle_dot_normalized, 'linewidth', 1, 'displayname', 'left leg angle dot');
                plot(time_marker, lleg_phase, 'linewidth', 1, 'displayname', 'left leg phase');
                
                % right leg
                figure; axes; hold on; title('right legswing')
%                 plot(time_marker, rleg_angle, 'linewidth', 1, 'displayname', 'right leg angle');
%                 plot(time_marker(right_leg_peak_locations), rleg_angle(right_leg_peak_locations), '+', 'linewidth', 2, 'color', color_peak);
                plot(time_marker, rleg_angle_normalized, 'linewidth', 1, 'displayname', 'right leg angle');
                plot(time_marker(right_leg_peak_locations), rleg_angle_normalized(right_leg_peak_locations), '+', 'linewidth', 2, 'color', color_peak);
                plot(time_marker, rleg_angle_dot_normalized, 'linewidth', 1, 'displayname', 'right leg angle dot');
                plot(time_marker, rleg_phase, 'linewidth', 1, 'displayname', 'right leg phase');
            end            
            
            
        end
    end
end













