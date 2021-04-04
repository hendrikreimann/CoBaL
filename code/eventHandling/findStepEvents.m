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
    addParameter(parser, 'VisualizationWindow', [])
    parse(parser, varargin{:})
    visualize = parser.Results.visualize;
    visualization_window = parser.Results.VisualizationWindow;

    % figure out folders
    if ~exist('analysis', 'dir')
        mkdir('analysis')
    end
    
    % load settings
    subject_settings = loadSettingsFromFile('subject');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            %% prepare
            % load data
            condition = condition_list{i_condition};
            [marker_trajectories, time_marker, sampling_rate_marker, marker_labels] = loadData(collection_date, subject_id, condition, i_trial, 'marker_trajectories');
            [cop_trajectories, time_forceplate, ~, ~, ~, ~] = loadData(collection_date, subject_id, condition, i_trial, 'total_forceplate_cop_world', 'optional');
            [left_foot_wrench_world, time_left_forceplate, ~, ~, ~, left_forceplate_available] = loadData(collection_date, subject_id, condition, i_trial, 'left_foot_wrench_world', 'optional');
            [right_foot_wrench_world, time_right_forceplate, ~, ~, ~, right_forceplate_available] = loadData(collection_date, subject_id, condition, i_trial, 'right_foot_wrench_world', 'optional');
            if left_forceplate_available & right_forceplate_available %#ok<AND2>
                left_fz_trajectory = left_foot_wrench_world(:, 3);
                right_fz_trajectory = right_foot_wrench_world(:, 3);
            end
            
            %% find events
            left_touchdown_times = [];
            left_pushoff_times = [];
            right_touchdown_times = [];
            right_pushoff_times = [];
            
            
            if any(strcmp(subject_settings.get('event_method', 1), 'cop_ap'))
                % get pushoff and touchdown indices
                cop_ap = cop_trajectories(:, 2);
                [~, pushoff_indices_both] = findpeaks(cop_ap, 'MinPeakProminence', subject_settings.get('peak_prominence_threshold'));
                [~, touchdown_indices_both] = findpeaks(-cop_ap, 'MinPeakProminence', subject_settings.get('peak_prominence_threshold'));
                
                LHEE_y_trajectory = getKinematics('heel', 'left', 'y');
                RHEE_y_trajectory = getKinematics('heel', 'right', 'y');
                
                % go through and assign each index to a foot
                pushoff_indices_left = [];
                pushoff_indices_right = [];
                for i_pushoff = 1 : length(pushoff_indices_both)
                    this_index_forceplate = pushoff_indices_both(i_pushoff);
                    this_time = time_forceplate(this_index_forceplate);
                    this_index_mocap = findClosestIndex(this_time, time_marker);
                    
                    left_heel_pos_ap_here = LHEE_y_trajectory(this_index_mocap);
                    right_heel_pos_ap_here = RHEE_y_trajectory(this_index_mocap);
                    if left_heel_pos_ap_here < right_heel_pos_ap_here
                        pushoff_indices_left = [pushoff_indices_left; this_index_forceplate]; %#ok<AGROW>
                        left_pushoff_times = [left_pushoff_times; this_time]; %#ok<AGROW>
                    end
                    if  right_heel_pos_ap_here < left_heel_pos_ap_here
                        pushoff_indices_right = [pushoff_indices_right; this_index_forceplate]; %#ok<AGROW>
                        right_pushoff_times = [right_pushoff_times; this_time]; %#ok<AGROW>
                    end
                    
                end
                touchdown_indices_left = [];
                touchdown_indices_right = [];
                for i_touchdown = 1 : length(touchdown_indices_both)
                    this_index_forceplate = touchdown_indices_both(i_touchdown);
                    this_time = time_forceplate(this_index_forceplate);
                    this_index_mocap = findClosestIndex(this_time, time_marker);
                    
                    left_heel_pos_ap_here = LHEE_y_trajectory(this_index_mocap);
                    right_heel_pos_ap_here = RHEE_y_trajectory(this_index_mocap);
                    if left_heel_pos_ap_here > right_heel_pos_ap_here
                        touchdown_indices_left = [touchdown_indices_left; this_index_forceplate]; %#ok<AGROW>
                        left_touchdown_times = [left_touchdown_times; this_time]; %#ok<AGROW>
                    end
                    if  right_heel_pos_ap_here > left_heel_pos_ap_here
                        touchdown_indices_right = [touchdown_indices_right; this_index_forceplate]; %#ok<AGROW>
                        right_touchdown_times = [right_touchdown_times; this_time]; %#ok<AGROW>
                    end
                end
                
                
                if visualize
                    if ~isempty(visualization_window)
                        window_frame_indices = findClosestIndex(visualization_window, time_forceplate);
                    else
                        window_frame_indices = [1, length(time_forceplate)];
                    end
                    window_indices = window_frame_indices(1) : window_frame_indices(2);
                    pushoff_indices_left_window = pushoff_indices_left(pushoff_indices_left>=window_frame_indices(1) & pushoff_indices_left<=window_frame_indices(2));
                    pushoff_indices_right_window = pushoff_indices_right(pushoff_indices_right>=window_frame_indices(1) & pushoff_indices_right<=window_frame_indices(2));
                    touchdown_indices_left_window = touchdown_indices_left(touchdown_indices_left>=window_frame_indices(1) & touchdown_indices_left<=window_frame_indices(2));
                    touchdown_indices_right_window = touchdown_indices_right(touchdown_indices_right>=window_frame_indices(1) & touchdown_indices_right<=window_frame_indices(2));
                    
                    figure; axes; hold on; title('CoP path on treadmill')
                    plot(cop_trajectories(window_indices, 1), cop_trajectories(window_indices, 2), 'color', 'b', 'DisplayName', 'CoP');
                    plot(cop_trajectories(pushoff_indices_left_window, 1), cop_trajectories(pushoff_indices_left_window, 2), '^', 'markersize', 8, 'color', 'g', 'markerfacecolor', 'g', 'DisplayName', 'left pushoff');
                    plot(cop_trajectories(pushoff_indices_right_window, 1), cop_trajectories(pushoff_indices_right_window, 2), '^', 'markersize', 8, 'color', 'r', 'markerfacecolor', 'r', 'DisplayName', 'right pushoff');
                    plot(cop_trajectories(touchdown_indices_left_window, 1), cop_trajectories(touchdown_indices_left_window, 2), 'v', 'markersize', 8, 'color', 'g', 'markerfacecolor', 'g', 'DisplayName', 'left heel-strike');
                    plot(cop_trajectories(touchdown_indices_right_window, 1), cop_trajectories(touchdown_indices_right_window, 2), 'v', 'markersize', 8, 'color', 'r', 'markerfacecolor', 'r', 'DisplayName', 'right heel-strike');
                    xlabel('ml position (m)'); ylabel('ap position (m)')

                    figure; axes_x = axes; hold on; title('CoP trajectories - left forceplate')
                    plot(time_forceplate(window_indices), cop_trajectories(window_indices, 1), 'linewidth', 2, 'color', 'b', 'DisplayName', 'CoP');
                    plot(time_forceplate(pushoff_indices_left_window), cop_trajectories(pushoff_indices_left_window, 1), '^', 'markersize', 8, 'color', 'g', 'markerfacecolor', 'g', 'DisplayName', 'left pushoff');
                    plot(time_forceplate(pushoff_indices_right_window), cop_trajectories(pushoff_indices_right_window, 1), '^', 'markersize', 8, 'color', 'r', 'markerfacecolor', 'r', 'DisplayName', 'right pushoff');
                    plot(time_forceplate(touchdown_indices_left_window), cop_trajectories(touchdown_indices_left_window, 1), 'v', 'markersize', 8, 'color', 'g', 'markerfacecolor', 'g', 'DisplayName', 'left heel-strike');
                    plot(time_forceplate(touchdown_indices_right_window), cop_trajectories(touchdown_indices_right_window, 1), 'v', 'markersize', 8, 'color', 'r', 'markerfacecolor', 'r', 'DisplayName', 'right heel-strike');
                    xlabel('time (s)'); ylabel('ml position (m)')

                    figure; axes_y = axes; hold on; title('CoP trajectories - right forceplate')
                    plot(time_forceplate(window_indices), cop_trajectories(window_indices, 2), 'linewidth', 2, 'color', 'b', 'DisplayName', 'CoP');
                    plot(time_forceplate(pushoff_indices_left_window), cop_trajectories(pushoff_indices_left_window, 2), '^', 'markersize', 8, 'color', 'g', 'markerfacecolor', 'g', 'DisplayName', 'left pushoff');
                    plot(time_forceplate(pushoff_indices_right_window), cop_trajectories(pushoff_indices_right_window, 2), '^', 'markersize', 8, 'color', 'r', 'markerfacecolor', 'r', 'DisplayName', 'right pushoff');
                    plot(time_forceplate(touchdown_indices_left_window), cop_trajectories(touchdown_indices_left_window, 2), 'v', 'markersize', 8, 'color', 'g', 'markerfacecolor', 'g', 'DisplayName', 'left heel-strike');
                    plot(time_forceplate(touchdown_indices_right_window), cop_trajectories(touchdown_indices_right_window, 2), 'v', 'markersize', 8, 'color', 'r', 'markerfacecolor', 'r', 'DisplayName', 'right heel-strike');
                    xlabel('time (s)'); ylabel('ap position (m)')

                    linkaxes([axes_x, axes_y], 'x');
                end
            end
                        
            %% find events for left foot
            if any(strcmp(subject_settings.get('left_touchdown_method', 1), 'heel_position_minima'))
                LHEE_z_trajectory = getKinematics('heel', 'left', 'z');
                % find touch down indices as negative peaks of the heel marker z-position
                [~, left_heel_peak_locations] = findpeaks(-LHEE_z_trajectory, 'MinPeakProminence', subject_settings.get('left_touchdown_peak_prominence_threshold'), 'MinPeakDistance', subject_settings.get('left_touchdown_peak_distance_threshold') * sampling_rate_marker);
                left_touchdown_indices_mocap = left_heel_peak_locations';
                left_touchdown_times = [left_touchdown_times; time_marker(left_touchdown_indices_mocap)]; %#ok<AGROW>
            end
            if any(strcmp(subject_settings.get('left_touchdown_method', 1), 'toe_position_minima'))
                LTOE_z_trajectory = getKinematics('toes', 'left', 'z');
                % find touch down indices as negative peaks of the heel marker z-position
                [~, left_toe_peak_locations] = findpeaks(-LTOE_z_trajectory, 'MinPeakProminence', subject_settings.get('left_touchdown_peak_prominence_threshold'), 'MinPeakDistance', subject_settings.get('left_touchdown_peak_distance_threshold') * sampling_rate_marker);
                left_touchdown_indices_mocap = left_toe_peak_locations';
                left_touchdown_times = [left_touchdown_times; time_marker(left_touchdown_indices_mocap)]; %#ok<AGROW>
            end
            if any(strcmp(subject_settings.get('left_touchdown_method', 1), 'toe_velocity_minima'))
                [~, LTOE_z_vel_trajectory] = getKinematics('toes', 'left', 'z');
                [~, left_toe_peak_locations] = findpeaks(-LTOE_z_vel_trajectory, 'MinPeakProminence', subject_settings.get('left_touchdown_peak_prominence_threshold'), 'MinPeakDistance', subject_settings.get('left_touchdown_peak_distance_threshold') * sampling_rate_marker);
                left_touchdown_indices_mocap = left_toe_peak_locations';
                left_touchdown_times = [left_touchdown_times; time_marker(left_touchdown_indices_mocap)]; %#ok<AGROW>
            end
            if any(strcmp(subject_settings.get('left_touchdown_method', 1), 'first_acceleration_peak_after_mid_swing'))
                [LHEE_z_trajectory, ~, LHEE_z_acc_trajectory] = getKinematics('heel', 'left', 'z');
                
                % find mid-swing as peaks of heel position
                [~, left_heel_midswing_locations] = findpeaks(LHEE_z_trajectory, 'MinPeakDistance', subject_settings.get('left_touchdown_peak_distance_threshold') * sampling_rate_marker);
                
                % find acceleration peaks
                [~, left_heel_acc_peak_locations] = findpeaks(LHEE_z_acc_trajectory, 'MinPeakProminence', subject_settings.get('left_touchdown_peak_prominence_threshold'));

                % identify acceleration peaks as touchdowns
                left_touchdown_indices_mocap = zeros(size(left_heel_midswing_locations));
                for i_step = 1 : length(left_heel_midswing_locations)
                    touchdown_index_index = left_heel_acc_peak_locations(find(left_heel_acc_peak_locations > left_heel_midswing_locations(i_step), 1, 'first'));
                    if ~isempty(touchdown_index_index)
                        left_touchdown_indices_mocap(i_step) = touchdown_index_index;
                    end
                end
                left_touchdown_indices_mocap(left_touchdown_indices_mocap==0) = [];
                left_touchdown_times = [left_touchdown_times; time_marker(left_touchdown_indices_mocap)]; %#ok<AGROW>
            end
            if any(strcmp(subject_settings.get('left_touchdown_method', 1), 'heel_ap_peak'))
                LHEE_y_trajectory = getKinematics('heel', 'left', 'y');
                
                % find velocity zeros
                [~, left_heel_peak_locations] = findpeaks(LHEE_y_trajectory, 'MinPeakProminence', subject_settings.get('left_touchdown_peak_prominence_threshold'));
                left_touchdown_indices_mocap = left_heel_peak_locations;
                left_touchdown_times = [left_touchdown_times; time_marker(left_heel_peak_locations)]; %#ok<AGROW>
            end
                
            if any(strcmp(subject_settings.get('left_pushoff_method', 1), 'first_velocity_peak_after_touchdown'))
                [~, LTOE_z_vel_trajectory] = getKinematics('toes', 'left', 'z');
                % for pushoff, find the first significant toes z-velocity peak after each touchdown
                [~, left_toes_vel_peak_locations] = findpeaks(LTOE_z_vel_trajectory, 'MinPeakProminence', subject_settings.get('left_pushoff_peak_prominence_threshold'));
                left_toes_vel_peak_locations = left_toes_vel_peak_locations';
                left_pushoff_indices_mocap = zeros(size(left_touchdown_indices_mocap));
                for i_touchdown = 1 : length(left_touchdown_indices_mocap)
                    pushoff_index_index = find(left_toes_vel_peak_locations > left_touchdown_indices_mocap(i_touchdown), 1, 'first');
                    if ~isempty(pushoff_index_index)
                        left_pushoff_indices_mocap(i_touchdown) = left_toes_vel_peak_locations(pushoff_index_index);
                    end
                end
                left_pushoff_indices_mocap(left_pushoff_indices_mocap==0) = [];
                left_pushoff_times = [left_pushoff_times; time_marker(left_pushoff_indices_mocap)]; %#ok<AGROW>
            end
            if any(strcmp(subject_settings.get('left_pushoff_method', 1), 'heel_acceleration_minima'))
                [~, ~, LHEE_z_acc_trajectory] = getKinematics('heel', 'left', 'z');
                % find touch down indices as negative peaks of the heel marker z-position
                [~, left_heel_acc_peak_locations] = findpeaks(-LHEE_z_acc_trajectory, 'MinPeakProminence', subject_settings.get('left_pushoff_vale_prominence_threshold'), 'MinPeakDistance', subject_settings.get('left_pushoff_vale_distance_threshold') * sampling_rate_marker);
                left_pushoff_indices_mocap = left_heel_acc_peak_locations';
                left_pushoff_times = [left_pushoff_times; time_marker(left_pushoff_indices_mocap)]; %#ok<AGROW>
            end
            if any(strcmp(subject_settings.get('left_pushoff_method', 1), 'forceplate_threshold'))
                left_pushoff_diff_forceplate = diff(sign(abs(left_fz_trajectory) - subject_settings.get('forceplate_load_threshold')));
                left_pushoff_times = [left_pushoff_times; time_left_forceplate(left_pushoff_diff_forceplate~=0)]; %#ok<AGROW>
            end
            if any(strcmp(subject_settings.get('left_pushoff_method', 1), 'toes_ap_vale'))
                LTOE_y_trajectory = getKinematics('toes', 'left', 'y');
                % find velocity zeros
                [~, left_toes_vale_locations] = findpeaks(-LTOE_y_trajectory, 'MinPeakProminence', subject_settings.get('left_pushoff_peak_prominence_threshold'));
                left_pushoff_times = [left_pushoff_times; time_marker(left_toes_vale_locations)]; %#ok<AGROW>
            end
            
            left_fullstance_times = [];
            if any(strcmp(subject_settings.get('left_fullstance_method', 1), 'first_zero_crossing_after_heelstrike'))
                % identify zero crossings of angle
                left_foot_angle_trajectory = calculateFootAngleKinematics('left');
                left_foot_angle_zero_crossing_indices = find(diff(sign(left_foot_angle_trajectory)));
                
                % identify acceleration peaks as touchdowns
                left_fullstance_indices_mocap = zeros(size(left_touchdown_times));
                for i_step = 1 : length(left_touchdown_times)
                    this_touchdown_time = left_touchdown_times(i_step);
                    zero_crossing_index = left_foot_angle_zero_crossing_indices(find((time_marker(left_foot_angle_zero_crossing_indices) > this_touchdown_time), 1, 'first'));
                    if ~isempty(zero_crossing_index)
                        left_fullstance_indices_mocap(i_step) = zero_crossing_index;
                    end
                end
                left_fullstance_indices_mocap(left_fullstance_indices_mocap==0) = [];
                left_fullstance_times = [left_fullstance_times; time_marker(left_fullstance_indices_mocap)]; %#ok<AGROW>
            end

            %% find events for right foot
            if any(strcmp(subject_settings.get('right_touchdown_method', 1), 'heel_position_minima'))
                RHEE_z_trajectory = getKinematics('heel', 'right', 'z');
                % find touch down indices as negative peaks of the heel marker z-position
                [~, right_heel_peak_locations] = findpeaks(-RHEE_z_trajectory, 'MinPeakProminence', subject_settings.get('right_touchdown_peak_prominence_threshold'), 'MinPeakDistance', subject_settings.get('right_touchdown_peak_distance_threshold'));
                right_touchdown_indices_mocap = right_heel_peak_locations';
                right_touchdown_times = [right_touchdown_times; time_marker(right_touchdown_indices_mocap)]; %#ok<AGROW>
            end
            if any(strcmp(subject_settings.get('right_touchdown_method', 1), 'toe_position_minima'))
                RTOE_z_trajectory = getKinematics('toes', 'right', 'z');
                % find touch down indices as negative peaks of the heel marker z-position
                [~, right_toe_peak_locations] = findpeaks(-RTOE_z_trajectory, 'MinPeakProminence', subject_settings.get('right_touchdown_peak_prominence_threshold'), 'MinPeakDistance', subject_settings.get('right_touchdown_peak_distance_threshold'));
                right_touchdown_indices_mocap = right_toe_peak_locations';
                right_touchdown_times = [right_touchdown_times; time_marker(right_touchdown_indices_mocap)]; %#ok<AGROW>
            end
            if any(strcmp(subject_settings.get('right_touchdown_method', 1), 'toe_velocity_minima'))
                [~, RTOE_z_vel_trajectory] = getKinematics('toes', 'right', 'z');
                [~, right_toe_peak_locations] = findpeaks(-RTOE_z_vel_trajectory, 'MinPeakProminence', subject_settings.get('right_touchdown_peak_prominence_threshold'), 'MinPeakDistance', subject_settings.get('right_touchdown_peak_distance_threshold') * sampling_rate_marker);
                right_touchdown_indices_mocap = right_toe_peak_locations';
                right_touchdown_times = [right_touchdown_times; time_marker(right_touchdown_indices_mocap)]; %#ok<AGROW>
            end
            if any(strcmp(subject_settings.get('right_touchdown_method', 1), 'first_acceleration_peak_after_mid_swing'))
                [RHEE_z_trajectory, ~, RHEE_z_acc_trajectory] = getKinematics('heel', 'right', 'z');
                
                % find mid-swing as peaks of heel position
                [~, right_heel_midswing_locations] = findpeaks(RHEE_z_trajectory, 'MinPeakDistance', subject_settings.get('right_touchdown_peak_distance_threshold') * sampling_rate_marker);
                
                % find acceleration peaks
                [~, right_heel_acc_peak_locations] = findpeaks(RHEE_z_acc_trajectory, 'MinPeakProminence', subject_settings.get('right_touchdown_peak_prominence_threshold'));

                % identify acceleration peaks as touchdowns
                right_touchdown_indices_mocap = zeros(size(right_heel_midswing_locations));
                for i_step = 1 : length(right_heel_midswing_locations)
                    touchdown_index_index = right_heel_acc_peak_locations(find(right_heel_acc_peak_locations > right_heel_midswing_locations(i_step), 1, 'first'));
                    if ~isempty(touchdown_index_index)
                        right_touchdown_indices_mocap(i_step) = touchdown_index_index;
                    end
                end
                right_touchdown_indices_mocap(right_touchdown_indices_mocap==0) = [];
                right_touchdown_times = [right_touchdown_times; time_marker(right_touchdown_indices_mocap)]; %#ok<AGROW>
            end
            if any(strcmp(subject_settings.get('right_touchdown_method', 1), 'heel_ap_peak'))
                RHEE_y_trajectory = getKinematics('heel', 'right', 'y');
                
                % find velocity zeros
                [~, right_heel_peak_locations] = findpeaks(RHEE_y_trajectory, 'MinPeakProminence', subject_settings.get('right_touchdown_peak_prominence_threshold'));
                right_touchdown_indices_mocap = right_heel_peak_locations;
                right_touchdown_times = [right_touchdown_times; time_marker(right_heel_peak_locations)]; %#ok<AGROW>
            end
  
            if any(strcmp(subject_settings.get('right_pushoff_method', 1), 'first_velocity_peak_after_touchdown'))
                [~, RTOE_z_vel_trajectory] = getKinematics('toes', 'right', 'z');
                % for pushoff, find the first significant toes z-velocity peak after each touchdown
                [~, right_toes_vel_peak_locations] = findpeaks(RTOE_z_vel_trajectory, 'MinPeakProminence', subject_settings.get('right_pushoff_peak_prominence_threshold'));
                right_toes_vel_peak_locations = right_toes_vel_peak_locations';
                right_pushoff_indices_mocap = zeros(size(right_touchdown_indices_mocap));
                for i_touchdown = 1 : length(right_touchdown_indices_mocap)
                    pushoff_index_index = find(right_toes_vel_peak_locations > right_touchdown_indices_mocap(i_touchdown), 1, 'first');
                    if ~isempty(pushoff_index_index)
                        right_pushoff_indices_mocap(i_touchdown) = right_toes_vel_peak_locations(pushoff_index_index);
                    end
                end
                right_pushoff_indices_mocap(right_pushoff_indices_mocap==0) = [];
                right_pushoff_times = [right_pushoff_times; time_marker(right_pushoff_indices_mocap)]; %#ok<AGROW>
            end
            if any(strcmp(subject_settings.get('right_pushoff_method', 1), 'heel_acceleration_minima'))
                [~, ~, RHEE_z_acc_trajectory] = getKinematics('heel', 'right', 'z');
                % find touch down indices as negative peaks of the heel marker z-position
                [~, right_heel_acc_peak_locations] = findpeaks(-RHEE_z_acc_trajectory, 'MinPeakProminence', subject_settings.get('right_pushoff_vale_prominence_threshold'), 'MinPeakDistance', subject_settings.get('right_pushoff_vale_distance_threshold') * sampling_rate_marker);
                right_pushoff_indices_mocap = right_heel_acc_peak_locations';
                right_pushoff_times = [right_pushoff_times; time_marker(right_pushoff_indices_mocap)]; %#ok<AGROW>
            end
            if any(strcmp(subject_settings.get('right_pushoff_method', 1), 'forceplate_threshold'))
                right_pushoff_diff_forceplate = diff(sign(abs(right_fz_trajectory) - subject_settings.get('forceplate_load_threshold')));
                right_pushoff_times = [right_pushoff_times; time_right_forceplate(right_pushoff_diff_forceplate~=0)]; %#ok<AGROW>
            end
            if any(strcmp(subject_settings.get('right_pushoff_method', 1), 'toes_ap_vale'))
                RTOE_y_trajectory = getKinematics('toes', 'right', 'y');
                % find velocity zeros
                [~, right_toes_vale_locations] = findpeaks(-RTOE_y_trajectory, 'MinPeakProminence', subject_settings.get('right_pushoff_peak_prominence_threshold'));
                right_pushoff_times = [right_pushoff_times; time_marker(right_toes_vale_locations)]; %#ok<AGROW>
            end

            right_fullstance_times = [];
            if any(strcmp(subject_settings.get('right_fullstance_method', 1), 'first_zero_crossing_after_heelstrike'))
                % identify zero crossings of angle
                right_foot_angle_trajectory = calculateFootAngleKinematics('right');
                right_foot_angle_zero_crossing_indices = find(diff(sign(right_foot_angle_trajectory)));
                
                % identify acceleration peaks as touchdowns
                right_fullstance_indices_mocap = zeros(size(right_touchdown_times));
                for i_step = 1 : length(right_touchdown_times)
                    this_touchdown_time = right_touchdown_times(i_step);
                    zero_crossing_index = right_foot_angle_zero_crossing_indices(find((time_marker(right_foot_angle_zero_crossing_indices) > this_touchdown_time), 1, 'first'));
                    if ~isempty(zero_crossing_index)
                        right_fullstance_indices_mocap(i_step) = zero_crossing_index;
                    end
                end
                right_fullstance_indices_mocap(right_fullstance_indices_mocap==0) = [];
                right_fullstance_times = [right_fullstance_times; time_marker(right_fullstance_indices_mocap)]; %#ok<AGROW>
            end
            
            % remove duplicates
            left_pushoff_times = unique(left_pushoff_times);
            left_touchdown_times = unique(left_touchdown_times);
            left_fullstance_times = unique(left_fullstance_times);
            right_pushoff_times = unique(right_pushoff_times);
            right_touchdown_times = unique(right_touchdown_times);
            right_fullstance_times = unique(right_fullstance_times);

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
            
            %% visualize
            color_heelstrike = [1 0 0];
            color_pushoff = [0 1 0];
%             color_peak = [0.5, 0.2, 1];

            force_scaler = 2e-4;
            vel_scaler = 0.1;
            acc_scaler = 0.01;



            if left_forceplate_available & right_forceplate_available %#ok<AND2>
                left_fz_trajectory_marker = spline(time_left_forceplate, left_fz_trajectory, time_marker);
                right_fz_trajectory_marker = spline(time_right_forceplate, right_fz_trajectory, time_marker);
            end
            if visualize
                
                if ~isempty(visualization_window)
                    window_frame_indices = findClosestIndex(visualization_window, time_marker);
                    if left_forceplate_available
                        window_frame_indices_forceplate = findClosestIndex(visualization_window, time_left_forceplate);
                    end
                else
                    window_frame_indices = [1, length(time_marker)];
                    if left_forceplate_available
                        window_frame_indices_forceplate = [1, length(time_left_forceplate)];
                    end
                end
                window_indices = window_frame_indices(1) : window_frame_indices(2);
                if left_forceplate_available
                    window_indices_forceplate = window_frame_indices_forceplate(1) : window_frame_indices_forceplate(2);
                end
                left_pushoff_indices_mocap_window = left_pushoff_indices_mocap(left_pushoff_indices_mocap>=window_frame_indices(1) & left_pushoff_indices_mocap<=window_frame_indices(2));
                right_pushoff_indices_mocap_window = right_pushoff_indices_mocap(right_pushoff_indices_mocap>=window_frame_indices(1) & right_pushoff_indices_mocap<=window_frame_indices(2));
                left_touchdown_indices_mocap_window = left_touchdown_indices_mocap(left_touchdown_indices_mocap>=window_frame_indices(1) & left_touchdown_indices_mocap<=window_frame_indices(2));
                right_touchdown_indices_mocap_window = right_touchdown_indices_mocap(right_touchdown_indices_mocap>=window_frame_indices(1) & right_touchdown_indices_mocap<=window_frame_indices(2));
                
                
                
                step_event_figures = zeros(1, 4);
                [LHEE_z_trajectory, ~, LHEE_z_acc_trajectory] = getKinematics('heel', 'left', 'z');
                [RHEE_z_trajectory, ~, RHEE_z_acc_trajectory] = getKinematics('heel', 'right', 'z');
                [LTOE_z_trajectory, LTOE_z_vel_trajectory] = getKinematics('toes', 'left', 'z');
                [RTOE_z_trajectory, RTOE_z_vel_trajectory] = getKinematics('toes', 'right', 'z');
                
                % left position
                step_event_figures(1) = figure; axes_left = axes; hold on; title('left foot marker positions')
                plot(time_marker(window_indices), LHEE_z_trajectory(window_indices), 'linewidth', 2, 'displayname', 'left heel vertical');
                plot(time_marker(window_indices), LTOE_z_trajectory(window_indices), 'linewidth', 2, 'displayname', 'left toes vertical');
%                 plot(time_marker(window_indices), left_foot_angle_trajectory(window_indices), 'linewidth', 2, 'displayname', 'left foot angle');
%                 plot(time_marker(window_indices), left_foot_angle_vel_trajectory*vel_scaler, 'linewidth', 2, 'displayname', 'left foot angle');
%                 plot(time_marker(window_indices), left_foot_angle_acc_trajectory*acc_scaler, 'linewidth', 2, 'displayname', 'left foot angle');
                plot(time_marker(left_touchdown_indices_mocap_window), LHEE_z_trajectory(left_touchdown_indices_mocap_window), 'v', 'color', 'g', 'markersize', 8, 'markerfacecolor', 'g', 'displayname', 'left touchdown');
                plot(time_marker(left_pushoff_indices_mocap_window), LTOE_z_trajectory(left_pushoff_indices_mocap_window), '^', 'color', 'g', 'markersize', 8, 'markerfacecolor', 'g', 'displayname', 'left pushoff');
%                 plot(time_marker(left_fullstance_indices_mocap), left_foot_angle_trajectory(left_fullstance_indices_mocap), '^', 'color', color_peak, 'displayname', 'left fullstance');
%                 plot(time_marker(left_leg_swing_onset_indices), LTOE_z_trajectory(left_leg_swing_onset_indices), '+', 'color', color_peak, 'displayname', 'left leg angle peaks');
                if left_forceplate_available
                    plot(time_left_forceplate(window_indices_forceplate), left_fz_trajectory(window_indices_forceplate)*force_scaler, 'linewidth', 2, 'displayname', 'left forceplate f_z')
                    h = plot(time_marker(left_touchdown_indices_mocap_window), left_fz_trajectory_marker(left_touchdown_indices_mocap_window)*force_scaler, 'v', 'markersize', 8, 'color', 'g', 'markerfacecolor', 'g');
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    h = plot(time_marker(left_pushoff_indices_mocap_window), left_fz_trajectory_marker(left_pushoff_indices_mocap_window)*force_scaler, '^', 'markersize', 8, 'color', 'g', 'markerfacecolor', 'g');
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                legend('toggle');

                % left derivatives
                step_event_figures(2) = figure; axes_left_derivatives = axes; hold on;  title('left foot marker derivatives (scaled for qualitative overview)')
                plot(time_marker(window_indices), LHEE_z_acc_trajectory(window_indices)*acc_scaler, 'linewidth', 2, 'displayname', 'left heel acceleration vertical');
                plot(time_marker(window_indices), LTOE_z_vel_trajectory(window_indices)*vel_scaler, 'linewidth', 2, 'displayname', 'left toes velocity vertical');
                plot(time_marker(left_touchdown_indices_mocap_window), LHEE_z_acc_trajectory(left_touchdown_indices_mocap_window)*acc_scaler, 'v', 'color', 'g', 'markersize', 8, 'markerfacecolor', 'g', 'displayname', 'left touchdown');
                plot(time_marker(left_pushoff_indices_mocap_window), LTOE_z_vel_trajectory(left_pushoff_indices_mocap_window)*vel_scaler, '^', 'color', 'g', 'markersize', 8, 'markerfacecolor', 'g', 'displayname', 'left pushoff');
                if left_forceplate_available
                    plot(time_left_forceplate(window_indices_forceplate), left_fz_trajectory(window_indices_forceplate)*force_scaler, 'linewidth', 2, 'displayname', 'left forceplate f_z')
                    h = plot(time_marker(left_touchdown_indices_mocap_window), left_fz_trajectory_marker(left_touchdown_indices_mocap_window)*force_scaler, 'v', 'markersize', 8, 'color', 'g', 'markerfacecolor', 'g');
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    h = plot(time_marker(left_pushoff_indices_mocap_window), left_fz_trajectory_marker(left_pushoff_indices_mocap_window)*force_scaler, '^', 'markersize', 8, 'color', 'g', 'markerfacecolor', 'g');
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                legend('toggle');

                % right position
                step_event_figures(3) = figure; axes_right = axes; hold on; title('right foot marker positions')
                plot(time_marker(window_indices), RHEE_z_trajectory(window_indices), 'linewidth', 2, 'displayname', 'right heel vertical');
                plot(time_marker(window_indices), RTOE_z_trajectory(window_indices), 'linewidth', 2, 'displayname', 'right toes vertical');
%                 plot(time_marker(window_indices), right_foot_angle_trajectory(window_indices), 'linewidth', 2, 'displayname', 'right foot angle');
%                 plot(time_marker(window_indices), right_foot_angle_vel_trajectory(window_indices)*vel_scaler, 'linewidth', 2, 'displayname', 'right foot angle');
%                 plot(time_marker(window_indices), right_foot_angle_acc_trajectory(window_indices)*acc_scaler, 'linewidth', 2, 'displayname', 'right foot angle');
                plot(time_marker(right_touchdown_indices_mocap_window), RHEE_z_trajectory(right_touchdown_indices_mocap_window), 'v', 'color', 'r', 'markersize', 8, 'markerfacecolor', 'r', 'displayname', 'right touchdown');
                plot(time_marker(right_pushoff_indices_mocap_window), RTOE_z_trajectory(right_pushoff_indices_mocap_window), '^', 'color', 'r', 'markersize', 8, 'markerfacecolor', 'r', 'displayname', 'right pushoff');
%                 plot(time_marker(right_fullstance_indices_mocap), right_foot_angle_trajectory(right_fullstance_indices_mocap), '^', 'color', color_peak, 'displayname', 'right fullstance');
                if right_forceplate_available
                    plot(time_right_forceplate(window_indices_forceplate), right_fz_trajectory(window_indices_forceplate)*force_scaler, 'linewidth', 2, 'displayname', 'right forceplate f_z')
                    h = plot(time_marker(right_touchdown_indices_mocap_window), right_fz_trajectory_marker(right_touchdown_indices_mocap_window)*force_scaler, 'v', 'markersize', 8, 'color', 'r', 'markerfacecolor', 'r');
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    h = plot(time_marker(right_pushoff_indices_mocap_window), right_fz_trajectory_marker(right_pushoff_indices_mocap_window)*force_scaler, '^', 'markersize', 8, 'color', 'r', 'markerfacecolor', 'r');
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                legend('toggle');

                % right derivatives
                step_event_figures(4) = figure; axes_right_derivatives = axes; hold on;  title('right foot marker derivatives (scaled for qualitative overview)') %#ok<NASGU>
                plot(time_marker(window_indices), RHEE_z_acc_trajectory(window_indices)*acc_scaler, 'linewidth', 2, 'displayname', 'right heel acceleration vertical');
                plot(time_marker(window_indices), RTOE_z_vel_trajectory(window_indices)*vel_scaler, 'linewidth', 2, 'displayname', 'right toes velocity vertical');
                plot(time_marker(right_touchdown_indices_mocap_window), RHEE_z_acc_trajectory(right_touchdown_indices_mocap_window)*acc_scaler, 'v', 'markersize', 8, 'color', 'r', 'markerfacecolor', 'r', 'displayname', 'right touchdown');
                plot(time_marker(right_pushoff_indices_mocap_window), RTOE_z_vel_trajectory(right_pushoff_indices_mocap_window)*vel_scaler, '^', 'markersize', 8, 'color', 'r', 'markerfacecolor', 'r', 'displayname', 'right pushoff');
                if right_forceplate_available
                    plot(time_right_forceplate(window_indices_forceplate), right_fz_trajectory(window_indices_forceplate)*force_scaler, 'linewidth', 2, 'displayname', 'right forceplate f_z')
                    h = plot(time_marker(right_touchdown_indices_mocap_window), right_fz_trajectory_marker(right_touchdown_indices_mocap_window)*force_scaler, 'v', 'markersize', 8, 'color', 'r', 'markerfacecolor', 'r');
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    h = plot(time_marker(right_pushoff_indices_mocap_window), right_fz_trajectory_marker(right_pushoff_indices_mocap_window)*force_scaler, '^', 'markersize', 8, 'color', 'r', 'markerfacecolor', 'r');
                    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                legend('toggle');

                linkaxes([axes_left axes_left_derivatives axes_right axes_right_derivatives], 'x')
                
%                 distFig(step_event_figures, 'rows', 2)
            end

            % change event variables to column vectors if necessary
            if isrow(left_pushoff_times)
                left_pushoff_times = left_pushoff_times';
            end
            if isrow(left_touchdown_times)
                left_touchdown_times = left_touchdown_times';
            end
            if isrow(right_pushoff_times)
                right_pushoff_times = right_pushoff_times';
            end
            if isrow(right_touchdown_times)
                right_touchdown_times = right_touchdown_times';
            end

            %% save
            % struct for saving
            variables_to_save = struct;
            
            event_data = ...
              { ...
                left_pushoff_times; ...
                left_touchdown_times; ...
                left_fullstance_times; ...
                right_pushoff_times; ...
                right_touchdown_times; ...
                right_fullstance_times; ...
              };
            event_labels = ...
              { ...
                'left_pushoff'; ...
                'left_touchdown'; ...
                'left_fullstance'; ...
                'right_pushoff'; ...
                'right_touchdown'; ...
                'right_fullstance'; ...
              };
            
            % add new variables to be saved
            variables_to_save.event_data = event_data;
            variables_to_save.event_labels = event_labels;
            
            step_events_file_name = ['analysis' filesep makeFileName(collection_date, subject_id, condition, i_trial, 'events.mat')];
            saveDataToFile(step_events_file_name, variables_to_save);

            disp(['Finding Step Events: condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' step_events_file_name]);
        end
    end
    
    function [pos, vel, acc] = getKinematics(landmark, side, component)
        % landmark options: heel, toes
        % side options: left, right
        % component options: x, y, z
        
        % figure out side
        if strcmp(side, 'left')
            marker_label = 'L';
        elseif strcmp(side, 'right')
            marker_label = 'R';
        else
            error(['Value "' side '" not recognized for side. Options are "left" or "right".']);
        end
        
        % figure out marker
        if strcmp(landmark, 'heel')
            marker_label = [marker_label 'HEE'];
        elseif strcmp(landmark, 'toes')
            marker_label = [marker_label 'TOE'];
        else
            error(['Value "' landmark '" not recognized for landmark. Options are "heel" or "toes".']);
        end
        
        % figure out component
        if nargin < 3
            component = 'all';
        end
        if strcmp(component, 'x')
            component_index = 1;
        elseif strcmp(component, 'y')
            component_index = 2;
        elseif strcmp(component, 'z')
            component_index = 3;
        elseif strcmp(component, 'all')
            component_index = 1:3;
        else
            error(['Value "' component '" not recognized for component. Options are "x", "y" or "z".']);
        end
        
        marker_trajectory = extractMarkerData(marker_trajectories, marker_labels, marker_label, 'trajectories');
        pos = marker_trajectory(:, component_index);
        
        if nargout > 1
            filter_order = 2;
            cutoff_frequency = 20; % cutoff frequency, in Hz
            [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate_marker/2));	% set filter parameters for butterworth filter: 2=order of filter;
            marker_velocity_trajectory = deriveByTime(nanfiltfilt(b, a, marker_trajectory), 1/sampling_rate_marker);
            vel = marker_velocity_trajectory(:, component_index);
        end
        if nargout > 2
            acc = deriveByTime(nanfiltfilt(b, a, vel), 1/sampling_rate_marker);
        end

    end

    function [angle, vel, acc] = calculateFootAngleKinematics(side)
        toes_trajectory = getKinematics('toes', side, 'all');
        heel_trajectory = getKinematics('heel', side, 'all');
        
        foot_vector = toes_trajectory - heel_trajectory;
        angle = zeros(size(heel_trajectory, 1), 1);
        for i_time = 1 : size(heel_trajectory, 1)
            angle(i_time) = atan2(foot_vector(i_time, 3), foot_vector(i_time, 2));
        end
        
        if nargout > 1
            filter_order = 2;
            cutoff_frequency = 20; % cutoff frequency, in Hz
            [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate_marker/2));	% set filter parameters for butterworth filter: 2=order of filter;
            vel = deriveByTime(nanfiltfilt(b, a, angle), 1/sampling_rate_marker);
        end
        if nargout > 2
            acc = deriveByTime(nanfiltfilt(b, a, vel), 1/sampling_rate_marker);
        end
       
        
        
    end
    
end














