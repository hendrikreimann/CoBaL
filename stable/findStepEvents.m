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
    settings = loadSettingsFile('subjectSettings.txt');
    
    
%     show_forceplate         = 1;

    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process

            %% prepare
            % load data
            condition = condition_list{i_condition};
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
            
            file_name_forceplate = ['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'forceplateTrajectories.mat')];
            if exist(file_name_forceplate, 'file')
                forceplate_data_available = true;
                loaded_forceplate_trajectories = load(file_name_forceplate);
                time_forceplate = loaded_forceplate_trajectories.time_forceplate;

                left_fz_trajectory = loaded_forceplate_trajectories.left_forceplate_wrench_world(:, 3);
                right_fz_trajectory = loaded_forceplate_trajectories.right_forceplate_wrench_world(:, 3);
                
                % this was for the VEPO lab, where the sides are inverted
%                 left_fz_trajectory = loaded_forceplate_trajectories.right_forceplate_wrench_world(:, 3); % left foot is on the right forceplate
%                 right_fz_trajectory = loaded_forceplate_trajectories.left_forceplate_wrench_world(:, 3); % right foot is on the left forceplate
            else
                forceplate_data_available = false;
            end

            % extract data
            LHEE_trajectory = extractMarkerTrajectories(marker_trajectories, marker_headers, 'LHEE');
            LHEE_z_trajectory = LHEE_trajectory(:, 3);
            LTOE_trajectory = extractMarkerTrajectories(marker_trajectories, marker_headers, 'LTOE');
            LTOE_z_trajectory = LTOE_trajectory(:, 3);
            RHEE_trajectory = extractMarkerTrajectories(marker_trajectories, marker_headers, 'RHEE');
            RHEE_z_trajectory = RHEE_trajectory(:, 3);
            RTOE_trajectory = extractMarkerTrajectories(marker_trajectories, marker_headers, 'RTOE');
            RTOE_z_trajectory = RTOE_trajectory(:, 3);
            
            % calculate derivatives
            filter_order = 2;
            cutoff_frequency = 20; % cutoff frequency, in Hz
            [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate_mocap/2));	% set filter parameters for butterworth filter: 2=order of filter;
            LHEE_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, LHEE_z_trajectory), 1/sampling_rate_mocap);
            RHEE_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, RHEE_z_trajectory), 1/sampling_rate_mocap);
            LHEE_z_acc_trajectory = deriveByTime(nanfiltfilt(b, a, LHEE_z_vel_trajectory), 1/sampling_rate_mocap);
            RHEE_z_acc_trajectory = deriveByTime(nanfiltfilt(b, a, RHEE_z_vel_trajectory), 1/sampling_rate_mocap);
            LTOE_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, LTOE_z_trajectory), 1/sampling_rate_mocap);
            RTOE_z_vel_trajectory = deriveByTime(nanfiltfilt(b, a, RTOE_z_trajectory), 1/sampling_rate_mocap);

            %% find events for left foot
            left_touchdown_times = [];
            if any(strcmp(settings.left_touchdown_method, 'heel_position_minima'))
                % find touch down indices as negative peaks of the heel marker z-position
                [~, left_heel_peak_locations] = findpeaks(-LHEE_z_trajectory, 'MinPeakProminence', settings.left_touchdown_peak_prominence_threshold, 'MinPeakDistance', settings.left_touchdown_peak_distance_threshold * sampling_rate_mocap);
                left_touchdown_indices_mocap = left_heel_peak_locations';
                left_touchdown_times = [left_touchdown_times; time_mocap(left_touchdown_indices_mocap)];
            end
            if any(strcmp(settings.left_touchdown_method, 'toe_position_minima'))
                % find touch down indices as negative peaks of the heel marker z-position
                [~, left_toe_peak_locations] = findpeaks(-LTOE_z_trajectory, 'MinPeakProminence', settings.left_touchdown_peak_prominence_threshold, 'MinPeakDistance', settings.left_touchdown_peak_distance_threshold * sampling_rate_mocap);
                left_touchdown_indices_mocap = left_toe_peak_locations';
                left_touchdown_times = [left_touchdown_times; time_mocap(left_touchdown_indices_mocap)];
            end
            if any(strcmp(settings.left_touchdown_method, 'toe_velocity_minima'))
                [~, left_toe_peak_locations] = findpeaks(-LTOE_z_vel_trajectory, 'MinPeakProminence', settings.left_touchdown_peak_prominence_threshold, 'MinPeakDistance', settings.left_touchdown_peak_distance_threshold * sampling_rate_mocap);
                left_touchdown_indices_mocap = left_toe_peak_locations';
                left_touchdown_times = [left_touchdown_times; time_mocap(left_touchdown_indices_mocap)];
            end
            if any(strcmp(settings.left_touchdown_method, 'first_acceleration_peak_after_mid_swing'))
                % left_touchdown_peak_distance_threshold: 0.8
                % left touchdown peak prominence threshold: 5
                
                % find mid-swing as peaks of heel position
                [~, left_heel_midswing_locations] = findpeaks(LHEE_z_trajectory, 'MinPeakDistance', settings.left_touchdown_peak_distance_threshold * sampling_rate_mocap);
                
                % find acceleration peaks
                [~, left_heel_acc_peak_locations] = findpeaks(LHEE_z_acc_trajectory, 'MinPeakProminence', settings.left_touchdown_peak_prominence_threshold);

                % identify acceleration peaks as touchdowns
                left_touchdown_indices_mocap = zeros(size(left_heel_midswing_locations));
                for i_step = 1 : length(left_heel_midswing_locations)
                    touchdown_index_index = left_heel_acc_peak_locations(find(left_heel_acc_peak_locations > left_heel_midswing_locations(i_step), 1, 'first'));
                    if ~isempty(touchdown_index_index)
                        left_touchdown_indices_mocap(i_step) = touchdown_index_index;
                    end
                end
                left_touchdown_indices_mocap(left_touchdown_indices_mocap==0) = [];
                left_touchdown_times = [left_touchdown_times; time_mocap(left_touchdown_indices_mocap)];
            end

            left_pushoff_times = [];
            if any(strcmp(settings.left_pushoff_method, 'first_velocity_peak_after_touchdown'))
                % for pushoff, find the first significant toes z-velocity peak after each touchdown
                [~, left_toes_vel_peak_locations] = findpeaks(LTOE_z_vel_trajectory, 'MinPeakProminence', settings.left_pushoff_peak_prominence_threshold);
                left_toes_vel_peak_locations = left_toes_vel_peak_locations';
                left_pushoff_indices_mocap = zeros(size(left_touchdown_indices_mocap));
                for i_touchdown = 1 : length(left_touchdown_indices_mocap)
                    pushoff_index_index = find(left_toes_vel_peak_locations > left_touchdown_indices_mocap(i_touchdown), 1, 'first');
                    if ~isempty(pushoff_index_index)
                        left_pushoff_indices_mocap(i_touchdown) = left_toes_vel_peak_locations(pushoff_index_index);
                    end
                end
                left_pushoff_indices_mocap(left_pushoff_indices_mocap==0) = [];
                left_pushoff_times = [left_pushoff_times; time_mocap(left_pushoff_indices_mocap)];
            end
            if any(strcmp(settings.left_pushoff_method, 'forceplate_threshold'))
                left_pushoff_diff_forceplate = diff(sign(left_fz_trajectory - settings.forceplate_load_threshold));
                left_pushoff_times = [left_pushoff_times; time_forceplate(left_pushoff_diff_forceplate~=0)];
            end

            %% find events for right foot
            right_touchdown_times = [];
            if any(strcmp(settings.right_touchdown_method, 'heel_position_minima'))
                % find touch down indices as negative peaks of the heel marker z-position
                [~, right_heel_peak_locations] = findpeaks(-RHEE_z_trajectory, 'MinPeakProminence', settings.right_touchdown_peak_prominence_threshold, 'MinPeakDistance', settings.right_touchdown_peak_distance_threshold);
                right_touchdown_indices_mocap = right_heel_peak_locations';
                right_touchdown_times = [right_touchdown_times; time_mocap(right_touchdown_indices_mocap)];
            end
            if any(strcmp(settings.right_touchdown_method, 'toe_position_minima'))
                % find touch down indices as negative peaks of the heel marker z-position
                [~, right_toe_peak_locations] = findpeaks(-RTOE_z_trajectory, 'MinPeakProminence', settings.right_touchdown_peak_prominence_threshold, 'MinPeakDistance', settings.right_touchdown_peak_distance_threshold);
                right_touchdown_indices_mocap = right_toe_peak_locations';
                right_touchdown_times = [right_touchdown_times; time_mocap(right_touchdown_indices_mocap)];
            end
            if any(strcmp(settings.right_touchdown_method, 'toe_velocity_minima'))
                [~, right_toe_peak_locations] = findpeaks(-RTOE_z_vel_trajectory, 'MinPeakProminence', settings.right_touchdown_peak_prominence_threshold, 'MinPeakDistance', settings.right_touchdown_peak_distance_threshold * sampling_rate_mocap);
                right_touchdown_indices_mocap = right_toe_peak_locations';
                right_touchdown_times = [right_touchdown_times; time_mocap(right_touchdown_indices_mocap)];
            end
            if any(strcmp(settings.right_touchdown_method, 'first_acceleration_peak_after_mid_swing'))
                % right_touchdown_peak_distance_threshold: 0.8
                % right touchdown peak prominence threshold: 5
                
                % find mid-swing as peaks of heel position
                [~, right_heel_midswing_locations] = findpeaks(RHEE_z_trajectory, 'MinPeakDistance', settings.right_touchdown_peak_distance_threshold * sampling_rate_mocap);
                
                % find acceleration peaks
                [~, right_heel_acc_peak_locations] = findpeaks(RHEE_z_acc_trajectory, 'MinPeakProminence', settings.right_touchdown_peak_prominence_threshold);

                % identify acceleration peaks as touchdowns
                right_touchdown_indices_mocap = zeros(size(right_heel_midswing_locations));
                for i_step = 1 : length(right_heel_midswing_locations)
                    touchdown_index_index = right_heel_acc_peak_locations(find(right_heel_acc_peak_locations > right_heel_midswing_locations(i_step), 1, 'first'));
                    if ~isempty(touchdown_index_index)
                        right_touchdown_indices_mocap(i_step) = touchdown_index_index;
                    end
                end
                right_touchdown_indices_mocap(right_touchdown_indices_mocap==0) = [];
                right_touchdown_times = [right_touchdown_times; time_mocap(right_touchdown_indices_mocap)];
            end
  
            right_pushoff_times = []; 
            if any(strcmp(settings.right_pushoff_method, 'first_velocity_peak_after_touchdown'))
                % for pushoff, find the first significant toes z-velocity peak after each touchdown
                [~, right_toes_vel_peak_locations] = findpeaks(RTOE_z_vel_trajectory, 'MinPeakProminence', settings.right_pushoff_peak_prominence_threshold);
                right_toes_vel_peak_locations = right_toes_vel_peak_locations';
                right_pushoff_indices_mocap = zeros(size(right_touchdown_indices_mocap));
                for i_touchdown = 1 : length(right_touchdown_indices_mocap)
                    pushoff_index_index = find(right_toes_vel_peak_locations > right_touchdown_indices_mocap(i_touchdown), 1, 'first');
                    if ~isempty(pushoff_index_index)
                        right_pushoff_indices_mocap(i_touchdown) = right_toes_vel_peak_locations(pushoff_index_index);
                    end
                end
                right_pushoff_indices_mocap(right_pushoff_indices_mocap==0) = [];
                right_pushoff_times = [right_pushoff_times; time_mocap(right_pushoff_indices_mocap)];
            end
            if any(strcmp(settings.right_pushoff_method, 'forceplate_threshold'))
                right_pushoff_diff_forceplate = diff(sign(right_fz_trajectory - settings.forceplate_load_threshold));
                right_pushoff_times = [right_pushoff_times; time_forceplate(right_pushoff_diff_forceplate~=0)];
            end

            % remove duplicates
            left_pushoff_times = unique(left_pushoff_times);
            left_touchdown_times = unique(left_touchdown_times);
            right_pushoff_times = unique(right_pushoff_times);
            right_touchdown_times = unique(right_touchdown_times);

            % determine indices
            left_touchdown_indices_mocap = zeros(size(left_touchdown_times));
            for i_index = 1 : length(left_touchdown_times)
                [~, index_mocap] = min(abs(time_mocap - left_touchdown_times(i_index)));
                left_touchdown_indices_mocap(i_index) = index_mocap;
            end
            left_pushoff_indices_mocap = zeros(size(left_pushoff_times));
            for i_index = 1 : length(left_pushoff_times)
                [~, index_mocap] = min(abs(time_mocap - left_pushoff_times(i_index)));
                left_pushoff_indices_mocap(i_index) = index_mocap;
            end
            right_touchdown_indices_mocap = zeros(size(right_touchdown_times));
            for i_index = 1 : length(right_touchdown_times)
                [~, index_mocap] = min(abs(time_mocap - right_touchdown_times(i_index)));
                right_touchdown_indices_mocap(i_index) = index_mocap;
            end
            right_pushoff_indices_mocap = zeros(size(right_pushoff_times));
            for i_index = 1 : length(right_pushoff_times)
                [~, index_mocap] = min(abs(time_mocap - right_pushoff_times(i_index)));
                right_pushoff_indices_mocap(i_index) = index_mocap;
            end

            %% visualize
            color_heelstrike = [1 0 0];
            color_pushoff = [0 1 0];

            force_scaler = 2e-4;
            vel_scaler = 1;
            acc_scaler = .2;



            if forceplate_data_available
                left_fz_trajectory_mocap = spline(time_forceplate, left_fz_trajectory, time_mocap);
                right_fz_trajectory_mocap = spline(time_forceplate, right_fz_trajectory, time_mocap);
            end
            if visualize
                % position
                
                % left events
                figure; axes_left = axes; hold on; title('left foot marker positions')
                plot(time_mocap, LHEE_z_trajectory, 'linewidth', 1);
                plot(time_mocap, LTOE_z_trajectory, 'linewidth', 1);
                plot(time_mocap(left_touchdown_indices_mocap), LHEE_z_trajectory(left_touchdown_indices_mocap), 'v', 'linewidth', 2, 'color', color_heelstrike);
                plot(time_mocap(left_pushoff_indices_mocap), LTOE_z_trajectory(left_pushoff_indices_mocap), '^', 'linewidth', 2, 'color', color_pushoff);
                if forceplate_data_available
                    plot(time_forceplate, left_fz_trajectory*force_scaler)
                    legend('heel', 'toes', 'touchdown', 'pushoff', 'fzl');
                    plot(time_mocap(left_touchdown_indices_mocap), left_fz_trajectory_mocap(left_touchdown_indices_mocap)*force_scaler, 'v', 'linewidth', 2, 'color', color_heelstrike);
                    plot(time_mocap(left_pushoff_indices_mocap), left_fz_trajectory_mocap(left_pushoff_indices_mocap)*force_scaler, '^', 'linewidth', 2, 'color', color_pushoff);
                else
                    legend('heel', 'toes', 'touchdown', 'pushoff');
                end

                figure; axes_left_derivatives = axes; hold on;  title('left foot marker derivatives')
        %         plot(time_mocap, LHEE_z_vel_trajectory, 'linewidth', 1);
                plot(time_mocap, LHEE_z_acc_trajectory*acc_scaler, 'linewidth', 1);
                plot(time_mocap, LTOE_z_vel_trajectory*vel_scaler, 'linewidth', 1);
        %         plot(time_mocap, LTOE_z_acc_trajectory, 'linewidth', 1);
        %         plot(time_mocap(left_touchdown_indices_mocap), LHEE_z_acc_trajectory(left_touchdown_indices_mocap), 'o', 'linewidth', 2);
                legend('heel acc', 'toes vel');
                plot(time_mocap(left_touchdown_indices_mocap), LHEE_z_acc_trajectory(left_touchdown_indices_mocap)*acc_scaler, 'v', 'linewidth', 2, 'color', color_heelstrike);
                plot(time_mocap(left_pushoff_indices_mocap), LTOE_z_vel_trajectory(left_pushoff_indices_mocap)*vel_scaler, '^', 'linewidth', 2, 'color', color_pushoff);

                linkaxes([axes_left, axes_left_derivatives], 'x')
        %         distFig('rows', 2)

                % right events
                figure; axes_right = axes; hold on; title('right foot marker positions')
                plot(time_mocap, RHEE_z_trajectory, 'linewidth', 1);
                plot(time_mocap, RTOE_z_trajectory, 'linewidth', 1);
                plot(time_mocap(right_touchdown_indices_mocap), RHEE_z_trajectory(right_touchdown_indices_mocap), 'v', 'linewidth', 2, 'color', color_heelstrike);
                plot(time_mocap(right_pushoff_indices_mocap), RTOE_z_trajectory(right_pushoff_indices_mocap), '^', 'linewidth', 2, 'color', color_pushoff);
                if forceplate_data_available
                    plot(time_forceplate, right_fz_trajectory*force_scaler)
                    legend('heel', 'toes', 'touchdown', 'pushoff', 'fzr');
                    plot(time_mocap(right_touchdown_indices_mocap), right_fz_trajectory_mocap(right_touchdown_indices_mocap)*force_scaler, 'v', 'linewidth', 2, 'color', color_heelstrike);
                    plot(time_mocap(right_pushoff_indices_mocap), right_fz_trajectory_mocap(right_pushoff_indices_mocap)*force_scaler, '^', 'linewidth', 2, 'color', color_pushoff);
                else
                    legend('heel', 'toes', 'touchdown', 'pushoff');
                end

                figure; axes_right_derivatives = axes; hold on; title('right foot marker derivatives')
        %         plot(time_mocap, RHEE_z_vel_trajectory, 'linewidth', 1);
                plot(time_mocap, RHEE_z_acc_trajectory*acc_scaler, 'linewidth', 1);
                plot(time_mocap, RTOE_z_vel_trajectory*vel_scaler, 'linewidth', 1);
        %         plot(time_mocap, RTOE_z_acc_trajectory, 'linewidth', 1);
                legend('heel acc', 'toes vel');
                plot(time_mocap(right_touchdown_indices_mocap), RHEE_z_acc_trajectory(right_touchdown_indices_mocap)*acc_scaler, 'v', 'linewidth', 2, 'color', color_heelstrike);
                plot(time_mocap(right_pushoff_indices_mocap), RTOE_z_vel_trajectory(right_pushoff_indices_mocap)*vel_scaler, '^', 'linewidth', 2, 'color', color_pushoff);

                linkaxes([axes_right, axes_right_derivatives], 'x')
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
            step_events_file_name = ['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'stepEvents')];
            save ...
              ( ...
                step_events_file_name, ...
                'left_pushoff_times', ...
                'left_touchdown_times', ...
                'right_pushoff_times', ...
                'right_touchdown_times' ...
              );

            disp(['Condition ' condition ', Trial ' num2str(i_trial) ' completed, saved as ' step_events_file_name]);
        end
    end
    if ~isempty(findobj('type','figure'))
        distFig('rows', 2)
    end
end















