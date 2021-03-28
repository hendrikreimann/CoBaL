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

function findEvents_MS(varargin)

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

    % load settings
    subject_settings = loadSettingsFromFile('subject');
    study_settings = loadSettingsFromFile('study');
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            event_data = struct;
            
            %% prepare
            % load data
            this_trial_type = condition_list{i_condition};
            [marker_trajectories, time_marker, sampling_rate_marker, marker_labels] = loadData(subject_settings.get('collection_date'), subject_settings.get('subject_id'), this_trial_type, i_trial, 'marker_trajectories');
            [treadmill_trajectories, time_treadmill, sampling_rate_treadmill, treadmill_labels] = loadData(subject_settings.get('collection_date'), subject_settings.get('subject_id'), this_trial_type, i_trial, 'treadmill_trajectories', 'optional');
            [left_foot_wrench_world, time_left_forceplate, ~, ~, ~, left_forceplate_available] = loadData(subject_settings.get('collection_date'), subject_settings.get('subject_id'), this_trial_type, i_trial, 'left_foot_wrench_world', 'optional');
            [right_foot_wrench_world, time_right_forceplate, ~, ~, ~, right_forceplate_available] = loadData(subject_settings.get('collection_date'), subject_settings.get('subject_id'), this_trial_type, i_trial, 'right_foot_wrench_world', 'optional');
            if left_forceplate_available & right_forceplate_available
                left_fz_trajectory = left_foot_wrench_world(:, 3);
                right_fz_trajectory = right_foot_wrench_world(:, 3);
            end
            
            % derive platform marker by time
            platform_marker_label = subject_settings.get('base_marker', true);
            if ~isempty(platform_marker_label)
                platform_trajectory = extractMarkerData(marker_trajectories, marker_labels, platform_marker_label);
                platform_trajectory = platform_trajectory(:, 1);

                % calculate derivatives
                filter_order = 4;
                cutoff_frequency = 20; % cutoff frequency, in Hz
                [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate_marker/2));	% set filter parameters for butterworth filter: 2=order of filter;
                platform_vel_trajectory = deriveByTime(nanfiltfilt(b, a, platform_trajectory), 1/sampling_rate_marker);
                platform_acc_trajectory = deriveByTime(nanfiltfilt(b, a, platform_vel_trajectory), 1/sampling_rate_marker);
            end

            % derive treadmill markers by time
            if ~isempty(treadmill_trajectories)
                % calculate derivatives
                filter_order = 2;
                cutoff_frequency = 20; % cutoff frequency, in Hz
                [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate_marker/2));	% set filter parameters for butterworth filter: 2=order of filter;
                treadmill_vel_trajectories = deriveByTime(nanfiltfilt(b, a, treadmill_trajectories), 1/sampling_rate_marker);
                treadmill_acc_trajectories = deriveByTime(nanfiltfilt(b, a, treadmill_vel_trajectories), 1/sampling_rate_marker);
            end

            shift_source = subject_settings.get('shift_source', true);
            if strcmp(shift_source, 'base_marker')
                shift_vel_trajectory = platform_vel_trajectory;
                shift_acc_trajectory = platform_acc_trajectory;
            elseif strcmp(shift_source, 'treadmill_trajectory')
                shift_vel_trajectory = treadmill_vel_trajectories;
                shift_acc_trajectory = treadmill_acc_trajectories;
            end
            

            %% find events
            platform_events = [];
            if any(strcmp(this_trial_type, study_settings.get('across_events_conditions', 1)))
%             if strcmp(condition_experimental(1:end-3), 'continuous')
                [~, platform_pos_peak_indices] = findpeaks(platform_trajectory, 'MinPeakProminence', subject_settings.get('platform_pos_peak_prominence_threshold'), 'MinPeakDistance', subject_settings.get('platform_pos_peak_distance_threshold') * sampling_rate_marker);
                [~, platform_neg_peak_indices] = findpeaks(-platform_trajectory, 'MinPeakProminence', subject_settings.get('platform_pos_peak_prominence_threshold'), 'MinPeakDistance', subject_settings.get('platform_pos_peak_distance_threshold') * sampling_rate_marker);
                
                platform_pos_peak_indices = platform_pos_peak_indices';
                platform_neg_peak_indices = platform_neg_peak_indices';

                platform_event_indices = [platform_pos_peak_indices platform_neg_peak_indices];
                platform_peak_indices = platform_pos_peak_indices;
                platform_peak_events = time_marker(platform_peak_indices);
                platform_vale_indices = platform_neg_peak_indices;
                platform_vale_events = time_marker(platform_vale_indices);
                
                event_data = ...
                  { ...
                    platform_peak_events; ...
                    platform_vale_events; ...
                  };
                event_labels = ...
                  { ...
                    'oscillation_peaks'; ...
                    'oscillation_vales'; ...
                  };
            end
            
            if any(strcmp(this_trial_type, study_settings.get('across_trials_conditions', 1)))
                distance_threshold = subject_settings.get('acc_peak_distance_threshold') * sampling_rate_marker;
                if distance_threshold > length(shift_acc_trajectory) - 2
                    distance_threshold = length(shift_acc_trajectory) - 2;
                end
                
                
                if strcmp(subject_settings.get('shift_method', 1), 'acceleration_peaks')
                    marker_acceleration_ap = shift_acc_trajectory(:, 1 : 3 : end);
                    number_of_shift_markers = size(marker_acceleration_ap, 2);
                    shift_start_indices = zeros(1, number_of_shift_markers);
                    shift_end_indices = zeros(1, number_of_shift_markers);
                    shift_start_times = zeros(1, number_of_shift_markers);
                    shift_end_times = zeros(1, number_of_shift_markers);
                    for i_marker = 1 : number_of_shift_markers
                        % use acceleration peaks
                        [~, acc_peak_indices] = findpeaks(marker_acceleration_ap(:, i_marker), 'MinPeakProminence', subject_settings.get('acc_peak_prominence_threshold'), 'MinPeakDistance', distance_threshold);
                        acc_peak_times = time_marker(acc_peak_indices);
                        [~, acc_vale_indices] = findpeaks(-marker_acceleration_ap(:, i_marker), 'MinPeakProminence', subject_settings.get('acc_peak_prominence_threshold'), 'MinPeakDistance', distance_threshold);
                        acc_vale_times = time_marker(acc_vale_indices);

                        % determine start and end
                        shift_start_candidate_index = min([acc_peak_indices; acc_vale_indices]);
                        shift_end_candidate_index = max([acc_peak_indices; acc_vale_indices]);
                        shift_start_times(i_marker) = time_marker(shift_start_candidate_index);
                        shift_end_times(i_marker) = time_marker(shift_end_candidate_index);

                    end
                    if std(shift_start_times) > subject_settings.get('marker_timing_tolerance')
                        warning(['Trial ' this_trial_type ' - ' num2str(i_trial) ': timing between markers exceeds tolerance'])
                    end
                    if std(shift_end_times) > subject_settings.get('marker_timing_tolerance')
                        warning(['Trial ' this_trial_type ' - ' num2str(i_trial) ': timing between markers exceeds tolerance'])
                    end

                    shift_start_time = mean(shift_start_times);
                    shift_end_time = mean(shift_end_times);
                end
                
                if strcmp(subject_settings.get('shift_method', 1), 'velocity_threshold')
                    % use velocity threshold crossing
                    marker_velocity_ap_abs = abs(shift_vel_trajectory(:, 1 : 3 : end));
                    
                    vel_threshold_breach = marker_velocity_ap_abs > subject_settings.get('shift_vel_threshold', 1);
                    shift_start_index = find(diff(vel_threshold_breach) > 0, 1);
                    shift_end_times = find(diff(vel_threshold_breach) < 0, 1);

                    shift_start_time = time_marker(shift_start_index);
                    shift_end_time = time_marker(shift_end_times);
                end

                if strcmp(subject_settings.get('shift_method', 1), 'acceleration_threshold')
                    % use absolute acceleration threshold crossing for start
                    marker_acceleration_ap = shift_acc_trajectory(:, 1 : 3 : end);
                    marker_acceleration_ap_abs = abs(marker_acceleration_ap);
                    acc_abs_threshold_breach = marker_acceleration_ap_abs > subject_settings.get('shift_acc_threshold', 1);
                    shift_start_index = find(diff(acc_abs_threshold_breach) > 0, 1, 'first');
                    
                    % determine if this was forward or backward
                    shift_direction = sign(shift_acc_trajectory(shift_start_index));
                    acc_threshold_breach = (- shift_direction * marker_acceleration_ap) > subject_settings.get('shift_acc_threshold', 1);
                    shift_end_index = find(diff(acc_threshold_breach) > 0, 1, 'first');
                    
                    shift_start_time = time_marker(shift_start_index);
                    shift_end_time = time_marker(shift_end_index);
                end
                
                event_data = ...
                  { ...
                    shift_start_time - 1; ...
                    shift_start_time - 0.2; ...
                    shift_start_time; ...
                    shift_start_time + 0.075; ...
                    shift_end_time; ...
                    shift_end_time + 1; ...
                    shift_end_time + 2; ...
                  };
                event_labels = ...
                  { ...
                    'perturbation_start_minus_1s'; ...
                    'perturbation_start_minus_200ms'; ...
                    'perturbation_start'; ...
                    'perturbation_start_plus_75ms'; ...
                    'perturbation_end'; ...
                    'perturbation_end_plus_1s'; ...
                    'perturbation_end_plus_2s'; ...
                  };
                
                % for visualization here
                shift_events = [shift_start_time, shift_end_time, shift_end_time + 1, shift_end_time + 2];
                shift_event_indices = findClosestIndex(shift_events, time_marker);
              
            end

            if any(strcmp(this_trial_type, study_settings.get('quiet_stance_conditions')))
                first_nanless_data_point = find(sum(isnan(marker_trajectories), 2) == 0, 1, 'first');
                last_nanless_data_point = find(sum(isnan(marker_trajectories), 2) == 0, 1, 'last');
                trial_start_time = time_marker(first_nanless_data_point);
                trial_end_time = time_marker(last_nanless_data_point);
                
                event_data = ...
                  { ...
                    trial_start_time; ...
                    trial_end_time; ...
                  };
                event_labels = ...
                  { ...
                    'trial_start_time'; ...
                    'trial_end_time'; ...
                  };
              
                % for visualization here
                platform_events = [trial_start_time, trial_end_time];
                platform_event_indices = findClosestIndex(platform_events, time_marker);
              
            end
            
            %% visualize
            color_heelstrike = [1 0 0];
            color_pushoff = [0 1 0];
            color_peak = [0.5, 0.2, 1];

            force_scaler = 2e-4;
            vel_scaler = 0.1;
            acc_scaler = 0.01;

            if visualize
                event_figures = zeros(1, 3);
                
                % position
                event_figures(1) = figure; axes_pos = axes; hold on; title('marker positions')
%                 plot(time_marker, platform_trajectory, 'linewidth', 1, 'displayname', 'shift pos');
                plot(time_marker, shift_vel_trajectory, 'linewidth', 1, 'displayname', 'shift vel');
                plot(time_marker, shift_acc_trajectory, 'linewidth', 1, 'displayname', 'shift acc');
                
                
%                 plot(time_marker, surround_trajectory, 'linewidth', 1, 'displayname', 'surround position');
                
%                 plot(time_marker(platform_event_indices), platform_trajectory(platform_event_indices), 'v', 'linewidth', 2, 'color', color_heelstrike, 'displayname', 'oscillation peaks');
                [~, shift_start_time_index] = min(abs(time_marker - shift_start_time));
                [~, shift_end_time_index] = min(abs(time_marker - shift_end_time));
%                 plot(shift_start_time, shift_vel_trajectory(shift_start_time_index), '>', 'linewidth', 2, 'color', color_heelstrike, 'displayname', 'start');
%                 plot(shift_end_time, shift_vel_trajectory(shift_end_time_index), '<', 'linewidth', 2, 'color', color_heelstrike, 'displayname', 'end');
                plot(shift_start_time, shift_acc_trajectory(shift_start_time_index), '>', 'linewidth', 2, 'color', color_heelstrike, 'displayname', 'start');
                plot(shift_end_time, shift_acc_trajectory(shift_end_time_index), '<', 'linewidth', 2, 'color', color_heelstrike, 'displayname', 'end');
                
%                 plot(time_marker(left_pushoff_indices_mocap), LTOE_z_trajectory(left_pushoff_indices_mocap), 'o', 'linewidth', 2, 'color', color_pushoff, 'displayname', 'left pushoff');
                legend('toggle');

%                 % velocities
%                 event_figures(2) = figure; axes_vel = axes; hold on;  title('marker velocities')
%                 plot(time_marker, platform_vel_trajectory, 'linewidth', 1, 'displayname', 'platform velocity');
% %                 plot(time_marker, surround_vel_trajectory, 'linewidth', 1, 'displayname', 'surround velocity');
% %                 plot(time_marker(left_touchdown_indices_mocap), platform_acc_trajectory(left_touchdown_indices_mocap)*acc_scaler, 'v', 'linewidth', 2, 'color', color_heelstrike, 'displayname', 'left touchdown');
% %                 plot(time_marker(left_pushoff_indices_mocap), LTOE_z_vel_trajectory(left_pushoff_indices_mocap)*vel_scaler, '^', 'linewidth', 2, 'color', color_pushoff, 'displayname', 'left pushoff');
%                 legend('toggle');
% 
%                 % accelerations
%                 event_figures(3) = figure; axes_acc = axes; hold on;  title('marker accelerations')
%                 plot(time_marker, platform_acc_trajectory, 'linewidth', 1, 'displayname', 'platform acceleration');
% %                 plot(time_marker, surround_acc_trajectory, 'linewidth', 1, 'displayname', 'surround acceleration');
%                 plot(time_marker(platform_event_indices), platform_acc_trajectory(platform_event_indices), 'v', 'linewidth', 2, 'color', color_heelstrike, 'displayname', 'oscillation peaks');
%                 legend('toggle');
%                 
%                 linkaxes([axes_pos axes_vel axes_acc], 'x')
            end


            %% save
            % struct for saving
            variables_to_save = struct;
            
            % add new variables to be saved
            variables_to_save.event_data = event_data;
            variables_to_save.event_labels = event_labels;
            
            step_events_file_name = ['analysis' filesep makeFileName(subject_settings.get('collection_date'), subject_settings.get('subject_id'), this_trial_type, i_trial, 'events')];
            saveDataToFile(step_events_file_name, variables_to_save);

            disp(['Finding Events: condition ' this_trial_type ', Trial ' num2str(i_trial) ' completed, saved as ' step_events_file_name]);
        end
    end
end















