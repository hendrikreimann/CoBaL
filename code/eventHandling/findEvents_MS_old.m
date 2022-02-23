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

function findEvents_MS_old(varargin)

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
    subject_settings = SettingsCustodian('subjectSettings.txt');
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            %% prepare
            % load data
            this_trial_type = condition_list{i_condition};
            [marker_trajectories, time_marker, sampling_rate_marker, marker_labels] = loadData(collection_date, subject_id, this_trial_type, i_trial, 'marker_trajectories');
            [left_foot_wrench_world, time_left_forceplate, ~, ~, ~, left_forceplate_available] = loadData(collection_date, subject_id, this_trial_type, i_trial, 'left_foot_wrench_world', 'optional');
            [right_foot_wrench_world, time_right_forceplate, ~, ~, ~, right_forceplate_available] = loadData(collection_date, subject_id, this_trial_type, i_trial, 'right_foot_wrench_world', 'optional');
            if left_forceplate_available & right_forceplate_available
                left_fz_trajectory = left_foot_wrench_world(:, 3);
                right_fz_trajectory = right_foot_wrench_world(:, 3);
            end
            
            platform_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'BackPlatform');
            platform_trajectory = platform_trajectory(:, 1);
            surround_trajectory = extractMarkerData(marker_trajectories, marker_labels, 'BackSurround');
            surround_trajectory = surround_trajectory(:, 1);
            
            % calculate derivatives
            filter_order = 2;
            cutoff_frequency = 20; % cutoff frequency, in Hz
            [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate_marker/2));	% set filter parameters for butterworth filter: 2=order of filter;
            platform_vel_trajectory = deriveByTime(nanfiltfilt(b, a, platform_trajectory), 1/sampling_rate_marker);
            surround_vel_trajectory = deriveByTime(nanfiltfilt(b, a, surround_trajectory), 1/sampling_rate_marker);
            platform_acc_trajectory = deriveByTime(nanfiltfilt(b, a, platform_vel_trajectory), 1/sampling_rate_marker);
            surround_acc_trajectory = deriveByTime(nanfiltfilt(b, a, surround_vel_trajectory), 1/sampling_rate_marker);
            

            %% find events
            platform_events = [];
            if strcmp(this_trial_type, 'continuousPerturbation')
                [~, platform_pos_peak_indices] = findpeaks(platform_trajectory, 'MinPeakProminence', subject_settings.get('platform_pos_peak_prominence_threshold'), 'MinPeakDistance', subject_settings.get('platform_pos_peak_distance_threshold') * sampling_rate_marker);
                [~, platform_neg_peak_indices] = findpeaks(-platform_trajectory, 'MinPeakProminence', subject_settings.get('platform_pos_peak_prominence_threshold'), 'MinPeakDistance', subject_settings.get('platform_pos_peak_distance_threshold') * sampling_rate_marker);
%                 [~, surround_pos_peak_indices] = findpeaks(surround_trajectory, 'MinPeakProminence', subject_settings.get('platform_pos_peak_prominence_threshold'), 'MinPeakDistance', subject_settings.get('platform_pos_peak_distance_threshold') * sampling_rate_marker);
                
                platform_pos_peak_indices = platform_pos_peak_indices';
                platform_neg_peak_indices = platform_neg_peak_indices';
%                 surround_pos_peak_indices = surround_pos_peak_indices';

                platform_event_indices = [platform_pos_peak_indices platform_neg_peak_indices];
                platform_peak_indices = platform_pos_peak_indices;
                platform_peak_events = time_marker(platform_peak_indices);
                platform_vale_indices = platform_neg_peak_indices;
                platform_vale_events = time_marker(platform_vale_indices);
                number_of_peaks = length(platform_peak_indices);
                number_of_vales = length(platform_vale_indices);

                % store
                event_data = cell(number_of_peaks + number_of_vales, 1);
                event_labels = cell(number_of_peaks + number_of_vales, 1);
                for i_peak = 1 : number_of_peaks
                    event_data{i_peak} = platform_peak_events(i_peak);
                    event_labels{i_peak} = ['peak_' zeroPrefixedIntegerString(i_peak, 2)];
                end
                for i_vale = 1 : number_of_vales
                    event_data{number_of_peaks + i_vale} = platform_vale_events(i_vale);
                    event_labels{number_of_peaks + i_vale} = ['vale_' zeroPrefixedIntegerString(i_vale, 2)];
                end
            end
            
            if length(this_trial_type) >= 4 && strcmp(this_trial_type(1:4), 'ramp')
                distance_threshold = subject_settings.get('platform_acc_peak_distance_threshold') * sampling_rate_marker;
                if distance_threshold > length(platform_acc_trajectory) - 2
                    distance_threshold = length(platform_acc_trajectory) - 2;
                end
                [~, platform_acc_peak_indices] = findpeaks(platform_acc_trajectory, 'MinPeakProminence', subject_settings.get('platform_acc_peak_prominence_threshold'), 'MinPeakDistance', distance_threshold);
                platform_acc_peak_times = time_marker(platform_acc_peak_indices);
                [~, platform_acc_vale_indices] = findpeaks(-platform_acc_trajectory, 'MinPeakProminence', subject_settings.get('platform_acc_peak_prominence_threshold'), 'MinPeakDistance', distance_threshold);
                platform_acc_vale_times = time_marker(platform_acc_vale_indices);
%                 platform_event_indices = [platform_acc_peak_indices; platform_acc_vale_indices];
                
                platform_shift_start_candidates = 2;
                platform_shift_end_candidates = 2 + [1.2 3.6 6 8.4 12] * 1/15; % platform shifted by 15cm/s
                ramp_times = subject_settings.get('ramp_times');
                
                platform_shift_start_time = 2;
                if length(this_trial_type) >= 7 && strcmp(this_trial_type(5:7), '012')
%                     platform_shift_end_time = platform_shift_start_time + 1.2 * 1/15;
                    platform_shift_end_time = platform_shift_start_time + ramp_times(1);
                elseif length(this_trial_type) >= 7 && strcmp(this_trial_type(5:7), '036')
%                     platform_shift_end_time = platform_shift_start_time + 3.6 * 1/15;
                    platform_shift_end_time = platform_shift_start_time + ramp_times(2);
                elseif length(this_trial_type) >= 7 && strcmp(this_trial_type(5:7), '060')
%                     platform_shift_end_time = platform_shift_start_time + 6 * 1/15;
                    platform_shift_end_time = platform_shift_start_time + ramp_times(3);
                elseif length(this_trial_type) >= 7 && strcmp(this_trial_type(5:7), '084')
%                     platform_shift_end_time = platform_shift_start_time + 7.4 * 1/15;
                    platform_shift_end_time = platform_shift_start_time + ramp_times(4);
                elseif length(this_trial_type) >= 7 && strcmp(this_trial_type(5:7), '120')
%                     platform_shift_end_time = platform_shift_start_time + 12 * 1/15;
                    platform_shift_end_time = platform_shift_start_time + ramp_times(5);
                else
                    platform_shift_end_time = platform_shift_start_time;
                end
                
%                 % determine start
%                 [~, start_index] = min(abs(platform_acc_vale_times - platform_shift_start_candidates));
%                 platform_shift_start_time = platform_shift_start_candidates(start_index);
%                 [~, end_index] = min(abs(platform_acc_peak_times - platform_shift_end_candidates));
%                 platform_shift_end_time = platform_shift_end_candidates(end_index);
                
                event_data = ...
                  { ...
                    platform_shift_start_time - 1; ...
                    platform_shift_start_time - 0.2; ...
                    platform_shift_start_time; ...
                    platform_shift_start_time + 0.075; ...
                    platform_shift_end_time; ...
                    platform_shift_end_time + 1; ...
                    platform_shift_end_time + 2; ...
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
              
                platform_events = [platform_shift_start_time, platform_shift_end_time, platform_shift_end_time + 1, platform_shift_end_time + 2];
                platform_event_indices = findClosestIndex(platform_events, time_marker);
              
            end

            if length(this_trial_type) >= 5 && strcmp(this_trial_type(1:5), 'quiet')
                trial_start_time = time_marker(1);
                trial_end_time = time_marker(end);
                
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
                plot(time_marker, platform_trajectory, 'linewidth', 1, 'displayname', 'platform position');
                plot(time_marker, surround_trajectory, 'linewidth', 1, 'displayname', 'surround position');
                
                plot(time_marker(platform_event_indices), platform_trajectory(platform_event_indices), 'v', 'linewidth', 2, 'color', color_heelstrike, 'displayname', 'oscillation peaks');
%                 plot(time_marker(left_pushoff_indices_mocap), LTOE_z_trajectory(left_pushoff_indices_mocap), 'o', 'linewidth', 2, 'color', color_pushoff, 'displayname', 'left pushoff');
                legend('toggle');

                % velocities
                event_figures(2) = figure; axes_vel = axes; hold on;  title('marker velocities')
                plot(time_marker, platform_vel_trajectory, 'linewidth', 1, 'displayname', 'platform velocity');
                plot(time_marker, surround_vel_trajectory, 'linewidth', 1, 'displayname', 'surround velocity');
%                 plot(time_marker(left_touchdown_indices_mocap), platform_acc_trajectory(left_touchdown_indices_mocap)*acc_scaler, 'v', 'linewidth', 2, 'color', color_heelstrike, 'displayname', 'left touchdown');
%                 plot(time_marker(left_pushoff_indices_mocap), LTOE_z_vel_trajectory(left_pushoff_indices_mocap)*vel_scaler, '^', 'linewidth', 2, 'color', color_pushoff, 'displayname', 'left pushoff');
                legend('toggle');

                % accelerations
                event_figures(3) = figure; axes_acc = axes; hold on;  title('marker accelerations')
                plot(time_marker, platform_acc_trajectory, 'linewidth', 1, 'displayname', 'platform acceleration');
                plot(time_marker, surround_acc_trajectory, 'linewidth', 1, 'displayname', 'surround acceleration');
                plot(time_marker(platform_event_indices), platform_acc_trajectory(platform_event_indices), 'v', 'linewidth', 2, 'color', color_heelstrike, 'displayname', 'oscillation peaks');
                legend('toggle');
                
                linkaxes([axes_pos axes_vel axes_acc], 'x')
            end


            %% save
            % struct for saving
            variables_to_save = struct;
            
            % add new variables to be saved
            variables_to_save.event_data = event_data;
            variables_to_save.event_labels = event_labels;
            
            step_events_file_name = ['analysis' filesep makeFileName(collection_date, subject_id, this_trial_type, i_trial, 'events')];
            saveDataToFile(step_events_file_name, variables_to_save);

            disp(['Finding Events: condition ' this_trial_type ', Trial ' num2str(i_trial) ' completed, saved as ' step_events_file_name]);
        end
    end
end















