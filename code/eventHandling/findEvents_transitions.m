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

function findEvents_transitions(varargin)

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
%             [cop_trajectories, time_forceplate, ~, ~, ~, ~] = loadData(collection_date, subject_id, condition, i_trial, 'total_forceplate_cop_world', 'optional');
            
            [forceplate_trajectories, time_forceplate, ~, labels_forceplate] = loadData(collection_date, subject_id, condition, i_trial, 'forceplate_trajectories', 'optional');
            cop_x = forceplate_trajectories(:, strcmp(labels_forceplate, 'copx'));
            cop_y = forceplate_trajectories(:, strcmp(labels_forceplate, 'copy'));
            
            [left_foot_wrench_world, time_left_forceplate, ~, ~, ~, left_forceplate_available] = loadData(collection_date, subject_id, condition, i_trial, 'left_foot_wrench_world', 'optional');
            [right_foot_wrench_world, time_right_forceplate, ~, ~, ~, right_forceplate_available] = loadData(collection_date, subject_id, condition, i_trial, 'right_foot_wrench_world', 'optional');
            if left_forceplate_available & right_forceplate_available %#ok<AND2>
                left_fz_trajectory = left_foot_wrench_world(:, 3);
                right_fz_trajectory = right_foot_wrench_world(:, 3);
            end
            
            %% find events
            
%             LHEE_z_trajectory = getKinematics('heel', 'left', 'z');
%             % find touch down indices as negative peaks of the heel marker z-position
%             [~, left_heel_peak_locations] = findpeaks(-LHEE_z_trajectory, 'MinPeakProminence', subject_settings.get('left_touchdown_peak_prominence_threshold'), 'MinPeakDistance', subject_settings.get('left_touchdown_peak_distance_threshold') * sampling_rate_marker);
%             left_touchdown_indices_mocap = left_heel_peak_locations';
%             left_touchdown_times = [left_touchdown_times; time_marker(left_touchdown_indices_mocap)]; %#ok<AGROW>
            
% ------------------------------------------------------------------------------------------------------------            
% the lines above to find left heel-strikes, delete them and add code to determine your own events here
% ------------------------------------------------------------------------------------------------------------            
                        
            


            %% save
            % struct for saving
            variables_to_save = struct;
            
% ------------------------------------------------------------------------------------------------------------            
% replace the event variable names and labels with yours here
% ------------------------------------------------------------------------------------------------------------            
            event_data = ...
              { ...
%                 left_pushoff_times; ...
%                 left_touchdown_times; ...
%                 left_fullstance_times; ...
%                 right_pushoff_times; ...
%                 right_touchdown_times; ...
%                 right_fullstance_times; ...
              };
            event_labels = ...
              { ...
%                 'left_pushoff'; ...
%                 'left_touchdown'; ...
%                 'left_fullstance'; ...
%                 'right_pushoff'; ...
%                 'right_touchdown'; ...
%                 'right_fullstance'; ...
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














