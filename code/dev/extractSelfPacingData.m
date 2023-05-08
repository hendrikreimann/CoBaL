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

function extractSelfPacingData(varargin)
    % parse arguments
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    trials_to_process_table = makeTrialsToProcessTable(condition_list, trial_number_list);

    %% prepare
    subject_settings = loadSettingsFromFile('subject');
    study_settings = loadSettingsFromFile('study');

    %% process
    for i_trial = 1 : height(trials_to_process_table)
        % create container for data
        trial_data = struct;
        
        % fill in basics
        trial_data.trial_type = trials_to_process_table{i_trial, 'trial type'};
        trial_data.trial_number = trials_to_process_table{i_trial, 'trial number'};

        % load belt speed data
        collection_date = subject_settings.get('collection_date');
        subject_id = subject_settings.get('subject_id');
        plc_data = load(['processed' filesep makeFileName(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'PLCData.mat')]);
        belt_velocity = (plc_data.belt_speed_left_trajectory + plc_data.belt_speed_right_trajectory) * 0.5;

        % rotate belt velocity according to forceplate frame
        v_belt = [0; -1; 0] * belt_velocity'; % default vector of belt movement
        R_fp_from_lab = reshape(study_settings.get('forceplate_rotation', true), 3, 3);
        v_belt_transformed = R_fp_from_lab * v_belt;
        % NOTE: this is tested properly only for flips of the y-axis

        % and integrate to belt translation
        belt_translation = cumtrapz(plc_data.time, v_belt_transformed');
        
        % load marker kinematics
        load_file_name = ['processed' filesep makeFileName(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'markerTrajectories.mat')];
        marker_data = load(load_file_name);
        
        % resample belt translation to marker time
        belt_translation_resampled = spline(plc_data.time, belt_translation', marker_data.time_mocap)';
        
        % add belt translation to ap-component of all data
        marker_data_belt = marker_data;
        for i_label = 1 : length(marker_data.marker_labels)
            this_label = marker_data.marker_labels{i_label};
            if ~isempty(this_label) && this_label(end) == 'x'
                marker_data_belt.marker_trajectories(:, i_label) = marker_data_belt.marker_trajectories(:, i_label) + belt_translation_resampled(:, 1);
            end
            if ~isempty(this_label) && this_label(end) == 'y'
                marker_data_belt.marker_trajectories(:, i_label) = marker_data_belt.marker_trajectories(:, i_label) + belt_translation_resampled(:, 2);
            end
            if ~isempty(this_label) && this_label(end) == 'z'
                marker_data_belt.marker_trajectories(:, i_label) = marker_data_belt.marker_trajectories(:, i_label) + belt_translation_resampled(:, 3);
            end
        end
        
        % save translated data
        marker_data_belt.markerBelt_trajectories = marker_data_belt.marker_trajectories;
        marker_data_belt = rmfield(marker_data_belt, 'marker_trajectories');
        save_folder = 'processed';
        save_file_name = makeFileName(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'markerBeltTrajectories.mat');
        save([save_folder filesep save_file_name], '-struct', 'marker_data_belt');
        addAvailableData('markerBelt_trajectories', 'time_mocap', 'sampling_rate_mocap', '_marker_labels', '_marker_directions', save_folder, save_file_name);
        
        disp(['processed ' load_file_name ' and translated to belt frame'])
    end
end

