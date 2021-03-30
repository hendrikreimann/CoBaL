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

% this function applies several basic processing steps to experimental data, e.g. filtering

% input: 
% Experimental data files generated by importAscii.m, in the subfolder "raw"
%
% output: 
% multiple files with processed data for each trial, in the subfolder "processed"


function preprocessMarkerData(varargin)
    % parse arguments
    [types_to_analyze, trials_to_analyze, types_to_exclude, trials_to_exclude] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'type', 'all')
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');

    % add excluded trials back in, because while we don't want to analyze them, we still want to pre-process them
    types_to_analyze = [types_to_analyze; types_to_exclude];
    trials_to_analyze = [trials_to_analyze; trials_to_exclude];
    
    % load static reference trial
    loaded_data = load(['raw' filesep makeFileName(collection_date, subject_id, subject_settings.get('static_reference_trial_type'), subject_settings.get('static_reference_trial_number'), 'markerTrajectoriesRaw.mat')]);
    marker_labels_reference = loaded_data.marker_labels;
    marker_directions_reference = loaded_data.marker_directions;
    marker_reference_trajectories = loaded_data.marker_raw_trajectories;
    i_time = 1;
    while any(isnan(marker_reference_trajectories(i_time, :)))
        i_time = i_time + 1;
    end
    marker_reference = marker_reference_trajectories(i_time, :);

    % load rigid body fill table
    marker_fill_table = subject_settings.getTable('marker_fill_table');

    data_dir = dir(['raw' filesep '*_markerTrajectoriesRaw.mat']);
    clear file_name_list;
    [file_name_list{1:length(data_dir)}] = deal(data_dir.name);
    number_of_files = length(file_name_list);
    for i_trial = 1 : number_of_files
        % load data
        raw_marker_file_name = file_name_list{i_trial};
        [date, subject_id, trial_type, trial_number] = getFileParameters(raw_marker_file_name);
        % does the caller want to process this file?
        if any(strcmp(trial_type, types_to_analyze))
            % condition is set to be processed, now check trial number
            trial_number_list_this_condition = trials_to_analyze{strcmp(trial_type, types_to_analyze)};
            if ismember(trial_number, trial_number_list_this_condition)
                % load data
                loaded_data = load(['raw' filesep raw_marker_file_name]);
                time_mocap = loaded_data.time_mocap;
                sampling_rate_mocap = loaded_data.sampling_rate_mocap;
                marker_labels = loaded_data.marker_labels;
                marker_directions = loaded_data.marker_directions;
                
                % figure out variable
                if isfield(loaded_data, 'marker_trajectories_raw')
                    raw_marker_trajectories = loaded_data.marker_trajectories_raw;
                elseif isfield(loaded_data, 'marker_raw_trajectories')
                    raw_marker_trajectories = loaded_data.marker_raw_trajectories;
                end

                if study_settings.get('filter_marker_data')
                    filter_order = 4;
                    cutoff_frequency = study_settings.get('marker_data_cutoff_frequency'); % in Hz
                    [b_marker, a_marker] = butter(filter_order, cutoff_frequency/(loaded_data.sampling_rate_mocap/2));
                    marker_trajectories = nanfiltfilt(b_marker, a_marker, raw_marker_trajectories);
                else
                    marker_trajectories = raw_marker_trajectories;
                end
                
                % perform rigid body fill if requested in subjectSettings.txt for this trial
                trial_type_match = strcmp(marker_fill_table.trial_type, trial_type);
                trial_number_match = strcmp(marker_fill_table.trial_number, num2str(trial_number));
                trial_match_indices = find(trial_type_match & trial_number_match);
                
                for i_index = 1 : length(trial_match_indices)
                    marker_to_fill = marker_fill_table.marker_to_fill{i_index};
                    source_marker_1 = marker_fill_table.marker_source_1{i_index};
                    source_marker_2 = marker_fill_table.marker_source_2{i_index};
                    source_marker_3 = marker_fill_table.marker_source_3{i_index};
                    
                    marker_trajectories = rigidBodyFillMarkerFromReference ...
                      ( ...
                        marker_trajectories, ...
                        marker_labels, ...
                        marker_directions, ...
                        marker_reference, ...
                        marker_labels_reference, ...
                        marker_directions_reference, ...
                        marker_to_fill, ...
                        source_marker_1, ...
                        source_marker_2, ...
                        source_marker_3, 1 ...
                      );
                end

                % compare marker labels to reference trial
                marker_labels_equal = 0;
                if length(loaded_data.marker_labels) == length(marker_labels_reference)
                    marker_labels_equal = 1;
                    % length is the same, now compare individual labels
                    for i_label = 1 : length(marker_labels_reference)
                        if ~strcmp(marker_labels{i_label}, marker_labels_reference{i_label})
                            marker_labels_equal = 0;
                        end
                    end
                end
                
                % re-order markers to match reference if necessary
                if ~marker_labels_equal
                    marker_trajectories_unsorted = marker_trajectories;
                    marker_labels_unsorted = loaded_data.marker_labels;
                    marker_trajectories = zeros(size(marker_trajectories, 1), length(marker_labels_reference)) * NaN;
                    marker_labels = marker_labels_reference;
                    marker_directions = marker_directions_reference;

                    for i_label = 1 : length(marker_labels_reference)
                        this_label = marker_labels{i_label};
                        this_label_index = find(strcmp(marker_labels_unsorted, this_label));
                        if ~isempty(this_label_index)
                            marker_trajectories(:, i_label) = marker_trajectories_unsorted(:, this_label_index);
                        end
                    end
                end

                % save
                save_folder = 'processed';
                save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'markerTrajectories.mat');
                save ...
                  ( ...
                    [save_folder filesep save_file_name], ...
                    'marker_trajectories', ...
                    'time_mocap', ...
                    'sampling_rate_mocap', ...
                    'marker_labels',  ...
                    'marker_directions' ...
                  );
                addAvailableData ...
                  ( ...
                    'marker_trajectories', ...
                    'time_mocap', ...
                    'sampling_rate_mocap', ...
                    '_marker_labels', ...
                    '_marker_directions', ...
                    save_folder, ...
                    save_file_name ...
                  );
%                     addAvailableData('marker_trajectories', 'time_mocap', 'sampling_rate_mocap', 'marker_labels', save_folder, save_file_name);
                disp(['processed ' raw_marker_file_name ' and saved as ' save_file_name])
            end
        end
    end
end
    
    