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


function preprocessRawData(varargin)
    % parse arguments
    [types_to_analyze, trials_to_analyze, types_to_exclude, trials_to_exclude] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'type', 'all')
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    type = parser.Results.type;

    % add excluded trials back in, because while we don't want to analyze them, we still want to pre-process them
    types_to_analyze = [types_to_analyze; types_to_exclude];
    trials_to_analyze = [trials_to_analyze; trials_to_exclude];
    
    % load settings
    study_settings = loadSettingsFromFile('study');
    
    %% emg
    if strcmp(type, 'emg') || strcmp(type, 'all')
        preprocessEmgData(varargin{:});
    end 
    
    %% forceplate data
    if strcmp(type, 'forceplate') || strcmp(type, 'all')
        preprocessForceplateData(varargin{:});
    end
    
    %% marker data
    if strcmp(type, 'marker') || strcmp(type, 'all')
        preprocessMarkerData(varargin{:});
    end

    %% treadmill marker data
    if strcmp(type, 'treadmill') || strcmp(type, 'all')
        
        data_dir = dir(['raw' filesep '*_treadmillTrajectoriesRaw.mat']);
        clear file_name_list;
        [file_name_list{1:length(data_dir)}] = deal(data_dir.name);
        number_of_files = length(file_name_list);
        for i_trial = 1 : number_of_files
            % load data
            raw_treadmill_file_name = file_name_list{i_trial};
            [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(raw_treadmill_file_name);
            % does the caller want to process this file?
            if any(strcmp(trial_type, types_to_analyze))
                % condition is set to be processed, now check trial number
                trial_number_list_this_condition = trials_to_analyze{strcmp(trial_type, types_to_analyze)};
                if ismember(trial_number, trial_number_list_this_condition)
                    load(['raw' filesep raw_treadmill_file_name]);

                    if study_settings.get('filter_marker_data')
                        filter_order = 4;
                        cutoff_frequency = study_settings.get('marker_data_cutoff_frequency'); % in Hz
                        [b_treadmill, a_treadmill] = butter(filter_order, cutoff_frequency/(sampling_rate_mocap/2));
                        treadmill_trajectories = nanfiltfilt(b_treadmill, a_treadmill, treadmill_trajectories_raw);
                    else
                        treadmill_trajectories = treadmill_trajectories_raw;
                    end
                    

                    % save
                    save_folder = 'processed';
                    save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'treadmillTrajectories.mat');
                    save ...
                      ( ...
                        [save_folder filesep save_file_name], ...
                        'treadmill_trajectories', ...
                        'time_mocap', ...
                        'sampling_rate_mocap', ...
                        'treadmill_labels',  ...
                        'treadmill_directions' ...
                      );
                    addAvailableData ...
                      ( ...
                        'treadmill_trajectories', ...
                        'time_mocap', ...
                        'sampling_rate_mocap', ...
                        '_treadmill_labels', ...
                        '_treadmill_directions', ...
                        save_folder, ...
                        save_file_name ...
                      );
%                     addAvailableData('treadmill_trajectories', 'time_mocap', 'sampling_rate_mocap', 'treadmill_labels', save_folder, save_file_name);
                    disp(['processed ' raw_treadmill_file_name ' and saved as ' save_file_name])
                end
            end
        end
    
    end
    
    


    
    
    
    
    
    
    
    
end
