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

% this script transform raw data from ascii into matlab format

% input:
% Experimental data files generated by e.g. Nexus, QTM, labview etc.
% Files are expected to be in subfolders named <subject code>_<source type>, e.g. XYZ_labview
% any source type is viable, but if it is not in the default list, it must be specified as a name-value pair,
% e.g. importAscii('sources', 'someSource') would look for data in the subfolder XYZ_someSource
%
% output:
% Files containing the same data in .mat format, with some additional information about where they came from.
% Output files will be saved to folders "raw" and "processed".

function importQTM(varargin)
parser = inputParser;
parser.KeepUnmatched = true;
sources_default = {'qtm', 'labview'};
addParameter(parser, 'sources', sources_default)
parse(parser, varargin{:})
sources = parser.Results.sources;

%% prepare
% set some parameters
millimeter_to_meter = 1e-3;
centimeter_to_meter = 1e-2;
milliseconds_to_seconds = 1e-3;
qtm_emg_scale = 1;

% initialize
total_number_of_trials_extracted_this_subject = 0;

% create folders if necessary
if ~exist('raw', 'dir')
    mkdir('raw')
end
if ~exist('processed', 'dir')
    mkdir('processed')
end
if ~exist('analysis', 'dir')
    mkdir('analysis')
end
current_path = pwd;
path_split = strsplit(current_path, filesep);
subject_code = path_split{end};

%% import data
for i_source = 1 : length(sources)
    source_dir = sources{i_source};
    if exist(source_dir, 'dir')
        
        % get list of files to import from this directory
        clear file_name_list_csv;
        data_dir_csv = dir([source_dir filesep '*csv']);
        [file_name_list_csv{1:length(data_dir_csv)}] = deal(data_dir_csv.name);
        
        clear file_name_list_mat;
        data_dir_mat = dir([source_dir filesep '*mat']);
        [file_name_list_mat{1:length(data_dir_mat)}] = deal(data_dir_mat.name);
        
        file_name_list = [file_name_list_mat file_name_list_csv];
        
        % go through files and import
        number_of_files = length(file_name_list);
        for i_file = 1 : number_of_files
            % file name stuff
            data_file_name = file_name_list{i_file};
            [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(data_file_name);
            if isempty(file_type)
                if data_file_name(end-2) == 'm'
                    file_type = 'qualisysData';
                else
                    file_type = 'unknown';
                end
            end
            
            %% qualisys
            if strcmp(file_type, 'qualisysData')
                % this is marker data from QTM
                % UNTESTED -- do stuff here (David)
                % note -- need to be able to handle multiple files, each
                %      with multiple trials
                
                data_source = 'qtm';
                varName = whos('-file', [source_dir, filesep, data_file_name]);
                temp_data = load([source_dir, filesep, data_file_name]);
                qtm_data = temp_data.(varName.name);
                
              

                if strcmp(trial_type, 'calibration')
                   % Markers -- only need/can access marker data from
                   % calibration file (no Analog?.. bc labview not running?)
                   
                   data_type = 'markers';
                   % deal with marker data
                   sampling_rate_mocap = qtm_data.FrameRate;
                   marker_labels = qtm_data.Trajectories.Labeled.Labels;
                   sampling_rate_mocap = qtm_data.FrameRate;
                   
                   tempMarkers = qtm_data.Trajectories.Labeled.Data(...
                      :,...
                      1:3,...
                      :) * millimeter_to_meter;
                  
                    marker_count = 1;
                    marker_trajectories_raw = [];
                    for i_marker = 1: size(tempMarkers,1)
                        this_marker = tempMarkers(i_marker,:,:);
                        marker_trajectories_raw(marker_count:marker_count+2,:) = reshape(this_marker, size(this_marker,2), size(this_marker,3)) * millimeter_to_meter; 
                        marker_count = marker_count + 3;
                    end
                    
                    marker_trajectories_raw = marker_trajectories_raw';

                  time_mocap = (1 : length(marker_trajectories_raw)) / sampling_rate_mocap;
                   if isrow(time_mocap)
                       time_mocap = time_mocap';
                   end
                   
                   % make directions
                    % NOTE: this defines directions and makes assumptions, make sure everything is right here
                    number_of_marker_trajectories = size(marker_trajectories_raw, 2);
                    marker_directions = cell(2, number_of_marker_trajectories);
                    [marker_directions{1, 1 : 3 : number_of_marker_trajectories}] = deal('right');
                    [marker_directions{2, 1 : 3 : number_of_marker_trajectories}] = deal('left');
                    [marker_directions{1, 2 : 3 : number_of_marker_trajectories}] = deal('forward');
                    [marker_directions{2, 2 : 3 : number_of_marker_trajectories}] = deal('backward');
                    [marker_directions{1, 3 : 3 : number_of_marker_trajectories}] = deal('up');
                    [marker_directions{2, 3 : 3 : number_of_marker_trajectories}] = deal('down');
                    
                  % save
                   save_folder = 'raw';
                   save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'markerTrajectoriesRaw.mat');
                   save ...
                       ( ...
                       [save_folder filesep save_file_name], ...
                       'marker_trajectories_raw', ...
                       'time_mocap', ...
                       'data_source', ...
                       'sampling_rate_mocap', ...
                       'marker_labels', ...
                       'marker_directions' ...
                       );
                   addAvailableData_new('marker_trajectories_raw', 'time_mocap', 'sampling_rate_mocap', 'marker_labels','_marker_directions', save_folder, save_file_name);
                   disp(['imported ' source_dir filesep data_file_name ' and saved as ' save_folder filesep save_file_name])
                                      
                else
                    % this is not calibration data, so break up into 2min chunks
                    analog_fs = qtm_data.Analog.Frequency;
                  
                    % Find trial start times

                    trigger_mask = contains(qtm_data.Analog.Labels, 'labview_sync');
                    trigger = qtm_data.Analog.Data(trigger_mask,:);

                    % normalise to the range [0,1] and round
                    trigger = round(trigger/max(trigger));
                    % find edges
                    trigger_edges = [diff(trigger),0];
                    analog_index = 1:length(qtm_data.Analog.Data);
                    % trials start on negative edge
                    start_indices = analog_index(trigger_edges==-1);
                    end_indices = analog_index(trigger_edges ==1);

                    if length(end_indices) ~= length(start_indices)
                        sprintf('ERROR: %d start indicies ~= %d end indices\nUSING 120 s trials\n',...
                            length(start_indices), length(end_indices))
                        end_indices = start_indices + 120*analog_fs;
                    end
                    
                    if end_indices(end) > length(qtm_data.Analog.Data)
                        end_indices(end) = length(qtm_data.Analog.Data);
                    end

                    %                  TO DO: Check number of start indices against number of
                    %                  LV files -- PROBLEM: multiple qtm files

                    %                 loop through trials within current file
                    
                    number_of_trials_in_this_qtm_file = length(start_indices);
                    
                    for i_trial_this_qtm_file = 1 : number_of_trials_in_this_qtm_file
                        
                        importing_trial_number = total_number_of_trials_extracted_this_subject + i_trial_this_qtm_file;

                        % EMG
                        %EMG data are in qtm_data.Analog.Data
                        data_type = 'emg';
                        emg_labels = qtm_data.Analog.Labels(14:29);
                        emg_trajectories_raw = qtm_data.Analog.Data(14:29,start_indices(i_trial_this_qtm_file):end_indices(i_trial_this_qtm_file))';
                        sampling_rate_emg = analog_fs;
                        number_of_samples = end_indices(i_trial_this_qtm_file) - start_indices(i_trial_this_qtm_file) + 1;
                        time_emg = (1 : number_of_samples) / sampling_rate_emg;
                        if isrow(time_emg)
                            time_emg = time_emg';
                        end
                        
                        % make directions
                        emg_directions = cell(2, length(emg_labels));
                        [emg_directions{1, :}] = deal('positive');
                        [emg_directions{2, :}] = deal('negative');
                        
                        % save emg data
                        save_folder = 'raw';
                        save_file_name = makeFileName(date, subject_id, trial_type, importing_trial_number, 'emgTrajectoriesRaw.mat');
                        save ...
                            ( ...
                            [save_folder filesep save_file_name], ...
                            'emg_trajectories_raw', ...
                            'time_emg', ...
                            'sampling_rate_emg', ...
                            'data_source', ...
                            'emg_labels', ...
                            'emg_directions' ...
                            );
                        addAvailableData_new('emg_trajectories_raw', 'time_emg', 'sampling_rate_emg', 'emg_labels', '_emg_directions', save_folder, save_file_name);



                        % Force data
                        %Force data are in qtm_data.Force(n).Force
                        data_type = 'forceplate';
                        if any(qtm_data.Force(1).Force)
                        forceplate_tajectories_Left = ...
	                          [ ...
    	                        qtm_data.Force(1).Force(:, start_indices(i_trial_this_qtm_file) : end_indices(i_trial_this_qtm_file))', ...
        	                    qtm_data.Force(1).Moment(:, start_indices(i_trial_this_qtm_file) : end_indices(i_trial_this_qtm_file))' ...
            	              ];
                	        forceplate_tajectories_Right = ...
                    	      [ ...
                        	    qtm_data.Force(2).Force(:, start_indices(i_trial_this_qtm_file) : end_indices(i_trial_this_qtm_file))', ...
                            	qtm_data.Force(2).Moment(:, start_indices(i_trial_this_qtm_file) : end_indices(i_trial_this_qtm_file))' ...
                          	  ];
	                        forceplate_trajectories_raw = [forceplate_tajectories_Left, forceplate_tajectories_Right];
                        else % currently taking volts... need to scale accordinginly
                            forceplate_trajectories_raw = [qtm_data.Analog.Data(1:12,start_indices(i_trial_this_qtm_file):end_indices(i_trial_this_qtm_file))]';
                        end
                        forceplate_labels = qtm_data.Analog.Labels(1:12);
                        forceplate_location_Acl = qtm_data.Force(1).ForcePlateLocation(4,:)/1000; % back left corner coordinates (m)
                        forceplate_location_Acr = qtm_data.Force(2).ForcePlateLocation(3,:)/1000; % back right corner coordinates (m)
                        sampling_rate_forceplate = analog_fs;
                        time_forceplate = (1 : number_of_samples) / sampling_rate_forceplate;
                        if isrow(time_forceplate)
                            time_forceplate = time_forceplate';
                        end
                        
                        % make directions
                        % NOTE: this defines directions and makes assumptions, make sure everything is right here
                        forceplate_directions = cell(2, length(forceplate_labels));
                        [forceplate_directions{1, [1 4 7 10]}] = deal('right');
                        [forceplate_directions{2, [1 4 7 10]}] = deal('left');
                        [forceplate_directions{1, [2 5 8 11]}] = deal('forward');
                        [forceplate_directions{2, [2 5 8 11]}] = deal('backward');
                        [forceplate_directions{1, [3 6 9 12]}] = deal('up');
                        [forceplate_directions{2, [3 6 9 12]}] = deal('down');
                        [forceplate_directions{1, [13 15 17]}] = deal('right');
                        [forceplate_directions{1, [14 16 18]}] = deal('left');

                        % save forceplate data
                        save_folder = 'raw';
                        save_file_name = makeFileName(date, subject_id, trial_type, importing_trial_number, 'forceplateTrajectoriesRaw.mat');
                        save ...
                            ( ...
                            [save_folder filesep save_file_name], ...
                            'forceplate_trajectories_raw', ...
                            'forceplate_labels', ...
                            'data_source', ...
                            'time_forceplate', ...
                            'sampling_rate_forceplate', ...
                            'forceplate_location_Acl', ...    
                            'forceplate_location_Acr', ...    
                            'forceplate_directions' ...
                            );
                        addAvailableData_new('forceplate_trajectories_raw', 'time_forceplate', 'sampling_rate_forceplate', 'forceplate_labels', '_forceplate_directions', save_folder, save_file_name);
                        
                       
                        % Markers
                        sampling_rate_mocap = qtm_data.FrameRate;
                        % align indices
                        start_indices_mocap = round(start_indices * ...
                            sampling_rate_mocap/analog_fs);
                        end_indices_mocap = round(end_indices * ...
                            sampling_rate_mocap/analog_fs);
                        number_of_frames = end_indices_mocap(i_trial_this_qtm_file) - start_indices_mocap(i_trial_this_qtm_file) + 1;
                        tempMarkers = qtm_data.Trajectories.Labeled.Data(...
                            :,...
                            1:3,...
                            start_indices_mocap(i_trial_this_qtm_file):end_indices_mocap(i_trial_this_qtm_file));
                     
                        data_type = 'markers';
                        % deal with marker data
                        
                        marker_labels = qtm_data.Trajectories.Labeled.Labels;
                        sampling_rate_mocap = qtm_data.FrameRate;
                        time_mocap = (1 : number_of_frames) / sampling_rate_mocap;
                        if isrow(time_mocap)
                           time_mocap = time_mocap';
                        end
                        
                        
                        
                        marker_count = 1;
                        marker_trajectories_raw = [];
                        for i_marker = 1: size(tempMarkers,1)
                            this_marker = tempMarkers(i_marker,:,:);
                            marker_trajectories_raw(marker_count:marker_count+2,:) = reshape(this_marker, size(this_marker,2), size(this_marker,3)) * millimeter_to_meter; 
                            marker_count = marker_count + 3;
                        end
                        marker_trajectories_raw = marker_trajectories_raw';
                         
                        % make directions
                        % NOTE: this defines directions and makes assumptions, make sure everything is right here
                        number_of_marker_trajectories = size(marker_trajectories_raw, 2);
                        marker_directions = cell(2, number_of_marker_trajectories);
                        [marker_directions{1, 1 : 3 : number_of_marker_trajectories}] = deal('right');
                        [marker_directions{2, 1 : 3 : number_of_marker_trajectories}] = deal('left');
                        [marker_directions{1, 2 : 3 : number_of_marker_trajectories}] = deal('forward');
                        [marker_directions{2, 2 : 3 : number_of_marker_trajectories}] = deal('backward');
                        [marker_directions{1, 3 : 3 : number_of_marker_trajectories}] = deal('up');
                        [marker_directions{2, 3 : 3 : number_of_marker_trajectories}] = deal('down');

                        
                        % save
                        save_folder = 'raw';
                        save_file_name = makeFileName(date, subject_id, trial_type, importing_trial_number, 'markerTrajectoriesRaw.mat');
                        save ...
                            ( ...
                            [save_folder filesep save_file_name], ...
                            'marker_trajectories_raw', ...
                            'time_mocap', ...
                            'data_source', ...
                            'sampling_rate_mocap', ...
                            'marker_labels', ...
                            'marker_directions' ...
                            );
                        addAvailableData_new('marker_trajectories_raw', 'time_mocap', 'sampling_rate_mocap', 'marker_labels', '_marker_directions', save_folder, save_file_name);

                        disp(['imported ' source_dir filesep data_file_name ' and saved as ' save_folder filesep save_file_name])
                    end
                    total_number_of_trials_extracted_this_subject = total_number_of_trials_extracted_this_subject + number_of_trials_in_this_qtm_file;
                end
                %% unknown
            elseif strcmp(file_type, 'unknown')
                disp(['FAILED to import ' data_file_name])
                
                %% labview
            else
                % assume this is labview data
                [imported_data, delimiter, nheaderlines] = importdata([source_dir filesep data_file_name], ',', 3);
                labview_trajectories = imported_data.data;
                
                % extract headers
                column_name_string = imported_data.textdata{1, 1};
                labview_header = strsplit(column_name_string, ',');
                number_of_data_columns = size(imported_data.textdata, 2);
                
                % extract data into properly named variables
                variables_to_save = struct();
                variables_to_save_list = {};
                for i_column = 1 : number_of_data_columns
                    variable_name = [strrep(labview_header{i_column}, ' ', '_'), '_trajectory'];
                    extract_string = ['variables_to_save.' variable_name ' = labview_trajectories(:, i_column);'];
                    eval(extract_string);
                    variables_to_save_list = [variables_to_save_list; variable_name];
                end
                
                % electrodes were inverted for subject STD (red was left, should be right), so correct this
                if strcmp(subject_code, 'STD')
                    figure; hold on;
                    plot(variables_to_save.GVS_out_trajectory);
                    variables_to_save.GVS_out_trajectory = -variables_to_save.GVS_out_trajectory;
                    plot(variables_to_save.GVS_out_trajectory);
                end
                
                % take special care of time, transform to seconds and rename according to file type
                if isrow(variables_to_save.time_trajectory)
                    variables_to_save.time_trajectory = variables_to_save.time_trajectory';
                end
                eval(['variables_to_save.time = variables_to_save.time_trajectory * milliseconds_to_seconds;']);
                variables_to_save.time = variables_to_save.time - variables_to_save.time(1);
                variables_to_save = rmfield(variables_to_save, 'time_trajectory');
                variables_to_save_list(strcmp(variables_to_save_list, 'time_trajectory')) = [];
                
                % add data source
                variables_to_save.data_source = 'labview'; % not tested yet
                
                % add sampling rate
                variables_to_save.sampling_rate = NaN;
                
                % save
                save_folder = 'processed';
                save_file_name = makeFileName(date, subject_id, trial_type, trial_number, file_type);
                save([save_folder filesep save_file_name], '-struct', 'variables_to_save');
                
                for i_variable = 1 : length(variables_to_save_list)
                    if ~checkDataAvailability(date, subject_id, trial_type, trial_number, variables_to_save_list{i_variable})
                        addAvailableData(variables_to_save_list{i_variable}, 'time', 'sampling_rate', variables_to_save_list{i_variable}, save_folder, save_file_name);
                    end
                end
                
                disp(['imported ' source_dir filesep data_file_name ' and saved as ' save_folder filesep save_file_name])
            end
            
        end
    end
    
end