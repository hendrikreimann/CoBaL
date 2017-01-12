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

% transform raw data from .csv or .tsv into matlab data

% set parameters
millimeter_to_meter = 1e-3;
centimeter_to_meter = 1e-2;
milliseconds_to_seconds = 1e-3;

% figure out folders
if ~exist('raw', 'dir')
    mkdir('raw')
end
if ~exist('processed', 'dir')
    mkdir('processed')
end
current_path = pwd;
path_split = strsplit(current_path, filesep);
subject_code = path_split{end};

% import data
potential_sources = {'device', 'devices', 'labview', 'marker', 'markers', 'ascii'};
% potential_sources = {'marker', 'markers', 'ascii'};
potential_sources = {'ascii'};
% potential_sources = {'devices'};
% potential_sources = {'markers'};
% potential_sources = {'labview'};
% potential_sources = {'neurocom'};
% potential_sources = {'neurocom', 'ascii'};
for i_source = 1 : length(potential_sources)
    source_dir = [subject_code '_' potential_sources{i_source}];
    if exist(source_dir, 'dir')
        % get list of files to import from this directory
        clear file_name_list_tsv;
        data_dir_tsv = dir([source_dir filesep '*tsv']);
        [file_name_list_tsv{1:length(data_dir_tsv)}] = deal(data_dir_tsv.name);

        clear file_name_list_csv;
        data_dir_csv = dir([source_dir filesep '*csv']);
        [file_name_list_csv{1:length(data_dir_csv)}] = deal(data_dir_csv.name);

        clear file_name_list_txt;
        data_dir_txt = dir([source_dir filesep '*txt']);
        [file_name_list_txt{1:length(data_dir_txt)}] = deal(data_dir_txt.name);
        
        file_name_list = [file_name_list_tsv file_name_list_csv file_name_list_txt];

        % go through files and import
        number_of_files = length(file_name_list);
        for i_file = 1 : number_of_files
            % file name stuff
            data_file_name = file_name_list{i_file};
            [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(data_file_name);

            % 
            if isempty(file_type)
                if data_file_name(end-2) == 't';
                    file_type = 'qualisysData';
                elseif data_file_name(end-2) == 'c';
                    file_type = 'nexus';
                else
                    file_type = 'unknown';
                end
            end

            if strcmp(file_type, 'nexus')
                import_more_data = true;
                number_of_header_lines = 5;
                data_source = 'nexus';
                while import_more_data
                    % import data
                    [imported_data, delimiter, number_of_header_lines_returned] = importdata([source_dir filesep data_file_name], ',', number_of_header_lines);
                    if isstruct(imported_data)
                        % extract info
                        data_class = imported_data.textdata{number_of_header_lines_returned-4, 1};
                        data_group = strsplit(imported_data.textdata{number_of_header_lines_returned-2, 1}, ',');
                        data_headers = strsplit(imported_data.textdata{number_of_header_lines_returned-1, 1}, ',');

                        number_of_samples = size(imported_data.data, 1);

                        if strcmp(data_class, 'Devices')
                            % deal with devices data
                            if strcmp(data_group{2}(1 : 17), 'Delsys Trigno EMG')
                                data_type = 'emg';
                                emg_headers = data_group(2 : end-1);
                                emg_trajectories_raw = imported_data.data(:, 3:end);
                                sampling_rate_emg = str2num(imported_data.textdata{number_of_header_lines_returned-3, 1});
                                time_emg = (1 : number_of_samples) / sampling_rate_emg;
                                % prune header
                                for i_column = 1 : length(emg_headers)
                                    emg_headers{i_column} = strrep(emg_headers{i_column}, 'Delsys Trigno EMG 1.2 - Sensor ', 'EMG');
                                end

                                % save emg data
                                matlab_data_file_name = ['raw' filesep makeFileName(date, subject_id, trial_type, trial_number, 'emgTrajectoriesRaw.mat')];
                                save ...
                                  ( ...
                                    matlab_data_file_name, ...
                                    'emg_trajectories_raw', ...
                                    'time_emg', ...
                                    'sampling_rate_emg', ...
                                    'data_source', ...
                                    'emg_headers' ...
                                  );

                            elseif strcmp(data_group{2}(1 : 6), 'Bertec')
                                data_type = 'forceplate';
                                forceplate_trajectories_raw = imported_data.data(:, 3:end);
                                forceplate_headers = data_headers(3 : end);
                                sampling_rate_forceplate = str2num(imported_data.textdata{number_of_header_lines_returned-3, 1});
                                time_forceplate = (1 : number_of_samples) / sampling_rate_forceplate;

                                % save forceplate data
                                matlab_data_file_name = ['raw' filesep makeFileName(date, subject_id, trial_type, trial_number, 'forceplateTrajectoriesRaw.mat')];
                                save ...
                                  ( ...
                                    matlab_data_file_name, ...
                                    'forceplate_trajectories_raw', ...
                                    'forceplate_headers', ...
                                    'data_source', ...
                                    'time_forceplate', ...
                                    'sampling_rate_forceplate' ...
                                  );
                            else
                                error(['data not recognized: file "' data_file_name])
                            end



                        elseif strcmp(data_class, 'Trajectories')
                            data_type = 'markers';

                            % deal with marker data
                            marker_trajectories = imported_data.data(:, 3:end) * millimeter_to_meter;
                            marker_headers_with_subject = data_group(2 : end-1);
                            sampling_rate_mocap = str2num(imported_data.textdata{number_of_header_lines_returned-3, 1});
                            time_mocap = (1 : number_of_samples) / sampling_rate_mocap;

                            % remove subject name from header strings
                            marker_headers = cell(size(marker_headers_with_subject));
                            for i_marker = 1 : length(marker_headers_with_subject)
                                marker_header = strsplit(marker_headers_with_subject{i_marker}, ':');
                                marker_headers{i_marker} = marker_header{2};
                            end

                            % save
        %                     matlab_data_file_name = [data_file_name(1 : end-4) '_markerTrajectoriesRaw.mat'];
                            matlab_data_file_name = ['raw' filesep makeFileName(date, subject_id, trial_type, trial_number, 'markerTrajectoriesRaw.mat')];
                            save ...
                              ( ...
                                matlab_data_file_name, ...
                                'marker_trajectories', ...
                                'time_mocap', ...
                                'data_source', ...
                                'sampling_rate_mocap', ...
                                'marker_headers' ...
                              );



                        else 
                            error(['unkown data type: ' data_class]); 
                        end

                    else
                        import_more_data = 0;
                    end

                    % prepare for next import
                    number_of_header_lines = number_of_header_lines + number_of_samples + 6;

                end
                disp(['imported ' source_dir filesep data_file_name])

            elseif strcmp(file_type, 'qualisysData')
                % this is marker data from QTM
                [imported_data, delimiter, nheaderlines] = importdata([source_dir filesep data_file_name], '\t', 10);
                
                data_type = 'markers';
                data_source = 'qtm';
                marker_headers_field = imported_data.textdata{10};
                marker_headers_strings = strsplit(marker_headers_field);
                marker_headers = marker_headers_strings(2 : end);
                
                marker_trajectories = imported_data.data(:, 3:end) * millimeter_to_meter;
                sampling_rate_field = imported_data.textdata{4, 1};
                sampling_rate_strings = strsplit(sampling_rate_field);
                sampling_rate_mocap = str2num(sampling_rate_strings{2});
                time_mocap = imported_data.data(:, 2);
                
                % set gaps to NaN instead of 0, which seems to be the incredibly annoying default output of QTM
                marker_trajectories(marker_trajectories==0) = NaN;

                % save marker data
                matlab_data_file_name = ['raw' filesep makeFileName(date, subject_id, trial_type, trial_number, 'markerTrajectoriesRaw.mat')];
                save ...
                  ( ...
                    matlab_data_file_name, ...
                    'marker_trajectories', ...
                    'time_mocap', ...
                    'data_source', ...
                    'sampling_rate_mocap', ...
                    'marker_headers' ...
                  );
                disp(['imported ' source_dir filesep data_file_name ' and saved as ' matlab_data_file_name])
            elseif strcmp(file_type, 'a')
                % this is analog data from QTM
                [imported_data, delimiter, nheaderlines] = importdata([source_dir filesep data_file_name], '\t', 13);
                
                data_type = 'emg';
                data_source = 'qtm';
                emg_headers_field = imported_data.textdata{9};
                emg_headers_strings = strsplit(emg_headers_field);
                emg_headers = emg_headers_strings(2 : end);
                for i_column = 1 : length(emg_headers)
                    emg_headers{i_column} = strrep(emg_headers{i_column}, 'CH', 'EMG');
                end
                
                emg_trajectories_raw = imported_data.data(:, 3:end);
                sampling_rate_field = imported_data.textdata{3, 1};
                sampling_rate_strings = strsplit(sampling_rate_field);
                sampling_rate_emg = str2num(sampling_rate_strings{2});
                time_emg = imported_data.data(:, 2);

                % save emg data
                matlab_data_file_name = ['raw' filesep makeFileName(date, subject_id, trial_type, trial_number, 'emgTrajectoriesRaw.mat')];
                save ...
                  ( ...
                    matlab_data_file_name, ...
                    'emg_trajectories_raw', ...
                    'time_emg', ...
                    'data_source', ...
                    'sampling_rate_emg', ...
                    'emg_headers' ...
                  );
                disp(['imported ' source_dir filesep data_file_name ' and saved as ' matlab_data_file_name])
                
            elseif strcmp(file_type, 'neurocomData')
                % this is data from the neurocom forceplate
                data_source = 'neurocom';
                number_of_header_lines = 30;
                number_of_columns = 19;
                sampling_rate_line_number = 22;
                total_time_line_number = 23;
                label_line_number = 30;
                format = repmat('%f', 1, number_of_columns);
                fid = fopen([source_dir filesep data_file_name]);
                
                header = cell(number_of_header_lines, 1);
                for i_line = 1 : number_of_header_lines
                    header{i_line} = fgetl(fid);
                end
                force_plate_data_cell = textscan(fid, format);
                fclose(fid);

                % extract information from header
                line_split = strsplit(header{sampling_rate_line_number}, ' ');
                sampling_rate_forceplate = str2num(line_split{end});
                line_split = strsplit(header{total_time_line_number}, ' ');
                total_time_rate = str2num(line_split{end});

                % extract and format time
                time_forceplate = force_plate_data_cell{1} * 1/sampling_rate_forceplate;
                
                % extract and format data
                forceplate_trajectories_raw = ...
                  [ ...
                    force_plate_data_cell{2}, ...
                    force_plate_data_cell{3}, ...
                    force_plate_data_cell{4}, ...
                    force_plate_data_cell{5} * centimeter_to_meter, ...
                    force_plate_data_cell{6} * centimeter_to_meter, ...
                    force_plate_data_cell{7} * centimeter_to_meter ...
                    force_plate_data_cell{8}, ...
                    force_plate_data_cell{9}, ...
                    force_plate_data_cell{10}, ...
                    force_plate_data_cell{11} * centimeter_to_meter, ...
                    force_plate_data_cell{12} * centimeter_to_meter, ...
                    force_plate_data_cell{13} * centimeter_to_meter ...
                    force_plate_data_cell{14} * centimeter_to_meter, ...
                    force_plate_data_cell{15} * centimeter_to_meter ...
                    force_plate_data_cell{16} * centimeter_to_meter, ...
                    force_plate_data_cell{17} * centimeter_to_meter ...
                    force_plate_data_cell{18} * centimeter_to_meter, ...
                    force_plate_data_cell{19} * centimeter_to_meter ...                    
                  ];
                
                % extract and format header
                line_split = strsplit(header{label_line_number}, ' ');
                forceplate_headers = cell(1, 18);
                for i_col = 2 : 19
                    cell_split = strsplit(line_split{i_col}, '.');
                    forceplate_headers{i_col-1} = cell_split{2};
                end
                
                % save forceplate data
                matlab_data_file_name = ['raw' filesep makeFileName(date, subject_id, trial_type, trial_number, 'forceplateTrajectoriesRaw.mat')];
                save ...
                  ( ...
                    matlab_data_file_name, ...
                    'forceplate_trajectories_raw', ...
                    'forceplate_headers', ...
                    'time_forceplate', ...
                    'data_source', ...
                    'sampling_rate_forceplate' ...
                  );
                disp(['imported ' source_dir filesep data_file_name])
                
                
                
                
                
                
            elseif strcmp(file_type, 'unknown')
                disp(['FAILED to import ' data_file_name])
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
                for i_column = 1 : number_of_data_columns
                    variable_name = [strrep(labview_header{i_column}, ' ', '_'), '_trajectory'];
                    extract_string = ['variables_to_save.' variable_name ' = labview_trajectories(:, i_column);'];
                    eval(extract_string);
                end
                
                % electrodes were inverted for subject STD (red was left, should be right), so correct this
                if strcmp(subject_code, 'STD')
                    figure; hold on;
                    plot(variables_to_save.GVS_out_trajectory);
                    variables_to_save.GVS_out_trajectory = -variables_to_save.GVS_out_trajectory;
                    plot(variables_to_save.GVS_out_trajectory);
                end
                

                % take special care of time, transform to seconds and rename according to file type
                eval(['variables_to_save.time_' file_type ' = variables_to_save.time_trajectory * milliseconds_to_seconds;']);
                variables_to_save = rmfield(variables_to_save, 'time_trajectory');
                
                % add data source
                variables_to_save.data_source = 'labview'; % not tested yet

                % save
                matlab_data_file_name = ['processed' filesep makeFileName(date, subject_id, trial_type, trial_number, file_type)];
                save(matlab_data_file_name, '-struct', 'variables_to_save');
                disp(['imported ' source_dir filesep data_file_name ' and saved as ' matlab_data_file_name])
            end


        end

        disp(['imported ' num2str(number_of_files) ' files'])        
        
        
        
        
        
        
    end
    
end





%%

% clear file_name_list;
% data_dir = dir('*.mat');
% [file_name_list{1:length(data_dir)}] = deal(data_dir.name);
% for i_file = 1 : length(file_name_list)
%     % file name stuff
%     data_file_name = file_name_list{i_file};
%     
%     movefile(data_file_name, ['../raw/' data_file_name]);
% end
% 
% disp(['moved imported files to ..' filesep 'raw'])













if false
    
    % is this marker data or force plate data?
    last_underscore = find(data_file_name == '_', 1, 'last');
    if strcmp(data_file_name(last_underscore+1 : end-4), 'labviewData')
        
        % import labview data
        [imported_data, delimiter, nheaderlines] = importdata(data_file_name, ',', 3);
        labview_trajectories = imported_data.data;
        
        % extract headers
        column_name_string = imported_data.textdata{1, 1};
        labview_header = strsplit(column_name_string, ',');
        number_of_data_columns = size(imported_data.textdata, 2);
        
        % extract data into properly named variables
        variables_to_save = struct();
        for i_column = 1 : number_of_data_columns
            variable_name = [strrep(labview_header{i_column}, ' ', '_'), '_trajectory'];
            extract_string = ['variables_to_save.' variable_name ' = labview_trajectories(:, i_column);'];
            eval(extract_string);
        end
        
        % take special care of time
        variables_to_save.time_labview = variables_to_save.time_trajectory * milliseconds_to_seconds;
        variables_to_save = rmfield(variables_to_save, 'time_trajectory');
        

        % save
        matlab_data_file_name = makeFileName(date, subject_id, 'walking', trial_number, 'labviewTrajectories');
        save(matlab_data_file_name, '-struct', 'variables_to_save');
            
    elseif strcmp(data_file_name(last_underscore+1 : end-4), 'armSenseData')
        % ignore for now
    else
        import_more_data = true;
        number_of_header_lines = 5;
        while import_more_data
            % import data
            [imported_data, delimiter, nheaderlines] = importdata(data_file_name, ',', number_of_header_lines);
            if isstruct(imported_data)
                % extract info
                data_class = imported_data.textdata{number_of_header_lines-4, 1};
                data_group = strsplit(imported_data.textdata{number_of_header_lines-2, 1}, ',');
                data_headers = strsplit(imported_data.textdata{number_of_header_lines-1, 1}, ',');

                number_of_samples = size(imported_data.data, 1);

                if strcmp(data_class, 'Devices')
                    % deal with devices data
                    if strcmp(data_group{2}(1 : 17), 'Delsys Trigno EMG')
%                     if strcmp(data_group{2}(1 : 26), 'Imported Delsys Trigno IMU')
                        data_type = 'emg';
                        emg_headers = data_group(2 : end-1);
                        emg_trajectories_raw = imported_data.data(:, 3:end);
                        sampling_rate_emg = str2num(imported_data.textdata{number_of_header_lines-3, 1});
                        time_emg = (1 : number_of_samples) / sampling_rate_emg;

                        % save emg data
                        matlab_data_file_name = makeFileName(date, subject_id, 'walking', trial_number, 'emgTrajectoriesRaw');
                        save ...
                          ( ...
                            matlab_data_file_name, ...
                            'emg_trajectories_raw', ...
                            'time_emg', ...
                            'sampling_rate_emg', ...
                            'emg_headers' ...
                          );
                    elseif strcmp(data_group{2}(1 : 26), 'Bertec Force Plate - Force')
                        data_type = 'forceplate';
                        forceplate_trajectories_raw = imported_data.data(:, 3:end);
                        forceplate_headers = data_headers(3 : end);
                        sampling_rate_forceplate = str2num(imported_data.textdata{number_of_header_lines-3, 1});
                        time_forceplate = (1 : number_of_samples) / sampling_rate_forceplate;

                        % save forceplate data
                        matlab_data_file_name = makeFileName(date, subject_id, 'walking', trial_number, 'forceplateTrajectoriesRaw');
                        save ...
                          ( ...
                            matlab_data_file_name, ...
                            'forceplate_trajectories_raw', ...
                            'forceplate_headers', ...
                            'time_forceplate', ...
                            'sampling_rate_forceplate' ...
                          );
                    else
                        error(['data not recognized: file "' data_file_name])
                    end

                    

                elseif strcmp(data_class, 'Trajectories')
                    data_type = 'markers';
                    
                    % deal with marker data
                    marker_trajectories = imported_data.data(:, 3:end) * millimeter_to_meter;
                    marker_headers_with_subject = data_group(2 : end-1);
                    sampling_rate_mocap = str2num(imported_data.textdata{number_of_header_lines-3, 1});
                    time_mocap = (1 : number_of_samples) / sampling_rate_mocap;
                    
                    % remove subject name from header strings
                    marker_headers = cell(size(marker_headers_with_subject));
                    for i_marker = 1 : length(marker_headers_with_subject)
                        marker_header = strsplit(marker_headers_with_subject{i_marker}, ':');
                        marker_headers{i_marker} = marker_header{2};
                    end

                    % save
                    matlab_data_file_name = [data_file_name(1 : end-4) '_markerTrajectories.mat'];
                    save ...
                      ( ...
                        matlab_data_file_name, ...
                        'marker_trajectories', ...
                        'time_mocap', ...
                        'sampling_rate_mocap', ...
                        'marker_headers' ...
                      );
                    
                else 
                    error(['unkown data type: ' data_class]); 
                end
                
            else
                import_more_data = 0;
            end

            % prepare for next import
            number_of_header_lines = number_of_header_lines + number_of_samples + 6;

        end
    end
    
    
    
    disp(['imported ' data_file_name])
end


% save([data_root directorySeparator 'markerData_raw.mat'], 'marker_trajectories_raw', 'file_name_list');


