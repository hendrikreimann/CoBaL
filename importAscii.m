% transform raw data from .csv or .tsv into matlab data

% set parameters
millimeter_to_meter = 1e-3; % millimeter to meter
milliseconds_to_seconds = 1e-3; % millimeter to meter

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
potential_sources = {'devices', 'labview', 'marker', 'ascii'};
% potential_sources = {'ascii'};
% potential_sources = {'devices'};
% potential_sources = {'markers'};
for i_source = 1 : length(potential_sources)
    source_dir = [subject_code '_' potential_sources{i_source}];
    if exist(source_dir, 'dir')
        % get list of files to import from this directory
        clear file_name_list;
        data_dir = dir([source_dir filesep '*sv']);
        [file_name_list{1:length(data_dir)}] = deal(data_dir.name);

        % go through files and import
        number_of_files = length(file_name_list);
        for i_file = 1 : number_of_files
            % file name stuff
            data_file_name = file_name_list{i_file};
            [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(data_file_name);

            % 
            if isempty(file_type)
                if data_file_name(end-2) == 't';
                    file_type = 'qualisys';
                elseif data_file_name(end-2) == 'c';
                    file_type = 'nexus';
                else
                    file_type = 'unknown';
                end
            end

            if strcmp(file_type, 'nexus')
                import_more_data = true;
                number_of_header_lines = 5;
                while import_more_data
                    % import data
                    [imported_data, delimiter, nheaderlines] = importdata([source_dir filesep data_file_name], ',', number_of_header_lines);
                    if isstruct(imported_data)
                        % extract info
                        data_class = imported_data.textdata{number_of_header_lines-4, 1};
                        data_group = strsplit(imported_data.textdata{number_of_header_lines-2, 1}, ',');
                        data_headers = strsplit(imported_data.textdata{number_of_header_lines-1, 1}, ',');

                        number_of_samples = size(imported_data.data, 1);

                        if strcmp(data_class, 'Devices')
                            % deal with devices data
                            if strcmp(data_group{2}(1 : 17), 'Delsys Trigno EMG')
                                data_type = 'emg';
                                emg_headers = data_group(2 : end-1);
                                emg_trajectories_raw = imported_data.data(:, 3:end);
                                sampling_rate_emg = str2num(imported_data.textdata{number_of_header_lines-3, 1});
                                time_emg = (1 : number_of_samples) / sampling_rate_emg;

                                % save emg data
                                matlab_data_file_name = ['raw' filesep makeFileName(date, subject_id, trial_type, trial_number, 'emgTrajectoriesRaw.mat')];
                                save ...
                                  ( ...
                                    matlab_data_file_name, ...
                                    'emg_trajectories_raw', ...
                                    'time_emg', ...
                                    'sampling_rate_emg', ...
                                    'emg_headers' ...
                                  );

                            elseif strcmp(data_group{2}(1 : 6), 'Bertec')
                                data_type = 'forceplate';
                                forceplate_trajectories_raw = imported_data.data(:, 3:end);
                                forceplate_headers = data_headers(3 : end);
                                sampling_rate_forceplate = str2num(imported_data.textdata{number_of_header_lines-3, 1});
                                time_forceplate = (1 : number_of_samples) / sampling_rate_forceplate;

                                % save forceplate data
        %                         matlab_data_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'forceplateTrajectoriesRaw');
                                matlab_data_file_name = ['raw' filesep makeFileName(date, subject_id, trial_type, trial_number, 'forceplateTrajectoriesRaw.mat')];
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
        %                     matlab_data_file_name = [data_file_name(1 : end-4) '_markerTrajectoriesRaw.mat'];
                            matlab_data_file_name = ['raw' filesep makeFileName(date, subject_id, trial_type, trial_number, 'markerTrajectoriesRaw.mat')];
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
                disp(['imported ' data_file_name])



            elseif strcmp(file_type, 'qualisys')
                % this is marker data from QTM
                [imported_data, delimiter, nheaderlines] = importdata([source_dir filesep data_file_name], '\t', 10);
                
                data_type = 'markers';
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
                    'sampling_rate_mocap', ...
                    'marker_headers' ...
                  );
                disp(['imported ' data_file_name ' and saved as ' matlab_data_file_name])
            elseif strcmp(file_type, 'a')
                % this is analog data from QTM
                [imported_data, delimiter, nheaderlines] = importdata([source_dir filesep data_file_name], '\t', 14);
                
                data_type = 'emg';
                emg_headers = imported_data.colheaders(3 : end);
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
                    'sampling_rate_emg', ...
                    'emg_headers' ...
                  );
                disp(['imported ' data_file_name ' and saved as ' matlab_data_file_name])
                
                

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

                % take special care of time, transform to seconds and rename according to file type
                eval(['variables_to_save.time_' file_type ' = variables_to_save.time_trajectory * milliseconds_to_seconds;']);
                variables_to_save = rmfield(variables_to_save, 'time_trajectory');

                % save
                matlab_data_file_name = ['processed' filesep makeFileName(date, subject_id, trial_type, trial_number, file_type)];
                save(matlab_data_file_name, '-struct', 'variables_to_save');
                disp(['imported ' data_file_name ' and saved as ' matlab_data_file_name])
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


