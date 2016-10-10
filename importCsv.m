% transform raw data from .csv into matlab data
transform_forceplate_data_to_Acw = 1;

current_directory = pwd;
data_dir = dir('*.csv');
clear file_name_list;
[file_name_list{1:length(data_dir)}] = deal(data_dir.name);

number_of_files = length(file_name_list);

millimeter_to_meter = 1e-3; % millimeter to meter
milliseconds_to_seconds = 1e-3; % millimeter to meter

for i_file = 1 : number_of_files
    % file name stuff
    csv_data_file_name = file_name_list{i_file};
    
    [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(csv_data_file_name);
    
    % is this marker data or force plate data?
    last_underscore = find(csv_data_file_name == '_', 1, 'last');
    if strcmp(csv_data_file_name(last_underscore+1 : end-4), 'labviewData')
        
        % import labview data
        [imported_data, delimiter, nheaderlines] = importdata(csv_data_file_name, ',', 3);
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
        
        % old code where the extraction and assignment was done by hand
%         variables_to_save.time_labview = labview_trajectories(:, 1) * milliseconds_to_seconds;
%         
%         variables_to_save.fxl_trajectory = labview_trajectories(:, 2);
%         variables_to_save.fyl_trajectory = labview_trajectories(:, 3);
%         variables_to_save.fzl_trajectory = labview_trajectories(:, 4);
%         variables_to_save.mxl_trajectory = labview_trajectories(:, 5);
%         variables_to_save.myl_trajectory = labview_trajectories(:, 6);
%         variables_to_save.mzl_trajectory = labview_trajectories(:, 7);
%         
%         variables_to_save.fxr_trajectory = labview_trajectories(:, 8);
%         variables_to_save.fyr_trajectory = labview_trajectories(:, 9);
%         variables_to_save.fzr_trajectory = labview_trajectories(:, 10);
%         variables_to_save.mxr_trajectory = labview_trajectories(:, 11);
%         variables_to_save.myr_trajectory = labview_trajectories(:, 12);
%         variables_to_save.mzr_trajectory = labview_trajectories(:, 13);
%         
%         variables_to_save.fx_trajectory = labview_trajectories(:, 14);
%         variables_to_save.fy_trajectory = labview_trajectories(:, 15);
%         variables_to_save.fz_trajectory = labview_trajectories(:, 16);
%         variables_to_save.mx_trajectory = labview_trajectories(:, 17);
%         variables_to_save.my_trajectory = labview_trajectories(:, 18);
%         variables_to_save.mz_trajectory = labview_trajectories(:, 19);
%         
%         variables_to_save.copxl_trajectory = labview_trajectories(:, 20);
%         variables_to_save.copyl_trajectory = labview_trajectories(:, 21);
%         variables_to_save.copxr_trajectory = labview_trajectories(:, 22);
%         variables_to_save.copyr_trajectory = labview_trajectories(:, 23);
%         variables_to_save.copx_trajectory = labview_trajectories(:, 24);
%         variables_to_save.copy_trajectory = labview_trajectories(:, 25);
%         
%         variables_to_save.belt_speed_left_trajectory = labview_trajectories(:, 26);
%         variables_to_save.belt_speed_right_trajectory = labview_trajectories(:, 27);
%         
%         
%         if size(labview_trajectories, 2) >= 29
%             variables_to_save.gvs_out_trajectory = labview_trajectories(:, 28);
%             variables_to_save.gvs_read_trajectory = labview_trajectories(:, 29);
%         else
%             variables_to_save.gvs_out_trajectory = zeros(size(variables_to_save.time_labview));
%             variables_to_save.gvs_read_trajectory = zeros(size(variables_to_save.time_labview));
%         end
%         
% %         visual_shift_ml_trajectory = forceplate_trajectories(:, 30);
% %         left_foot_state = forceplate_trajectories(:, 31);
% %         right_foot_state = forceplate_trajectories(:, 32);
% %         stimulus_foot_state = forceplate_trajectories(:, 33);
% %         heel_strike_count = forceplate_trajectories(:, 34);
%         
%         % replacement for the five lines above for legacy file format without visual shift trajectory
%         if size(labview_trajectories, 2) >= 34
%             variables_to_save.visual_shift_ml_trajectory = [];
%             variables_to_save.left_foot_state = labview_trajectories(:, 30);
%             variables_to_save.right_foot_state = labview_trajectories(:, 31);
%             variables_to_save.stimulus_foot_state = labview_trajectories(:, 32);
%             variables_to_save.heel_strike_count = labview_trajectories(:, 33);
%         else
%             variables_to_save.visual_shift_ml_trajectory = zeros(size(variables_to_save.time_labview));
%             variables_to_save.left_foot_state = zeros(size(variables_to_save.time_labview));
%             variables_to_save.right_foot_state = zeros(size(variables_to_save.time_labview));
%             variables_to_save.stimulus_foot_state = zeros(size(variables_to_save.time_labview));
%             variables_to_save.heel_strike_count = zeros(size(variables_to_save.time_labview));
%             
%         end
        
        % save
        matlab_data_file_name = makeFileName(date, subject_id, 'walking', trial_number, 'labviewTrajectories');
        save(matlab_data_file_name, '-struct', 'variables_to_save');
            
    elseif strcmp(csv_data_file_name(last_underscore+1 : end-4), 'armSenseData')
        % ignore for now
    else
        import_more_data = true;
        number_of_header_lines = 5;
        while import_more_data
            % import data
            [imported_data, delimiter, nheaderlines] = importdata(csv_data_file_name, ',', number_of_header_lines);
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
                        error(['data not recognized: file "' csv_data_file_name])
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
                    matlab_data_file_name = [csv_data_file_name(1 : end-4) '_markerTrajectories.mat'];
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
    
    
    
    disp(['imported ' csv_data_file_name])
end
disp(['imported ' num2str(number_of_files) ' files'])

% save([data_root directorySeparator 'markerData_raw.mat'], 'marker_trajectories_raw', 'file_name_list');


