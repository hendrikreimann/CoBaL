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
    if strcmp(csv_data_file_name(last_underscore+1 : end-4), 'forcePlateData')
        
        % this is labview data, fix this later
        
        % import force plate data
        [imported_data, delimiter, nheaderlines] = importdata(csv_data_file_name, ',', 3);
        forceplate_trajectories = imported_data.data;
        
        % extract headers
        column_name_string = imported_data.textdata{1, 1};
        number_of_columns = size(imported_data.textdata, 2);
        
        time_labview = forceplate_trajectories(:, 1) * milliseconds_to_seconds;
        
        fxl_trajectory = forceplate_trajectories(:, 2);
        fyl_trajectory = forceplate_trajectories(:, 3);
        fzl_trajectory = forceplate_trajectories(:, 4);
        mxl_trajectory = forceplate_trajectories(:, 5);
        myl_trajectory = forceplate_trajectories(:, 6);
        mzl_trajectory = forceplate_trajectories(:, 7);
        
        fxr_trajectory = forceplate_trajectories(:, 8);
        fyr_trajectory = forceplate_trajectories(:, 9);
        fzr_trajectory = forceplate_trajectories(:, 10);
        mxr_trajectory = forceplate_trajectories(:, 11);
        myr_trajectory = forceplate_trajectories(:, 12);
        mzr_trajectory = forceplate_trajectories(:, 13);
        
        fx_trajectory = forceplate_trajectories(:, 14);
        fy_trajectory = forceplate_trajectories(:, 15);
        fz_trajectory = forceplate_trajectories(:, 16);
        mx_trajectory = forceplate_trajectories(:, 17);
        my_trajectory = forceplate_trajectories(:, 18);
        mz_trajectory = forceplate_trajectories(:, 19);
        
        copxl_trajectory = forceplate_trajectories(:, 20);
        copyl_trajectory = forceplate_trajectories(:, 21);
        copxr_trajectory = forceplate_trajectories(:, 22);
        copyr_trajectory = forceplate_trajectories(:, 23);
        copx_trajectory = forceplate_trajectories(:, 24);
        copy_trajectory = forceplate_trajectories(:, 25);
        
        belt_speed_left_trajectory = forceplate_trajectories(:, 26);
        belt_speed_right_trajectory = forceplate_trajectories(:, 27);
        
        gvs_out_trajectory = forceplate_trajectories(:, 28);
        gvs_read_trajectory = forceplate_trajectories(:, 29);
        
%         visual_shift_ml_trajectory = forceplate_trajectories(:, 30);
%         left_foot_state = forceplate_trajectories(:, 31);
%         right_foot_state = forceplate_trajectories(:, 32);
%         stimulus_foot_state = forceplate_trajectories(:, 33);
%         heel_strike_count = forceplate_trajectories(:, 34);
        
        % replacement for the five lines above for legacy file format without viusal shift trajectory
        visual_shift_ml_trajectory = [];
        left_foot_state = forceplate_trajectories(:, 30);
        right_foot_state = forceplate_trajectories(:, 31);
        stimulus_foot_state = forceplate_trajectories(:, 32);
        heel_strike_count = forceplate_trajectories(:, 33);
        
        % save
        matlab_data_file_name = makeFileName(date, subject_id, 'walking', trial_number, 'labviewTrajectories');
        
%             'left_forceplate_wrench_Acw', ...
%             'left_forceplate_cop_Acw', ...
%             'right_forceplate_wrench_Acw', ...
%             'right_forceplate_cop_Acw', ...
        save ...
          ( ...
            matlab_data_file_name, ...
            'time_labview', ...
            'fxl_trajectory', ...
            'fxl_trajectory', ...
            'fyl_trajectory', ...
            'fzl_trajectory', ...
            'mxl_trajectory', ...
            'myl_trajectory', ...
            'mzl_trajectory', ...
            'fxr_trajectory', ...
            'fyr_trajectory', ...
            'fzr_trajectory', ...
            'mxr_trajectory', ...
            'myr_trajectory', ...
            'mzr_trajectory', ...
            'fx_trajectory', ...
            'fy_trajectory', ...
            'fz_trajectory', ...
            'mx_trajectory', ...
            'my_trajectory', ...
            'mz_trajectory', ...
            'copxl_trajectory', ...
            'copyl_trajectory', ...
            'copxr_trajectory', ...
            'copyr_trajectory', ...
            'copx_trajectory', ...
            'copy_trajectory', ...
            'belt_speed_left_trajectory', ...
            'belt_speed_right_trajectory', ...
            'gvs_out_trajectory', ...
            'gvs_read_trajectory', ...
            'visual_shift_ml_trajectory', ...
            'left_foot_state', ...
            'right_foot_state', ...
            'stimulus_foot_state', ...
            'heel_strike_count' ...
          );
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
                    elseif strcmp(data_group{2}(1 : 27), 'Bertec Force Plates - Force')
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


