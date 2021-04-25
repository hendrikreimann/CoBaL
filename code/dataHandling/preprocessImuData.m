function preprocessImuData(varargin)

    parser = inputParser;
    parser.KeepUnmatched = true;
    parse(parser, varargin{:})

    data_source_dir = 'raw';


    clear file_name_list;
    data_dir_mat = dir([data_source_dir filesep '*_gyroscopetrajectoriesraw.mat']);
    [file_name_list{1:length(data_dir_mat)}] = deal(data_dir_mat.name);

    % go through files and import
    number_of_files = length(file_name_list);
    for i_file = 1 : number_of_files
        % file name stuff
        data_file_name = file_name_list{i_file};
        disp(['extracting data from ' data_file_name])
        [date, subject_id, trial_type, trial_number] = getFileParameters(data_file_name);

        % Loading .mat

        loaded_data = load(['raw' filesep data_file_name]);
        data_type = 'gyroscope';

        %% Separate each variable %%

        % LEFT WRIST
        left_wrist_x_gyro = loaded_data.left_wrist_gyro_raw_trajectories(:,1);
        left_wrist_y_gyro = loaded_data.left_wrist_gyro_raw_trajectories(:,2);
        left_wrist_z_gyro = loaded_data.left_wrist_gyro_raw_trajectories(:,3);

        % RIGHT WRIST
        right_wrist_x_gyro = loaded_data.right_wrist_gyro_raw_trajectories(:,1);
        right_wrist_y_gyro = loaded_data.right_wrist_gyro_raw_trajectories(:,2);
        right_wrist_z_gyro = loaded_data.right_wrist_gyro_raw_trajectories(:,3);

        % LEFT FOOT
        left_foot_x_gyro = loaded_data.left_foot_gyro_raw_trajectories(:,1);
        left_foot_y_gyro = loaded_data.left_foot_gyro_raw_trajectories(:,2);
        left_foot_z_gyro = loaded_data.left_foot_gyro_raw_trajectories(:,3);

        % RIGHT FOOT
        right_foot_x_gyro = loaded_data.right_foot_gyro_raw_trajectories(:,1);
        right_foot_y_gyro = loaded_data.right_foot_gyro_raw_trajectories(:,2);
        right_foot_z_gyro = loaded_data.right_foot_gyro_raw_trajectories(:,3);

        % LUMBAR
        lumbar_x_gyro = loaded_data.lumbar_gyro_raw_trajectories(:,1);
        lumbar_y_gyro = loaded_data.lumbar_gyro_raw_trajectories(:,2);
        lumbar_z_gyro = loaded_data.lumbar_gyro_raw_trajectories(:,3);

        % STERNUM
        sternum_x_gyro = loaded_data.sternum_gyro_raw_trajectories(:,1);
        sternum_y_gyro = loaded_data.sternum_gyro_raw_trajectories(:,2);
        sternum_z_gyro = loaded_data.sternum_gyro_raw_trajectories(:,3);

        % auxiliary
        sensor_sampling_rate = loaded_data.sensor_sampling_rate;
        left_wrist_time_sensor = loaded_data.left_wrist_time_sensor;
        right_wrist_time_sensor = loaded_data.right_wrist_time_sensor;
        left_foot_time_sensor = loaded_data.left_foot_time_sensor;
        right_foot_time_sensor = loaded_data.right_foot_time_sensor;
        lumbar_time_sensor = loaded_data.lumbar_time_sensor;
        sternum_time_sensor = loaded_data.sternum_time_sensor;

        % save .mat with left x,y,z and right x,y,z per condition. 
        save_folder = 'processed';
        save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'gyroscopetrajectories.mat');
        save ...
          ( ...
            [save_folder filesep save_file_name], ...
            'left_wrist_x_gyro', ...
            'left_wrist_y_gyro', ...
            'left_wrist_z_gyro', ...
            'right_wrist_x_gyro', ...
            'right_wrist_y_gyro', ...
            'right_wrist_z_gyro',...
            'left_foot_x_gyro',...
            'left_foot_y_gyro',...
            'left_foot_z_gyro',...
            'right_foot_x_gyro',...
            'right_foot_y_gyro',...
            'right_foot_z_gyro',...
            'lumbar_x_gyro',...
            'lumbar_y_gyro', ...
            'lumbar_z_gyro',...
            'sternum_x_gyro',...
            'sternum_y_gyro',...
            'sternum_z_gyro', ...
            'trial_type', ...
            'data_type', ...
            'sensor_sampling_rate', ...
            'left_wrist_time_sensor', ...
            'right_wrist_time_sensor', ...
            'left_foot_time_sensor', ...
            'right_foot_time_sensor', ...
            'lumbar_time_sensor', ...
            'sternum_time_sensor' ...
          );
        addAvailableData('left_wrist_x_gyro', 'left_wrist_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('left_wrist_y_gyro', 'left_wrist_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('left_wrist_z_gyro', 'left_wrist_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('right_wrist_x_gyro', 'right_wrist_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('right_wrist_y_gyro', 'right_wrist_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('right_wrist_z_gyro', 'right_wrist_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('left_foot_x_gyro', 'left_foot_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('left_foot_y_gyro', 'left_foot_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('left_foot_z_gyro', 'left_foot_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('right_foot_x_gyro', 'right_foot_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('right_foot_y_gyro', 'right_foot_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('right_foot_z_gyro', 'right_foot_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('lumbar_x_gyro', 'lumbar_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('lumbar_y_gyro', 'lumbar_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('lumbar_z_gyro', 'lumbar_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('sternum_x_gyro', 'sternum_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('sternum_y_gyro', 'sternum_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('sternum_z_gyro', 'sternum_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);

    end


    % import accelerometer data

    clear file_name_list;
    data_dir_mat = dir([data_source_dir filesep '*_accelerometertrajectoriesraw.mat']);
    [file_name_list{1:length(data_dir_mat)}] = deal(data_dir_mat.name);

    % go through files and import
    number_of_files = length(file_name_list);
    for i_file = 1 : number_of_files
        % file name stuff
        data_file_name = file_name_list{i_file};
        disp(['extracting data from ' data_file_name])
        [date, subject_id, trial_type, trial_number] = getFileParameters(data_file_name);

        % loading .mat
        loaded_data = load(['raw' filesep data_file_name]);
        data_type = 'accelerometer';

        % Separate each variable

        % LEFT WRIST
        left_wrist_x_accel = loaded_data.left_wrist_accel_raw_trajectories(:,1);
        left_wrist_y_accel = loaded_data.left_wrist_accel_raw_trajectories(:,2);
        left_wrist_z_accel = loaded_data.left_wrist_accel_raw_trajectories(:,3);

        % RIGHT WRIST
        right_wrist_x_accel = loaded_data.right_wrist_accel_raw_trajectories(:,1);
        right_wrist_y_accel = loaded_data.right_wrist_accel_raw_trajectories(:,2);
        right_wrist_z_accel = loaded_data.right_wrist_accel_raw_trajectories(:,3);

        % LEFT FOOT
        left_foot_x_accel = loaded_data.left_foot_accel_raw_trajectories(:,1);
        left_foot_y_accel = loaded_data.left_foot_accel_raw_trajectories(:,2);
        left_foot_z_accel = loaded_data.left_foot_accel_raw_trajectories(:,3);

        % RIGHT FOOT
        right_foot_x_accel = loaded_data.right_foot_accel_raw_trajectories(:,1);
        right_foot_y_accel = loaded_data.right_foot_accel_raw_trajectories(:,2);
        right_foot_z_accel = loaded_data.right_foot_accel_raw_trajectories(:,3);

        % LUMBAR
        lumbar_x_accel = loaded_data.lumbar_accel_raw_trajectories(:,1);
        lumbar_y_accel = loaded_data.lumbar_accel_raw_trajectories(:,2);
        lumbar_z_accel = loaded_data.lumbar_accel_raw_trajectories(:,3);

        % STERNUM
        sternum_x_accel = loaded_data.sternum_accel_raw_trajectories(:,1);
        sternum_y_accel = loaded_data.sternum_accel_raw_trajectories(:,2);
        sternum_z_accel = loaded_data.sternum_accel_raw_trajectories(:,3);

        % auxiliary
        sensor_sampling_rate = loaded_data.sensor_sampling_rate;
        left_wrist_time_sensor = loaded_data.left_wrist_time_sensor;
        right_wrist_time_sensor = loaded_data.right_wrist_time_sensor;
        left_foot_time_sensor = loaded_data.left_foot_time_sensor;
        right_foot_time_sensor = loaded_data.right_foot_time_sensor;
        lumbar_time_sensor = loaded_data.lumbar_time_sensor;
        sternum_time_sensor = loaded_data.sternum_time_sensor;

        %%% SAVE  
        % Saves .mat with left x,y,z and right x,y,z per condition. 
        save_folder = 'processed';
        save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'accelerometertrajectories.mat');
        save ...
          ( ...
            [save_folder filesep save_file_name], ...
            'left_wrist_x_accel', ...
            'left_wrist_y_accel', ...
            'left_wrist_z_accel', ...
            'right_wrist_x_accel', ...
            'right_wrist_y_accel', ...
            'right_wrist_z_accel',...
            'left_foot_x_accel',...
            'left_foot_y_accel',...
            'left_foot_z_accel',...
            'right_foot_x_accel',...
            'right_foot_y_accel',...
            'right_foot_z_accel',...
            'lumbar_x_accel',...
            'lumbar_y_accel', ...
            'lumbar_z_accel',...
            'sternum_x_accel',...
            'sternum_y_accel',...
            'sternum_z_accel', ...
            'trial_type', ...
            'data_type', ...
            'sensor_sampling_rate', ...
            'left_wrist_time_sensor', ...
            'right_wrist_time_sensor', ...
            'left_foot_time_sensor', ...
            'right_foot_time_sensor', ...
            'lumbar_time_sensor', ...
            'sternum_time_sensor' ...
          );
        addAvailableData('left_wrist_x_accel', 'left_wrist_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('left_wrist_y_accel', 'left_wrist_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('left_wrist_z_accel', 'left_wrist_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('right_wrist_x_accel', 'right_wrist_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('right_wrist_y_accel', 'right_wrist_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('right_wrist_z_accel', 'right_wrist_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('left_foot_x_accel', 'left_foot_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('left_foot_y_accel', 'left_foot_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('left_foot_z_accel', 'left_foot_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('right_foot_x_accel', 'right_foot_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('right_foot_y_accel', 'right_foot_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('right_foot_z_accel', 'right_foot_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('lumbar_x_accel', 'lumbar_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('lumbar_y_accel', 'lumbar_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('lumbar_z_accel', 'lumbar_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('sternum_x_accel', 'sternum_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('sternum_y_accel', 'sternum_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
        addAvailableData('sternum_z_accel', 'sternum_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);

    end
end
  
  
 