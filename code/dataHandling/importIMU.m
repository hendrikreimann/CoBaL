

%%% Imports .h5 files to CoBal format
%%% 
%%% To run need to rename the files in WalkH5 folder. Must manually look at
%%% walkCSV to see which trial is which. Labeling should follow our usual
%%% naming convention 
    % [date_subjID_trialtype_trialnumber (000X)]
    % run from subject folder...

%%% Sensor Correspondence  %%%
% Right_Wrist = 5995;
% Left_wrist = 5996;
% Right_Foot = 6048;
% Left_Foot = 6049;
% Sternum = 5988;
% Lumbar = 6016;

% Gyro data is in rad/sec
% Commented out the raw trajectories in available Variables for now.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function importIMU (varargin)

parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'APDM', true)
    addParameter(parser, 'Gyroscope', false)
    addParameter(parser, 'Accelerometer', false)
    parse(parser, varargin{:})
   
    APDM = parser.Results.APDM;

   %%% where the .h5 files are
    walk_source_dir = 'WalkH5';

 if parser.Results.APDM
        clear file_name_list;
        data_dir_mat = dir([walk_source_dir filesep '*.h5']);
        [file_name_list{1:length(data_dir_mat)}] = deal(data_dir_mat.name);

        % go through files and import
        number_of_files = length(file_name_list);
        importing_trial_number = 1;      %%% can make this update each time but its just 1 trial per conditon
        for i_file = 1 : number_of_files
    %     for i_file = 3 : number_of_files % XXX to speed things up during bug-fixing
            % file name stuff
            data_file_name = file_name_list{i_file};
            disp(['Importing data from ' data_file_name])
            [date, subject_id, trial_type, trial_number] = getFileParameters(data_file_name);

    
            % create folders if necessary
            if ~directoryExists('raw')
                mkdir('raw')
            end
            if ~directoryExists('processed')
                mkdir('processed')
            end
            if ~directoryExists('analysis')
                mkdir('analysis')
            end
            current_path = pwd;
            path_split = strsplit(current_path, filesep);
            subject_code = path_split{end};

            %%% Sensor Info
       
            %%% Displays all displays the entire HDF5 file's metadata. 
            % h5info(filename);
            % h5disp(filename);
            %%% Stores .h5 file info as file info.            
            %file_info = h5info(filename);

            %%%Sample rate is 128. Found in /sensors/number/configuration
            sensor_sampling_rate = 128;


            %%% Puts data in a struct  
            filename = data_file_name;

            data = struct;
            data.RightWristGyroData = h5read([walk_source_dir filesep filename],'/Sensors/5995/Gyroscope');
            data.RightWristGyroTime = h5read([walk_source_dir filesep filename],'/Sensors/5995/Time');       
            data.LeftWristGyroData = h5read([walk_source_dir filesep filename],'/Sensors/5996/Gyroscope');
            data.LeftWristGyroTime = h5read([walk_source_dir filesep filename],'/Sensors/5996/Time');
            data.RightFootGyroData = h5read([walk_source_dir filesep filename],'/Sensors/6048/Gyroscope');
            data.RightFootGyroTime = h5read([walk_source_dir filesep filename],'/Sensors/6048/Time');
            data.LeftFootGyroData = h5read([walk_source_dir filesep filename],'/Sensors/6049/Gyroscope');
            data.LeftFootGyroTime = h5read([walk_source_dir filesep filename],'/Sensors/6049/Time');
            data.LumbarGyroData = h5read([walk_source_dir filesep filename],'/Sensors/6016/Gyroscope');
            data.LumbarGyroTime = h5read([walk_source_dir filesep filename],'/Sensors/6016/Time');
            data.SternumGyroData = h5read([walk_source_dir filesep filename],'/Sensors/5988/Gyroscope');
            data.SternumGyroTime = h5read([walk_source_dir filesep filename],'/Sensors/5988/Time');

            %%% Saves above as .mat file if needed  %%%
            %save_file_name = '[NeedName]SensorInfo.mat';
             %           save(save_file_name, '-struct', 'data');

            %%% Get sensor trajectories and puts them in same matrix dimensions as time. X Y Z columns  


            right_wrist_gyro_raw_trajectories = data.RightWristGyroData';
            left_wrist_gyro_raw_trajectories = data.LeftWristGyroData';
            right_foot_gyro_raw_trajectories = data.RightFootGyroData';
            left_foot_gyro_raw_trajectories = data.LeftFootGyroData';
            lumbar_gyro_raw_trajectories = data.LumbarGyroData';
            sternum_gyro_raw_trajectories = data.SternumGyroData';



            %%% Time units microseconds. Found in /sensors/number/ dataset 'Time' 
            %%% Get time_sensor  %%%
            microseconds_to_seconds = 1e-6;       

            % Make time into double
            right_wrist_time_sensor = double(data.RightWristGyroTime).* microseconds_to_seconds;
            right_wrist_time_sensor = right_wrist_time_sensor - right_wrist_time_sensor(1);

            left_wrist_time_sensor = double(data.LeftWristGyroTime).* microseconds_to_seconds;
            left_wrist_time_sensor = left_wrist_time_sensor - left_wrist_time_sensor(1);

            right_foot_time_sensor = double(data.RightFootGyroTime).* microseconds_to_seconds;
            right_foot_time_sensor = right_foot_time_sensor - right_foot_time_sensor(1);

            left_foot_time_sensor = double(data.LeftFootGyroTime).* microseconds_to_seconds;
            left_foot_time_sensor = left_foot_time_sensor - left_foot_time_sensor(1);

            lumbar_time_sensor = double(data.LumbarGyroTime).* microseconds_to_seconds;
            lumbar_time_sensor = lumbar_time_sensor - lumbar_time_sensor(1);

            sternum_time_sensor = double(data.SternumGyroTime).* microseconds_to_seconds;
            sternum_time_sensor = sternum_time_sensor - sternum_time_sensor(1);



            %%% Get Sensor labels   %%%
            sensor_labels = {'x', 'y', 'z'};

            %%% Get Sensor_Directions  %%%   
            %%%% NEED TO DOUBLE CHECK THESE %%% Placeholder for now
            sensor_directions = cell(2, length(sensor_labels));
            sensor_directions{1, 1} = 'right - (placeholder)';
            sensor_directions{2, 1} = 'left - (placeholder)';
            sensor_directions{1, 2} = 'forward - (placeholder)';
            sensor_directions{2, 2} = 'backward - (placeholder)';
            sensor_directions{1, 3} = 'up - (placeholder)';
            sensor_directions{2, 3} = 'down - (placeholder)';



            %%% SAVE
            save_this_trial = 1;
            if save_this_trial
                save_folder = 'raw';

                save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'gyroscopeTrajectoriesRaw.mat');

                save ...
                    ( ...
                    [save_folder filesep save_file_name], ...
                    'right_wrist_gyro_raw_trajectories', ...
                    'left_wrist_gyro_raw_trajectories', ...
                    'right_foot_gyro_raw_trajectories', ...
                    'left_foot_gyro_raw_trajectories', ...
                    'lumbar_gyro_raw_trajectories', ...
                    'sternum_gyro_raw_trajectories', ...
                    'right_wrist_time_sensor', ...
                    'left_wrist_time_sensor', ...
                    'right_foot_time_sensor', ...
                    'left_foot_time_sensor', ...
                    'lumbar_time_sensor', ...
                    'sternum_time_sensor', ...
                    'sensor_sampling_rate', ...
                    'sensor_labels', ...
                    'sensor_directions' ...
                    );
                addAvailableData('right_wrist_gyro_raw_trajectories', 'right_wrist_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
                addAvailableData('left_wrist_gyro_raw_trajectories', 'left_wrist_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
                addAvailableData('right_foot_gyro_raw_trajectories', 'right_foot_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
                addAvailableData('left_foot_gyro_raw_trajectories', 'left_foot_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
                addAvailableData('lumbar_gyro_raw_trajectories', 'lumbar_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
                addAvailableData('sternum_gyro_raw_trajectories', 'sternum_time_sensor', 'sensor_sampling_rate', 'sensor_labels', 'sensor_directions', save_folder, save_file_name);
            end
            

            %%% Puts accelerometer data in a struct %%% 
            data = struct;
            data.RightWristAccelData = h5read([walk_source_dir filesep filename],'/Sensors/5995/Accelerometer');
            data.RightWristAccelTime = h5read([walk_source_dir filesep filename],'/Sensors/5995/Time');       
            data.LeftWristAccelData = h5read([walk_source_dir filesep filename],'/Sensors/5996/Accelerometer');
            data.LeftWristAccelTime = h5read([walk_source_dir filesep filename],'/Sensors/5996/Time');
            data.RightFootAccelData = h5read([walk_source_dir filesep filename],'/Sensors/6048/Accelerometer');
            data.RightFootAccelTime = h5read([walk_source_dir filesep filename],'/Sensors/6048/Time');
            data.LeftFootAccelData = h5read([walk_source_dir filesep filename],'/Sensors/6049/Accelerometer');
            data.LeftFootAccelTime = h5read([walk_source_dir filesep filename],'/Sensors/6049/Time');
            data.LumbarAccelData = h5read([walk_source_dir filesep filename],'/Sensors/6016/Accelerometer');
            data.LumbarAccelTime = h5read([walk_source_dir filesep filename],'/Sensors/6016/Time');
            data.SternumAccelData = h5read([walk_source_dir filesep filename],'/Sensors/5988/Accelerometer');
            data.SternumAccelTime = h5read([walk_source_dir filesep filename],'/Sensors/5988/Time');

            %%% Acceleration Trajectories  %%%
            right_wrist_accel_raw_trajectories = data.RightWristAccelData';
            left_wrist_accel_raw_trajectories = data.LeftWristAccelData';
            right_foot_accel_raw_trajectories = data.RightFootAccelData';
            left_foot_accel_raw_trajectories = data.LeftFootAccelData';
            lumbar_accel_raw_trajectories = data.LumbarAccelData';
            sternum_accel_raw_trajectories = data.SternumAccelData';

    
            microseconds_to_seconds = 1e-6;       

            % Make time into double
            right_wrist_time_sensor = double(data.RightWristAccelTime).* microseconds_to_seconds;
            right_wrist_time_sensor = right_wrist_time_sensor - right_wrist_time_sensor(1);

            left_wrist_time_sensor = double(data.LeftWristAccelTime).* microseconds_to_seconds;
            left_wrist_time_sensor = left_wrist_time_sensor - left_wrist_time_sensor(1);

            right_foot_time_sensor = double(data.RightFootAccelTime).* microseconds_to_seconds;
            right_foot_time_sensor = right_foot_time_sensor - right_foot_time_sensor(1);

            left_foot_time_sensor = double(data.LeftFootAccelTime).* microseconds_to_seconds;
            left_foot_time_sensor = left_foot_time_sensor - left_foot_time_sensor(1);

            lumbar_time_sensor = double(data.LumbarAccelTime).* microseconds_to_seconds;
            lumbar_time_sensor = lumbar_time_sensor - lumbar_time_sensor(1);

            sternum_time_sensor = double(data.SternumAccelTime).* microseconds_to_seconds;
            sternum_time_sensor = sternum_time_sensor - sternum_time_sensor(1);



            %%% Get Sensor labels   %%%
            sensor_labels = {'x', 'y', 'z'};

            %%% Get Sensor_Directions  %%%   
            sensor_directions = cell(2, length(sensor_labels));
            [sensor_directions{1, [1 4 7 10]}] = deal('right');
            [sensor_directions{2, [1 4 7 10]}] = deal('left');
            [sensor_directions{1, [2 5 8 11]}] = deal('forward');
            [sensor_directions{2, [2 5 8 11]}] = deal('backward');
            [sensor_directions{1, [3 6 9 12]}] = deal('up');
            [sensor_directions{2, [3 6 9 12]}] = deal('down');

            %%% Save 
            save_this_trial = 1;
            if save_this_trial
                save_folder = 'raw';

                save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'accelerometerTrajectoriesRaw.mat');

                save ...
                    ( ...
                    [save_folder filesep save_file_name], ...
                    'right_wrist_accel_raw_trajectories', ...
                    'left_wrist_accel_raw_trajectories', ...
                    'right_foot_accel_raw_trajectories', ...
                    'left_foot_accel_raw_trajectories', ...
                    'lumbar_accel_raw_trajectories', ...
                    'sternum_accel_raw_trajectories', ...
                    'right_wrist_time_sensor', ...
                    'left_wrist_time_sensor', ...
                    'left_foot_time_sensor', ...
                    'right_foot_time_sensor', ...
                    'lumbar_time_sensor', ...
                    'sternum_time_sensor', ...
                    'sensor_sampling_rate', ...
                    'sensor_labels', ...
                    'sensor_directions' ...
                    );
            
                addAvailableData('right_wrist_accel_raw_trajectories', 'right_wrist_time_sensor', 'sensor_sampling_rate', '_sensor_labels', '_sensor_directions', save_folder, save_file_name);
                addAvailableData('left_wrist_accel_raw_trajectories', 'left_wrist_time_sensor', 'sensor_sampling_rate', '_sensor_labels', '_sensor_directions', save_folder, save_file_name);
                addAvailableData('right_foot_accel_raw_trajectories', 'right_foot_time_sensor', 'sensor_sampling_rate', '_sensor_labels', '_sensor_directions', save_folder, save_file_name);
                addAvailableData('left_foot_accel_raw_trajectories', 'left_foot_time_sensor', 'sensor_sampling_rate', '_sensor_labels', '_sensor_directions', save_folder, save_file_name);
                addAvailableData('lumbar_accel_raw_trajectories', 'lumbar_time_sensor', 'sensor_sampling_rate', '_sensor_labels', '_sensor_directions', save_folder, save_file_name);
                addAvailableData('sternum_accel_raw_trajectories', 'sternum_time_sensor', 'sensor_sampling_rate', '_sensor_labels', '_sensor_directions', save_folder, save_file_name);
            end
        end
    end
end
    
            



            