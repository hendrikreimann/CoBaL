% transform raw data from .csv into matlab data
transform_force_plate_data_to_Acw = 1;

current_directory = pwd;
data_dir = dir('*.csv');
clear file_name_list;
[file_name_list{1:length(data_dir)}] = deal(data_dir.name);

number_of_files = length(file_name_list);

millimeter_to_meter = 1e-3; % millimeter to meter
milliseconds_to_seconds = 1e-3; % millimeter to meter

for i_trial = 1 : number_of_files
    % file name stuff
    csv_data_file_name = file_name_list{i_trial};
    
    % is this marker data or force plate data?
    last_underscore = find(csv_data_file_name == '_', 1, 'last');
    if strcmp(csv_data_file_name(last_underscore+1 : end-4), 'forcePlateData')
        % import force plate data
        [imported_data, delimiter, nheaderlines] = importdata(csv_data_file_name, ',', 1);
        force_plate_trajectories = imported_data.data;
        
        time_force_plate = force_plate_trajectories(:, 1) * milliseconds_to_seconds;
        
        fxl_trajectory = force_plate_trajectories(:, 2);
        fyl_trajectory = force_plate_trajectories(:, 3);
        fzl_trajectory = force_plate_trajectories(:, 4);
        mxl_trajectory = force_plate_trajectories(:, 5);
        myl_trajectory = force_plate_trajectories(:, 6);
        mzl_trajectory = force_plate_trajectories(:, 7);
        
        fxr_trajectory = force_plate_trajectories(:, 8);
        fyr_trajectory = force_plate_trajectories(:, 9);
        fzr_trajectory = force_plate_trajectories(:, 10);
        mxr_trajectory = force_plate_trajectories(:, 11);
        myr_trajectory = force_plate_trajectories(:, 12);
        mzr_trajectory = force_plate_trajectories(:, 13);
        
        fx_trajectory = force_plate_trajectories(:, 14);
        fy_trajectory = force_plate_trajectories(:, 15);
        fz_trajectory = force_plate_trajectories(:, 16);
        mx_trajectory = force_plate_trajectories(:, 17);
        my_trajectory = force_plate_trajectories(:, 18);
        mz_trajectory = force_plate_trajectories(:, 19);
        
        copxl_trajectory = force_plate_trajectories(:, 20);
        copyl_trajectory = force_plate_trajectories(:, 21);
        copxr_trajectory = force_plate_trajectories(:, 22);
        copyr_trajectory = force_plate_trajectories(:, 23);
        copx_trajectory = force_plate_trajectories(:, 24);
        copy_trajectory = force_plate_trajectories(:, 25);
        
        belt_speed_left_trajectory = force_plate_trajectories(:, 26);
        belt_speed_right_trajectory = force_plate_trajectories(:, 27);
        
        stim_sent_trajectory = force_plate_trajectories(:, 28);
        stim_read_trajectory = force_plate_trajectories(:, 29);
        
        left_foot_state = force_plate_trajectories(:, 30);
        right_foot_state = force_plate_trajectories(:, 31);
        stimulus_foot_state = force_plate_trajectories(:, 32);
        heel_strike_count = force_plate_trajectories(:, 33);
        
        if transform_force_plate_data_to_Acw
            % extract and process data
            left_force_plate_wrench_Acl = [fxl_trajectory fyl_trajectory fzl_trajectory mxl_trajectory myl_trajectory mzl_trajectory];
            left_force_plate_cop_Acl = [copxl_trajectory copyl_trajectory zeros(size(copxl_trajectory))];
            right_force_plate_wrench_Acr = [fxr_trajectory fyr_trajectory fzr_trajectory mxr_trajectory myr_trajectory mzr_trajectory];
            right_force_plate_cop_Acr = [copxr_trajectory copyr_trajectory zeros(size(copxr_trajectory))];

            % define forceplate rotation and translation
            Acr_to_Acw_rotation = [-1 0 0; 0 1 0; 0 0 -1];
            Acr_to_Acw_translation = [0.5588; 0; 0];
            Acr_to_Acw_trafo = [Acr_to_Acw_rotation Acr_to_Acw_translation; 0 0 0 1];
            Acl_to_Acw_rotation = [-1 0 0; 0 1 0; 0 0 -1];
            Acl_to_Acw_translation = [-0.5588; 0; 0];
            Acl_to_Acw_trafo = [Acl_to_Acw_rotation Acl_to_Acw_translation; 0 0 0 1];
            Acr_to_Acw_adjoint = rigidToAdjointTransformation(Acr_to_Acw_trafo);
            Acl_to_Acw_adjoint = rigidToAdjointTransformation(Acl_to_Acw_trafo);

            % transform
            left_force_plate_wrench_Acw = (Acl_to_Acw_adjoint' * left_force_plate_wrench_Acl')';
            left_force_plate_cop_Acw = (eye(2, 4) * Acl_to_Acw_trafo * [left_force_plate_cop_Acl ones(size(left_force_plate_cop_Acl, 1), 1)]')';
            right_force_plate_wrench_Acw = (Acr_to_Acw_adjoint' * right_force_plate_wrench_Acr')';
            right_force_plate_cop_Acw = (eye(2, 4) * Acr_to_Acw_trafo * [right_force_plate_cop_Acr ones(size(right_force_plate_cop_Acr, 1), 1)]')';
            
            % re-zero CoP for low loads
            left_force_plate_low_load_indicator = copxl_trajectory == 0;
            left_force_plate_cop_Acw(left_force_plate_low_load_indicator, :) = 0;
            right_force_plate_low_load_indicator = copxr_trajectory == 0;
            right_force_plate_cop_Acw(right_force_plate_low_load_indicator, :) = 0;
        end        
        
        
        % save
        matlab_data_file_name = [csv_data_file_name(1 : end-4) '.mat'];
        save ...
          ( ...
            matlab_data_file_name, ...
            'time_force_plate', ...
            'left_force_plate_wrench_Acw', ...
            'left_force_plate_cop_Acw', ...
            'right_force_plate_wrench_Acw', ...
            'right_force_plate_cop_Acw', ...
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
            'stim_sent_trajectory', ...
            'stim_read_trajectory', ...
            'left_foot_state', ...
            'right_foot_state', ...
            'stimulus_foot_state', ...
            'heel_strike_count' ...
          );
    else
        % import marker data
        [A, delimiter, nheaderlines] = importdata(csv_data_file_name, ',', 5);
        marker_trajectories_raw = A.data(:, 3:end) * millimeter_to_meter;
        
        sampling_rate = str2num(A.textdata{2, 1});
        time_mocap = A.data(:, 1) * sampling_rate^(-1); % can I find the sampling rate automatically?
        
        % save
        matlab_data_file_name = [csv_data_file_name(1 : end-4) '_rawTrajectories.mat'];
        save(matlab_data_file_name, 'marker_trajectories_raw', 'time_mocap', 'sampling_rate');
    end
    
    
    
    disp(['imported ' csv_data_file_name])
end
disp(['imported ' num2str(number_of_files) ' files'])

% save([data_root directorySeparator 'markerData_raw.mat'], 'marker_trajectories_raw', 'file_name_list');


