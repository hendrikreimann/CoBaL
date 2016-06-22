do_plots                = 0;
map_sensors_manually    = 1;

current_directory = pwd;
data_dir = dir('*_emgTrajectoriesRaw.mat');
clear file_name_list;
[file_name_list{1:length(data_dir)}] = deal(data_dir.name);

load('subjectInfo.mat');

number_of_files = length(file_name_list);
for i_trial = 1 : number_of_files
    % load data
    raw_emg_file_name = file_name_list{i_trial};
    [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(raw_emg_file_name);
    load(raw_emg_file_name);
    
    % define filters
    % initial low pass filter
    filter_order_low = 4;
    cutoff_frequency_low = 500; % in Hz
    [b_low, a_low] = butter(filter_order_low, cutoff_frequency_low/(sampling_rate_emg/2), 'low');
    % high pass filter at 20 hz to get rid of DC offset
    filter_order_high = 4;
    cutoff_frequency_high = 20; % in Hz
    [b_high, a_high] = butter(filter_order_high, cutoff_frequency_high/(sampling_rate_emg/2), 'high');
    % low pass filter below 10 Hz -- aggressive smoothing after rectification
    filter_order_final = 4;
    cutoff_frequency_final = 10; % in Hz
    [b_final, a_final] = butter(filter_order_final, cutoff_frequency_high/(sampling_rate_emg/2), 'low');
    
    % filter
    emg_trajectories_filtered_lowpass = filtfilt(b_low, a_low, emg_trajectories_raw);
    emg_trajectories_filtered_highpass = filtfilt(b_high, a_high, emg_trajectories_filtered_lowpass);
    
    % rectify signal
    emg_trajectories_rectified = abs(emg_trajectories_filtered_highpass);
    
    % filter again
    emg_trajectories = filtfilt(b_final, a_final, emg_trajectories_rectified);
    
    if map_sensors_manually
        lglutmed_sensor_index = 7;
        ltibiant_sensor_index = 8;
        lperolng_sensor_index = 9;
        rglutmed_sensor_index = 1;
        rtibiant_sensor_index = 2;
        rperolng_sensor_index = 3;
        
        emg_headers{lglutmed_sensor_index} = 'LGLUTMED';
        emg_headers{ltibiant_sensor_index} = 'LTIBIANT';
        emg_headers{lperolng_sensor_index} = 'LPEROLNG';
        emg_headers{rglutmed_sensor_index} = 'RGLUTMED';
        emg_headers{rtibiant_sensor_index} = 'RTIBIANT';
        emg_headers{rperolng_sensor_index} = 'RPEROLNG';
    end
    
    % save
    save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'emgTrajectories');
    save ...
      ( ...
        save_file_name, ...
        'emg_trajectories', ...
        'time_emg', ...
        'sampling_rate_emg', ...
        'emg_headers' ...
      );
    disp(['filtered and saved as ' save_file_name])
    
    % delete raw data file
%     delete(raw_emg_file_name);
    
    % visualize
    if do_plots
        i_channel = 1;
        figure; axes; hold on
        plot(time_emg, emg_trajectories_raw(:, i_channel));
        plot(time_emg, emg_trajectories_filtered_lowpass(:, i_channel));
        plot(time_emg, emg_trajectories_filtered_highpass(:, i_channel));
        plot(time_emg, emg_trajectories_rectified(:, i_channel));
        plot(time_emg, emg_trajectories(:, i_channel), 'linewidth', 2);
    end    
    
    
    
    
    
end