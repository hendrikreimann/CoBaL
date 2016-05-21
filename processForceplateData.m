do_plots = 0;

current_directory = pwd;
data_dir = dir('*_forceplateTrajectoriesRaw.mat');
clear file_name_list;
[file_name_list{1:length(data_dir)}] = deal(data_dir.name);

load('subjectInfo.mat');



number_of_files = length(file_name_list);
for i_trial = 1 : number_of_files
% for i_trial = 1
    % load data
    raw_forceplate_file_name = file_name_list{i_trial};
    [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(raw_forceplate_file_name);
    load(raw_forceplate_file_name);
    
    % filter
    filter_order_low = 4;
    cutoff_frequency_low = 50; % in Hz
    [b_lowpass, a_lowpass] = butter(filter_order_low, cutoff_frequency_low/(sampling_rate_forceplate/2), 'low');
    forceplate_trajectories_filtered = filtfilt(b_lowpass, a_lowpass, forceplate_trajectories_raw);
    
    % extract
    fxl_trajectory = forceplate_trajectories_filtered(:, 1);
    fyl_trajectory = forceplate_trajectories_filtered(:, 2);
    fzl_trajectory = forceplate_trajectories_filtered(:, 3);
    mxl_trajectory = forceplate_trajectories_filtered(:, 4);
    myl_trajectory = forceplate_trajectories_filtered(:, 5);
    mzl_trajectory = forceplate_trajectories_filtered(:, 6);
    fxr_trajectory = forceplate_trajectories_filtered(:, 7);
    fyr_trajectory = forceplate_trajectories_filtered(:, 8);
    fzr_trajectory = forceplate_trajectories_filtered(:, 9);
    mxr_trajectory = forceplate_trajectories_filtered(:, 10);
    myr_trajectory = forceplate_trajectories_filtered(:, 11);
    mzr_trajectory = forceplate_trajectories_filtered(:, 12);
    
    % calculate CoP
    copxl_trajectory = - myl_trajectory ./ fzl_trajectory;
    copyl_trajectory = mxl_trajectory ./ fzl_trajectory;
    copxr_trajectory = - myr_trajectory ./ fzr_trajectory;
    copyr_trajectory = mxr_trajectory ./ fzr_trajectory;
    
    % zero CoP for no contact times
    fz_threshold = 60;
    copxl_trajectory(fzl_trajectory < fz_threshold) = 0;
    copyl_trajectory(fzl_trajectory < fz_threshold) = 0;
    copxr_trajectory(fzr_trajectory < fz_threshold) = 0;
    copyr_trajectory(fzr_trajectory < fz_threshold) = 0;
    
    % save
    save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'forceplateTrajectories');
    save ...
      ( ...
        save_file_name, ...
        'forceplate_trajectories', ...
        'fxl_trajectory', ...
        'fyl_trajectory', ...
        'fzl_trajectory', ...
        'mxl_trajectory', ...
        'myl_trajectory', ...
        'mzl_trajectory', ...
        'copxl_trajectory', ...
        'copyl_trajectory', ...
        'fxr_trajectory', ...
        'fyr_trajectory', ...
        'fzr_trajectory', ...
        'mxr_trajectory', ...
        'myr_trajectory', ...
        'mzr_trajectory', ...
        'copxr_trajectory', ...
        'copyr_trajectory', ...
        'time_forceplate', ...
        'sampling_rate_forceplate' ...
      );
    disp(['filtered and saved as ' save_file_name])
    
    
    % visualize
    if do_plots
        i_channel = 1;
        figure; axes; hold on
        plot(time_forceplate, forceplate_trajectories_raw(:, i_channel));
        plot(time_forceplate, forceplate_trajectories_filtered_lowpass(:, i_channel));
        plot(time_forceplate, forceplate_trajectories_filtered_highpass(:, i_channel));
        plot(time_forceplate, forceplate_trajectories_rectified(:, i_channel));
        plot(time_forceplate, forceplate_trajectories_filtered(:, i_channel), 'linewidth', 2);
    end    
    
    
    
    
    
end

compareForcePlateData