subject = 'AJH';
data_root = '/Users/reimajbi/Drive_UD/GVS';
trial_to_process = 1;
load([data_root filesep subject filesep 'subjectInfo.mat']);
loaded_data = load([data_root filesep subject filesep 'raw' filesep makeFileName(date, subject_id, 'walking', trial_to_process, 'emgTrajectoriesRaw.mat')]);
subject_settings = SettingsCustodian([data_root filesep subject filesep 'subjectSettings.txt']);
emg_period_estimate = subject_settings.get('emg_period_estimate');
emg_offset_x_estimates_by_trial = subject_settings.get('emg_offset_x_estimates_by_trial');
emg_offset_y_estimates_by_channel = subject_settings.get('emg_offset_y_estimates_by_channel');
emg_amplitude_y_estimates_by_channel = subject_settings.get('emg_amplitude_y_estimates_by_channel');
time_emg = loaded_data.time_emg;

filter_order = 4;
cutoff_frequency_low = 3; % in Hz
[b_low, a_low] = butter(filter_order, cutoff_frequency_low/(loaded_data.sampling_rate_emg/2), 'low');


channels_to_process = 1:5;
for i_channel = 1 : length(channels_to_process)
    this_channel = channels_to_process(i_channel);
    data_raw = loaded_data.emg_trajectories_raw(:, channels_to_process(i_channel));
    
    this_offset_time = emg_offset_x_estimates_by_trial(trial_to_process);
    this_offset_gain = emg_offset_y_estimates_by_channel(this_channel);
    this_offset_amplitude = emg_amplitude_y_estimates_by_channel(this_channel);
    
    % method one: direct reweight
    gain_estimate = generateSinusoid(loaded_data.time_emg, emg_period_estimate, this_offset_time, this_offset_gain, this_offset_amplitude);
    data_for_weighting = data_raw;
    data_reweighted_one = data_for_weighting .* gain_estimate.^(-1);
    
    % method two:
    % filter
    filter_order = 4;
    cutoff_frequency_low = 1; % in Hz
    [b_low, a_low] = butter(filter_order, cutoff_frequency_low/(loaded_data.sampling_rate_emg/2), 'low');

    new_indices = (14 : 27 : length(time_emg))';
    new_data_from_min = zeros(size(new_indices));
    new_data_from_max = zeros(size(new_indices));

    for i_index = 1 : length(new_indices)
        index_center = new_indices(i_index);
        index_range = (-13 : 13) + index_center;
        index_range(index_range > length(data_raw)) = [];
        new_data_from_min(i_index) = min(data_raw(index_range));
        new_data_from_max(i_index) = max(data_raw(index_range));
    end


    raw_data_average = (new_data_from_max - new_data_from_min)/2;
    filtered_average = filtfilt(b_low, a_low, raw_data_average);

    % treat filtered oscillation as gain and remove
    cleaned_time = time_emg(new_indices);
    cleaned_data = raw_data_average .* filtered_average.^(-1);

    % upsample
    time_padded = [time_emg(1); cleaned_time; time_emg(end)];
    data_padded = [cleaned_data(1); cleaned_data; cleaned_data(end)];

    if any(~isnan(data_padded))
        data_upsampled = spline(time_padded, data_padded, time_emg);
    else
        data_upsampled = raw_data;
    end
    data_reweighted_two = data_upsampled;
    
    
    figure; hold on
    plot(loaded_data.time_emg, data_raw)
    plot(loaded_data.time_emg, gain_estimate, 'linewidth', 2)
    
    figure; hold on
    plot(loaded_data.time_emg, data_raw * 1e6)
    plot(loaded_data.time_emg, data_reweighted_one)
    plot(loaded_data.time_emg, data_reweighted_two)
end