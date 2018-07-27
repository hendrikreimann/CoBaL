
data_root = '/Users/reimajbi/Drive_UD/GVS';
subject = 'AJH';
trial_to_process = 1;
load([data_root filesep subject filesep 'subjectInfo.mat']);
loaded_data = load([data_root filesep subject filesep 'raw' filesep makeFileName(date, subject_id, 'walking', trial_to_process, 'emgTrajectoriesRaw.mat')]);

raw_data = emg_trajectories_raw(:, 5);

new_indices = (14 : 27 : length(time_emg))';
new_data_from_min = zeros(size(new_indices));
new_data_from_max = zeros(size(new_indices));

for i_index = 1 : length(new_indices)
    index_center = new_indices(i_index);
    index_range = (-13 : 13) + index_center;
    index_range(index_range > length(raw_data)) = [];
    new_data_from_min(i_index) = min(raw_data(index_range));
    new_data_from_max(i_index) = max(raw_data(index_range));
end

raw_data_average = (new_data_from_max - new_data_from_min)/2;

% filter
filter_order = 4;
cutoff_frequency_low = 1; % in Hz
[b_low, a_low] = butter(filter_order, cutoff_frequency_low/(sampling_rate_emg/2), 'low');
filtered_average = filtfilt(b_low, a_low, raw_data_average);

% treat filtered oscillation as gain and remove
cleaned_time = time_emg(new_indices);
cleaned_data = raw_data_average .* filtered_average.^(-1);

time_padded = [time_emg(1); cleaned_time; time_emg(end)];
data_padded = [cleaned_data(1); cleaned_data; cleaned_data(end)];
data_upsampled = spline(time_padded, data_padded, time_emg);


filter_order_final = 4;
cutoff_frequency_final = 3; % in Hz
[b_final, a_final] = butter(filter_order_final, cutoff_frequency_final/(sampling_rate_emg/2), 'low');
data_final = filtfilt(b_final, a_final, data_upsampled);


figure; hold on
plot(raw_data)
plot(new_indices, -new_data_from_min)
plot(new_indices, new_data_from_max)
plot(new_indices, (new_data_from_max - new_data_from_min)/2)
plot(new_indices, filtered_average, 'linewidth', 3)


figure; hold on;
plot(cleaned_time, cleaned_data)
plot(time_emg, data_upsampled)
plot(time_emg, data_final, 'linewidth', 2)


