
collect_data = 0;
fit_sinusoids = 1;
calculate_averages = 1;
print_results = 1;
visualize = 1;

data_root = '/Users/reimajbi/Drive_UD/GVS';
subject = 'AJH';

%% collect data
channels_to_analyze = 1 : 12;
number_of_channels = 16;
number_of_trials = 20;
number_of_time_points = 120 * 2000;
number_of_channels_to_analyze = length(channels_to_analyze);
if collect_data
    load([data_root filesep subject filesep 'subjectInfo.mat']);
    collected_data = cell(number_of_channels, 1);
    for i_trial = 1 : number_of_trials
        loaded_data = load([data_root filesep subject filesep 'raw' filesep makeFileName(date, subject_id, 'walking', i_trial, 'emgTrajectoriesRaw.mat')]);
        for i_channel = 1 : number_of_channels
            collected_data{i_channel} = [collected_data{i_channel} loaded_data.emg_trajectories_raw(:, i_channel)];
        end
    end
    time_emg = loaded_data.time_emg;
end

%% calculate averages of windowed min and max
if fit_sinusoids
    % define filter
    filter_order = 4;
    cutoff_frequency_low = 3; % in Hz
    [b_low, a_low] = butter(filter_order, cutoff_frequency_low/(loaded_data.sampling_rate_emg/2), 'low');
    
    
    new_indices = (14 : 27 : length(time_emg))';
    
    % find time offset
    channels_to_analyze_for_time_offset = [5 10];
    number_of_channels_to_analyze_for_time_offset = length(channels_to_analyze_for_time_offset);
    period_data = zeros(number_of_channels_to_analyze, number_of_trials);
    offset_x_data = zeros(number_of_channels_to_analyze, number_of_trials);
    offset_y_data = zeros(number_of_channels_to_analyze, number_of_trials);
    amplitude_data = zeros(number_of_channels_to_analyze, number_of_trials);
    for i_channel = 1 : length(channels_to_analyze_for_time_offset)
        for i_trial = 1 : number_of_trials
            raw_data = collected_data{channels_to_analyze_for_time_offset(i_channel)}(:, i_trial);
            min_data = zeros(size(new_indices));
            max_data = zeros(size(new_indices));

            for i_index = 1 : length(new_indices)
                index_center = new_indices(i_index);
                index_range = (-13 : 13) + index_center;
                index_range(index_range > length(raw_data)) = [];
                min_data(i_index) = min(raw_data(index_range));
                max_data(i_index) = max(raw_data(index_range));
            end
            raw_data_average = (max_data - min_data)/2;
            
            % filter
            filtered_average = filtfilt(b_low, a_low, raw_data_average);
            data_to_fit = raw_data_average;
%             fit_data = filtered_average;

            % fit
            time = time_emg(new_indices);
            numeric_scale = 1e6;
            temporal_scale = 1/70;
            [~, max_index] = max(filtered_average);
            time_offset_guess = time(max_index) - 70/4;
%             [fit_data, fit_period, fit_offset_x, fit_offset_y, fit_amplitude] = sinusoidFit(time, raw_data_average*numeric_scale, 70, 0, 0, 1);
            [fit_data, fit_period, fit_offset_x, fit_offset_y, fit_amplitude] = ...
                sinusoidFit(time, data_to_fit*numeric_scale, 70, time_offset_guess, 0, range(data_to_fit)*numeric_scale);
            
            period_data(channels_to_analyze_for_time_offset(i_channel), i_trial) = fit_period;
            offset_x_data(channels_to_analyze_for_time_offset(i_channel), i_trial) = fit_offset_x;
            offset_y_data(channels_to_analyze_for_time_offset(i_channel), i_trial) = fit_offset_y * 1/numeric_scale;
            amplitude_data(channels_to_analyze_for_time_offset(i_channel), i_trial) = fit_amplitude * 1/numeric_scale;
            
        end
        
    end
    for i_trial = 1 : number_of_trials
        this_trial_offsets = offset_x_data(channels_to_analyze_for_time_offset, i_trial);
        this_trial_center = this_trial_offsets(1);
        while any(this_trial_offsets > this_trial_center + 35)
            this_trial_offsets(this_trial_offsets > this_trial_center + 35) = this_trial_offsets(this_trial_offsets > this_trial_center + 35) - 70;
        end
        while any(this_trial_offsets < this_trial_center - 35)
            this_trial_offsets(this_trial_offsets < this_trial_center - 35) = this_trial_offsets(this_trial_offsets < this_trial_center - 35) + 70;
        end
        offset_x_data(channels_to_analyze_for_time_offset, i_trial) = this_trial_offsets;
    end
    offset_x_estimates_by_trial = mean(offset_x_data(channels_to_analyze_for_time_offset, :), 1);
    
    % fit all channels
    channels_to_analyze_for_time_offset = [5 10];
    average_data = cell(length(channels_to_analyze), 1);
    period_data = zeros(number_of_channels_to_analyze, number_of_trials);
%     offset_x_data = zeros(number_of_channels_to_analyze, number_of_trials);
    offset_y_data = zeros(number_of_channels_to_analyze, number_of_trials);
    amplitude_data = zeros(number_of_channels_to_analyze, number_of_trials);
    for i_channel = 1 : length(channels_to_analyze)
        for i_trial = 1 : number_of_trials
            raw_data = collected_data{i_channel}(:, i_trial);
            min_data = zeros(size(new_indices));
            max_data = zeros(size(new_indices));

            for i_index = 1 : length(new_indices)
                index_center = new_indices(i_index);
                index_range = (-13 : 13) + index_center;
                index_range(index_range > length(raw_data)) = [];
                min_data(i_index) = min(raw_data(index_range));
                max_data(i_index) = max(raw_data(index_range));
            end
            raw_data_average = (max_data - min_data)/2;
            average_data{i_channel} = [average_data{i_channel} raw_data_average];
            
            % filter
            filtered_average = filtfilt(b_low, a_low, raw_data_average);

            % fit
            time = time_emg(new_indices);
            numeric_scale = 1e6;
            temporal_scale = 1/70;
            [~, max_index] = max(filtered_average);
            time_offset_guess = time(max_index) - 70/4;
%             [fit_data, fit_period, fit_offset_x, fit_offset_y, fit_amplitude] = sinusoidFit(time, raw_data_average*numeric_scale, 70, 0, 0, 1);
            [fit_data, fit_period, fit_offset_x, fit_offset_y, fit_amplitude] = ...
                sinusoidFit_constrainedOffset(time, filtered_average*numeric_scale, 70, offset_x_estimates_by_trial(i_trial), 0, range(filtered_average)*numeric_scale);
            
            period_data(i_channel, i_trial) = fit_period;
%             offset_x_data(i_channel, i_trial) = fit_offset_x;
            offset_y_data(i_channel, i_trial) = fit_offset_y * 1/numeric_scale;
            amplitude_data(i_channel, i_trial) = fit_amplitude * 1/numeric_scale;
            
%             figure; hold on;
%             plot(time, raw_data_average)
%             plot(time, filtered_average, 'linewidth', 2)
%             plot(time, fit_data * 1/numeric_scale, 'linewidth', 2)
        end    
    %     figure; plot(period_data{i_channel, 1}); title('period')
    %     figure; plot(offset_x_data{i_channel, 1}); title('offset_x_data')
    %     figure; plot(offset_y_data{i_channel, 1}); title('offset_y_data')
    %     figure; plot(amplitude_data{i_channel, 1}); title('amplitude_data')
    end    
end

%% calculate averages
if calculate_averages
    period_estimate = median(reshape(period_data, 1, numel(period_data)));
    offset_y_estimates_by_channel = median(offset_y_data, 2);
    amplitude_y_estimates_by_channel = median(amplitude_data, 2);
    
    % clamp time-offset to median \pm 35
    for i_trial = 1 : number_of_trials
        this_trial_offsets = offset_x_data(:, i_trial);
%         this_trial_center = median(this_trial_offsets);
        this_trial_center = 0;
        while any(this_trial_offsets > this_trial_center + 35)
            this_trial_offsets(this_trial_offsets > this_trial_center + 35) = this_trial_offsets(this_trial_offsets > this_trial_center + 35) - 70;
        end
        while any(this_trial_offsets < this_trial_center - 35)
            this_trial_offsets(this_trial_offsets < this_trial_center - 35) = this_trial_offsets(this_trial_offsets < this_trial_center - 35) + 70;
        end
    end
%     offset_x_data(:, i_trial) = this_trial_offsets;
%     offset_x_estimates_by_trial = median(offset_x_data, 1);
    
end

%% print results
if print_results
    fprintf('\n')
    fprintf([subject '\n']);
    
    % period
    fprintf('emg_period_estimate: ');
    fprintf(num2str(period_estimate));
    fprintf('\n')
    
    % offset_x
    fprintf('emg_offset_x_estimates_by_trial: ');
    for i_trial = 1 : number_of_trials-1
        fprintf(num2str(offset_x_estimates_by_trial(i_trial)));
        fprintf(', ')
    end
    fprintf(num2str(offset_x_estimates_by_trial(end)));
    fprintf('\n')

    % offset_y
    fprintf('emg_offset_y_estimates_by_channel: ');
    for i_channel = 1 : number_of_channels_to_analyze-1
        fprintf(num2str(offset_y_estimates_by_channel(i_channel)));
        fprintf(', ')
    end
    fprintf(num2str(offset_y_estimates_by_channel(end)));
    fprintf('\n')
    
    % amplitude_y
    fprintf('emg_amplitude_y_estimates_by_channel: ');
    for i_channel = 1 : number_of_channels_to_analyze-1
        fprintf(num2str(amplitude_y_estimates_by_channel(i_channel)));
        fprintf(', ')
    end
    fprintf(num2str(amplitude_y_estimates_by_channel(end)));
    fprintf('\n')
    
    fprintf('\n')
end

%% visualize
if visualize
    face_color = [1 1 1]*0.2;
    face_alpha = 0.2;

%     % period should be a single constant
%     figure; surf(period_data);  title([subject ' - period'])
%     patch('xdata', [1 number_of_trials number_of_trials 1], 'ydata', [1 1 number_of_channels_to_analyze number_of_channels_to_analyze], 'zdata', [1 1 1 1]*period_estimate, 'FaceColor', face_color, 'FaceAlpha', face_alpha)

    % time offset for specific channels
    
%     offset_x_data_clean(abs(offset_x_data_clean) > 50) = NaN;
    figure; hold on;
    for i_channel = 1 : length(channels_to_analyze_for_time_offset)
        plot(1:number_of_trials, offset_x_data(channels_to_analyze_for_time_offset(i_channel), :))
    end
    
%     surf(offset_x_data_clean); title([subject ' - x-offset']); xlabel('trial'); ylabel('channel');
%     for ...
%         plot3(i_trial*[1 1], [1 number_of_channels_to_analyze], offset_x_estimates_by_trial(i_trial)*[1 1], 'linewidth', 3, 'color', 'r');
%     end

    % amplitude offset should vary by channel, but not by trial
    figure; hold on
    surf(offset_y_data); title([subject ' - y-offset']); xlabel('trial'); ylabel('channel')
    for i_channel = 1 : number_of_channels_to_analyze
        plot3([1 number_of_trials], [1 1]*i_channel, offset_y_estimates_by_channel(i_channel)*[1 1], 'linewidth', 3, 'color', 'r');
    end

    figure; hold on
    surf(amplitude_data); title([subject ' - amplitude']); xlabel('trial'); ylabel('channel')
    for i_channel = 1 : number_of_channels_to_analyze
        plot3([1 number_of_trials], [1 1]*i_channel, amplitude_y_estimates_by_channel(i_channel)*[1 1], 'linewidth', 3, 'color', 'r');
    end
end










