function onset_time = estimateOnsetTime(data_mean, data_cinv, time)

    sampling_rate = 1/median(diff(time));

    exclude_range_start = 0;
    exclude_range_end = 0;
    filter_order = 2;
    cutoff_frequency = 10; % cutoff frequency, in Hz
    [filter_b, filter_a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;

    data_here_filtered = filtfilt(filter_b, filter_a, data_mean);
    data_here_dot = deriveByTime(data_here_filtered, 1/sampling_rate);
    data_here_dot_pruned = [zeros(exclude_range_start, 1); data_here_dot(exclude_range_start+1 : end-exclude_range_end); zeros(exclude_range_end, 1)];


    exclude_range_start = 40;
    cinv_lower = data_mean - data_cinv;
    cinv_upper = data_mean + data_cinv;
    cinv_excludes_zero = (cinv_lower > 0) | cinv_upper < 0;
    % error here
    critical_index = find(cinv_excludes_zero' & (1:length(time)>exclude_range_start), 1, 'first');
    % error here
    if isempty(critical_index) || critical_index >= 98
        onset_time = 0;
        return;
    end
    
    data_of_interest = data_mean(critical_index : end);
    data_dot_of_interest = data_here_dot(critical_index : end);
    if max(data_of_interest) > abs(min(data_of_interest))
        [~, peak_index_local] = findpeaks(data_dot_of_interest);
        if isempty(peak_index_local)
            [~, peak_index_local] = max(data_dot_of_interest);
        end
    else
        [~, peak_index_local] = findpeaks(-data_dot_of_interest);
        if isempty(peak_index_local)
            [~, peak_index_local] = max(-data_dot_of_interest);
        end
    end
    peak_index_localPeak = critical_index - 1 + peak_index_local(1);
    data_here_slope_at_localPeak = data_here_dot(peak_index_localPeak);
    data_here_slope_intersect_localPeak =  - time(peak_index_localPeak)*data_here_slope_at_localPeak + data_here_filtered(peak_index_localPeak);
    onset_time = - data_here_slope_intersect_localPeak / data_here_slope_at_localPeak;

    % report
%     disp(['Variable "' this_variable ' - onset time: ' num2str(onset_time_here_localPeak)]);
    % save(['stats_' this_variable '.mat'], 'h_results_stimRight', 'h_results_stimLeft');

    % visualize onset estimate
    %     time = 1 : 100;
    figure; hold on;
    shadedErrorBar(time, data_mean, data_cinv);
    plot(time, time*0, 'k');
    plot(time, data_here_filtered);
    plot(time, data_here_dot);
    plot(time(peak_index_localPeak), data_mean(peak_index_localPeak), 'o');
    plot(time(peak_index_localPeak), data_here_dot(peak_index_localPeak), 'o');
    plot(time, time*data_here_slope_at_localPeak + data_here_slope_intersect_localPeak)
    plot(onset_time, 0, 'o')

end