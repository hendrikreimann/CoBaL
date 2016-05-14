do_plots = 0;

current_directory = pwd;
data_dir = dir('*_markerTrajectoriesRaw.mat');
clear file_name_list;
[file_name_list{1:length(data_dir)}] = deal(data_dir.name);

load('subjectInfo.mat');

filter_order = 2;
cutoff_frequency = 20; % cutoff frequency, in Hz

gap_length_limit = 0.15; % maximum of allowable gap being splined over, in seconds

number_of_files = length(file_name_list);
for i_trial = 1 : number_of_files
    % load data
    raw_marker_file_name = file_name_list{i_trial};
    [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(raw_marker_file_name);
    load(raw_marker_file_name);
    
    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate_mocap/2));	% set filter parameters for butterworth filter: 2=order of filter;
    
    % process data
    marker_trajectories = zeros(size(marker_trajectories_raw));

    if do_plots
        figure; axes; hold on; title(['trial ' trial_id]);
    end
    time = (1/sampling_rate_mocap : 1/sampling_rate_mocap : size(marker_trajectories_raw, 1) / sampling_rate_mocap)';
    number_of_time_steps = length(time);
    time_steps = 1 : number_of_time_steps;
    gap_length_limit_time_steps = ceil(gap_length_limit * sampling_rate_mocap);
    gap_template = ones(gap_length_limit_time_steps, 1);

    for i_column = 1 : size(marker_trajectories_raw, 2)
        column_raw = marker_trajectories_raw(:, i_column);

        % find the gaps
        gaps = isnan(column_raw);
        gap_start_indices = [];
        gap_end_indices = [];
        first_data_point = find(~isnan(column_raw), 1);
        if isempty(first_data_point)
            first_data_point = number_of_time_steps+1;
        end
        for i_time = first_data_point : number_of_time_steps
            if isnan(column_raw(i_time))
                % check if this is the start of a gap
                if (i_time == 1) || (~isnan(column_raw(i_time-1)))
                    gap_start_indices = [gap_start_indices; i_time]; %#ok<*AGROW>
                end
                % check if this is the end of a gap
                if (i_time == number_of_time_steps) || (~isnan(column_raw(i_time+1)))
                    gap_end_indices = [gap_end_indices; i_time];
                end

            end
        end

        % find the big gaps
        big_gaps = [];
        for i_gap = 1 : length(gap_start_indices)
            % check if this is a long gap
            if ((gap_end_indices(i_gap) - gap_start_indices(i_gap)) > gap_length_limit_time_steps)
                big_gaps = [big_gaps; [gap_start_indices(i_gap) gap_end_indices(i_gap)]];
            end
        end

        % find the pieces of data without big gaps
        data_stretches_without_big_gaps = [];
        if isempty(big_gaps)
            data_stretches_without_big_gaps = [first_data_point number_of_time_steps];
        else
            if  big_gaps(1, 1) > 1
                data_stretches_without_big_gaps = [first_data_point big_gaps(1, 1)-1];
            end
            for i_gap = 1 : size(big_gaps, 1) - 1
                % add the stretch after this big gaps to the list
                data_stretches_without_big_gaps = [data_stretches_without_big_gaps; big_gaps(i_gap, 2)+1 big_gaps(i_gap+1, 1)-1];
            end
            if ~isempty(big_gaps) && big_gaps(end, 2) < number_of_time_steps
                data_stretches_without_big_gaps = [data_stretches_without_big_gaps; big_gaps(end, 2)+1 number_of_time_steps];
            end
            % TODO: check if the limit cases are treated correctly

            % make indicator for plotting
            big_gap_indicator = zeros(size(column_raw));
            for i_gap = 1 : size(big_gaps, 1)
                big_gap_indicator((time_steps >= big_gaps(i_gap, 1)) & (time_steps <= big_gaps(i_gap, 2))) = 1;
            end
            big_gap_indicator(big_gap_indicator == 0) = NaN;
        end

        % fill in small holes in the data by spline
        if numel(column_raw) - numel(find(isnan(column_raw))) > 2
            [~, column_splined] = evalc('csaps(time, column_raw, 1, time);');
            % filter the stretches of data without big holes
            column_filtered_by_stretch = zeros(size(column_raw)) * NaN;
            for i_stretch = 1 : size(data_stretches_without_big_gaps, 1)
                data_stretch = column_splined(data_stretches_without_big_gaps(i_stretch, 1) : data_stretches_without_big_gaps(i_stretch, 2));
                if length(data_stretch) > 3*filter_order
                    data_stretch_filtered = filtfilt(b, a, data_stretch);
                    column_filtered_by_stretch(data_stretches_without_big_gaps(i_stretch, 1) : data_stretches_without_big_gaps(i_stretch, 2)) = data_stretch_filtered;
                end
            end
        else
            column_filtered_by_stretch = column_raw;
        end

        marker_trajectories(:, i_column) = column_filtered_by_stretch;

        if do_plots
            plot(time, column_raw, 'b-', 'linewidth', 5);
            plot(time, column_splined, 'r-', 'linewidth', 2);
            plot(time, column_filtered_by_stretch, 'c-', 'linewidth', 2);
            legend('raw data', 'splined', 'splined and filtered')
        end
        
    end

    save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'markerTrajectories');
    save(save_file_name, 'marker_trajectories', 'time_mocap', 'sampling_rate_mocap');
    disp(['filtered and saved as ' save_file_name])
    
    
    
    
    
    
    
    
    
    
end

disp(['filtered ' num2str(number_of_files) ' files'])
