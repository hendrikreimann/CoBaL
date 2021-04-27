% apply the specified filter for each stretch of not-NaN data that is long enough
function y = nanfiltfilt(b, a, x)
% function y = nanfiltfilt(x, filterdesign)
    
    y = zeros(size(x));
    number_of_time_steps = size(x, 1);
    for i_column = 1 : size(x, 2)
        column_raw = x(:, i_column);

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
        
        gaps = [gap_start_indices gap_end_indices];

        % find the pieces of data without gaps
        data_stretches_without_gaps = [];
        if isempty(gaps)
            data_stretches_without_gaps = [first_data_point number_of_time_steps];
        else
            if  gaps(1, 1) > 1
                data_stretches_without_gaps = [first_data_point gaps(1, 1)-1];
            end
            for i_gap = 1 : size(gaps, 1) - 1
                % add the stretch after this big gaps to the list
                data_stretches_without_gaps = [data_stretches_without_gaps; gaps(i_gap, 2)+1 gaps(i_gap+1, 1)-1];
            end
            if ~isempty(gaps) && gaps(end, 2) < number_of_time_steps
                data_stretches_without_gaps = [data_stretches_without_gaps; gaps(end, 2)+1 number_of_time_steps];
            end
            % TODO: check if the limit cases are treated correctly

        end
        
        % filter the stretches without gaps
        filter_order = length(a) - 1;
        column_filtered_by_stretch = zeros(size(column_raw)) * NaN;
        for i_stretch = 1 : size(data_stretches_without_gaps, 1)
            data_stretch = column_raw(data_stretches_without_gaps(i_stretch, 1) : data_stretches_without_gaps(i_stretch, 2));
            % standard single filter
%             if length(data_stretch) > 3*filter_order
%                 data_stretch_filtered = filtfilt(b, a, data_stretch);
%                 column_filtered_by_stretch(data_stretches_without_gaps(i_stretch, 1) : data_stretches_without_gaps(i_stretch, 2)) = data_stretch_filtered;
%             end
            % call filtfilt once forward and once backward on this stretch
            if length(data_stretch) > 3*filter_order
                data_stretch_filtered_forward = filtfilt(b, a, data_stretch);
                data_stretch_filtered_backward = flip(filtfilt(b, a, flip(data_stretch)));
                midpoint = round(length(data_stretch) / 2);
                data_stretch_filtered = zeros(size(data_stretch));
                data_stretch_filtered(1 : midpoint) = data_stretch_filtered_forward(1 : midpoint);
                data_stretch_filtered(midpoint+1 : end) = data_stretch_filtered_backward(midpoint+1 : end);
                
%                 hold on; 
%                 plot(data_stretch, 'linewidth', 6); 
%                 plot(data_stretch_filtered_forward, 'linewidth', 4); 
%                 plot(data_stretch_filtered_backward, 'linewidth', 4); 
%                 plot(data_stretch_filtered, 'linewidth', 2)
                
                column_filtered_by_stretch(data_stretches_without_gaps(i_stretch, 1) : data_stretches_without_gaps(i_stretch, 2)) = data_stretch_filtered;
            end

% plot(data_stretch); hold on; plot(data_stretch_filtered)
        end
        
        y(:, i_column) = column_filtered_by_stretch;
    end









end