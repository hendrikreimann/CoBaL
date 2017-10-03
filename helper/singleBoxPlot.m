%     This file is part of the CoBaL code base
%     Copyright (C) 2017 Hendrik Reimann <hendrikreimann@gmail.com>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

function singleBoxPlot(target_axes_handle, abscissa, data, color, label, show_outliers, width)
    % these should be name-value pair arguments later
    if nargin < 6
        show_outliers = true;
    end
    if nargin < 7
        width = 0.8;
    end
    
    

    if length(data) == 1
        bar ...
          ( ...
            abscissa, ...
            data, ...
            width, ...
            'FaceColor', color, ...
            'parent', target_axes_handle, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off' ...
          );
        return
    end
    if iscolumn(data)
        data = data';
    end
    
    % extract data
    data_median = median(data);
    data_quartile_1 = prctile(data, 25);
    data_quartile_3 = prctile(data, 75);
    data_iqr = iqr(data);
    data_upper_inner_fence = data_quartile_3 + 1.5*data_iqr;
    data_lower_inner_fence = data_quartile_1 - 1.5*data_iqr;
    data_upper_adjacent = max([data(data<=data_upper_inner_fence) -inf]);
    data_lower_adjacent = min([data(data>=data_lower_inner_fence) inf]);
    outliers = data(data>data_upper_inner_fence | data<data_lower_inner_fence);
    
    % plot
    box_x_data = [abscissa-width/2 abscissa+width/2 abscissa+width/2 abscissa-width/2 abscissa-width/2];
    box_y_data = [data_quartile_1 data_quartile_1 data_quartile_3 data_quartile_3 data_quartile_1];
    patch ...
      ( ...
        box_x_data, ...
        box_y_data, ...
        color, ...
        'parent', target_axes_handle, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off' ...
      );
    plot(target_axes_handle, abscissa + width*[-0.5 0.5], [data_median data_median], 'color', 'k', 'HandleVisibility', 'off'); % median
    plot(target_axes_handle, [abscissa abscissa], [data_quartile_3 data_upper_adjacent], 'k--', 'HandleVisibility', 'off'); % upper range
    plot(target_axes_handle, [abscissa abscissa], [data_lower_adjacent data_quartile_1], 'k--', 'HandleVisibility', 'off'); % lower range
    plot(target_axes_handle, abscissa+width*[-0.25 0.25], [data_lower_adjacent data_lower_adjacent], 'k-', 'HandleVisibility', 'off'); % max
    plot(target_axes_handle, abscissa+width*[-0.25 0.25], [data_upper_adjacent data_upper_adjacent], 'k-', 'HandleVisibility', 'off'); % min
    if show_outliers
        plot(target_axes_handle, abscissa * ones(size(outliers)), outliers, '+', 'color', [1; 1; 1] * 0.7, 'HandleVisibility', 'off');
    end
    
    % labels
    if ~isempty(label)
        xtick = get(target_axes_handle, 'xtick');
        if ~ismember(abscissa, xtick)
            xtick = sort([xtick, abscissa]);
            set(target_axes_handle, 'xtick', xtick);
        end
        xticklabels = get(target_axes_handle, 'xticklabel');
        xticklabels{xtick == abscissa} = label;
        set(target_axes_handle, 'xticklabel', xticklabels);
    end
end