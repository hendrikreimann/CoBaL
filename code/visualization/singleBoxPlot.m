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

function singleBoxPlot(data, varargin)
    %% parse input
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'abscissa', 1)
    addParameter(parser, 'axes', gca)
    addParameter(parser, 'width', 0.8)
    addParameter(parser, 'facecolor', [0.8, 0.8, 0.8])
    addParameter(parser, 'edgecolor', [0.5, 0.5, 0.5])
    addParameter(parser, 'meancolor', [0.5, 0.5, 0.5])
    addParameter(parser, 'markercolor', [0.5, 0.5, 0.5])
    addParameter(parser, 'mean_linewidth', 1)
    addParameter(parser, 'plot_mean', true)
    addParameter(parser, 'mediancolor', [0.7, 0.2, 0.07])
    addParameter(parser, 'median_linewidth', 1)
    addParameter(parser, 'plot_median', true)
    addParameter(parser, 'facealpha', 1)
    addParameter(parser, 'xlabel', [])
    addParameter(parser, 'bandwidth', '')
    addParameter(parser, 'show_outliers', true)
    addParameter(parser, 'show_data', true)
    parse(parser, varargin{:})
    abscissa = parser.Results.abscissa;
    axes_handle = parser.Results.axes;
    width = parser.Results.width;
    facecolor = parser.Results.facecolor;
    edgecolor = parser.Results.edgecolor;
    meancolor = parser.Results.meancolor;
    markercolor = parser.Results.markercolor;
    mean_linewidth = parser.Results.mean_linewidth;
    plot_mean = parser.Results.plot_mean;
    mediancolor = parser.Results.mediancolor;
    median_linewidth = parser.Results.median_linewidth;
    plot_median = parser.Results.plot_median;
    facealpha = parser.Results.facealpha;
    xlabel = parser.Results.xlabel;
    bandwidth = parser.Results.bandwidth;
    show_outliers = parser.Results.show_outliers;    
    show_data = parser.Results.show_data;

    if length(data) == 1
        bar ...
          ( ...
            abscissa, ...
            data, ...
            width, ...
            'FaceColor', facecolor, ...
            'parent', axes_handle, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off' ...
          );
        return
    end
    if iscolumn(data)
        data = data';
    end
    
    % extract data
    data_median = nanmedian(data);
    data_mean = nanmean(data);
    data_quartile_1 = prctile(data, 25);
    data_quartile_3 = prctile(data, 75);
    data_iqr = iqr(data);
    data_upper_inner_fence = data_quartile_3 + 1.5*data_iqr;
    data_lower_inner_fence = data_quartile_1 - 1.5*data_iqr;
    data_upper_adjacent = max([data(data<=data_upper_inner_fence) -inf]);
    data_lower_adjacent = min([data(data>=data_lower_inner_fence) inf]);
    outliers = data(data>data_upper_inner_fence | data<data_lower_inner_fence);
    
    % plot
    axes_hold = ishold(axes_handle);
    hold(axes_handle, 'on');
    box_x_data = [abscissa-width/2 abscissa+width/2 abscissa+width/2 abscissa-width/2 abscissa-width/2];
    box_y_data = [data_quartile_1 data_quartile_1 data_quartile_3 data_quartile_3 data_quartile_1];
    patch ...
      ( ...
        'xdata', box_x_data, ...
        'ydata', box_y_data, ...
        'facecolor', facecolor, ...
        'parent', axes_handle, ...
        'EdgeColor', edgecolor, ...
        'HandleVisibility', 'off' ...
      );
    plot(axes_handle, abscissa + width*[-0.5 0.5], [data_median data_median], 'color', 'k', 'HandleVisibility', 'off'); % median
    plot(axes_handle, abscissa + width*[-0.5 0.5], [data_mean data_mean], 'color', [1 1 1]*0.5, 'HandleVisibility', 'off'); % mean
    plot(axes_handle, [abscissa abscissa], [data_quartile_3 data_upper_adjacent], 'k--', 'HandleVisibility', 'off'); % upper range
    plot(axes_handle, [abscissa abscissa], [data_lower_adjacent data_quartile_1], 'k--', 'HandleVisibility', 'off'); % lower range
    plot(axes_handle, abscissa+width*[-0.25 0.25], [data_lower_adjacent data_lower_adjacent], 'k-', 'HandleVisibility', 'off'); % max
    plot(axes_handle, abscissa+width*[-0.25 0.25], [data_upper_adjacent data_upper_adjacent], 'k-', 'HandleVisibility', 'off'); % min
    if show_outliers
        plot(axes_handle, abscissa * ones(size(outliers)), outliers, '+', 'color', [1; 1; 1] * 0.7, 'HandleVisibility', 'off');
    end
    
    if show_data
        scatter(axes_handle, ones(size(data))*abscissa, data, 'jitter', 1, 'jitterAmount', width/2, 'MarkerFaceColor', markercolor, 'MarkerEdgeColor', 'none')
    end
    
    
    if ~axes_hold
        hold(axes_handle, 'off');
    end
    
    % labels
    if ~isempty(xlabel)
        xtick = get(axes_handle, 'xtick');
        if ~ismember(abscissa, xtick)
            xtick = sort([xtick, abscissa]);
            set(axes_handle, 'xtick', xtick);
        end
        xticklabels = get(axes_handle, 'xticklabel');
        xticklabels{xtick == abscissa} = label;
        set(axes_handle, 'xticklabel', xticklabels);
    end
end