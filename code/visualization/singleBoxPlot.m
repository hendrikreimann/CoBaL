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
    addParameter(parser, 'FaceColor', [0.8, 0.8, 0.8])
    addParameter(parser, 'EdgeColor', [0.5, 0.5, 0.5])
    addParameter(parser, 'MarkerColor', [0.5, 0.5, 0.5])
    addParameter(parser, 'EdgeLinewidth', 1)
    addParameter(parser, 'WiskColor', [0.7, 0.2, 0.07])
    addParameter(parser, 'WiskLinewidth', 1)
    addParameter(parser, 'FaceAlpha', 1)
    
    addParameter(parser, 'PlotMean', true)
    addParameter(parser, 'MeanLinewidth', 1)
    addParameter(parser, 'MeanColor', [0.5, 0.5, 0.5])
    
    addParameter(parser, 'PlotMedian', true)
    addParameter(parser, 'MedianColor', [0.7, 0.2, 0.07])
    addParameter(parser, 'MedianLinewidth', 1)
    
    addParameter(parser, 'xlabel', [])
    addParameter(parser, 'BandWidth', '')
    
    addParameter(parser, 'ShowOutliers', true)
    addParameter(parser, 'ShowData', false)
    parse(parser, varargin{:})
    abscissa = parser.Results.abscissa;
    axes_handle = parser.Results.axes;
    width = parser.Results.width;
    facecolor = parser.Results.FaceColor;
    edgecolor = parser.Results.EdgeColor;
    meancolor = parser.Results.MeanColor;
    markercolor = parser.Results.MarkerColor;
    mean_linewidth = parser.Results.MeanLinewidth;
    edge_linewidth = parser.Results.EdgeLinewidth;
    plot_mean = parser.Results.PlotMean;
    median_color = parser.Results.MedianColor;
    median_linewidth = parser.Results.MedianLinewidth;
    wiskcolor = parser.Results.WiskColor;
    wisk_linewidth = parser.Results.WiskLinewidth;
    plot_median = parser.Results.PlotMedian;
    facealpha = parser.Results.FaceAlpha;
    xlabel = parser.Results.xlabel;
    bandwidth = parser.Results.BandWidth;
    show_outliers = parser.Results.ShowOutliers;    
    show_data = parser.Results.ShowData;

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
        'LineWidth', edge_linewidth, ...
        'HandleVisibility', 'off' ...
      );
    if plot_median
        plot(axes_handle, abscissa + width*[-0.5 0.5], [data_median data_median], 'color', median_color, 'linewidth', median_linewidth, 'HandleVisibility', 'off'); % median
    end
    if plot_mean
        plot(axes_handle, abscissa + width*[-0.5 0.5], [data_mean data_mean], 'color', meancolor, 'linewidth', mean_linewidth, 'HandleVisibility', 'off'); % mean
    end
    plot(axes_handle, [abscissa abscissa], [data_quartile_3 data_upper_adjacent], '--', 'color', wiskcolor, 'linewidth', wisk_linewidth, 'HandleVisibility', 'off'); % upper range
    plot(axes_handle, [abscissa abscissa], [data_lower_adjacent data_quartile_1], '--', 'color', wiskcolor, 'linewidth', wisk_linewidth, 'HandleVisibility', 'off'); % lower range
    plot(axes_handle, abscissa+width*[-0.25 0.25], [data_lower_adjacent data_lower_adjacent], '-', 'color', wiskcolor, 'linewidth', wisk_linewidth, 'HandleVisibility', 'off'); % max
    plot(axes_handle, abscissa+width*[-0.25 0.25], [data_upper_adjacent data_upper_adjacent], '-', 'color', wiskcolor, 'linewidth', wisk_linewidth, 'HandleVisibility', 'off'); % min
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
        xticklabels{xtick == abscissa} = xlabel;
        set(axes_handle, 'xticklabel', xticklabels);
    end
end