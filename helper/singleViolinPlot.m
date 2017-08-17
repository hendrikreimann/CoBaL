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

% function singleViolinPlot(target_axes_handle, abscissa, data, color, label, show_outliers)
function singleViolinPlot(data, varargin)
    %% parse input
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'abscissa', 1)
    addParameter(parser, 'axes', gca)
    addParameter(parser, 'width', 0.4)
    addParameter(parser, 'facecolor', [0.8, 0.8, 0.8])
    addParameter(parser, 'edgecolor', [0.5, 0.5, 0.5])
    addParameter(parser, 'meancolor', [0.5, 0.5, 0.5])
    addParameter(parser, 'plot_mean', true)
    addParameter(parser, 'mediancolor', [0.7, 0.2, 0.07])
    addParameter(parser, 'plot_median', true)
    addParameter(parser, 'facealpha', 1)
    addParameter(parser, 'xlabel', '')
    addParameter(parser, 'bandwidth', '')
    addParameter(parser, 'show_outliers', true)
    parse(parser, varargin{:})
    abscissa = parser.Results.abscissa;
    axes_handle = parser.Results.axes;
    width = parser.Results.width;
    facecolor = parser.Results.facecolor;
    edgecolor = parser.Results.edgecolor;
    meancolor = parser.Results.meancolor;
    plot_mean = parser.Results.plot_mean;
    mediancolor = parser.Results.mediancolor;
    plot_median = parser.Results.plot_median;
    facealpha = parser.Results.facealpha;
    xlabel = parser.Results.xlabel;
    bandwidth = parser.Results.bandwidth;
    show_outliers = parser.Results.show_outliers;

    data_median = nanmedian(data);
    data_mean = nanmean(data);
    
    if ~show_outliers
        data_quartile_1 = prctile(data, 25);
        data_quartile_3 = prctile(data, 75);
        data_iqr = iqr(data);
        data_upper_inner_fence = data_quartile_3 + 1.5*data_iqr;
        data_lower_inner_fence = data_quartile_1 - 1.5*data_iqr;
        data_upper_adjacent = max([data(data<=data_upper_inner_fence) -inf]);
        data_lower_adjacent = min([data(data>=data_lower_inner_fence) inf]);
        outliers = data(data>data_upper_inner_fence | data<data_lower_inner_fence);
        data = data(data<=data_upper_inner_fence & data>=data_lower_inner_fence);
    end    
    
	% Calculate the kernel density
    if ~isempty(bandwidth)
        [data_density, data_range, bandwidth_used]=ksdensity(data, 'bandwidth', bandwidth);
    end
    if isempty(bandwidth)
        [data_density, data_range, bandwidth_used]=ksdensity(data);
    end
    density_normalized = data_density/max(data_density)*width; %normalize

    % plot density estimate
    x_values = abscissa + [-flip(density_normalized) density_normalized];
    y_values = [flip(data_range) data_range];
    patch ...
      ( ...
        'Parent', axes_handle, ...
        'xdata', x_values, ...
        'ydata', y_values, ...
        'facecolor', facecolor, ...
        'FaceAlpha', facealpha, ...
        'EdgeColor', edgecolor ...
      );
    
    % plot mean
    if plot_mean
        density_at_mean = interp1(data_range, density_normalized, data_mean);
        plot ...
          ( ...
            abscissa + [-density_at_mean density_at_mean], ...
            [data_mean data_mean], ...
            'Parent', axes_handle, ...
            'color', meancolor ...
          );
    end
    
    % plot median
    if plot_median
        density_at_median = interp1(data_range, density_normalized, data_median);
        plot ...
          ( ...
            abscissa + [-density_at_median density_at_median], ...
            [data_median data_median], ...
            'Parent', axes_handle, ...
            'color', mediancolor ...
          );
    end    

    % labels
    xtick = get(axes_handle, 'xtick');
    xticklabels = get(axes_handle, 'xticklabel');
    xticklabels{xtick == abscissa} = xlabel;
    set(axes_handle, 'xticklabel', xticklabels);




end