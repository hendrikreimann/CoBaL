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

    function singleScatterPlot(data, varargin)
    %% parse input
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'abscissa', 1)
    addParameter(parser, 'axes', gca)
    addParameter(parser, 'width', 0.8)
    addParameter(parser, 'color', [0.8, 0.8, 0.8])
    addParameter(parser, 'marker', 'x')
    addParameter(parser, 'DataMarkerSize', 18)
    addParameter(parser, 'linewidth', 1)
    addParameter(parser, 'xlabel', [])
    addParameter(parser, 'seed', 0)
    parse(parser, varargin{:})
    abscissa = parser.Results.abscissa;
    axes_handle = parser.Results.axes;
    width = parser.Results.width;
    color = parser.Results.color;
    marker = parser.Results.marker;
    data_marker_size = parser.Results.DataMarkerSize;
    linewidth = parser.Results.linewidth;
    xlabel = parser.Results.xlabel;
    seed = parser.Results.seed;
        
    [data_density, data_range]=ksdensity(data);
    density_normalized = data_density/max(data_density); %normalize

    old_stream = RandStream.getGlobalStream;
    new_stream = RandStream.create('mrg32k3a', 'seed', seed);
    RandStream.setGlobalStream(new_stream);
    jitter = (rand(size(data)) * width) - width/2;

    % normalize jitter by spread
    for i_point = 1 : length(data)
        % calculate position in range
        this_point = data(i_point);
        this_point_position_in_range = (this_point - data_range(1)) / (data_range(end) - data_range(1));
        this_point_spread = density_normalized(round(this_point_position_in_range*100));
        jitter(i_point) = jitter(i_point) * this_point_spread;
    end

    scatterplot = scatter(axes_handle, ones(size(data))*abscissa, data, data_marker_size, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none');
    set(scatterplot, 'xdata', ones(size(data))*abscissa + jitter)
%         uistack(scatterplot, 'bottom');
    RandStream.setGlobalStream(old_stream);
        
%     % plot
%     plot ...
%       ( ...
%         x_values, ...
%         data, ...
%         'Parent', axes_handle, ...
%         'marker', marker, ...
%         'markersize', markersize, ...
%         'linewidth', linewidth*2, ...
%         'linestyle', 'none', ...
%         'color', color ...
%       );        

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