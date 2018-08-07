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
    addParameter(parser, 'markersize', 8)
    addParameter(parser, 'linewidth', 1)
    addParameter(parser, 'plot_mean', true)
    addParameter(parser, 'xlabel', [])
    addParameter(parser, 'bandwidth', '')
    parse(parser, varargin{:})
    abscissa = parser.Results.abscissa;
    axes_handle = parser.Results.axes;
    width = parser.Results.width;
    color = parser.Results.color;
    marker = parser.Results.marker;
    markersize = parser.Results.markersize;
    linewidth = parser.Results.linewidth;
    plot_mean = parser.Results.plot_mean;
    xlabel = parser.Results.xlabel;
    bandwidth = parser.Results.bandwidth;
    
    % scatter x-values
    x_values = abscissa + (rand(numel(data), 1) - 0.5) * width;
    
    
    % plot
    plot ...
      ( ...
        x_values, ...
        data, ...
        'Parent', axes_handle, ...
        'marker', marker, ...
        'markersize', markersize, ...
        'linewidth', linewidth*2, ...
        'linestyle', 'none', ...
        'color', color ...
      );
    
    % plot mean
    if plot_mean
        plot ...
          ( ...
            abscissa, ...
            mean(data), ...
            'Parent', axes_handle, ...
            'marker', marker, ...
            'markersize', markersize*3, ...
            'linewidth', linewidth*5, ...
            'linestyle', 'none', ...
            'color', color ...
          );
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