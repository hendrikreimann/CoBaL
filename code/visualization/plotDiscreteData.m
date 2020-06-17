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

function plotDiscreteData(data, varargin)
    settings = parseSettings(varargin{:});
    
    % extract statistics
    stats = extractStatistics(settings, data);
    
    % remember hold status (treat the whole thing created here as one plot)
    axes_hold = ishold(settings.axes);
    hold(settings.axes, 'on');
    
    % plot spread
    if settings.spread.show
        if strcmp(settings.spread.style, 'box')
            plotBox(settings, stats)
        elseif strcmp(settings.spread.style, 'violin')
            plotViolin(settings, stats)
        else
            warning(['Plotting spread with style "' settings.spread.style '" not supported, please change in settings'])
        end
        
    end
    
    % plot mean
    if settings.mean.show
        plotMark(settings, stats, 'mean')
    end
    
    % plot median
    if settings.median.show
        plotMark(settings, stats, 'median')
    end
    
    % plot individual data
    if settings.individual_data.show
        plotIndividualData(settings, data, stats)
    end
    
    % set label
    if ~isempty(settings.label)
        xtick = get(settings.axes, 'xtick');
        if ~ismember(settings.abscissa, xtick)
            xtick = sort([xtick, settings.abscissa]);
            set(settings.axes, 'xtick', xtick);
        end
        xticklabels = get(settings.axes, 'xticklabel');
        xticklabels{xtick == settings.abscissa} = settings.label;
        set(settings.axes, 'xticklabel', xticklabels);
    end
    
    % return hold to previous state
    if ~axes_hold
        hold(settings.axes, 'off');
    end
    
    
    
    
end

function settings = parseSettings(varargin)
    % set up parser
    parser = inputParser;
    parser.KeepUnmatched = true;
    
    % where to plot
    addParameter(parser, 'abscissa', 1)
    addParameter(parser, 'axes', gca)
    addParameter(parser, 'width', 0.8)
    addParameter(parser, 'label', [])
    
    % general color settings
    addParameter(parser, 'color', [])
    
    % individual data
    addParameter(parser, 'ShowIndividualData', 1)
    addParameter(parser, 'IndividualDataMarkerStyle', 'o')
    addParameter(parser, 'IndividualDataMarkerSize', 18)
    addParameter(parser, 'IndividualDataColor', [0.5, 0.5, 0.5])
    addParameter(parser, 'seed', 0)
    
    % mean
    addParameter(parser, 'ShowMean', true)
    addParameter(parser, 'MeanStyle', 'd')
    addParameter(parser, 'MeanMarkerSize', 18)
    addParameter(parser, 'MeanColor', [0.5, 0.5, 0.5])
    addParameter(parser, 'MeanLinewidth', 3)
    
    % median
    addParameter(parser, 'ShowMedian', true)
    addParameter(parser, 'MedianStyle', 'd')
    addParameter(parser, 'MedianMarkerSize', 18)
    addParameter(parser, 'MedianColor', [0.5, 0.5, 0.5])
    addParameter(parser, 'MedianLinewidth', 3)
    
    % spread
    addParameter(parser, 'ShowSpread', 1)
    addParameter(parser, 'SpreadStyle', 'box')
    addParameter(parser, 'SpreadFaceColor', [0.8, 0.8, 0.8])
    addParameter(parser, 'SpreadFaceAlpha', 1)
    addParameter(parser, 'SpreadEdgeColor', [0.5, 0.5, 0.5])
    addParameter(parser, 'SpreadEdgeLinewidth', 1)
    addParameter(parser, 'SpreadWiskColor', [0.5, 0.5, 0.5])
    addParameter(parser, 'SpreadWiskLinewidth', 1)
    addParameter(parser, 'SpreadWiskLinestyle', '-')
    addParameter(parser, 'SpreadDensityEstimateBandwidth', [])
        
    % do the actual parsing
    parse(parser, varargin{:})
    
    % store results in settings struct
    settings = struct;
    settings.abscissa = parser.Results.abscissa;
    settings.axes = parser.Results.axes;
    settings.width = parser.Results.width;
    settings.label = parser.Results.label;
    settings.color = parser.Results.color;
    
    settings.individual_data.show = parser.Results.ShowIndividualData;
    settings.individual_data.marker_style = parser.Results.IndividualDataMarkerStyle;
    settings.individual_data.marker_size = parser.Results.IndividualDataMarkerSize;
    settings.individual_data.color = parser.Results.IndividualDataColor;
    settings.seed = parser.Results.seed;
    
    settings.mean.show = parser.Results.ShowMean;
    settings.mean.style = parser.Results.MeanStyle;
    settings.mean.marker_size = parser.Results.MeanMarkerSize;
    settings.mean.color = parser.Results.MeanColor;
    settings.mean.linewidth = parser.Results.MeanLinewidth;
    
    settings.median.show = parser.Results.ShowMedian;
    settings.median.style = parser.Results.MedianStyle;
    settings.median.marker_size = parser.Results.MedianMarkerSize;
    settings.median.color = parser.Results.MedianColor;
    settings.median.linewidth = parser.Results.MedianLinewidth;
    
    settings.spread.show = parser.Results.ShowSpread;
    settings.spread.style = parser.Results.SpreadStyle;
    settings.spread.face_color = parser.Results.SpreadFaceColor;
    settings.spread.face_alpha = parser.Results.SpreadFaceAlpha;
    settings.spread.edge_color = parser.Results.SpreadEdgeColor;
    settings.spread.edge_linewidth = parser.Results.SpreadEdgeLinewidth;
    settings.spread.wisk_color = parser.Results.SpreadWiskColor;
    settings.spread.wisk_linewidth = parser.Results.SpreadWiskLinewidth;
    settings.spread.wisk_linestyle = parser.Results.SpreadWiskLinestyle;
    settings.spread.density_estimate_bandwidth = parser.Results.SpreadDensityEstimateBandwidth;
    
    % apply general color setting if there was one
    if ~isempty(settings.color) && any(strcmp(parser.UsingDefaults, 'MeanColor')) 
        settings.mean.color = lightenColor(settings.color, 0.5);
    end
    if ~isempty(settings.color) && any(strcmp(parser.UsingDefaults, 'MedianColor')) 
        settings.median.color = lightenColor(settings.color, 0.5);
    end
    if ~isempty(settings.color) && any(strcmp(parser.UsingDefaults, 'IndividualDataColor')) 
        settings.individual_data.color = lightenColor(settings.color, 0.5);
    end
    if ~isempty(settings.color) && any(strcmp(parser.UsingDefaults, 'SpreadFaceColor')) 
        settings.spread.face_color = settings.color;
    end
    
end

function stats = extractStatistics(settings, data)
    stats = struct;
    stats.median = nanmedian(data);
    stats.mean = nanmean(data);
    stats.quartile_1 = prctile(data, 25);
    stats.quartile_3 = prctile(data, 75);
    stats.inter_quartile_range = iqr(data);
    stats.upper_inner_fence = stats.quartile_3 + 1.5*stats.inter_quartile_range;
    stats.lower_inner_fence = stats.quartile_1 - 1.5*stats.inter_quartile_range;
    stats.upper_adjacent = max([data(data<=stats.upper_inner_fence) -inf]);
    stats.lower_adjacent = min([data(data>=stats.lower_inner_fence) inf]);
    
	% Calculate the kernel density
    if ~isempty(settings.spread.density_estimate_bandwidth)
        [stats.data_density, stats.data_range] = ksdensity(data, 'bandwidth', bandwidth);
    else
        [stats.data_density, stats.data_range] = ksdensity(data);
    end
    stats.density_normalized = stats.data_density/max(stats.data_density) * settings.width/2; %normalize
end

function plotBox(settings, stats)
    % box
    box_x_data = settings.abscissa + [-settings.width/2 settings.width/2 settings.width/2 -settings.width/2 -settings.width/2];
    box_y_data = [stats.quartile_1 stats.quartile_1 stats.quartile_3 stats.quartile_3 stats.quartile_1];
    patch ...
      ( ...
        'xdata', box_x_data, ...
        'ydata', box_y_data, ...
        'parent', settings.axes, ...
        'FaceColor', settings.spread.face_color, ...
        'EdgeColor', settings.spread.edge_color, ...
        'LineWidth', settings.spread.edge_linewidth, ...
        'HandleVisibility', 'off' ...
      );
  
    % wisks
    plot ...
      ( ...
        settings.axes, ...
        settings.abscissa*[1 1], ...
        [stats.quartile_3 stats.upper_adjacent], ...
        'linestyle', settings.spread.wisk_linestyle, ...
        'color', settings.spread.wisk_color, ...
        'linewidth', settings.spread.wisk_linewidth, ...
        'HandleVisibility', 'off' ...
      ); % upper range
    plot ...
      ( ...
        settings.axes, ...
    	settings.abscissa*[1 1], ...
        [stats.lower_adjacent stats.quartile_1], ...
        'linestyle', settings.spread.wisk_linestyle, ...
        'color', settings.spread.wisk_color, ...
        'linewidth', settings.spread.wisk_linewidth, ...
        'HandleVisibility', 'off' ...
      ); % lower range
    plot ...
      ( ...
        settings.axes, ...
    	settings.abscissa+settings.width*[-0.25 0.25], ...
        [stats.lower_adjacent stats.lower_adjacent], ...
        'linestyle', settings.spread.wisk_linestyle, ...
        'color', settings.spread.wisk_color, ...
        'linewidth', settings.spread.wisk_linewidth, ...
        'HandleVisibility', 'off' ...
        ); % max
    plot ...
      ( ...
        settings.axes, ...
    	settings.abscissa+settings.width*[-0.25 0.25], ...
        [stats.upper_adjacent stats.upper_adjacent], ...
        'linestyle', settings.spread.wisk_linestyle, ...
        'color', settings.spread.wisk_color, ...
        'linewidth', settings.spread.wisk_linewidth, ...
        'HandleVisibility', 'off' ...
      ); % min
end

function plotViolin(settings, stats)

    % plot density estimate
    x_values = settings.abscissa + [-flip(stats.density_normalized) stats.density_normalized];
    y_values = [flip(stats.data_range) stats.data_range];
    patch ...
      ( ...
        'xdata', x_values, ...
        'ydata', y_values, ...
        'parent', settings.axes, ...
        'FaceColor', settings.spread.face_color, ...
        'EdgeColor', settings.spread.edge_color, ...
        'LineWidth', settings.spread.edge_linewidth, ...
        'HandleVisibility', 'off' ...
      );
  
end

function plotMark(settings, stats, location)
    ydata = stats.(location);
    style = settings.(location).style;
    
    if strcmp(style, 'line')
        if settings.spread.show && strcmp(settings.spread.style, 'violin')
            % use width of the violin at the desired location
            width = interp1(stats.data_range, stats.density_normalized, ydata);
        else
            width = settings.width * 0.5;
        end
        
        % plot line
        plot ...
          ( ...
            settings.abscissa + [-width width], ...
            [ydata ydata], ...
            'Parent', settings.axes, ...
            'linewidth', settings.(location).linewidth, ...
            'color', settings.(location).color ...
          );
    elseif strcmp(style, 'bar')
        box_x_data = settings.abscissa + [-settings.width/2 settings.width/2 settings.width/2 -settings.width/2 -settings.width/2];
        box_y_data = [0 0 ydata ydata 0];
        patch ...
          ( ...
            box_x_data, ...
            box_y_data, ...
            settings.(location).color, ...
            'parent', settings.axes, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off' ...
          );
        
    else
        % plot marker
        plot ...
          ( ...
            settings.abscissa, ...
            ydata, ...
            'Parent', settings.axes, ...
            'Marker', style, ...
            'MarkerSize', settings.(location).marker_size, ...
            'MarkerFaceColor', settings.(location).color, ...
            'MarkerEdgeColor', 'none' ...
          );
    end    
end

function plotIndividualData(settings, data, stats)
    % create plot
    scatterplot = ...
        scatter ...
          ( ...
            settings.axes, ...
            ones(size(data))*settings.abscissa, ...
            data, ...
            settings.individual_data.marker_size, ...
            'MarkerFaceColor', settings.individual_data.color, ...
            'MarkerEdgeColor', 'none' ...
          );
      
    % create jitter
    old_stream = RandStream.getGlobalStream;
    new_stream = RandStream.create('mrg32k3a', 'seed', settings.seed);
    RandStream.setGlobalStream(new_stream);
    jitter = 2*rand(size(data)) - 1;
    RandStream.setGlobalStream(old_stream);
        
    % normalize jitter by spread
    for i_point = 1 : length(data)
        % calculate position in range
        this_point = data(i_point);
        this_point_position_in_range = (this_point - stats.data_range(1)) / (stats.data_range(end) - stats.data_range(1));
        this_point_spread = stats.density_normalized(round(this_point_position_in_range*100));
        jitter(i_point) = jitter(i_point) * this_point_spread;
    end
    
    % apply jitter
    set(scatterplot, 'xdata', ones(size(data))*settings.abscissa + jitter)



end










