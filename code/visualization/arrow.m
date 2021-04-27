function [arrowline, arrowhead] = arrow(varargin)
    % parse input
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'linewidth', 1)
    addParameter(parser, 'color', [0 0 0])
    addParameter(parser, 'start', 0)
    addParameter(parser, 'end', 0)
    addParameter(parser, 'height', 0)
    addParameter(parser, 'width', 0)
    parse(parser, varargin{:})

    startpoint = parser.Results.start;
    endpoint = parser.Results.end;
    linewidth = parser.Results.linewidth;
    color = parser.Results.color;
    height = parser.Results.height;
    width = parser.Results.width;
    
    % calculate default height and width
    if height == 0
        height = 0.05 * norm(endpoint - startpoint);
    end
    if width == 0
        width = height / 2;
    end

    % calculate triangle points
    vector = (endpoint - startpoint) * 1 / norm(endpoint - startpoint);
    normal = [vector(2), -vector(1)];
    corner_one = endpoint - vector * height + normal * width/2;
    corner_two = endpoint - vector * height - normal * width/2;
    line_end = endpoint - vector * height;
    
    % draw line
    arrowline = plot([startpoint(1), line_end(1)], [startpoint(2), line_end(2)], 'linewidth', linewidth, 'color', color);

    xdata = [endpoint(1), corner_one(1), corner_two(1)];
    ydata = [endpoint(2), corner_one(2), corner_two(2)];
    arrowhead = patch('XData', xdata, 'YData', ydata, 'facecolor', color, 'edgecolor', 'none');
    
end