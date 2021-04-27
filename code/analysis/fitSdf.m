%     This file is part of the CoBaL code base
%     Copyright (C) 2019 Hendrik Reimann <hendrikreimann@gmail.com>
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

function [sdf_fit, breakpoint, slope_short_term, slope_long_term] = fitSdf(sdf_data, sdf_time)
    % normalize sdf data
    normalization_factor = sdf_data(end);
    sdf_data_normalized = sdf_data * 1/normalization_factor;


    % initial guess for breakpoint is 2 seconds
    breakpoint_init = 2;
    if sdf_time(end) < 2
        breakpoint_init = sdf_time(end)/2
    end
    [~, breakpoint_index] = min(abs(sdf_time - breakpoint_init));
    sdf_breakpoint = sdf_data_normalized(breakpoint_index);
    sdf_end = sdf_data_normalized(end);
    
    % extract data for short term and long term segments
    time_short_term_init = sdf_time(1 : breakpoint_index);
    data_short_term = sdf_data_normalized(1 : breakpoint_index);
    time_long_term_init = sdf_time(breakpoint_index : end);
    data_long_term = sdf_data_normalized(breakpoint_index+1 : end);
    
    % determine initial guesses for lines
    slope_short_term_init = sdf_breakpoint / sdf_time(breakpoint_index);
    
    coefficients = polyfit([sdf_time(breakpoint_index), sdf_time(end)], [sdf_breakpoint, sdf_end], 1);
    slope_long_term_init = coefficients(1);
    offset_long_term = coefficients(2);    
    
    line_short_term_init = time_short_term_init*slope_short_term_init;
    line_long_term_init = time_long_term_init*slope_long_term_init + offset_long_term;
    
    % set fit options
    options = optimset ...
    ( ...
        'GradObj', 'off', ...
        'Display','off', ...
        'LargeScale', 'off', ...
        'DerivativeCheck', 'on' ...
    );

    % fit
    sdf_parameters_init = [breakpoint_init, slope_short_term_init, slope_long_term_init];
    sdf_parameters_fit = fminunc(@sdfDataDifference, sdf_parameters_init, options);
    
    % unpack parameters
    breakpoint_fit = sdf_parameters_fit(1);
    slope_short_term_fit = sdf_parameters_fit(2);
    slope_long_term_fit = sdf_parameters_fit(3);

    % time vectors
    [~, breakpoint_index_fit] = find(sdf_time < breakpoint_fit, 1, 'last');
    time_short_term_fit = [sdf_time(1 : breakpoint_index_fit) breakpoint_fit];
    time_long_term_fit = [breakpoint_fit sdf_time(breakpoint_index_fit+1 : end)];

    % determine offset
    sdf_breakpoint_fit = slope_short_term_fit * breakpoint_fit;
    offset_long_term_fit = sdf_breakpoint_fit - slope_long_term_fit * breakpoint_fit;

    % construct fits
    line_short_term_fit = time_short_term_fit * slope_short_term_fit;
    line_long_term_fit = time_long_term_fit * slope_long_term_fit + offset_long_term_fit;
    sdf_fit = [line_short_term_fit(1:end-1) line_long_term_fit(2:end)];


    figure; axes; hold on
    plot(sdf_time, sdf_data_normalized, 'linewidth', 3);
    plot(time_short_term_init, line_short_term_init, 'linewidth', 3, 'color', [1 0.2 0]);
    plot(time_long_term_init, line_long_term_init, 'linewidth', 3, 'color', [1 0.2 0]);
    plot(time_short_term_fit, line_short_term_fit, 'linewidth', 3, 'color', [0.2 1 0]);
    plot(time_long_term_fit, line_long_term_fit, 'linewidth', 3, 'color', [0.2 1 0]);
    plot(sdf_time, sdf_fit, 'linewidth', 1);

    % re-normalize
    sdf_fit = sdf_fit * normalization_factor;
    breakpoint = breakpoint_fit;
    slope_short_term = slope_short_term_fit * normalization_factor;
    slope_long_term = slope_long_term_fit * normalization_factor;
    
    
    function f = sdfDataDifference(sdf_parameters)
        % unpack parameters
        breakpoint_here = sdf_parameters(1);
        slope_short_term_here = sdf_parameters(2);
        slope_long_term_here = sdf_parameters(3);
        
        % time vectors
        [~, breakpoint_index_here] = find(sdf_time < breakpoint_here, 1, 'last');
        time_short_term_here = [sdf_time(1 : breakpoint_index_here) breakpoint_here];
        time_long_term_here = [breakpoint_here sdf_time(breakpoint_index_here+1 : end)];
        
        % determine offset
        sdf_breakpoint_here = slope_short_term_here * breakpoint_here;
        offset_long_term_here = sdf_breakpoint_here - slope_long_term_here * breakpoint_here;
        
        % construct fits
        line_short_term_here = time_short_term_here * slope_short_term_here;
        line_long_term_here = time_long_term_here * slope_long_term_here + offset_long_term_here;
        sdf_fit_here = [line_short_term_here(1:end-1) line_long_term_here(2:end)];
        
        % calculate error
        difference = sdf_data_normalized - sdf_fit_here;
        f = sum(difference.^2);
        
%         disp(num2str(f))
        
%         figure; axes; hold on
%         plot(sdf_time, sdf_data_normalized, 'linewidth', 3);
%         plot(time_short_term_here, line_short_term_here, 'linewidth', 3);
%         plot(time_long_term_here, line_long_term_here, 'linewidth', 3);
%         plot(sdf_time, sdf_fit_here, 'linewidth', 1);
%         plot(sdf_time, difference, 'linewidth', 2);
        

    end
    
    
end

