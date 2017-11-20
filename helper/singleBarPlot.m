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

function singleBoxPlot(target_axes_handle, abscissa, data, color, label, width)
    if nargin < 6
        width = 0.8;
    end
    data_mean = mean(data);
    data_std_error = std(data)/ sqrt( length( data ));
      
%     bar ...
%        ( ...
%         abscissa, ...
%         data_mean, ...
%         width, ...
%         'FaceColor', color, ...
%         'parent', target_axes_handle, ...
%         'EdgeColor', 'none', ...
%         'HandleVisibility', 'off' ...
%         );
%     return

%     bar_x_data = [abscissa-width/2 abscissa+width/2 abscissa+width/2 abscissa-width/2 abscissa-width/2];
%     bar_mean_data = [data_quartile_1 data_quartile_1 data_quartile_3 data_quartile_3 data_quartile_1];
%     
%       patch ...
%       ( ...
%         bar_x_data, ...
%         data_mean, ...
%         color, ...
%         'parent', target_axes_handle, ...
%         'EdgeColor', 'none', ...
%         'HandleVisibility', 'off' ...
%       );
    axes(target_axes_handle)
    bar(abscissa, data_mean, 'FaceColor', color) % main bar
    errorbar(abscissa, data_mean, data_std_error, 'k'); % error bars
  
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
