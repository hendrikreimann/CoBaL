%     This file is part of the CoBaL code base
%     Copyright (C) 2017-18 Hendrik Reimann <hendrikreimann@gmail.com>
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

function [scaled_abscissa, band_limits] = createScaledAbscissa(band_scales, number_of_time_steps_normalized)
    number_of_bands = length(band_scales);
    scaled_abscissa = [];
    band_limits = 0;
    for i_band = 1 : number_of_bands
        scaled_abscissa_this_band = linspace(0, band_scales(i_band), number_of_time_steps_normalized)';
        if i_band > 1
            % start time of this band is end time of the last band, so remove the duplicate point
            scaled_abscissa_this_band = scaled_abscissa_this_band + scaled_abscissa(end);
            scaled_abscissa_this_band = scaled_abscissa_this_band(2:end);
        end
        scaled_abscissa = [scaled_abscissa; scaled_abscissa_this_band]; %#ok<AGROW>
        band_limits = [band_limits scaled_abscissa(end)]; %#ok<AGROW>
    end
end
