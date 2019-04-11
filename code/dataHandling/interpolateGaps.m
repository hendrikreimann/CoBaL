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

function filled_trajectory = interpolateGaps(trajectory, maximal_gap_number_of_indices)
    number_of_time_steps = length(trajectory);
    if any(isnan(trajectory))
        % determine start of first fillable gap
        this_index = 1;
        while (this_index < number_of_time_steps) && any(isnan(trajectory(this_index : end)))
            % find start and end of next gap
            this_gap_start_index = this_index + find(isnan(trajectory(this_index : end)), 1) - 1;
            this_gap_end_index = this_gap_start_index;
            while this_gap_end_index < number_of_time_steps && isnan(trajectory(this_gap_end_index+1))
                this_gap_end_index = this_gap_end_index+1;
            end

            % fill gap if not too long and not at the start or end
            gap_size = this_gap_end_index - this_gap_start_index + 1;
            if (gap_size < maximal_gap_number_of_indices) && (this_gap_end_index < number_of_time_steps) && (this_gap_start_index > 1)
                trajectory(this_gap_start_index-1 : this_gap_end_index+1) = ...
                    linspace(trajectory(this_gap_start_index-1), trajectory(this_gap_end_index+1), gap_size+2);
            end
            
            % move on
            this_index = this_gap_end_index + 1;
        end
    end
    filled_trajectory = trajectory;
end
