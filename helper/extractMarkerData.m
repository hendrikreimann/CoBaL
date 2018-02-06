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

function data = extractMarkerData(marker_trajectories, marker_labels, specified_marker, output)
    if nargin < 4
        output = 'trajectories';
    end
    index_x = find(strcmp(marker_labels, [specified_marker '_x']));
    index_y = find(strcmp(marker_labels, [specified_marker '_y']));
    index_z = find(strcmp(marker_labels, [specified_marker '_z']));
    if strcmp(output, 'indices')
        data = [index_x index_y index_z];
    end
    if strcmp(output, 'trajectories')
        trajectory_x = marker_trajectories(:, index_x);
        trajectory_y = marker_trajectories(:, index_y);
        trajectory_z = marker_trajectories(:, index_z);
        data = [trajectory_x trajectory_y trajectory_z];
    end
end