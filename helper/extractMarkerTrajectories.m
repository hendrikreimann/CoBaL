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

function extracted_trajectory = extractMarkerTrajectories(marker_trajectories, marker_labels, specified_marker)
    warning('This function is depcrecated. Use "extractMarkerData()" instead');

    marker_number = find(strcmp(marker_labels, specified_marker));
    markers_indices = reshape([(marker_number - 1) * 3 + 1; (marker_number - 1) * 3 + 2; (marker_number - 1) * 3 + 3], 1, length(marker_number)*3);
    extracted_trajectory = marker_trajectories(:, markers_indices);
end