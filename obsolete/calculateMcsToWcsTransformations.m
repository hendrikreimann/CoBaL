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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% compare the kinematic tree against the kinematic chain


function transformations = calculateMcsToWcsTransformations(marker_positions, marker_headers, markers_by_segment)

    % find marker numbers
    number_of_segments = size(markers_by_segment, 1);
    number_of_markers_per_segment = size(markers_by_segment, 2);
    marker_numbers = zeros(size(markers_by_segment));
    for i_segment = 1 : number_of_segments
        for i_marker = 1 : number_of_markers_per_segment
            marker_numbers(i_segment, i_marker) = find(strcmp(marker_headers, markers_by_segment(i_segment, i_marker)));
        end
    end

    % find marker indices
    marker_indices_x = zeros(size(marker_numbers));
    marker_indices_y = zeros(size(marker_numbers));
    marker_indices_z = zeros(size(marker_numbers));
    for i_segment = 1 : number_of_segments
        for i_marker = 1 : number_of_markers_per_segment
            marker_number = marker_numbers(i_segment, i_marker);
            marker_indices_x(i_segment, i_marker) = (marker_number - 1) * 3 + 1;
            marker_indices_y(i_segment, i_marker) = (marker_number - 1) * 3 + 2;
            marker_indices_z(i_segment, i_marker) = (marker_number - 1) * 3 + 3;
        end
    end
    
    % calculate transformations
    transformations = cell(number_of_segments, 1);
    for i_segment = 1 : number_of_segments
        % extract marker positions
        marker_positions_this_segment = zeros(number_of_markers_per_segment, 3);
        for i_marker = 1 : number_of_markers_per_segment
            marker_positions_this_segment(i_marker, :) = marker_positions([marker_indices_x(i_segment, i_marker) marker_indices_y(i_segment, i_marker) marker_indices_z(i_segment, i_marker)]);
        end
   
        % define coordinate frame
        p = marker_positions_this_segment(1, :)';
        u1 = marker_positions_this_segment(2, :)' - p;
        u2 = marker_positions_this_segment(3, :)' - p;
        u3 = cross(u1, u2);
        
        R = orthonormalizeBasis([u1, u2, u3]);
        
        A = [u1 u2 u3];
        Q = zeros(3, number_of_markers_per_segment);
        R_qr = zeros(size(number_of_markers_per_segment, number_of_markers_per_segment));
        for i_col = 1 : number_of_markers_per_segment
            v = A(:, i_col);
            for i_row = 1 : i_col-1
                R_qr(i_row, i_col) = Q(:, i_row)' * A(:, i_col);
                v = v - R_qr(i_row, i_col) * Q(:, i_row);
            end
            R_qr(i_col, i_col) = norm(v);
            Q(:, i_col) = v/R_qr(i_col, i_col);
        end
        R = Q;
        
        
        transformations{i_segment} = [R p; [0 0 0 1]];
        
    end    
end