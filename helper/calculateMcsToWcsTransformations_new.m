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


function transformations = calculateMcsToWcsTransformations_new(marker_positions, marker_headers, segment_labels, subject_settings)

    
    transformations = cell(length(segment_labels), 1);
    for i_segment = 1 : length(segment_labels)
        this_segment_label = segment_labels{i_segment};

        if strcmp(this_segment_label, 'PELVIS')
            LASI = extractMarkerTrajectories(marker_positions, marker_headers, 'LASI')';
            RASI = extractMarkerTrajectories(marker_positions, marker_headers, 'RASI')';
            LPSI = extractMarkerTrajectories(marker_positions, marker_headers, 'LPSI')';
            RPSI = extractMarkerTrajectories(marker_positions, marker_headers, 'RPSI')';
            
            pelvis_base_point = mean([LASI, RASI, LPSI, RPSI], 2);
            pelvis_left = (LASI + LPSI) / 2;
            pelvis_right = (RASI + RPSI) / 2;
            pelvis_left_to_right = pelvis_right - pelvis_left;
            MASIS = (LASI + RASI) / 2;
            MPSIS = (LPSI + RPSI) / 2;
            MPSIS_to_MASIS = MASIS - MPSIS;
            
            positions_this_segment = [pelvis_base_point, pelvis_base_point+pelvis_left_to_right, pelvis_base_point+MPSIS_to_MASIS]';
            
        else % define three markers and use these to define the basis
            this_segment_markers = subject_settings.get(['markers_' this_segment_label]);

            
%             if strcmp(this_segment_label, 'HEAD')
%                 this_segment_markers = {'RFHD', 'LFHD', 'RBHD'};
%             elseif strcmp(this_segment_label, 'TORSO')
%                 this_segment_markers = {'C7', 'CLAV', 'T10'};
%             elseif strcmp(this_segment_label, 'LUPPERARM')
%                 this_segment_markers = {'LSHOULDERCOR', 'LELBOWCOR', 'LELB'};
%             elseif strcmp(this_segment_label, 'RUPPERARM')
%                 this_segment_markers = {'RSHOULDERCOR', 'RELBOWCOR', 'RELB'};
            if strcmp(this_segment_label, 'LFOREARM')
                this_segment_markers = {'LFRA', 'LWRA', 'LWRB'};
            elseif strcmp(this_segment_label, 'RFOREARM')
                this_segment_markers = {'RFRA', 'RWRA', 'RWRB'};
%             elseif strcmp(this_segment_label, 'LHAND')
%                 this_segment_markers = {'LWRB', 'LWRA', 'LFIN'}; % this is not good because these markers are not actually on a rigid body. Check later if I can change that
%             elseif strcmp(this_segment_label, 'RHAND')
%                 this_segment_markers = {'RWRB', 'RWRA', 'RFIN'};
%             elseif strcmp(this_segment_label, 'LTHIGH')
%                 this_segment_markers = {'LTHI', 'LHIPCOR', 'LKNE'};
%             elseif strcmp(this_segment_label, 'RTHIGH')
%                 this_segment_markers = {'RTHI', 'RHIPCOR', 'RKNE'};
%                 
%             elseif strcmp(this_segment_label, 'LSHANK')
%                 this_segment_markers = {'LKNEECOR', 'LKNE', 'LANK'};
%             elseif strcmp(this_segment_label, 'RSHANK')
%                 this_segment_markers = {'RKNEECOR', 'RKNE', 'RANK'};
%             elseif strcmp(this_segment_label, 'LFOOT')
% %                 this_segment_markers = {'LANK', 'LHEE', 'LTOE'};
%                 this_segment_markers = {'LHEE', 'LTOE', 'LTOEL'};
%             elseif strcmp(this_segment_label, 'RFOOT')
% %                 this_segment_markers = {'RANK', 'RHEE', 'RTOE'};
%                 this_segment_markers = {'RHEE', 'RTOE', 'RTOEL'};
            end                
            
            % find marker numbers
            number_of_markers_this_segment = length(this_segment_markers);
            marker_numbers = zeros(1, number_of_markers_this_segment);
            for i_marker = 1 : number_of_markers_this_segment
                marker_numbers(i_marker) = find(strcmp(marker_headers, this_segment_markers(i_marker)));
            end

            % find marker indices
            marker_indices_x = zeros(1, length(marker_numbers));
            marker_indices_y = zeros(1, length(marker_numbers));
            marker_indices_z = zeros(1, length(marker_numbers));
            for i_marker = 1 : number_of_markers_this_segment
                marker_number = marker_numbers(i_marker);
                marker_indices_x(i_marker) = (marker_number - 1) * 3 + 1;
                marker_indices_y(i_marker) = (marker_number - 1) * 3 + 2;
                marker_indices_z(i_marker) = (marker_number - 1) * 3 + 3;
            end

            % extract marker positions
            positions_this_segment = zeros(number_of_markers_this_segment, 3);
            for i_marker = 1 : number_of_markers_this_segment
                positions_this_segment(i_marker, :) = marker_positions([marker_indices_x(i_marker) marker_indices_y(i_marker) marker_indices_z(i_marker)]);
            end

        end
        
        % define coordinate frame based on positions
        p = positions_this_segment(1, :)';
        u1 = positions_this_segment(2, :)' - p;
        u2 = positions_this_segment(3, :)' - p;
        u3 = cross(u1, u2);
        R = orthogonalizeBasis([u1 u2 u3]);
        
        % define transformation
        transformations{i_segment} = [R p; [0 0 0 1]];

        
    end



end