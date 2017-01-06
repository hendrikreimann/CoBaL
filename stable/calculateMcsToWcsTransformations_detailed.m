
function transformations = calculateMcsToWcsTransformations_detailed(marker_positions, marker_headers, segment_labels)

    
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
            number_of_markers_this_segment = 3;
            
        else % define three markers and use these to define the basis
            if strcmp(this_segment_label, 'HEAD')
                this_segment_markers = {'RFHD', 'LFHD', 'RBHD'};
            elseif strcmp(this_segment_label, 'TORSO')
                this_segment_markers = {'C7', 'CLAV', 'T10'};
            elseif strcmp(this_segment_label, 'LUPPERARM')
                this_segment_markers = {'LSHOULDERCOR', 'LELBOWCOR', 'LELB'};
            elseif strcmp(this_segment_label, 'RUPPERARM')
                this_segment_markers = {'RSHOULDERCOR', 'RELBOWCOR', 'RELB'};
            elseif strcmp(this_segment_label, 'LFOREARM')
                this_segment_markers = {'LELB', 'LWRA', 'LWRB'};
            elseif strcmp(this_segment_label, 'RFOREARM')
                this_segment_markers = {'RELB', 'RWRA', 'RWRB'};
            elseif strcmp(this_segment_label, 'LHAND')
                this_segment_markers = {'LWRA', 'LWRB', 'LFIN'};
            elseif strcmp(this_segment_label, 'RHAND')
                this_segment_markers = {'RWRA', 'RWRB', 'RFIN'};
            elseif strcmp(this_segment_label, 'LTHIGH')
                this_segment_markers = {'LHIPCOR', 'LKNEECOR', 'LKNE'};
            elseif strcmp(this_segment_label, 'RTHIGH')
                this_segment_markers = {'RHIPCOR', 'RKNEECOR', 'RKNE'};
            elseif strcmp(this_segment_label, 'LSHANK')
                this_segment_markers = {'LKNEECOR', 'LKNE', 'LANK'};
            elseif strcmp(this_segment_label, 'RSHANK')
                this_segment_markers = {'RKNEECOR', 'RKNE', 'RANK'};
            elseif strcmp(this_segment_label, 'LFOOT')
                this_segment_markers = {'LANK', 'LHEE', 'LTOE'};
            elseif strcmp(this_segment_label, 'RFOOT')
                this_segment_markers = {'RANK', 'RHEE', 'RTOE'};
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