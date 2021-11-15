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

% this function transforms larger blocks of data, such as marker coordinates, between coordinate frames, keeping track
% of the labels and directions

function [transformed_trajectories, transformed_labels, transformed_directions] = transformSpatialData(source_trajectories, source_labels, source_directions, rigid_transformation, study_settings)
    transformed_trajectories = zeros(size(source_trajectories));
    transformed_labels = cell(size(source_labels));
    transformed_directions = cell(size(source_directions));
    
    
    i_trajectory = 1;
    while i_trajectory <= length(source_labels)
        if i_trajectory <= length(source_labels)-2 && isSpatialLabelTriplet(source_labels(i_trajectory : i_trajectory+2))
            % this is the start of a spatial coordinate triplet, so transform
            this_spatial_trajectory_source = source_trajectories(:, i_trajectory : i_trajectory+2);
            this_label_source = source_labels(:, i_trajectory : i_trajectory+2);
            this_direction_source = source_directions(:, i_trajectory : i_trajectory+2);
            
            % transform trajectories
            this_spatial_trajectory_transformed = (eye(3, 4) * rigid_transformation * [this_spatial_trajectory_source'; ones(1, size(this_spatial_trajectory_source, 1))])';
            transformed_trajectories(:, i_trajectory : i_trajectory+2) = this_spatial_trajectory_transformed;
            
            % transform labels
            label_root = this_label_source{1}(1:end-1);
            transformed_labels{i_trajectory} = [label_root 'x'];
            transformed_labels{i_trajectory+1} = [label_root 'y'];
            transformed_labels{i_trajectory+2} = [label_root 'z'];
            
            % transform directions
            if any(any(~((rigid_transformation(1:3, 1:3) == 1) | (rigid_transformation(1:3, 1:3) == 0))))
                transformed_directions(:, i_trajectory) = {'TBD', 'TBD'};
                transformed_directions(:, i_trajectory+1) = {'TBD', 'TBD'};
                transformed_directions(:, i_trajectory+2) = {'TBD', 'TBD'};
                warning(['Directions for spatial variable ' label_root ' could not be transformed automatically, using placeholder'])
            else
                transformed_directions(:, i_trajectory) = this_direction_source(:, logical(rigid_transformation(1, :)));
                transformed_directions(:, i_trajectory+1) = this_direction_source(:, logical(rigid_transformation(2, :)));
                transformed_directions(:, i_trajectory+2) = this_direction_source(:, logical(rigid_transformation(3, :)));
            end            
            i_trajectory = i_trajectory + 3;
        else
            % simply copy over
            transformed_trajectories(:, i_trajectory) = source_trajectories(:, i_trajectory);
            transformed_labels(i_trajectory) = source_labels(i_trajectory);
            transformed_directions(:, i_trajectory) = source_directions(:, i_trajectory);
            i_trajectory = i_trajectory + 1;
        end
    end

    % go through directions and direction of cartesian basis vectors with labels from study settings
    x_pos_label = study_settings.get('direction_x_pos');
    x_neg_label = study_settings.get('direction_x_neg');
    y_pos_label = study_settings.get('direction_y_pos');
    y_neg_label = study_settings.get('direction_y_neg');
    z_pos_label = study_settings.get('direction_z_pos');
    z_neg_label = study_settings.get('direction_z_neg');
    
    for i_column = 1 : size(transformed_directions, 2)
        this_direction_pos = transformed_directions{1, i_column};
        this_direction_neg = transformed_directions{2, i_column};
        if strcmp(this_direction_pos, 'x+') && strcmp(this_direction_neg, 'x-')
            transformed_directions{1, i_column} = x_pos_label;
            transformed_directions{2, i_column} = x_neg_label;
        end
        if strcmp(this_direction_pos, 'y+') && strcmp(this_direction_neg, 'y-')
            transformed_directions{1, i_column} = y_pos_label;
            transformed_directions{2, i_column} = y_neg_label;
        end
        if strcmp(this_direction_pos, 'z+') && strcmp(this_direction_neg, 'z-')
            transformed_directions{1, i_column} = z_pos_label;
            transformed_directions{2, i_column} = z_neg_label;
        end
        
    end
    
    
end