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

function directions_opensim = createDirectionsForOpensimData(labels_opensim, directions_settings)
    directions_map = directions_settings.get('specific_directions');
    generic_directions = directions_settings.get('generic_directions');
    directions_opensim = cell(2, length(labels_opensim));

    i_label = 1;
    while i_label <= length(labels_opensim)
        this_label = labels_opensim{i_label};
        
        % check if this is the start of a spatial coordinate triplet
        if i_label <= length(labels_opensim)-2 && isSpatialLabelTriplet(labels_opensim(i_label : i_label+2))
            directions_opensim(:, i_label : i_label+2) = generic_directions';
            i_label = i_label + 3;
        else
            if any(strcmp(directions_map(:, 1), this_label))
                this_label_index = strcmp(directions_map(:, 1), this_label);
                this_label_directions = directions_map(this_label_index, 2:3);
            else
                this_label_directions = {'direction not specified', 'direction not specified'};
            end
            directions_opensim(:, i_label) = this_label_directions;
            i_label = i_label + 1;
        end
    end
end