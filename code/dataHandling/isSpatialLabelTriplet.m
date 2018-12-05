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

function result = isSpatialLabelTriplet(labels)
    label_1 = labels{1};
    label_2 = labels{2};
    label_3 = labels{3};
    
    result = false;
    if strcmp(label_1(1:end-1), label_2(1:end-1)) && strcmp(label_1(1:end-1), label_3(1:end-1)) && strcmp(label_2(1:end-1), label_3(1:end-1))
        if label_1(end) == 'x' && label_2(end) == 'y' && label_3(end) == 'z'
            result = true;
        end
        if label_1(end) == 'X' && label_2(end) == 'Y' && label_3(end) == 'Z'
            result = true;
        end
    end
end