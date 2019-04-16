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

function [data, names, directions] = addOrOverwriteResultsData(data, names, directions, new_data, new_name, new_directions)
    index_in_existing_data = find(strcmp(names, new_name));
    if isempty(index_in_existing_data)
        data = [data; new_data];
        names = [names; new_name];
        directions = [directions; new_directions];
    else
        data{index_in_existing_data} = new_data;
        directions(index_in_existing_data,:) = new_directions;
        % TODO: check whether directions are different
    end
end
