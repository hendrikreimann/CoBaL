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

function existing_data_struct = addOrReplaceResultsData(existing_data_struct, new_data_struct)
    new_data = new_data_struct.data;
    new_name = new_data_struct.name;
    new_directions = new_data_struct.directions;

    index_in_existing_data = find(strcmp(existing_data_struct.names, new_name));
    if isempty(index_in_existing_data)
        existing_data_struct.data = [existing_data_struct.data; new_data];
        existing_data_struct.names = [existing_data_struct.names; new_name];
        existing_data_struct.directions = [existing_data_struct.directions; new_directions];
    else
        existing_data_struct.data{index_in_existing_data} = new_data;
        existing_data_struct.directions(index_in_existing_data,:) = new_directions;
        % TODO: check whether directions are different
    end
end
