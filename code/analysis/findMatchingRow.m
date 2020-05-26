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

function row_index = findMatchingRow(cell_array, row_to_match)
    matches = false(size(cell_array));
    for i_column = 1 : size(cell_array, 2)
        entry_to_match = row_to_match{i_column};
        match_this_factor = strcmp(cell_array(:, i_column), entry_to_match);
        matches(:, i_column) = match_this_factor;
    end
    row_index = find(all(matches, 2));

end
