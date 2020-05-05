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

function data = addOrReplaceResultsData(data, new_data, data_label)
    % create the data fields if they don't exist yet
    if ~isfield(data, [data_label '_data_session'])
        data.([data_label '_data_session']) = {};
        data.([data_label '_directions_session']) = {};
        data.([data_label '_names_session']) = {};
    end

    index_in_existing_data = find(strcmp(data.([data_label '_names_session']), new_data.name), 1, 'first');
    if isempty(index_in_existing_data)
        data.([data_label '_data_session']) = [data.([data_label '_data_session']); new_data.data];
        data.([data_label '_names_session']) = [data.([data_label '_names_session']); new_data.name];
        data.([data_label '_directions_session']) = [data.([data_label '_directions_session']); new_data.directions];
    else
        % update data
        data.([data_label '_data_session']){index_in_existing_data} = new_data.data;

        % compare directions and warn if they don't match
        old_directions = data.([data_label '_directions_session'])(index_in_existing_data, :);
        if ~isequal(old_directions, new_data.directions)
            warning(['Updating data for variable ' new_data.name ', but directions do not match.'])
        end
    end
end
