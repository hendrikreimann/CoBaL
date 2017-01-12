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

function clarified_string = clarifyConditionString(condition_string)
    proxy_string = condition_string;
    proxy_string = strrep(proxy_string, 'LEFT', 'stance foot LEFT');
    proxy_string = strrep(proxy_string, 'RIGHT', 'stance foot RIGHT');
    proxy_string = strrep(proxy_string, 'POSITIVE', 'illusion LEFT');
    proxy_string = strrep(proxy_string, 'NEGATIVE', 'illusion RIGHT');
    clarified_string = proxy_string;
end