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

function clarified_string = conditionStringToFilename(condition_string)
    proxy_string = condition_string;
    proxy_string = strrep(proxy_string, 'STANCE_LEFT', 'stanceL');
    proxy_string = strrep(proxy_string, 'STANCE_RIGHT', 'stanceR');
    proxy_string = strrep(proxy_string, 'ILLUSION_LEFT', 'illuL');
    proxy_string = strrep(proxy_string, 'ILLUSION_RIGHT', 'illuR');
    proxy_string = strrep(proxy_string, 'ONE', 'step1');
    proxy_string = strrep(proxy_string, 'TWO', 'step2');
    proxy_string = strrep(proxy_string, 'THREE', 'step3');
    proxy_string = strrep(proxy_string, 'FOUR', 'step4');
    clarified_string = proxy_string;
end