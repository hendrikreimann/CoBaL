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

function fileName = makeFileName(date, subject_id, trial_type, trial_number, file_type)
    % transform to required format if necessary
    if isstring(trial_type)
        trial_type = char(trial_type);
    end

    % catch funny case
    if trial_type(1) == '.'
        fileName = [date '_' subject_id trial_type];
    else
        fileName = [date '_' subject_id '_' trial_type];
    end
    
    if nargin > 3
        if ~isnumeric(trial_number)
            trial_number = str2num(trial_number);
        end
        trial_number = zeroPrefixedIntegerString(trial_number, 3);
        fileName = [fileName '_' trial_number];
    end
    if nargin > 4
        fileName = [fileName '_' file_type];
    end
return
