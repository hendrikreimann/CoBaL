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

function [date, subject_id, trial_type, trial_number, file_type, success] = getFileParameters(fileName)
    % remove file ending if necessary
    file_split = strsplit(fileName, '.');

    % initialize
    date = [];
    subject_id = [];
    trial_type = [];
    trial_number = [];
    file_type = '';
    success = false;
    
    % split at underscores
    elements = strsplit(file_split{1}, '_');
    if length(elements) > 3
        date = elements{1};
        subject_id = elements{2};
        trial_type = elements{3};
        trial_number = str2double(elements{4});
        success = true;
    else
        
    end
    
    if length(elements) > 4
        file_type = elements{5};
    else
        file_type = '';
    end

end
