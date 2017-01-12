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

function [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(fileName)
    % remove file ending if necessary
    file_split = strsplit(fileName, '.');

    % split at underscores
    elements = strsplit(file_split{1}, '_');
    date = elements{1};
    subject_id = elements{2};
    trial_type = elements{3};
    trial_number = str2double(elements{4});
    if length(elements) > 4
        file_type = elements{5};
    else
        file_type = '';
    end

%     underscores = find(fileName == '_' | fileName == '-'); % Location of underscores
%     date = fileName(1:underscores(1)-1); % date
%     subject_id = fileName(underscores(1)+1:underscores(2)-1);
%     trial_type = fileName(underscores(2)+1:underscores(3)-1);
%     trial_number = str2double(fileName(underscores(3)+1:underscores(4)-1));
%     file_type = fileName(underscores(4)+1:length(fileName)-4);
return
