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

function settings = loadSettingsFromFile(source)
    % determine file name
    if isequal(source(end-3:end), '.txt')
        % source was specified as a file name, so just use that
        settings_file_name = source;
    else
        % source was specified as a general type, so use default
        settings_file_name = [source 'Settings.txt'];
    end

    % look for specified file, either here or one folder up
    if exist(settings_file_name, 'file')
        settings_file = settings_file_name;
    end    
    if exist(['..' filesep settings_file_name], 'file')
        settings_file = ['..' filesep settings_file_name];
    end
    
    % load settings from this file
    settings = SettingsCustodian(settings_file);
end

