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

classdef SettingsCustodian < handle
    properties
        settings_file = '';
        settings_struct = struct;
    end
    methods
        function this = SettingsCustodian(settings_file)
            this.settings_file = settings_file;
            this.settings_struct = loadSettingsFile(settings_file);
            
        end
        function data = get(this, property_name)
            if isfield(this.settings_struct, property_name)
                eval(['data = this.settings_struct.' property_name ';']);
            else
                data = [];
            end
        end
    
        function settings_names = getAllSettingsNames(this)
            settings_names = fieldnames(this.settings_struct);
        end
        
    end
end