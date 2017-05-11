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
            % set to default
            data = this.getDefaultSetting(property_name);
            used_default = true;
            
            % try reading from settings file
            if isfield(this.settings_struct, property_name)
                eval(['data = this.settings_struct.' property_name ';']);
                used_default = false;
            end
            
            
            
            if used_default
                data_string = data;
                if isempty(data_string)
                    disp(['Setting not found in file ' this.settings_file ' - ' property_name]);
                else
                    if isnumeric(data_string)
                        data_string = num2str(data_string);
                    end
                    disp(['Setting not found in file ' this.settings_file ', using default - ' property_name ': ' data_string]);
                end
            end
        end
        
        function default_data = getDefaultSetting(this, property_name)
            % general default is empty set
            default_data = [];
            
            % apply specific default values
            if strcmp(property_name, 'data_stretch_padding')
                default_data = 0;
            end
            if strcmp(property_name, 'left_fullstance_method')
                default_data = 'first_zero_crossing_after_heelstrike';
            end
            if strcmp(property_name, 'right_fullstance_method')
                default_data = 'first_zero_crossing_after_heelstrike';
            end
            if strcmp(property_name, 'left_armswing_peak_prominence_threshold')
                default_data = 10;
            end
            if strcmp(property_name, 'left_armswing_peak_distance_threshold')
                default_data = 0.25;
            end
            if strcmp(property_name, 'right_armswing_peak_prominence_threshold')
                default_data = 10;
            end
            if strcmp(property_name, 'right_armswing_peak_distance_threshold')
                default_data = 0.25;
            end
            if strcmp(property_name, 'left_legswing_peak_prominence_threshold')
                default_data = 10;
            end
            if strcmp(property_name, 'left_legswing_peak_distance_threshold')
                default_data = 0.25;
            end
            if strcmp(property_name, 'right_legswing_peak_prominence_threshold')
                default_data = 10;
            end
            if strcmp(property_name, 'right_legswing_peak_distance_threshold')
                default_data = 0.25;
            end
            if strcmp(property_name, 'plot_control')
                default_data = 1;
            end
            if strcmp(property_name, 'plot_response')
                default_data = 0;
            end
            if strcmp(property_name, 'apply_forceplate_offset')
                default_data = 0;
            end
            if strcmp(property_name, 'figure_settings_file')
                default_data = 'eventGuiFigureSettings.mat';
            end
            

            
            
            
            
            
            
        end
    
        function settings_names = getAllSettingsNames(this)
            settings_names = fieldnames(this.settings_struct);
        end
        
    end
end