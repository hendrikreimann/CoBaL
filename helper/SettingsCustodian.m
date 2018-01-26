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
        property_not_found_list = {};
        verbose = true;
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
                data_string = [];
                % transform this entry into a string
                if iscell(data)
                    for i_cell = 1 : length(data)
                        this_cell = data{i_cell};
                        % make sure this is a string
                        if isnumeric(this_cell)
                            this_cell_string = num2str(this_cell);
                        else
                            this_cell_string = this_cell;
                        end
                        
                        % append to entry string
                        if isempty(data_string)
                            data_string = this_cell_string;
                        else
                            data_string = [data_string ', ' this_cell_string];
                        end
                    end
                end
                
                % data
                report_string = [];
                if isempty(data_string)
                    report_string = ['Setting not found in file ' this.settings_file ' - ' property_name];
                else
                    if isnumeric(data_string)
                        data_string = num2str(data_string);
                    end
                    report_string = ['Setting not found in file ' this.settings_file ', using default - ' property_name ': ' data_string];
                end
                
                if ~any(strcmp(property_name, this.property_not_found_list))
                    if this.verbose
                        disp(report_string)
                    end
                    this.property_not_found_list = [this.property_not_found_list; property_name];
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
            if strcmp(property_name, 'show_outliers')
                default_data = 1;
            end
            if strcmp(property_name, 'discrete_data_plot_style')
                default_data = 'box';
            end
            if strcmp(property_name, 'apply_forceplate_offset')
                default_data = 0;
            end
            if strcmp(property_name, 'figure_settings_file')
                default_data = 'eventGuiFigureSettings.mat';
            end
            if strcmp(property_name, 'markers_HEAD')
                default_data = {'RFHD', 'LFHD', 'RBHD'};
            end
            if strcmp(property_name, 'markers_TORSO')
                default_data = {'C7', 'CLAV', 'T10'};
            end
            if strcmp(property_name, 'markers_LUPPERARM')
                default_data = {'LSHOULDERCOR', 'LELBOWCOR', 'LELB'};
            end
            if strcmp(property_name, 'markers_RUPPERARM')
                default_data = {'RSHOULDERCOR', 'RELBOWCOR', 'RELB'};
            end
            if strcmp(property_name, 'markers_LFOREARM')
                default_data = {'LFRM', 'LWRA', 'LWRB'};
            end
            if strcmp(property_name, 'markers_RFOREARM')
                default_data = {'RFRM', 'RWRA', 'RWRB'};
            end
            if strcmp(property_name, 'markers_LHAND')
                default_data = {'LWRB', 'LWRA', 'LFIN'};
            end
            if strcmp(property_name, 'markers_RHAND')
                default_data = {'RWRB', 'RWRA', 'RFIN'};
            end
            if strcmp(property_name, 'markers_LTHIGH')
                default_data = {'LTHI', 'LHIPCOR', 'LKNE'};
            end
            if strcmp(property_name, 'markers_RTHIGH')
                default_data = {'RTHI', 'RHIPCOR', 'RKNE'};
            end
            if strcmp(property_name, 'markers_LSHANK')
                default_data = {'LKNEECOR', 'LKNE', 'LANK'};
            end
            if strcmp(property_name, 'markers_RSHANK')
                default_data = {'RKNEECOR', 'RKNE', 'RANK'};
            end
            if strcmp(property_name, 'markers_LFOOT')
                default_data = {'LHEE', 'LTOE', 'LTOEL'};
            end
            if strcmp(property_name, 'markers_RFOOT')
                default_data = {'RHEE', 'RTOE', 'RTOEL'};
            end

            if strcmp(property_name, 'acceptable_number_of_zeros_per_stretch')
                default_data = 5;
            end












            
            
            
            
            
            
        end
    
        function settings_names = getAllSettingsNames(this)
            settings_names = fieldnames(this.settings_struct);
        end
        
    end
end