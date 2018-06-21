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
        
        force_string_list = ...
          { ...
            'preferred_level_order'; ...
          }
        force_cell_list = ...
          {
            'band_labels'; ...
          }
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
            
            % try reading from settings struct
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
            
            % force string
            if any(strcmp(this.force_string_list, property_name))
                if ~iscell(data)
                    data_cell = cell(size(data));
                    for i_entry = 1 : numel(data)
                        this_entry = data(i_entry);
                        if isnumeric(this_entry)
                            data_cell{i_entry} = num2str(this_entry);
                        else
                            data_cell{i_entry} = this_entry;
                        end
                        
                    end
                    data = data_cell;
                end
            end
            
            % force cell
            if any(strcmp(this.force_cell_list, property_name))
                if ~iscell(data)
                    data = {data};
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
            if strcmp(property_name, 'mark_pushoff')
                default_data = 0;
            end
            if strcmp(property_name, 'mark_bands')
                default_data = 0;
            end
            if strcmp(property_name, 'merge_bands')
                default_data = 0;
            end
            if strcmp(property_name, 'xtick_label_rotation')
                default_data = 0;
            end
            if strcmp(property_name, 'group_bands_within_conditions')
                default_data = 0;
            end
            if strcmp(property_name, 'emg_time_offset')
                default_data = 0;
            end
            if strcmp(property_name, 'trial_types_to_ignore')
                default_data = {'calibration', 'emg'};
            end
            
            
            
            
            if strcmp(property_name, 'marker_to_segment_map')
                default_data = ...
                  { ...
                    'LFHD', '26'; ...
                    'RFHD', '26'; ...
                    'LBHD', '26'; ...
                    'RBHD', '26'; ...
                    'C7', '23'; ...
                    'T10', '23'; ...
                    'CLAV', '23'; ...
                    'STRN', '23'; ...
                    'RBAK', '23'; ...
                    'LSHO', '23'; ...
                    'LUPA', '29'; ...
                    'LELB', '29'; ...
                    'LFRM', '31'; ...
                    'LWRA', '31'; ...
                    'LWRB', '31'; ...
                    'LFIN', '32'; ...
                    'RSHO', '23'; ...
                    'RUPA', '35'; ...
                    'RELB', '35'; ...
                    'RFRM', '37'; ...
                    'RWRA', '37'; ...
                    'RWRB', '37'; ...
                    'RFIN', '38'; ...
                    'LASI', '6'; ...
                    'RASI', '6'; ...
                    'LPSI', '6'; ...
                    'RPSI', '6'; ...
                    'LTHI', '9'; ...
                    'LTHIA', '9'; ...
                    'LKNE', '9'; ...
                    'LTIB', '11'; ...
                    'LTIBA', '11'; ...
                    'LANK', '11'; ...
                    'LHEE', '13'; ...
                    'LTOE', '13'; ...
                    'LTOEL', '13'; ...
                    'RTHI', '16'; ...
                    'RTHIA', '16'; ...
                    'RKNE', '16'; ...
                    'RTIB', '18'; ...
                    'RTIBA', '18'; ...
                    'RANK', '18'; ...
                    'RHEE', '20'; ...
                    'RTOE', '20'; ...
                    'RTOEL', '20'; ...
                  };
            end
            if strcmp(property_name, 'marker_to_color_map')
                default_data = ...
                  { ...
                    'LFHD', '1', '0', '0'; ...
                    'RFHD', '0', '1', '0'; ...
                    'LBHD', '1', '0', '0'; ...
                    'RBHD', '0', '1', '0'; ...
                    'C7', '0', '0', '1'; ...
                    'T10', '0', '0', '1'; ...
                    'CLAV', '0', '0', '1'; ...
                    'STRN', '0', '0', '1'; ...
                    'RBAK', '0', '0', '1'; ...
                    'LSHO', '1', '0', '0'; ...
                    'LUPA', '1', '0', '0'; ...
                    'LELB', '1', '0', '0'; ...
                    'LFRM', '1', '0', '0'; ...
                    'LWRA', '1', '0', '0'; ...
                    'LWRB', '1', '0', '0'; ...
                    'LFIN', '1', '0', '0'; ...
                    'RSHO', '0', '1', '0'; ...
                    'RUPA', '0', '1', '0'; ...
                    'RELB', '0', '1', '0'; ...
                    'RFRM', '0', '1', '0'; ...
                    'RWRA', '0', '1', '0'; ...
                    'RWRB', '0', '1', '0'; ...
                    'RFIN', '0', '1', '0'; ...
                    'LASI', '1', '0', '0'; ...
                    'RASI', '0', '1', '0'; ...
                    'LPSI', '1', '0', '0'; ...
                    'RPSI', '0', '1', '0'; ...
                    'LTHI', '1', '0', '0'; ...
                    'LTHIA', '1', '0', '0'; ...
                    'LKNE', '1', '0', '0'; ...
                    'LTIB', '1', '0', '0'; ...
                    'LTIBA', '1', '0', '0'; ...
                    'LANK', '1', '0', '0'; ...
                    'LHEE', '1', '0', '0'; ...
                    'LTOE', '1', '0', '0'; ...
                    'LTOEL', '1', '0', '0'; ...
                    'RTHI', '0', '1', '0'; ...
                    'RTHIA', '0', '1', '0'; ...
                    'RKNE', '0', '1', '0'; ...
                    'RTIB', '0', '1', '0'; ...
                    'RTIBA', '0', '1', '0'; ...
                    'RANK', '0', '1', '0'; ...
                    'RHEE', '0', '1', '0'; ...
                    'RTOE', '0', '1', '0'; ...
                    'RTOEL', '0', '1', '0'; ...
                  };
            end






            
            
            
            
            
            
        end
    
        function settings_names = getAllSettingsNames(this)
            settings_names = fieldnames(this.settings_struct);
        end
        
    end
end