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
            'collection_date'; ...
          }
        force_cell_list = ...
          {
            'band_labels'; ...
            'trial_types_to_ignore'; ...
            'analog_data_to_import'; ...
            'emg_data_to_import'; ...
          }
    end
    methods
        function this = SettingsCustodian(settings_file)
            this.settings_file = settings_file;
            this.loadSettingsFile(settings_file);
        end
        
        function loadSettingsFile(this, filename)
            this.settings_struct = struct();

            % open file
            if ~exist(filename, 'file')
                error(['Failed to load "' filename '".'])
            end

            fileID = fopen(filename, 'r');

            % read text to cell
            text_line = fgetl(fileID);
            text_cell = {};
            while ischar(text_line)
                text_cell = [text_cell; text_line]; %#ok<AGROW>
                text_line = fgetl(fileID);
            end
            fclose(fileID);

            % prune lines
            lines_to_prune = false(size(text_cell, 1), 1);
            for i_line = 1 : size(text_cell, 1)
                this_line = text_cell{i_line};

                % remove initial white space
                while ~isempty(this_line) && (this_line(1) == ' ' || double(this_line(1)) == 9)
                    this_line(1) = [];
                end
                text_cell{i_line} = this_line; %#ok<AGROW>

                % remove comments
                if length(this_line) > 1 && strcmp(this_line(1:2), '//')
                    lines_to_prune(i_line) = true;
                end

                % flag lines consisting only of white space
                if all(ismember(this_line, ' ') | double(this_line) == 9)
                    lines_to_prune(i_line) = true;
                end


            end
            text_cell(lines_to_prune) = [];

            % extract data and store in settings struct
            while ~isempty(text_cell)
                text_cell = this.parseNextBlock(text_cell);
            end
        end
        
        function [text_cell] = parseNextBlock(this, text_cell) %#ok<INUSL>
            % get first line of remaining text
            text_line = text_cell{1};

            if (length(text_cell) > 1) && (text_cell{2}(1) == '{')
                % this is the beginning of a block
                line_split = strsplit(text_line, ':');
                variable_name = strrep(line_split{1}, ' ', '_');

                % get data
                block_end_line_index = find(strcmp(text_cell, '}'), 1, 'first');
                variable_data_lines = text_cell(3 : block_end_line_index-1);
                variable_value = {};
                for i_line = 1 : length(variable_data_lines)
                    this_line_text = variable_data_lines{i_line};
                    % remove leading white space
                    while ~isempty(this_line_text) && (this_line_text(1) == ' ' || double(this_line_text(1)) == 9)
                        this_line_text(1) = [];
                    end
                    this_line_text = strrep(this_line_text, ', ', ',');
                    % break up by commas
                    this_line_cell = strsplit(this_line_text, ',');
                    % deal with array entries
                    for i_entry = 1 : length(this_line_cell)
                        entry_to_process = this_line_cell{i_entry};
                        % remove leading spaces
                        while ~isempty(entry_to_process) && (entry_to_process(1) == ' ' || double(entry_to_process(1)) == 9)
                            entry_to_process(1) = [];
                        end
                        % remove spaces at end
                        while ~isempty(entry_to_process) && (entry_to_process(end) == ' ' || double(entry_to_process(1)) == 9)
                            entry_to_process(end) = [];
                        end

                        % are there spaces, but all non-spaces numbers?
                        if any(entry_to_process == ' ') && ~strcmp(entry_to_process, '-') && all(ismember(strrep(entry_to_process, ' ', ''), '0123456789-.'))
                            % treat this as a numerical array
                            this_entry_cell = strsplit(entry_to_process, ' ');
                            this_entry_mat = cell2mat(cellfun(@str2num, this_entry_cell, 'un', 0));
                            this_line_cell{i_entry} = this_entry_mat;
                        else
                            % put back original entry, minus leading white space
                            this_line_cell{i_entry} = entry_to_process;
                        end
                    end
                    variable_value(i_line, :) = this_line_cell; %#ok<AGROW>
                end

                % try to transform this into a double array
                variable_array = zeros(size(variable_value)) * NaN;
                for i_row = 1 : size(variable_value, 1)
                    for i_col = 1 : size(variable_value, 2)
                        entry_to_process = variable_value{i_row, i_col};
                        % are these exclusively numbers (and decimal points and minus signs)?
                        if ~strcmp(entry_to_process, '-') && all(ismember(entry_to_process, '0123456789-.'))
                            variable_array(i_row, i_col) = str2num(entry_to_process); %#ok<ST2NM>
                        end
                    end
                end
                if ~any(isnan(variable_array))
                    variable_value = variable_array;
                end

                % store in struct
                this.settings_struct.(variable_name) = variable_value;

                % remove parsed line
                text_cell = text_cell(block_end_line_index+1:end);

                return
            end

            % parse single entry
            line_split = strsplit(text_line, ':');
            variable_name = strrep(line_split{1}, ' ', '_');
            variable_value_string = line_split{2};
            while ~isempty(variable_value_string) && (variable_value_string(1) == ' ' || double(variable_value_string(1)) == 9)
                variable_value_string(1) = [];
            end
            variable_value_string = strrep(variable_value_string, ', ', ',');
            variable_value_cell = strsplit(variable_value_string, ',');
            if length(variable_value_cell) == 1
                variable_value = variable_value_cell{1};

                % try to transform to a single double
                if ~isempty(str2num(variable_value)) %#ok<ST2NM>
                    variable_value = str2num(variable_value); %#ok<ST2NM,NASGU>
                end
            else
                variable_value = variable_value_cell;

                % try to transform to a double array
                variable_value_array = zeros(size(variable_value)) * NaN;
                for i_entry = 1 : length(variable_value)
                    entry_to_process = variable_value{i_entry};
                    entry_to_process = strrep(entry_to_process, 'step', ''); % remove the word "step" to avoid output when trying to transform a string containing it to a number
                    if ~isempty(str2num(entry_to_process)) %#ok<ST2NM>
                        variable_value_array(i_entry) = str2num(entry_to_process); %#ok<ST2NM>
                    end
                end
                if ~any(isnan(variable_value_array))
                    variable_value = variable_value_array;
                end
            end

            % add to settings
            try
                this.settings_struct.(variable_name) = variable_value;
            catch error
                disp(['Variable name causing error: ' variable_name])
                throw(error)
            end

            % remove parsed line
            text_cell = text_cell(2:end);
        end        
        
        function data = get(this, property_name, optional)
            if nargin < 3
                optional = false;
            end
            
            % set to default
            data = this.getDefaultSetting(property_name);
            used_default = true;
            
            % try reading from settings struct
            if isfield(this.settings_struct, property_name)
                data = this.settings_struct.(property_name);
                used_default = false;
            end

            if used_default && ~optional
                error(['Required setting "' property_name '" missing in file ' this.settings_file])
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
                            data_string = [data_string ', ' this_cell_string]; %#ok<AGROW>
                        end
                    end
                end
                
                % data
                report_string = []; %#ok<NASGU>
                if isempty(data_string)
                    report_string = ['Setting not found in file ' this.settings_file ' - ' property_name];
                else
                    if isnumeric(data_string)
                        data_string = num2str(data_string);
                    end
                    report_string = ['Setting not found in file ' this.settings_file ', using default - ' property_name ': ' data_string];
                end
                
                if ~any(strcmp(property_name, this.property_not_found_list))
                    if this.verbose && ~optional
                        disp(report_string)
                    end
                    this.property_not_found_list = [this.property_not_found_list; property_name];
                end
            end
            
            % force string
            if any(strcmp(this.force_string_list, property_name))
                data_old = data;
                
                % is this an individual number?
                if isnumeric(data_old)
                    data = num2str(data_old);
                end
                % not a single number, but not a cell?
                if ~isnumeric(data_old) && ~iscell(data_old)
                    data_cell = cell(size(data_old));
                    for i_entry = 1 : numel(data_old)
                        this_entry = data_old(i_entry);
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
                if ~iscell(data) & ~isempty(data) %#ok<AND2>
                    data = {data};
                end
            end
            
            % check for emptyness
            if strcmp(data, '~')
                data = [];
            end
        end
        
        function table_data = getTable(this, table_name, optional)
            if nargin < 3
                optional = false;
            end
            
            table_header_label = [table_name '_header'];
            
            table_header = this.get(table_header_label, optional);
            table_body = this.get(table_name, optional);
            
            if isempty(table_header)
                table_data = table;
            else
                if isempty(table_body)
                    variable_types = cell(1, length(table_header));
                    variable_types(:) = {'string'};
                    table_data = table('size', [0, length(table_header)], 'VariableTypes', variable_types, 'VariableNames', table_header);
                else
                    table_data = cell2table(table_body, 'VariableNames', table_header);
                end
            end
            
        end
        
        function answer = isfield(this, property_name)
            answer = isfield(this.settings_struct, property_name);
        end
        function default_data = getDefaultSetting(this, property_name) %#ok<INUSL>
            % general default is empty set
            default_data = [];
            
            % apply specific default values
            if strcmp(property_name, 'data_stretch_padding')
                default_data = 0;
            end
            if strcmp(property_name, 'time_to_nearest_heelstrike_before_trigger_threshold')
                default_data = 0.10;
            end
            if strcmp(property_name, 'time_to_nearest_heelstrike_after_trigger_threshold')
                default_data = 0.3;
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
            if strcmp(property_name, 'show_individual_data')
                default_data = 1;
            end
            if strcmp(property_name, 'show_average_data')
                default_data = 1;
            end
            if strcmp(property_name, 'show_spread_data')
                default_data = 1;
            end
            if strcmp(property_name, 'show_single_data_points')
                default_data = 0;
            end
            
            if strcmp(property_name, 'edge_color')
                default_data = [0.4 0.4 0.4];
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
            if strcmp(property_name, 'qtm_import_mode')
                default_data = 'events';
            end
            if strcmp(property_name, 'trial_types_to_ignore')
                default_data = {'calibration', 'emg'};
            end
            if strcmp(property_name, 'force_plates_to_import')
                default_data = [1, 2];
            end
            if strcmp(property_name, 'emg_cutoff_frequency_low')
            
                default_data = 500;
            end
            
            if strcmp(property_name, 'variables_to_plot_header')
                default_data = {'variable_name', 'variable_type', 'variable_label', 'y_axis_label', 'save_file_string', 'x_axis_lower_limit', 'x_axis_upper_limit', 'y_axis_lower_limit', 'y_axis_upper_limit'};
            end
            if strcmp(property_name, 'variables_to_plot_discrete_header')
                default_data = {'variable_name', 'variable_type', 'variable_label', 'y-axis_label', 'save_file_string', 'x_axis_lower limit', 'x-axis_upper_limit', 'y_axis_lower_limit', 'y_axis_upper_limit'};
            end
            if strcmp(property_name, 'variables_to_plot_continuous_header')
                default_data = {'variable_name', 'variable_type', 'variable_label', 'y_axis_label', 'save_file_string', 'x_axis_lower_limit', 'x_axis_upper_limit', 'y_axis_lower_limit', 'y_axis_upper_limit'};
            end
            if strcmp(property_name, 'convert_to_mm')
                default_data = false;
            end
            if strcmp(property_name, 'inverse_kinematics_source')
                default_data = 'opensim';
            end


            if strcmp(property_name, 'analysis_table')
                default_data = ...
                  { ...
                    'invert by condition','inversion_variables','inversion_variables_header'; ...
                    'select from multiple variables by condition','selection_variables','selection_variables_header'; ...
                    'integrate over time','integration_variables','integration_variables_header'; ...
                    'integrate over time','analysis_variables_from_integration','analysis_variables_from_integration_header'; ...
                    'take value at point given by percentage within band','band_percent_variables','band_percent_variables_header'; ...
                    'take value at point given by percentage within band','analysis_variables_from_band_percent','analysis_variables_from_band_percent_header'; ...
                    'calculate mean over time','analysis_variables_from_mean', 'analysis_variables_from_mean_header'; ...
                    'calculate rms over time','analysis_variables_from_rms', 'analysis_variables_from_rms_header'; ...
                    'take value at point given by absolute time within band', 'analysis_variables_from_time_point', 'analysis_variables_from_time_point_header'; ...
                    'take value at band end', 'analysis_variables_from_band_end', 'analysis_variables_from_band_end_header'; ...
                    'take extremum within whole band', 'analysis_variables_from_extrema', 'analysis_variables_from_extrema_header'; ...
                    'take extremum over range', 'analysis_variables_from_extrema_range', 'analysis_variables_from_extrema_range_header'; ...
                    'combine two variables', 'combination_variables', 'combination_variables_header'; ...
                  };
             end
            if strcmp(property_name, 'analysis_table_header')
                default_data = {'action','settings_table','settings_table_header'};
            end
            if strcmp(property_name, 'inversion_variables_header')
                default_data = {'new_variable_name', 'source_variable_name', 'source_type', 'relevant_condition', 'information_table', 'direction_label_positive', 'direction_label_negative'};
            end
            if strcmp(property_name, 'selection_variables_header')
                default_data = {'new_variable_name','source_type','relevant_condition','information_table'};
            end
            if strcmp(property_name, 'selection_variables_header')
                default_data = {'new_variable_name','source_variable_name','source_variable_type','start','start_variable_type','end','end_variable_type'};
            end
            if strcmp(property_name, 'analysis_variables_from_band_end_header')
                default_data = {'new_variable_name','source_variable_name','source_type'};
            end
            if strcmp(property_name, 'analysis_variables_from_extrema_header')
                default_data = {'new_variable_name','source_variable_name','source_type','extremum_type'};
            end
            if strcmp(property_name, 'analysis_variables_from_extrema_range_header')
                default_data = {'new_variable_name','source_variable_name','source_type','extremum_type'};
            end
            
            if strcmp(property_name, 'emg_import_map_header')
                default_data = {'label_in_qtm_file','label_in_cobal'};
            end
            if strcmp(property_name, 'emg_import_map')
                default_data = {};
            end
            if strcmp(property_name, 'filter_order_belt')
                default_data = 4;
            end
            if strcmp(property_name, 'filter_cutoff_belt')
                default_data = 10;
            end
            if strcmp(property_name, 'labview_resampling_rate')
                default_data = 100;
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