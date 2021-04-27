%     This file is part of the CoBaL code base
%     Copyright (C) 2020 Hendrik Reimann <hendrikreimann@gmail.com>
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

classdef stretchGui < handle
    properties
        data_custodian;
        settings_custodian;
        control_figure;
        figure_settings;
        figures_list;
        axes_list;
        line_handle_list;
        point_handle_list;
        number_of_figures;
        
        current_stretch;
        current_time;
        
        % controls
        current_stretch_edit
        time_slider
        
        canvas_on_screen_width;
        canvas_on_screen_height;
    end
    
    methods
        function this = stretchGui(varargin)
            % parse arguments
            parser = inputParser;
            parser.KeepUnmatched = true;
            addParameter(parser, 'settings', 'stretchGui')
            addParameter(parser, 'root', pwd)
            addParameter(parser, 'source', '')
            addParameter(parser, 'subjects', [])
            parse(parser, varargin{:})
            settings_file = parser.Results.settings;
            root = parser.Results.root;
            source= parser.Results.source;
            subjects = parser.Results.subjects;
            
            % load stuff
            this.settings_custodian = loadSettingsFromFile(settings_file);
            this.parseFigureSettings();
            this.data_custodian = StretchDataCustodian(root, source, subjects);
            
            % create figures
            this.createFigures();
            this.createController();
            this.setStretch(1);
        end
        function parseFigureSettings(this)
            figure_header = this.settings_custodian.get('figures_header');
            figures_table = this.settings_custodian.get('figures');

            this.number_of_figures = size(figures_table, 1);
            this.figure_settings = cell(this.number_of_figures, 1);
            for i_figure = 1 : this.number_of_figures
                figure_settings_here = struct;
                % store data
                figure_settings_here.data = this.settings_custodian.get(figures_table{i_figure, strcmp(figure_header, 'data')});
                figure_settings_here.data_header = this.settings_custodian.get(figures_table{i_figure, strcmp(figure_header, 'data header')});
                % store individual entries in info table
                info = this.settings_custodian.get(figures_table{i_figure, strcmp(figure_header, 'info')});
                for i_entry = 1 : size(info, 1)
                    figure_settings_here.(info{i_entry, 1}) = info{i_entry, 2};
                end
                % go through expected entries and enter defaults if needed
                if ~isfield(figure_settings_here, 'title')
                    figure_settings_here.title = 'untitled figure';
                end
                if ~isfield(figure_settings_here, 'origin_horizontal')
                    figure_settings_here.origin_horizontal = 0.25;
                end
                if ~isfield(figure_settings_here, 'origin_vertical')
                    figure_settings_here.origin_vertical = 0.25;
                end
                if ~isfield(figure_settings_here, 'width')
                    figure_settings_here.width = 0.5;
                end
                if ~isfield(figure_settings_here, 'height')
                    figure_settings_here.height = 0.5;
                end

                % transform entries to numbers if needed
                if ~isnumeric(figure_settings_here.origin_horizontal)
                    figure_settings_here.origin_horizontal = str2double(figure_settings_here.origin_horizontal);
                end
                if ~isnumeric(figure_settings_here.origin_vertical)
                    figure_settings_here.origin_vertical = str2double(figure_settings_here.origin_vertical);
                end
                if ~isnumeric(figure_settings_here.width)
                    figure_settings_here.width = str2double(figure_settings_here.width);
                end
                if ~isnumeric(figure_settings_here.height)
                    figure_settings_here.height = str2double(figure_settings_here.height);
                end
                if isfield(figure_settings_here, 'x_lim_lower') && ~isnumeric(figure_settings_here.x_lim_lower)
                    figure_settings_here.x_lim_lower = str2double(figure_settings_here.x_lim_lower);
                end
                if isfield(figure_settings_here, 'x_lim_upper') && ~isnumeric(figure_settings_here.x_lim_upper)
                    figure_settings_here.x_lim_upper = str2double(figure_settings_here.x_lim_upper);
                end
                if isfield(figure_settings_here, 'y_lim_lower') && ~isnumeric(figure_settings_here.y_lim_lower)
                    figure_settings_here.y_lim_lower = str2double(figure_settings_here.y_lim_lower);
                end
                if isfield(figure_settings_here, 'y_lim_upper') && ~isnumeric(figure_settings_here.y_lim_upper)
                    figure_settings_here.y_lim_upper = str2double(figure_settings_here.y_lim_upper);
                end


                this.figure_settings{i_figure} = figure_settings_here;
            end
        end
        function createController(this)
            screen_size = get(0,'ScreenSize');
            figure_height = 600;
            figure_width = 420;
%             this.canvas_on_screen_width = screen_size(3) - figure_width;
%             this.canvas_on_screen_height = screen_size(4) - 24; % 24 is the size of the menu bar on the Mac, check for windows systems later
            this.control_figure = uifigure('position', [screen_size(3)-figure_width screen_size(3)-figure_height figure_width figure_height], 'Units', 'pixels', 'KeyPressFcn', @this.processKeyPress);
            
            stretch_panel_height = 60;
            stretch_panel = uipanel(this.control_figure, 'Title', 'Stretch', 'FontSize', 12, 'BackgroundColor', 'white', 'Units', 'pixels', 'Position', [5, figure_height-stretch_panel_height-5, figure_width-10, stretch_panel_height]);
            uilabel(stretch_panel, 'text', 'Selected stretch:', 'Position', [5, stretch_panel_height-40, 180, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
            this.current_stretch_edit = uieditfield(stretch_panel, 'BackgroundColor', 'white', 'Position', [90, stretch_panel_height-40, 40, 20], 'Value', '0');
            this.time_slider = uislider(stretch_panel, 'Position', [5, stretch_panel_height-55, figure_width-20, 3]);

            normalized_time = this.data_custodian.getNormalizedTime();
            set(this.time_slider, 'limit', [1, normalized_time(end)], 'value', 1, 'ValueChangingFcn',@(sld,event) this.timeSliderMoving(event));
        end
        function createFigures(this)
            % create containers
            this.number_of_figures = size(this.figure_settings, 1);
            this.figures_list = cell(this.number_of_figures, 1);
            this.axes_list = cell(this.number_of_figures, 1);
            this.line_handle_list = cell(this.number_of_figures, 1);
            
            % define canvas
            screen_size = get(0,'ScreenSize');
            this.canvas_on_screen_width = screen_size(3);
            this.canvas_on_screen_height = screen_size(4) - 24; % 24 is the size of the menu bar on the Mac, check for windows systems later
            
            % create figures
            for i_figure = 1 : size(this.figure_settings, 1)
                % create a figure with these settings
                figure_settings_here = this.figure_settings{i_figure};

                % determine position
                origin_x = figure_settings_here.origin_horizontal * this.canvas_on_screen_width;
                origin_y = figure_settings_here.origin_vertical * this.canvas_on_screen_height;
                width = figure_settings_here.width * this.canvas_on_screen_width;
                height = figure_settings_here.height * this.canvas_on_screen_height - 72;
                figure_position = [origin_x origin_y width height];

                % create figures and axes
                new_figure = figure('Position', figure_position);
                new_axes = axes();
                hold on;
                this.figures_list{i_figure} = new_figure;
                this.axes_list{i_figure} = new_axes;
                
                % set axis scale
                if isfield(figure_settings_here, 'x_lim_lower')
                    set(new_axes, 'XLimMode', 'manual');
                    new_axes.XLim(1) = figure_settings_here.x_lim_lower;
                end
                if isfield(figure_settings_here, 'x_lim_upper')
                    set(new_axes, 'XLimMode', 'manual');
                    new_axes.XLim(2) = figure_settings_here.x_lim_upper;
                end
                if isfield(figure_settings_here, 'y_lim_lower')
                    set(new_axes, 'YLimMode', 'manual');
                    new_axes.YLim(1) = figure_settings_here.y_lim_lower;
                end
                if isfield(figure_settings_here, 'y_lim_upper')
                    set(new_axes, 'YLimMode', 'manual');
                    new_axes.YLim(2) = figure_settings_here.y_lim_upper;
                end

                % create empty plots
                data_to_plot = figure_settings_here.data;
                data_header = figure_settings_here.data_header;
                line_handles = cell(size(data_to_plot, 1), 1);
                point_handles = cell(size(data_to_plot, 1), 1);
                for i_plot = 1 : size(data_to_plot, 1)
                    color_R = str2double(data_to_plot{i_plot, strcmp(data_header, 'color R')});
                    color_G = str2double(data_to_plot{i_plot, strcmp(data_header, 'color G')});
                    color_B = str2double(data_to_plot{i_plot, strcmp(data_header, 'color B')});
                    line_handles{i_plot} = plot(0, 0, 'color', [color_R, color_G, color_B]);
                    point_handles{i_plot} = plot(0, 0, 'o', 'color', [color_R, color_G, color_B]);
                end
                this.line_handle_list{i_figure} = line_handles;
                this.point_handle_list{i_figure} = point_handles;
            end
            
        end
        function setCurrentTime(this, current_time)
            this.current_time = current_time;
            this.updateCurrentTimePlots();
        end
        function setStretch(this, stretch_index)
            this.current_stretch = stretch_index;
            for i_figure = 1 : this.number_of_figures
                figure_settings_here = this.figure_settings{i_figure};
                line_handles_here = this.line_handle_list{i_figure};
                for i_variable = 1 : size(figure_settings_here.data, 1)
                    % get data
                    x_data_name = figure_settings_here.data{i_variable, strcmp(figure_settings_here.data_header, 'x-data name')};
                    x_data_type = figure_settings_here.data{i_variable, strcmp(figure_settings_here.data_header, 'x-data type')};
                    y_data_name = figure_settings_here.data{i_variable, strcmp(figure_settings_here.data_header, 'y-data name')};
                    y_data_type = figure_settings_here.data{i_variable, strcmp(figure_settings_here.data_header, 'y-data type')};
                    x_data = this.data_custodian.getData(x_data_name, x_data_type);
                    y_data = this.data_custodian.getData(y_data_name, y_data_type);
                    x_data_this_stretch = x_data.variable_data(:, stretch_index);
                    y_data_this_stretch = y_data.variable_data(:, stretch_index);
                    
                    % update plot
                    this_line = line_handles_here{i_variable};
                    set(this_line, 'xdata', x_data_this_stretch, 'ydata', y_data_this_stretch);
                end
            end
            
            set(this.current_stretch_edit, 'Value', num2str(stretch_index));
            this.updateCurrentTimePlots();
        end
        function updateCurrentTimePlots(this)
            
            for i_figure = 1 : this.number_of_figures
                figure_settings_here = this.figure_settings{i_figure};
                point_handles_here = this.point_handle_list{i_figure};
                for i_variable = 1 : size(figure_settings_here.data, 1)
                    % get data
                    x_data_name = figure_settings_here.data{i_variable, strcmp(figure_settings_here.data_header, 'x-data name')};
                    x_data_type = figure_settings_here.data{i_variable, strcmp(figure_settings_here.data_header, 'x-data type')};
                    y_data_name = figure_settings_here.data{i_variable, strcmp(figure_settings_here.data_header, 'y-data name')};
                    y_data_type = figure_settings_here.data{i_variable, strcmp(figure_settings_here.data_header, 'y-data type')};
                    x_data = this.data_custodian.getData(x_data_name, x_data_type);
                    y_data = this.data_custodian.getData(y_data_name, y_data_type);
                    x_data_this_stretch = x_data.variable_data(:, this.current_stretch);
                    y_data_this_stretch = y_data.variable_data(:, this.current_stretch);
                    
                    % update plot
                    this_point = point_handles_here{i_variable};
                    set(this_point, 'xdata', x_data_this_stretch(this.current_time), 'ydata', y_data_this_stretch(this.current_time));
                end
            end
            
        end
        
        % callbacks
        function processKeyPress(this, ~, eventdata)
            if strcmp(eventdata.Key, 'a') && isempty(eventdata.Modifier)
                this.setStretch(this.current_stretch - 1);
            end
            if strcmp(eventdata.Key, 'd') && isempty(eventdata.Modifier)
                this.setStretch(this.current_stretch + 1);
            end
        end
        function timeSliderMoving(this, sender, ~)
            this.setCurrentTime(round(sender.Value));
        end
    end
    
end








