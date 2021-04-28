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

classdef eventFigure < handle
    properties
        title;
        
        main_figure;
        main_axes;
        
        controller;
        data_custodian;
        event_data;
        
        data_plots;
        event_plots;
        stretch_patches;
        selected_event_plot;
        selected_time_plot;
        ignore_marker_plot;
        problem_marker_plots;
        data_plot_offsets;
        data_plot_scale_factors;
        has_legend_entries = false;

        time_extension_steps = [0.1 0.2 0.5 1 2 5 10 20 40 60 120];
        
        patch_color = [0 0 0];
        patch_alpha = 0.05;
        problem_color = [0.5 0 1];
        problem_alpha = 0.5;
    end
    methods
        function this = eventFigure(figure_settings, controller, data_custodian, eventData)
            % determine position
            canvas_width = controller.canvas_on_screen_width;
            canvas_height = controller.canvas_on_screen_height;
            origin_x = figure_settings.origin_horizontal * canvas_width;
            origin_y = figure_settings.origin_vertical * canvas_height;
            width = figure_settings.width * canvas_width;
            height = figure_settings.height * canvas_height - 72;
            
            figure_position = [origin_x origin_y width height];
            
            % create figure, axes and helper plots
            this.main_figure = figure('Position', figure_position, 'KeyPressFcn', @controller.processKeyPress);
            this.main_axes = axes('ButtonDownFcn', @this.stepEventFigureClicked);
            this.title = figure_settings.title;
            title(figure_settings.title);
            hold on;
            this.selected_event_plot = plot(0, 0, 'o', 'markersize', 15, 'linewidth', 3, 'color', [1 0.5 0], 'visible', 'off', 'HandleVisibility', 'off');
%             this.selected_time_plot = plot(0, 0, '+', 'markersize', 25, 'linewidth', 1, 'color', [1 1 1]*0.7, 'ButtonDownFcn', @this.stepEventFigureClicked, 'HandleVisibility', 'off');
            this.selected_time_plot = plot([0 0], [0 0], 'linewidth', 1, 'color', [1 1 1]*0.7, 'ButtonDownFcn', @this.stepEventFigureClicked, 'HandleVisibility', 'off', 'visible', 'off');
            this.ignore_marker_plot = plot(0, 0, 'x', 'markersize', 25, 'linewidth', 2, 'color', [1 0 0.5]*0.7, 'ButtonDownFcn', @this.stepEventFigureClicked, 'HandleVisibility', 'off', 'visible', 'off');
            this.problem_marker_plots = {};

            % register with controller
            if strcmp(controller.figureSelectionBox.String, '<no figure>')
                controller.figureSelectionBox.String = figure_settings.title;
                controller.figureSelectionBox.UserData = {this};
            else
                controller.figureSelectionBox.String = [controller.figureSelectionBox.String; {figure_settings.title}];
                controller.figureSelectionBox.UserData = [controller.figureSelectionBox.UserData; {this}];
            end
            
            % store handles
            this.controller = controller;
            this.data_custodian = data_custodian;
            this.event_data = eventData;
            
        end
        function addDataPlot(this, data_label, legend_label, color, scale_factor, offset, linewidth)
            % define defaults if not all properties are specified
            if nargin < 3
                color = rand(1, 3);
            end
            if nargin < 4
                scale_factor = 1;
            end
            if nargin < 5
                offset = 0;
            end
            
            % create the plot
            time = this.data_custodian.getTimeData(data_label);
            data = this.data_custodian.getBasicVariableData(data_label);
            new_plot = plot ...
              ( ...
                time, ...
                (data + offset) * scale_factor, ...
                'color', color, ...
                'linewidth', linewidth, ...
                'ButtonDownFcn', @this.stepEventFigureClicked ...
              );
          
            % set properties related to legend
            if strcmp(legend_label, '~')
                set(new_plot, 'HandleVisibility', 'off')
            else
                set(new_plot, 'DisplayName', strrep(legend_label, '_', ' '))
                this.has_legend_entries = true;
            end
            new_plot.UserData = data_label;
            this.data_plots{length(this.data_plots)+1} = new_plot;
            
            % store offset and scale factors
            this.data_plot_offsets(length(this.data_plot_offsets)+1) = offset;
            this.data_plot_scale_factors(length(this.data_plot_scale_factors)+1) = scale_factor;
            
        end
        function addEventPlot(this, data_label, event_label, legend_label, color, marker, marker_size)
            % define defaults if not all properties are specified
            if nargin < 3
                color = rand(1, 3);
            end
            if nargin < 4
                marker = 'o';
            end
            if nargin < 5
                marker_size = 6;
            end
            data_plot_handle = this.getDataPlot(data_label);
            
            % create the plot with placeholder data
            new_plot = plot ...
              ( ...
                this.main_axes, ...
                0, 0, 'o', ...
                'linewidth', 2, ...
                'color', color, ...
                'marker', marker, ...
                'markersize', marker_size, ...
                'ButtonDownFcn', @this.stepEventFigureClicked ...
              );
          
            % set properties related to legend
            if strcmp(legend_label, '~')
                set(new_plot, 'HandleVisibility', 'off')
            else
                set(new_plot, 'DisplayName', strrep(legend_label, '_', ' '))
                this.has_legend_entries = true;
            end
            new_plot.UserData = {data_plot_handle, event_label};
            this.event_plots{length(this.event_plots)+1} = new_plot;
            
            % update to show actual data
            this.updateEventPlots();
        end
        function data_plot = getDataPlot(this, data_label)
            i_plot = 1;
            while i_plot <= length(this.data_plots) && ~strcmp(this.data_plots{i_plot}.UserData, data_label)
                i_plot = i_plot+1;
            end
            if i_plot > length(this.data_plots)
                i_plot = 1;
                warning(['Label ' data_label ' not found, defaulting to ' this.data_plots{i_plot}.UserData])
            end
            data_plot = this.data_plots{i_plot};
        end
        
        function stepEventFigureClicked(this, sender, eventdata) %#ok<INUSD>
            % get click coordinates in pixels
            current_point_pixel = get(this.main_figure, 'CurrentPoint');
            current_point_pixel = current_point_pixel(1,1:2);
            
            if strcmp(this.controller.getEventMode, 'select')
                [event_label, event_time] = this.getClosestEvent(current_point_pixel);
                this.controller.setSelectedEvent(event_label, event_time);
            elseif strcmp(this.controller.getEventMode, 'add')
                point_time = this.getClosestDataPoint(current_point_pixel);
                this.controller.addEvent(point_time);
            end
                
        end
        function [event_label, event_time, distance] = getClosestEvent(this, point_pixel)
            % determine closest event for each type
            number_of_event_types = length(this.event_plots);
            number_of_problems = length(this.problem_marker_plots);
            candidate_distances = zeros(1, number_of_event_types + number_of_problems);
            candidate_event_indices = zeros(1, number_of_event_types + number_of_problems);
            for i_type = 1 : number_of_event_types
                [candidate_distances(i_type), candidate_event_indices(i_type)] = this.calculatePointToCurvePixelDistance(this.event_plots{i_type}, point_pixel);
            end
%             [candidate_distances(length(this.event_plots)+1), candidate_event_indices(length(this.event_plots)+1)] = this.calculatePointToCurvePixelDistance(this.ignore_marker_plot, point_pixel);
            for i_problem = 1 : number_of_problems
                this_problem_marker = this.problem_marker_plots{i_problem};
                [candidate_distances(number_of_event_types + i_problem), candidate_event_indices(number_of_event_types + i_problem)] ...
                    = this.calculatePointToCurvePixelDistance(this_problem_marker, point_pixel);
                
            end
            
            % find the one with minimal distance among these candidates
            [distance, type_index] = min(candidate_distances);
%             if type_index == length(this.event_plots)+1
%                 event_label = 'ignore_times';
%                 event_index = candidate_event_indices(length(this.event_plots)+1);
%                 event_times = this.event_data.getEventTimes(event_label);
%                 event_time = event_times(event_index);
%             else
%                 event_label = this.event_plots{type_index}.UserData{2};
%                 event_index = candidate_event_indices(type_index);
%                 event_times = this.event_data.getEventTimes(event_label);
%                 event_time = event_times(event_index);
%             end            
            if type_index <= number_of_event_types
                event_label = this.event_plots{type_index}.UserData{2};
                event_index = candidate_event_indices(type_index);
                event_times = this.event_data.getEventTimes(event_label);
                event_time = event_times(event_index);
            else
                event_label = 'problem';
                problem_index = type_index - number_of_event_types;
                event_time = this.event_data.problem_table.start_time(problem_index);
            end            
        end
        function point_time = getClosestDataPoint(this, point_pixel)
            candidate_distances = zeros(1, length(this.data_plots));
            candidate_time_indices = zeros(1, length(this.data_plots));
            for i_plot = 1 : length(this.data_plots)
                [candidate_distances(i_plot), candidate_time_indices(i_plot)] = this.calculatePointToCurvePixelDistance(this.data_plots{i_plot}, point_pixel);
            end
            
            % find the one with minimal distance among these candidates
            [~, plot_index] = min(candidate_distances);
            point_time = this.data_plots{plot_index}.XData(candidate_time_indices(plot_index));
        end
        function [distance, point_index] = calculatePointToCurvePixelDistance(this, curve, point_pixel)
            % transform curve to pixel coordinates
            axes_origin_pixel = getpixelposition(this.main_axes);
            xlim = get(this.main_axes, 'xlim');
            ylim = get(this.main_axes, 'ylim');
            
            curve_x_pixel = axes_origin_pixel(1) + axes_origin_pixel(3) * (curve.XData-xlim(1))/(xlim(2)-xlim(1));
            curve_y_pixel = axes_origin_pixel(2) + axes_origin_pixel(4) * (curve.YData-ylim(1))/(ylim(2)-ylim(1));
            curve_data_pixel = [curve_x_pixel; curve_y_pixel];
            
            
            
            % find closest point on curve
            difference_vectors = curve_data_pixel - repmat(point_pixel', 1, size(curve_data_pixel, 2));
            distance_squared = sum((difference_vectors.^2), 1);
            [distance, point_index] = min(distance_squared);
            
            if isempty(distance)
                distance = inf;
                point_index = -1;
            end
        end
        function toggleLegend(this)
            if this.has_legend_entries
                legend(this.main_axes, 'toggle')
            end
        end
        
        function setting_struct = getSetting(this)
            setting_struct = struct();
            setting_struct.title = this.main_axes.Title.String;
            setting_struct.position = this.main_figure.Position;
            
            setting_struct.data_labels = cell(size(this.data_plots));
            for i_data_plot = 1 : length(this.data_plots)
                setting_struct.data_labels{i_data_plot} = this.data_plots{i_data_plot}.UserData;
            end
            
        end
        function applySettings(this, setting_struct)
            this.main_figure.Position = setting_struct.position;
%             this.main_axes.Title.String = setting_struct.title;
            
        end
        function updateEventPlots(this)
            % loop through all event plots
            for i_plot = 1 : length(this.event_plots)
                % get handle of parent data plot
                data_plot_handle = this.event_plots{i_plot}.UserData{1};
                
                % get event label
                event_label = this.event_plots{i_plot}.UserData{2};
                
                % get data
                time = get(data_plot_handle, 'xdata');
                data = get(data_plot_handle, 'ydata');
                event_time = this.event_data.getEventTimes(event_label);
                event_data = interp1(time, data, event_time); %#ok<PROP>
                
                % update
                set(this.event_plots{i_plot}, 'xdata', event_time, 'ydata', event_data); %#ok<PROP>
            end
            this.updateIgnoreMarkerPlot();
            this.updateProblemMarkers();
        end
        function updateDataPlots(this)
            % loop through all data plots
            for i_plot = 1 : length(this.data_plots)
                % get plot handle, modifiers and label
                data_plot_handle = this.data_plots{i_plot};
                data_label = this.data_plots{i_plot}.UserData;
                offset = this.data_plot_offsets(i_plot);
                scale_factor = this.data_plot_scale_factors(i_plot);
                
                % get data
                time = this.data_custodian.getTimeData(data_label);
                data = this.data_custodian.getBasicVariableData(data_label);
                
                % update
                set(data_plot_handle, 'xdata', time, 'ydata', (data + offset) * scale_factor);
            end
        end
        function updateStretchPatches(this)
            % delete old stretch plots
            delete(this.stretch_patches);
            this.stretch_patches = [];
            
            ylimits = get(this.main_axes, 'ylim');
            
            % determine blocks of stretches
            stretch_start_times = this.event_data.stretch_start_times;
            stretch_end_times = this.event_data.stretch_end_times;
            [block_start_times, sort_indices] = sort(stretch_start_times);
            block_end_times = stretch_end_times(sort_indices);
            
            index = 2;
            while index <= length(block_end_times)
                while index <= length(block_end_times) && (block_start_times(index) == block_end_times(index - 1))
                    block_start_times(index) = [];
                    block_end_times(index-1) = [];
                end
                index = index + 1;
            end            
            
            
            for i_block = 1 : length(block_start_times)
                stretch_start = block_start_times(i_block);
                stretch_end = block_end_times(i_block);
                
                patch_x = [stretch_start stretch_end stretch_end stretch_start];
                patch_y = [ylimits(1) ylimits(1) ylimits(2) ylimits(2)];
                patch_handle = ...
                    patch ...
                      ( ...
                        patch_x, ...
                        patch_y, ...
                        this.patch_color, ...
                        'parent', this.main_axes, ...
                        'EdgeColor', lightenColor(this.patch_color, 0.8), ...
                        'FaceAlpha', this.patch_alpha, ...
                        'HandleVisibility', 'off', ...
                        'ButtonDownFcn', @this.stepEventFigureClicked ...
                      ); 
                uistack(patch_handle, 'bottom')
                this.stretch_patches = [this.stretch_patches, patch_handle];
            end
        end
        function updateIgnoreMarkerPlot(this)
            ignore_times = this.event_data.ignore_times;
            set(this.ignore_marker_plot, 'xdata', ignore_times, 'ydata', zeros(size(ignore_times)));
        end
        function updateProblemMarkers(this)
            % clear out old problem markers
            for i_plot = 1 : numel(this.problem_marker_plots)
                delete(this.problem_marker_plots{i_plot});
            end
            this.problem_marker_plots = {};
            
            % create new plots
            problems = this.event_data.problem_table;
            problems_added_to_legend = false;
            ylimits = get(this.main_axes, 'ylim');
            for i_problem = 1 : size(problems, 1)
                this_start_time = problems.start_time(i_problem);
                this_end_time = problems.end_time(i_problem);
                this_reason = problems.reason(i_problem);
                
                if this_start_time == this_end_time
                    % draw single marker
                    this.problem_marker_plots{i_problem} = plot ...
                      ( ...
                        this_start_time, ...
                        mean(ylimits), ...
                        'x', ...
                        'color', this.problem_color, ...
                        'parent', this.main_axes, ...
                        'linewidth', 2, ...
                        'markersize', 12, ...
                        'ButtonDownFcn', @this.stepEventFigureClicked ...
                      );
                    if ~problems_added_to_legend
                        set(this.problem_marker_plots{i_problem}, 'DisplayName', 'problems');
                        problems_added_to_legend = true;
                    else
                        set(this.problem_marker_plots{i_problem}, 'HandleVisibility', 'off');
                    end


                else
                    
                    this.problem_marker_plots{i_problem} = plot ...
                      ( ...
                        [this_start_time, this_end_time], ...
                        mean(ylimits) * [1 1], ...
                        '-', ...
                        'color', this.problem_color, ...
                        'parent', this.main_axes, ...
                        'linewidth', 8, ...
                        'ButtonDownFcn', @this.stepEventFigureClicked ...
                      );
                    if ~problems_added_to_legend
                        set(this.problem_marker_plots{i_problem}, 'DisplayName', 'problems');
                        problems_added_to_legend = true;
                    else
                        set(this.problem_marker_plots{i_problem}, 'HandleVisibility', 'off');
                    end
                    
                end
            end
            
            
        end
        function updateSelectedEventPlot(this)
            % go through all event plots and check whether they are of the selected event
            selected_event_plot_x_data = [];
            selected_event_plot_y_data = [];
            for i_plot = 1 : length(this.event_plots)
                if strcmp(this.controller.event_data.selected_event_label, this.event_plots{i_plot}.UserData{2})
                    % event type of this one and the selected is a match
                    selected_event_plot_x_data = [selected_event_plot_x_data this.controller.event_data.selected_event_time]; %#ok<AGROW>
                    y_data_point = interp1(this.event_plots{i_plot}.UserData{1}.XData, this.event_plots{i_plot}.UserData{1}.YData, this.controller.event_data.selected_event_time);
                    selected_event_plot_y_data = [selected_event_plot_y_data y_data_point]; %#ok<AGROW>
                end
            end
            if strcmp(this.controller.event_data.selected_event_label, 'ignore_times')
                % event type of this one and the selected is a match
                selected_event_plot_x_data = [selected_event_plot_x_data this.controller.event_data.selected_event_time];
                y_data_point = 0;
                selected_event_plot_y_data = [selected_event_plot_y_data y_data_point];
            end
            if strcmp(this.controller.event_data.selected_event_label, 'problem')
                % event type of this one and the selected is a match
                selected_event_plot_x_data = this.controller.event_data.selected_event_time;
                selected_event_plot_y_data = mean(this.main_axes.YLim);
            end
            set(this.selected_event_plot, 'xdata', selected_event_plot_x_data, 'ydata', selected_event_plot_y_data);
            if isempty(selected_event_plot_x_data)
                set(this.selected_event_plot, 'visible', 'off');
            else
                set(this.selected_event_plot, 'visible', 'on');
            end
        end
        function updateSelectedTimePlot(this)
            if isempty(this.event_data.selected_time)
                set(this.selected_time_plot, 'visible', 'off');
            else
                set(this.selected_time_plot, 'xdata', this.event_data.selected_time * [1 1]);
                
                % scale properly by hiding plot, then adjust y-data of plot and
                % showing again
                set(this.selected_time_plot, 'visible', 'off');
                set(this.main_axes, 'YLimMode',  'auto');
                ylimits = this.main_axes.YLim;
                set(this.selected_time_plot, 'ydata', ylimits);
                
                set(this.selected_time_plot, 'visible', 'on');
            end
            
            
        end
        function updateTimeWindow(this, new_time_window)
            if isempty(this.event_data.selected_time)
                % just update xlimits
                set(this.main_axes, 'xlim', new_time_window);
            else
                % hide selected time marker, then update, then update
                % selected time marker
                
                set(this.selected_time_plot, 'visible', 'off');
                set(this.main_axes, 'xlim', new_time_window);
                ylimits = this.main_axes.YLim;
                set(this.selected_time_plot, 'ydata', ylimits);
                set(this.selected_time_plot, 'visible', 'on');
            end
            
        end
    end
    
end