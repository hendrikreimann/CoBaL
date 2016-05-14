classdef stepEventFigure < handle;
    properties
        main_figure;
        main_axes;
        
        controller;
        trial_data;
        event_data;
        
        data_plots;
        event_plots;
        selected_event_plot;

    end
    methods
        function this = stepEventFigure(figureTitle, controller, trialData, eventData)
           
           this.main_figure = figure('KeyPressFcn', @controller.processKeyPress);
           this.main_axes = axes('ButtonDownFcn', @this.ViewerClickCallback);
           title(figureTitle);
           hold on;
           this.selected_event_plot = plot(0, 0, 'o', 'markersize', 15, 'linewidth', 3, 'color', [1 0.5 0], 'visible', 'off');
           
           this.controller = controller;
           this.trial_data = trialData;
           this.event_data = eventData;
           
        end
        function addDataPlot(this, data_label, color)
            if nargin < 3
                color = rand(1, 3);
            end
            new_plot = plot ...
              ( ...
                this.trial_data.getTime(data_label), ...
                this.trial_data.getData(data_label), ...
                'color', color, ...
                'ButtonDownFcn', @this.ViewerClickCallback ...
              );
            new_plot.UserData = data_label;
            this.data_plots{length(this.data_plots)+1} = new_plot;
        end
        function addEventPlot(this, data_label, event_label, color, marker)
            if nargin < 3
                color = rand(1, 3);
            end
            if nargin < 4
                marker = 'o';
            end
            data_plot_handle = this.getDataPlot(data_label);
            
            new_plot = plot ...
              ( ...
                this.main_axes, ...
                0, 0, 'o', ...
                'linewidth', 2, ...
                'color', color, ...
                'marker', marker, ...
                'ButtonDownFcn', @this.ViewerClickCallback ...
              );
            new_plot.UserData = {data_plot_handle, event_label};
            this.event_plots{length(this.event_plots)+1} = new_plot;
            
            this.updateEventPlots();
        end
        function data_plot = getDataPlot(this, data_label)
            i_plot = 1;
            while i_plot < length(this.data_plots) && ~strcmp(this.data_plots{i_plot}.UserData, data_label)
                i_plot = i_plot+1;
            end
            data_plot = this.data_plots{i_plot};
%             eval(['data = this.' data_label ';']);
        end
        
        function ViewerClickCallback(this, sender, eventdata)
            % get click coordinates in pixels
            current_point_pixel = get(this.main_figure, 'CurrentPoint');
            current_point_pixel = current_point_pixel(1,1:2);
            
            [event_label, event_time, distance] = this.getClosestEvent(current_point_pixel);
            
            
            this.controller.setSelectedEvent(event_label, event_time);
                
        end
        function [event_label, event_time, distance] = getClosestEvent(this, point_pixel)
            
            % determine closest event for each type
            candidate_distances = zeros(1, length(this.event_plots));
            candidate_event_indices = zeros(1, length(this.event_plots));
            for i_type = 1 : length(this.event_plots)
                [candidate_distances(i_type), candidate_event_indices(i_type)] = this.calculatePointToCurvePixelDistance(this.event_plots{i_type}, point_pixel);
            end
            
            % find the one with minimal distance among these candidates
            [distance, type_index] = min(candidate_distances);
            event_label = this.event_plots{type_index}.UserData{2};
            event_index = candidate_event_indices(type_index);
            event_times = this.event_data.getEventTimes(event_label);
            event_time = event_times(event_index);
            
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
            for i_plot = 1 : length(this.event_plots)
                data_plot_handle = this.event_plots{i_plot}.UserData{1};
                event_label = this.event_plots{i_plot}.UserData{2};
                time = get(data_plot_handle, 'xdata');
                data = get(data_plot_handle, 'ydata');
                event_time = this.event_data.getEventTimes(event_label);
                event_data = interp1(time, data, event_time);
                set(this.event_plots{i_plot}, 'xdata', event_time, 'ydata', event_data);
            end
        end
        function updateSelectedEventPlot(this)
            % go through all event plots and check whether they are of the selected event
            selected_event_plot_x_data = [];
            selected_event_plot_y_data = [];
            for i_plot = 1 : length(this.event_plots)
                if strcmp(this.controller.selected_event_label, this.event_plots{i_plot}.UserData{2})
                    % event type of this one and the selected is a match
                    selected_event_plot_x_data = [selected_event_plot_x_data this.controller.selected_event_time];
                    y_data_point = interp1(this.event_plots{i_plot}.UserData{1}.XData, this.event_plots{i_plot}.UserData{1}.YData, this.controller.selected_event_time);
                    selected_event_plot_y_data = [selected_event_plot_y_data y_data_point];
                end
            end
            set(this.selected_event_plot, 'xdata', selected_event_plot_x_data, 'ydata', selected_event_plot_y_data);
            set(this.selected_event_plot, 'visible', 'on');
            
        end
    end
    
end