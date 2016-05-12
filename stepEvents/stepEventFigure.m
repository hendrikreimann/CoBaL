classdef stepEventFigure < handle;
    properties
        main_figure;
        main_axes;
        
        trial_data;
        event_data;
        
        data_plots;
        event_plots;

    end
    methods
        function this = stepEventFigure(figureTitle, trialData, eventData)
           
           this.main_figure = figure('ButtonDownFcn', @this.ViewerClickCallback);
           this.main_axes = axes('ButtonDownFcn', @this.ViewerClickCallback);
           title(figureTitle);
           hold on;
           
           this.trial_data = trialData;
           this.event_data = eventData;
           
        end
        function addDataPlot(this, data_label)
            new_plot = plot ...
              ( ...
                this.trial_data.getTime(data_label), ...
                this.trial_data.getData(data_label), ...
                'ButtonDownFcn', @this.ViewerClickCallback ...
              );
            new_plot.UserData = data_label;
            this.data_plots{length(this.data_plots)+1} = new_plot;
        end
        function addEventPlot(this, data_label, event_label)
            data_plot_handle = this.getDataPlot(data_label);
            
            new_plot = plot ...
              ( ...
                this.main_axes, ...
                0, 0, 'o', ...
                'linewidth', 2, ...
                'ButtonDownFcn', @this.ViewerClickCallback ...
              );
            new_plot.UserData = {data_plot_handle, event_label};
            this.event_plots{length(this.event_plots)+1} = new_plot;
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
            % get click coordinates
%             pixel_coordinates = get(viewer_figure, 'CurrentPoint');
%             pixel_coordinates = pixel_coordinates(1,1:2);
% 
%             % take action
%             if strcmp(current_action_state, 'add liftoff event')
%                 % find closest point on any curve to the click
%                 [curve_identifier, index] = pixelCoordinatesToClosestCurvePoint(pixel_coordinates(1, 1:2), sender);
% 
%                 % add the event
%                 addEvent(curve_identifier(1), curve_identifier(2), 1, index);
%                 createContactTrajectories();
%                 updatePlots();
%                 set(viewer_figure, 'Pointer', 'arrow');
%                 set(add_liftoff_event_button, 'BackgroundColor', [.94 .94 .94]);
%                 current_action_state = 'default';
%             elseif strcmp(current_action_state, 'add touchdown event')
%                 % find closest point on any curve to the click
%                 [curve_identifier, index] = pixelCoordinatesToClosestCurvePoint(pixel_coordinates(1, 1:2), sender);
% 
%                 % add the event
%                 addEvent(curve_identifier(1), curve_identifier(2), 2, index);
%                 createContactTrajectories();
%                 updatePlots();
%                 set(viewer_figure, 'Pointer', 'arrow');
%                 set(add_touchdown_event_button, 'BackgroundColor', [.94 .94 .94]);
%                 current_action_state = 'default';
%             elseif strcmp(current_action_state, 'default')
%                 current_event_identifier = pixelCoordinatesToClosestEventIdentifier(pixel_coordinates(1, 1:2), sender);
%                 updateCurrentEventIdentifierPlot();
%             else
%                 error('illegal action state')
%             end
disp('blubclick')
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
            this.main_axes.Title.String = setting_struct.title;
            
        end
        function updateEventPlots(this)
            for i_plot = 1 : length(this.event_plots)
                data_plot_handle = this.event_plots{i_plot}.UserData{1};
                event_label = this.event_plots{i_plot}.UserData{2};
                time = get(data_plot_handle, 'xdata');
                data = get(data_plot_handle, 'ydata');
                event_time = this.event_data.getData(event_label);
                event_data = interp1(time, data, event_time);
                set(this.event_plots{i_plot}, 'xdata', event_time, 'ydata', event_data);
            end
        end
    end
    
end