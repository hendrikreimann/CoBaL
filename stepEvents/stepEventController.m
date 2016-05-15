classdef stepEventController < handle
    properties
        trial_data;
        event_data;
        
        selected_event_label;
        selected_event_time;
        
        event_time_normal_step = 0.005;
        event_time_large_step = 0.050;
        
        control_figure;
        saveFigureSettingsButton;
        loadFigureSettingsButton;
        findEventsButton;
        saveEventsButton;
        figureSelectionBox;
        heel_pos_peak_width;
        toes_vel_peak_width;
    end
    methods
        function controller = stepEventController(trial_data, event_data)
            controller.trial_data = trial_data;
            controller.event_data = event_data;

            figure_height = 600;
            figure_width = 420;
            controller.control_figure = figure('position', [1600 300 figure_width figure_height], 'Units', 'pixels', 'KeyPressFcn', @controller.processKeyPress);

            % figure control
            figures_panel_height = 100;
            figures_panel = uipanel(controller.control_figure, 'Title', 'Figure Control', 'FontSize', 12, 'BackgroundColor', 'white', 'Units', 'pixels', 'Position', [5, figure_height-figures_panel_height-5, figure_width-10, figures_panel_height]);
            controller.figureSelectionBox = uicontrol(figures_panel, 'Style', 'popup', 'String', '<no figure>', 'Position', [5, figures_panel_height-40, 395, 20], 'Fontsize', 12, 'HorizontalAlignment', 'left');

            controller.saveFigureSettingsButton = uicontrol(figures_panel, 'Style', 'pushbutton', 'Position', [5, figures_panel_height-100, 130, 60], 'Fontsize', 12, 'String', 'Save Figure Settings', 'callback', @controller.saveFigureSettings, 'HitTest', 'off');
            controller.loadFigureSettingsButton = uicontrol(figures_panel, 'Style', 'pushbutton', 'Position', [140, figures_panel_height-100, 130, 60], 'Fontsize', 12, 'String', 'Load Figure Settings', 'callback', @controller.loadFigureSettings, 'HitTest', 'off');

            % event controls
            events_panel_height = 255;
            events_panel = uipanel(controller.control_figure, 'Title', 'Events Control', 'FontSize', 12, 'BackgroundColor', 'white', 'Units', 'pixels', 'Position', [5, figure_height-figures_panel_height-events_panel_height-5, figure_width-10, events_panel_height]);
            controller.findEventsButton = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [5, events_panel_height-75, 130, 60], 'Fontsize', 12, 'String', 'Find Events', 'callback', @controller.findEvents, 'HitTest', 'off');
            controller.saveEventsButton = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [140, events_panel_height-75, 130, 60], 'Fontsize', 12, 'String', 'Save Events', 'callback', @event_data.saveEvents, 'HitTest', 'off');

            uicontrol(events_panel, 'Style', 'text', 'string', 'Heel pos peak prominence (m):', 'Position', [5, events_panel_height-100, 190, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
            controller.heel_pos_peak_width = uicontrol(events_panel, 'Style', 'edit', 'BackgroundColor', 'white', 'Position', [150, events_panel_height-97, 40, 20], 'String', '0.05');
            uicontrol(events_panel, 'Style', 'text', 'string', 'Toes vel peak prominence (m):', 'Position', [5, events_panel_height-120, 190, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
            controller.toes_vel_peak_width = uicontrol(events_panel, 'Style', 'edit', 'BackgroundColor', 'white', 'Position', [150, events_panel_height-117, 40, 20], 'String', '0.05');
        end
        
        
        function setSelectedEvent(this, event_label, event_time)
            if nargin == 1
                event_label = 'left_touchdown';
                event_times = this.event_data.getEventTimes(event_label);
                event_time = event_times(1);
            end
            this.selected_event_label = event_label;
            this.selected_event_time = event_time;
            
            this.updateSelectedEventPlots();
        end
        function updateSelectedEventPlots(this)
            for i_figure = 1 : length(this.figureSelectionBox.String)
                this.figureSelectionBox.UserData{i_figure}.updateSelectedEventPlot();
            end
        end
        function updateEventPlots(this)
            for i_figure = 1 : length(this.figureSelectionBox.String)
                this.figureSelectionBox.UserData{i_figure}.updateEventPlots();
            end
        end
    
        
        function processKeyPress(this, sender, eventdata)
            if strcmp(eventdata.Key, 'a') && isempty(eventdata.Modifier)
                this.selectPreviousEvent();
                this.updateSelectedEventPlots();
            elseif strcmp(eventdata.Key, 'd') && isempty(eventdata.Modifier)
                this.selectNextEvent();
                this.updateSelectedEventPlots();
            elseif strcmp(eventdata.Key, 'w') && isempty(eventdata.Modifier)
                this.updateTimeWindow('zoom in');
            elseif strcmp(eventdata.Key, 'x') && isempty(eventdata.Modifier)
                this.updateTimeWindow('zoom out');
            elseif strcmp(eventdata.Key, 's') && isempty(eventdata.Modifier)
                this.updateTimeWindow('center');
            elseif strcmp(eventdata.Key, 'q') && isempty(eventdata.Modifier)
                this.updateTimeWindow('previous');
            elseif strcmp(eventdata.Key, 'e') && isempty(eventdata.Modifier)
                this.updateTimeWindow('next');
            elseif strcmp(eventdata.Key, 'z') || strcmp(eventdata.Key, 'c')
                this.moveSelectedEvent(sender, eventdata);
            elseif strcmp(eventdata.Key, 'x') && strcmp(eventdata.Modifier, 'command')
                % this won't work on windows systems, adapt that!
                this.quit();
            end
        end
        function saveFigureSettings(this, sender, eventdata)
            figure_settings = cell(1, length(this.figureSelectionBox));
            for i_figure = 1 : length(this.figureSelectionBox.String)
                figure_settings{i_figure} = this.figureSelectionBox.UserData{i_figure}.getSetting();
            end
            control_figure_setting = struct();
            control_figure_setting.position = this.control_figure.Position;

            save_file = '/Users/reimajbi/Library/Application Support/stepEventFinder/figureSettings.mat';
            save(save_file, 'figure_settings', 'control_figure_setting');
        end
        function loadFigureSettings(this, sender, eventdata)
            load_file = '/Users/reimajbi/Library/Application Support/stepEventFinder/figureSettings.mat';

            load(load_file, 'figure_settings', 'control_figure_setting')
            for i_figure = 1 : length(figure_settings)
    %             if length(step_event_figures) < i_figure
    %                 step_event_figures{i_figure} = createStepEventFigure();
    %             end
                this.figureSelectionBox.UserData{i_figure}.applySettings(figure_settings{i_figure});
            end
            this.control_figure.Position = control_figure_setting.position;
        end
        function findEvents(this, sender, eventdata)
            % this should eventually be moved to WalkingEventData
            
            % left touchdown
            [~, left_touchdown_indices_mocap] = findpeaks(-this.trial_data.getData('left_heel_z_pos'), 'MinPeakProminence', str2num(this.heel_pos_peak_width.String));
            time = this.trial_data.getTime('left_heel_z_pos');
            left_touchdown_times = time(left_touchdown_indices_mocap);
            this.event_data.setEventTimes(left_touchdown_times, 'left_touchdown');

            % left pushoff
            [~, left_pushoff_indices_mocap] = findpeaks(this.trial_data.getData('left_toes_z_vel'), 'MinPeakProminence', str2num(this.toes_vel_peak_width.String));
            time = this.trial_data.getTime('left_toes_z_vel');
            left_pushoff_times = time(left_pushoff_indices_mocap);
            this.event_data.setEventTimes(left_pushoff_times, 'left_pushoff');

            % right touchdown
            [~, right_touchdown_indices_mocap] = findpeaks(-this.trial_data.getData('right_heel_z_pos'), 'MinPeakProminence', str2num(this.heel_pos_peak_width.String));
            time = this.trial_data.getTime('right_heel_z_pos');
            right_touchdown_times = time(right_touchdown_indices_mocap);
            this.event_data.setEventTimes(right_touchdown_times, 'right_touchdown');

            % right pushoff
            [~, right_pushoff_indices_mocap] = findpeaks(this.trial_data.getData('right_toes_z_vel'), 'MinPeakProminence', str2num(this.toes_vel_peak_width.String));
            time = this.trial_data.getTime('right_toes_z_vel');
            right_pushoff_times = time(right_pushoff_indices_mocap);
            this.event_data.setEventTimes(right_pushoff_times, 'right_pushoff');

            this.updateEventPlots();


        end
        
        function selectNextEvent(this)
            % find index of currently selected event
            currently_selected_event_time = this.selected_event_time;
            event_data_of_current_type = this.event_data.getEventTimes(this.selected_event_label);
            event_data_of_current_type_after_currently_selected = event_data_of_current_type(event_data_of_current_type > currently_selected_event_time);
            
            if isempty(event_data_of_current_type_after_currently_selected)
                % the selected event was the last of this type, so go to next type
                this.selected_event_label = this.event_data.getNextEventTypeLabel(this.selected_event_label);
                event_data_of_current_type = this.event_data.getEventTimes(this.selected_event_label);
                this.selected_event_time = event_data_of_current_type(1);
            else
                % select next one
                this.selected_event_time = event_data_of_current_type_after_currently_selected(1);
            end
        end
        function selectPreviousEvent(this)
            % find index of currently selected event
            currently_selected_event_time = this.selected_event_time;
            event_data_of_current_type = this.event_data.getEventTimes(this.selected_event_label);
            event_data_of_current_type_before_currently_selected = event_data_of_current_type(event_data_of_current_type < currently_selected_event_time);
            
            if isempty(event_data_of_current_type_before_currently_selected)
                % the selected event was the last of this type, so go to next type
                this.selected_event_label = this.event_data.getPreviousEventTypeLabel(this.selected_event_label);
                event_data_of_current_type = this.event_data.getEventTimes(this.selected_event_label);
                this.selected_event_time = event_data_of_current_type(end);
            else
                % select next one
                this.selected_event_time = event_data_of_current_type_before_currently_selected(end);
            end
        end
        function moveSelectedEvent(this, sender, eventdata)
            selected_event_time_current = this.selected_event_time;
            
            % step forward or backward
            if strcmp(eventdata.Key, 'z') && isempty(eventdata.Modifier)
                selected_event_time_new = selected_event_time_current - this.event_time_normal_step;
            elseif strcmp(eventdata.Key, 'z') && strcmp(eventdata.Modifier, 'shift')
                selected_event_time_new = selected_event_time_current - this.event_time_large_step;
            elseif strcmp(eventdata.Key, 'c') && isempty(eventdata.Modifier)
                selected_event_time_new = selected_event_time_current + this.event_time_normal_step;
            elseif strcmp(eventdata.Key, 'c') && strcmp(eventdata.Modifier, 'shift')
                selected_event_time_new = selected_event_time_current + this.event_time_large_step;
            else
                % modifier not recognized, don't do anything
                return;
            end
            
            % apply limits
            this.selected_event_time = this.event_data.updateEventTime(this.selected_event_label, selected_event_time_current, selected_event_time_new);
            
            % update plots
            this.updateSelectedEventPlots();
            this.updateEventPlots();
        end
        function updateTimeWindow(this, mode)
            % get current time extension and calculate new value
            epsilon = 0.01; % to avoid numerical inconsistencies, should be around the largest delta t in the data
            current_range = this.figureSelectionBox.UserData{1}.main_axes.XLim;
            current_extension = current_range(2) - current_range(1);
            current_center = mean(current_range);
            if strcmp(mode, 'zoom in')
                new_extension = max(this.figureSelectionBox.UserData{1}.time_extension_steps(this.figureSelectionBox.UserData{1}.time_extension_steps < current_extension - epsilon));
                new_center = current_center;
            elseif strcmp(mode, 'zoom out')
                new_extension = min(this.figureSelectionBox.UserData{1}.time_extension_steps(this.figureSelectionBox.UserData{1}.time_extension_steps > current_extension + epsilon));
                new_center = current_center;
            elseif strcmp(mode, 'center')
                new_extension = current_extension;
                new_center = this.selected_event_time;
            elseif strcmp(mode, 'previous')
                new_extension = current_extension;
                new_center = current_center - current_extension;
            elseif strcmp(mode, 'next')
                new_extension = current_extension;
                new_center = current_center + current_extension;
            end
            if isempty(new_extension)
                return
            end
            
            % calculate new time range
            new_x_lim = new_center + [-0.5 0.5]*new_extension;
            if new_x_lim(1) < 0
                new_x_lim = new_x_lim - new_x_lim(1);
            end
            if new_x_lim(2) > this.trial_data.recording_time
                new_x_lim = new_x_lim - (new_x_lim(2) - this.trial_data.recording_time)
            end
            set(this.figureSelectionBox.UserData{1}.main_axes, 'xlim', new_x_lim);
            
        end
        
        function quit(this)
            for i_figure = 1 : length(this.figureSelectionBox.String)
                try
                    close(this.figureSelectionBox.UserData{i_figure}.main_figure);
                catch exception
                    if ~strcmp(exception.identifier, 'MATLAB:close:InvalidFigureHandle')
                        rethrow(exception)
                    end
                end
            end
            close(this.control_figure);
        end
    end
    
end