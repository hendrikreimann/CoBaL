classdef stepEventController < handle
    properties
        trial_data;
        event_data;
        
        event_time_normal_step = 0.005;
        event_time_large_step = 0.050;
        
        control_figure;
        saveFigureSettingsButton;
        loadFigureSettingsButton;
        findEventsButton;
        saveEventsButton;
        addEventButtons;
        
        figureSelectionBox;
        heel_pos_peak_width;
        toes_vel_peak_width;
        
        color_selected = [1 0.5 0];
        color_normal = [0 0 0];
    end
    methods
        function this = stepEventController(trial_data, event_data)
            this.trial_data = trial_data;
            this.event_data = event_data;

            figure_height = 600;
            figure_width = 420;
            this.control_figure = figure('position', [1600 300 figure_width figure_height], 'Units', 'pixels', 'KeyPressFcn', @this.processKeyPress);

            % figure control
            figures_panel_height = 100;
            figures_panel = uipanel(this.control_figure, 'Title', 'Figure Control', 'FontSize', 12, 'BackgroundColor', 'white', 'Units', 'pixels', 'Position', [5, figure_height-figures_panel_height-5, figure_width-10, figures_panel_height]);
            this.figureSelectionBox = uicontrol(figures_panel, 'Style', 'popup', 'String', '<no figure>', 'Position', [5, figures_panel_height-40, 395, 20], 'Fontsize', 12, 'HorizontalAlignment', 'left');

            this.saveFigureSettingsButton = uicontrol(figures_panel, 'Style', 'pushbutton', 'Position', [5, figures_panel_height-100, 130, 60], 'Fontsize', 12, 'String', 'Save Figure Settings', 'callback', @this.saveFigureSettings);
            this.loadFigureSettingsButton = uicontrol(figures_panel, 'Style', 'pushbutton', 'Position', [140, figures_panel_height-100, 130, 60], 'Fontsize', 12, 'String', 'Load Figure Settings', 'callback', @this.loadFigureSettings);

            % event controls
            events_panel_height = 255;
            events_panel = uipanel(this.control_figure, 'Title', 'Events Control', 'FontSize', 12, 'BackgroundColor', 'white', 'Units', 'pixels', 'Position', [5, figure_height-figures_panel_height-events_panel_height-5, figure_width-10, events_panel_height]);
            this.findEventsButton = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [5, events_panel_height-75, 130, 60], 'Fontsize', 12, 'String', 'Find Events', 'callback', @this.findEvents);
            this.saveEventsButton = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [140, events_panel_height-75, 130, 60], 'Fontsize', 12, 'String', 'Save Events', 'callback', @event_data.saveEvents);
            
            this.addEventButtons(1) = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [5, events_panel_height-135, 100, 60], 'Fontsize', 12, 'String', '<html><center>Add Left<br></center>Pushoff', 'callback', @this.addEventButtonPressed, 'UserData', 'left_pushoff');
            this.addEventButtons(2) = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [105, events_panel_height-135, 100, 60], 'Fontsize', 12, 'String', '<html><center>Add Left<br></center>Touchdown', 'callback', @this.addEventButtonPressed, 'UserData', 'left_touchdown');
            this.addEventButtons(3) = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [205, events_panel_height-135, 100, 60], 'Fontsize', 12, 'String', '<html><center>Add Right<br></center>Pushoff', 'callback', @this.addEventButtonPressed, 'UserData', 'right_pushoff');
            this.addEventButtons(4)= uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [305, events_panel_height-135, 100, 60], 'Fontsize', 12, 'String', '<html><center>Add Right<br></center>Touchdown', 'callback', @this.addEventButtonPressed, 'UserData', 'right_touchdown');

            uicontrol(events_panel, 'Style', 'text', 'string', 'Heel pos peak prominence (m):', 'Position', [5, events_panel_height-160, 190, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
            this.heel_pos_peak_width = uicontrol(events_panel, 'Style', 'edit', 'BackgroundColor', 'white', 'Position', [150, events_panel_height-157, 40, 20], 'String', '0.05');
            uicontrol(events_panel, 'Style', 'text', 'string', 'Toes vel peak prominence (m):', 'Position', [5, events_panel_height-180, 190, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
            this.toes_vel_peak_width = uicontrol(events_panel, 'Style', 'edit', 'BackgroundColor', 'white', 'Position', [150, events_panel_height-177, 40, 20], 'String', '0.05');
        end
        
        
        function setSelectedEvent(this, event_label, event_time)
            if nargin == 1
                event_label = 'left_touchdown';
                event_times = this.event_data.getEventTimes(event_label);
                event_time = event_times(1);
            end
            this.event_data.selected_event_label = event_label;
            this.event_data.selected_event_time = event_time;
            
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
        function event_mode = getEventMode(this)
            event_mode = 'select';
            for i_button = 1 : length(this.addEventButtons)
                if get(this.addEventButtons(i_button), 'ForegroundColor') == this.color_selected;
                    event_mode = 'add';
                end
            end

        end
        
        function processKeyPress(this, sender, eventdata)
            if strcmp(eventdata.Key, 'a') && isempty(eventdata.Modifier)
                this.event_data.selectPreviousEvent();
                this.updateSelectedEventPlots();
            elseif strcmp(eventdata.Key, 'd') && isempty(eventdata.Modifier)
                this.event_data.selectNextEvent();
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
            elseif strcmp(eventdata.Key, 'delete') || strcmp(eventdata.Key, 'backspace')
                this.event_data.updateEventTime(this.event_data.selected_event_label, this.event_data.selected_event_time, NaN)
                this.updateEventPlots();
                this.updateSelectedEventPlots();
            elseif strcmp(eventdata.Key, 'escape')
                set(this.addLeftPushoffButton, 'ForegroundColor', [0 0 0]);
                set(this.addLeftTouchdownButton, 'ForegroundColor', [0 0 0]);
                set(this.addRightPushoffButton, 'ForegroundColor', [0 0 0]);
                set(this.addRightTouchdownButton, 'ForegroundColor', [0 0 0]);
                
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
        
        function deleteEvent(this, sender, eventdata)
            % this should eventually be moved to WalkingEventData
            
            selected_event_time_current = this.selected_event_time;
            
            this.selected_event_time = this.event_data.updateEventTime(this.selected_event_label, selected_event_time_current, NaN);
            
        end
        function addEvent(this, event_time)
            for i_button = 1 : length(this.addEventButtons)
                if get(this.addEventButtons(i_button), 'ForegroundColor') == this.color_selected
                    % this button was pressed to add an event
                    event_label = get(this.addEventButtons(i_button), 'UserData');
                    this.event_data.addEventTime(event_time, event_label);
                    this.updateEventPlots();
                    this.updateSelectedEventPlots();
                    set(this.addEventButtons(i_button), 'ForegroundColor', this.color_normal);
                end
            end
            
        end
        function addEventButtonPressed(this, sender, eventdata)
            % check if this button was pressed
            if get(sender, 'ForegroundColor') == this.color_normal
                % set all buttons to normal
                for i_button = 1 : length(this.addEventButtons)
                    set(this.addEventButtons(i_button), 'ForegroundColor', this.color_normal);
                end
                % set this button to pressed
                set(sender, 'ForegroundColor', this.color_selected);
                
                % change mouse cursors to crosshair
                for i_figure = 1 : length(this.figureSelectionBox.String)
                    set(this.figureSelectionBox.UserData{i_figure}.main_figure, 'Pointer', 'crosshair');
                end
            elseif get(sender, 'ForegroundColor') == this.color_selected
                set(sender, 'ForegroundColor', this.color_normal);
                % change mouse cursors to normal
                for i_figure = 1 : length(this.figureSelectionBox.String)
                    set(this.figureSelectionBox.UserData{i_figure}.main_figure, 'Pointer', 'arrow');
                end
            end            
            
        end
        function moveSelectedEvent(this, sender, eventdata)
            % this should eventually be moved to WalkingEventData
            
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
            
            % update
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