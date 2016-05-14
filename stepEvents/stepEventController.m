classdef stepEventController < handle
    properties
        trial_data;
        event_data;
        
        selected_event_label;
        selected_event_time;
        
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

            controller.saveFigureSettingsButton = uicontrol(figures_panel, 'Style', 'pushbutton', 'Position', [5, figures_panel_height-100, 130, 60], 'Fontsize', 12, 'String', 'Save Figure Settings');
            controller.loadFigureSettingsButton = uicontrol(figures_panel, 'Style', 'pushbutton', 'Position', [140, figures_panel_height-100, 130, 60], 'Fontsize', 12, 'String', 'Load Figure Settings');

            % event controls
            events_panel_height = 255;
            events_panel = uipanel(controller.control_figure, 'Title', 'Events Control', 'FontSize', 12, 'BackgroundColor', 'white', 'Units', 'pixels', 'Position', [5, figure_height-figures_panel_height-events_panel_height-5, figure_width-10, events_panel_height]);
            controller.findEventsButton = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [5, events_panel_height-75, 130, 60], 'Fontsize', 12, 'String', 'Find Events');
            controller.saveEventsButton = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [140, events_panel_height-75, 130, 60], 'Fontsize', 12, 'String', 'Save Events');

            uicontrol(events_panel, 'Style', 'text', 'string', 'Heel pos peak prominence (m):', 'Position', [5, events_panel_height-100, 190, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
            controller.heel_pos_peak_width = uicontrol(events_panel, 'Style', 'edit', 'BackgroundColor', 'white', 'Position', [150, events_panel_height-97, 40, 20], 'String', '0.05');
            uicontrol(events_panel, 'Style', 'text', 'string', 'Toes vel peak prominence (m):', 'Position', [5, events_panel_height-120, 190, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
            controller.toes_vel_peak_width = uicontrol(events_panel, 'Style', 'edit', 'BackgroundColor', 'white', 'Position', [150, events_panel_height-117, 40, 20], 'String', '0.05');
        end
        function setSelectedEvent(this, event_label, event_time)
            this.selected_event_label = event_label;
            this.selected_event_time = event_time;
            
            this.updateSelectedEventPlots();
        end
        function updateSelectedEventPlots(this)
            for i_figure = 1 : length(this.figureSelectionBox.String)
                this.figureSelectionBox.UserData{i_figure}.updateSelectedEventPlot();
            end
        end
        
        function processKeyPress(this, sender, eventdata)
            if strcmp(eventdata.Key, 'a')
                this.selectPreviousEvent();
                this.updateSelectedEventPlots();
            elseif strcmp(eventdata.Key, 'd')
                this.selectNextEvent();
                this.updateSelectedEventPlots();
            end
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
        
    end
    
end