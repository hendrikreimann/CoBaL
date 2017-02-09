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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% compare the kinematic tree against the kinematic chain

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
        addIgnoreMarkerButton;
        lockAddEventsButton
        findStretchesButton;
        saveStretchesButton;
        automateStretchFindingButton;
        previousTrialButton
        nextTrialButton
        quitButton
        
        condition_label
        trial_number_label
        
        figureSelectionBox;
        selected_time_edit;
        
        scene_figure = [];
        
        kinematic_tree_controller = [];
        
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
            events_panel_height = 260;
            events_panel = uipanel(this.control_figure, 'Title', 'Events Control', 'FontSize', 12, 'BackgroundColor', 'white', 'Units', 'pixels', 'Position', [5, figure_height-figures_panel_height-events_panel_height-5, figure_width-10, events_panel_height]);
            this.findEventsButton = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [5, events_panel_height-75, 130, 60], 'Fontsize', 12, 'String', 'Find Events', 'callback', @this.findEvents);
            this.saveEventsButton = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [140, events_panel_height-75, 130, 60], 'Fontsize', 12, 'String', 'Save Events', 'callback', @event_data.saveEvents);
            this.lockAddEventsButton = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [275, events_panel_height-75, 130, 60], 'Fontsize', 12, 'String', 'Lock Add', 'callback', @this.lockAddEvents);
            
            this.addEventButtons(1) = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [5, events_panel_height-135, 100, 60], 'Fontsize', 12, 'String', '<html><center>Add Left<br></center>Pushoff', 'callback', @this.addEventButtonPressed, 'UserData', 'left_pushoff');
            this.addEventButtons(2) = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [105, events_panel_height-135, 100, 60], 'Fontsize', 12, 'String', '<html><center>Add Left<br></center>Touchdown', 'callback', @this.addEventButtonPressed, 'UserData', 'left_touchdown');
            this.addEventButtons(3) = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [205, events_panel_height-135, 100, 60], 'Fontsize', 12, 'String', '<html><center>Add Right<br></center>Pushoff', 'callback', @this.addEventButtonPressed, 'UserData', 'right_pushoff');
            this.addEventButtons(4)= uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [305, events_panel_height-135, 100, 60], 'Fontsize', 12, 'String', '<html><center>Add Right<br></center>Touchdown', 'callback', @this.addEventButtonPressed, 'UserData', 'right_touchdown');

            this.addEventButtons(5) = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [5, events_panel_height-195, 130, 60], 'Fontsize', 12, 'String', 'Add Ignore Marker', 'callback', @this.addEventButtonPressed, 'UserData', 'ignore');
            
            this.findStretchesButton = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [5, events_panel_height-255, 130, 60], 'Fontsize', 12, 'String', 'Find Stretches', 'callback', @this.findStretches);
            this.saveStretchesButton = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [140, events_panel_height-255, 130, 60], 'Fontsize', 12, 'String', 'Save Stretches', 'callback', @this.saveStretches);
            this.automateStretchFindingButton = uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [275, events_panel_height-255, 130, 60], 'Fontsize', 12, 'String', 'Auto Save', 'callback', @this.automateStretchFinding);
            
            % figure control
            files_panel_height = 95;
            files_panel = uipanel(this.control_figure, 'Title', 'Files', 'FontSize', 12, 'BackgroundColor', 'white', 'Units', 'pixels', 'Position', [5, figure_height-figures_panel_height-events_panel_height-files_panel_height-5, figure_width-10, files_panel_height]);
            uicontrol(files_panel, 'Style', 'text', 'Position', [5, files_panel_height-32, 70, 15], 'Fontsize', 12, 'String', 'Condition:', 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
            this.condition_label = uicontrol(files_panel, 'Style', 'text', 'Position', [75, files_panel_height-32, 80, 15], 'Fontsize', 12, 'String', this.trial_data.condition, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
            uicontrol(files_panel, 'Style', 'text', 'Position', [195, files_panel_height-32, 50, 15], 'Fontsize', 12, 'String', 'Trial:', 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
            this.trial_number_label = uicontrol(files_panel, 'Style', 'text', 'Position', [235, files_panel_height-32, 80, 15], 'Fontsize', 12, 'String', num2str(this.trial_data.trial_number), 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
            this.previousTrialButton = uicontrol(files_panel, 'Style', 'pushbutton', 'Position', [5, figures_panel_height-100, 130, 60], 'Fontsize', 12, 'String', 'Previous Trial', 'callback', @this.loadPreviousTrial);
            this.nextTrialButton = uicontrol(files_panel, 'Style', 'pushbutton', 'Position', [140, figures_panel_height-100, 130, 60], 'Fontsize', 12, 'String', 'Next Trial', 'callback', @this.loadNextTrial);
            this.quitButton = uicontrol(files_panel, 'Style', 'pushbutton', 'Position', [275, figures_panel_height-100, 130, 60], 'Fontsize', 12, 'String', 'Quit', 'callback', @this.quit);
            
            % scene control
            scene_panel_height = 60;
            scene_panel = uipanel(this.control_figure, 'Title', 'Files', 'FontSize', 12, 'BackgroundColor', 'white', 'Units', 'pixels', 'Position', [5, figure_height-figures_panel_height-events_panel_height-files_panel_height-scene_panel_height-5, figure_width-10, scene_panel_height]);
            uicontrol(scene_panel, 'Style', 'text', 'string', 'Selected time:', 'Position', [5, scene_panel_height-40, 190, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
            this.selected_time_edit = uicontrol(scene_panel, 'Style', 'edit', 'BackgroundColor', 'white', 'Position', [80, scene_panel_height-37, 40, 20], 'String', '0');
        end
        
        
        function setSelectedEvent(this, event_label, event_time)
            if nargin == 1
                event_label = 'left_touchdown';
                event_times = this.event_data.getEventTimes(event_label);
                event_time = event_times(1);
            end
            this.event_data.selected_event_label = event_label;
            this.event_data.selected_event_time = event_time;
            this.trial_data.selected_time = event_time;
            
            this.updateSelectedEventPlots();
            this.updateSelectedTime();
        end
        function updateSelectedEventPlots(this)
            for i_figure = 1 : size(this.figureSelectionBox.String, 1)
                this.figureSelectionBox.UserData{i_figure}.updateSelectedEventPlot();
            end
        end
        function updateSelectedTime(this)
            % determine index to display
            [~, index_mocap] = min(abs(this.trial_data.time_marker - this.trial_data.selected_time));
            
            % update step event figures
            for i_figure = 1 : size(this.figureSelectionBox.String, 1)
                this.figureSelectionBox.UserData{i_figure}.updateSelectedTimePlot();
            end
            
            % extract marker data
            marker_data = this.trial_data.marker_positions(index_mocap, :);
            
            % extract joint center data
            if ~isempty(this.trial_data.joint_center_positions)
                joint_center_data = this.trial_data.joint_center_positions(index_mocap, :);
            else 
                joint_center_data = [];
            end                
            
            % extract CoM data
            if ~isempty(this.trial_data.com_positions)
                com_data = this.trial_data.com_positions(index_mocap, :);
            else 
                com_data = [];
            end
            
            % extract joint angle data and update
            if ~isempty(this.trial_data.joint_angles)
                joint_angle_data = this.trial_data.joint_angles(index_mocap, :);
            else 
                joint_angle_data = [];
            end
            
            % update scene figure
            this.scene_figure.update([marker_data joint_center_data com_data]);
            
            % update kinematic chain stick figure
            if ~isempty(this.kinematic_tree_controller)
                this.kinematic_tree_controller.kinematicTree.jointAngles = joint_angle_data';
                this.kinematic_tree_controller.kinematicTree.updateConfiguration();
                this.kinematic_tree_controller.update();
                
                this.kinematic_tree_controller.updateRecordedMarkerPlots([marker_data joint_center_data]);
            end
            
            this.selected_time_edit.String = num2str(this.trial_data.selected_time);
        end
        function updateEventPlots(this)
            for i_figure = 1 : size(this.figureSelectionBox.String, 1)
                this.figureSelectionBox.UserData{i_figure}.updateEventPlots();
            end
        end
        function updateDataPlots(this)
            for i_figure = 1 : size(this.figureSelectionBox.String, 1)
                this.figureSelectionBox.UserData{i_figure}.updateDataPlots();
            end
        end
        function updateStretchPatches(this)
            for i_figure = 1 : size(this.figureSelectionBox.String, 1)
                this.figureSelectionBox.UserData{i_figure}.updateStretchPatches();
            end
        end            
        function updateTrialLabels(this)
            this.condition_label.String = this.trial_data.condition;
            this.trial_number_label.String = num2str(this.trial_data.trial_number);
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
                this.updateSelectedTime();
            end
            if strcmp(eventdata.Key, 'd') && isempty(eventdata.Modifier)
                this.event_data.selectNextEvent();
                this.updateSelectedEventPlots();
                this.updateSelectedTime();
            end
            if strcmp(eventdata.Key, 'w') && isempty(eventdata.Modifier)
                this.updateTimeWindow('zoom in');
            end
            if strcmp(eventdata.Key, 'x') && isempty(eventdata.Modifier)
                this.updateTimeWindow('zoom out');
            end
            if strcmp(eventdata.Key, 's') && isempty(eventdata.Modifier)
                this.updateTimeWindow('center');
            end
            if strcmp(eventdata.Key, 'q') && isempty(eventdata.Modifier)
                this.updateTimeWindow('previous');
            end
            if strcmp(eventdata.Key, 'e') && isempty(eventdata.Modifier)
                this.updateTimeWindow('next');
            end
            if strcmp(eventdata.Key, 'leftarrow')
                if isempty(eventdata.Modifier)
                    this.trial_data.stepSelectedTime('back')
                elseif strcmp(eventdata.Modifier, 'shift')
                    this.trial_data.stepSelectedTime('back', 5)
                elseif strcmp(eventdata.Modifier, 'alt')
                    this.trial_data.stepSelectedTime('back', 25)
                end
                this.updateSelectedTime();
            end
            if strcmp(eventdata.Key, 'rightarrow')
                if isempty(eventdata.Modifier)
                    this.trial_data.stepSelectedTime('forward')
                elseif strcmp(eventdata.Modifier, 'shift')
                    this.trial_data.stepSelectedTime('forward', 5)
                elseif strcmp(eventdata.Modifier, 'alt')
                    this.trial_data.stepSelectedTime('forward', 25)
                end
                this.updateSelectedTime();
            end
            if strcmp(eventdata.Key, 'z') || strcmp(eventdata.Key, 'c')
                this.moveSelectedEvent(sender, eventdata);
            end
            if strcmp(eventdata.Modifier, 'command') & strcmp(eventdata.Key, 'x')
                % this won't work on windows systems, adapt that!
                this.quit();
            end
            if strcmp(eventdata.Key, 'delete') || strcmp(eventdata.Key, 'backspace')
                this.event_data.updateEventTime(this.event_data.selected_event_label, this.event_data.selected_event_time, NaN);
                this.updateEventPlots();
                this.updateSelectedEventPlots();
                this.updateSelectedTime();
                
                if get(this.automateStretchFindingButton, 'ForegroundColor') == this.color_selected
                    % automatic processing requested
                    this.event_data.saveEvents;
                    this.findStretches(sender, eventdata);
                end            
            end  
            if strcmp(eventdata.Key, 'escape')
                set(this.addLeftPushoffButton, 'ForegroundColor', [0 0 0]);
                set(this.addLeftTouchdownButton, 'ForegroundColor', [0 0 0]);
                set(this.addRightPushoffButton, 'ForegroundColor', [0 0 0]);
                set(this.addRightTouchdownButton, 'ForegroundColor', [0 0 0]);
                
            end
        end
        function saveFigureSettings(this, sender, eventdata)
            
            figure_settings = cell(1, length(this.figureSelectionBox));
            for i_figure = 1 : size(this.figureSelectionBox.String, 1)
                figure_settings{i_figure} = this.figureSelectionBox.UserData{i_figure}.getSetting();
            end
            
            scene_figure_setting = struct();
            scene_figure_setting.position = this.scene_figure.scene_figure.Position;
            
            control_figure_setting = struct();
            control_figure_setting.position = this.control_figure.Position;
            
            kinematic_tree_figure_setting = struct();
            if ~isempty(this.kinematic_tree_controller)
                kinematic_tree_figure_setting.position = this.kinematic_tree_controller.sceneFigure.Position;
            end
            settings_file = [getUserSettingsPath filesep 'eventGuiFigureSettings.mat'];
            save(settings_file, 'figure_settings', 'control_figure_setting', 'scene_figure_setting', 'kinematic_tree_figure_setting');
        end
        function loadFigureSettings(this, sender, eventdata)
%             settings_file = '/Users/reimajbi/Library/Application Support/stepEventFinder/figureSettings.mat';
            settings_file = [getUserSettingsPath filesep 'eventGuiFigureSettings.mat'];

            if exist(settings_file, 'file')
                load(settings_file, 'figure_settings', 'control_figure_setting', 'scene_figure_setting', 'kinematic_tree_figure_setting')
                for i_figure = 1 : length(figure_settings)
        %             if length(step_event_figures) < i_figure
        %                 step_event_figures{i_figure} = createStepEventFigure();
        %             end
                    this.figureSelectionBox.UserData{i_figure}.applySettings(figure_settings{i_figure});
                end
                this.control_figure.Position = control_figure_setting.position;
                this.scene_figure.scene_figure.Position = scene_figure_setting.position;
                if ~isempty(this.kinematic_tree_controller) && any(strcmp(fieldnames(kinematic_tree_figure_setting), 'position'))
                    this.kinematic_tree_controller.sceneFigure.Position = kinematic_tree_figure_setting.position;
                end
            end
        end
        function findEvents(this, sender, eventdata)
            % find events
            findStepEvents('condition', this.trial_data.condition, 'trials', this.trial_data.trial_number);
            
            % load results from events file
            this.event_data.loadEvents;

            % update plots
            this.updateEventPlots;
            this.event_data.selectNextEvent;
            this.updateSelectedEventPlots;
            this.updateSelectedTime
            
            this.updateEventPlots();
            
            if get(this.automateStretchFindingButton, 'ForegroundColor') == this.color_selected
                % automatic processing requested
                this.findStretches(sender, eventdata);
            end            
        end
        
        function findStretches(this, sender, eventdata)
            % find stretches
            findRelevantDataStretches('condition', this.trial_data.condition, 'trials', this.trial_data.trial_number);
            
            % load results from stretch file
            this.event_data.loadStretches;
            
            % update plots
            this.updateStretchPatches;
        end
        function saveStretches(this, sender, eventdata)
            disp('Pushing the "saveStretches" button does nothing, will be removed')
        end
        function automateStretchFinding(this, sender, eventdata)
            if get(sender, 'ForegroundColor') == this.color_normal
                set(sender, 'ForegroundColor', this.color_selected);
                
                % process changes that have been made but might not have been saved
                this.event_data.saveEvents;
                this.findStretches(sender, eventdata)
            else
                set(sender, 'ForegroundColor', this.color_normal);
            end
            
        end
        
        function addEvent(this, event_time)
            for i_button = 1 : length(this.addEventButtons)
                if get(this.addEventButtons(i_button), 'ForegroundColor') == this.color_selected
                    % this button was pressed to add an event
                    event_label = get(this.addEventButtons(i_button), 'UserData');
                    this.event_data.addEventTime(event_time, event_label);
                    this.updateEventPlots();
                    this.updateSelectedEventPlots();
                    if get(this.lockAddEventsButton, 'ForegroundColor') == this.color_normal
                        set(this.addEventButtons(i_button), 'ForegroundColor', this.color_normal);
                        % change mouse cursors to normal
                        for i_figure = 1 : size(this.figureSelectionBox.String, 1)
                            set(this.figureSelectionBox.UserData{i_figure}.main_figure, 'Pointer', 'arrow');
                        end
                    end
                end
            end
            if get(this.automateStretchFindingButton, 'ForegroundColor') == this.color_selected
                % automatic processing requested
                this.event_data.saveEvents;
                this.findStretches();
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
                for i_figure = 1 : size(this.figureSelectionBox.String, 1)
                    set(this.figureSelectionBox.UserData{i_figure}.main_figure, 'Pointer', 'crosshair');
                end
            elseif get(sender, 'ForegroundColor') == this.color_selected
                set(sender, 'ForegroundColor', this.color_normal);
                % change mouse cursors to normal
                for i_figure = 1 : size(this.figureSelectionBox.String, 1)
                    set(this.figureSelectionBox.UserData{i_figure}.main_figure, 'Pointer', 'arrow');
                end
            end
        end
        function lockAddEvents(this, sender, eventdata)
            if get(sender, 'ForegroundColor') == this.color_normal
                set(sender, 'ForegroundColor', this.color_selected);
            else
                set(sender, 'ForegroundColor', this.color_normal);
            end
            
        end
        function moveSelectedEvent(this, sender, eventdata)
            % this should eventually be moved to WalkingEventData
            
            selected_event_time_current = this.event_data.selected_event_time;
            
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
            this.event_data.selected_event_time = this.event_data.updateEventTime(this.event_data.selected_event_label, selected_event_time_current, selected_event_time_new);
            
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
                new_center = this.event_data.selected_event_time;
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
                new_x_lim = new_x_lim - (new_x_lim(2) - this.trial_data.recording_time);
            end
            
            for i_figure = 1 : size(this.figureSelectionBox.String, 1)
                set(this.figureSelectionBox.UserData{i_figure}.main_axes, 'xlim', new_x_lim);
            end        
        end
        
        function loadPreviousTrial(this, sender, eventdata)
            current_condition = this.trial_data.condition;
            current_trial_number = this.trial_data.trial_number;
            
            [condition_list, trial_number_list] = parseTrialArguments();
            condition_index = find(strcmp(condition_list, current_condition));
            if current_trial_number == trial_number_list{condition_index}(1)
                new_condition = condition_list{condition_index - 1};
                new_trial_number = trial_number_list{condition_index - 1}(end);
            else
                new_condition = this.trial_data.condition;
                trial_index = find(trial_number_list{condition_index} == current_trial_number);
                new_trial_number = trial_number_list{condition_index}(trial_index - 1);
            end
            this.trial_data.condition = new_condition;
            this.trial_data.trial_number = new_trial_number;
            this.trial_data.loadMarkerTrajectories;
            this.trial_data.loadForceplateTrajectories;
            this.event_data.loadEvents;
            this.event_data.loadStretches;
            
            this.updateDataPlots;
            this.updateEventPlots;
            this.updateStretchPatches;
            this.event_data.selectNextEvent;
            this.updateSelectedEventPlots;
            this.updateTrialLabels;
        end
        function loadNextTrial(this, sender, eventdata)
            current_condition = this.trial_data.condition;
            current_trial_number = this.trial_data.trial_number;
            
            [condition_list, trial_number_list] = parseTrialArguments();
            condition_index = find(strcmp(condition_list, current_condition));
            if current_trial_number == trial_number_list{condition_index}(end)
                new_condition = condition_list{condition_index + 1};
                new_trial_number = trial_number_list{condition_index + 1}(1);
            else
                new_condition = this.trial_data.condition;
                trial_index = find(trial_number_list{condition_index} == current_trial_number);
                new_trial_number = trial_number_list{condition_index}(trial_index + 1);
            end
            this.trial_data.condition = new_condition;
            this.trial_data.trial_number = new_trial_number;
            this.trial_data.loadMarkerTrajectories;
            this.trial_data.loadForceplateTrajectories;
            this.event_data.loadEvents;
            this.event_data.loadStretches;

            this.updateDataPlots;
            this.updateEventPlots;
            this.updateStretchPatches;
            this.event_data.selectNextEvent;
            this.updateSelectedEventPlots;
            this.updateTrialLabels;
            
            
            
            this.updateSelectedTime
        end
        function quit(this, sender, eventdata)
            try
                close(this.scene_figure.scene_figure);
            catch exception
                if ~strcmp(exception.identifier, 'MATLAB:close:InvalidFigureHandle')
                    rethrow(exception)
                end
            end
            for i_figure = 1 : size(this.figureSelectionBox.String, 1)
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