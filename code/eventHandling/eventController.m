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

classdef eventController < handle
    properties
        data_custodian;
        event_data;
%         figure_settings_file;
        
        event_time_normal_step = 0.005;
        event_time_large_step = 0.050;
        
        control_figure;
        saveFigureSettingsButton;
        loadFigureSettingsButton;
        toggleLegendsButton;
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
        event_time_update_button;
        event_type_text;
        
        scene_figure = [];
        kinematic_tree_stick_figure = [];
        kinematic_tree_controller = [];
        
        color_selected = [1 0.5 0];
        color_normal = [0 0 0];
        
        canvas_on_screen_width;
        canvas_on_screen_height;
    end
    methods
        function this = eventController(data_custodian, event_data)
            this.data_custodian = data_custodian;
            this.event_data = event_data;
            
            screen_size = get(0,'ScreenSize');
            figure_height = 600;
            figure_width = 420;
            this.canvas_on_screen_width = screen_size(3) - figure_width;
            this.canvas_on_screen_height = screen_size(4) - 24; % 24 is the size of the menu bar on the Mac, check for windows systems later
            if ismac
                this.control_figure = figure ...
                  ( ...
                    'position', [screen_size(3)-figure_width screen_size(3)-figure_height figure_width figure_height], ...
                    'Units', 'pixels', ...
                    'KeyPressFcn', @this.processKeyPress ...
                  );
            elseif ispc
                this.control_figure = figure('position', [0 0 figure_width figure_height], 'Units', 'pixels', 'KeyPressFcn', @this.processKeyPress);
            end
            % figure control
            figures_panel_height = 100;
            figures_panel = uipanel ...
              ( ...
                this.control_figure, ...
                'Title', 'Figure Control', ...
                'FontSize', 12, ...
                'BackgroundColor', 'white', ...
                'Units', 'pixels', ...
                'Position', [5, figure_height-figures_panel_height-5, figure_width-10, figures_panel_height] ...
              );
            this.figureSelectionBox = uicontrol ...
              ( ...
                figures_panel, ...
                'Style', 'popup', ...
                'String', '<no figure>', ...
                'Position', [5, figures_panel_height-40, 395, 20], ...
                'Fontsize', 12, ...
                'HorizontalAlignment', 'left', ...
                'KeyPressFcn', @this.processKeyPress ...
              );

            this.toggleLegendsButton = uicontrol ...
              ( ...
                figures_panel, ...
                'Style', 'pushbutton', ...
                'Position', [275, figures_panel_height-100, 130, 60], ...
                'Fontsize', 12, ...
                'String', 'Toggle Legends', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'callback', @this.toggleLegends ...
              );

            % event controls
            events_panel_height = 260;
            events_panel = uipanel ...
              ( ...
                this.control_figure, ...
                'Title', 'Events Control', ...
                'FontSize', 12, ...
                'BackgroundColor', 'white', ...
                'Units', 'pixels', ...
                'Position', [5, figure_height-figures_panel_height-events_panel_height-5, figure_width-10, events_panel_height] ...
              );
            this.findEventsButton = uicontrol ...
              ( ...
                events_panel, ...
                'Style', 'pushbutton', ...
                'Position', [5, events_panel_height-75, 130, 60], 'Fontsize', 12, 'String', 'Find Events', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'callback', @this.findEvents ...
              );
            this.saveEventsButton = uicontrol ...
              ( ...
                events_panel, ...
                'Style', 'pushbutton', ...
                'Position', [140, events_panel_height-75, 130, 60], ...
                'Fontsize', 12, ...
                'String', 'Save Events', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'callback', @event_data.saveEvents ...
              );
            this.lockAddEventsButton = uicontrol ...
              ( ...
                events_panel, ...
                'Style', 'pushbutton', ...
                'Position', [275, events_panel_height-75, 130, 60], ...
                'Fontsize', 12, ...
                'String', 'Lock Add', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'callback', @this.lockAddEvents ...
              );
            
            this.addEventButtons(1) = uicontrol ...
              ( ...
                events_panel, ...
                'Style', 'pushbutton', ...
                'Position', [5, events_panel_height-135, 100, 60], ...
                'Fontsize', 12, ...
                'String', '<html><center>Add Left<br></center>Pushoff', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'callback', @this.addEventButtonPressed, ...
                'UserData', 'left_pushoff' ...
              );
            this.addEventButtons(2) = uicontrol ...
              ( ...
                events_panel, ...
                'Style', 'pushbutton', ...
                'Position', [105, events_panel_height-135, 100, 60], ...
                'Fontsize', 12, ...
                'String', '<html><center>Add Left<br></center>Touchdown', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'callback', @this.addEventButtonPressed, ...
                'UserData', 'left_touchdown' ...
              );
            this.addEventButtons(3) = uicontrol ...
              ( ...
                events_panel, ...
                'Style', 'pushbutton', ...
                'Position', [205, events_panel_height-135, 100, 60], ...
                'Fontsize', 12, ...
                'String', '<html><center>Add Right<br></center>Pushoff', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'callback', @this.addEventButtonPressed, ...
                'UserData', 'right_pushoff' ...
              );
            this.addEventButtons(4)= uicontrol ...
              ( ...
                events_panel, ...
                'Style', 'pushbutton', ...
                'Position', [305, events_panel_height-135, 100, 60], ...
                'Fontsize', 12, ...
                'String', '<html><center>Add Right<br></center>Touchdown', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'callback', @this.addEventButtonPressed, ...
                'UserData', 'right_touchdown' ...
              );

            this.addEventButtons(5) = uicontrol ...
              ( ...
                events_panel, ...
                'Style', 'pushbutton', ...
                'Position', [5, events_panel_height-195, 130, 60], ...
                'Fontsize', 12, ...
                'String', 'Add Problem Marker', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'callback', @this.addEventButtonPressed, ...
                'UserData', 'problem' ...
              );
            
            this.findStretchesButton = uicontrol ...
              ( ...
                events_panel, ...
                'Style', 'pushbutton', ...
                'Position', [5, events_panel_height-255, 130, 60], ...
                'Fontsize', 12, ...
                'String', 'Find Stretches', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'callback', @this.findStretches ...
              );
            this.automateStretchFindingButton = uicontrol ...
              ( ...
                events_panel, ...
                'Style', 'pushbutton', ...
                'Position', [275, events_panel_height-255, 130, 60], ...
                'Fontsize', 12, ...
                'String', '<html><center>Auto-Find<br></center>Stretches', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'callback', @this.automateStretchFinding, ...
                'ForegroundColor', this.color_selected ...
              );
            
            % figure control
            files_panel_height = 95;
            files_panel = uipanel ...
              ( ...
                this.control_figure, ...
                'Title', 'Files', ...
                'FontSize', 12, ...
                'BackgroundColor', 'white', ...
                'Units', 'pixels', ...
                'Position', [5, figure_height-figures_panel_height-events_panel_height-files_panel_height-5, figure_width-10, files_panel_height] ...
              );
            uicontrol ...
              ( ... 
                files_panel, ...
                'Style', 'text', ...
                'Position', [5, files_panel_height-32, 70, 15], ...
                'Fontsize', 12, ...
                'String', 'Condition:', ...
                'HorizontalAlignment', 'left', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'BackgroundColor', 'white' ...
              );
            this.condition_label = uicontrol ...
              ( ...
                files_panel, ...
                'Style', 'text', ...
                'Position', [75, files_panel_height-32, 80, 15], ...
                'Fontsize', 12, ...
                'String', this.data_custodian.trial_type, ...
                'HorizontalAlignment', 'left', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'BackgroundColor', 'white' ...
              );
            uicontrol ...
              ( ...
                files_panel, ...
                'Style', 'text', ...
                'Position', [195, files_panel_height-32, 50, 15], ...
                'Fontsize', 12, ...
                'String', 'Trial:', ...
                'HorizontalAlignment', 'left', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'BackgroundColor', 'white' ...
              );
            this.trial_number_label = uicontrol ...
              ( ...
                files_panel, ...
                'Style', 'text', ...
                'Position', [235, files_panel_height-32, 80, 15], ...
                'Fontsize', 12, ...
                'String', num2str(this.data_custodian.trial_number), ...
                'HorizontalAlignment', 'left', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'BackgroundColor', 'white' ...
              );
            this.previousTrialButton = uicontrol ...
              ( ...
                files_panel, ...
                'Style', 'pushbutton', ...
                'Position', [5, figures_panel_height-100, 130, 60], ...
                'Fontsize', 12, ...
                'String', 'Previous Trial', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'callback', @this.loadPreviousTrial ...
              );
            this.nextTrialButton = uicontrol ...
              ( ...
                files_panel, ...
                'Style', 'pushbutton', ...
                'Position', [140, figures_panel_height-100, 130, 60], ...
                'Fontsize', 12, ...
                'String', 'Next Trial', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'callback', @this.loadNextTrial ...
              );
            this.quitButton = uicontrol ...
              ( ...
                files_panel, ...
                'Style', 'pushbutton', ...
                'Position', [275, figures_panel_height-100, 130, 60], ...
                'Fontsize', 12, ...
                'String', 'Quit', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'callback', @this.quit ...
              );
            
            % scene control
            scene_panel_height = 70;
            scene_panel = uipanel ...
              ( ...
                this.control_figure, ...
                'Title', 'Selection', ...
                'FontSize', 12, ...
                'BackgroundColor', 'white', ...
                'Units', 'pixels', ...
                'Position', [5, figure_height-figures_panel_height-events_panel_height-files_panel_height-scene_panel_height-5, figure_width-10, scene_panel_height] ...
              );
            uicontrol ...
              ( ...
                scene_panel, ...
                'Style', 'text', ...
                'string', 'Selected time:', ...
                'Position', [5, scene_panel_height-40, 190, 20], ...
                'Fontsize', 10, ...
                'HorizontalAlignment', 'left', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'BackgroundColor', 'white' ...
              );
            this.selected_time_edit = uicontrol ...
              ( ...
                scene_panel, ...
                'Style', 'edit', ...
                'BackgroundColor', 'white', ...
                'Position', [80, scene_panel_height-37, 40, 20], ...
                'KeyPressFcn', @this.processKeyPress, ...
                'String', '0' ...
              );
            this.event_time_update_button = uicontrol ...
              ( ...
                scene_panel, ...
                'Style', 'pushbutton', ...
                'Position', [130, scene_panel_height-37, 60, 20], ...
                'Fontsize', 12, ...
                'String', 'update', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'callback', @this.updateEventTime ...
              );
            this.event_type_text = uicontrol ...
              ( ...
                scene_panel, ...
                'Style', 'text', ...
                'string', 'Selected event type:', ...
                'Position', [5, scene_panel_height-65, 190, 20], ...
                'Fontsize', 10, ...
                'HorizontalAlignment', 'left', ...
                'KeyPressFcn', @this.processKeyPress, ...
                'BackgroundColor', 'white' ...
              );
        end
        
        function updateEventTime(this, varargin)
            selected_event_time_new = str2double(this.selected_time_edit.String);

            selected_event_time_current = this.event_data.selected_event_time;
            this.event_data.selected_event_time = this.event_data.updateEventTime(this.event_data.selected_event_label, selected_event_time_current, selected_event_time_new);
            
            % update plots
            this.updateSelectedEventPlots();
            this.updateSelectedTime(this.event_data.selected_event_time);
            this.updateEventPlots();
        end
        function setSelectedEvent(this, event_label, event_time)
            if nargin == 1
                event_label = 'left_touchdown';
                event_times = this.event_data.getEventTimes(event_label);
                event_time = event_times(1);
            end
            this.event_data.selected_event_label = event_label;
            this.event_data.selected_event_time = event_time;
            this.event_data.selected_time = event_time;
            
            
            set(this.event_type_text, 'string', ['Selected event type: ' event_label])
            
            this.updateSelectedEventPlots();
            this.updateSelectedTime();
        end
        function updateSelectedEventPlots(this)
            for i_figure = 1 : size(this.figureSelectionBox.String, 1)
                this.figureSelectionBox.UserData{i_figure}.updateSelectedEventPlot();
            end
        end
        function updateSelectedTime(this, new_time)
            if nargin > 1
                this.event_data.selected_time = new_time;
            end
            
            % update step event figures
            for i_figure = 1 : size(this.figureSelectionBox.String, 1)
                this.figureSelectionBox.UserData{i_figure}.updateSelectedTimePlot();
            end
            
            % update scene figures
            if ~isempty(this.scene_figure)
                
                % determine index to display
                [~, index_mocap] = min(abs(this.data_custodian.getTimeData('marker_trajectories') - this.event_data.selected_time));

                % extract marker data
                marker_trajectories = this.data_custodian.getBasicVariableData('marker_trajectories');
                marker_data = marker_trajectories(index_mocap, :);
%                 joint_center_trajectories = this.data_custodian.getBasicVariableData('joint_center_trajectories');
%                 joint_center_data = joint_center_trajectories(index_mocap, :);
%                 com_trajectories = this.data_custodian.getBasicVariableData('com_trajectories');
%                 com_data = com_trajectories(index_mocap, :);
                % TODO: deal with cases where we have only marker data and no kinematic data yet
            
%             
%             % extract joint angle data and update
%             if ~isempty(this.trial_data.joint_angles)
%                 joint_angle_data = this.trial_data.joint_angles(index_mocap, :);
%             else 
%                 joint_angle_data = [];
%             end
%             
%             % update scene figure
%                 this.scene_figure.update([marker_data joint_center_data com_data]);
%             
%             % update kinematic chain stick figure
%             if ~isempty(this.kinematic_tree_controller)
%                 this.kinematic_tree_stick_figure.kinematicTree.jointAngles = joint_angle_data';
%                 this.kinematic_tree_stick_figure.kinematicTree.updateConfiguration();
%                 this.kinematic_tree_stick_figure.update();
%                 
%                 this.kinematic_tree_stick_figure.updateRecordedMarkerPlots([marker_data joint_center_data]);
%             end
            
            end            
            
            this.selected_time_edit.String = num2str(this.event_data.selected_time);
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
            this.condition_label.String = this.data_custodian.trial_type;
            this.trial_number_label.String = num2str(this.data_custodian.trial_number);
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
            if strcmp(eventdata.Key, 'r') && isempty(eventdata.Modifier)
                this.loadPreviousTrial(sender, eventdata)
            end
            if strcmp(eventdata.Key, 't') && isempty(eventdata.Modifier)
                this.loadNextTrial(sender, eventdata)
            end
            if strcmp(eventdata.Key, 'leftarrow')
                if isempty(eventdata.Modifier)
                    this.event_data.stepSelectedTime('back')
                elseif strcmp(eventdata.Modifier, 'shift')
                    this.event_data.stepSelectedTime('back', 5)
                elseif strcmp(eventdata.Modifier, 'alt')
                    this.event_data.stepSelectedTime('back', 25)
                end
                this.updateSelectedTime();
            end
            if strcmp(eventdata.Key, 'rightarrow')
                if isempty(eventdata.Modifier)
                    this.event_data.stepSelectedTime('forward')
                elseif strcmp(eventdata.Modifier, 'shift')
                    this.event_data.stepSelectedTime('forward', 5)
                elseif strcmp(eventdata.Modifier, 'alt')
                    this.event_data.stepSelectedTime('forward', 25)
                end
                this.updateSelectedTime();
            end
            if strcmp(eventdata.Key, 'z') || strcmp(eventdata.Key, 'c')
                this.moveSelectedEvent(sender, eventdata);
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
        function saveFigureSettings(this, sender, eventdata) %#ok<INUSD>
            % get figure settings from trajectory figures
            figure_settings = cell(1, length(this.figureSelectionBox));
            for i_figure = 1 : size(this.figureSelectionBox.String, 1)
                figure_settings{i_figure} = this.figureSelectionBox.UserData{i_figure}.getSetting();
            end
            
            % get figure settings from scene figure
            scene_figure_setting = struct();
            scene_figure_setting.position = this.scene_figure.scene_figure.Position; %#ok<STRNU>
            
            % get figure settings from control figure
            control_figure_setting = struct();
            control_figure_setting.position = this.control_figure.Position; %#ok<STRNU>
            
            % get figure settings from kinematic tree figure
            kinematic_tree_figure_setting = struct();
            if ~isempty(this.kinematic_tree_stick_figure)
                kinematic_tree_figure_setting.position = this.kinematic_tree_stick_figure.sceneFigure.Position; %#ok<STRNU>
            end
            
            % save settings to file
%             settings_file = [getUserSettingsPath filesep this.figure_settings_file];
%             save(settings_file, 'figure_settings', 'control_figure_setting', 'scene_figure_setting', 'kinematic_tree_figure_setting');
        end
        function loadFigureSettings(this, sender, eventdata) %#ok<INUSD>
%             settings_file = [getUserSettingsPath filesep this.figure_settings_file];
% 
%             if exist(settings_file, 'file')
%                 % load settings
%                 load(settings_file, 'figure_settings', 'control_figure_setting', 'scene_figure_setting', 'kinematic_tree_figure_setting')
%                 
%                 % apply for controller and stick figure
%                 this.control_figure.Position = control_figure_setting.position;
%                 if ~isempty(this.scene_figure)
%                     this.scene_figure.scene_figure.Position = scene_figure_setting.position;
%                 end
%                 if ~isempty(this.kinematic_tree_stick_figure) && any(strcmp(fieldnames(kinematic_tree_figure_setting), 'position'))
%                     this.kinematic_tree_stick_figure.sceneFigure.Position = kinematic_tree_figure_setting.position;
%                 end
%                 
%                 % apply for trajectory figure
%                 for i_figure = 1 : length(figure_settings) %#ok<USENS>
%                     this_figure_settings = figure_settings{i_figure};
%                     this_figure_title = this_figure_settings.title;
%                     
%                     % cycle through available figures and look for a match
%                     for j_figure = 1 : length(this.figureSelectionBox.UserData)
%                         if strcmp(this.figureSelectionBox.UserData{j_figure}.title, this_figure_title)
%                             this.figureSelectionBox.UserData{j_figure}.applySettings(this_figure_settings);
%                         end
%                     end
%                 end
%             end
        end
        function toggleLegends(this, sender, eventdata) %#ok<INUSD>
            for i_figure = 1 : size(this.figureSelectionBox.String, 1)
                this.figureSelectionBox.UserData{i_figure}.toggleLegend;
            end
        end
        
        function findEvents(this, sender, eventdata)
            % find events
            if strcmp(this.data_custodian.study_settings.get('study_type', 1), 'MS')
                findEvents_MS('condition', this.data_custodian.trial_type, 'trials', this.data_custodian.trial_number);
            else
                findStepEvents('condition', this.data_custodian.trial_type, 'trials', this.data_custodian.trial_number);
            end
            
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
        
        function findStretches(this, sender, eventdata) %#ok<INUSD>
            if strcmp(this.data_custodian.study_settings.get('study_type', 1), 'MS')
                return
            end
            
            % find stretches
            determineStretchesToAnalyze('condition', this.data_custodian.trial_type, 'trials', this.data_custodian.trial_number);
            
            
            % load results from stretch file
            this.event_data.loadStretches;
            
            % update plots
            this.updateStretchPatches;
        end
        function saveStretches(this, sender, eventdata) %#ok<INUSD>
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
        function addEventButtonPressed(this, sender, eventdata) %#ok<INUSD>
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
        function lockAddEvents(this, sender, eventdata) %#ok<INUSD>
            if get(sender, 'ForegroundColor') == this.color_normal
                set(sender, 'ForegroundColor', this.color_selected);
            else
                set(sender, 'ForegroundColor', this.color_normal);
            end
            
        end
        function moveSelectedEvent(this, sender, eventdata) %#ok<INUSL>
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
            this.updateSelectedTime(this.event_data.selected_event_time);
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
            if new_x_lim(1) < this.data_custodian.getRecordingTimeStart()
                new_x_lim = new_x_lim - (new_x_lim(1) - this.data_custodian.getRecordingTimeStart());
            end
            if new_x_lim(2) > this.data_custodian.getRecordingTimeEnd()
                new_x_lim = new_x_lim - (new_x_lim(2) - this.data_custodian.getRecordingTimeEnd());
            end
            
            for i_figure = 1 : size(this.figureSelectionBox.String, 1)
                this.figureSelectionBox.UserData{i_figure}.updateTimeWindow(new_x_lim);
%                 set(this.figureSelectionBox.UserData{i_figure}.main_axes, 'xlim', new_x_lim);
            end        
        end
        
        function loadPreviousTrial(this, sender, eventdata) %#ok<INUSD>
            current_condition = this.data_custodian.trial_type;
            current_trial_number = this.data_custodian.trial_number;
            
            [condition_list, trial_number_list] = parseTrialArguments();
            condition_index = find(strcmp(condition_list, current_condition));
            if current_trial_number == trial_number_list{condition_index}(1)
                % first trial of this type
                if condition_index == 1
                    % first type, so just stay at the current trial
                    new_condition = current_condition;
                    new_trial_number = current_trial_number;
                else
                    % go to next trial type, first trial
                    new_condition = condition_list{condition_index - 1};
                    new_trial_number = trial_number_list{condition_index - 1}(end);
                end                    
            else
                % go to previous trial of this type
                new_condition = this.data_custodian.trial_type;
                trial_index = find(trial_number_list{condition_index} == current_trial_number);
                new_trial_number = trial_number_list{condition_index}(trial_index - 1);
            end
            
            % load new data
            this.data_custodian.prepareBasicVariables(new_condition, new_trial_number)
            this.event_data.loadEvents;
            this.event_data.loadStretches;
            
            % update everything
            this.updateDataPlots;
            this.updateEventPlots;
            this.updateStretchPatches;
            this.event_data.selectClosestEvent;
            this.updateSelectedEventPlots;
            this.updateTrialLabels;
            this.updateSelectedTime;
        end
        function loadNextTrial(this, sender, eventdata) %#ok<INUSD>
            % figure out new type/number based on current values
            current_condition = this.data_custodian.trial_type;
            current_trial_number = this.data_custodian.trial_number;
            
            [condition_list, trial_number_list] = parseTrialArguments();
            condition_index = find(strcmp(condition_list, current_condition));
            if current_trial_number == trial_number_list{condition_index}(end)
                % last trial of this type
                if condition_index == length(condition_list)
                    % last type, so just stay at the current trial
                    new_condition = current_condition;
                    new_trial_number = current_trial_number;
                else
                    % go to next trial type, first trial
                    new_condition = condition_list{condition_index + 1};
                    new_trial_number = trial_number_list{condition_index + 1}(1);
                end
            else
                % go to next trial of this type
                new_condition = this.data_custodian.trial_type;
                trial_index = find(trial_number_list{condition_index} == current_trial_number);
                new_trial_number = trial_number_list{condition_index}(trial_index + 1);
            end
            
            % load new data
            this.data_custodian.prepareBasicVariables(new_condition, new_trial_number)
            this.event_data.loadEvents;
            this.event_data.loadStretches;

            % update everything
            this.updateDataPlots;
            this.updateEventPlots;
            this.updateStretchPatches;
            this.event_data.selectClosestEvent;
            this.updateSelectedEventPlots;
            this.updateTrialLabels;
            this.updateSelectedTime;
        end
        function quit(this, sender, eventdata) %#ok<INUSD>
            try
                close(this.scene_figure.scene_figure);
            catch exception
                if strcmp(exception.identifier, 'MATLAB:close:InvalidFigureHandle')
                    % do nothing
                elseif strcmp(exception.identifier, 'MATLAB:structRefFromNonStruct')
                    % do nothing
                else
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
            try
                close(this.kinematic_tree_stick_figure.sceneFigure);
            catch exception
                if ~ ...
                     ( ...
                       strcmp(exception.identifier, 'MATLAB:close:InvalidFigureHandle') ...
                       || strcmp(exception.identifier, 'MATLAB:structRefFromNonStruct') ...
                     )
                    rethrow(exception)
                end
            end
            try
                close(this.kinematic_tree_controller);
            catch exception
                if ~strcmp(exception.identifier, 'MATLAB:close:InvalidFigureHandle')
                    rethrow(exception)
                end
            end
            close(this.control_figure);
        end
    end
    
end