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

function stepEventGui(varargin)

    % parse arguments
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'color_scheme', 'intra')
    parse(parser, varargin{:})
    color_scheme = parser.Results.color_scheme;

    condition = condition_list{1};
    trial_to_process = trial_number_list{1}(1);
    
    
    %% load data
    trial_data = WalkingTrialData(pwd, condition, trial_to_process);
    event_data = WalkingEventData(trial_data);
    load('subjectModel.mat');
    
    % load settings
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    
    gui_settings_file = '';
    if exist(['..' filesep 'eventGuiSettings.txt'], 'file')
        gui_settings_file = ['..' filesep 'eventGuiSettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'eventGuiSettings.txt'], 'file')
        gui_settings_file = ['..' filesep '..' filesep 'eventGuiSettings.txt'];
    end
    gui_settings = SettingsCustodian(gui_settings_file);
    
    % init gui
    controller = stepEventController(trial_data, event_data);
    
    % show stick figure
    scene_bound = ...
      [ ...
        study_settings.get('scene_bound_x_min') study_settings.get('scene_bound_x_max'); ...
        study_settings.get('scene_bound_y_min') study_settings.get('scene_bound_y_max'); ...
        study_settings.get('scene_bound_z_min') study_settings.get('scene_bound_z_max') ...
      ];
    positions = trial_data.marker_positions(1, :);
    headers = trial_data.marker_labels;
    if trial_data.kinematic_data_available
        positions = [positions trial_data.joint_center_positions(1, :)];
        headers = [headers trial_data.joint_center_labels(1, :)];
    end
    if trial_data.com_data_available
        positions = [positions trial_data.com_positions(1, :)];
        headers = [headers trial_data.com_labels(1, :)];
    end
    
    controller.scene_figure = stickFigure(positions, headers, scene_bound);
    controller.scene_figure.setColors('extended plug-in gait');
    controller.scene_figure.addLines('extended plug-in gait');

    if trial_data.joint_angle_data_available
        positions = trial_data.marker_positions(1, :);
        headers = trial_data.marker_labels;
        if ~isempty(trial_data.joint_center_positions)
            positions = [positions trial_data.joint_center_positions(1, :)];
            headers = [headers trial_data.joint_center_labels(1, :)];
        end
        controller.kinematic_tree_controller = KinematicTreeController(kinematic_tree, scene_bound, 'none', positions);
        Link = linkprop([controller.kinematic_tree_controller.sceneAxes controller.scene_figure.scene_axes], {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'}); 
        setappdata(gcf, 'StoreTheLink', Link);
    end
    
    % link perspectives of the stick figures
%     Link = linkprop([controller.kinematic_tree_controller.sceneAxes controller.scene_figure.scene_axes], {'CameraUpVector', 'CameraPosition', 'CameraTarget'}); 
    
    % get list of figures from settings
    figure_list = getFiguresListFromSettings(gui_settings);
    for i_figure = 1 : size(figure_list, 1)
        % create a figure with these settings
        new_figure = stepEventFigure(figure_list{i_figure, 2}, controller, trial_data, event_data);
        
        plot_list = gui_settings.get(figure_list{i_figure, 1});
        for i_plot = 1 : size(plot_list, 1)
            if strcmp(plot_list{i_plot, 1}, 'data')
                variable_label = plot_list{i_plot, 2};
                scale_factor = str2num(plot_list{i_plot, 3});
                offset = str2num(plot_list{i_plot, 4});
                color = [str2num(plot_list{i_plot, 5}) str2num(plot_list{i_plot, 6}) str2num(plot_list{i_plot, 7})];
                new_figure.addDataPlot(variable_label, color, scale_factor, offset);
            end
            if strcmp(plot_list{i_plot, 1}, 'event')
                event_label = plot_list{i_plot, 2};
                variable_label = plot_list{i_plot, 3};
                marker = plot_list{i_plot, 4};
                color = [str2num(plot_list{i_plot, 5}) str2num(plot_list{i_plot, 6}) str2num(plot_list{i_plot, 7})];
                new_figure.addEventPlot(variable_label, event_label, color, marker);
            end
                
        end
        
        
    end
    
    % load settings
    controller.loadFigureSettings();
    controller.updateStretchPatches();
    
    % select event (first left touchdown is default

%     event_label = 'left_touchdown';
    event_label = 'left_pushoff';
    event_times = controller.event_data.getEventTimes(event_label);
    if isempty(event_times)
        event_time = 0;
    else
        event_time = event_times(1);
    end
    controller.setSelectedEvent(event_label, event_time);









end

function figures_list = getFiguresListFromSettings(settings)
    field_names = settings.getAllSettingsNames;
    figures_list = {};
    for i_field = 1 : length(field_names)
        this_field_name = field_names{i_field};
        if length(this_field_name) > 6 && strcmp(this_field_name(end-5:end), 'Figure')
            figure_label = strrep(this_field_name(1:end-6), '_', ' ');
            figures_list = [figures_list; {this_field_name, figure_label}]; %#ok<AGROW>
        end
    end
    
end

