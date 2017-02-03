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
    
    %% set plot stuff
    if strcmp(color_scheme, 'side')
        color_left_heel = [0.5 0.8 0];
        color_left_toes = [0 0.8 0.5];
        color_left_fz = [0.2 0.8 0.2];
        color_left_touchdown = [0.0 0.5 0.8];
        color_left_pushoff = [0.0 0.5 0.8];
        marker_left_touchdown = 'v';
        marker_left_pushoff = '^';

        color_right_heel = [0.8 0.5 0];
        color_right_toes = [0.8 0 0.5];
        color_right_fz = [0.8 0.2 0.2];
        color_right_touchdown = [0.5 0 0.8];
        color_right_pushoff = [0.5 0 0.8];
        marker_right_touchdown = 'v';
        marker_right_pushoff = '^';

        scale_factor_heel = 6;
        scale_factor_toes = 20;
        scale_factor_fz = 1/600;
        offset_heel = - 0.07;
        offset_toes = - 0.05;
    elseif strcmp(color_scheme, 'intra')
        color_left_heel = [232 26 75]*1/255;
        color_left_toes = [57 181 74]*1/255;
        color_left_fz = [134 50 140]*1/255;
        color_left_touchdown = [0 131 202]*1/255;
        color_left_pushoff = [241 90 34]*1/255;

        color_right_heel = [232 26 75]*1/255;
        color_right_toes = [57 181 74]*1/255;
        color_right_fz = [134 50 140]*1/255;
        color_right_touchdown = [0 131 202]*1/255;
        color_right_pushoff = [241 90 34]*1/255;
        
        color_stimulus_state = [255 198 11]*1/255;
        
        marker_right_touchdown = 'v';
        marker_right_pushoff = '^';
        marker_left_touchdown = 'v';
        marker_left_pushoff = '^';
        
        scale_factor_heel = 1;
        scale_factor_toes = 1;
        scale_factor_fz = 1/2000;
        scale_factor_stimulus_state = 1/8;
        offset_heel = - 0.03;
        offset_toes = - 0.0;
    end
    
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
    study_settings = loadSettingsFile(study_settings_file);
    
    % init gui
    controller = stepEventController(trial_data, event_data);
    
    % show stick figure
    scene_bound = ...
      [ ...
        study_settings.scene_bound_x_min study_settings.scene_bound_x_max; ...
        study_settings.scene_bound_y_min study_settings.scene_bound_y_max; ...
        study_settings.scene_bound_z_min study_settings.scene_bound_z_max ...
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
    
    % create position figures
    left_foot_pos_figure = stepEventFigure('Positions Left', controller, trial_data, event_data);
    left_foot_pos_figure.addDataPlot('left_heel_z_pos', color_left_heel, scale_factor_heel, offset_heel);
    left_foot_pos_figure.addDataPlot('left_toes_z_pos', color_left_toes, scale_factor_toes, offset_toes);
    
%     step_event_figure.addDataPlot('left_foot_fz', color_left_fz, scale_factor_fz);
%     left_foot_pos_figure.addDataPlot('right_fz', color_right_fz, scale_factor_fz);
    left_foot_pos_figure.addEventPlot('left_heel_z_pos', 'left_touchdown', color_left_touchdown, marker_left_touchdown);
    left_foot_pos_figure.addEventPlot('left_toes_z_pos', 'left_pushoff', color_left_pushoff, marker_left_pushoff);
%     step_event_figure.addEventPlot('left_fz', 'left_pushoff', color_left_pushoff, marker_left_pushoff);
%     step_event_figure.addEventPlot('left_fz', 'left_touchdown', color_left_touchdown, marker_left_touchdown);
%     step_event_figure.addDataPlot('stimulus_state', color_stimulus_state, scale_factor_stimulus_state);
%     step_event_figure.addEventPlot('left_fz', 'left_pushoff', color_left_pushoff, marker_left_pushoff);
%     step_event_figure.addEventPlot('left_fz', 'left_touchdown', color_left_touchdown, marker_left_touchdown);
    
    right_foot_pos_figure = stepEventFigure('Positions Right', controller, trial_data, event_data);
    right_foot_pos_figure.addDataPlot('right_heel_z_pos', color_right_heel, scale_factor_heel, offset_heel);
    right_foot_pos_figure.addDataPlot('right_toes_z_pos', color_right_toes, scale_factor_toes, offset_toes);
%     step_event_figure.addDataPlot('right_fz', color_right_fz, scale_factor_fz);
%     right_foot_pos_figure.addDataPlot('left_fz', color_left_fz, scale_factor_fz);
    right_foot_pos_figure.addEventPlot('right_heel_z_pos', 'right_touchdown', color_right_touchdown, marker_right_touchdown);
    right_foot_pos_figure.addEventPlot('right_toes_z_pos', 'right_pushoff', color_right_pushoff, marker_right_pushoff);
%     step_event_figure.addEventPlot('right_fz', 'right_pushoff', color_right_pushoff, marker_right_pushoff);
%     step_event_figure.addEventPlot('right_fz', 'right_touchdown', color_right_touchdown, marker_right_touchdown);
%     step_event_figure.addDataPlot('stimulus_state', color_stimulus_state, scale_factor_stimulus_state);

    % create velocity figures
    left_foot_vel_figure = stepEventFigure('Velocities Left', controller, trial_data, event_data);
%     step_event_figure.addDataPlot('left_fz', color_left_fz, scale_factor_fz*5);
    
%     step_event_figure.addEventPlot('left_heel_y_vel', 'left_touchdown', color_left_touchdown, marker_left_touchdown);
%     step_event_figure.addEventPlot('left_heel_z_vel', 'left_pushoff', color_left_pushoff, marker_left_pushoff);
%     step_event_figure.addEventPlot('right_heel_y_vel', 'right_touchdown', color_right_touchdown, marker_right_touchdown);
%     step_event_figure.addEventPlot('right_heel_z_vel', 'right_pushoff', color_right_pushoff, marker_right_pushoff);
    
    left_foot_vel_figure.addDataPlot('left_toes_z_vel', color_left_toes);
    left_foot_vel_figure.addEventPlot('left_toes_z_vel', 'left_pushoff', color_left_pushoff, marker_left_pushoff);
%     step_event_figure.addEventPlot('left_toes_z_vel', 'left_touchdown', color_left_touchdown, marker_left_touchdown);
    
    right_foot_vel_figure = stepEventFigure('Velocities Right', controller, trial_data, event_data);
%     step_event_figure.addDataPlot('right_heel_z_vel', color_right_heel);
    right_foot_vel_figure.addDataPlot('right_toes_z_vel', color_right_toes);
    right_foot_vel_figure.addEventPlot('right_toes_z_vel', 'right_pushoff', color_right_pushoff, marker_right_pushoff);
%     step_event_figure.addEventPlot('right_toes_z_vel', 'right_touchdown', color_right_touchdown, marker_right_touchdown);
    
    % create acceleration figure
    left_foot_acc_figure = stepEventFigure('Acceleration Left', controller, trial_data, event_data);
    left_foot_acc_figure.addDataPlot('left_heel_z_acc', color_left_heel);
    left_foot_acc_figure.addEventPlot('left_heel_z_acc', 'left_touchdown', color_left_touchdown, marker_left_touchdown);
%     step_event_figure.addDataPlot('left_toes_z_acc', color_left_toes);
%     step_event_figure.addEventPlot('left_toes_z_acc', 'left_pushoff', color_left_pushoff, marker_left_pushoff);
    
    right_foot_acc_figure = stepEventFigure('Acceleration Right', controller, trial_data, event_data);
    right_foot_acc_figure.addDataPlot('right_heel_z_acc', color_right_heel);
    right_foot_acc_figure.addEventPlot('right_heel_z_acc', 'right_touchdown', color_right_touchdown, marker_right_touchdown)
%     step_event_figure.addDataPlot('right_toes_z_acc', color_right_toes);
%     step_event_figure.addEventPlot('right_toes_z_acc', 'right_pushoff', color_right_pushoff, marker_right_pushoff);
    
    step_event_axes = zeros(size(controller.figureSelectionBox.String));
    for i_figure = 1 : length((controller.figureSelectionBox.String))
        step_event_axes(i_figure) = controller.figureSelectionBox.UserData{i_figure}.main_axes;
    end
    linkaxes(step_event_axes, 'x');
    
    % load settings
    controller.loadFigureSettings();
    controller.updateStretchPatches();
    
%     controller.findEvents();

    % select event (first left touchdown is default

    event_label = 'left_touchdown';
    event_label = 'left_pushoff';
    event_times = controller.event_data.getEventTimes(event_label);
    event_time = event_times(1);
    controller.setSelectedEvent(event_label, event_time);









end