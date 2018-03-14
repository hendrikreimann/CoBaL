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

function eventGui(varargin)

    %% parse arguments
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'settings', 'eventGuiSettings.txt')
    parse(parser, varargin{:})
    settings_file = parser.Results.settings;

    condition = condition_list{1};
    trial_to_process = trial_number_list{1}(1);
    
    %% load
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
    if exist(['..' filesep settings_file], 'file')
        gui_settings_file = ['..' filesep settings_file];
    end    
    if exist(['..' filesep '..' filesep settings_file], 'file')
        gui_settings_file = ['..' filesep '..' filesep settings_file];
    end
    gui_settings = SettingsCustodian(gui_settings_file);
    
    % load data
    figure_list = getFiguresListFromSettings(gui_settings);
    variables_list = getVariableListFromSettings(gui_settings, figure_list);
    data_custodian = WalkingDataCustodian(variables_list);
    data_custodian.prepareBasicVariables(condition, trial_to_process);
    
%     trial_data = WalkingTrialData(pwd, condition, trial_to_process);
%     event_data = WalkingEventData(trial_data);
    event_data = [];
    
    %% init gui
    figure_settings_file = gui_settings.get('figure_settings_file');
    controller = eventController(data_custodian, event_data, figure_settings_file);
    
    %% stick figure and kinematic tree figure
    scene_bound = ...
      [ ...
        study_settings.get('scene_bound_x_min') study_settings.get('scene_bound_x_max'); ...
        study_settings.get('scene_bound_y_min') study_settings.get('scene_bound_y_max'); ...
        study_settings.get('scene_bound_z_min') study_settings.get('scene_bound_z_max') ...
      ];
  
    
  
  
  
%     positions = trial_data.marker_positions(1, :);
%     headers = trial_data.marker_labels;
%     if trial_data.kinematic_data_available
%         positions = [positions trial_data.joint_center_positions(1, :)];
%         headers = [headers trial_data.joint_center_labels(1, :)];
%     end
%     if trial_data.com_data_available
%         positions = [positions trial_data.com_positions(1, :)];
%         headers = [headers trial_data.com_labels(1, :)];
%     end
%     
%     % stick figure
%     controller.scene_figure = stickFigure(positions, headers, scene_bound);
%     controller.scene_figure.setColors('extended plug-in gait');
%     controller.scene_figure.addLines('extended plug-in gait');
%     
%     % kinematic tree figure
%     if trial_data.joint_angle_data_available
%         positions = trial_data.marker_positions(1, :);
%         headers = trial_data.marker_labels;
%         if ~isempty(trial_data.joint_center_positions)
%             positions = [positions trial_data.joint_center_positions(1, :)];
%             headers = [headers trial_data.joint_center_labels(1, :)];
%         end
%         [controller.kinematic_tree_stick_figure, controller.kinematic_tree_controller] = KinematicTreeController(kinematic_tree, scene_bound, 'none', positions);
%         % link perspectives of the stick figures
%         Link = linkprop([controller.kinematic_tree_stick_figure.sceneAxes controller.scene_figure.scene_axes], {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'}); 
%         setappdata(gcf, 'StoreTheLink', Link);
%     end
    
    %% trajectory figures
    % get list of figures from settings
    for i_figure = 1 : size(figure_list, 1)
        % create a figure with these settings
        new_figure = eventFigure(figure_list{i_figure, 2}, controller, data_custodian, event_data);
        
        plot_list = gui_settings.get(figure_list{i_figure, 1});
        for i_plot = 1 : size(plot_list, 1)
            if strcmp(plot_list{i_plot, 1}, 'data')
                variable_label = plot_list{i_plot, 2};
                legend_label = plot_list{i_plot, 3};
                scale_factor = str2num(plot_list{i_plot, 4});
                offset = str2num(plot_list{i_plot, 5});
                color = [str2num(plot_list{i_plot, 6}) str2num(plot_list{i_plot, 7}) str2num(plot_list{i_plot, 8})];
                linewidth = str2num(plot_list{i_plot, 9});
                new_figure.addDataPlot(variable_label, legend_label, color, scale_factor, offset, linewidth);
            end
            if strcmp(plot_list{i_plot, 1}, 'event')
                event_label = plot_list{i_plot, 2};
                variable_label = plot_list{i_plot, 3};
                legend_label = plot_list{i_plot, 4};
                marker = plot_list{i_plot, 5};
                color = [str2num(plot_list{i_plot, 6}) str2num(plot_list{i_plot, 7}) str2num(plot_list{i_plot, 8})];
                marker_size = str2num(plot_list{i_plot, 9});
                new_figure.addEventPlot(variable_label, event_label, legend_label, color, marker, marker_size);
            end
                
        end
        
        
    end
    
    %% initialize figures
    % load settings
    controller.loadFigureSettings();
%     controller.updateStretchPatches();
    
%     % select event (first left touchdown is default
%     event_label = 'left_pushoff';
% %     event_label = 'left_touchdown';
%     event_times = controller.event_data.getEventTimes(event_label);
%     if isempty(event_times)
%         event_time = 0;
%     else
%         event_time = event_times(1);
%     end
%     controller.setSelectedEvent(event_label, event_time);

end

%% nested function
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

function variable_names = getVariableListFromSettings(settings, figure_list)
    variable_names = {};
    for i_figure = 1 : length(figure_list)
        this_figure_settings = settings.get(figure_list{i_figure});
        for i_variable = 1 : size(this_figure_settings, 1);
            this_entry_type = this_figure_settings{i_variable, 1};
            this_variable_name = this_figure_settings{i_variable, 2};
            if strcmp(this_entry_type, 'data') && ~any(strcmp(variable_names, this_variable_name))
                variable_names = [variable_names; this_variable_name];
            end
        end
    end

end


