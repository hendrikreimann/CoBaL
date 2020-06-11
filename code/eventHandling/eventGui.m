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
    addParameter(parser, 'settings', 'eventGui')
    parse(parser, varargin{:})
    settings_file = parser.Results.settings;

    condition = condition_list{1};
    trial_to_process = trial_number_list{1}(1);
    
    %% load
    gui_settings = loadSettingsFromFile(settings_file); % TODO: only using default, not specified settings file
    
    % load data
    if gui_settings.isfield('figures')
        % settings are specified in new way, so use that
        figure_settings = getFigureSettings(gui_settings);
    else
        % settings are specified in old way
        warning('Figure settings are in old format, please update.')
        figure_list = getFiguresListFromSettings(gui_settings);
        
        % transform into new format
        figure_settings = cell(length(figure_list), 1);
        for i_figure = 1 : length(figure_list)
            this_figure_settings = struct;
            
            % store data
            this_figure_settings.data = gui_settings.get(figure_list{i_figure, 1});
            % store individual entries in info table
            this_figure_settings.title = figure_list{i_figure, 2};
            this_figure_settings.origin_horizontal = 25;
            this_figure_settings.origin_vertical = 25;
            this_figure_settings.width = 50;
            this_figure_settings.height = 50;

            figure_settings{i_figure} = this_figure_settings; %#ok<AGROW>
        end
        
    end
    variables_list = getVariableListFromSettings(figure_settings);
    
    data_custodian = WalkingDataCustodian(variables_list);
    data_custodian.prepareBasicVariables(condition, trial_to_process);
    
    event_data = eventData(data_custodian);
    
    %% init gui
%     figure_settings_file = gui_settings.get('figure_settings_file');
    controller = eventController(data_custodian, event_data);
    
    %% stick figure and kinematic tree figure
  
    if gui_settings.get('show_simple_stick_figure', 1)
        scene_bound = ...
          [ ...
            gui_settings.get('scene_bound_x_min') gui_settings.get('scene_bound_x_max'); ...
            gui_settings.get('scene_bound_y_min') gui_settings.get('scene_bound_y_max'); ...
            gui_settings.get('scene_bound_z_min') gui_settings.get('scene_bound_z_max') ...
          ];
      
        marker_trajectories = data_custodian.getBasicVariableData('marker_trajectories');
        joint_center_trajectories = data_custodian.getBasicVariableData('joint_center_trajectories');
        com_trajectories = data_custodian.getBasicVariableData('com_trajectories');
        % TODO: deal with cases where we have only marker data and no kinematic data yet

        marker_labels = data_custodian.basic_variable_labels.marker_trajectories;
        joint_center_labels = data_custodian.basic_variable_labels.joint_center_trajectories;
        com_labels = data_custodian.basic_variable_labels.com_trajectories;

        positions = [marker_trajectories(1, :), joint_center_trajectories(1, :), com_trajectories(1, :)];
        labels = [marker_labels joint_center_labels com_labels];

        % stick figure
        controller.scene_figure = stickFigure(positions, labels, scene_bound);
        controller.scene_figure.setColors('extended plug-in gait');
        controller.scene_figure.addLines('extended plug-in gait');
    end    

    
    %% trajectory figures
    % get list of figures from settings
    for i_figure = 1 : size(figure_settings, 1)
        % create a figure with these settings
        this_figure_settings = figure_settings{i_figure};
        new_figure = eventFigure(this_figure_settings, controller, data_custodian, event_data);
        
        % TODO: this adds the data to the figure, should be moved into the
        % constructor of the figure
        plot_list = this_figure_settings.data;
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
    controller.updateStretchPatches();
    
    % select event (first left touchdown is default
    if ~isempty(controller.event_data.event_labels)
        event_label = controller.event_data.event_labels{1};
    else
        event_label = '';
    end
    event_times = controller.event_data.getEventTimes(event_label);
    if isempty(event_times)
        event_time = 0;
    else
        event_time = event_times(1);
    end
%     controller.setSelectedEvent(event_label, event_time);
    controller.setSelectedEvent([], []);

end

%% nested functions

function figure_settings = getFigureSettings(gui_settings)
    figure_header = gui_settings.get('figures_header');
    figures_table = gui_settings.get('figures');
    
    number_of_figures = size(figures_table, 1);
    figure_settings = cell(number_of_figures, 1);
    for i_figure = 1 : number_of_figures
        this_figure_settings = struct;
        % store data
        this_figure_settings.data = gui_settings.get(figures_table{i_figure, strcmp(figure_header, 'data')});
        % store individual entries in info table
        info = gui_settings.get(figures_table{i_figure, strcmp(figure_header, 'info')});
        for i_entry = 1 : size(info, 1)
            this_figure_settings.(info{i_entry, 1}) = info{i_entry, 2};
        end
        % go through expected entries and enter defaults if needed
        if ~isfield(this_figure_settings, 'title')
            this_figure_settings.title = 'untitled figure';
        end
        if ~isfield(this_figure_settings, 'origin_horizontal')
            this_figure_settings.origin_horizontal = 25;
        end
        if ~isfield(this_figure_settings, 'origin_vertical')
            this_figure_settings.origin_vertical = 25;
        end
        if ~isfield(this_figure_settings, 'width')
            this_figure_settings.width = 50;
        end
        if ~isfield(this_figure_settings, 'height')
            this_figure_settings.height = 50;
        end
        
        % transform entries to numbers if needed
        if ~isnumeric(this_figure_settings.origin_horizontal)
            this_figure_settings.origin_horizontal = str2double(this_figure_settings.origin_horizontal);
        end
        if ~isnumeric(this_figure_settings.origin_vertical)
            this_figure_settings.origin_vertical = str2double(this_figure_settings.origin_vertical);
        end
        if ~isnumeric(this_figure_settings.width)
            this_figure_settings.width = str2double(this_figure_settings.width);
        end
        if ~isnumeric(this_figure_settings.height)
            this_figure_settings.height = str2double(this_figure_settings.height);
        end
        
        figure_settings{i_figure} = this_figure_settings;
    end
end

function variable_names = getVariableListFromSettings(figure_settings)
    % get variables for trajectory figures
    variable_names = {};
    for i_figure = 1 : size(figure_settings, 1)
        this_figure_settings = figure_settings{i_figure};
        this_figure_data_table = this_figure_settings.data;
        for i_variable = 1 : size(this_figure_data_table, 1);
            this_entry_type = this_figure_data_table{i_variable, 1};
            this_variable_name = this_figure_data_table{i_variable, 2};
            if strcmp(this_entry_type, 'data') && ~any(strcmp(variable_names, this_variable_name))
                variable_names = [variable_names; this_variable_name];
            end
        end
    end
    

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



