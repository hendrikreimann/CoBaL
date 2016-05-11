classdef stepEventFigure;
    properties
        main_figure;
        main_axes;
        
        left_heel_z_pos_plot;
        left_heel_z_vel_plot;
        left_heel_z_acc_plot;
        left_toes_z_pos_plot;
        left_toes_z_vel_plot;
        left_toes_z_acc_plot;
        right_heel_z_pos_plot;
        right_heel_z_vel_plot;
        right_heel_z_acc_plot;
        right_toes_z_pos_plot;
        right_toes_z_vel_plot;
        right_toes_z_acc_plot;
        
        left_touchdown_plots;
        left_pushoff_plots;
        right_touchdown_plots;
        right_pushoff_plots;
    end
    methods
        function this = stepEventFigure(varargin)
           parser = inputParser;
           default_trial_data = createTrialDataStruct;
           default_event_data = createEventDataStruct;
           default_color_struct = createColorStruct;
           default_title = 'Data Figure';

           addOptional(parser, 'trialdata', default_trial_data);
           addOptional(parser, 'eventdata', default_event_data);
           addOptional(parser, 'colors', default_color_struct);
           addOptional(parser, 'title', default_title);

           parse(parser, varargin{:});
           
           trial_data = parser.Results.trialdata;
           event_data = parser.Results.eventdata;
           color_struct = parser.Results.colors;
           
           this.main_figure = figure;
           this.main_axes = axes;
           title(parser.Results.title);
           hold on;
           
           this.left_heel_z_pos_plot = plot(trial_data.time_mocap, trial_data.left_heel_z_pos_trajectory, 'color', color_struct.left_heel_pos);
           this.left_heel_z_vel_plot = plot(trial_data.time_mocap, trial_data.left_heel_z_vel_trajectory, 'color', color_struct.left_heel_vel);
           this.left_heel_z_acc_plot = plot(trial_data.time_mocap, trial_data.left_heel_z_acc_trajectory, 'color', color_struct.left_heel_acc);
           
           this.left_toes_z_pos_plot = plot(trial_data.time_mocap, trial_data.left_toes_z_pos_trajectory, 'color', color_struct.left_toes_pos);
           this.left_toes_z_vel_plot = plot(trial_data.time_mocap, trial_data.left_toes_z_vel_trajectory, 'color', color_struct.left_toes_vel);
           this.left_toes_z_acc_plot = plot(trial_data.time_mocap, trial_data.left_toes_z_acc_trajectory, 'color', color_struct.left_toes_acc);
           
           this.right_heel_z_pos_plot = plot(trial_data.time_mocap, trial_data.right_heel_z_pos_trajectory, 'color', color_struct.right_heel_pos);
           this.right_heel_z_vel_plot = plot(trial_data.time_mocap, trial_data.right_heel_z_vel_trajectory, 'color', color_struct.right_heel_vel);
           this.right_heel_z_acc_plot = plot(trial_data.time_mocap, trial_data.right_heel_z_acc_trajectory, 'color', color_struct.right_heel_acc);
           
           this.right_toes_z_pos_plot = plot(trial_data.time_mocap, trial_data.right_toes_z_pos_trajectory, 'color', color_struct.right_toes_pos);
           this.right_toes_z_vel_plot = plot(trial_data.time_mocap, trial_data.right_toes_z_vel_trajectory, 'color', color_struct.right_toes_vel);
           this.right_toes_z_acc_plot = plot(trial_data.time_mocap, trial_data.right_toes_z_acc_trajectory, 'color', color_struct.right_toes_acc);
           
           this.left_touchdown_plots = plot(0, 0, 'v', 'linewidth', 2, 'color', color_struct.left_touchdown);
           this.left_pushoff_plots = plot([], [], 'v', 'linewidth', 2, 'color', color_struct.left_pushoff);
           this.right_touchdown_plots = plot([], [], 'v', 'linewidth', 2, 'color', color_struct.right_touchdown);
           this.right_pushoff_plots = plot([], [], 'v', 'linewidth', 2, 'color', color_struct.right_pushoff);
           
           this.updateEventPlots(event_data);
        end
        function setting_struct = getSetting(this)
            setting_struct = struct();
            setting_struct.title = this.main_axes.Title.String;
            setting_struct.position = this.main_figure.Position;
            
            setting_struct.left_heel_z_pos_plot_visibility = get(this.left_heel_z_pos_plot, 'Visible');
            setting_struct.left_heel_z_vel_plot_visibility = get(this.left_heel_z_vel_plot, 'Visible');
            setting_struct.left_heel_z_acc_plot_visibility = get(this.left_heel_z_acc_plot, 'Visible');
            setting_struct.left_toes_z_pos_plot_visibility = get(this.left_toes_z_pos_plot, 'Visible');
            setting_struct.left_toes_z_vel_plot_visibility = get(this.left_toes_z_vel_plot, 'Visible');
            setting_struct.left_toes_z_acc_plot_visibility = get(this.left_toes_z_acc_plot, 'Visible');
            setting_struct.right_heel_z_pos_plot_visibility = get(this.right_heel_z_pos_plot, 'Visible');
            setting_struct.right_heel_z_vel_plot_visibility = get(this.right_heel_z_vel_plot, 'Visible');
            setting_struct.right_heel_z_acc_plot_visibility = get(this.right_heel_z_acc_plot, 'Visible');
            setting_struct.right_toes_z_pos_plot_visibility = get(this.right_toes_z_pos_plot, 'Visible');
            setting_struct.right_toes_z_vel_plot_visibility = get(this.right_toes_z_vel_plot, 'Visible');
            setting_struct.right_toes_z_acc_plot_visibility = get(this.right_toes_z_acc_plot, 'Visible');
            
            setting_struct.left_heel_z_pos_plot_color = get(this.left_heel_z_pos_plot, 'color');
            setting_struct.left_heel_z_vel_plot_color = get(this.left_heel_z_vel_plot, 'Color');
            setting_struct.left_heel_z_acc_plot_color = get(this.left_heel_z_acc_plot, 'Color');
            setting_struct.left_toes_z_pos_plot_color = get(this.left_toes_z_pos_plot, 'Color');
            setting_struct.left_toes_z_vel_plot_color = get(this.left_toes_z_vel_plot, 'Color');
            setting_struct.left_toes_z_acc_plot_color = get(this.left_toes_z_acc_plot, 'Color');
            setting_struct.right_heel_z_pos_plot_color = get(this.right_heel_z_pos_plot, 'Color');
            setting_struct.right_heel_z_vel_plot_color = get(this.right_heel_z_vel_plot, 'Color');
            setting_struct.right_heel_z_acc_plot_color = get(this.right_heel_z_acc_plot, 'Color');
            setting_struct.right_toes_z_pos_plot_color = get(this.right_toes_z_pos_plot, 'Color');
            setting_struct.right_toes_z_vel_plot_color = get(this.right_toes_z_vel_plot, 'Color');
            setting_struct.right_toes_z_acc_plot_color = get(this.right_toes_z_acc_plot, 'Color');
            
        end
        function applySettings(this, setting_struct)
            this.main_figure.Position = setting_struct.position;
            this.main_axes.Title.String = setting_struct.title;
            
            set(this.left_heel_z_pos_plot, 'Visible', setting_struct.left_heel_z_pos_plot_visibility);
            set(this.left_heel_z_vel_plot, 'Visible', setting_struct.left_heel_z_vel_plot_visibility);
            set(this.left_heel_z_acc_plot, 'Visible', setting_struct.left_heel_z_acc_plot_visibility);
            set(this.left_toes_z_pos_plot, 'Visible', setting_struct.left_toes_z_pos_plot_visibility);
            set(this.left_toes_z_vel_plot, 'Visible', setting_struct.left_toes_z_vel_plot_visibility);
            set(this.left_toes_z_acc_plot, 'Visible', setting_struct.left_toes_z_acc_plot_visibility);
            set(this.right_heel_z_pos_plot, 'Visible', setting_struct.right_heel_z_pos_plot_visibility);
            set(this.right_heel_z_vel_plot, 'Visible', setting_struct.right_heel_z_vel_plot_visibility);
            set(this.right_heel_z_acc_plot, 'Visible', setting_struct.right_heel_z_acc_plot_visibility);
            set(this.right_toes_z_pos_plot, 'Visible', setting_struct.right_toes_z_pos_plot_visibility);
            set(this.right_toes_z_vel_plot, 'Visible', setting_struct.right_toes_z_vel_plot_visibility);
            set(this.right_toes_z_acc_plot, 'Visible', setting_struct.right_toes_z_acc_plot_visibility);
            
            set(this.left_heel_z_pos_plot, 'Color', setting_struct.left_heel_z_pos_plot_color);
            set(this.left_heel_z_vel_plot, 'Color', setting_struct.left_heel_z_vel_plot_color);
            set(this.left_heel_z_acc_plot, 'Color', setting_struct.left_heel_z_acc_plot_color);
            set(this.left_toes_z_pos_plot, 'Color', setting_struct.left_toes_z_pos_plot_color);
            set(this.left_toes_z_vel_plot, 'Color', setting_struct.left_toes_z_vel_plot_color);
            set(this.left_toes_z_acc_plot, 'Color', setting_struct.left_toes_z_acc_plot_color);
            set(this.right_heel_z_pos_plot, 'Color', setting_struct.right_heel_z_pos_plot_color);
            set(this.right_heel_z_vel_plot, 'Color', setting_struct.right_heel_z_vel_plot_color);
            set(this.right_heel_z_acc_plot, 'Color', setting_struct.right_heel_z_acc_plot_color);
            set(this.right_toes_z_pos_plot, 'Color', setting_struct.right_toes_z_pos_plot_color);
            set(this.right_toes_z_vel_plot, 'Color', setting_struct.right_toes_z_vel_plot_color);
            set(this.right_toes_z_acc_plot, 'Color', setting_struct.right_toes_z_acc_plot_color);
        end
        function updateEventPlots(this, event_data)
            time_mocap = get(this.left_heel_z_pos_plot, 'xdata');
            set(this.left_touchdown_plots, 'xdata', time_mocap(event_data.left_touchdown), 'ydata', zeros(1, length(event_data.left_touchdown)));
        end
    end
    
end