function stepEventGui(dataDirectory, trialToProcess)
    if nargin < 2
        trialToProcess = 1;
    end
    if nargin < 1
        dataDirectory = pwd;
    end
    
    % load data    
    trial_data = loadTrialData;
    
    % load or create event data
    event_data = createEventDataStruct();
    
    % init gui
    color_struct = createColorStruct();

    controller = createController();
    step_event_figures{1} = createStepEventFigure('Positions');
    setCurrentFigureToDefault();
    step_event_figures{2} = createStepEventFigure('Derivatives');
    setCurrentFigureToDefault();
    
%     loadFigureSettings();



    %% GUI
    function controller = createController()
        controller = struct();

        figure_height = 600;
        figure_width = 420;
        controller.control_figure = figure('position', [1600 300 figure_width figure_height], 'Units', 'pixels');
        
        % figure control
        figure_panel_height = 255;
        figure_panel = uipanel(controller.control_figure, 'Title', 'Figure Control', 'FontSize', 12, 'BackgroundColor', 'white', 'Units', 'pixels', 'Position', [5, figure_height-figure_panel_height-5, figure_width-10, figure_panel_height]);
        
    	controller.left_heel_pos_box = uicontrol(figure_panel, 'Style', 'checkbox', 'Value', 0, 'Callback', @updatePlotVisibility, 'Position', [5, figure_panel_height - 40, 20, 20], 'Value', 1, 'BackgroundColor', 'white');
        uicontrol(figure_panel, 'Style', 'text', 'string', 'left heel pos', 'Position', [30, figure_panel_height-43, 150, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
    	controller.left_heel_vel_box = uicontrol(figure_panel, 'Style', 'checkbox', 'Value', 0, 'Callback', @updatePlotVisibility, 'Position', [5, figure_panel_height - 60, 20, 20], 'Value', 1, 'BackgroundColor', 'white');
        uicontrol(figure_panel, 'Style', 'text', 'string', 'left heel vel', 'Position', [30, figure_panel_height-63, 150, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
    	controller.left_heel_acc_box = uicontrol(figure_panel, 'Style', 'checkbox', 'Value', 0, 'Callback', @updatePlotVisibility, 'Position', [5, figure_panel_height - 80, 20, 20], 'Value', 1, 'BackgroundColor', 'white');
        uicontrol(figure_panel, 'Style', 'text', 'string', 'left heel acc', 'Position', [30, figure_panel_height-83, 150, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');

    	controller.left_toes_pos_box = uicontrol(figure_panel, 'Style', 'checkbox', 'Value', 0, 'Callback', @updatePlotVisibility, 'Position', [105, figure_panel_height - 40, 20, 20], 'Value', 1, 'BackgroundColor', 'white');
        uicontrol(figure_panel, 'Style', 'text', 'string', 'left toes pos', 'Position', [130, figure_panel_height-43, 150, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
    	controller.left_toes_vel_box = uicontrol(figure_panel, 'Style', 'checkbox', 'Value', 0, 'Callback', @updatePlotVisibility, 'Position', [105, figure_panel_height - 60, 20, 20], 'Value', 1, 'BackgroundColor', 'white');
        uicontrol(figure_panel, 'Style', 'text', 'string', 'left toes vel', 'Position', [130, figure_panel_height-63, 150, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
    	controller.left_toes_acc_box = uicontrol(figure_panel, 'Style', 'checkbox', 'Value', 0, 'Callback', @updatePlotVisibility, 'Position', [105, figure_panel_height - 80, 20, 20], 'Value', 1, 'BackgroundColor', 'white');
        uicontrol(figure_panel, 'Style', 'text', 'string', 'left toes acc', 'Position', [130, figure_panel_height-83, 150, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');

    	controller.right_heel_pos_box = uicontrol(figure_panel, 'Style', 'checkbox', 'Value', 0, 'Callback', @updatePlotVisibility, 'Position', [205, figure_panel_height - 40, 20, 20], 'Value', 1, 'BackgroundColor', 'white');
        uicontrol(figure_panel, 'Style', 'text', 'string', 'right heel pos', 'Position', [230, figure_panel_height-43, 150, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
    	controller.right_heel_vel_box = uicontrol(figure_panel, 'Style', 'checkbox', 'Value', 0, 'Callback', @updatePlotVisibility, 'Position', [205, figure_panel_height - 60, 20, 20], 'Value', 1, 'BackgroundColor', 'white');
        uicontrol(figure_panel, 'Style', 'text', 'string', 'right heel vel', 'Position', [230, figure_panel_height-63, 150, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
    	controller.right_heel_acc_box = uicontrol(figure_panel, 'Style', 'checkbox', 'Value', 0, 'Callback', @updatePlotVisibility, 'Position', [205, figure_panel_height - 80, 20, 20], 'Value', 1, 'BackgroundColor', 'white');
        uicontrol(figure_panel, 'Style', 'text', 'string', 'right heel acc', 'Position', [230, figure_panel_height-83, 150, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');

    	controller.right_toes_pos_box = uicontrol(figure_panel, 'Style', 'checkbox', 'Value', 0, 'Callback', @updatePlotVisibility, 'Position', [305, figure_panel_height - 40, 20, 20], 'Value', 1, 'BackgroundColor', 'white');
        uicontrol(figure_panel, 'Style', 'text', 'string', 'right toes pos', 'Position', [330, figure_panel_height-43, 150, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
    	controller.right_toes_vel_box = uicontrol(figure_panel, 'Style', 'checkbox', 'Value', 0, 'Callback', @updatePlotVisibility, 'Position', [305, figure_panel_height - 60, 20, 20], 'Value', 1, 'BackgroundColor', 'white');
        uicontrol(figure_panel, 'Style', 'text', 'string', 'right toes vel', 'Position', [330, figure_panel_height-63, 150, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
    	controller.right_toes_acc_box = uicontrol(figure_panel, 'Style', 'checkbox', 'Value', 0, 'Callback', @updatePlotVisibility, 'Position', [305, figure_panel_height - 80, 20, 20], 'Value', 1, 'BackgroundColor', 'white');
        uicontrol(figure_panel, 'Style', 'text', 'string', 'right toes acc', 'Position', [330, figure_panel_height-83, 150, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');

        controller.figureSelectionBox = uicontrol(figure_panel, 'Style', 'listbox', 'Position', [5, figure_panel_height-185, 395, 100], 'Fontsize', 12, 'HorizontalAlignment', 'left', 'callback', @updateVisibilityCheckBoxes);
        uicontrol(figure_panel, 'Style', 'pushbutton', 'Position', [5, figure_panel_height-250, 130, 60], 'Fontsize', 12, 'String', 'Save Figure Settings', 'callback', @saveFigureSettings);
        uicontrol(figure_panel, 'Style', 'pushbutton', 'Position', [140, figure_panel_height-250, 130, 60], 'Fontsize', 12, 'String', 'Load Figure Settings', 'callback', @loadFigureSettings);
        
        % event controls
        events_panel_height = 255;
        events_panel = uipanel(controller.control_figure, 'Title', 'Events Control', 'FontSize', 12, 'BackgroundColor', 'white', 'Units', 'pixels', 'Position', [5, figure_height-figure_panel_height-events_panel_height-5, figure_width-10, events_panel_height]);
        uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [5, events_panel_height-75, 130, 60], 'Fontsize', 12, 'String', 'Find Events', 'callback', @findEvents);
        uicontrol(events_panel, 'Style', 'pushbutton', 'Position', [140, events_panel_height-75, 130, 60], 'Fontsize', 12, 'String', 'Save Events', 'callback', @saveEvents);
        
        uicontrol(events_panel, 'Style', 'text', 'string', 'Heel peak mrominence (m):', 'Position', [5, events_panel_height-100, 180, 20], 'Fontsize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
        controller.heel_pos_peak_width = uicontrol(events_panel, 'Style', 'edit', 'BackgroundColor', 'white', 'Position', [150, events_panel_height-97, 40, 20], 'String', '0.05');
    end
    function step_event_figure = createStepEventFigure(title)
        step_event_figure = stepEventFigure ...
          ( ...
            'trialdata', trial_data, ...
            'colors', color_struct, ...
            'eventdata', event_data, ...
            'title', title ...
          );
        
        controller.figureSelectionBox.String = [controller.figureSelectionBox.String; {title}];
    end
    function setCurrentFigureToDefault()
        controller.left_heel_pos_box.Value = 1;
        controller.left_heel_vel_box.Value = 0;
        controller.left_heel_acc_box.Value = 0;

        controller.left_toes_pos_box.Value = 1;
        controller.left_toes_vel_box.Value = 0;
        controller.left_toes_acc_box.Value = 0;
        
        controller.right_heel_pos_box.Value = 1;
        controller.right_heel_vel_box.Value = 0;
        controller.right_heel_acc_box.Value = 0;
        
        controller.right_toes_pos_box.Value = 1;
        controller.right_toes_vel_box.Value = 0;
        controller.right_toes_acc_box.Value = 0;
        
        updatePlotVisibility();
    end
    function updatePlotVisibility(varargin)
        figure_index = controller.figureSelectionBox.Value;

        % left heel
        if controller.left_heel_pos_box.Value
            set(step_event_figures{figure_index}.left_heel_z_pos_plot, 'visible', 'on')
        else
            set(step_event_figures{figure_index}.left_heel_z_pos_plot, 'visible', 'off')
        end
        if controller.left_heel_vel_box.Value
            set(step_event_figures{figure_index}.left_heel_z_vel_plot, 'visible', 'on')
        else
            set(step_event_figures{figure_index}.left_heel_z_vel_plot, 'visible', 'off')
        end
        if controller.left_heel_acc_box.Value
            set(step_event_figures{figure_index}.left_heel_z_acc_plot, 'visible', 'on')
        else
            set(step_event_figures{figure_index}.left_heel_z_acc_plot, 'visible', 'off')
        end

        % left toes
        if controller.left_toes_pos_box.Value
            set(step_event_figures{figure_index}.left_toes_z_pos_plot, 'visible', 'on')
        else
            set(step_event_figures{figure_index}.left_toes_z_pos_plot, 'visible', 'off')
        end
        if controller.left_toes_vel_box.Value
            set(step_event_figures{figure_index}.left_toes_z_vel_plot, 'visible', 'on')
        else
            set(step_event_figures{figure_index}.left_toes_z_vel_plot, 'visible', 'off')
        end
        if controller.left_toes_acc_box.Value
            set(step_event_figures{figure_index}.left_toes_z_acc_plot, 'visible', 'on')
        else
            set(step_event_figures{figure_index}.left_toes_z_acc_plot, 'visible', 'off')
        end

        % right heel
        if controller.right_heel_pos_box.Value
            set(step_event_figures{figure_index}.right_heel_z_pos_plot, 'visible', 'on')
        else
            set(step_event_figures{figure_index}.right_heel_z_pos_plot, 'visible', 'off')
        end
        if controller.right_heel_vel_box.Value
            set(step_event_figures{figure_index}.right_heel_z_vel_plot, 'visible', 'on')
        else
            set(step_event_figures{figure_index}.right_heel_z_vel_plot, 'visible', 'off')
        end
        if controller.right_heel_acc_box.Value
            set(step_event_figures{figure_index}.right_heel_z_acc_plot, 'visible', 'on')
        else
            set(step_event_figures{figure_index}.right_heel_z_acc_plot, 'visible', 'off')
        end

        % right toes
        if controller.right_toes_pos_box.Value
            set(step_event_figures{figure_index}.right_toes_z_pos_plot, 'visible', 'on')
        else
            set(step_event_figures{figure_index}.right_toes_z_pos_plot, 'visible', 'off')
        end
        if controller.right_toes_vel_box.Value
            set(step_event_figures{figure_index}.right_toes_z_vel_plot, 'visible', 'on')
        else
            set(step_event_figures{figure_index}.right_toes_z_vel_plot, 'visible', 'off')
        end
        if controller.right_toes_acc_box.Value
            set(step_event_figures{figure_index}.right_toes_z_acc_plot, 'visible', 'on')
        else
            set(step_event_figures{figure_index}.right_toes_z_acc_plot, 'visible', 'off')
        end
    end
    function updateVisibilityCheckBoxes(varargin)
        figure_index = controller.figureSelectionBox.Value;

        % left heel
        if strcmp(get(step_event_figures{figure_index}.left_heel_z_pos_plot, 'visible'), 'on')
            controller.left_heel_pos_box.Value = 1;
        else
            controller.left_heel_pos_box.Value = 0;
        end
        if strcmp(get(step_event_figures{figure_index}.left_heel_z_vel_plot, 'visible'), 'on')
            controller.left_heel_vel_box.Value = 1;
        else
            controller.left_heel_vel_box.Value = 0;
        end
        if strcmp(get(step_event_figures{figure_index}.left_heel_z_acc_plot, 'visible'), 'on')
            controller.left_heel_acc_box.Value = 1;
        else
            controller.left_heel_acc_box.Value = 0;
        end

        % left toes
        if strcmp(get(step_event_figures{figure_index}.left_toes_z_pos_plot, 'visible'), 'on')
            controller.left_toes_pos_box.Value = 1;
        else
            controller.left_toes_pos_box.Value = 0;
        end
        if strcmp(get(step_event_figures{figure_index}.left_toes_z_vel_plot, 'visible'), 'on')
            controller.left_toes_vel_box.Value = 1;
        else
            controller.left_toes_vel_box.Value = 0;
        end
        if strcmp(get(step_event_figures{figure_index}.left_toes_z_acc_plot, 'visible'), 'on')
            controller.left_toes_acc_box.Value = 1;
        else
            controller.left_toes_acc_box.Value = 0;
        end

        % right heel
        if strcmp(get(step_event_figures{figure_index}.right_heel_z_pos_plot, 'visible'), 'on')
            controller.right_heel_pos_box.Value = 1;
        else
            controller.right_heel_pos_box.Value = 0;
        end
        if strcmp(get(step_event_figures{figure_index}.right_heel_z_vel_plot, 'visible'), 'on')
            controller.right_heel_vel_box.Value = 1;
        else
            controller.right_heel_vel_box.Value = 0;
        end
        if strcmp(get(step_event_figures{figure_index}.right_heel_z_acc_plot, 'visible'), 'on')
            controller.right_heel_acc_box.Value = 1;
        else
            controller.right_heel_acc_box.Value = 0;
        end

        % right toes
        if strcmp(get(step_event_figures{figure_index}.right_toes_z_pos_plot, 'visible'), 'on')
            controller.right_toes_pos_box.Value = 1;
        else
            controller.right_toes_pos_box.Value = 0;
        end
        if strcmp(get(step_event_figures{figure_index}.right_toes_z_vel_plot, 'visible'), 'on')
            controller.right_toes_vel_box.Value = 1;
        else
            controller.right_toes_vel_box.Value = 0;
        end
        if strcmp(get(step_event_figures{figure_index}.right_toes_z_acc_plot, 'visible'), 'on')
            controller.right_toes_acc_box.Value = 1;
        else
            controller.right_toes_acc_box.Value = 0;
        end
    end
    function saveFigureSettings(varargin)
        figure_settings = cell(1, length(controller.figureSelectionBox));
        for i_figure = 1 : length(controller.figureSelectionBox.String)
            figure_settings{i_figure} = step_event_figures{i_figure}.getSetting();
        end
        control_figure_setting = struct();
        control_figure_setting.position = controller.control_figure.Position;
        
        save_file = '/Users/reimajbi/Library/Application Support/stepEventFinder/figureSettings.mat';
        save(save_file, 'figure_settings', 'control_figure_setting');
    end
    function loadFigureSettings(varargin)
%         load_file = '/Users/reimajbi/Dropbox/CoBaL/figureSettings.mat';
        load_file = '/Users/reimajbi/Library/Application Support/stepEventFinder/figureSettings.mat';
        
        load(load_file, 'figure_settings', 'control_figure_setting')
        for i_figure = 1 : length(figure_settings)
            if length(step_event_figures) < i_figure
                step_event_figures{i_figure} = createStepEventFigure();
            end
            step_event_figures{i_figure}.applySettings(figure_settings{i_figure});
        end
        controller.control_figure.Position = control_figure_setting.position;
    end
    function updateEventPlots()
        for i_figure = 1 : length(step_event_figures)
            step_event_figures{i_figure}.updateEventPlots(event_data);
        end
    end

    %% events
    function findEvents(varargin)
        min_peak_prominence = str2num(controller.heel_pos_peak_width.String);

        % left
        [~, left_touchdown_indices_mocap] = findpeaks(-trial_data.left_heel_z_pos_trajectory, 'MinPeakProminence', min_peak_prominence);

        % right
        [~, right_touchdown_indices_mocap] = findpeaks(-trial_data.right_heel_z_pos_trajectory, 'MinPeakProminence', min_peak_prominence);

        event_data.left_touchdown = left_touchdown_indices_mocap;
        event_data.right_touchdown = right_touchdown_indices_mocap;

        updateEventPlots();


%         if strcmp(method_pushoff, 'first_velocity_peak')
%             % for pushoff, find the first significant toes z-velocity peak after each heelstrike
%             min_peak_prominence = 0.4;
%             [~, left_toes_vel_peak_locations] = findpeaks(left_toes_marker_z_vel_trajectory, 'MinPeakProminence', min_peak_prominence);
%             left_pushoff_indices_mocap = zeros(size(left_touchdown_indices_mocap));
%             for i_touchdown = 1 : length(left_touchdown_indices_mocap)
%                 pushoff_index_index = find(left_toes_vel_peak_locations > left_touchdown_indices_mocap(i_touchdown), 1, 'first');
%                 if ~isempty(pushoff_index_index)
%                     left_pushoff_indices_mocap(i_touchdown) = left_toes_vel_peak_locations(pushoff_index_index);
%                 end
%             end
%             left_pushoff_indices_mocap(left_pushoff_indices_mocap==0) = [];
% 
%             min_peak_prominence = 0.4;
%             [~, right_toes_vel_peak_locations] = findpeaks(right_toes_marker_z_vel_trajectory, 'MinPeakProminence', min_peak_prominence);
%             right_pushoff_indices_mocap = zeros(size(right_touchdown_indices_mocap));
%             for i_touchdown = 1 : length(right_touchdown_indices_mocap)
%                 pushoff_index_index = find(right_toes_vel_peak_locations > right_touchdown_indices_mocap(i_touchdown), 1, 'first');
%                 if ~isempty(pushoff_index_index)
%                     right_pushoff_indices_mocap(i_touchdown) = right_toes_vel_peak_locations(pushoff_index_index);
%                 end
%             end
%             right_pushoff_indices_mocap(right_pushoff_indices_mocap==0) = [];
%         end
        
    end
    function saveEvents(varargin)
        
    end

    %% file system interaction
    function data = loadTrialData
        
        loaded_subject_info = load([dataDirectory filesep 'subjectInfo.mat']);
        loaded_marker_trajectories = load([dataDirectory filesep makeFileName(loaded_subject_info.date, loaded_subject_info.subject_id, 'walking', trialToProcess, 'markerTrajectories')]);
        loaded_force_plate_data = load([dataDirectory filesep makeFileName(loaded_subject_info.date, loaded_subject_info.subject_id, 'walking', trialToProcess, 'forcePlateData')]);
        
        sampling_rate_mocap = loaded_marker_trajectories.sampling_rate_mocap;
        
        % extract data
        left_heel_marker = 34;
        left_toes_marker = 35;
        right_heel_marker = 42;
        right_toes_marker = 43;
        left_heel_marker_indices = reshape([(left_heel_marker - 1) * 3 + 1; (left_heel_marker - 1) * 3 + 2; (left_heel_marker - 1) * 3 + 3], 1, length(left_heel_marker)*3);
        left_toes_marker_indices = reshape([(left_toes_marker - 1) * 3 + 1; (left_toes_marker - 1) * 3 + 2; (left_toes_marker - 1) * 3 + 3], 1, length(left_toes_marker)*3);
        right_heel_marker_indices = reshape([(right_heel_marker - 1) * 3 + 1; (right_heel_marker - 1) * 3 + 2; (right_heel_marker - 1) * 3 + 3], 1, length(right_heel_marker)*3);
        right_toes_marker_indices = reshape([(right_toes_marker - 1) * 3 + 1; (right_toes_marker - 1) * 3 + 2; (right_toes_marker - 1) * 3 + 3], 1, length(right_toes_marker)*3);
        left_heel_z_trajectory = loaded_marker_trajectories.marker_trajectories(:, left_heel_marker_indices(3));
        left_toes_z_trajectory = loaded_marker_trajectories.marker_trajectories(:, left_toes_marker_indices(3));
        right_heel_z_trajectory = loaded_marker_trajectories.marker_trajectories(:, right_heel_marker_indices(3));
        right_toes_z_trajectory = loaded_marker_trajectories.marker_trajectories(:, right_toes_marker_indices(3));

        % calculate derivatives
        filter_order = 2;
        cutoff_frequency = 20; % cutoff frequency, in Hz
        [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate_mocap/2));	% set filter parameters for butterworth filter: 2=order of filter;
        left_heel_z_vel_trajectory = deriveByTime(filtfilt(b, a, left_heel_z_trajectory), 1/sampling_rate_mocap);
        right_heel_z_vel_trajectory = deriveByTime(filtfilt(b, a, right_heel_z_trajectory), 1/sampling_rate_mocap);
        left_heel_z_acc_trajectory = deriveByTime(filtfilt(b, a, left_heel_z_vel_trajectory), 1/sampling_rate_mocap);
        right_heel_z_acc_trajectory = deriveByTime(filtfilt(b, a, right_heel_z_vel_trajectory), 1/sampling_rate_mocap);
        left_toes_z_vel_trajectory = deriveByTime(filtfilt(b, a, left_toes_z_trajectory), 1/sampling_rate_mocap);
        right_toes_z_vel_trajectory = deriveByTime(filtfilt(b, a, right_toes_z_trajectory), 1/sampling_rate_mocap);
        left_toes_z_acc_trajectory = deriveByTime(filtfilt(b, a, left_toes_z_vel_trajectory), 1/sampling_rate_mocap);
        right_toes_z_acc_trajectory = deriveByTime(filtfilt(b, a, right_toes_z_vel_trajectory), 1/sampling_rate_mocap);        
        
        % package
        data = createTrialDataStruct;
        data.time_mocap = loaded_marker_trajectories.time_mocap;
        data.sampling_rate_mocap = sampling_rate_mocap;
        
        data.left_heel_z_pos_trajectory = left_heel_z_trajectory;
        data.left_heel_z_vel_trajectory = left_heel_z_vel_trajectory;
        data.left_heel_z_acc_trajectory = left_heel_z_acc_trajectory;
        
        data.right_heel_z_pos_trajectory = right_heel_z_trajectory;
        data.right_heel_z_vel_trajectory = right_heel_z_vel_trajectory;
        data.right_heel_z_acc_trajectory = right_heel_z_acc_trajectory;
        
        data.left_toes_z_pos_trajectory = left_toes_z_trajectory;
        data.left_toes_z_vel_trajectory = left_toes_z_vel_trajectory;
        data.left_toes_z_acc_trajectory = left_toes_z_acc_trajectory;
        
        data.right_toes_z_pos_trajectory = right_toes_z_trajectory;
        data.right_toes_z_vel_trajectory = right_toes_z_vel_trajectory;
        data.right_toes_z_acc_trajectory = right_toes_z_acc_trajectory;
        
    end





end