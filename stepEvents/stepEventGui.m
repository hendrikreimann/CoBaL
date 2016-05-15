function stepEventGui(dataDirectory, trialToProcess)
    if nargin < 2
        trialToProcess = 11;
    end
    if nargin < 1
        dataDirectory = pwd;
    end
    
    %% set colors
    color_left_heel = [0.3 0.7 0];
    color_left_toes = [0 0.7 0.3];
    color_left_touchdown = [0.0 0.3 0.7];
    color_left_pushoff = [0.0 0.3 0.7];
    marker_left_touchdown = 'v';
    marker_left_pushoff = '^';

    color_right_heel = [0.7 0.3 0];
    color_right_toes = [0.7 0 0.3];
    color_right_touchdown = [0.3 0 0.7];
    color_right_pushoff = [0.3 0 0.7];
    marker_right_touchdown = 'v';
    marker_right_pushoff = '^';
    
    
    
    %% load data
    trial_data = WalkingTrialData(dataDirectory, trialToProcess);
    
    test_time = trial_data.getTime('left_toes_z_pos');
    test_data = trial_data.getData('left_heel_z_pos');
    
    % load or create event data
    event_data = WalkingEventData(trial_data);
    
    % init gui
    controller = createController;
    
    % create position figures
    step_event_figure = createStepEventFigure('Positions Left');
    step_event_figure.addDataPlot('left_heel_z_pos', color_left_heel);
    step_event_figure.addDataPlot('left_toes_z_pos', color_left_toes);
    step_event_figure.addEventPlot('left_heel_z_pos', 'left_touchdown', color_left_touchdown, marker_left_touchdown);
    step_event_figure.addEventPlot('left_toes_z_pos', 'left_pushoff', color_left_pushoff, marker_left_pushoff);
    
    step_event_figure = createStepEventFigure('Positions Right');
    step_event_figure.addDataPlot('right_heel_z_pos', color_right_heel);
    step_event_figure.addDataPlot('right_toes_z_pos', color_right_toes);
    step_event_figure.addEventPlot('right_heel_z_pos', 'right_touchdown', color_right_touchdown, marker_right_touchdown);
    step_event_figure.addEventPlot('right_toes_z_pos', 'right_pushoff', color_right_pushoff, marker_right_pushoff);

    % create velocity figures
    step_event_figure = createStepEventFigure('Velocities Left');
%     step_event_figure.addDataPlot('left_heel_z_vel', color_left_heel);
    step_event_figure.addDataPlot('left_toes_z_vel', color_left_toes);
    step_event_figure.addEventPlot('left_toes_z_vel', 'left_pushoff', color_left_pushoff, marker_left_pushoff);
    
    step_event_figure = createStepEventFigure('Velocities Right');
%     step_event_figure.addDataPlot('right_heel_z_vel', color_right_heel);
    step_event_figure.addDataPlot('right_toes_z_vel', color_right_toes);
    step_event_figure.addEventPlot('right_toes_z_vel', 'right_pushoff', color_right_pushoff, marker_right_pushoff);
    
    % create acceleration figure
    step_event_figure = createStepEventFigure('Acceleration Left');
    step_event_figure.addDataPlot('left_heel_z_acc', color_left_heel);
%     step_event_figure.addDataPlot('left_toes_z_acc', color_left_toes);
    step_event_figure.addEventPlot('left_heel_z_acc', 'left_touchdown', color_left_touchdown, marker_left_touchdown);
    
    step_event_figure = createStepEventFigure('Acceleration Right');
    step_event_figure.addDataPlot('right_heel_z_acc', color_right_heel);
%     step_event_figure.addDataPlot('right_toes_z_acc', color_right_toes);
    step_event_figure.addEventPlot('right_heel_z_acc', 'right_touchdown', color_right_touchdown, marker_right_touchdown)
    
    step_event_axes = zeros(size(controller.figureSelectionBox.String));
    for i_figure = 1 : length((controller.figureSelectionBox.String))
        step_event_axes(i_figure) = controller.figureSelectionBox.UserData{i_figure}.main_axes;
    end
    linkaxes(step_event_axes, 'x');
    
    % load settings
    loadFigureSettings();
    findEvents();
    controller.setSelectedEvent();

    %% GUI
    function controller = createController()
        controller = stepEventController(trial_data, event_data);
        set(controller.saveFigureSettingsButton, 'callback', @saveFigureSettings);
        set(controller.loadFigureSettingsButton, 'callback', @loadFigureSettings);
        set(controller.findEventsButton, 'callback', @findEvents);
        set(controller.saveEventsButton, 'callback', @saveEvents);
        
        % move all these callbacks to the controller
    end
    function step_event_figure = createStepEventFigure(title)
        step_event_figure = stepEventFigure(title, controller, trial_data, event_data);
        if strcmp(controller.figureSelectionBox.String, '<no figure>')
            controller.figureSelectionBox.String = title;
            controller.figureSelectionBox.UserData = {step_event_figure};
        else
            controller.figureSelectionBox.String = [controller.figureSelectionBox.String; {title}];
            controller.figureSelectionBox.UserData = [controller.figureSelectionBox.UserData; {step_event_figure}];
        end
    end

    function saveFigureSettings(varargin)
        figure_settings = cell(1, length(controller.figureSelectionBox));
        for i_figure = 1 : length(controller.figureSelectionBox.String)
            figure_settings{i_figure} = controller.figureSelectionBox.UserData{i_figure}.getSetting();
        end
        control_figure_setting = struct();
        control_figure_setting.position = controller.control_figure.Position;
        
        save_file = '/Users/reimajbi/Library/Application Support/stepEventFinder/figureSettings.mat';
        save(save_file, 'figure_settings', 'control_figure_setting');
    end
    function loadFigureSettings(varargin)
        load_file = '/Users/reimajbi/Library/Application Support/stepEventFinder/figureSettings.mat';
        
        load(load_file, 'figure_settings', 'control_figure_setting')
        for i_figure = 1 : length(figure_settings)
%             if length(step_event_figures) < i_figure
%                 step_event_figures{i_figure} = createStepEventFigure();
%             end
            controller.figureSelectionBox.UserData{i_figure}.applySettings(figure_settings{i_figure});
        end
        controller.control_figure.Position = control_figure_setting.position;
    end

    %% events
    function findEvents(varargin)
        % left touchdown
        [~, left_touchdown_indices_mocap] = findpeaks(-trial_data.getData('left_heel_z_pos'), 'MinPeakProminence', str2num(controller.heel_pos_peak_width.String));
        time = trial_data.getTime('left_heel_z_pos');
        left_touchdown_times = time(left_touchdown_indices_mocap);
        event_data.setEventTimes(left_touchdown_times, 'left_touchdown');
        
        % left pushoff
        [~, left_pushoff_indices_mocap] = findpeaks(trial_data.getData('left_toes_z_vel'), 'MinPeakProminence', str2num(controller.toes_vel_peak_width.String));
        time = trial_data.getTime('left_toes_z_vel');
        left_pushoff_times = time(left_pushoff_indices_mocap);
        event_data.setEventTimes(left_pushoff_times, 'left_pushoff');
        
        % right touchdown
        [~, right_touchdown_indices_mocap] = findpeaks(-trial_data.getData('right_heel_z_pos'), 'MinPeakProminence', str2num(controller.heel_pos_peak_width.String));
        time = trial_data.getTime('right_heel_z_pos');
        right_touchdown_times = time(right_touchdown_indices_mocap);
        event_data.setEventTimes(right_touchdown_times, 'right_touchdown');

        % right pushoff
        [~, right_pushoff_indices_mocap] = findpeaks(trial_data.getData('right_toes_z_vel'), 'MinPeakProminence', str2num(controller.toes_vel_peak_width.String));
        time = trial_data.getTime('right_toes_z_vel');
        right_pushoff_times = time(right_pushoff_indices_mocap);
        event_data.setEventTimes(right_pushoff_times, 'right_pushoff');

        controller.updateEventPlots();


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