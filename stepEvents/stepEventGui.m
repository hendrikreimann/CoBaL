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
    event_data = WalkingEventData(trial_data);
    
    % init gui
    controller = stepEventController(trial_data, event_data);
    
    % create position figures
    step_event_figure = stepEventFigure('Positions Left', controller, trial_data, event_data);
    step_event_figure.addDataPlot('left_heel_z_pos', color_left_heel);
    step_event_figure.addDataPlot('left_toes_z_pos', color_left_toes);
    step_event_figure.addEventPlot('left_heel_z_pos', 'left_touchdown', color_left_touchdown, marker_left_touchdown);
    step_event_figure.addEventPlot('left_toes_z_pos', 'left_pushoff', color_left_pushoff, marker_left_pushoff);
    
    step_event_figure = stepEventFigure('Positions Right', controller, trial_data, event_data);
    step_event_figure.addDataPlot('right_heel_z_pos', color_right_heel);
    step_event_figure.addDataPlot('right_toes_z_pos', color_right_toes);
    step_event_figure.addEventPlot('right_heel_z_pos', 'right_touchdown', color_right_touchdown, marker_right_touchdown);
    step_event_figure.addEventPlot('right_toes_z_pos', 'right_pushoff', color_right_pushoff, marker_right_pushoff);

    % create velocity figures
    step_event_figure = stepEventFigure('Velocities Left', controller, trial_data, event_data);
%     step_event_figure.addDataPlot('left_heel_z_vel', color_left_heel);
    step_event_figure.addDataPlot('left_toes_z_vel', color_left_toes);
    step_event_figure.addEventPlot('left_toes_z_vel', 'left_pushoff', color_left_pushoff, marker_left_pushoff);
    
    step_event_figure = stepEventFigure('Velocities Right', controller, trial_data, event_data);
%     step_event_figure.addDataPlot('right_heel_z_vel', color_right_heel);
    step_event_figure.addDataPlot('right_toes_z_vel', color_right_toes);
    step_event_figure.addEventPlot('right_toes_z_vel', 'right_pushoff', color_right_pushoff, marker_right_pushoff);
    
    % create acceleration figure
    step_event_figure = stepEventFigure('Acceleration Left', controller, trial_data, event_data);
    step_event_figure.addDataPlot('left_heel_z_acc', color_left_heel);
%     step_event_figure.addDataPlot('left_toes_z_acc', color_left_toes);
    step_event_figure.addEventPlot('left_heel_z_acc', 'left_touchdown', color_left_touchdown, marker_left_touchdown);
    
    step_event_figure = stepEventFigure('Acceleration Right', controller, trial_data, event_data);
    step_event_figure.addDataPlot('right_heel_z_acc', color_right_heel);
%     step_event_figure.addDataPlot('right_toes_z_acc', color_right_toes);
    step_event_figure.addEventPlot('right_heel_z_acc', 'right_touchdown', color_right_touchdown, marker_right_touchdown)
    
    step_event_axes = zeros(size(controller.figureSelectionBox.String));
    for i_figure = 1 : length((controller.figureSelectionBox.String))
        step_event_axes(i_figure) = controller.figureSelectionBox.UserData{i_figure}.main_axes;
    end
    linkaxes(step_event_axes, 'x');
    
    % load settings
    controller.loadFigureSettings();
    controller.findEvents();
    controller.setSelectedEvent();









end