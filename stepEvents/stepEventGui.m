function stepEventGui(varargin)
    [condition_list, trial_number_list] = parseTrialArguments(varargin);
    condition = condition_list{1};
    trial_to_process = trial_number_list{1}(1);
    
    if nargin < 3
        color_scheme = 'intra';
    else
        color_scheme = varargin{3};
    end

    
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
        offset_heel = - 0.0;
        offset_toes = - 0.0;
    end
    
    %% load data
    trial_data = WalkingTrialData(pwd, condition, trial_to_process);
    event_data = WalkingEventData(trial_data);
    
    % init gui
    controller = stepEventController(trial_data, event_data);
    
    % create position figures
    step_event_figure = stepEventFigure('Positions Left', controller, trial_data, event_data);
    step_event_figure.addDataPlot('left_heel_z_pos', color_left_heel, scale_factor_heel, offset_heel);
    step_event_figure.addDataPlot('left_toes_z_pos', color_left_toes, scale_factor_toes, offset_toes);
%     step_event_figure.addDataPlot('left_fz', color_left_fz, scale_factor_fz);
    step_event_figure.addEventPlot('left_heel_z_pos', 'left_touchdown', color_left_touchdown, marker_left_touchdown);
    step_event_figure.addEventPlot('left_toes_z_pos', 'left_pushoff', color_left_pushoff, marker_left_pushoff);
%     step_event_figure.addEventPlot('left_fz', 'left_pushoff', color_left_pushoff, marker_left_pushoff);
%     step_event_figure.addEventPlot('left_fz', 'left_touchdown', color_left_touchdown, marker_left_touchdown);
%     step_event_figure.addDataPlot('stimulus_state', color_stimulus_state, scale_factor_stimulus_state);
    
    step_event_figure = stepEventFigure('Positions Right', controller, trial_data, event_data);
    step_event_figure.addDataPlot('right_heel_z_pos', color_right_heel, scale_factor_heel, offset_heel);
    step_event_figure.addDataPlot('right_toes_z_pos', color_right_toes, scale_factor_toes, offset_toes);
%     step_event_figure.addDataPlot('right_fz', color_right_fz, scale_factor_fz);
    step_event_figure.addEventPlot('right_heel_z_pos', 'right_touchdown', color_right_touchdown, marker_right_touchdown);
    step_event_figure.addEventPlot('right_toes_z_pos', 'right_pushoff', color_right_pushoff, marker_right_pushoff);
%     step_event_figure.addEventPlot('right_fz', 'right_pushoff', color_right_pushoff, marker_right_pushoff);
%     step_event_figure.addEventPlot('right_fz', 'right_touchdown', color_right_touchdown, marker_right_touchdown);
%     step_event_figure.addDataPlot('stimulus_state', color_stimulus_state, scale_factor_stimulus_state);

    % create velocity figures
    step_event_figure = stepEventFigure('Velocities Left', controller, trial_data, event_data);
%     step_event_figure.addDataPlot('left_heel_z_vel', color_left_heel);
    step_event_figure.addDataPlot('left_toes_z_vel', color_left_toes);
    step_event_figure.addEventPlot('left_toes_z_vel', 'left_pushoff', color_left_pushoff, marker_left_pushoff);
    step_event_figure.addEventPlot('left_toes_z_vel', 'left_touchdown', color_left_touchdown, marker_left_touchdown);
    
    step_event_figure = stepEventFigure('Velocities Right', controller, trial_data, event_data);
%     step_event_figure.addDataPlot('right_heel_z_vel', color_right_heel);
    step_event_figure.addDataPlot('right_toes_z_vel', color_right_toes);
    step_event_figure.addEventPlot('right_toes_z_vel', 'right_pushoff', color_right_pushoff, marker_right_pushoff);
    step_event_figure.addEventPlot('right_toes_z_vel', 'right_touchdown', color_right_touchdown, marker_right_touchdown);
    
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
    
    
%     controller.findEvents();
    controller.setSelectedEvent();









end