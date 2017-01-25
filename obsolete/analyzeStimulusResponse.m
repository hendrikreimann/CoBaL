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

function analyzeStimulusResponse(varargin)
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    
    save_data                           = 1;
    normalize_emg_mean                  = 1;
    normalize_emg_peak                  = 0;
    
    visualize_com                       = 0;
    % setup analysis pattern - this should be done by editing a file
    process_data_balance                = 1;
    process_data_armswing               = 1;
    process_data_forceplate             = 1;
    process_data_emg                    = 1;
    process_data_angles                 = 0;
    process_data_torques                = 0;
    
    wait_times = [0 0.150 0.450];
    wait_time_labels = {'0ms', '150ms', '450ms'};

    condition_labels = {'stance foot', 'perturbation', 'delay', 'index'};
    
    %% define conditions
    % order of condition is STANCE_FOOT, PERTURBATION, DELAY, INDEX
    % STANCE_FOOT: RIGHT, LEFT
    % PERTURBATION: POSITIVE, NEGATIVE
    % DELAY: 0, 150, 450
    % INDEX: ONE, TWO
    %
    % this should be loaded from a file


    % for phase-dependent GVS
    conditions_control = ...
      {
        'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'walking'; ...
        'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'CONTROL', 'walking'; ...
      };

    conditions_to_analyze = ...
      {
        'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'ONE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_RIGHT', '150ms', 'ONE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_RIGHT', '450ms', 'ONE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'ONE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_LEFT', '150ms', 'ONE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_LEFT', '450ms', 'ONE', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'TWO', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_RIGHT', '150ms', 'TWO', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_RIGHT', '450ms', 'TWO', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'TWO', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_LEFT', '150ms', 'TWO', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_LEFT', '450ms', 'TWO', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'THREE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_RIGHT', '150ms', 'THREE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_RIGHT', '450ms', 'THREE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'THREE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_LEFT', '150ms', 'THREE', 'walking'; ...
        'STANCE_RIGHT', 'ILLUSION_LEFT', '450ms', 'THREE', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'FOUR', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_RIGHT', '150ms', 'FOUR', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_RIGHT', '450ms', 'FOUR', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'FOUR', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_LEFT', '150ms', 'FOUR', 'walking'; ...
        'STANCE_LEFT', 'ILLUSION_LEFT', '450ms', 'FOUR', 'walking'; ...

      };

%     % for vision
%     conditions_control = ...
%       {
%         'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'walking'; ...
%         'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'CONTROL', 'walking'; ...
%       };
% 
%     conditions_to_analyze = ...
%       {
%         'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'ONE', 'walking'; ...
%         'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'ONE', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'ONE', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'ONE', 'walking'; ...
%         'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'TWO', 'walking'; ...
%         'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'TWO', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'TWO', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'TWO', 'walking'; ...
%         'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'THREE', 'walking'; ...
%         'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'THREE', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'THREE', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'THREE', 'walking'; ...
%         'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'FOUR', 'walking'; ...
%         'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'FOUR', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'FOUR', 'walking'; ...
%         'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'FOUR', 'walking'; ...
%       };
    
    number_of_conditions_control = size(conditions_control, 1);
    number_of_conditions_to_analyze = size(conditions_to_analyze, 1);
    second_to_minute = 1/60;

    number_of_time_steps_normalized = 100;
    
    % determine which control condition applies to which perturbation condition
    applicable_control_condition_indices = zeros(number_of_conditions_to_analyze, 1);
    for i_condition = 1 : number_of_conditions_to_analyze
        if strcmp(conditions_to_analyze(i_condition, 1), 'STANCE_LEFT')
            applicable_control_condition_indices(i_condition) = 1;
        elseif strcmp(conditions_to_analyze(i_condition, 1), 'STANCE_RIGHT')
            applicable_control_condition_indices(i_condition) = 2;
        end
    end

    %% extract data
    %
    % This section loads the trial data from disk and extracts the stretches of interest.
    % Stretch start and end points are defined depending upon experimental conditions.
    % Stretch data is time-normalized and stored in arrays with size (number_of_time_steps_normalized x number_of_stretches)
    % named <variable name>_normalized_all
    %

    load('subjectInfo.mat', 'date', 'subject_id');
    
    % initialize containers
%     stretch_length_indices_forceplate = [];
    condition_stance_foot_list_all = {};
    condition_perturbation_list_all = {};
    condition_delay_list_all = {};
    condition_index_list_all = {};
    condition_experimental_list_all = {};
    origin_trial_list_all = [];
    origin_start_time_list_all = [];
    origin_end_time_list_all = [];
    step_times_all = [];
    pushoff_times_all = [];
    cadence_all = [];
    stim_start_time_relative_to_stretch_all = [];    
    
    if process_data_balance
        step_width_all = [];
        step_length_all = [];
        step_speed_all = [];
        
        com_x_pos_normalized_all = [];
        com_y_pos_normalized_all = [];
        
        c7_x_pos_normalized_all = [];
        c7_y_pos_normalized_all = [];
        c7_z_pos_normalized_all = [];
        lpsi_x_pos_normalized_all = [];
        lpsi_y_pos_normalized_all = [];
        lpsi_z_pos_normalized_all = [];
        rpsi_x_pos_normalized_all = [];
        rpsi_y_pos_normalized_all = [];
        rpsi_z_pos_normalized_all = [];
        lheel_x_pos_normalized_all = [];
        lheel_y_pos_normalized_all = [];
        lheel_z_pos_normalized_all = [];
        rheel_x_pos_normalized_all = [];
        rheel_y_pos_normalized_all = [];
        rheel_z_pos_normalized_all = [];
        
        lheel_x_pos_stancefoot_normalized_all = [];
        rheel_x_pos_stancefoot_normalized_all = [];
        lheel_x_pos_mpsis_normalized_all = [];
        rheel_x_pos_mpsis_normalized_all = [];
        lheel_x_pos_com_normalized_all = [];
        rheel_x_pos_com_normalized_all = [];
        
        foot_placement_world_all = [];
        foot_placement_stancefoot_all = [];
        foot_placement_mpsis_all = [];
        foot_placement_com_all = [];
        
        % angles
        lleg_angle_ml_normalized_all = [];
        rleg_angle_ml_normalized_all = [];
        trunk_angle_ap_normalized_all = [];
        trunk_angle_ml_normalized_all = [];
        
    end
    if process_data_armswing
        lelb_x_pos_normalized_all = [];
        lelb_y_pos_normalized_all = [];
        lelb_z_pos_normalized_all = [];
        lwra_x_pos_normalized_all = [];
        lwra_y_pos_normalized_all = [];
        lwra_z_pos_normalized_all = [];
        lwrb_x_pos_normalized_all = [];
        lwrb_y_pos_normalized_all = [];
        lwrb_z_pos_normalized_all = [];
        relb_x_pos_normalized_all = [];
        relb_y_pos_normalized_all = [];
        relb_z_pos_normalized_all = [];
        rwra_x_pos_normalized_all = [];
        rwra_y_pos_normalized_all = [];
        rwra_z_pos_normalized_all = [];
        rwrb_x_pos_normalized_all = [];
        rwrb_y_pos_normalized_all = [];
        rwrb_z_pos_normalized_all = [];
        
        % angles
        linclination_normalized_all = [];
        rinclination_normalized_all = [];
    end
    if process_data_forceplate
        cop_x_normalized_all = [];
        cop_x_stancefoot_normalized_all = [];
        cop_x_mpsis_normalized_all = [];
        cop_x_com_normalized_all = [];
        f_x_normalized_all = [];
        f_z_normalized_all = [];
        m_y_normalized_all = [];
    end
    if process_data_emg
        lglutmed_normalized_all = [];
        ltibiant_normalized_all = [];
        lgastroc_normalized_all = [];
        lperolng_normalized_all = [];
        rglutmed_normalized_all = [];
        rtibiant_normalized_all = [];
        rgastroc_normalized_all = [];
        rperolng_normalized_all = [];
        ltnsrflt_normalized_all = [];
        rtnsrflt_normalized_all = [];
        lglutmed_max_all = [];
        ltibiant_max_all = [];
        lgastroc_max_all = [];
        lperolng_max_all = [];
        rglutmed_max_all = [];
        rtibiant_max_all = [];
        rgastroc_max_all = [];
        rperolng_max_all = [];
        ltnsrflt_max_all = [];
        rtnsrflt_max_all = [];
    end
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            % load data
            condition = condition_list{i_condition};
            load(['analysis' filesep makeFileName(date, subject_id, condition, i_trial, 'relevantDataStretches')]);
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'kinematicTrajectories')]);
            number_of_stretches_trial = length(condition_stance_foot_list);

            condition_stance_foot_list_trial = condition_stance_foot_list;
            condition_perturbation_list_trial = condition_perturbation_list;
            condition_delay_list_trial = condition_delay_list;
            condition_index_list_trial = condition_index_list;
            condition_experimental_list_trial = condition_experimental_list;
            origin_trial_list_trial = zeros(number_of_stretches_trial, 1);
            origin_start_time_list_trial = zeros(number_of_stretches_trial, 1);
            origin_end_time_list_trial = zeros(number_of_stretches_trial, 1);
            step_times_trial = zeros(1, number_of_stretches_trial);
            pushoff_times_trial = zeros(1, number_of_stretches_trial);
            cadence_trial = zeros(1, number_of_stretches_trial);
            
            if process_data_balance

                % define markers and indices
                c7_marker = find(strcmp(marker_headers, 'C7'));
                lpsi_marker = find(strcmp(marker_headers, 'LPSI'));
                rpsi_marker = find(strcmp(marker_headers, 'RPSI'));
                lheel_marker = find(strcmp(marker_headers, 'LHEE'));
                rheel_marker = find(strcmp(marker_headers, 'RHEE'));
                com_body = find(strcmp(com_labels, 'BODYCOM'));

                c7_marker_indices = reshape([(c7_marker - 1) * 3 + 1; (c7_marker - 1) * 3 + 2; (c7_marker - 1) * 3 + 3], 1, length(c7_marker)*3);
                lpsi_marker_indices = reshape([(lpsi_marker - 1) * 3 + 1; (lpsi_marker - 1) * 3 + 2; (lpsi_marker - 1) * 3 + 3], 1, length(lpsi_marker)*3);
                rpsi_marker_indices = reshape([(rpsi_marker - 1) * 3 + 1; (rpsi_marker - 1) * 3 + 2; (rpsi_marker - 1) * 3 + 3], 1, length(rpsi_marker)*3);
                lheel_marker_indices = reshape([(lheel_marker - 1) * 3 + 1; (lheel_marker - 1) * 3 + 2; (lheel_marker - 1) * 3 + 3], 1, length(lheel_marker)*3);
                rheel_marker_indices = reshape([(rheel_marker - 1) * 3 + 1; (rheel_marker - 1) * 3 + 2; (rheel_marker - 1) * 3 + 3], 1, length(rheel_marker)*3);
                com_body_indices = reshape([(com_body - 1) * 3 + 1; (com_body - 1) * 3 + 2; (com_body - 1) * 3 + 3], 1, length(com_body) *3);
                
                % rename relevant trajectories
                c7_x_pos_trajectory = marker_trajectories(:, c7_marker_indices(1));
                c7_y_pos_trajectory = marker_trajectories(:, c7_marker_indices(2));
                c7_z_pos_trajectory = marker_trajectories(:, c7_marker_indices(3));
                lpsi_x_pos_trajectory = marker_trajectories(:, lpsi_marker_indices(1));
                lpsi_y_pos_trajectory = marker_trajectories(:, lpsi_marker_indices(2));
                lpsi_z_pos_trajectory = marker_trajectories(:, lpsi_marker_indices(3));
                rpsi_x_pos_trajectory = marker_trajectories(:, rpsi_marker_indices(1));
                rpsi_y_pos_trajectory = marker_trajectories(:, rpsi_marker_indices(2));
                rpsi_z_pos_trajectory = marker_trajectories(:, rpsi_marker_indices(3));
                lheel_x_pos_trajectory = marker_trajectories(:, lheel_marker_indices(1));
                lheel_y_pos_trajectory = marker_trajectories(:, lheel_marker_indices(2));
                lheel_z_pos_trajectory = marker_trajectories(:, lheel_marker_indices(3));
                rheel_x_pos_trajectory = marker_trajectories(:, rheel_marker_indices(1));
                rheel_y_pos_trajectory = marker_trajectories(:, rheel_marker_indices(2));
                rheel_z_pos_trajectory = marker_trajectories(:, rheel_marker_indices(3));
                com_x_pos_trajectory = com_trajectories(:, com_body_indices(1));
                com_y_pos_trajectory = com_trajectories(:, com_body_indices(2));
                com_z_pos_trajectory = com_trajectories(:, com_body_indices(3));
                
                % initialize containers
                step_width_trial = zeros(1, number_of_stretches_trial);
                step_length_trial = zeros(1, number_of_stretches_trial);
                step_speed_trial = zeros(1, number_of_stretches_trial);
                c7_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                c7_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                c7_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lpsi_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lpsi_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lpsi_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rpsi_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rpsi_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rpsi_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lheel_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lheel_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lheel_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rheel_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rheel_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rheel_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lheel_x_pos_stancefoot_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rheel_x_pos_stancefoot_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lheel_x_pos_mpsis_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rheel_x_pos_mpsis_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lheel_x_pos_com_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rheel_x_pos_com_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lleg_angle_ml_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rleg_angle_ml_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                trunk_angle_ap_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                trunk_angle_ml_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                foot_placement_world_trial = zeros(1, number_of_stretches_trial);
                foot_placement_stancefoot_trial = zeros(1, number_of_stretches_trial);
                foot_placement_mpsis_trial = zeros(1, number_of_stretches_trial);
                foot_placement_com_trial = zeros(1, number_of_stretches_trial);
                
            end                     
            if process_data_armswing
                
                % define markers and indices
                lelb_marker = find(strcmp(marker_headers, 'LELB'));
                lwra_marker = find(strcmp(marker_headers, 'LWRA'));
                lwrb_marker = find(strcmp(marker_headers, 'LWRB'));
                relb_marker = find(strcmp(marker_headers, 'RELB'));
                rwra_marker = find(strcmp(marker_headers, 'RWRA'));
                rwrb_marker = find(strcmp(marker_headers, 'RWRB'));
                
                lelb_marker_indices = reshape([(lelb_marker - 1) * 3 + 1; (lelb_marker - 1) * 3 + 2; (lelb_marker - 1) * 3 + 3], 1, length(lelb_marker)*3);
                lwra_marker_indices = reshape([(lwra_marker - 1) * 3 + 1; (lwra_marker - 1) * 3 + 2; (lwra_marker - 1) * 3 + 3], 1, length(lwra_marker)*3);
                lwrb_marker_indices = reshape([(lwrb_marker - 1) * 3 + 1; (lwrb_marker - 1) * 3 + 2; (lwrb_marker - 1) * 3 + 3], 1, length(lwrb_marker)*3);
                relb_marker_indices = reshape([(relb_marker - 1) * 3 + 1; (relb_marker - 1) * 3 + 2; (relb_marker - 1) * 3 + 3], 1, length(relb_marker)*3);
                rwra_marker_indices = reshape([(rwra_marker - 1) * 3 + 1; (rwra_marker - 1) * 3 + 2; (rwra_marker - 1) * 3 + 3], 1, length(rwra_marker)*3);
                rwrb_marker_indices = reshape([(rwrb_marker - 1) * 3 + 1; (rwrb_marker - 1) * 3 + 2; (rwrb_marker - 1) * 3 + 3], 1, length(rwrb_marker)*3);
                
                % rename relevant trajectories
                lelb_x_pos_trajectory = marker_trajectories(:, lelb_marker_indices(1));
                lelb_y_pos_trajectory = marker_trajectories(:, lelb_marker_indices(2));
                lelb_z_pos_trajectory = marker_trajectories(:, lelb_marker_indices(3));
                lwra_x_pos_trajectory = marker_trajectories(:, lwra_marker_indices(1));
                lwra_y_pos_trajectory = marker_trajectories(:, lwra_marker_indices(2));
                lwra_z_pos_trajectory = marker_trajectories(:, lwra_marker_indices(3));
                lwrb_x_pos_trajectory = marker_trajectories(:, lwrb_marker_indices(1));
                lwrb_y_pos_trajectory = marker_trajectories(:, lwrb_marker_indices(2));
                lwrb_z_pos_trajectory = marker_trajectories(:, lwrb_marker_indices(3));
                relb_x_pos_trajectory = marker_trajectories(:, relb_marker_indices(1));
                relb_y_pos_trajectory = marker_trajectories(:, relb_marker_indices(2));
                relb_z_pos_trajectory = marker_trajectories(:, relb_marker_indices(3));
                rwra_x_pos_trajectory = marker_trajectories(:, rwra_marker_indices(1));
                rwra_y_pos_trajectory = marker_trajectories(:, rwra_marker_indices(2));
                rwra_z_pos_trajectory = marker_trajectories(:, rwra_marker_indices(3));
                rwrb_x_pos_trajectory = marker_trajectories(:, rwrb_marker_indices(1));
                rwrb_y_pos_trajectory = marker_trajectories(:, rwrb_marker_indices(2));
                rwrb_z_pos_trajectory = marker_trajectories(:, rwrb_marker_indices(3));

                % initialize containers
                lelb_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lelb_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lelb_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lwra_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lwra_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lwra_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lwrb_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lwrb_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lwrb_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                relb_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                relb_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                relb_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rwra_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rwra_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rwra_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rwrb_x_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rwrb_y_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rwrb_z_pos_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                
                linclination_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rinclination_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            
            end
            if process_data_forceplate
                % XXX this should be replaced with a system where the subject info file stores the available data for each trial
                forceplate_file_name = ['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'forceplateTrajectories') '.mat'];
                if exist(forceplate_file_name)
                    load(forceplate_file_name);
                    % initialize containers
                    force_plate_data_available = 1;
                else
                    force_plate_data_available = 0;
                    disp(['Warning: forceplate data not available for trial ' num2str(i_trial)]);
                end
                cop_x_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                cop_x_stancefoot_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                cop_x_mpsis_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                cop_x_com_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                f_x_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                f_z_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                m_y_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            end
            if process_data_emg
                emg_file_name = ['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'emgTrajectories') '.mat'];
                if exist(emg_file_name)
                    load(emg_file_name);
                    emg_data_available = 1;
                else
                    emg_data_available = 0;
                end
                % rename relevant trajectories
                lglutmed_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LGLUTMED'));
                ltibiant_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LTIBIANT'));
                lgastroc_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LGASTROC'));
                lperolng_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LPEROLNG'));
                rglutmed_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RGLUTMED'));
                rtibiant_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RTIBIANT'));
                rgastroc_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RGASTROC'));
                rperolng_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RPEROLNG'));
                ltnsrflt_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LTNSRFLT'));
                rtnsrflt_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RTNSRFLT'));
                lglutmed_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                ltibiant_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lgastroc_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                lperolng_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rglutmed_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rtibiant_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rgastroc_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rperolng_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                ltnsrflt_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
                rtnsrflt_normalized_trial = zeros(number_of_time_steps_normalized, number_of_stretches_trial);
            end
            
            % go through each stretch to extract data, time-normalize it and store it in lists
            for i_stretch = 1 : number_of_stretches_trial
                % time
                [~, start_index_mocap] = min(abs(time_mocap - stretch_start_times(i_stretch)));
                [~, pushoff_index_mocap] = min(abs(time_mocap - stretch_pushoff_times(i_stretch)));
                [~, end_index_mocap] = min(abs(time_mocap - stretch_end_times(i_stretch)));
                step_times_trial(i_stretch) = time_mocap(end_index_mocap) - time_mocap(start_index_mocap);
                pushoff_times_trial(i_stretch) = time_mocap(pushoff_index_mocap) - time_mocap(start_index_mocap);
                cadence_trial(i_stretch) = (step_times_trial(i_stretch) * second_to_minute)^(-1);
            
                % balance data
                if process_data_balance
                    % define times
                    time_extracted_mocap = time_mocap(start_index_mocap : end_index_mocap);
                    
                    % extract
                    c7_x_pos_extracted_stretch = c7_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    c7_y_pos_extracted_stretch = c7_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    c7_z_pos_extracted_stretch = c7_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    lpsi_x_pos_extracted_stretch = lpsi_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    lpsi_y_pos_extracted_stretch = lpsi_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    lpsi_z_pos_extracted_stretch = lpsi_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    rpsi_x_pos_extracted_stretch = rpsi_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    rpsi_y_pos_extracted_stretch = rpsi_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    rpsi_z_pos_extracted_stretch = rpsi_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    lheel_x_pos_extracted_stretch = lheel_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    lheel_y_pos_extracted_stretch = lheel_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    lheel_z_pos_extracted_stretch = lheel_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    rheel_x_pos_extracted_stretch = rheel_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    rheel_y_pos_extracted_stretch = rheel_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    rheel_z_pos_extracted_stretch = rheel_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    com_x_pos_extracted_stretch = com_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    com_y_pos_extracted_stretch = com_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    
                    % calculate vectors
                    mpsi_x_pos_extracted_stretch = (lpsi_x_pos_extracted_stretch + rpsi_x_pos_extracted_stretch) * 0.5;
                    mpsi_y_pos_extracted_stretch = (lpsi_y_pos_extracted_stretch + rpsi_y_pos_extracted_stretch) * 0.5;
                    mpsi_z_pos_extracted_stretch = (lpsi_z_pos_extracted_stretch + rpsi_z_pos_extracted_stretch) * 0.5;
                    lleg_vector_x = mpsi_x_pos_extracted_stretch - lheel_x_pos_extracted_stretch;
                    lleg_vector_y = mpsi_y_pos_extracted_stretch - lheel_y_pos_extracted_stretch;
                    lleg_vector_z = mpsi_z_pos_extracted_stretch - lheel_z_pos_extracted_stretch;
                    rleg_vector_x = mpsi_x_pos_extracted_stretch - rheel_x_pos_extracted_stretch;
                    rleg_vector_y = mpsi_y_pos_extracted_stretch - rheel_y_pos_extracted_stretch;
                    rleg_vector_z = mpsi_z_pos_extracted_stretch - rheel_z_pos_extracted_stretch;
                    trunk_vector_x = c7_x_pos_extracted_stretch - mpsi_x_pos_extracted_stretch;
                    trunk_vector_y = c7_y_pos_extracted_stretch - mpsi_y_pos_extracted_stretch;
                    trunk_vector_z = c7_z_pos_extracted_stretch - mpsi_z_pos_extracted_stretch;

                    % calculate angles
                    trunk_angle_ap_extracted_stretch = rad2deg(atan2(trunk_vector_y, trunk_vector_z));
                    trunk_angle_ml_extracted_stretch = rad2deg(atan2(trunk_vector_x, trunk_vector_z));
                    lleg_angle_ml_extracted_stretch = rad2deg(atan2(lleg_vector_x, lleg_vector_z));
                    rleg_angle_ml_extracted_stretch = rad2deg(atan2(rleg_vector_x, rleg_vector_z));
                    
                    % define stance foot heel as spatial point of reference
                    if strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_RIGHT')
                        stance_foot_heel_x_initial = rheel_x_pos_extracted_stretch(1);
                        stance_foot_heel_y_initial = rheel_y_pos_extracted_stretch(1);
                    elseif strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_LEFT')
                        stance_foot_heel_x_initial = lheel_x_pos_extracted_stretch(1);
                        stance_foot_heel_y_initial = lheel_y_pos_extracted_stretch(1);
                    else
                        error('stance condition should be either "STANCE_LEFT or STANCE_RIGHT"');
                    end
                  
                    % calculate step width and length
                    if strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_RIGHT')
                        step_width_trial(i_stretch) = abs(lheel_x_pos_extracted_stretch(end) - rheel_x_pos_extracted_stretch(1));
                        step_length_trial(i_stretch) = lheel_y_pos_extracted_stretch(end) - rheel_y_pos_extracted_stretch(1);
                        foot_placement_world_trial(i_stretch) = lheel_x_pos_extracted_stretch(end);
                        foot_placement_stancefoot_trial(i_stretch) = lheel_x_pos_extracted_stretch(end) - rheel_x_pos_extracted_stretch(1);
                        foot_placement_mpsis_trial(i_stretch) = lheel_x_pos_extracted_stretch(end) - mpsi_x_pos_extracted_stretch(end);
                        foot_placement_com_trial(i_stretch) = lheel_x_pos_extracted_stretch(end) - com_x_pos_extracted_stretch(end);
                    elseif strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_LEFT')
                        step_width_trial(i_stretch) = abs(rheel_x_pos_extracted_stretch(end) - lheel_x_pos_extracted_stretch(1));
                        step_length_trial(i_stretch) = rheel_y_pos_extracted_stretch(end) - lheel_y_pos_extracted_stretch(1);
                        foot_placement_world_trial(i_stretch) = rheel_x_pos_extracted_stretch(end);
                        foot_placement_stancefoot_trial(i_stretch) = rheel_x_pos_extracted_stretch(end) - lheel_x_pos_extracted_stretch(1);
                        foot_placement_mpsis_trial(i_stretch) = rheel_x_pos_extracted_stretch(end) - mpsi_x_pos_extracted_stretch(end);
                        foot_placement_com_trial(i_stretch) = rheel_x_pos_extracted_stretch(end) - com_x_pos_extracted_stretch(end);
                    end

                    % normalize mocap data in time
                    time_normalized_mocap = linspace(time_extracted_mocap(1), time_extracted_mocap(end), number_of_time_steps_normalized);
                    lpsi_x_pos_normalized_stretch = spline(time_extracted_mocap, lpsi_x_pos_extracted_stretch, time_normalized_mocap);
                    lheel_x_pos_normalized_stretch = spline(time_extracted_mocap, lheel_x_pos_extracted_stretch, time_normalized_mocap);
                    lheel_y_pos_normalized_stretch = spline(time_extracted_mocap, lheel_y_pos_extracted_stretch, time_normalized_mocap);
                    rpsi_x_pos_normalized_stretch = spline(time_extracted_mocap, rpsi_x_pos_extracted_stretch, time_normalized_mocap);
                    rheel_x_pos_normalized_stretch = spline(time_extracted_mocap, rheel_x_pos_extracted_stretch, time_normalized_mocap);
                    rheel_y_pos_normalized_stretch = spline(time_extracted_mocap, rheel_y_pos_extracted_stretch, time_normalized_mocap);
                    mpsi_x_pos_normalized_stretch = spline(time_extracted_mocap, mpsi_x_pos_extracted_stretch, time_normalized_mocap);
                    mpsi_y_pos_normalized_stretch = spline(time_extracted_mocap, mpsi_y_pos_extracted_stretch, time_normalized_mocap);
                    com_x_pos_normalized_stretch = spline(time_extracted_mocap, com_x_pos_extracted_stretch, time_normalized_mocap);
                    com_y_pos_normalized_stretch = spline(time_extracted_mocap, com_y_pos_extracted_stretch, time_normalized_mocap);

                    trunk_angle_ap_normalized_stretch = spline(time_extracted_mocap, trunk_angle_ap_extracted_stretch, time_normalized_mocap);
                    trunk_angle_ml_normalized_stretch = spline(time_extracted_mocap, trunk_angle_ml_extracted_stretch, time_normalized_mocap);
                    lleg_angle_ml_normalized_stretch = spline(time_extracted_mocap, lleg_angle_ml_extracted_stretch, time_normalized_mocap);
                    rleg_angle_ml_normalized_stretch = spline(time_extracted_mocap, rleg_angle_ml_extracted_stretch, time_normalized_mocap);
                    
                    % store
                    com_x_pos_normalized_trial(:, i_stretch) = com_x_pos_normalized_stretch;
                    lpsi_x_pos_normalized_trial(:, i_stretch) = lpsi_x_pos_normalized_stretch;
                    rpsi_x_pos_normalized_trial(:, i_stretch) = rpsi_x_pos_normalized_stretch;
                    lheel_x_pos_normalized_trial(:, i_stretch) = lheel_x_pos_normalized_stretch;
                    rheel_x_pos_normalized_trial(:, i_stretch) = rheel_x_pos_normalized_stretch;
                    lheel_y_pos_normalized_trial(:, i_stretch) = lheel_y_pos_normalized_stretch;
                    rheel_y_pos_normalized_trial(:, i_stretch) = rheel_y_pos_normalized_stretch;
                    trunk_angle_ap_normalized_trial(:, i_stretch) = trunk_angle_ap_normalized_stretch;
                    trunk_angle_ml_normalized_trial(:, i_stretch) = trunk_angle_ml_normalized_stretch;
                    lleg_angle_ml_normalized_trial(:, i_stretch) = lleg_angle_ml_normalized_stretch;
                    rleg_angle_ml_normalized_trial(:, i_stretch) = rleg_angle_ml_normalized_stretch;
                    
                    lheel_x_pos_stancefoot_normalized_trial(:, i_stretch) = lheel_x_pos_normalized_stretch - stance_foot_heel_x_initial;
                    rheel_x_pos_stancefoot_normalized_trial(:, i_stretch) = rheel_x_pos_normalized_stretch - stance_foot_heel_x_initial;
                    lheel_x_pos_mpsis_normalized_trial(:, i_stretch) = lheel_x_pos_normalized_stretch - mpsi_x_pos_normalized_stretch;
                    rheel_x_pos_mpsis_normalized_trial(:, i_stretch) = rheel_x_pos_normalized_stretch - mpsi_x_pos_normalized_stretch;
                    lheel_x_pos_com_normalized_trial(:, i_stretch) = lheel_x_pos_normalized_stretch - com_x_pos_normalized_stretch;
                    rheel_x_pos_com_normalized_trial(:, i_stretch) = rheel_x_pos_normalized_stretch - com_x_pos_normalized_stretch;

                    % time and origin
                    origin_trial_list_trial(i_stretch) = i_trial;
                    origin_start_time_list_trial(i_stretch) = time_mocap(start_index_mocap);
                    origin_end_time_list_trial(i_stretch) = time_mocap(end_index_mocap);
                    
                    % calculate step time and speed
                    step_speed_trial(i_stretch) = step_length_trial(i_stretch) / step_times_trial(i_stretch);
                end
                
                % armswing data
                if process_data_armswing
                    % define times
                    [~, start_index_mocap] = min(abs(time_mocap - stretch_start_times(i_stretch)));
                    [~, end_index_mocap] = min(abs(time_mocap - stretch_end_times(i_stretch)));
                    time_extracted_mocap = time_mocap(start_index_mocap : end_index_mocap);
                    
                    % extract
                    lelb_x_pos_extracted_stretch = lelb_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    lelb_y_pos_extracted_stretch = lelb_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    lelb_z_pos_extracted_stretch = lelb_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    lwra_x_pos_extracted_stretch = lwra_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    lwra_y_pos_extracted_stretch = lwra_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    lwra_z_pos_extracted_stretch = lwra_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    lwrb_x_pos_extracted_stretch = lwrb_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    lwrb_y_pos_extracted_stretch = lwrb_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    lwrb_z_pos_extracted_stretch = lwrb_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    relb_x_pos_extracted_stretch = relb_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    relb_y_pos_extracted_stretch = relb_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    relb_z_pos_extracted_stretch = relb_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    rwra_x_pos_extracted_stretch = rwra_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    rwra_y_pos_extracted_stretch = rwra_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    rwra_z_pos_extracted_stretch = rwra_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    rwrb_x_pos_extracted_stretch = rwrb_x_pos_trajectory(start_index_mocap : end_index_mocap);
                    rwrb_y_pos_extracted_stretch = rwrb_y_pos_trajectory(start_index_mocap : end_index_mocap);
                    rwrb_z_pos_extracted_stretch = rwrb_z_pos_trajectory(start_index_mocap : end_index_mocap);
                    
                    % calculate vectors
                    lwrm_x_pos = (lwra_x_pos_extracted_stretch + lwrb_x_pos_extracted_stretch) * 0.5;
                    lwrm_y_pos = (lwra_y_pos_extracted_stretch + lwrb_y_pos_extracted_stretch) * 0.5;
                    lwrm_z_pos = (lwra_z_pos_extracted_stretch + lwrb_z_pos_extracted_stretch) * 0.5;
                    larm_vector_x = lelb_x_pos_extracted_stretch - lwrm_x_pos;
                    larm_vector_y = lelb_y_pos_extracted_stretch - lwrm_y_pos;
                    larm_vector_z = lelb_z_pos_extracted_stretch - lwrm_z_pos;
                    rwrm_x_pos = (rwra_x_pos_extracted_stretch + rwrb_x_pos_extracted_stretch) * 0.5;
                    rwrm_y_pos = (rwra_y_pos_extracted_stretch + rwrb_y_pos_extracted_stretch) * 0.5;
                    rwrm_z_pos = (rwra_z_pos_extracted_stretch + rwrb_z_pos_extracted_stretch) * 0.5;
                    rarm_vector_x = relb_x_pos_extracted_stretch - rwrm_x_pos;
                    rarm_vector_y = relb_y_pos_extracted_stretch - rwrm_y_pos;
                    rarm_vector_z = relb_z_pos_extracted_stretch - rwrm_z_pos;
                    
                    % calculate angles
                    vertical = [0; 0; 1];
                    linclination_extracted_stretch = zeros(length(time_extracted_mocap), 1);
                    rinclination_extracted_stretch = zeros(length(time_extracted_mocap), 1);
                    for i_time = 1 : length(time_extracted_mocap)
                        left_arm_vector = normVector([larm_vector_x(i_time); larm_vector_y(i_time); larm_vector_z(i_time)]);
                        left_arm_vector_projected_length = norm(left_arm_vector(1:2));
                        linclination_extracted_stretch(i_time) = rad2deg(atan2(left_arm_vector_projected_length, left_arm_vector(3)));
                        
                        right_arm_vector = normVector([rarm_vector_x(i_time); rarm_vector_y(i_time); rarm_vector_z(i_time)]);
                        right_arm_vector_projected_length = norm(right_arm_vector(1:2));
                        rinclination_extracted_stretch(i_time) = rad2deg(atan2(right_arm_vector_projected_length, right_arm_vector(3)));
                        
                    end
                    
                    % define stance foot heel as spatial point of reference
                    if strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_RIGHT')
                        stance_foot_heel_x_initial = rheel_x_pos_extracted_stretch(1);
                        stance_foot_heel_y_initial = rheel_y_pos_extracted_stretch(1);
                    elseif strcmp(condition_stance_foot_list_trial{i_stretch}, 'STANCE_LEFT')
                        stance_foot_heel_x_initial = lheel_x_pos_extracted_stretch(1);
                        stance_foot_heel_y_initial = lheel_y_pos_extracted_stretch(1);
                    else
                        error('stance condition should be either "STANCE_LEFT or STANCE_RIGHT"');
                    end
                  
                    % normalize data in time
                    time_normalized_mocap = linspace(time_extracted_mocap(1), time_extracted_mocap(end), number_of_time_steps_normalized);
                    linclination_normalized_stretch = spline(time_extracted_mocap, linclination_extracted_stretch, time_normalized_mocap);
                    rinclination_normalized_stretch = spline(time_extracted_mocap, rinclination_extracted_stretch, time_normalized_mocap);
                    
                    % store
                    linclination_normalized_trial(:, i_stretch) = linclination_normalized_stretch;
                    rinclination_normalized_trial(:, i_stretch) = rinclination_normalized_stretch;
                end

                % forceplate data
                if process_data_forceplate
                    if force_plate_data_available
                        % define times
                        [~, start_index_forceplate] = min(abs(time_forceplate - stretch_start_times(i_stretch)));
                        [~, end_index_forceplate] = min(abs(time_forceplate - stretch_end_times(i_stretch)));
                        time_extracted_forceplate = time_forceplate(start_index_forceplate : end_index_forceplate);

                        % extract
                        cop_x_extracted_stretch = total_forceplate_cop_world(start_index_forceplate : end_index_forceplate, 1);
                        f_x_extracted_stretch = total_forceplate_wrench_world(start_index_forceplate : end_index_forceplate, 1);
                        f_z_extracted_stretch = total_forceplate_wrench_world(start_index_forceplate : end_index_forceplate, 3);
                        m_y_extracted_stretch = total_forceplate_wrench_world(start_index_forceplate : end_index_forceplate, 5);

                        % normalize
                        time_normalized_forceplate = linspace(time_extracted_forceplate(1), time_extracted_forceplate(end), number_of_time_steps_normalized);
                        cop_x_normalized_stretch = spline(time_extracted_forceplate, cop_x_extracted_stretch, time_normalized_forceplate);
                        f_x_normalized_stretch = spline(time_extracted_forceplate, f_x_extracted_stretch, time_normalized_forceplate);
                        f_z_normalized_stretch = spline(time_extracted_forceplate, f_z_extracted_stretch, time_normalized_forceplate);
                        m_y_normalized_stretch = spline(time_extracted_forceplate, m_y_extracted_stretch, time_normalized_forceplate);

                        % use stance foot heel as spatial reference and store
                        cop_x_normalized_trial(:, i_stretch) = cop_x_normalized_stretch;
                        cop_x_stancefoot_normalized_trial(:, i_stretch) = cop_x_normalized_stretch - stance_foot_heel_x_initial;
                        cop_x_mpsis_normalized_trial(:, i_stretch) = cop_x_normalized_stretch - mpsi_x_pos_normalized_stretch;
                        cop_x_com_normalized_trial(:, i_stretch) = cop_x_normalized_stretch - com_x_pos_normalized_stretch;
                        f_x_normalized_trial(:, i_stretch) = f_x_normalized_stretch;
                        f_z_normalized_trial(:, i_stretch) = f_z_normalized_stretch;
                        m_y_normalized_trial(:, i_stretch) = m_y_normalized_stretch;
                    end
                end
                if process_data_emg
                    if emg_data_available
                        % define times
                        [~, start_index_emg] = min(abs(time_emg - stretch_start_times(i_stretch)));
                        [~, end_index_emg] = min(abs(time_emg - stretch_end_times(i_stretch)));
                        time_extracted_emg = time_emg(start_index_emg : end_index_emg);

                        % extract
                        lglutmed_extracted_stretch = lglutmed_trajectory(start_index_emg : end_index_emg, 1);
                        ltibiant_extracted_stretch = ltibiant_trajectory(start_index_emg : end_index_emg, 1);
                        lgastroc_extracted_stretch = lgastroc_trajectory(start_index_emg : end_index_emg, 1);
                        lperolng_extracted_stretch = lperolng_trajectory(start_index_emg : end_index_emg, 1);
                        rglutmed_extracted_stretch = rglutmed_trajectory(start_index_emg : end_index_emg, 1);
                        rtibiant_extracted_stretch = rtibiant_trajectory(start_index_emg : end_index_emg, 1);
                        rgastroc_extracted_stretch = rgastroc_trajectory(start_index_emg : end_index_emg, 1);
                        rperolng_extracted_stretch = rperolng_trajectory(start_index_emg : end_index_emg, 1);
                        ltnsrflt_extracted_stretch = ltnsrflt_trajectory(start_index_emg : end_index_emg, 1);
                        rtnsrflt_extracted_stretch = rtnsrflt_trajectory(start_index_emg : end_index_emg, 1);

                        % normalize
                        time_normalized_emg = linspace(time_extracted_emg(1), time_extracted_emg(end), number_of_time_steps_normalized);
                        lglutmed_normalized_stretch = spline(time_extracted_emg, lglutmed_extracted_stretch, time_normalized_emg);
                        ltibiant_normalized_stretch = spline(time_extracted_emg, ltibiant_extracted_stretch, time_normalized_emg);
                        lgastroc_normalized_stretch = spline(time_extracted_emg, lgastroc_extracted_stretch, time_normalized_emg);
                        lperolng_normalized_stretch = spline(time_extracted_emg, lperolng_extracted_stretch, time_normalized_emg);
                        rglutmed_normalized_stretch = spline(time_extracted_emg, rglutmed_extracted_stretch, time_normalized_emg);
                        rtibiant_normalized_stretch = spline(time_extracted_emg, rtibiant_extracted_stretch, time_normalized_emg);
                        rgastroc_normalized_stretch = spline(time_extracted_emg, rgastroc_extracted_stretch, time_normalized_emg);
                        rperolng_normalized_stretch = spline(time_extracted_emg, rperolng_extracted_stretch, time_normalized_emg);
                        ltnsrflt_normalized_stretch = spline(time_extracted_emg, ltnsrflt_extracted_stretch, time_normalized_emg);
                        rtnsrflt_normalized_stretch = spline(time_extracted_emg, rtnsrflt_extracted_stretch, time_normalized_emg);

                        % store
                        lglutmed_normalized_trial(:, i_stretch) = lglutmed_normalized_stretch;
                        ltibiant_normalized_trial(:, i_stretch) = ltibiant_normalized_stretch;
                        lgastroc_normalized_trial(:, i_stretch) = lgastroc_normalized_stretch;
                        lperolng_normalized_trial(:, i_stretch) = lperolng_normalized_stretch;
                        rglutmed_normalized_trial(:, i_stretch) = rglutmed_normalized_stretch;
                        rtibiant_normalized_trial(:, i_stretch) = rtibiant_normalized_stretch;
                        rgastroc_normalized_trial(:, i_stretch) = rgastroc_normalized_stretch;
                        rperolng_normalized_trial(:, i_stretch) = rperolng_normalized_stretch;
                        ltnsrflt_normalized_trial(:, i_stretch) = ltnsrflt_normalized_stretch;
                        rtnsrflt_normalized_trial(:, i_stretch) = rtnsrflt_normalized_stretch;
                    end
                end
            end
                
            % append trial containers to total containers
            condition_stance_foot_list_all = [condition_stance_foot_list_all; condition_stance_foot_list_trial];
            condition_perturbation_list_all = [condition_perturbation_list_all; condition_perturbation_list_trial];
            condition_delay_list_all = [condition_delay_list_all; condition_delay_list_trial];
            condition_index_list_all = [condition_index_list_all; condition_index_list_trial];
            condition_experimental_list_all = [condition_experimental_list_all; condition_experimental_list_trial];
            origin_trial_list_all = [origin_trial_list_all; origin_trial_list_trial];
            origin_start_time_list_all = [origin_start_time_list_all; origin_start_time_list_trial];
            origin_end_time_list_all = [origin_end_time_list_all; origin_end_time_list_trial];
            step_times_all = [step_times_all step_times_trial];
            pushoff_times_all = [pushoff_times_all pushoff_times_trial];
            cadence_all = [cadence_all cadence_trial];

            if process_data_balance
                step_width_all = [step_width_all step_width_trial];
                step_length_all = [step_length_all step_length_trial];
                step_speed_all = [step_speed_all step_speed_trial];
                
                com_x_pos_normalized_all = [com_x_pos_normalized_all com_x_pos_normalized_trial];
                lpsi_x_pos_normalized_all = [lpsi_x_pos_normalized_all lpsi_x_pos_normalized_trial];
                rpsi_x_pos_normalized_all = [rpsi_x_pos_normalized_all rpsi_x_pos_normalized_trial];
                lheel_x_pos_normalized_all = [lheel_x_pos_normalized_all lheel_x_pos_normalized_trial];
                rheel_x_pos_normalized_all = [rheel_x_pos_normalized_all rheel_x_pos_normalized_trial];
                lheel_y_pos_normalized_all = [lheel_y_pos_normalized_all lheel_y_pos_normalized_trial];
                rheel_y_pos_normalized_all = [rheel_y_pos_normalized_all rheel_y_pos_normalized_trial];
                
                lheel_x_pos_stancefoot_normalized_all = [lheel_x_pos_stancefoot_normalized_all lheel_x_pos_stancefoot_normalized_trial];
                rheel_x_pos_stancefoot_normalized_all = [rheel_x_pos_stancefoot_normalized_all rheel_x_pos_stancefoot_normalized_trial];
                lheel_x_pos_mpsis_normalized_all = [lheel_x_pos_mpsis_normalized_all lheel_x_pos_mpsis_normalized_trial];
                rheel_x_pos_mpsis_normalized_all = [rheel_x_pos_mpsis_normalized_all rheel_x_pos_mpsis_normalized_trial];
                lheel_x_pos_com_normalized_all = [lheel_x_pos_com_normalized_all lheel_x_pos_com_normalized_trial];
                rheel_x_pos_com_normalized_all = [rheel_x_pos_com_normalized_all rheel_x_pos_com_normalized_trial];
                
                trunk_angle_ap_normalized_all = [trunk_angle_ap_normalized_all trunk_angle_ap_normalized_trial];
                trunk_angle_ml_normalized_all = [trunk_angle_ml_normalized_all trunk_angle_ml_normalized_trial];
                lleg_angle_ml_normalized_all = [lleg_angle_ml_normalized_all lleg_angle_ml_normalized_trial];
                rleg_angle_ml_normalized_all = [rleg_angle_ml_normalized_all rleg_angle_ml_normalized_trial];
                
                foot_placement_world_all = [foot_placement_world_all foot_placement_world_trial];
                foot_placement_stancefoot_all = [foot_placement_stancefoot_all foot_placement_stancefoot_trial];
                foot_placement_mpsis_all = [foot_placement_mpsis_all foot_placement_mpsis_trial];
                foot_placement_com_all = [foot_placement_com_all foot_placement_com_trial];
            end

            if process_data_armswing
                linclination_normalized_all = [linclination_normalized_all linclination_normalized_trial];
                rinclination_normalized_all = [rinclination_normalized_all rinclination_normalized_trial];
            end
                
            if process_data_forceplate
                cop_x_normalized_all = [cop_x_normalized_all cop_x_normalized_trial];
                cop_x_stancefoot_normalized_all = [cop_x_stancefoot_normalized_all cop_x_stancefoot_normalized_trial];
                cop_x_mpsis_normalized_all = [cop_x_mpsis_normalized_all cop_x_mpsis_normalized_trial];
                cop_x_com_normalized_all = [cop_x_com_normalized_all cop_x_com_normalized_trial];
                f_x_normalized_all = [f_x_normalized_all f_x_normalized_trial];
                f_z_normalized_all = [f_z_normalized_all f_z_normalized_trial];
                m_y_normalized_all = [m_y_normalized_all m_y_normalized_trial];
            end

            if process_data_emg
                lglutmed_normalized_all = [lglutmed_normalized_all lglutmed_normalized_trial];
                ltibiant_normalized_all = [ltibiant_normalized_all ltibiant_normalized_trial];
                lgastroc_normalized_all = [lgastroc_normalized_all lgastroc_normalized_trial];
                lperolng_normalized_all = [lperolng_normalized_all lperolng_normalized_trial];
                rglutmed_normalized_all = [rglutmed_normalized_all rglutmed_normalized_trial];
                rtibiant_normalized_all = [rtibiant_normalized_all rtibiant_normalized_trial];
                rgastroc_normalized_all = [rgastroc_normalized_all rgastroc_normalized_trial];
                rperolng_normalized_all = [rperolng_normalized_all rperolng_normalized_trial];                
                ltnsrflt_normalized_all = [ltnsrflt_normalized_all ltnsrflt_normalized_trial];                
                rtnsrflt_normalized_all = [rtnsrflt_normalized_all rtnsrflt_normalized_trial];                
            end
            
            disp(['Condition ' condition ', Trial ' num2str(i_trial) ' done extracted ' num2str(length(step_times_trial)) ' stretches']);
        end
    end
    number_of_stretches = length(step_times_all);
    
    
    %% extract conditions
    %
    % Compile information about which stretch belongs to which condition
    %
    % allocate indicators
    conditions_control_indicators = false(number_of_stretches, number_of_conditions_control);
    conditions_to_analyze_indicators = false(number_of_stretches, number_of_conditions_to_analyze);

    % extract control condition indicators
    for i_condition = 1 : number_of_conditions_control
        stance_foot_indicator = strcmp(condition_stance_foot_list_all, conditions_control(i_condition, 1));
        perturbation_indicator = strcmp(condition_perturbation_list_all, conditions_control(i_condition, 2));
        delay_indicator = strcmp(condition_delay_list_all, conditions_control(i_condition, 3));
        index_indicator = strcmp(condition_index_list_all, conditions_control(i_condition, 4));

        this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator;
        conditions_control_indicators(:, i_condition) = this_condition_indicator;
    end

    % extract indicators for conditions to analyze
    for i_condition = 1 : number_of_conditions_to_analyze
        stance_foot_indicator = strcmp(condition_stance_foot_list_all, conditions_to_analyze(i_condition, 1));
        perturbation_indicator = strcmp(condition_perturbation_list_all, conditions_to_analyze(i_condition, 2));
        delay_indicator = strcmp(condition_delay_list_all, conditions_to_analyze(i_condition, 3));
        index_indicator = strcmp(condition_index_list_all, conditions_to_analyze(i_condition, 4));

        this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator;
        conditions_to_analyze_indicators(:, i_condition) = this_condition_indicator;
    end

    % give feedback about number of trials per condition
    trials_per_condition_control = sum(conditions_control_indicators)';
    conditions_control_with_number = conditions_control;
    for i_condition = 1 : number_of_conditions_control
        conditions_control_with_number{i_condition, 5} = num2str(trials_per_condition_control(i_condition));
    end
    conditions_control_with_labels = [condition_labels 'number of stretches'; conditions_control_with_number];

    trials_per_condition_to_analyze = sum(conditions_to_analyze_indicators)';
    conditions_to_analyze_with_number = conditions_to_analyze;
    for i_condition = 1 : number_of_conditions_to_analyze
        conditions_to_analyze_with_number{i_condition, 5} = num2str(trials_per_condition_to_analyze(i_condition));
    end
    conditions_to_analyze_with_labels = [condition_labels 'number of stretches'; conditions_to_analyze_with_number];

    disp('Control conditions:')
    disp(conditions_control_with_labels);

    disp('Conditions to analyze:')
    disp(conditions_to_analyze_with_labels);

    disp(['Number of control stretches: ' num2str(sum(trials_per_condition_control))]);
    disp(['Number of stimulus stretches: ' num2str(sum(trials_per_condition_to_analyze))]);
    disp(['Number of unassigned stretches: ' num2str(number_of_stretches - sum(trials_per_condition_control) - sum(trials_per_condition_to_analyze))]);

    %% Normalize EMG
    if (normalize_emg_mean || normalize_emg_peak) && (process_data_emg)
        if normalize_emg_mean && normalize_emg_peak
            error('Choose only ONE normalization method at top menu')
        end
        % LOAD BASELINE and ASSIGN
        i_trial = 1;
        condition = 'emg';
        load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'emgTrajectories')]);
        lglutmed_baseline_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LGLUTMED'));
        ltibiant_baseline_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LTIBIANT'));
        lgastroc_baseline_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LGASTROC'));
        lperolng_baseline_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LPEROLNG'));
        rglutmed_baseline_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RGLUTMED'));
        rtibiant_baseline_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RTIBIANT'));
        rgastroc_baseline_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RGASTROC'));
        rperolng_baseline_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RPEROLNG'));
        ltnsrflt_baseline_trajectory = emg_trajectories(:, strcmp(emg_headers, 'LTNSRFLT'));
        rtnsrflt_baseline_trajectory = emg_trajectories(:, strcmp(emg_headers, 'RTNSRFLT'));
        
        % NORMALIZED EMG LOW
        lglutmed_emg_low = mean(lglutmed_baseline_trajectory,1);
        ltibiant_emg_low = mean(ltibiant_baseline_trajectory,1);
        lgastroc_emg_low = mean(lgastroc_baseline_trajectory,1);
        lperolng_emg_low = mean(lperolng_baseline_trajectory,1);
        rglutmed_emg_low = mean(rglutmed_baseline_trajectory,1);
        rtibiant_emg_low = mean(rtibiant_baseline_trajectory,1);
        rgastroc_emg_low = mean(rgastroc_baseline_trajectory,1);
        rperolng_emg_low = mean(rperolng_baseline_trajectory,1);
        ltnsrflt_emg_low = mean(ltnsrflt_baseline_trajectory,1);
        rtnsrflt_emg_low = mean(rtnsrflt_baseline_trajectory,1);
        
        % CONTROL MEAN LEFT AND RIGHT
        % Need to find control indicator and stance foot setup and proper averaging across
        % all time series
        lglutmed_control_mean_stanceleft = mean(lglutmed_normalized_all(:,conditions_control_indicators(:,1)),2); 
        ltibiant_control_mean_stanceleft = mean(ltibiant_normalized_all(:,conditions_control_indicators(:,1)),2);
        lgastroc_control_mean_stanceleft = mean(lgastroc_normalized_all(:,conditions_control_indicators(:,1)),2);
        lperolng_control_mean_stanceleft = mean(lperolng_normalized_all(:,conditions_control_indicators(:,1)),2);
        rglutmed_control_mean_stanceleft = mean(rglutmed_normalized_all(:,conditions_control_indicators(:,1)),2);
        rtibiant_control_mean_stanceleft = mean(rtibiant_normalized_all(:,conditions_control_indicators(:,1)),2);
        rgastroc_control_mean_stanceleft = mean(rgastroc_normalized_all(:,conditions_control_indicators(:,1)),2);
        rperolng_control_mean_stanceleft = mean(rperolng_normalized_all(:,conditions_control_indicators(:,1)),2);  
        ltnsrflt_control_mean_stanceleft = mean(ltnsrflt_normalized_all(:,conditions_control_indicators(:,1)),2);  
        rtnsrflt_control_mean_stanceleft = mean(rtnsrflt_normalized_all(:,conditions_control_indicators(:,1)),2);  

        lglutmed_control_mean_stanceright = mean(lglutmed_normalized_all(:,conditions_control_indicators(:,2)),2); 
        ltibiant_control_mean_stanceright = mean(ltibiant_normalized_all(:,conditions_control_indicators(:,2)),2);
        lgastroc_control_mean_stanceright = mean(lgastroc_normalized_all(:,conditions_control_indicators(:,2)),2);
        lperolng_control_mean_stanceright = mean(lperolng_normalized_all(:,conditions_control_indicators(:,2)),2);
        rglutmed_control_mean_stanceright = mean(rglutmed_normalized_all(:,conditions_control_indicators(:,2)),2);
        rtibiant_control_mean_stanceright = mean(rtibiant_normalized_all(:,conditions_control_indicators(:,2)),2);
        rgastroc_control_mean_stanceright = mean(rgastroc_normalized_all(:,conditions_control_indicators(:,2)),2);
        rperolng_control_mean_stanceright = mean(rperolng_normalized_all(:,conditions_control_indicators(:,2)),2);
        ltnsrflt_control_mean_stanceright = mean(ltnsrflt_normalized_all(:,conditions_control_indicators(:,2)),2);
        rtnsrflt_control_mean_stanceright = mean(rtnsrflt_normalized_all(:,conditions_control_indicators(:,2)),2);

        if normalize_emg_mean
            % NORMALIZED EMG HIGH LEFT AND RIGHT
            % Average across time
            lglutmed_emg_high_stanceleft = mean(lglutmed_control_mean_stanceleft,1);
            ltibiant_emg_high_stanceleft = mean(ltibiant_control_mean_stanceleft,1);
            lgastroc_emg_high_stanceleft = mean(lgastroc_control_mean_stanceleft,1);
            lperolng_emg_high_stanceleft = mean(lperolng_control_mean_stanceleft,1);
            rglutmed_emg_high_stanceleft = mean(rglutmed_control_mean_stanceleft,1);
            rtibiant_emg_high_stanceleft = mean(rtibiant_control_mean_stanceleft,1);
            rgastroc_emg_high_stanceleft = mean(rgastroc_control_mean_stanceleft,1);
            rperolng_emg_high_stanceleft = mean(rperolng_control_mean_stanceleft,1);
            ltnsrflt_emg_high_stanceleft = mean(ltnsrflt_control_mean_stanceleft,1);
            rtnsrflt_emg_high_stanceleft = mean(rtnsrflt_control_mean_stanceleft,1);

            lglutmed_emg_high_stanceright = mean(lglutmed_control_mean_stanceright,1);
            ltibiant_emg_high_stanceright = mean(ltibiant_control_mean_stanceright,1);
            lgastroc_emg_high_stanceright = mean(lgastroc_control_mean_stanceright,1);
            lperolng_emg_high_stanceright = mean(lperolng_control_mean_stanceright,1);
            rglutmed_emg_high_stanceright = mean(rglutmed_control_mean_stanceright,1);
            rtibiant_emg_high_stanceright = mean(rtibiant_control_mean_stanceright,1);
            rgastroc_emg_high_stanceright = mean(rgastroc_control_mean_stanceright,1);
            rperolng_emg_high_stanceright = mean(rperolng_control_mean_stanceright,1);
            ltnsrflt_emg_high_stanceright = mean(ltnsrflt_control_mean_stanceright,1);
            rtnsrflt_emg_high_stanceright = mean(rtnsrflt_control_mean_stanceright,1);

            % NORMALIZED EMG HIGH
            % Avg two values for each muscle
            lglutmed_emg_high = mean([lglutmed_emg_high_stanceleft, lglutmed_emg_high_stanceright]);
            ltibiant_emg_high = mean([ltibiant_emg_high_stanceleft, ltibiant_emg_high_stanceright]);
            lgastroc_emg_high = mean([lgastroc_emg_high_stanceleft, lgastroc_emg_high_stanceright]);
            lperolng_emg_high = mean([lperolng_emg_high_stanceleft, lperolng_emg_high_stanceright]);
            rglutmed_emg_high = mean([rglutmed_emg_high_stanceleft, rglutmed_emg_high_stanceright]);
            rtibiant_emg_high = mean([rtibiant_emg_high_stanceleft, rtibiant_emg_high_stanceright]);
            rgastroc_emg_high = mean([rgastroc_emg_high_stanceleft, rgastroc_emg_high_stanceright]);
            rperolng_emg_high = mean([rperolng_emg_high_stanceleft, rperolng_emg_high_stanceright]);
            ltnsrflt_emg_high = mean([ltnsrflt_emg_high_stanceleft, ltnsrflt_emg_high_stanceright]);
            rtnsrflt_emg_high = mean([rtnsrflt_emg_high_stanceleft, rtnsrflt_emg_high_stanceright]);
        end
       
        if normalize_emg_peak
            % NORMALIZED EMG HIGH LEFT AND RIGHT
            % Average across time
            lglutmed_emg_high_stanceleft = max(lglutmed_control_mean_stanceleft);
            ltibiant_emg_high_stanceleft = max(ltibiant_control_mean_stanceleft);
            lgastroc_emg_high_stanceleft = max(lgastroc_control_mean_stanceleft);
            lperolng_emg_high_stanceleft = max(lperolng_control_mean_stanceleft);
            rglutmed_emg_high_stanceleft = max(rglutmed_control_mean_stanceleft);
            rtibiant_emg_high_stanceleft = max(rtibiant_control_mean_stanceleft);
            rgastroc_emg_high_stanceleft = max(rgastroc_control_mean_stanceleft);
            rperolng_emg_high_stanceleft = max(rperolng_control_mean_stanceleft);
            ltnsrflt_emg_high_stanceleft = max(ltnsrflt_control_mean_stanceleft);
            rtnsrflt_emg_high_stanceleft = max(rtnsrflt_control_mean_stanceleft);

            lglutmed_emg_high_stanceright = max(lglutmed_control_mean_stanceright);
            ltibiant_emg_high_stanceright = max(ltibiant_control_mean_stanceright);
            lgastroc_emg_high_stanceright = max(lgastroc_control_mean_stanceright);
            lperolng_emg_high_stanceright = max(lperolng_control_mean_stanceright);
            rglutmed_emg_high_stanceright = max(rglutmed_control_mean_stanceright);
            rtibiant_emg_high_stanceright = max(rtibiant_control_mean_stanceright);
            rgastroc_emg_high_stanceright = max(rgastroc_control_mean_stanceright);
            rperolng_emg_high_stanceright = max(rperolng_control_mean_stanceright);
            ltnsrflt_emg_high_stanceright = max(ltnsrflt_control_mean_stanceright);
            rtnsrflt_emg_high_stanceright = max(rtnsrflt_control_mean_stanceright);

            % NORMALIZED EMG HIGH
            % Avg two values for each muscle
            lglutmed_emg_high = max([lglutmed_emg_high_stanceleft, lglutmed_emg_high_stanceright]);
            ltibiant_emg_high = max([ltibiant_emg_high_stanceleft, ltibiant_emg_high_stanceright]);
            lgastroc_emg_high = max([lgastroc_emg_high_stanceleft, lgastroc_emg_high_stanceright]);
            lperolng_emg_high = max([lperolng_emg_high_stanceleft, lperolng_emg_high_stanceright]);
            rglutmed_emg_high = max([rglutmed_emg_high_stanceleft, rglutmed_emg_high_stanceright]);
            rtibiant_emg_high = max([rtibiant_emg_high_stanceleft, rtibiant_emg_high_stanceright]);
            rgastroc_emg_high = max([rgastroc_emg_high_stanceleft, rgastroc_emg_high_stanceright]);
            rperolng_emg_high = max([rperolng_emg_high_stanceleft, rperolng_emg_high_stanceright]);
            ltnsrflt_emg_high = max([ltnsrflt_emg_high_stanceleft, ltnsrflt_emg_high_stanceright]);
            rtnsrflt_emg_high = max([rtnsrflt_emg_high_stanceleft, rtnsrflt_emg_high_stanceright]);

        end
            % RESCALE TRAJECTORIES
            lglutmed_normalized_rescaled_all = ((lglutmed_normalized_all - lglutmed_emg_low)/(lglutmed_emg_high - lglutmed_emg_low));
            ltibiant_normalized_rescaled_all = ((ltibiant_normalized_all - ltibiant_emg_low)/(ltibiant_emg_high - ltibiant_emg_low));
            lgastroc_normalized_rescaled_all = ((lgastroc_normalized_all - lgastroc_emg_low)/(lgastroc_emg_high - lgastroc_emg_low));
            lperolng_normalized_rescaled_all = ((lperolng_normalized_all - lperolng_emg_low)/(lperolng_emg_high - lperolng_emg_low));
            rglutmed_normalized_rescaled_all = ((rglutmed_normalized_all - rglutmed_emg_low)/(rglutmed_emg_high - rglutmed_emg_low));
            rtibiant_normalized_rescaled_all = ((rtibiant_normalized_all - rtibiant_emg_low)/(rtibiant_emg_high - rtibiant_emg_low));
            rgastroc_normalized_rescaled_all = ((rgastroc_normalized_all - rgastroc_emg_low)/(rgastroc_emg_high - rgastroc_emg_low));
            rperolng_normalized_rescaled_all = ((rperolng_normalized_all - rperolng_emg_low)/(rperolng_emg_high - rperolng_emg_low));
            ltnsrflt_normalized_rescaled_all = ((ltnsrflt_normalized_all - ltnsrflt_emg_low)/(ltnsrflt_emg_high - ltnsrflt_emg_low));
            rtnsrflt_normalized_rescaled_all = ((rtnsrflt_normalized_all - rtnsrflt_emg_low)/(rtnsrflt_emg_high - rtnsrflt_emg_low));
    end
    
    
    %% calculate responses
    if process_data_balance
        % calculate control means
        step_width_control_means = zeros(1, number_of_conditions_control);
        foot_placement_world_control_means = zeros(1, number_of_conditions_control);
        foot_placement_stancefoot_control_means = zeros(1, number_of_conditions_control);
        foot_placement_mpsis_control_means = zeros(1, number_of_conditions_control);
        foot_placement_com_control_means = zeros(1, number_of_conditions_control);
        com_x_pos_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lheel_x_pos_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rheel_x_pos_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lheel_y_pos_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rheel_y_pos_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lheel_x_pos_stancefoot_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rheel_x_pos_stancefoot_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lheel_x_pos_mpsis_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rheel_x_pos_mpsis_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lheel_x_pos_com_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rheel_x_pos_com_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        trunk_angle_ml_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lleg_angle_ml_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rleg_angle_ml_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            step_width_control_means(:, i_condition) = mean(step_width_all(:, condition_indicator), 2);
            foot_placement_world_control_means(:, i_condition) = mean(foot_placement_world_all(:, condition_indicator), 2);
            foot_placement_stancefoot_control_means(:, i_condition) = mean(foot_placement_stancefoot_all(:, condition_indicator), 2);
            foot_placement_mpsis_control_means(:, i_condition) = mean(foot_placement_mpsis_all(:, condition_indicator), 2);
            foot_placement_com_control_means(:, i_condition) = mean(foot_placement_com_all(:, condition_indicator), 2);
            com_x_pos_control_means(:, i_condition) = mean(com_x_pos_normalized_all(:, condition_indicator), 2);
            lheel_x_pos_control_means(:, i_condition) = mean(lheel_x_pos_normalized_all(:, condition_indicator), 2);
            rheel_x_pos_control_means(:, i_condition) = mean(rheel_x_pos_normalized_all(:, condition_indicator), 2);
            lheel_y_pos_control_means(:, i_condition) = mean(lheel_y_pos_normalized_all(:, condition_indicator), 2);
            rheel_y_pos_control_means(:, i_condition) = mean(rheel_y_pos_normalized_all(:, condition_indicator), 2);
            lheel_x_pos_stancefoot_control_means(:, i_condition) = mean(lheel_x_pos_stancefoot_normalized_all(:, condition_indicator), 2);
            rheel_x_pos_stancefoot_control_means(:, i_condition) = mean(rheel_x_pos_stancefoot_normalized_all(:, condition_indicator), 2);
            lheel_x_pos_mpsis_control_means(:, i_condition) = mean(lheel_x_pos_mpsis_normalized_all(:, condition_indicator), 2);
            rheel_x_pos_mpsis_control_means(:, i_condition) = mean(rheel_x_pos_mpsis_normalized_all(:, condition_indicator), 2);  
            lheel_x_pos_com_control_means(:, i_condition) = mean(lheel_x_pos_com_normalized_all(:, condition_indicator), 2);
            rheel_x_pos_com_control_means(:, i_condition) = mean(rheel_x_pos_com_normalized_all(:, condition_indicator), 2);
            trunk_angle_ml_control_means(:, i_condition) = mean(trunk_angle_ml_normalized_all(:, condition_indicator), 2);
            lleg_angle_ml_control_means(:, i_condition) = mean(lleg_angle_ml_normalized_all(:, condition_indicator), 2);
            rleg_angle_ml_control_means(:, i_condition) = mean(rleg_angle_ml_normalized_all(:, condition_indicator), 2);
        end
        
        % calculate stimulus responses
        step_width_response = zeros(1, length(step_width_all));
        foot_placement_world_response = zeros(1, length(step_width_all));
        foot_placement_stancefoot_response = zeros(1, length(step_width_all));
        foot_placement_mpsis_response = zeros(1, length(step_width_all));
        foot_placement_com_response = zeros(1, length(step_width_all));
        com_x_pos_response = zeros(size(lheel_x_pos_normalized_all));
        lheel_x_pos_response = zeros(size(lheel_x_pos_normalized_all));
        rheel_x_pos_response = zeros(size(rheel_x_pos_normalized_all));
        lheel_y_pos_response = zeros(size(lheel_y_pos_normalized_all));
        rheel_y_pos_response = zeros(size(rheel_y_pos_normalized_all));
        lheel_x_pos_stancefoot_response = zeros(size(lheel_x_pos_normalized_all));
        rheel_x_pos_stancefoot_response = zeros(size(rheel_x_pos_normalized_all));
        lheel_x_pos_mpsis_response = zeros(size(lheel_x_pos_normalized_all));
        rheel_x_pos_mpsis_response = zeros(size(rheel_x_pos_normalized_all));
        lheel_x_pos_com_response = zeros(size(lheel_x_pos_normalized_all));
        rheel_x_pos_com_response = zeros(size(rheel_x_pos_normalized_all));
        trunk_angle_ml_response = zeros(size(trunk_angle_ml_normalized_all));
        lleg_angle_ml_response = zeros(size(lleg_angle_ml_normalized_all));
        rleg_angle_ml_response = zeros(size(rleg_angle_ml_normalized_all));
        
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            step_width_response(condition_indicator) = step_width_all(condition_indicator) - step_width_control_means(:, i_condition);
            foot_placement_world_response(condition_indicator) = foot_placement_world_all(condition_indicator) - foot_placement_world_control_means(:, i_condition);
            foot_placement_stancefoot_response(condition_indicator) = foot_placement_stancefoot_all(condition_indicator) - foot_placement_stancefoot_control_means(:, i_condition);
            foot_placement_mpsis_response(condition_indicator) = foot_placement_mpsis_all(condition_indicator) - foot_placement_mpsis_control_means(:, i_condition);
            foot_placement_com_response(condition_indicator) = foot_placement_com_all(condition_indicator) - foot_placement_com_control_means(:, i_condition);
            com_x_pos_response(:, condition_indicator) = com_x_pos_normalized_all(:, condition_indicator) - repmat(com_x_pos_control_means(:, i_condition), 1, sum(condition_indicator));
            lheel_x_pos_response(:, condition_indicator) = lheel_x_pos_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_control_means(:, i_condition), 1, sum(condition_indicator));
            rheel_x_pos_response(:, condition_indicator) = rheel_x_pos_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_control_means(:, i_condition), 1, sum(condition_indicator));
            lheel_y_pos_response(:, condition_indicator) = lheel_y_pos_normalized_all(:, condition_indicator) - repmat(lheel_y_pos_control_means(:, i_condition), 1, sum(condition_indicator));
            rheel_y_pos_response(:, condition_indicator) = rheel_y_pos_normalized_all(:, condition_indicator) - repmat(rheel_y_pos_control_means(:, i_condition), 1, sum(condition_indicator));
            lheel_x_pos_stancefoot_response(:, condition_indicator) = lheel_x_pos_stancefoot_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_stancefoot_control_means(:, i_condition), 1, sum(condition_indicator));
            rheel_x_pos_stancefoot_response(:, condition_indicator) = rheel_x_pos_stancefoot_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_stancefoot_control_means(:, i_condition), 1, sum(condition_indicator));
            lheel_x_pos_mpsis_response(:, condition_indicator) = lheel_x_pos_mpsis_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_mpsis_control_means(:, i_condition), 1, sum(condition_indicator));
            rheel_x_pos_mpsis_response(:, condition_indicator) = rheel_x_pos_mpsis_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_mpsis_control_means(:, i_condition), 1, sum(condition_indicator));
            lheel_x_pos_com_response(:, condition_indicator) = lheel_x_pos_com_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_com_control_means(:, i_condition), 1, sum(condition_indicator));
            rheel_x_pos_com_response(:, condition_indicator) = rheel_x_pos_com_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_com_control_means(:, i_condition), 1, sum(condition_indicator));       
            trunk_angle_ml_response(:, condition_indicator) = trunk_angle_ml_normalized_all(:, condition_indicator) - repmat(trunk_angle_ml_control_means(:, i_condition), 1, sum(condition_indicator));
            lleg_angle_ml_response(:, condition_indicator) = lleg_angle_ml_normalized_all(:, condition_indicator) - repmat(lleg_angle_ml_control_means(:, i_condition), 1, sum(condition_indicator));
            rleg_angle_ml_response(:, condition_indicator) = rleg_angle_ml_normalized_all(:, condition_indicator) - repmat(rleg_angle_ml_control_means(:, i_condition), 1, sum(condition_indicator));
        end
        for i_condition = 1 : number_of_conditions_to_analyze
            condition_indicator = conditions_to_analyze_indicators(:, i_condition);
            step_width_response(condition_indicator) = step_width_all(condition_indicator) - step_width_control_means(:, applicable_control_condition_indices(i_condition));
            foot_placement_world_response(condition_indicator) = foot_placement_world_all(condition_indicator) - foot_placement_world_control_means(:, applicable_control_condition_indices(i_condition));
            foot_placement_stancefoot_response(condition_indicator) = foot_placement_stancefoot_all(condition_indicator) - foot_placement_stancefoot_control_means(:, applicable_control_condition_indices(i_condition));
            foot_placement_mpsis_response(condition_indicator) = foot_placement_mpsis_all(condition_indicator) - foot_placement_mpsis_control_means(:, applicable_control_condition_indices(i_condition));
            foot_placement_com_response(condition_indicator) = foot_placement_com_all(condition_indicator) - foot_placement_com_control_means(:, applicable_control_condition_indices(i_condition));
            com_x_pos_response(:, condition_indicator) = com_x_pos_normalized_all(:, condition_indicator) - repmat(com_x_pos_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lheel_x_pos_response(:, condition_indicator) = lheel_x_pos_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rheel_x_pos_response(:, condition_indicator) = rheel_x_pos_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lheel_y_pos_response(:, condition_indicator) = lheel_y_pos_normalized_all(:, condition_indicator) - repmat(lheel_y_pos_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rheel_y_pos_response(:, condition_indicator) = rheel_y_pos_normalized_all(:, condition_indicator) - repmat(rheel_y_pos_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lheel_x_pos_stancefoot_response(:, condition_indicator) = lheel_x_pos_stancefoot_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_stancefoot_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rheel_x_pos_stancefoot_response(:, condition_indicator) = rheel_x_pos_stancefoot_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_stancefoot_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lheel_x_pos_mpsis_response(:, condition_indicator) = lheel_x_pos_mpsis_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_mpsis_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rheel_x_pos_mpsis_response(:, condition_indicator) = rheel_x_pos_mpsis_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_mpsis_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lheel_x_pos_com_response(:, condition_indicator) = lheel_x_pos_com_normalized_all(:, condition_indicator) - repmat(lheel_x_pos_com_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rheel_x_pos_com_response(:, condition_indicator) = rheel_x_pos_com_normalized_all(:, condition_indicator) - repmat(rheel_x_pos_com_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            trunk_angle_ml_response(:, condition_indicator) = trunk_angle_ml_normalized_all(:, condition_indicator) - repmat(trunk_angle_ml_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lleg_angle_ml_response(:, condition_indicator) = lleg_angle_ml_normalized_all(:, condition_indicator) - repmat(lleg_angle_ml_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rleg_angle_ml_response(:, condition_indicator) = rleg_angle_ml_normalized_all(:, condition_indicator) - repmat(rleg_angle_ml_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
        end
        
    end
    
    if process_data_forceplate
        % calculate control means
        cop_x_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        cop_x_stancefoot_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        cop_x_mpsis_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        cop_x_com_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        f_x_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        f_z_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        m_y_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        % removed these for now, as we're currently not using them. They're easy to add in later
%         lcop_x_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         rcop_x_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         fxl_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         fzl_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         myl_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         fxr_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         fzr_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
%         myr_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            cop_x_control_means(:, i_condition) = mean(cop_x_normalized_all(:, condition_indicator), 2);
            cop_x_stancefoot_control_means(:, i_condition) = mean(cop_x_stancefoot_normalized_all(:, condition_indicator), 2);
            cop_x_mpsis_control_means(:, i_condition) = mean(cop_x_mpsis_normalized_all(:, condition_indicator), 2);
            cop_x_com_control_means(:, i_condition) = mean(cop_x_com_normalized_all(:, condition_indicator), 2);
            f_x_control_means(:, i_condition) = mean(f_x_normalized_all(:, condition_indicator), 2);
            f_z_control_means(:, i_condition) = mean(f_z_normalized_all(:, condition_indicator), 2);
            m_y_control_means(:, i_condition) = mean(m_y_normalized_all(:, condition_indicator), 2);
%             lcop_x_control_means(:, i_condition) = mean(lcop_x_normalized_all(:, condition_indicator), 2);
%             rcop_x_control_means(:, i_condition) = mean(rcop_x_normalized_all(:, condition_indicator), 2);
%             fxl_control_means(:, i_condition) = mean(fxl_normalized_all(:, condition_indicator), 2);
%             fzl_control_means(:, i_condition) = mean(fzl_normalized_all(:, condition_indicator), 2);
%             myl_control_means(:, i_condition) = mean(myl_normalized_all(:, condition_indicator), 2);
%             fxr_control_means(:, i_condition) = mean(fxr_normalized_all(:, condition_indicator), 2);
%             fzr_control_means(:, i_condition) = mean(fzr_normalized_all(:, condition_indicator), 2);
%             myr_control_means(:, i_condition) = mean(myr_normalized_all(:, condition_indicator), 2);
        end
        
        % calculate stimulus responses
        cop_x_response = zeros(size(cop_x_normalized_all));
        cop_x_stancefoot_response = zeros(size(cop_x_normalized_all));
        cop_x_mpsis_response = zeros(size(cop_x_normalized_all));
        cop_x_com_response = zeros(size(cop_x_normalized_all));
        f_x_response = zeros(size(f_x_normalized_all));
        f_z_response = zeros(size(f_z_normalized_all));
        m_y_response = zeros(size(m_y_normalized_all));
%         lcop_x_response = zeros(size(lcop_x_normalized_all));
%         rcop_x_response = zeros(size(lcop_x_normalized_all));
%         fxl_response = zeros(size(fxl_normalized_all));
%         fzl_response = zeros(size(fzl_normalized_all));
%         myl_response = zeros(size(myl_normalized_all));
%         fxr_response = zeros(size(fxl_normalized_all));
%         fzr_response = zeros(size(fzl_normalized_all));
%         myr_response = zeros(size(myl_normalized_all));
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            cop_x_response(:, condition_indicator) = cop_x_normalized_all(:, condition_indicator) - repmat(cop_x_control_means(:, i_condition), 1, sum(condition_indicator));
            cop_x_stancefoot_response(:, condition_indicator) = cop_x_stancefoot_normalized_all(:, condition_indicator) - repmat(cop_x_stancefoot_control_means(:, i_condition), 1, sum(condition_indicator));
            cop_x_mpsis_response(:, condition_indicator) = cop_x_mpsis_normalized_all(:, condition_indicator) - repmat(cop_x_mpsis_control_means(:, i_condition), 1, sum(condition_indicator));
            cop_x_com_response(:, condition_indicator) = cop_x_com_normalized_all(:, condition_indicator) - repmat(cop_x_com_control_means(:, i_condition), 1, sum(condition_indicator));
            f_x_response(:, condition_indicator) = f_x_normalized_all(:, condition_indicator) - repmat(f_x_control_means(:, i_condition), 1, sum(condition_indicator));
            f_z_response(:, condition_indicator) = f_z_normalized_all(:, condition_indicator) - repmat(f_z_control_means(:, i_condition), 1, sum(condition_indicator));
            m_y_response(:, condition_indicator) = m_y_normalized_all(:, condition_indicator) - repmat(m_y_control_means(:, i_condition), 1, sum(condition_indicator));
%             fxl_response(:, condition_indicator) = fxl_normalized_all(:, condition_indicator) - repmat(fxl_control_means(:, i_condition), 1, sum(condition_indicator));
%             fzl_response(:, condition_indicator) = fzl_normalized_all(:, condition_indicator) - repmat(fzl_control_means(:, i_condition), 1, sum(condition_indicator));
%             myl_response(:, condition_indicator) = myl_normalized_all(:, condition_indicator) - repmat(myl_control_means(:, i_condition), 1, sum(condition_indicator));
%             lcop_x_response(:, condition_indicator) = lcop_x_normalized_all(:, condition_indicator) - repmat(lcop_x_control_means(:, i_condition), 1, sum(condition_indicator));
%             fxr_response(:, condition_indicator) = fxr_normalized_all(:, condition_indicator) - repmat(fxr_control_means(:, i_condition), 1, sum(condition_indicator));
%             fzr_response(:, condition_indicator) = fzr_normalized_all(:, condition_indicator) - repmat(fzr_control_means(:, i_condition), 1, sum(condition_indicator));
%             myr_response(:, condition_indicator) = myr_normalized_all(:, condition_indicator) - repmat(myr_control_means(:, i_condition), 1, sum(condition_indicator));
%             rcop_x_response(:, condition_indicator) = rcop_x_normalized_all(:, condition_indicator) - repmat(rcop_x_control_means(:, i_condition), 1, sum(condition_indicator));
        end
        for i_condition = 1 : number_of_conditions_to_analyze
            condition_indicator = conditions_to_analyze_indicators(:, i_condition);
            cop_x_response(:, condition_indicator) = cop_x_normalized_all(:, condition_indicator) - repmat(cop_x_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            cop_x_stancefoot_response(:, condition_indicator) = cop_x_stancefoot_normalized_all(:, condition_indicator) - repmat(cop_x_stancefoot_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            cop_x_mpsis_response(:, condition_indicator) = cop_x_mpsis_normalized_all(:, condition_indicator) - repmat(cop_x_mpsis_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            cop_x_com_response(:, condition_indicator) = cop_x_com_normalized_all(:, condition_indicator) - repmat(cop_x_com_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            f_x_response(:, condition_indicator) = f_x_normalized_all(:, condition_indicator) - repmat(f_x_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            f_z_response(:, condition_indicator) = f_z_normalized_all(:, condition_indicator) - repmat(f_z_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            m_y_response(:, condition_indicator) = m_y_normalized_all(:, condition_indicator) - repmat(m_y_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             fxl_response(:, condition_indicator) = fxl_normalized_all(:, condition_indicator) - repmat(fxl_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             fzl_response(:, condition_indicator) = fzl_normalized_all(:, condition_indicator) - repmat(fzl_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             myl_response(:, condition_indicator) = myl_normalized_all(:, condition_indicator) - repmat(myl_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             lcop_x_response(:, condition_indicator) = lcop_x_normalized_all(:, condition_indicator) - repmat(lcop_x_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             fxr_response(:, condition_indicator) = fxr_normalized_all(:, condition_indicator) - repmat(fxr_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             fzr_response(:, condition_indicator) = fzr_normalized_all(:, condition_indicator) - repmat(fzr_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             myr_response(:, condition_indicator) = myr_normalized_all(:, condition_indicator) - repmat(myr_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             rcop_x_response(:, condition_indicator) = rcop_x_normalized_all(:, condition_indicator) - repmat(rcop_x_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
        end
        
    end
    
    if process_data_emg
        % calculate control means
        lglutmed_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        ltibiant_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lgastroc_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        lperolng_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rglutmed_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rtibiant_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rgastroc_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rperolng_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        ltnsrflt_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        rtnsrflt_control_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        
        if normalize_emg_mean || normalize_emg_peak
            lglutmed_control_rescaled_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
            ltibiant_control_rescaled_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
            lgastroc_control_rescaled_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
            lperolng_control_rescaled_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
            rglutmed_control_rescaled_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
            rtibiant_control_rescaled_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
            rgastroc_control_rescaled_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
            rperolng_control_rescaled_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
            ltnsrflt_control_rescaled_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
            rtnsrflt_control_rescaled_means = zeros(number_of_time_steps_normalized, number_of_conditions_control);
        end
        
        for i_condition = 1 : number_of_conditions_control
            condition_indicator = conditions_control_indicators(:, i_condition);
            lglutmed_control_means(:, i_condition) = mean(lglutmed_normalized_all(:, condition_indicator), 2);
            ltibiant_control_means(:, i_condition) = mean(ltibiant_normalized_all(:, condition_indicator), 2);
            lgastroc_control_means(:, i_condition) = mean(lgastroc_normalized_all(:, condition_indicator), 2);
            lperolng_control_means(:, i_condition) = mean(lperolng_normalized_all(:, condition_indicator), 2);
            rglutmed_control_means(:, i_condition) = mean(rglutmed_normalized_all(:, condition_indicator), 2);
            rtibiant_control_means(:, i_condition) = mean(rtibiant_normalized_all(:, condition_indicator), 2);
            rgastroc_control_means(:, i_condition) = mean(rgastroc_normalized_all(:, condition_indicator), 2);
            rperolng_control_means(:, i_condition) = mean(rperolng_normalized_all(:, condition_indicator), 2);
            ltnsrflt_control_means(:, i_condition) = mean(ltnsrflt_normalized_all(:, condition_indicator), 2);
            rtnsrflt_control_means(:, i_condition) = mean(rtnsrflt_normalized_all(:, condition_indicator), 2);
            if normalize_emg_mean || normalize_emg_peak
                condition_indicator = conditions_control_indicators(:, i_condition);
                lglutmed_control_rescaled_means(:, i_condition) = mean(lglutmed_normalized_rescaled_all(:, condition_indicator), 2);
                ltibiant_control_rescaled_means(:, i_condition) = mean(ltibiant_normalized_rescaled_all(:, condition_indicator), 2);
                lgastroc_control_rescaled_means(:, i_condition) = mean(lgastroc_normalized_rescaled_all(:, condition_indicator), 2);
                lperolng_control_rescaled_means(:, i_condition) = mean(lperolng_normalized_rescaled_all(:, condition_indicator), 2);
                rglutmed_control_rescaled_means(:, i_condition) = mean(rglutmed_normalized_rescaled_all(:, condition_indicator), 2);
                rtibiant_control_rescaled_means(:, i_condition) = mean(rtibiant_normalized_rescaled_all(:, condition_indicator), 2);
                rgastroc_control_rescaled_means(:, i_condition) = mean(rgastroc_normalized_rescaled_all(:, condition_indicator), 2);
                rperolng_control_rescaled_means(:, i_condition) = mean(rperolng_normalized_rescaled_all(:, condition_indicator), 2);
                ltnsrflt_control_rescaled_means(:, i_condition) = mean(ltnsrflt_normalized_rescaled_all(:, condition_indicator), 2);
                rtnsrflt_control_rescaled_means(:, i_condition) = mean(rtnsrflt_normalized_rescaled_all(:, condition_indicator), 2);
            end
        end
        
        % calculate stimulus responses
        lglutmed_response = zeros(size(lglutmed_normalized_all));
        ltibiant_response = zeros(size(ltibiant_normalized_all));
        lgastroc_response = zeros(size(lgastroc_normalized_all));
        lperolng_response = zeros(size(lperolng_normalized_all));
        rglutmed_response = zeros(size(rglutmed_normalized_all));
        rtibiant_response = zeros(size(rtibiant_normalized_all));
        rgastroc_response = zeros(size(rgastroc_normalized_all));
        rperolng_response = zeros(size(rperolng_normalized_all));
        ltnsrflt_response = zeros(size(ltnsrflt_normalized_all));
        rtnsrflt_response = zeros(size(rtnsrflt_normalized_all));
        if normalize_emg_mean || normalize_emg_peak
            lglutmed_rescaled_response = zeros(size(lglutmed_normalized_all));
            ltibiant_rescaled_response = zeros(size(ltibiant_normalized_all));
            lgastroc_rescaled_response = zeros(size(lgastroc_normalized_all));
            lperolng_rescaled_response = zeros(size(lperolng_normalized_all));
            rglutmed_rescaled_response = zeros(size(rglutmed_normalized_all));
            rtibiant_rescaled_response = zeros(size(rtibiant_normalized_all));
            rgastroc_rescaled_response = zeros(size(rgastroc_normalized_all));
            rperolng_rescaled_response = zeros(size(rperolng_normalized_all));
            ltnsrflt_rescaled_response = zeros(size(ltnsrflt_normalized_all));
            rtnsrflt_rescaled_response = zeros(size(rtnsrflt_normalized_all));
        end
        
        for i_condition = 1 : number_of_conditions_to_analyze
            condition_indicator = conditions_to_analyze_indicators(:, i_condition);
            lglutmed_response(:, condition_indicator) = lglutmed_normalized_all(:, condition_indicator) - repmat(lglutmed_control_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            ltibiant_response(:, condition_indicator) = ltibiant_normalized_all(:, condition_indicator) - repmat(ltibiant_control_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lgastroc_response(:, condition_indicator) = lgastroc_normalized_all(:, condition_indicator) - repmat(lgastroc_control_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            lperolng_response(:, condition_indicator) = lperolng_normalized_all(:, condition_indicator) - repmat(lperolng_control_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rglutmed_response(:, condition_indicator) = rglutmed_normalized_all(:, condition_indicator) - repmat(rglutmed_control_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rtibiant_response(:, condition_indicator) = rtibiant_normalized_all(:, condition_indicator) - repmat(rtibiant_control_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rgastroc_response(:, condition_indicator) = rgastroc_normalized_all(:, condition_indicator) - repmat(rgastroc_control_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rperolng_response(:, condition_indicator) = rperolng_normalized_all(:, condition_indicator) - repmat(rperolng_control_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            ltnsrflt_response(:, condition_indicator) = ltnsrflt_normalized_all(:, condition_indicator) - repmat(ltnsrflt_control_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            rtnsrflt_response(:, condition_indicator) = rtnsrflt_normalized_all(:, condition_indicator) - repmat(rtnsrflt_control_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            if normalize_emg_mean || normalize_emg_peak   
                 condition_indicator = conditions_to_analyze_indicators(:, i_condition);
                 lglutmed_rescaled_response(:, condition_indicator) = lglutmed_normalized_rescaled_all(:, condition_indicator) - repmat(lglutmed_control_rescaled_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
                 ltibiant_rescaled_response(:, condition_indicator) = ltibiant_normalized_rescaled_all(:, condition_indicator) - repmat(ltibiant_control_rescaled_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
                 lgastroc_rescaled_response(:, condition_indicator) = lgastroc_normalized_rescaled_all(:, condition_indicator) - repmat(lgastroc_control_rescaled_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
                 lperolng_rescaled_response(:, condition_indicator) = lperolng_normalized_rescaled_all(:, condition_indicator) - repmat(lperolng_control_rescaled_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
                 rglutmed_rescaled_response(:, condition_indicator) = rglutmed_normalized_rescaled_all(:, condition_indicator) - repmat(rglutmed_control_rescaled_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
                 rtibiant_rescaled_response(:, condition_indicator) = rtibiant_normalized_rescaled_all(:, condition_indicator) - repmat(rtibiant_control_rescaled_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
                 rgastroc_rescaled_response(:, condition_indicator) = rgastroc_normalized_rescaled_all(:, condition_indicator) - repmat(rgastroc_control_rescaled_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
                 rperolng_rescaled_response(:, condition_indicator) = rperolng_normalized_rescaled_all(:, condition_indicator) - repmat(rperolng_control_rescaled_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
                 ltnsrflt_rescaled_response(:, condition_indicator) = ltnsrflt_normalized_rescaled_all(:, condition_indicator) - repmat(ltnsrflt_control_rescaled_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
                 rtnsrflt_rescaled_response(:, condition_indicator) = rtnsrflt_normalized_rescaled_all(:, condition_indicator) - repmat(rtnsrflt_control_rescaled_means(:,applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
            end
        end
    end
     
%         for i_condition = 1 : number_of_conditions_to_analyze
%             condition_indicator = conditions_to_analyze_indicators(:, i_condition);
%             lglutmed_response(:, condition_indicator) = lglutmed_normalized_all(:, condition_indicator) - repmat(lglutmed_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             ltibiant_response(:, condition_indicator) = ltibiant_normalized_all(:, condition_indicator) - repmat(ltibiant_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             lgastroc_response(:, condition_indicator) = lgastroc_normalized_all(:, condition_indicator) - repmat(lgastroc_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             lperolng_response(:, condition_indicator) = lperolng_normalized_all(:, condition_indicator) - repmat(lperolng_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             rglutmed_response(:, condition_indicator) = rglutmed_normalized_all(:, condition_indicator) - repmat(rglutmed_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             rtibiant_response(:, condition_indicator) = rtibiant_normalized_all(:, condition_indicator) - repmat(rtibiant_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             rgastroc_response(:, condition_indicator) = rgastroc_normalized_all(:, condition_indicator) - repmat(rgastroc_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%             rperolng_response(:, condition_indicator) = rperolng_normalized_all(:, condition_indicator) - repmat(rperolng_control_means(:, applicable_control_condition_indices(i_condition)), 1, sum(condition_indicator));
%         end    

%     end
    
    
    %% save data
    if save_data

        save ...
          ( ...
            ['analysis' filesep makeFileName(date, subject_id, 'resultsConditions')], ...
            'number_of_time_steps_normalized', ...
            'origin_trial_list_all', ...
            'origin_start_time_list_all', ...
            'origin_end_time_list_all', ...
            'conditions_control', ...
            'conditions_to_analyze', ...
            'condition_stance_foot_list_all', ...
            'condition_perturbation_list_all', ...
            'condition_delay_list_all', ...
            'condition_index_list_all', ...
            'condition_experimental_list_all', ...
            'conditions_control_indicators', ...
            'conditions_to_analyze_indicators', ...
            'step_times_all', ...
            'pushoff_times_all', ...
            'cadence_all' ...
          );

        if process_data_balance
            save ...
              ( ...
                ['analysis' filesep makeFileName(date, subject_id, 'resultsBalance')], ...
                'step_width_all', ...
                'step_length_all', ...
                'step_speed_all', ...
                'foot_placement_world_all', ...
                'foot_placement_stancefoot_all', ...
                'foot_placement_mpsis_all', ...
                'foot_placement_com_all', ...
                'com_x_pos_normalized_all', ...
                'lheel_x_pos_normalized_all', ...
                'rheel_x_pos_normalized_all', ...
                'lheel_y_pos_normalized_all', ...
                'rheel_y_pos_normalized_all', ...
                'lheel_x_pos_stancefoot_normalized_all', ...
                'rheel_x_pos_stancefoot_normalized_all', ...
                'lheel_x_pos_mpsis_normalized_all', ...
                'rheel_x_pos_mpsis_normalized_all', ...
                'lheel_x_pos_com_normalized_all',...
                'rheel_x_pos_com_normalized_all',...
                'trunk_angle_ap_normalized_all', ...
                'trunk_angle_ml_normalized_all', ...
                'lleg_angle_ml_normalized_all', ...
                'rleg_angle_ml_normalized_all', ...
                'step_width_response', ...
                'com_x_pos_response', ...
                'foot_placement_world_response', ...
                'foot_placement_stancefoot_response', ...
                'foot_placement_mpsis_response', ...
                'foot_placement_com_response', ...
                'lheel_x_pos_response', ...
                'rheel_x_pos_response', ...
                'lheel_y_pos_response', ...
                'rheel_y_pos_response', ...
                'lheel_x_pos_stancefoot_response', ...
                'rheel_x_pos_stancefoot_response', ...
                'lheel_x_pos_mpsis_response', ...
                'rheel_x_pos_mpsis_response', ...
                'lheel_x_pos_com_response', ...
                'rheel_x_pos_com_response', ...
                'trunk_angle_ml_response', ...
                'lleg_angle_ml_response', ...
                'rleg_angle_ml_response' ...
              );
        end

        if process_data_armswing
            save ...
              ( ...
                ['analysis' filesep makeFileName(date, subject_id, 'resultsArmswing')], ...
                'linclination_normalized_all', ...
                'rinclination_normalized_all' ...
              );
        end

        if process_data_forceplate
            save ...
              ( ...
                ['analysis' filesep makeFileName(date, subject_id, 'resultsForceplate')], ...
                'cop_x_normalized_all', ...
                'cop_x_response', ...
                'cop_x_stancefoot_normalized_all', ...
                'cop_x_stancefoot_response', ...
                'cop_x_mpsis_normalized_all', ...
                'cop_x_mpsis_response', ...
                'cop_x_com_normalized_all', ...
                'cop_x_com_response', ...
                'f_x_normalized_all', ...
                'f_x_response', ...
                'f_z_normalized_all', ...
                'f_z_response', ...
                'm_y_normalized_all', ...
                'm_y_response' ...
              );
        end
        
        if process_data_emg
            save ...
              ( ...
                ['analysis' filesep makeFileName(date, subject_id, 'resultsEmg')], ...
                'lglutmed_normalized_all', ...
                'ltibiant_normalized_all', ...
                'lgastroc_normalized_all', ...
                'lperolng_normalized_all', ...
                'rglutmed_normalized_all', ...
                'rtibiant_normalized_all', ...
                'rgastroc_normalized_all', ...
                'rperolng_normalized_all', ...
                'ltnsrflt_normalized_all', ...
                'rtnsrflt_normalized_all', ...
                'lglutmed_normalized_rescaled_all', ...
                'ltibiant_normalized_rescaled_all', ...
                'lgastroc_normalized_rescaled_all', ...
                'lperolng_normalized_rescaled_all', ...
                'rglutmed_normalized_rescaled_all', ...
                'rtibiant_normalized_rescaled_all', ...
                'rgastroc_normalized_rescaled_all', ...
                'rperolng_normalized_rescaled_all', ...
                'ltnsrflt_normalized_rescaled_all', ...
                'rtnsrflt_normalized_rescaled_all', ...
                'lglutmed_response', ...
                'ltibiant_response', ...
                'lgastroc_response', ...
                'lperolng_response', ...
                'rglutmed_response', ...
                'rtibiant_response', ...
                'rgastroc_response', ...
                'rperolng_response', ...
                'ltnsrflt_response', ...
                'rtnsrflt_response', ...
                'lglutmed_rescaled_response', ...
                'ltibiant_rescaled_response', ...
                'lgastroc_rescaled_response', ...
                'lperolng_rescaled_response', ...
                'rglutmed_rescaled_response', ...
                'rtibiant_rescaled_response', ...
                'rgastroc_rescaled_response', ...
                'rperolng_rescaled_response', ...
                'ltnsrflt_rescaled_response', ...
                'rtnsrflt_rescaled_response' ...
                )
        end    
    end
    
    

