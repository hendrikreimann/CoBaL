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

% input

% description

% there are different levels of variables that this object takes care of
% - basic variables can be loaded or calculated directly from loaded variables at each point in time
% - stretch variables are variables defined for each stretch. That can be a single number or a trajectory with
% normalized time. 

% order of processing is
% Construction:
% - variables are determined based on the list of variables_to_analyze, loaded from the study_settings
% For each trial
% 1. Prepare data: loads and calculates basic variables
% 2. <apply stretches and normalize in time>
% 3. calculate stretch variables

% to add a new variable, you need to do the following things:
%
% 1. Provide information about all variables this new one depends upon. Go to the determineVariables function and add a
% new conditional statement "if this.isVariableToAnalyze('<your new variable name>') ... end". Within the
% statement, add the basic variables and stretch variables that the new variable depends upon by calling addBasicVariable 
% and addStretchVariable. The order matters here, make sure for each variable, all dependencies have already been added, 
% i.e. don't call "this.addBasicVariable('lheel_y')" before "this.addBasicVariable('marker_trajectories')", because 
% lheel_y depends upon the marker_trajectories.
%
% 2. If you added a new basic variable, but is not already available, i.e. one of the outputs of the
% preprocess scripts, you need to specify how to calculate it. Go to prepareBasicVariables and add a new conditional
% statement "if strcmp(variable_name, '<your new variable name>')". Within the statement, calculate your new
% variable, then add it to the variable_data. Also add the corresponding time vector to the time_data.
%
% 3. If the variable you want to add is a stretch variable, you need to specify how to calculate it. Go to 
% calculateStretchVariables and add a new conditional statement "if strcmp(variable_name, '<your new variable name>')".
% Within the statement, calculate the value of your new variable for a single stretch, depending upon
% this_stretch_start_time and this_stretch_end_time. Assign that value to the variable "stretch_data".


classdef WalkingDataCustodian < handle
    properties
        data_directory = [];
        date = [];
        subject_id = [];
        trial_type = [];
        trial_number = [];
        variables_to_analyze = {};
        basic_variable_names = {};
        basic_variable_load_failures = {};
        stretch_variable_names = {};
        stretch_variable_directions = {};
        
        % these are private and should only be accessed using get functions
        number_of_time_steps_normalized;
        basic_variable_data;
        basic_variable_labels;
        basic_variable_directions;
        stretch_variable_data;
        time_data;
        
        subject_info;
        subject_settings;
        study_settings;
        emg_normalization_values;
        emg_normalization_labels;
    end
    
    methods
        % constructor
        function this = WalkingDataCustodian(variables_to_analyze, data_directory)
            if nargin < 2
                data_directory = pwd;
            end
            this.subject_info = load([data_directory filesep 'subjectInfo.mat']);
            this.study_settings = loadSettingsFromFile('study', data_directory);
            this.subject_settings = loadSettingsFromFile('subject', data_directory);
            
            % load this information from the subjects.mat and studySettings.txt files
            this.data_directory = data_directory;
            this.date = this.subject_settings.get('collection_date');
            this.subject_id = this.subject_settings.get('subject_id');

            % prepare things for EMG normalization
            emg_normalization_file_name = [data_directory filesep 'analysis' filesep makeFileName(this.date, this.subject_id, 'emgNormalization.mat')];
            if exist(emg_normalization_file_name, 'file')
                emg_normalization_data = load(emg_normalization_file_name);
                this.emg_normalization_values = emg_normalization_data.emg_normalization_values;
                this.emg_normalization_labels = emg_normalization_data.emg_variable_names;
            end
            
            % prepare list of variables to analyze
            if nargin < 1
                variables_to_analyze = this.study_settings.get('stretch_variables');
            end
            this.variables_to_analyze = variables_to_analyze;
            this.number_of_time_steps_normalized = this.study_settings.get('number_of_time_steps_normalized');
            this.determineVariables();
        end
        
        % initialization
        function addBasicVariable(this, variable_name)
            if ~any(strcmp(this.basic_variable_names, variable_name))
                this.basic_variable_names = [this.basic_variable_names; variable_name];
            end
        end
        function addStretchVariable(this, variable_name)
            if ~any(strcmp(this.stretch_variable_names, variable_name))
                this.stretch_variable_names = [this.stretch_variable_names; variable_name];
                this.stretch_variable_directions = [this.stretch_variable_directions; {'TBD', 'TBD'}];
            end
        end
        function determineVariables(this)
            
            % go through list of variables to analyze and check if they are elementary variables
            for i_variable = 1 : length(this.variables_to_analyze)
                this_variable_name = this.variables_to_analyze{i_variable};
                
                % check if this is a compound name, listing a loaded variable and a label
                if any(this_variable_name==':')
                    this_variable_split = strsplit(this_variable_name, ':');
                    this_variable_type = this_variable_split{1};
                    this_variable_label = this_variable_split{2};
                    
                    if strcmp(this_variable_type, 'elementary') || strcmp(this_variable_type, 'derivative')
                        % not a true compound, but label is elementary variable
                        this.addBasicVariable(this_variable_label)
                        this.addStretchVariable(this_variable_label)
                    else
                        % this is a true compound
                        this.addBasicVariable([this_variable_type '_trajectories'])
                        this.addStretchVariable(this_variable_name)
                    end
                end

                % use an initial value as a reference
                if strcmp(this_variable_name, 'lheel_from_mpsis_initial_x')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('lheel_from_mpsis_initial_x')                
                end
                if strcmp(this_variable_name, 'rheel_from_mpsis_initial_x')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('rheel_from_mpsis_initial_x')
                end
                if strcmp(this_variable_name, 'mpsis_from_mpsis_initial_x')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('mpsis_from_mpsis_initial_x')
                end            
                if strcmp(this_variable_name, 'lheel_from_com_initial_x')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('com_position_trajectories')
                    this.addStretchVariable('lheel_from_com_initial_x')                
                end
                if strcmp(this_variable_name, 'rheel_from_com_initial_x')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('com_position_trajectories')
                    this.addStretchVariable('rheel_from_com_initial_x')
                end
                if strcmp(this_variable_name, 'com_from_com_initial_x')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('com_position_trajectories')
                    this.addStretchVariable('com_from_com_initial_x')
                end

                % balance stuff
                if strcmp(this_variable_name, 'step_placement_x')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('step_placement_x')
                end
                if strcmp(this_variable_name, 'xcom_x')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('com_position_trajectories')
                    this.addBasicVariable('com_velocity_trajectories')
                    this.addStretchVariable('xcom_x')
                end
                if strcmp(this_variable_name, 'xcom_y')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('com_position_trajectories')
                    this.addBasicVariable('com_velocity_trajectories')
                    this.addStretchVariable('xcom_y')
                end
                if strcmp(this_variable_name, 'mos_x')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('com_position_trajectories')
                    this.addBasicVariable('com_velocity_trajectories')
                    this.addStretchVariable('xcom_x')
                    this.addStretchVariable('mos_x')
                end
                if strcmp(this_variable_name, 'mos_y')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('com_position_trajectories')
                    this.addBasicVariable('com_velocity_trajectories')
                    this.addStretchVariable('xcom_y')
                    this.addStretchVariable('mos_y')
                end
                if strcmp(this_variable_name, 'xcom_mpsis_x')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('derivative:LPSI_x_vel')
                    this.addBasicVariable('derivative:RPSI_x_vel')
                    this.addStretchVariable('xcom_mpsis_x')
                end
                if strcmp(this_variable_name, 'xcom_mpsis_y')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('derivative:LPSI_y_vel')
                    this.addBasicVariable('derivative:RPSI_y_vel')
                    this.addStretchVariable('xcom_mpsis_y')
                end
                if strcmp(this_variable_name, 'mos_mpsis_x')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('derivative:LPSI_x_vel')
                    this.addBasicVariable('derivative:RPSI_x_vel')
                    this.addStretchVariable('xcom_mpsis_x')
                    this.addStretchVariable('mos_mpsis_x')
                end
                if strcmp(this_variable_name, 'mos_mpsis_y')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('derivative:LPSI_y_vel')
                    this.addBasicVariable('derivative:RPSI_y_vel')
                    this.addStretchVariable('xcom_mpsis_y')
                    this.addStretchVariable('mos_mpsis_y')
                end
                if strcmp(this_variable_name, 'com_rough_x')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('com_rough_x')
                end
                if strcmp(this_variable_name, 'com_rough_y')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('com_rough_y')
                end  
                if strcmp(this_variable_name, 'com_rough_z')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('com_rough_z')
                end
                if strcmp(this_variable_name, 'com_rough_x_vel')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('derivative:LPSI_x_vel')
                    this.addBasicVariable('derivative:RPSI_x_vel')
                    this.addBasicVariable('derivative:LASI_x_vel')
                    this.addBasicVariable('derivative:RASI_x_vel')
                    this.addBasicVariable('derivative:LKNE_x_vel')
                    this.addBasicVariable('derivative:RKNE_x_vel')
                    this.addBasicVariable('derivative:LANK_x_vel')
                    this.addBasicVariable('derivative:RANK_x_vel')
                    this.addBasicVariable('derivative:LTOE_x_vel')
                    this.addBasicVariable('derivative:RTOE_x_vel')
                    this.addBasicVariable('derivative:LFHD_x_vel')
                    this.addBasicVariable('derivative:RFHD_x_vel')
                    this.addBasicVariable('derivative:LBHD_x_vel')
                    this.addBasicVariable('derivative:RBHD_x_vel')
                    this.addBasicVariable('derivative:C7_x_vel')
                    this.addStretchVariable('com_rough_x_vel')
                end  
                if strcmp(this_variable_name, 'com_rough_y_vel')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('derivative:LPSI_y_vel')
                    this.addBasicVariable('derivative:RPSI_y_vel')
                    this.addBasicVariable('derivative:LASI_y_vel')
                    this.addBasicVariable('derivative:RASI_y_vel')
                    this.addBasicVariable('derivative:LKNE_y_vel')
                    this.addBasicVariable('derivative:RKNE_y_vel')
                    this.addBasicVariable('derivative:LANK_y_vel')
                    this.addBasicVariable('derivative:RANK_y_vel')
                    this.addBasicVariable('derivative:LTOE_y_vel')
                    this.addBasicVariable('derivative:RTOE_y_vel')
                    this.addBasicVariable('derivative:LFHD_y_vel')
                    this.addBasicVariable('derivative:RFHD_y_vel')
                    this.addBasicVariable('derivative:LBHD_y_vel')
                    this.addBasicVariable('derivative:RBHD_y_vel')
                    this.addBasicVariable('derivative:C7_y_vel')
                    this.addStretchVariable('com_rough_y_vel')
                end 
                if strcmp(this_variable_name, 'xcom_rough_x')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('derivative:LPSI_x_vel')
                    this.addBasicVariable('derivative:RPSI_x_vel')
                    this.addBasicVariable('derivative:LASI_x_vel')
                    this.addBasicVariable('derivative:RASI_x_vel')
                    this.addBasicVariable('derivative:LKNE_x_vel')
                    this.addBasicVariable('derivative:RKNE_x_vel')
                    this.addBasicVariable('derivative:LANK_x_vel')
                    this.addBasicVariable('derivative:RANK_x_vel')
                    this.addBasicVariable('derivative:LTOE_x_vel')
                    this.addBasicVariable('derivative:RTOE_x_vel')
                    this.addBasicVariable('derivative:LFHD_x_vel')
                    this.addBasicVariable('derivative:RFHD_x_vel')
                    this.addBasicVariable('derivative:LBHD_x_vel')
                    this.addBasicVariable('derivative:RBHD_x_vel')
                    this.addBasicVariable('derivative:C7_x_vel')
                    this.addStretchVariable('com_rough_x')
                    this.addStretchVariable('com_rough_y')
                    this.addStretchVariable('com_rough_z')
                    this.addStretchVariable('com_rough_x_vel')
                    this.addStretchVariable('xcom_rough_x')
                end
                if strcmp(this_variable_name, 'xcom_rough_y')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('derivative:LPSI_y_vel')
                    this.addBasicVariable('derivative:RPSI_y_vel')
                    this.addBasicVariable('derivative:LASI_y_vel')
                    this.addBasicVariable('derivative:RASI_y_vel')
                    this.addBasicVariable('derivative:LKNE_y_vel')
                    this.addBasicVariable('derivative:RKNE_y_vel')
                    this.addBasicVariable('derivative:LANK_y_vel')
                    this.addBasicVariable('derivative:RANK_y_vel')
                    this.addBasicVariable('derivative:LTOE_y_vel')
                    this.addBasicVariable('derivative:RTOE_y_vel')
                    this.addBasicVariable('derivative:LFHD_y_vel')
                    this.addBasicVariable('derivative:RFHD_y_vel')
                    this.addBasicVariable('derivative:LBHD_y_vel')
                    this.addBasicVariable('derivative:RBHD_y_vel')
                    this.addBasicVariable('derivative:C7_y_vel')
                    this.addStretchVariable('com_rough_x')
                    this.addStretchVariable('com_rough_y')
                    this.addStretchVariable('com_rough_z')
                    this.addStretchVariable('com_rough_y_vel')
                    this.addStretchVariable('xcom_rough_y')
                end
                if strcmp(this_variable_name, 'mos_rough_x')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('derivative:LPSI_x_vel')
                    this.addBasicVariable('derivative:RPSI_x_vel')
                    this.addBasicVariable('derivative:LASI_x_vel')
                    this.addBasicVariable('derivative:RASI_x_vel')
                    this.addBasicVariable('derivative:LKNE_x_vel')
                    this.addBasicVariable('derivative:RKNE_x_vel')
                    this.addBasicVariable('derivative:LANK_x_vel')
                    this.addBasicVariable('derivative:RANK_x_vel')
                    this.addBasicVariable('derivative:LTOE_x_vel')
                    this.addBasicVariable('derivative:RTOE_x_vel')
                    this.addBasicVariable('derivative:LFHD_x_vel')
                    this.addBasicVariable('derivative:RFHD_x_vel')
                    this.addBasicVariable('derivative:LBHD_x_vel')
                    this.addBasicVariable('derivative:RBHD_x_vel')
                    this.addBasicVariable('derivative:C7_x_vel')
                    this.addStretchVariable('com_rough_x')
                    this.addStretchVariable('com_rough_y')
                    this.addStretchVariable('com_rough_z')
                    this.addStretchVariable('com_rough_x_vel')
                    this.addStretchVariable('com_rough_x_vel')
                    this.addStretchVariable('xcom_rough_x')
                    this.addStretchVariable('mos_rough_x')
                end
                if strcmp(this_variable_name, 'mos_rough_y')
                    this.addBasicVariable('marker_trajectories')
                    this.addBasicVariable('derivative:LPSI_y_vel')
                    this.addBasicVariable('derivative:RPSI_y_vel')
                    this.addBasicVariable('derivative:LASI_y_vel')
                    this.addBasicVariable('derivative:RASI_y_vel')
                    this.addBasicVariable('derivative:LKNE_y_vel')
                    this.addBasicVariable('derivative:RKNE_y_vel')
                    this.addBasicVariable('derivative:LANK_y_vel')
                    this.addBasicVariable('derivative:RANK_y_vel')
                    this.addBasicVariable('derivative:LTOE_y_vel')
                    this.addBasicVariable('derivative:RTOE_y_vel')
                    this.addBasicVariable('derivative:LFHD_y_vel')
                    this.addBasicVariable('derivative:RFHD_y_vel')
                    this.addBasicVariable('derivative:LBHD_y_vel')
                    this.addBasicVariable('derivative:RBHD_y_vel')
                    this.addBasicVariable('derivative:C7_y_vel')
                    this.addStretchVariable('com_rough_x')
                    this.addStretchVariable('com_rough_y')
                    this.addStretchVariable('com_rough_z')
                    this.addStretchVariable('com_rough_y_vel')
                    this.addStretchVariable('xcom_rough_y')
                    this.addStretchVariable('mos_rough_y')
                end

                % gait parameters
                if strcmp(this_variable_name, 'step_length')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('step_length')
                end
                if strcmp(this_variable_name, 'step_width')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('step_width')
                end
                if strcmp(this_variable_name, 'step_time')
                    this.addStretchVariable('step_time')
                end
                if strcmp(this_variable_name, 'pushoff_time')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('step_time')
                    this.addStretchVariable('pushoff_time')
                end
                if strcmp(this_variable_name, 'midstance_index')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('step_time')
                    this.addStretchVariable('pushoff_time')
                    this.addStretchVariable('midstance_index')
                end
                if strcmp(this_variable_name, 'midstance_time')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('step_time')
                    this.addStretchVariable('pushoff_time')
                    this.addStretchVariable('midstance_time')
                end
                if strcmp(this_variable_name, 'midswing_time')
                    this.addStretchVariable('step_time')
                    this.addStretchVariable('pushoff_time')
                    this.addStretchVariable('midswing_time')
                end
                if strcmp(this_variable_name, 'cadence')
                    this.addStretchVariable('step_time')
                    this.addStretchVariable('cadence')
                end
                if strcmp(this_variable_name, 'velocity')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('step_length')
                    this.addStretchVariable('step_time')
                    this.addStretchVariable('velocity')
                end
                if strcmp(this_variable_name, 'band_duration')
                    this.addStretchVariable('band_duration')
                end

                % segment angles and lengths
                if strcmp(this_variable_name, 'leg_length_l')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('leg_length_l')
                end
                if strcmp(this_variable_name, 'leg_length_r')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('leg_length_r')
                end
                if strcmp(this_variable_name, 'leg_angle_ap_l')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('leg_angle_ap_l')
                end
                if strcmp(this_variable_name, 'leg_angle_ap_r')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('leg_angle_ap_r')
                end
                if strcmp(this_variable_name, 'head_angle_ap')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('head_angle_ap')
                end
                if strcmp(this_variable_name, 'head_angle_ml')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('head_angle_ml')
                end
                if strcmp(this_variable_name, 'trunk_angle_ap')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('trunk_angle_ap')
                end
                if strcmp(this_variable_name, 'trunk_angle_ml')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('trunk_angle_ml')
                end
                if strcmp(this_variable_name, 'left_foot_angle_ap')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('left_foot_angle_ap')
                end
                if strcmp(this_variable_name, 'left_foot_angle_ml')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('left_foot_angle_ml')
                end
                if strcmp(this_variable_name, 'left_foot_angle_yaw')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('left_foot_angle_yaw')
                end
                if strcmp(this_variable_name, 'right_foot_angle_ap')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('right_foot_angle_ap')
                end
                if strcmp(this_variable_name, 'right_foot_angle_ml')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('right_foot_angle_ml')
                end            
                if strcmp(this_variable_name, 'right_foot_angle_yaw')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('right_foot_angle_yaw')
                end

                % clearance above obstacle
                if strcmp(this_variable_name, 'heel_clearance')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('heel_clearance')
                end
                if strcmp(this_variable_name, 'toes_clearance')
                    this.addBasicVariable('marker_trajectories')
                    this.addStretchVariable('toes_clearance')
                end
                
                % is this a normal variable? did we add it successfully
                if ~any(this_variable_name==':') && ~any(strcmp(this_variable_name, this.stretch_variable_names))
                    % no, we haven't added a stretch variable, so add this as a basic and a stretch variable (unless it's stimulus_response, which is special)
                    if ~strcmp(this_variable_name, 'stimulus_response_x')
                        this.addBasicVariable(this_variable_name)
                        this.addStretchVariable(this_variable_name)
                    end
                end
            end
        end
        
        % interface
        function result = isBasicVariable(this, variable_name)
            result = any(strcmp(this.basic_variable_names, variable_name));
        end
        function T = getRecordingTimeStart(this)
            T = inf;
            for i_variable = 1 : size(this.basic_variable_names)
                this_variable_name = this.basic_variable_names{i_variable};
                this_variable_time = this.time_data.(this_variable_name);
                T = min([this_variable_time(1), T]);
            end
        end
        function T = getRecordingTimeEnd(this)
            T = inf;
            for i_variable = 1 : size(this.basic_variable_names)
                this_variable_name = this.basic_variable_names{i_variable};
                this_variable_time = this.time_data.(this_variable_name);
                T = min([this_variable_time(end), T]);
            end
        end
        function time_data = getTimeData(this, variable_name)
            % check if this is a sub-variable
            if any(variable_name==':')
                this_variable_split = strsplit(variable_name, ':');
                this_variable_type = this_variable_split{1};
                this_variable_label = this_variable_split{2};
                if strcmp(this_variable_type, 'elementary') || strcmp(this_variable_type, 'derivative')
                    % not a true compound, but label is elementary variable
                    name_to_use = this_variable_label;
                else                    
                    % this is a true compound
                    name_to_use = [this_variable_type '_trajectories'];
                end
            else
                name_to_use = variable_name;
            end
            
            time_data = this.time_data.(name_to_use);
        end
        function [variable_data, variable_directions] = getBasicVariableData(this, variable_name)
            % unpack name
            if any(variable_name==':')
                this_variable_split = strsplit(variable_name, ':');
                this_variable_type = this_variable_split{1};
                this_variable_label = this_variable_split{2};
                
                if strcmp(this_variable_type, 'elementary') || strcmp(this_variable_type, 'derivative')
                    % not a true compound, but label is elementary variable
                    variable_data = this.basic_variable_data.(this_variable_label);
                    if nargout > 1
                        variable_directions = this.basic_variable_directions.(this_variable_label);
                    end
                else                    
                    % this is a true compound
                    name_to_use = [this_variable_type '_trajectories'];
                    trajectory_data = this.basic_variable_data.(name_to_use);
                    trajectory_labels = this.basic_variable_labels.(name_to_use);
                    variable_data = trajectory_data(:, strcmp(trajectory_labels, this_variable_label));
                    
                    if nargout > 1
                        direction_data = this.basic_variable_directions.(name_to_use);
                        variable_directions = direction_data(:, strcmp(trajectory_labels, this_variable_label));
                    end
                end
            else
                variable_data = this.basic_variable_data.(variable_name);
                if nargout > 1
                    variable_directions = this.basic_variable_directions.(variable_name);
                end
            end
            
        end
        
        function prepareBasicVariables(this, trial_type, trial_number, variables_to_prepare)
            if nargin < 4
                variables_to_prepare = this.basic_variable_names;
            end
            data_to_remove_header = this.subject_settings.get('data_to_remove_header', 1);
            data_to_remove = this.subject_settings.get('data_to_remove', 1);
            if ~isempty(data_to_remove)
                this_trial_rows = strcmp(data_to_remove(:, strcmp(data_to_remove_header, 'trial_number')), num2str(trial_number));
                data_to_remove_this_trial = data_to_remove(this_trial_rows, 1);
            else
                data_to_remove_this_trial = [];
            end
            
            % clear out old data
            this.basic_variable_data = struct;
            this.basic_variable_labels = struct;
            this.basic_variable_directions = struct;
            this.stretch_variable_data = struct;
            this.time_data = struct;
            this.trial_type = trial_type;
            this.trial_number = trial_number;

            % load basic variables
            for i_variable = 1 : length(variables_to_prepare)
                variable_name = variables_to_prepare{i_variable};
                
                % try loading
                [data, time, sampling_rate, labels, directions, success] = loadData(this.date, this.subject_id, trial_type, trial_number, variable_name, 'optional'); %#ok<ASGLU>
                
                % store
                if success
                    this.basic_variable_data.(variable_name) = data;
                    this.time_data.(variable_name) = time;
                    this.basic_variable_labels.(variable_name) = labels;
                    this.basic_variable_directions.(variable_name) = directions;
                end
                
                % check if this is a compound name, listing a loaded variable and a label
                if any(variable_name==':')
                    this_variable_split = strsplit(variable_name, ':');
                    this_variable_type = this_variable_split{1};
                    this_variable_label = this_variable_split{2};
                    if strcmp(this_variable_type, 'elementary') || strcmp(this_variable_type, 'derivative')
                        % not a true compound, but label is elementary variable
                        [data, time, sampling_rate, labels, directions, success] = loadData(this.date, this.subject_id, trial_type, trial_number, this_variable_label, 'optional'); %#ok<ASGLU>
                        if success
                            this.basic_variable_data.(this_variable_label) = data;
                            this.time_data.(this_variable_label) = time;
                            this.basic_variable_labels.(this_variable_label) = labels;
                            this.basic_variable_directions.(this_variable_label) = directions;
                        end
                    else                    
                        % this is a true compound
                        if any(strcmp(this.basic_variable_labels.([this_variable_type '_trajectories']), this_variable_label))
                            % found this label in the data, so all is good
                            success = 1;
                        end
                    end
                end

                % remove data flagged in subject settings
                if any(strcmp(data_to_remove_this_trial, variable_name))
                    if any(variable_name==':')
                        this_variable_split = strsplit(variable_name, ':');
                        this_variable_type = this_variable_split{1};
                        this_variable_label = this_variable_split{2};
                        column = strcmp(this.basic_variable_labels.([this_variable_type '_trajectories']), this_variable_label);
                        this.basic_variable_data.([this_variable_type '_trajectories'])(:, column) = NaN;
                    else
                        this.basic_variable_data.(variable_name)(:) = NaN;
                    end
                end
                
                % report and store failure
                if ~success && ~any(strcmp(variable_name, this.basic_variable_load_failures))
                    disp(['Warning: tried to process variable ' variable_name ' but failed to load data.'])
                    this.basic_variable_load_failures = [this.basic_variable_load_failures; variable_name];
                end
                
                
            end
        end
        function stretch_variables = calculateStretchVariables(this, stretch_times, stance_foot_data, relevant_condition_data, variables_to_calculate)
            if nargin < 6
                variables_to_calculate = this.stretch_variable_names;
            end
            
            number_of_stretch_variables = length(variables_to_calculate);
            number_of_stretches = size(stretch_times, 1);
            number_of_bands = size(stance_foot_data, 2);
            stretch_variables = cell(number_of_stretch_variables, 1);
            
            for i_variable = 1 : number_of_stretch_variables
                variable_name = variables_to_calculate{i_variable};
                
                this.registerStretchVariableDirections(variable_name);
                
                % extract and normalize data from stretches
                for i_stretch = 1 : number_of_stretches
                    % create containers
                    stretch_data = [];
                    
                    % extract time
                    this_stretch_times = stretch_times(i_stretch, :);
                    
                    % calculate normalized stretch data for the basic variables
                    if this.isBasicVariable(variable_name) || any(variable_name==':')
                        if ~any(strcmp(this.basic_variable_load_failures, variable_name))
                            stretch_data = this.getTimeNormalizedData(variable_name, this_stretch_times);
                        else
                            stretch_data = NaN;
                        end
                        
                    end
                    
                    % calculate stretch variables that are not basic variables or need special attention
                    % REMINDER: if you add a variable here, make sure to also add it below in registerStretchVariableDirections
                    
                    % use an initial value as a reference
                    if strcmp(variable_name, 'lheel_from_mpsis_initial_x')
                        LHEE_x = this.getTimeNormalizedData('marker:LHEE_x', this_stretch_times);
                        LPSI_x = this.getTimeNormalizedData('marker:LPSI_x', this_stretch_times);
                        RPSI_x = this.getTimeNormalizedData('marker:RPSI_x', this_stretch_times);
                        mpsis_x = (LPSI_x + RPSI_x) * 0.5;
                        stretch_data = LHEE_x - mpsis_x(1);
                    end
                    if strcmp(variable_name, 'rheel_from_mpsis_initial_x')
                        RHEE_x = this.getTimeNormalizedData('marker:RHEE_x', this_stretch_times);
                        LPSI_x = this.getTimeNormalizedData('marker:LPSI_x', this_stretch_times);
                        RPSI_x = this.getTimeNormalizedData('marker:RPSI_x', this_stretch_times);
                        mpsis_x = (LPSI_x + RPSI_x) * 0.5;
                        stretch_data = RHEE_x - mpsis_x(1);
                    end
                    if strcmp(variable_name, 'mpsis_from_mpsis_initial_x')
                        LPSI_x = this.getTimeNormalizedData('marker:LPSI_x', this_stretch_times);
                        RPSI_x = this.getTimeNormalizedData('marker:RPSI_x', this_stretch_times);
                        mpsis_x = (LPSI_x + RPSI_x) * 0.5;
                        stretch_data = mpsis_x - mpsis_x(1);
                    end
                    if strcmp(variable_name, 'lheel_from_com_initial_x')
                        LHEE_x = this.getTimeNormalizedData('marker:LHEE_x', this_stretch_times);
                        com_x =  this.getTimeNormalizedData('com_position:center_of_mass_x', this_stretch_times);
                        stretch_data = LHEE_x - com_x(1);
                    end
                    if strcmp(variable_name, 'rheel_from_com_initial_x')
                        RHEE_x = this.getTimeNormalizedData('marker:RHEE_x', this_stretch_times);
                        com_x =  this.getTimeNormalizedData('com_position:center_of_mass_x', this_stretch_times);
                        stretch_data = RHEE_x - com_x(1);
                    end
                    if strcmp(variable_name, 'com_from_com_initial_x')
                        com_x =  this.getTimeNormalizedData('com_position:center_of_mass_x', this_stretch_times);
                        stretch_data = com_x - com_x(1);
                    end
                    
                    % balance stuff
                    if strcmp(variable_name, 'step_placement_x')
                        LHEE_x = this.getTimeNormalizedData('marker:LHEE_x', this_stretch_times);
                        RHEE_x = this.getTimeNormalizedData('marker:RHEE_x', this_stretch_times);
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            [~, band_end_indices] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stretch_data(i_band) = LHEE_x(band_end_indices) - RHEE_x(band_end_indices);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stretch_data(i_band) = RHEE_x(band_end_indices) - LHEE_x(band_end_indices);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stretch_data(i_band) = NaN;
                            end
                        end
                    end
                    if strcmp(variable_name, 'xcom_x')
                        % get CoM position
                        com_x =  this.getTimeNormalizedData('com_position:center_of_mass_x', this_stretch_times);
                        com_y =  this.getTimeNormalizedData('com_position:center_of_mass_y', this_stretch_times);
                        com_z =  this.getTimeNormalizedData('com_position:center_of_mass_z', this_stretch_times);
                        
                        % get instantaneous leg length
                        LANK_x = this.getTimeNormalizedData('marker:LANK_x', this_stretch_times);
                        RANK_x = this.getTimeNormalizedData('marker:RANK_x', this_stretch_times);
                        LANK_y = this.getTimeNormalizedData('marker:LANK_y', this_stretch_times);
                        RANK_y = this.getTimeNormalizedData('marker:RANK_y', this_stretch_times);
                        LANK_z = this.getTimeNormalizedData('marker:LANK_z', this_stretch_times);
                        RANK_z = this.getTimeNormalizedData('marker:RANK_z', this_stretch_times);
                        
                        stance_ankle_x_data = zeros(size(com_x));
                        stance_ankle_y_data = zeros(size(com_x));
                        stance_ankle_z_data = zeros(size(com_x));
                        for i_band = number_of_bands : -1 : 1
                            % going backward makes a difference for the
                            % junction points between two steps. We want to
                            % use data from the earlier step for BoS, so we
                            % go backward
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stance_ankle_x_data(band_start_index : band_end_index) = RANK_x(band_start_index : band_end_index);
                                stance_ankle_y_data(band_start_index : band_end_index) = RANK_y(band_start_index : band_end_index);
                                stance_ankle_z_data(band_start_index : band_end_index) = RANK_z(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stance_ankle_x_data(band_start_index : band_end_index) = LANK_x(band_start_index : band_end_index);
                                stance_ankle_y_data(band_start_index : band_end_index) = LANK_y(band_start_index : band_end_index);
                                stance_ankle_z_data(band_start_index : band_end_index) = LANK_z(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stance_ankle_x_data(band_start_index : band_end_index) = NaN;
                                stance_ankle_y_data(band_start_index : band_end_index) = NaN;
                                stance_ankle_z_data(band_start_index : band_end_index) = NaN;
                            end
                        end
                        leg_vector_x = com_x - stance_ankle_x_data;
                        leg_vector_y = com_y - stance_ankle_y_data;
                        leg_vector_z = com_z - stance_ankle_z_data;
                        leg_length_data = (leg_vector_x.^2 + leg_vector_y.^2 + leg_vector_z.^2).^(0.5);
                        
                        % calculate XCoM
                        com_x_vel =  this.getTimeNormalizedData('com_velocity:center_of_mass_x', this_stretch_times);
                        omega_0 = (9.81 * leg_length_data.^(-1)).^(0.5);
                       
                        stretch_data = com_x + omega_0.^(-1) .* com_x_vel;
                    end
                    if strcmp(variable_name, 'xcom_y')
                        % get CoM position
                        com_x =  this.getTimeNormalizedData('com_position:center_of_mass_x', this_stretch_times);
                        com_y =  this.getTimeNormalizedData('com_position:center_of_mass_y', this_stretch_times);
                        com_z =  this.getTimeNormalizedData('com_position:center_of_mass_z', this_stretch_times);
                        
                        % get instantaneous leg length
                        LANK_x = this.getTimeNormalizedData('marker:LANK_x', this_stretch_times);
                        RANK_x = this.getTimeNormalizedData('marker:RANK_x', this_stretch_times);
                        LANK_y = this.getTimeNormalizedData('marker:LANK_y', this_stretch_times);
                        RANK_y = this.getTimeNormalizedData('marker:RANK_y', this_stretch_times);
                        LANK_z = this.getTimeNormalizedData('marker:LANK_z', this_stretch_times);
                        RANK_z = this.getTimeNormalizedData('marker:RANK_z', this_stretch_times);
                        stance_ankle_x_data = zeros(size(com_x));
                        stance_ankle_y_data = zeros(size(com_x));
                        stance_ankle_z_data = zeros(size(com_x));
                        for i_band = number_of_bands : -1 : 1
                            % going backward makes a difference for the
                            % junction points between two steps. We want to
                            % use data from the earlier step for BoS, so we
                            % go backward
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stance_ankle_x_data(band_start_index : band_end_index) = RANK_x(band_start_index : band_end_index);
                                stance_ankle_y_data(band_start_index : band_end_index) = RANK_y(band_start_index : band_end_index);
                                stance_ankle_z_data(band_start_index : band_end_index) = RANK_z(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stance_ankle_x_data(band_start_index : band_end_index) = LANK_x(band_start_index : band_end_index);
                                stance_ankle_y_data(band_start_index : band_end_index) = LANK_y(band_start_index : band_end_index);
                                stance_ankle_z_data(band_start_index : band_end_index) = LANK_z(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stance_ankle_x_data(band_start_index : band_end_index) = NaN;
                                stance_ankle_y_data(band_start_index : band_end_index) = NaN;
                                stance_ankle_z_data(band_start_index : band_end_index) = NaN;
                            end
                        end
                        leg_vector_x = com_x - stance_ankle_x_data;
                        leg_vector_y = com_y - stance_ankle_y_data;
                        leg_vector_z = com_z - stance_ankle_z_data;
                        leg_length_data = (leg_vector_x.^2 + leg_vector_y.^2 + leg_vector_z.^2).^(0.5);
                        
                        % calculate XCoM
                        com_y_vel =  this.getTimeNormalizedData('com_velocity:center_of_mass_y', this_stretch_times);
                        omega_0 = (9.81 * leg_length_data.^(-1)).^(0.5);
                       
                        stretch_data = com_y + omega_0.^(-1) .* com_y_vel;                    
                    end
                    if strcmp(variable_name, 'mos_x')
                        % first calculate base of support
                        LTOEL_x = this.getTimeNormalizedData('marker:LTOEL_x', this_stretch_times);
                        RTOEL_x = this.getTimeNormalizedData('marker:RTOEL_x', this_stretch_times);
                        bos_x_data = zeros(size(LTOEL_x));
                        for i_band = 1 : number_of_bands
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                bos_x_data(band_start_index : band_end_index) = RTOEL_x(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                bos_x_data(band_start_index : band_end_index) = LTOEL_x(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                bos_x_data(band_start_index : band_end_index) = NaN;
                            end
                        end
                        
                        % set BoS to NaN at junction points between two steps
                        for i_band = 2 : number_of_bands
                            band_start_index = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            bos_x_data(band_start_index) = NaN;
                        end
                        
                        % now calculate XCoM - BoS
                        xcom_x = stretch_variables{strcmp(this.stretch_variable_names, 'xcom_x')}(:, i_stretch);
                        stretch_data = xcom_x - bos_x_data;
                    end
                    if strcmp(variable_name, 'mos_y')
                        % first calculate base of support
                        LTOEL_y = this.getTimeNormalizedData('marker:LTOEL_y', this_stretch_times);
                        RTOEL_y = this.getTimeNormalizedData('marker:RTOEL_y', this_stretch_times);
                        bos_y_data = zeros(size(LTOEL_y));
                        for i_band = 1 : number_of_bands
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                bos_y_data(band_start_index : band_end_index) = RTOEL_y(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                bos_y_data(band_start_index : band_end_index) = LTOEL_y(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                bos_y_data(band_start_index : band_end_index) = NaN;
                            end
                        end
                        
                        % set BoS to NaN at junction points between two steps
                        for i_band = 2 : number_of_bands
                            band_start_index = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            bos_y_data(band_start_index) = NaN;
                        end
                        
                        % now calculate XCoM - BoS
                        xcom_y = stretch_variables{strcmp(this.stretch_variable_names, 'xcom_y')}(:, i_stretch);
                        stretch_data = xcom_y - bos_y_data;
                    end
                    if strcmp(variable_name, 'xcom_mpsis_x')
                        % get mpsis position
                        LPSI_x = this.getTimeNormalizedData('marker:LPSI_x', this_stretch_times);
                        RPSI_x = this.getTimeNormalizedData('marker:RPSI_x', this_stretch_times);
                        LPSI_y = this.getTimeNormalizedData('marker:LPSI_y', this_stretch_times);
                        RPSI_y = this.getTimeNormalizedData('marker:RPSI_y', this_stretch_times);
                        LPSI_z = this.getTimeNormalizedData('marker:LPSI_z', this_stretch_times);
                        RPSI_z = this.getTimeNormalizedData('marker:RPSI_z', this_stretch_times);
                        
                        mpsis_x = (LPSI_x + RPSI_x) * 0.5;
                        mpsis_y = (LPSI_y + RPSI_y) * 0.5;
                        mpsis_z = (LPSI_z + RPSI_z) * 0.5;
                        
                        % get instantaneous leg length
                        LANK_x = this.getTimeNormalizedData('marker:LANK_x', this_stretch_times);
                        RANK_x = this.getTimeNormalizedData('marker:RANK_x', this_stretch_times);
                        LANK_y = this.getTimeNormalizedData('marker:LANK_y', this_stretch_times);
                        RANK_y = this.getTimeNormalizedData('marker:RANK_y', this_stretch_times);
                        LANK_z = this.getTimeNormalizedData('marker:LANK_z', this_stretch_times);
                        RANK_z = this.getTimeNormalizedData('marker:RANK_z', this_stretch_times);
                        stance_ankle_x_data = zeros(size(mpsis_x));
                        stance_ankle_y_data = zeros(size(mpsis_x));
                        stance_ankle_z_data = zeros(size(mpsis_x));
                        for i_band = number_of_bands : -1 : 1
                            % going backward makes a difference for the
                            % junction points between two steps. We want to
                            % use data from the earlier step for BoS, so we
                            % go backward
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stance_ankle_x_data(band_start_index : band_end_index) = RANK_x(band_start_index : band_end_index);
                                stance_ankle_y_data(band_start_index : band_end_index) = RANK_y(band_start_index : band_end_index);
                                stance_ankle_z_data(band_start_index : band_end_index) = RANK_z(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stance_ankle_x_data(band_start_index : band_end_index) = LANK_x(band_start_index : band_end_index);
                                stance_ankle_y_data(band_start_index : band_end_index) = LANK_y(band_start_index : band_end_index);
                                stance_ankle_z_data(band_start_index : band_end_index) = LANK_z(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stance_ankle_x_data(band_start_index : band_end_index) = NaN;
                                stance_ankle_y_data(band_start_index : band_end_index) = NaN;
                                stance_ankle_z_data(band_start_index : band_end_index) = NaN;
                            end
                        end
                        leg_vector_x = mpsis_x - stance_ankle_x_data;
                        leg_vector_y = mpsis_y - stance_ankle_y_data;
                        leg_vector_z = mpsis_z - stance_ankle_z_data;
                        leg_length_data = (leg_vector_x.^2 + leg_vector_y.^2 + leg_vector_z.^2).^(0.5);
                        
                        % calculate XCoM_mpsis_x
                        LPSI_x_vel = this.getTimeNormalizedData('derivative:LPSI_x_vel', this_stretch_times);
                        RPSI_x_vel = this.getTimeNormalizedData('derivative:RPSI_x_vel', this_stretch_times);
                        mpsis_x_vel = (LPSI_x_vel + RPSI_x_vel) * 0.5;
                        omega_0 = (9.81 * leg_length_data.^(-1)).^(0.5);
                       
                        stretch_data = mpsis_x + omega_0.^(-1) .* mpsis_x_vel;
                    end
                    if strcmp(variable_name, 'xcom_mpsis_y')
                        % get mpsis position
                        LPSI_x = this.getTimeNormalizedData('marker:LPSI_x', this_stretch_times);
                        RPSI_x = this.getTimeNormalizedData('marker:RPSI_x', this_stretch_times);
                        LPSI_y = this.getTimeNormalizedData('marker:LPSI_y', this_stretch_times);
                        RPSI_y = this.getTimeNormalizedData('marker:RPSI_y', this_stretch_times);
                        LPSI_z = this.getTimeNormalizedData('marker:LPSI_z', this_stretch_times);
                        RPSI_z = this.getTimeNormalizedData('marker:RPSI_z', this_stretch_times);
                        
                        mpsis_x = (LPSI_x + RPSI_x) * 0.5;
                        mpsis_y = (LPSI_y + RPSI_y) * 0.5;
                        mpsis_z = (LPSI_z + RPSI_z) * 0.5;
                        
                        % get instantaneous leg length
                        LANK_x = this.getTimeNormalizedData('marker:LANK_x', this_stretch_times);
                        RANK_x = this.getTimeNormalizedData('marker:RANK_x', this_stretch_times);
                        LANK_y = this.getTimeNormalizedData('marker:LANK_y', this_stretch_times);
                        RANK_y = this.getTimeNormalizedData('marker:RANK_y', this_stretch_times);
                        LANK_z = this.getTimeNormalizedData('marker:LANK_z', this_stretch_times);
                        RANK_z = this.getTimeNormalizedData('marker:RANK_z', this_stretch_times);
                        stance_ankle_x_data = zeros(size(mpsis_x));
                        stance_ankle_y_data = zeros(size(mpsis_x));
                        stance_ankle_z_data = zeros(size(mpsis_x));
                        for i_band = number_of_bands : -1 : 1
                            % going backward makes a difference for the
                            % junction points between two steps. We want to
                            % use data from the earlier step for BoS, so we
                            % go backward
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stance_ankle_x_data(band_start_index : band_end_index) = RANK_x(band_start_index : band_end_index);
                                stance_ankle_y_data(band_start_index : band_end_index) = RANK_y(band_start_index : band_end_index);
                                stance_ankle_z_data(band_start_index : band_end_index) = RANK_z(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stance_ankle_x_data(band_start_index : band_end_index) = LANK_x(band_start_index : band_end_index);
                                stance_ankle_y_data(band_start_index : band_end_index) = LANK_y(band_start_index : band_end_index);
                                stance_ankle_z_data(band_start_index : band_end_index) = LANK_z(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stance_ankle_x_data(band_start_index : band_end_index) = NaN;
                                stance_ankle_y_data(band_start_index : band_end_index) = NaN;
                                stance_ankle_z_data(band_start_index : band_end_index) = NaN;
                            end
                        end
                        leg_vector_x = mpsis_x - stance_ankle_x_data;
                        leg_vector_y = mpsis_y - stance_ankle_y_data;
                        leg_vector_z = mpsis_z - stance_ankle_z_data;
                        leg_length_data = (leg_vector_x.^2 + leg_vector_y.^2 + leg_vector_z.^2).^(0.5);
                        
                        % calculate XCoM_mpsis_y
                        LPSI_y_vel = this.getTimeNormalizedData('derivative:LPSI_y_vel', this_stretch_times);
                        RPSI_y_vel = this.getTimeNormalizedData('derivative:RPSI_y_vel', this_stretch_times);
                        mpsis_y_vel = (LPSI_y_vel + RPSI_y_vel) * 0.5;
                        omega_0 = (9.81 * leg_length_data.^(-1)).^(0.5);
                       
                        stretch_data = mpsis_y + omega_0.^(-1) .* mpsis_y_vel;                    
                    end
                    if strcmp(variable_name, 'mos_mpsis_x')
                        % first calculate base of support
                        LTOEL_x = this.getTimeNormalizedData('marker:LTOEL_x', this_stretch_times);
                        RTOEL_x = this.getTimeNormalizedData('marker:RTOEL_x', this_stretch_times);
                        bos_x_data = zeros(size(LTOEL_x));
                        for i_band = 1 : number_of_bands
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                bos_x_data(band_start_index : band_end_index) = RTOEL_x(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                bos_x_data(band_start_index : band_end_index) = LTOEL_x(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                bos_x_data(band_start_index : band_end_index) = NaN;
                            end
                        end
                        
                        % set BoS to NaN at junction points between two steps
                        for i_band = 2 : number_of_bands
                            band_start_index = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            bos_x_data(band_start_index) = NaN;
                        end
                        
                        % now calculate XCoM - BoS
                        xcom_mpsis_x = stretch_variables{strcmp(this.stretch_variable_names, 'xcom_mpsis_x')}(:, i_stretch);
                        stretch_data = xcom_mpsis_x - bos_x_data;
                    end
                    if strcmp(variable_name, 'mos_mpsis_y')
                        % first calculate base of support
                        LTOEL_y = this.getTimeNormalizedData('marker:LTOEL_y', this_stretch_times);
                        RTOEL_y = this.getTimeNormalizedData('marker:RTOEL_y', this_stretch_times);
                        bos_y_data = zeros(size(LTOEL_y));
                        for i_band = 1 : number_of_bands
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                bos_y_data(band_start_index : band_end_index) = RTOEL_y(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                bos_y_data(band_start_index : band_end_index) = LTOEL_y(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                bos_y_data(band_start_index : band_end_index) = NaN;
                            end
                        end
                        
                        % set BoS to NaN at junction points between two steps
                        for i_band = 2 : number_of_bands
                            band_start_index = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            bos_y_data(band_start_index) = NaN;
                        end
                        
                        % now calculate XCoM - BoS
                        xcom_mpsis_y = stretch_variables{strcmp(this.stretch_variable_names, 'xcom_mpsis_y')}(:, i_stretch);
                        stretch_data = xcom_mpsis_y - bos_y_data;
                    end
                    if strcmp(variable_name, 'com_rough_x')
                        % grab required marker trajectories
                        LASI_trajectory = this.getTimeNormalizedData('marker:LASI_x', this_stretch_times);
                        RASI_trajectory = this.getTimeNormalizedData('marker:RASI_x', this_stretch_times);
                        LPSI_trajectory = this.getTimeNormalizedData('marker:LPSI_x', this_stretch_times);
                        RPSI_trajectory = this.getTimeNormalizedData('marker:RPSI_x', this_stretch_times);
                        LKNE_trajectory = this.getTimeNormalizedData('marker:LKNE_x', this_stretch_times);
                        RKNE_trajectory = this.getTimeNormalizedData('marker:RKNE_x', this_stretch_times);
                        LANK_trajectory = this.getTimeNormalizedData('marker:LANK_x', this_stretch_times);
                        RANK_trajectory = this.getTimeNormalizedData('marker:RANK_x', this_stretch_times);
                        LTOE_trajectory = this.getTimeNormalizedData('marker:LTOE_x', this_stretch_times);
                        RTOE_trajectory = this.getTimeNormalizedData('marker:RTOE_x', this_stretch_times);
                        C7_trajectory = this.getTimeNormalizedData('marker:C7_x', this_stretch_times);
                        LFHD_trajectory = this.getTimeNormalizedData('marker:LFHD_x', this_stretch_times);
                        RFHD_trajectory = this.getTimeNormalizedData('marker:RFHD_x', this_stretch_times);
                        LBHD_trajectory = this.getTimeNormalizedData('marker:LBHD_x', this_stretch_times);
                        RBHD_trajectory = this.getTimeNormalizedData('marker:RBHD_x', this_stretch_times);

                        % calculate joint center trajectories
                        hip_left_center = (LASI_trajectory + LPSI_trajectory) * 0.5;
                        hip_right_center = (RASI_trajectory + RPSI_trajectory) * 0.5;
                        head_tip = (LFHD_trajectory + RFHD_trajectory + LBHD_trajectory + RBHD_trajectory) * 0.25;

                        % calculate segment center trajectories
                        pelvis_center = (hip_left_center + hip_right_center) * 0.5;
                        trunk_center = (pelvis_center + C7_trajectory) * 0.5;
                        head_center = (C7_trajectory + head_tip) * 0.5;
                        thigh_left_center = (hip_left_center + LKNE_trajectory) * 0.5;
                        thigh_right_center = (hip_right_center + RKNE_trajectory) * 0.5;
                        shank_left_center = (LKNE_trajectory + LANK_trajectory) * 0.5;
                        shank_right_center = (RKNE_trajectory + RANK_trajectory) * 0.5;
                        foot_left_center = (LANK_trajectory + LTOE_trajectory) * 0.5;
                        foot_right_center = (RANK_trajectory + RTOE_trajectory) * 0.5;

                        % define weight factors                
                        % according to R. Dumas , L. Cheze, J.-P. Verriest: "Adjustments to McConville et al. and Young et al. body
                        % segment inertial parameters", Journal of Biomechanics 40 (2007) 543?553
                        head_mass_factor      = 0.067;
                        torso_mass_factor     = 0.333;
                        arm_mass_factor       = 0.024;
                        forearm_mass_factor   = 0.017;
                        hand_mass_factor      = 0.006;
                        pelvis_mass_factor    = 0.142;
                        thigh_mass_factor     = 0.123;
                        shank_mass_factor     = 0.048;
                        foot_mass_factor      = 0.012;

                        w_head = head_mass_factor;
                        w_hat = torso_mass_factor + 2*arm_mass_factor + 2*forearm_mass_factor + 2*hand_mass_factor;
                        w_pelvis = pelvis_mass_factor;
                        w_thigh = thigh_mass_factor;
                        w_shank = shank_mass_factor;
                        w_foot = foot_mass_factor;

                        % calculate CoM as weighted sum
                        com_trajectory = ...
                            w_head * head_center ...
                            + w_hat * trunk_center ...
                            + w_pelvis * pelvis_center ...
                            + w_thigh * thigh_left_center ...
                            + w_thigh * thigh_right_center ...
                            + w_shank * shank_left_center ...
                            + w_shank * shank_right_center ...
                            + w_foot * foot_left_center ...
                            + w_foot * foot_right_center ...
                            ;
                        stretch_data = com_trajectory;
                    end
                    if strcmp(variable_name, 'com_rough_y')
                        % grab required marker trajectories
                        LASI_trajectory = this.getTimeNormalizedData('marker:LASI_y', this_stretch_times);
                        RASI_trajectory = this.getTimeNormalizedData('marker:RASI_y', this_stretch_times);
                        LPSI_trajectory = this.getTimeNormalizedData('marker:LPSI_y', this_stretch_times);
                        RPSI_trajectory = this.getTimeNormalizedData('marker:RPSI_y', this_stretch_times);
                        LKNE_trajectory = this.getTimeNormalizedData('marker:LKNE_y', this_stretch_times);
                        RKNE_trajectory = this.getTimeNormalizedData('marker:RKNE_y', this_stretch_times);
                        LANK_trajectory = this.getTimeNormalizedData('marker:LANK_y', this_stretch_times);
                        RANK_trajectory = this.getTimeNormalizedData('marker:RANK_y', this_stretch_times);
                        LTOE_trajectory = this.getTimeNormalizedData('marker:LTOE_y', this_stretch_times);
                        RTOE_trajectory = this.getTimeNormalizedData('marker:RTOE_y', this_stretch_times);
                        C7_trajectory = this.getTimeNormalizedData('marker:C7_y', this_stretch_times);
                        LFHD_trajectory = this.getTimeNormalizedData('marker:LFHD_y', this_stretch_times);
                        RFHD_trajectory = this.getTimeNormalizedData('marker:RFHD_y', this_stretch_times);
                        LBHD_trajectory = this.getTimeNormalizedData('marker:LBHD_y', this_stretch_times);
                        RBHD_trajectory = this.getTimeNormalizedData('marker:RBHD_y', this_stretch_times);

                        % calculate joint center trajectories
                        hip_left_center = (LASI_trajectory + LPSI_trajectory) * 0.5;
                        hip_right_center = (RASI_trajectory + RPSI_trajectory) * 0.5;
                        head_tip = (LFHD_trajectory + RFHD_trajectory + LBHD_trajectory + RBHD_trajectory) * 0.25;

                        % calculate segment center trajectories
                        pelvis_center = (hip_left_center + hip_right_center) * 0.5;
                        trunk_center = (pelvis_center + C7_trajectory) * 0.5;
                        head_center = (C7_trajectory + head_tip) * 0.5;
                        thigh_left_center = (hip_left_center + LKNE_trajectory) * 0.5;
                        thigh_right_center = (hip_right_center + RKNE_trajectory) * 0.5;
                        shank_left_center = (LKNE_trajectory + LANK_trajectory) * 0.5;
                        shank_right_center = (RKNE_trajectory + RANK_trajectory) * 0.5;
                        foot_left_center = (LANK_trajectory + LTOE_trajectory) * 0.5;
                        foot_right_center = (RANK_trajectory + RTOE_trajectory) * 0.5;

                        % define weight factors                
                        % according to R. Dumas , L. Cheze, J.-P. Verriest: "Adjustments to McConville et al. and Young et al. body
                        % segment inertial parameters", Journal of Biomechanics 40 (2007) 543?553
                        head_mass_factor      = 0.067;
                        torso_mass_factor     = 0.333;
                        arm_mass_factor       = 0.024;
                        forearm_mass_factor   = 0.017;
                        hand_mass_factor      = 0.006;
                        pelvis_mass_factor    = 0.142;
                        thigh_mass_factor     = 0.123;
                        shank_mass_factor     = 0.048;
                        foot_mass_factor      = 0.012;

                        w_head = head_mass_factor;
                        w_hat = torso_mass_factor + 2*arm_mass_factor + 2*forearm_mass_factor + 2*hand_mass_factor;
                        w_pelvis = pelvis_mass_factor;
                        w_thigh = thigh_mass_factor;
                        w_shank = shank_mass_factor;
                        w_foot = foot_mass_factor;

                        % calculate CoM as weighted sum
                        com_trajectory = ...
                            w_head * head_center ...
                            + w_hat * trunk_center ...
                            + w_pelvis * pelvis_center ...
                            + w_thigh * thigh_left_center ...
                            + w_thigh * thigh_right_center ...
                            + w_shank * shank_left_center ...
                            + w_shank * shank_right_center ...
                            + w_foot * foot_left_center ...
                            + w_foot * foot_right_center ...
                            ;
                        stretch_data = com_trajectory;
                    end
                    if strcmp(variable_name, 'com_rough_z')
                        % grab required marker trajectories
                        LASI_trajectory = this.getTimeNormalizedData('marker:LASI_z', this_stretch_times);
                        RASI_trajectory = this.getTimeNormalizedData('marker:RASI_z', this_stretch_times);
                        LPSI_trajectory = this.getTimeNormalizedData('marker:LPSI_z', this_stretch_times);
                        RPSI_trajectory = this.getTimeNormalizedData('marker:RPSI_z', this_stretch_times);
                        LKNE_trajectory = this.getTimeNormalizedData('marker:LKNE_z', this_stretch_times);
                        RKNE_trajectory = this.getTimeNormalizedData('marker:RKNE_z', this_stretch_times);
                        LANK_trajectory = this.getTimeNormalizedData('marker:LANK_z', this_stretch_times);
                        RANK_trajectory = this.getTimeNormalizedData('marker:RANK_z', this_stretch_times);
                        LTOE_trajectory = this.getTimeNormalizedData('marker:LTOE_z', this_stretch_times);
                        RTOE_trajectory = this.getTimeNormalizedData('marker:RTOE_z', this_stretch_times);
                        C7_trajectory = this.getTimeNormalizedData('marker:C7_z', this_stretch_times);
                        LFHD_trajectory = this.getTimeNormalizedData('marker:LFHD_z', this_stretch_times);
                        RFHD_trajectory = this.getTimeNormalizedData('marker:RFHD_z', this_stretch_times);
                        LBHD_trajectory = this.getTimeNormalizedData('marker:LBHD_z', this_stretch_times);
                        RBHD_trajectory = this.getTimeNormalizedData('marker:RBHD_z', this_stretch_times);

                        % calculate joint center trajectories
                        hip_left_center = (LASI_trajectory + LPSI_trajectory) * 0.5;
                        hip_right_center = (RASI_trajectory + RPSI_trajectory) * 0.5;
                        head_tip = (LFHD_trajectory + RFHD_trajectory + LBHD_trajectory + RBHD_trajectory) * 0.25;

                        % calculate segment center trajectories
                        pelvis_center = (hip_left_center + hip_right_center) * 0.5;
                        trunk_center = (pelvis_center + C7_trajectory) * 0.5;
                        head_center = (C7_trajectory + head_tip) * 0.5;
                        thigh_left_center = (hip_left_center + LKNE_trajectory) * 0.5;
                        thigh_right_center = (hip_right_center + RKNE_trajectory) * 0.5;
                        shank_left_center = (LKNE_trajectory + LANK_trajectory) * 0.5;
                        shank_right_center = (RKNE_trajectory + RANK_trajectory) * 0.5;
                        foot_left_center = (LANK_trajectory + LTOE_trajectory) * 0.5;
                        foot_right_center = (RANK_trajectory + RTOE_trajectory) * 0.5;

                        % define weight factors                
                        % according to R. Dumas , L. Cheze, J.-P. Verriest: "Adjustments to McConville et al. and Young et al. body
                        % segment inertial parameters", Journal of Biomechanics 40 (2007) 543?553
                        head_mass_factor      = 0.067;
                        torso_mass_factor     = 0.333;
                        arm_mass_factor       = 0.024;
                        forearm_mass_factor   = 0.017;
                        hand_mass_factor      = 0.006;
                        pelvis_mass_factor    = 0.142;
                        thigh_mass_factor     = 0.123;
                        shank_mass_factor     = 0.048;
                        foot_mass_factor      = 0.012;

                        w_head = head_mass_factor;
                        w_hat = torso_mass_factor + 2*arm_mass_factor + 2*forearm_mass_factor + 2*hand_mass_factor;
                        w_pelvis = pelvis_mass_factor;
                        w_thigh = thigh_mass_factor;
                        w_shank = shank_mass_factor;
                        w_foot = foot_mass_factor;

                        % calculate CoM as weighted sum
                        com_trajectory = ...
                            w_head * head_center ...
                            + w_hat * trunk_center ...
                            + w_pelvis * pelvis_center ...
                            + w_thigh * thigh_left_center ...
                            + w_thigh * thigh_right_center ...
                            + w_shank * shank_left_center ...
                            + w_shank * shank_right_center ...
                            + w_foot * foot_left_center ...
                            + w_foot * foot_right_center ...
                            ;
                        stretch_data = com_trajectory;
                    end
                    if strcmp(variable_name, 'com_rough_x_vel')
                        % grab required marker trajectories
                        LASI_trajectory = this.getTimeNormalizedData('derivative:LASI_x_vel', this_stretch_times);
                        RASI_trajectory = this.getTimeNormalizedData('derivative:RASI_x_vel', this_stretch_times);
                        LPSI_trajectory = this.getTimeNormalizedData('derivative:LPSI_x_vel', this_stretch_times);
                        RPSI_trajectory = this.getTimeNormalizedData('derivative:RPSI_x_vel', this_stretch_times);
                        LKNE_trajectory = this.getTimeNormalizedData('derivative:LKNE_x_vel', this_stretch_times);
                        RKNE_trajectory = this.getTimeNormalizedData('derivative:RKNE_x_vel', this_stretch_times);
                        LANK_trajectory = this.getTimeNormalizedData('derivative:LANK_x_vel', this_stretch_times);
                        RANK_trajectory = this.getTimeNormalizedData('derivative:RANK_x_vel', this_stretch_times);
                        LTOE_trajectory = this.getTimeNormalizedData('derivative:LTOE_x_vel', this_stretch_times);
                        RTOE_trajectory = this.getTimeNormalizedData('derivative:RTOE_x_vel', this_stretch_times);
                        C7_trajectory = this.getTimeNormalizedData('derivative:C7_x_vel', this_stretch_times);
                        LFHD_trajectory = this.getTimeNormalizedData('derivative:LFHD_x_vel', this_stretch_times);
                        RFHD_trajectory = this.getTimeNormalizedData('derivative:RFHD_x_vel', this_stretch_times);
                        LBHD_trajectory = this.getTimeNormalizedData('derivative:LBHD_x_vel', this_stretch_times);
                        RBHD_trajectory = this.getTimeNormalizedData('derivative:RBHD_x_vel', this_stretch_times);

                        % calculate joint center trajectories
                        hip_left_center = (LASI_trajectory + LPSI_trajectory) * 0.5;
                        hip_right_center = (RASI_trajectory + RPSI_trajectory) * 0.5;
                        head_tip = (LFHD_trajectory + RFHD_trajectory + LBHD_trajectory + RBHD_trajectory) * 0.25;

                        % calculate segment center trajectories
                        pelvis_center = (hip_left_center + hip_right_center) * 0.5;
                        trunk_center = (pelvis_center + C7_trajectory) * 0.5;
                        head_center = (C7_trajectory + head_tip) * 0.5;
                        thigh_left_center = (hip_left_center + LKNE_trajectory) * 0.5;
                        thigh_right_center = (hip_right_center + RKNE_trajectory) * 0.5;
                        shank_left_center = (LKNE_trajectory + LANK_trajectory) * 0.5;
                        shank_right_center = (RKNE_trajectory + RANK_trajectory) * 0.5;
                        foot_left_center = (LANK_trajectory + LTOE_trajectory) * 0.5;
                        foot_right_center = (RANK_trajectory + RTOE_trajectory) * 0.5;

                        % define weight factors                
                        % according to R. Dumas , L. Cheze, J.-P. Verriest: "Adjustments to McConville et al. and Young et al. body
                        % segment inertial parameters", Journal of Biomechanics 40 (2007) 543?553
                        head_mass_factor      = 0.067;
                        torso_mass_factor     = 0.333;
                        arm_mass_factor       = 0.024;
                        forearm_mass_factor   = 0.017;
                        hand_mass_factor      = 0.006;
                        pelvis_mass_factor    = 0.142;
                        thigh_mass_factor     = 0.123;
                        shank_mass_factor     = 0.048;
                        foot_mass_factor      = 0.012;

                        w_head = head_mass_factor;
                        w_hat = torso_mass_factor + 2*arm_mass_factor + 2*forearm_mass_factor + 2*hand_mass_factor;
                        w_pelvis = pelvis_mass_factor;
                        w_thigh = thigh_mass_factor;
                        w_shank = shank_mass_factor;
                        w_foot = foot_mass_factor;

                        % calculate CoM as weighted sum
                        com_trajectory = ...
                            w_head * head_center ...
                            + w_hat * trunk_center ...
                            + w_pelvis * pelvis_center ...
                            + w_thigh * thigh_left_center ...
                            + w_thigh * thigh_right_center ...
                            + w_shank * shank_left_center ...
                            + w_shank * shank_right_center ...
                            + w_foot * foot_left_center ...
                            + w_foot * foot_right_center ...
                            ;
                        stretch_data = com_trajectory;
                    end
                    if strcmp(variable_name, 'com_rough_y_vel')
                        % grab required marker trajectories
                        LASI_trajectory = this.getTimeNormalizedData('derivative:LASI_y_vel', this_stretch_times);
                        RASI_trajectory = this.getTimeNormalizedData('derivative:RASI_y_vel', this_stretch_times);
                        LPSI_trajectory = this.getTimeNormalizedData('derivative:LPSI_y_vel', this_stretch_times);
                        RPSI_trajectory = this.getTimeNormalizedData('derivative:RPSI_y_vel', this_stretch_times);
                        LKNE_trajectory = this.getTimeNormalizedData('derivative:LKNE_y_vel', this_stretch_times);
                        RKNE_trajectory = this.getTimeNormalizedData('derivative:RKNE_y_vel', this_stretch_times);
                        LANK_trajectory = this.getTimeNormalizedData('derivative:LANK_y_vel', this_stretch_times);
                        RANK_trajectory = this.getTimeNormalizedData('derivative:RANK_y_vel', this_stretch_times);
                        LTOE_trajectory = this.getTimeNormalizedData('derivative:LTOE_y_vel', this_stretch_times);
                        RTOE_trajectory = this.getTimeNormalizedData('derivative:RTOE_y_vel', this_stretch_times);
                        C7_trajectory = this.getTimeNormalizedData('derivative:C7_y_vel', this_stretch_times);
                        LFHD_trajectory = this.getTimeNormalizedData('derivative:LFHD_y_vel', this_stretch_times);
                        RFHD_trajectory = this.getTimeNormalizedData('derivative:RFHD_y_vel', this_stretch_times);
                        LBHD_trajectory = this.getTimeNormalizedData('derivative:LBHD_y_vel', this_stretch_times);
                        RBHD_trajectory = this.getTimeNormalizedData('derivative:RBHD_y_vel', this_stretch_times);

                        % calculate joint center trajectories
                        hip_left_center = (LASI_trajectory + LPSI_trajectory) * 0.5;
                        hip_right_center = (RASI_trajectory + RPSI_trajectory) * 0.5;
                        head_tip = (LFHD_trajectory + RFHD_trajectory + LBHD_trajectory + RBHD_trajectory) * 0.25;

                        % calculate segment center trajectories
                        pelvis_center = (hip_left_center + hip_right_center) * 0.5;
                        trunk_center = (pelvis_center + C7_trajectory) * 0.5;
                        head_center = (C7_trajectory + head_tip) * 0.5;
                        thigh_left_center = (hip_left_center + LKNE_trajectory) * 0.5;
                        thigh_right_center = (hip_right_center + RKNE_trajectory) * 0.5;
                        shank_left_center = (LKNE_trajectory + LANK_trajectory) * 0.5;
                        shank_right_center = (RKNE_trajectory + RANK_trajectory) * 0.5;
                        foot_left_center = (LANK_trajectory + LTOE_trajectory) * 0.5;
                        foot_right_center = (RANK_trajectory + RTOE_trajectory) * 0.5;

                        % define weight factors                
                        % according to R. Dumas , L. Cheze, J.-P. Verriest: "Adjustments to McConville et al. and Young et al. body
                        % segment inertial parameters", Journal of Biomechanics 40 (2007) 543?553
                        head_mass_factor      = 0.067;
                        torso_mass_factor     = 0.333;
                        arm_mass_factor       = 0.024;
                        forearm_mass_factor   = 0.017;
                        hand_mass_factor      = 0.006;
                        pelvis_mass_factor    = 0.142;
                        thigh_mass_factor     = 0.123;
                        shank_mass_factor     = 0.048;
                        foot_mass_factor      = 0.012;

                        w_head = head_mass_factor;
                        w_hat = torso_mass_factor + 2*arm_mass_factor + 2*forearm_mass_factor + 2*hand_mass_factor;
                        w_pelvis = pelvis_mass_factor;
                        w_thigh = thigh_mass_factor;
                        w_shank = shank_mass_factor;
                        w_foot = foot_mass_factor;

                        % calculate CoM as weighted sum
                        com_trajectory = ...
                            w_head * head_center ...
                            + w_hat * trunk_center ...
                            + w_pelvis * pelvis_center ...
                            + w_thigh * thigh_left_center ...
                            + w_thigh * thigh_right_center ...
                            + w_shank * shank_left_center ...
                            + w_shank * shank_right_center ...
                            + w_foot * foot_left_center ...
                            + w_foot * foot_right_center ...
                            ;
                        stretch_data = com_trajectory;
                    end
                    if strcmp(variable_name, 'xcom_rough_x')
                        % get com position
                        com_rough_x = stretch_variables{strcmp(this.stretch_variable_names, 'com_rough_x')}(:, i_stretch);
                        com_rough_y = stretch_variables{strcmp(this.stretch_variable_names, 'com_rough_y')}(:, i_stretch);
                        com_rough_z = stretch_variables{strcmp(this.stretch_variable_names, 'com_rough_z')}(:, i_stretch);
                        
                        % get instantaneous leg length
                        LANK_x = this.getTimeNormalizedData('marker:LANK_x', this_stretch_times);
                        RANK_x = this.getTimeNormalizedData('marker:RANK_x', this_stretch_times);
                        LANK_y = this.getTimeNormalizedData('marker:LANK_y', this_stretch_times);
                        RANK_y = this.getTimeNormalizedData('marker:RANK_y', this_stretch_times);
                        LANK_z = this.getTimeNormalizedData('marker:LANK_z', this_stretch_times);
                        RANK_z = this.getTimeNormalizedData('marker:RANK_z', this_stretch_times);
                        stance_ankle_x_data = zeros(size(com_rough_x));
                        stance_ankle_y_data = zeros(size(com_rough_x));
                        stance_ankle_z_data = zeros(size(com_rough_x));
                        for i_band = number_of_bands : -1 : 1
                            % going backward makes a difference for the
                            % junction points between two steps. We want to
                            % use data from the earlier step for BoS, so we
                            % go backward
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stance_ankle_x_data(band_start_index : band_end_index) = RANK_x(band_start_index : band_end_index);
                                stance_ankle_y_data(band_start_index : band_end_index) = RANK_y(band_start_index : band_end_index);
                                stance_ankle_z_data(band_start_index : band_end_index) = RANK_z(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stance_ankle_x_data(band_start_index : band_end_index) = LANK_x(band_start_index : band_end_index);
                                stance_ankle_y_data(band_start_index : band_end_index) = LANK_y(band_start_index : band_end_index);
                                stance_ankle_z_data(band_start_index : band_end_index) = LANK_z(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stance_ankle_x_data(band_start_index : band_end_index) = NaN;
                                stance_ankle_y_data(band_start_index : band_end_index) = NaN;
                                stance_ankle_z_data(band_start_index : band_end_index) = NaN;
                            end
                        end
                        leg_vector_x = com_rough_x - stance_ankle_x_data;
                        leg_vector_y = com_rough_y - stance_ankle_y_data;
                        leg_vector_z = com_rough_z - stance_ankle_z_data;
                        leg_length_data = (leg_vector_x.^2 + leg_vector_y.^2 + leg_vector_z.^2).^(0.5);
                        
                        % calculate XCoM_rough_x
                        com_rough_x_vel = stretch_variables{strcmp(this.stretch_variable_names, 'com_rough_x_vel')}(:, i_stretch);
                        omega_0 = (9.81 * leg_length_data.^(-1)).^(0.5);
                       
                        stretch_data = com_rough_x + omega_0.^(-1) .* com_rough_x_vel;
                    end
                    if strcmp(variable_name, 'xcom_rough_y')
                        % get mpsis position
                        com_rough_x = stretch_variables{strcmp(this.stretch_variable_names, 'com_rough_x')}(:, i_stretch);
                        com_rough_y = stretch_variables{strcmp(this.stretch_variable_names, 'com_rough_y')}(:, i_stretch);
                        com_rough_z = stretch_variables{strcmp(this.stretch_variable_names, 'com_rough_z')}(:, i_stretch);
                        
                        % get instantaneous leg length
                        LANK_x = this.getTimeNormalizedData('marker:LANK_x', this_stretch_times);
                        RANK_x = this.getTimeNormalizedData('marker:RANK_x', this_stretch_times);
                        LANK_y = this.getTimeNormalizedData('marker:LANK_y', this_stretch_times);
                        RANK_y = this.getTimeNormalizedData('marker:RANK_y', this_stretch_times);
                        LANK_z = this.getTimeNormalizedData('marker:LANK_z', this_stretch_times);
                        RANK_z = this.getTimeNormalizedData('marker:RANK_z', this_stretch_times);
                        stance_ankle_x_data = zeros(size(com_rough_y));
                        stance_ankle_y_data = zeros(size(com_rough_y));
                        stance_ankle_z_data = zeros(size(com_rough_y));
                        for i_band = number_of_bands : -1 : 1
                            % going backward makes a difference for the
                            % junction points between two steps. We want to
                            % use data from the earlier step for BoS, so we
                            % go backward
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stance_ankle_x_data(band_start_index : band_end_index) = RANK_x(band_start_index : band_end_index);
                                stance_ankle_y_data(band_start_index : band_end_index) = RANK_y(band_start_index : band_end_index);
                                stance_ankle_z_data(band_start_index : band_end_index) = RANK_z(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stance_ankle_x_data(band_start_index : band_end_index) = LANK_x(band_start_index : band_end_index);
                                stance_ankle_y_data(band_start_index : band_end_index) = LANK_y(band_start_index : band_end_index);
                                stance_ankle_z_data(band_start_index : band_end_index) = LANK_z(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stance_ankle_x_data(band_start_index : band_end_index) = NaN;
                                stance_ankle_y_data(band_start_index : band_end_index) = NaN;
                                stance_ankle_z_data(band_start_index : band_end_index) = NaN;
                            end
                        end
                        leg_vector_x = com_rough_x - stance_ankle_x_data;
                        leg_vector_y = com_rough_y - stance_ankle_y_data;
                        leg_vector_z = com_rough_z - stance_ankle_z_data;
                        leg_length_data = (leg_vector_x.^2 + leg_vector_y.^2 + leg_vector_z.^2).^(0.5);
                        
                        % calculate xcom_rough_y
                        com_rough_y_vel = stretch_variables{strcmp(this.stretch_variable_names, 'com_rough_y_vel')}(:, i_stretch);
                        omega_0 = (9.81 * leg_length_data.^(-1)).^(0.5);
                       
                        stretch_data = com_rough_y + omega_0.^(-1) .* com_rough_y_vel;
                    end
                    if strcmp(variable_name, 'mos_rough_x')
                        % first calculate base of support
                        LTOEL_x = this.getTimeNormalizedData('marker:LTOEL_x', this_stretch_times);
                        RTOEL_x = this.getTimeNormalizedData('marker:RTOEL_x', this_stretch_times);
                        bos_x_data = zeros(size(LTOEL_x));
                        for i_band = 1 : number_of_bands
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                bos_x_data(band_start_index : band_end_index) = RTOEL_x(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                bos_x_data(band_start_index : band_end_index) = LTOEL_x(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                bos_x_data(band_start_index : band_end_index) = NaN;
                            end
                        end
                        
                        % set BoS to NaN at junction points between two steps
                        for i_band = 2 : number_of_bands
                            band_start_index = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            bos_x_data(band_start_index) = NaN;
                        end
                        
                        % now calculate XCoM - BoS
                        xcom_rough_x = stretch_variables{strcmp(this.stretch_variable_names, 'xcom_rough_x')}(:, i_stretch);
                        stretch_data = xcom_rough_x - bos_x_data;
                    end  
                    if strcmp(variable_name, 'mos_rough_y')
                        % first calculate base of support
                        LTOEL_y = this.getTimeNormalizedData('marker:LTOEL_y', this_stretch_times);
                        RTOEL_y = this.getTimeNormalizedData('marker:RTOEL_y', this_stretch_times);
                        bos_y_data = zeros(size(LTOEL_y));
                        for i_band = 1 : number_of_bands
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                bos_y_data(band_start_index : band_end_index) = RTOEL_y(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                bos_y_data(band_start_index : band_end_index) = LTOEL_y(band_start_index : band_end_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                bos_y_data(band_start_index : band_end_index) = NaN;
                            end
                        end
                        
                        % set BoS to NaN at junction points between two steps
                        for i_band = 2 : number_of_bands
                            band_start_index = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            bos_y_data(band_start_index) = NaN;
                        end
                        
                        % now calculate XCoM - BoS
                        xcom_rough_y = stretch_variables{strcmp(this.stretch_variable_names, 'xcom_rough_y')}(:, i_stretch);
                        stretch_data = xcom_rough_y - bos_y_data;
                    end
                    
                    % gait parameters
                    if strcmp(variable_name, 'step_length')
                        LHEE_y = this.getTimeNormalizedData('marker:LHEE_y', this_stretch_times);
                        RHEE_y = this.getTimeNormalizedData('marker:RHEE_y', this_stretch_times);
                        
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            [~, band_end_indices] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stretch_data(i_band) = LHEE_y(band_end_indices) - RHEE_y(band_end_indices);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stretch_data(i_band) = RHEE_y(band_end_indices) - LHEE_y(band_end_indices);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stretch_data(i_band) = NaN;
                            end
                        end
                    end
                    if strcmp(variable_name, 'step_width')
                        LHEE_x = this.getTimeNormalizedData('marker:LHEE_x', this_stretch_times);
                        RHEE_x = this.getTimeNormalizedData('marker:RHEE_x', this_stretch_times);
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            [~, band_end_indices] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stretch_data(i_band) = abs(LHEE_x(band_end_indices) - RHEE_x(band_end_indices));
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stretch_data(i_band) = abs(RHEE_x(band_end_indices) - LHEE_x(band_end_indices));
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stretch_data(i_band) = NaN;
                            end
                        end
                    end
                    if strcmp(variable_name, 'step_time')
                        stretch_data = diff(this_stretch_times)';
                    end
                    if strcmp(variable_name, 'pushoff_time')
                        % load events
                        event_data = load(['analysis' filesep makeFileName(this.date, this.subject_id, this.trial_type, this.trial_number, 'events.mat')]);
                        left_pushoff_times = event_data.event_data{strcmp(event_data.event_labels, 'left_pushoff')};
                        right_pushoff_times = event_data.event_data{strcmp(event_data.event_labels, 'right_pushoff')};
                        
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                % find first left push-off after band start
                                band_start_time = this_stretch_times(i_band);
                                band_end_time = this_stretch_times(i_band+1);
                                this_pushoff_time = min(left_pushoff_times(left_pushoff_times >= band_start_time));
                                if this_pushoff_time >= band_end_time
                                    this_pushoff_time = band_start_time;
                                end
                                stretch_data(i_band) = this_pushoff_time - band_start_time;
                                
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                band_start_time = this_stretch_times(i_band);
                                band_end_time = this_stretch_times(i_band+1);
                                this_pushoff_time = min(right_pushoff_times(right_pushoff_times >= band_start_time));
                                if this_pushoff_time >= band_end_time
                                    this_pushoff_time = band_start_time;
                                end
                                stretch_data(i_band) = this_pushoff_time - band_start_time;
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stretch_data(i_band) = NaN;
                            end
                        end
                    end
                    if strcmp(variable_name, 'midstance_index')
                        LANK_y = this.getTimeNormalizedData('marker:LANK_y', this_stretch_times);
                        RANK_y = this.getTimeNormalizedData('marker:RANK_y', this_stretch_times);
                        RASI_y = this.getTimeNormalizedData('marker:RASI_y', this_stretch_times);
                        LASI_y = this.getTimeNormalizedData('marker:LASI_y', this_stretch_times);
                        RPSI_y = this.getTimeNormalizedData('marker:RPSI_y', this_stretch_times);
                        LPSI_y = this.getTimeNormalizedData('marker:LPSI_y', this_stretch_times);
                        pelvis_y = (RASI_y + LASI_y + RPSI_y + LPSI_y) * 0.25;
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stance_ankle_this_band = RANK_y(band_start_index : band_end_index);
                                pelvis_this_band = pelvis_y(band_start_index : band_end_index);
                                [~, zero_crossing_index] = min(abs(stance_ankle_this_band - pelvis_this_band));
                                stretch_data(i_band) = band_start_index - 1 + zero_crossing_index;
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stance_ankle_this_band = LANK_y(band_start_index : band_end_index);
                                pelvis_this_band = pelvis_y(band_start_index : band_end_index);
                                [~, zero_crossing_index] = min(abs(stance_ankle_this_band - pelvis_this_band));
                                stretch_data(i_band) = band_start_index - 1 + zero_crossing_index;
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stretch_data(i_band) = NaN;
                            end
                        end
                        
                    end
                    if strcmp(variable_name, 'midstance_time')
                        LANK_y = this.getTimeNormalizedData('marker:LANK_y', this_stretch_times);
                        RANK_y = this.getTimeNormalizedData('marker:RANK_y', this_stretch_times);
                        RASI_y = this.getTimeNormalizedData('marker:RASI_y', this_stretch_times);
                        LASI_y = this.getTimeNormalizedData('marker:LASI_y', this_stretch_times);
                        RPSI_y = this.getTimeNormalizedData('marker:RPSI_y', this_stretch_times);
                        LPSI_y = this.getTimeNormalizedData('marker:LPSI_y', this_stretch_times);
                        pelvis_y = (RASI_y + LASI_y + RPSI_y + LPSI_y) * 0.25;
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            band_start_time = this_stretch_times(i_band);
                            band_end_time = this_stretch_times(i_band+1);
                            band_time = linspace(band_start_time, band_end_time, band_end_index-band_start_index+1);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stance_ankle_this_band = RANK_y(band_start_index : band_end_index);
                                pelvis_this_band = pelvis_y(band_start_index : band_end_index);
                                [~, zero_crossing_index] = min(abs(stance_ankle_this_band - pelvis_this_band));
                                stretch_data(i_band) = band_time(zero_crossing_index) - band_start_time;
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stance_ankle_this_band = LANK_y(band_start_index : band_end_index);
                                pelvis_this_band = pelvis_y(band_start_index : band_end_index);
                                [~, zero_crossing_index] = min(abs(stance_ankle_this_band - pelvis_this_band));
                                stretch_data(i_band) = band_time(zero_crossing_index) - band_start_time;
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stretch_data(i_band) = NaN;
                            end
                        end
                    end
                    if strcmp(variable_name, 'midswing_time')
                        LTOE_y = this.getTimeNormalizedData('marker:LTOE_y', this_stretch_times);
                        RTOE_y = this.getTimeNormalizedData('marker:RTOE_y', this_stretch_times);
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            band_start_time = this_stretch_times(i_band);
                            band_end_time = this_stretch_times(i_band+1);
                            band_time = linspace(band_start_time, band_end_time, band_end_index-band_start_index+1);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stance_toes_this_band = RTOE_y(band_start_index : band_end_index);
                                swing_toes_this_band = LTOE_y(band_start_index : band_end_index);
                                [~, zero_crossing_index] = min(abs(stance_toes_this_band - swing_toes_this_band));
                                stretch_data(i_band) = band_time(zero_crossing_index) - band_start_time;
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stance_toes_this_band = LTOE_y(band_start_index : band_end_index);
                                swing_toes_this_band = RTOE_y(band_start_index : band_end_index);
                                [~, zero_crossing_index] = min(abs(stance_toes_this_band - swing_toes_this_band));
                                stretch_data(i_band) = band_time(zero_crossing_index) - band_start_time;
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stretch_data(i_band) = NaN;
                            end
                        end
                    end
                    if strcmp(variable_name, 'cadence')
                        second_to_minute = 1/60;
                        step_time = stretch_variables{strcmp(this.stretch_variable_names, 'step_time')}(:, i_stretch);
                        stretch_data = (step_time .* second_to_minute).^(-1);
                    end
                    if strcmp(variable_name, 'velocity')
                        step_time = stretch_variables{strcmp(this.stretch_variable_names, 'step_time')}(:, i_stretch);
                        step_length = stretch_variables{strcmp(this.stretch_variable_names, 'step_length')}(:, i_stretch);
                        stretch_data = step_length ./ step_time;
                    end
                    if strcmp(variable_name, 'band_duration')
                        stretch_data = diff(this_stretch_times)';
                    end
                    
                    % segment angles and lengths
                    if strcmp(variable_name, 'leg_length_l')
                        LANK_x = this.getTimeNormalizedData('marker:LANK_x', this_stretch_times);
                        LANK_y = this.getTimeNormalizedData('marker:LANK_y', this_stretch_times);
                        LANK_z = this.getTimeNormalizedData('marker:LANK_z', this_stretch_times);

                        LASI_x = this.getTimeNormalizedData('marker:LASI_x', this_stretch_times);
                        LASI_y = this.getTimeNormalizedData('marker:LASI_y', this_stretch_times);
                        LASI_z = this.getTimeNormalizedData('marker:LASI_z', this_stretch_times);
                        LPSI_x = this.getTimeNormalizedData('marker:LPSI_x', this_stretch_times);
                        LPSI_y = this.getTimeNormalizedData('marker:LPSI_y', this_stretch_times);
                        LPSI_z = this.getTimeNormalizedData('marker:LPSI_z', this_stretch_times);
                        
                        pelvis_midpoint_x = (LASI_x + LPSI_x) * 0.5;
                        pelvis_midpoint_y = (LASI_y + LPSI_y) * 0.5;
                        pelvis_midpoint_z = (LASI_z + LPSI_z) * 0.5;
                        
                        leg_vector_x = pelvis_midpoint_x - LANK_x;
                        leg_vector_y = pelvis_midpoint_y - LANK_y;
                        leg_vector_z = pelvis_midpoint_z - LANK_z;
                        
                        stretch_data = (leg_vector_x.^2 + leg_vector_y.^2 + leg_vector_z.^2).^(0.5);
                    end
                    if strcmp(variable_name, 'leg_length_r')
                        RANK_x = this.getTimeNormalizedData('marker:RANK_x', this_stretch_times);
                        RANK_y = this.getTimeNormalizedData('marker:RANK_y', this_stretch_times);
                        RANK_z = this.getTimeNormalizedData('marker:RANK_z', this_stretch_times);

                        RASI_x = this.getTimeNormalizedData('marker:RASI_x', this_stretch_times);
                        RASI_y = this.getTimeNormalizedData('marker:RASI_y', this_stretch_times);
                        RASI_z = this.getTimeNormalizedData('marker:RASI_z', this_stretch_times);
                        RPSI_x = this.getTimeNormalizedData('marker:RPSI_x', this_stretch_times);
                        RPSI_y = this.getTimeNormalizedData('marker:RPSI_y', this_stretch_times);
                        RPSI_z = this.getTimeNormalizedData('marker:RPSI_z', this_stretch_times);
                        
                        pelvis_midpoint_x = (RASI_x + RPSI_x) * 0.5;
                        pelvis_midpoint_y = (RASI_y + RPSI_y) * 0.5;
                        pelvis_midpoint_z = (RASI_z + RPSI_z) * 0.5;
                        
                        leg_vector_x = pelvis_midpoint_x - RANK_x;
                        leg_vector_y = pelvis_midpoint_y - RANK_y;
                        leg_vector_z = pelvis_midpoint_z - RANK_z;
                        
                        stretch_data = (leg_vector_x.^2 + leg_vector_y.^2 + leg_vector_z.^2).^(0.5);
                    end
                    if strcmp(variable_name, 'leg_angle_ap_l')
                        LANK_y = this.getTimeNormalizedData('marker:LANK_y', this_stretch_times);
                        LANK_z = this.getTimeNormalizedData('marker:LANK_z', this_stretch_times);

                        LASI_y = this.getTimeNormalizedData('marker:LASI_y', this_stretch_times);
                        LASI_z = this.getTimeNormalizedData('marker:LASI_z', this_stretch_times);
                        LPSI_y = this.getTimeNormalizedData('marker:LPSI_y', this_stretch_times);
                        LPSI_z = this.getTimeNormalizedData('marker:LPSI_z', this_stretch_times);
                        
                        pelvis_midpoint_y = (LASI_y + LPSI_y) * 0.5;
                        pelvis_midpoint_z = (LASI_z + LPSI_z) * 0.5;
                        
                        leg_vector_y = pelvis_midpoint_y - LANK_y;
                        leg_vector_z = pelvis_midpoint_z - LANK_z;
                        
                        stretch_data = atan2(leg_vector_y, leg_vector_z);
                    end
                    if strcmp(variable_name, 'leg_angle_ap_r')
                        RANK_y = this.getTimeNormalizedData('marker:RANK_y', this_stretch_times);
                        RANK_z = this.getTimeNormalizedData('marker:RANK_z', this_stretch_times);

                        RASI_y = this.getTimeNormalizedData('marker:RASI_y', this_stretch_times);
                        RASI_z = this.getTimeNormalizedData('marker:RASI_z', this_stretch_times);
                        RPSI_y = this.getTimeNormalizedData('marker:RPSI_y', this_stretch_times);
                        RPSI_z = this.getTimeNormalizedData('marker:RPSI_z', this_stretch_times);
                        
                        pelvis_midpoint_y = (RASI_y + RPSI_y) * 0.5;
                        pelvis_midpoint_z = (RASI_z + RPSI_z) * 0.5;
                        
                        leg_vector_y = pelvis_midpoint_y - RANK_y;
                        leg_vector_z = pelvis_midpoint_z - RANK_z;
                        
                        stretch_data = atan2(leg_vector_y, leg_vector_z);
                    end
                    if strcmp(variable_name, 'head_angle_ap')
                        % calculate angle trajectory
                        C7_y_trajectory = this.getTimeNormalizedData('marker:C7_y', this_stretch_times);
                        LBHD_y_trajectory = this.getTimeNormalizedData('marker:LBHD_y', this_stretch_times);
                        RBHD_y_trajectory = this.getTimeNormalizedData('marker:RBHD_y', this_stretch_times);
                        C7_z_trajectory = this.getTimeNormalizedData('marker:C7_z', this_stretch_times);
                        LBHD_z_trajectory = this.getTimeNormalizedData('marker:LBHD_z', this_stretch_times);
                        RBHD_z_trajectory = this.getTimeNormalizedData('marker:RBHD_z', this_stretch_times);
                        MBHD_y_trajectory = (LBHD_y_trajectory + RBHD_y_trajectory) * 0.5;
                        MBHD_z_trajectory = (LBHD_z_trajectory + RBHD_z_trajectory) * 0.5;
                        head_vector_y = MBHD_y_trajectory - C7_y_trajectory;
                        head_vector_z = MBHD_z_trajectory - C7_z_trajectory;
                        stretch_data = atan2(head_vector_y, head_vector_z);
                    end
                    if strcmp(variable_name, 'head_angle_ml')
                        % calculate angle trajectory
                        C7_x_trajectory = this.getTimeNormalizedData('marker:C7_x', this_stretch_times);
                        LBHD_x_trajectory = this.getTimeNormalizedData('marker:LBHD_x', this_stretch_times);
                        RBHD_x_trajectory = this.getTimeNormalizedData('marker:RBHD_x', this_stretch_times);
                        C7_z_trajectory = this.getTimeNormalizedData('marker:C7_z', this_stretch_times);
                        LBHD_z_trajectory = this.getTimeNormalizedData('marker:LBHD_z', this_stretch_times);
                        RBHD_z_trajectory = this.getTimeNormalizedData('marker:RBHD_z', this_stretch_times);
                        MBHD_x_trajectory = (LBHD_x_trajectory + RBHD_x_trajectory) * 0.5;
                        MBHD_z_trajectory = (LBHD_z_trajectory + RBHD_z_trajectory) * 0.5;
                        head_vector_x = MBHD_x_trajectory - C7_x_trajectory;
                        head_vector_z = MBHD_z_trajectory - C7_z_trajectory;
                        stretch_data = atan2(head_vector_x, head_vector_z);
                    end
                    if strcmp(variable_name, 'trunk_angle_ap')
                        % calculate angle trajectory
                        C7_y_trajectory = this.getTimeNormalizedData('marker:C7_y', this_stretch_times);
                        LPSI_y_trajectory = this.getTimeNormalizedData('marker:LPSI_y', this_stretch_times);
                        RPSI_y_trajectory = this.getTimeNormalizedData('marker:RPSI_y', this_stretch_times);
                        C7_z_trajectory = this.getTimeNormalizedData('marker:C7_z', this_stretch_times);
                        LPSI_z_trajectory = this.getTimeNormalizedData('marker:LPSI_z', this_stretch_times);
                        RPSI_z_trajectory = this.getTimeNormalizedData('marker:RPSI_z', this_stretch_times);
                        MPSI_y_trajectory = (LPSI_y_trajectory + RPSI_y_trajectory) * 0.5;
                        MPSI_z_trajectory = (LPSI_z_trajectory + RPSI_z_trajectory) * 0.5;
                        trunk_vector_y = C7_y_trajectory - MPSI_y_trajectory;
                        trunk_vector_z = C7_z_trajectory - MPSI_z_trajectory;
                        stretch_data = atan2(trunk_vector_y, trunk_vector_z);
                    end
                    if strcmp(variable_name, 'trunk_angle_ml')
                        % calculate angle trajectory
                        C7_x_trajectory = this.getTimeNormalizedData('marker:C7_x', this_stretch_times);
                        LPSI_x_trajectory = this.getTimeNormalizedData('marker:LPSI_x', this_stretch_times);
                        RPSI_x_trajectory = this.getTimeNormalizedData('marker:RPSI_x', this_stretch_times);
                        C7_z_trajectory = this.getTimeNormalizedData('marker:C7_z', this_stretch_times);
                        LPSI_z_trajectory = this.getTimeNormalizedData('marker:LPSI_z', this_stretch_times);
                        RPSI_z_trajectory = this.getTimeNormalizedData('marker:RPSI_z', this_stretch_times);
                        MPSI_x_trajectory = (LPSI_x_trajectory + RPSI_x_trajectory) * 0.5;
                        MPSI_z_trajectory = (LPSI_z_trajectory + RPSI_z_trajectory) * 0.5;
                        trunk_vector_x = C7_x_trajectory - MPSI_x_trajectory;
                        trunk_vector_z = C7_z_trajectory - MPSI_z_trajectory;
                        stretch_data = atan2(trunk_vector_x, trunk_vector_z);
                    end
                    
                    if strcmp(variable_name, 'left_foot_angle_ap')
                        % calculate angle trajectory
                        LTOE_x_trajectory = this.getTimeNormalizedData('marker:LTOE_x', this_stretch_times);
                        LTOEL_x_trajectory = this.getTimeNormalizedData('marker:LTOEL_x', this_stretch_times);
                        LHEE_x_trajectory = this.getTimeNormalizedData('marker:LHEE_x', this_stretch_times);
                        LTOE_y_trajectory = this.getTimeNormalizedData('marker:LTOE_y', this_stretch_times);
                        LTOEL_y_trajectory = this.getTimeNormalizedData('marker:LTOEL_y', this_stretch_times);
                        LHEE_y_trajectory = this.getTimeNormalizedData('marker:LHEE_y', this_stretch_times);
                        LTOE_z_trajectory = this.getTimeNormalizedData('marker:LTOE_z', this_stretch_times);
                        LTOEL_z_trajectory = this.getTimeNormalizedData('marker:LTOEL_z', this_stretch_times);
                        LHEE_z_trajectory = this.getTimeNormalizedData('marker:LHEE_z', this_stretch_times);
                        
                        LTOEM_x_trajectory = (LTOE_x_trajectory + LTOEL_x_trajectory) * 0.5;
                        LTOEM_y_trajectory = (LTOE_y_trajectory + LTOEL_y_trajectory) * 0.5;
                        LTOEM_z_trajectory = (LTOE_z_trajectory + LTOEL_z_trajectory) * 0.5;
                        
                        foot_vector_x = LTOEM_x_trajectory - LHEE_x_trajectory;
                        foot_vector_y = LTOEM_y_trajectory - LHEE_y_trajectory;
                        foot_vector_z = LTOEM_z_trajectory - LHEE_z_trajectory;
                        foot_vector_xy = (foot_vector_x.^2 + foot_vector_y.^2).^(0.5);
                        stretch_data = atan2(foot_vector_z, foot_vector_xy);
                    end
                    if strcmp(variable_name, 'left_foot_angle_ml')
                        % calculate angle trajectory
                        LTOE_x_trajectory = this.getTimeNormalizedData('marker:LTOE_x', this_stretch_times);
                        LTOEL_x_trajectory = this.getTimeNormalizedData('marker:LTOEL_x', this_stretch_times);
                        LTOE_y_trajectory = this.getTimeNormalizedData('marker:LTOE_y', this_stretch_times);
                        LTOEL_y_trajectory = this.getTimeNormalizedData('marker:LTOEL_y', this_stretch_times);
                        LTOE_z_trajectory = this.getTimeNormalizedData('marker:LTOE_z', this_stretch_times);
                        LTOEL_z_trajectory = this.getTimeNormalizedData('marker:LTOEL_z', this_stretch_times);
                        
                        foot_vector_x = LTOEL_x_trajectory - LTOE_x_trajectory;
                        foot_vector_y = LTOEL_y_trajectory - LTOE_y_trajectory;
                        foot_vector_z = LTOEL_z_trajectory - LTOE_z_trajectory;
                        foot_vector_xy = (foot_vector_x.^2 + foot_vector_y.^2).^(0.5);
                        stretch_data = atan2(foot_vector_z, foot_vector_xy);
                    end
                    if strcmp(variable_name, 'left_foot_angle_yaw')
                        % calculate angle trajectory
                        LTOE_x_trajectory = this.getTimeNormalizedData('marker:LTOE_x', this_stretch_times);
                        LTOEL_x_trajectory = this.getTimeNormalizedData('marker:LTOEL_x', this_stretch_times);
                        LHEE_x_trajectory = this.getTimeNormalizedData('marker:LHEE_x', this_stretch_times);
                        LTOE_y_trajectory = this.getTimeNormalizedData('marker:LTOE_y', this_stretch_times);
                        LTOEL_y_trajectory = this.getTimeNormalizedData('marker:LTOEL_y', this_stretch_times);
                        LHEE_y_trajectory = this.getTimeNormalizedData('marker:LHEE_y', this_stretch_times);
                        
                        LTOEM_x_trajectory = (LTOE_x_trajectory + LTOEL_x_trajectory) * 0.5;
                        LTOEM_y_trajectory = (LTOE_y_trajectory + LTOEL_y_trajectory) * 0.5;
                        
                        foot_vector_x = LTOEM_x_trajectory - LHEE_x_trajectory;
                        foot_vector_y = LTOEM_y_trajectory - LHEE_y_trajectory;
                        stretch_data = atan2(foot_vector_x, foot_vector_y);
                    end
                    if strcmp(variable_name, 'right_foot_angle_ap')
                        % calculate angle trajectory
                        RTOE_x_trajectory = this.getTimeNormalizedData('marker:RTOE_x', this_stretch_times);
                        RTOEL_x_trajectory = this.getTimeNormalizedData('marker:RTOEL_x', this_stretch_times);
                        RHEE_x_trajectory = this.getTimeNormalizedData('marker:RHEE_x', this_stretch_times);
                        RTOE_y_trajectory = this.getTimeNormalizedData('marker:RTOE_y', this_stretch_times);
                        RTOEL_y_trajectory = this.getTimeNormalizedData('marker:RTOEL_y', this_stretch_times);
                        RHEE_y_trajectory = this.getTimeNormalizedData('marker:RHEE_y', this_stretch_times);
                        RTOE_z_trajectory = this.getTimeNormalizedData('marker:RTOE_z', this_stretch_times);
                        RTOEL_z_trajectory = this.getTimeNormalizedData('marker:RTOEL_z', this_stretch_times);
                        RHEE_z_trajectory = this.getTimeNormalizedData('marker:RHEE_z', this_stretch_times);
                        
                        RTOEM_x_trajectory = (RTOE_x_trajectory + RTOEL_x_trajectory) * 0.5;
                        RTOEM_y_trajectory = (RTOE_y_trajectory + RTOEL_y_trajectory) * 0.5;
                        RTOEM_z_trajectory = (RTOE_z_trajectory + RTOEL_z_trajectory) * 0.5;
                        
                        foot_vector_x = RTOEM_x_trajectory - RHEE_x_trajectory;
                        foot_vector_y = RTOEM_y_trajectory - RHEE_y_trajectory;
                        foot_vector_z = RTOEM_z_trajectory - RHEE_z_trajectory;
                        foot_vector_xy = (foot_vector_x.^2 + foot_vector_y.^2).^(0.5);
                        stretch_data = atan2(foot_vector_z, foot_vector_xy);
                    end
                    if strcmp(variable_name, 'right_foot_angle_ml')
                        % calculate angle trajectory
                        RTOE_x_trajectory = this.getTimeNormalizedData('marker:RTOE_x', this_stretch_times);
                        RTOEL_x_trajectory = this.getTimeNormalizedData('marker:RTOEL_x', this_stretch_times);
                        RTOE_y_trajectory = this.getTimeNormalizedData('marker:RTOE_y', this_stretch_times);
                        RTOEL_y_trajectory = this.getTimeNormalizedData('marker:RTOEL_y', this_stretch_times);
                        RTOE_z_trajectory = this.getTimeNormalizedData('marker:RTOE_z', this_stretch_times);
                        RTOEL_z_trajectory = this.getTimeNormalizedData('marker:RTOEL_z', this_stretch_times);
                        
                        foot_vector_x = -(RTOEL_x_trajectory - RTOE_x_trajectory);
                        foot_vector_y = -(RTOEL_y_trajectory - RTOE_y_trajectory);
                        foot_vector_z = -(RTOEL_z_trajectory - RTOE_z_trajectory);
                        foot_vector_xy = (foot_vector_x.^2 + foot_vector_y.^2).^(0.5);
                        stretch_data = atan2(foot_vector_z, foot_vector_xy);
                    end
                    if strcmp(variable_name, 'right_foot_angle_yaw')
                        % calculate angle trajectory
                        RTOE_x_trajectory = this.getTimeNormalizedData('marker:RTOE_x', this_stretch_times);
                        RTOEL_x_trajectory = this.getTimeNormalizedData('marker:RTOEL_x', this_stretch_times);
                        RHEE_x_trajectory = this.getTimeNormalizedData('marker:RHEE_x', this_stretch_times);
                        RTOE_y_trajectory = this.getTimeNormalizedData('marker:RTOE_y', this_stretch_times);
                        RTOEL_y_trajectory = this.getTimeNormalizedData('marker:RTOEL_y', this_stretch_times);
                        RHEE_y_trajectory = this.getTimeNormalizedData('marker:RHEE_y', this_stretch_times);
                        
                        RTOEM_x_trajectory = (RTOE_x_trajectory + RTOEL_x_trajectory) * 0.5;
                        RTOEM_y_trajectory = (RTOE_y_trajectory + RTOEL_y_trajectory) * 0.5;
                        
                        foot_vector_x = RTOEM_x_trajectory - RHEE_x_trajectory;
                        foot_vector_y = RTOEM_y_trajectory - RHEE_y_trajectory;
                        stretch_data = atan2(foot_vector_x, foot_vector_y);
                    end
                    
                    % clearance above obstacle
                    if strcmp(variable_name, 'heel_clearance')
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stretch_data(i_band) = NaN;
                            else
                                % get relevant data
                                if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                    swing_marker_y_complete = this.getBasicVariableData('marker:LHEE_y');
                                    swing_marker_z_complete = this.getBasicVariableData('marker:LHEE_z');
                                    variable_time = this.getTimeData('marker:LHEE_y');
                                end
                                if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                    swing_marker_y_complete = this.getBasicVariableData('marker:RHEE_y');
                                    swing_marker_z_complete = this.getBasicVariableData('marker:RHEE_z');
                                    variable_time = this.getTimeData('marker:RHEE_y');
                                end
                                this_band_start_time = this_stretch_times(i_band);
                                this_band_end_time = this_stretch_times(i_band+1);
                                [~, start_index] = min(abs(variable_time - this_band_start_time));
                                [~, end_index] = min(abs(variable_time - this_band_end_time));
                                swing_heel_marker_y = swing_marker_y_complete(start_index : end_index);
                                swing_heel_marker_z = swing_marker_z_complete(start_index : end_index);

                                % find time step of obstacle crossing
                                obstacle_pos_y = NaN;
                                if strcmp(relevant_condition_data{i_stretch}, 'OBS_NEAR')
                                    obstacle_pos_y = this.subject_settings.get('near_distance');
                                end
                                if strcmp(relevant_condition_data{i_stretch}, 'OBS_FAR')
                                    obstacle_pos_y = this.subject_settings.get('far_distance');
                                end
                                [~, crossing_time_step] = min(abs(swing_heel_marker_y - obstacle_pos_y));

                                % extract swing foot heel marker
                                if strcmp(relevant_condition_data{i_stretch}, 'OBS_NO')
                                    stretch_data(i_band) = NaN;
                                else
                                    obstacle_height = this.subject_settings.get('obstacle_height');
                                    stretch_data(i_band) = swing_heel_marker_z(crossing_time_step) - obstacle_height;
                                end
                            end
                        end
                    end
                    if strcmp(variable_name, 'toes_clearance')
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stretch_data(i_band) = NaN;
                            else
                                % get relevant data
                                if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                    swing_marker_y_complete = this.getBasicVariableData('marker:LTOE_y');
                                    swing_marker_z_complete = this.getBasicVariableData('marker:LTOE_z');
                                    variable_time = this.getTimeData('marker:LTOE_y');
                                end
                                if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                    swing_marker_y_complete = this.getBasicVariableData('marker:RTOE_y');
                                    swing_marker_z_complete = this.getBasicVariableData('marker:RTOE_z');
                                    variable_time = this.getTimeData('marker:RTOE_y');
                                end
                                this_band_start_time = this_stretch_times(i_band);
                                this_band_end_time = this_stretch_times(i_band+1);
                                [~, start_index] = min(abs(variable_time - this_band_start_time));
                                [~, end_index] = min(abs(variable_time - this_band_end_time));
                                swing_toes_marker_y = swing_marker_y_complete(start_index : end_index);
                                swing_toes_marker_z = swing_marker_z_complete(start_index : end_index);

                                % find time step of obstacle crossing
                                obstacle_pos_y = NaN;
                                if strcmp(relevant_condition_data{i_stretch}, 'OBS_NEAR')
                                    obstacle_pos_y = this.subject_settings.get('toe_marker') + this.subject_settings.get('near_distance');
                                end
                                if strcmp(relevant_condition_data{i_stretch}, 'OBS_FAR')
                                    obstacle_pos_y = this.subject_settings.get('toe_marker') + this.subject_settings.get('far_distance');
                                end
                                [~, crossing_time_step] = min(abs(swing_toes_marker_y - obstacle_pos_y));

                                % extract swing foot toes marker
                                if strcmp(relevant_condition_data{i_stretch}, 'OBS_NO')
                                    stretch_data(i_band) = NaN;
                                else
                                    obstacle_height = this.subject_settings.get('obstacle_height');
                                    stretch_data(i_band) = swing_toes_marker_z(crossing_time_step) - obstacle_height;
                                end
                            end
                        end
                    end
                    
                    % store in cell
                    stretch_variables{i_variable} = [stretch_variables{i_variable} stretch_data];
                end
                
            end
        end
        function data_normalized = getTimeNormalizedData(this, variable_name, band_times)
            % extract data
            try
                variable_time = this.getTimeData(variable_name);
                variable_data = this.getBasicVariableData(variable_name);
                
                band_time_indices = zeros(size(band_times));
                for i_band_time = 1 : length(band_times)
                    [~, time_index] = min(abs(variable_time - band_times(i_band_time)));
                    band_time_indices(i_band_time) = time_index;
                end
                
                time_extracted = variable_time(band_time_indices(1) : band_time_indices(end));
                data_extracted = variable_data(band_time_indices(1) : band_time_indices(end));
                band_time_indices_local = band_time_indices - band_time_indices(1) + 1;
%                 if any(isnan(data_extracted))
%                     exception = MException('CoBaL:NaN', 'Data contains NaN values.');
%                     throw(exception)
%                 end
            catch error
                disp(['Error while processing variable ''' variable_name ''''])
                throw(error)
            end
                
            % normalize data in time
            if ~isempty(time_extracted) && ~any(isnan(data_extracted))
                % create normalized time
                number_of_bands = length(band_time_indices_local) - 1;
                time_normalized = [];
                for i_band = 1 : number_of_bands
                    time_normalized_this_band = linspace(time_extracted(band_time_indices_local(i_band)), time_extracted(band_time_indices_local(i_band+1)), this.number_of_time_steps_normalized)';
                    if i_band > 1
                        % start time of this band is end time of the last band, so remove the duplicate point
                        time_normalized_this_band = time_normalized_this_band(2:end);
                    end
                    time_normalized = [time_normalized; time_normalized_this_band]; %#ok<AGROW>
                end
                
                % time-normalize data
                data_normalized = spline(time_extracted, data_extracted, time_normalized);
            else
                number_of_bands = length(band_time_indices_local) - 1;
                time_normalized = [];
                for i_band = 1 : number_of_bands
                    time_normalized_this_band = linspace(time_extracted(band_time_indices_local(i_band)), time_extracted(band_time_indices_local(i_band+1)), this.number_of_time_steps_normalized)';
                    if i_band > 1
                        % start time of this band is end time of the last band, so remove the duplicate point
                        time_normalized_this_band = time_normalized_this_band(2:end);
                    end
                    time_normalized = [time_normalized; time_normalized_this_band]; %#ok<AGROW>
                end
                
                % time-normalize data
                data_normalized = time_normalized * NaN;
            end
        end
        
        function registerStretchVariableDirections(this, variable_name)
            this_variable_index = strcmp(variable_name, this.stretch_variable_names);
            
            % determine directions
            if this.isBasicVariable(variable_name)
                stretch_directions_new = this.basic_variable_directions.(variable_name);
            end

            % check if this is a compound name, listing a loaded variable and a label
            if any(variable_name==':')
                this_variable_split = strsplit(variable_name, ':');
                this_variable_type = this_variable_split{1};
                this_variable_label = this_variable_split{2};
                
                this_type_labels = this.basic_variable_labels.([this_variable_type '_trajectories']);
                this_type_directions = this.basic_variable_directions.([this_variable_type '_trajectories']);
                
                stretch_directions_new = this_type_directions(:, strcmp(this_type_labels, this_variable_label));
            end

            % use an initial value as a reference
            if strcmp(variable_name, 'lheel_from_mpsis_initial_x')
                [~, LHEE_x_directions] = this.getBasicVariableData('marker:LHEE_x');
                [~, LPSI_x_directions] = this.getBasicVariableData('marker:LPSI_x');
                [~, RPSI_x_directions] = this.getBasicVariableData('marker:RPSI_x');
                if ~strcmp(LHEE_x_directions{1}, LPSI_x_directions{1})
                    error('LHEE_x and LPSI_x directions are different from each other')
                end
                if ~strcmp(LHEE_x_directions{2}, LPSI_x_directions{2})
                    error('LHEE_x and LPSI_x directions are different from each other')
                end
                if ~strcmp(LHEE_x_directions{1}, RPSI_x_directions{1})
                    error('LHEE_x and RPSI_x directions are different from each other')
                end
                if ~strcmp(LHEE_x_directions{2}, RPSI_x_directions{2})
                    error('LHEE_x and RPSI_x directions are different from each other')
                end
                stretch_directions_new = LHEE_x_directions;
            end
            if strcmp(variable_name, 'rheel_from_mpsis_initial_x')
                [~, RHEE_x_directions] = this.getBasicVariableData('marker:RHEE_x');
                [~, LPSI_x_directions] = this.getBasicVariableData('marker:LPSI_x');
                [~, RPSI_x_directions] = this.getBasicVariableData('marker:RPSI_x');
                if ~strcmp(RHEE_x_directions{1}, LPSI_x_directions{1})
                    error('RHEE_x and LPSI_x directions are different from each other')
                end
                if ~strcmp(RHEE_x_directions{2}, LPSI_x_directions{2})
                    error('RHEE_x and LPSI_x directions are different from each other')
                end
                if ~strcmp(RHEE_x_directions{1}, RPSI_x_directions{1})
                    error('RHEE_x and RPSI_x directions are different from each other')
                end
                if ~strcmp(RHEE_x_directions{2}, RPSI_x_directions{2})
                    error('RHEE_x and RPSI_x directions are different from each other')
                end
                stretch_directions_new = RHEE_x_directions;
            end
            if strcmp(variable_name, 'mpsis_from_mpsis_initial_x')
                [~, LPSI_x_directions] = this.getBasicVariableData('marker:LPSI_x');
                [~, RPSI_x_directions] = this.getBasicVariableData('marker:RPSI_x');
                if ~strcmp(LPSI_x_directions{1}, RPSI_x_directions{1})
                    error('LPSI_x and RPSI_x directions are different from each other')
                end
                if ~strcmp(LPSI_x_directions{2}, RPSI_x_directions{2})
                    error('LPSI_x and RPSI_x directions are different from each other')
                end
                stretch_directions_new = LPSI_x_directions;
            end
            if strcmp(variable_name, 'lheel_from_com_initial_x')
                [~, LHEE_x_directions] = this.getBasicVariableData('marker:LHEE_x');
                [~, com_x_directions] = this.getBasicVariableData('com_position:center_of_mass_x');
                if ~strcmp(LHEE_x_directions{1}, com_x_directions{1})
                    error('LHEE_x and com_x directions are different from each other')
                end
                if ~strcmp(LHEE_x_directions{2}, com_x_directions{2})
                    error('LHEE_x and com_x directions are different from each other')
                end
                stretch_directions_new = LHEE_x_directions;
            end
            if strcmp(variable_name, 'rheel_from_com_initial_x')
                [~, RHEE_x_directions] = this.getBasicVariableData('marker:RHEE_x');
                [~, com_x_directions] = this.getBasicVariableData('com_position:center_of_mass_x');
                if ~strcmp(RHEE_x_directions{1}, com_x_directions{1})
                    error('RHEE_x and com_x directions are different from each other')
                end
                if ~strcmp(RHEE_x_directions{2}, com_x_directions{2})
                    error('RHEE_x and com_x directions are different from each other')
                end
                stretch_directions_new = RHEE_x_directions;
            end
            if strcmp(variable_name, 'com_from_com_initial_x')
                [~, com_x_directions] = this.getBasicVariableData('com_position:center_of_mass_x');
                stretch_directions_new = com_x_directions;
            end
            
            % balance stuff
            if strcmp(variable_name, 'step_placement_x')
                [~, LHEE_x_directions] = this.getBasicVariableData('marker:LHEE_x');
                [~, RHEE_x_directions] = this.getBasicVariableData('marker:RHEE_x');
                if ~strcmp(LHEE_x_directions{1}, RHEE_x_directions{1})
                    error('LHEE_x and RHEE_x directions are different from each other')
                end
                if ~strcmp(LHEE_x_directions{2}, RHEE_x_directions{2})
                    error('LHEE_x and RHEE_x directions are different from each other')
                end
                stretch_directions_new = LHEE_x_directions;
            end
            if strcmp(variable_name, 'xcom_x')
                [~, com_x_directions] = this.getBasicVariableData('com_position:center_of_mass_x');
                stretch_directions_new = com_x_directions;
            end
            if strcmp(variable_name, 'xcom_y')
                [~, com_y_directions] = this.getBasicVariableData('com_position:center_of_mass_y');
                stretch_directions_new = com_y_directions;
            end
            if strcmp(variable_name, 'mos_x')
                [~, com_x_directions] = this.getBasicVariableData('com_position:center_of_mass_x');
                stretch_directions_new = com_x_directions;
            end
            if strcmp(variable_name, 'mos_y')
                [~, com_y_directions] = this.getBasicVariableData('com_position:center_of_mass_y');
                stretch_directions_new = com_y_directions;
            end
            if strcmp(variable_name, 'xcom_mpsis_x')
                [~, LPSI_x_directions] = this.getBasicVariableData('marker:LPSI_x');
                [~, RPSI_x_directions] = this.getBasicVariableData('marker:RPSI_x');
                if ~strcmp(LPSI_x_directions{1}, RPSI_x_directions{1})
                    error('LPSI_x and RPSI_x directions are different from each other')
                end
                if ~strcmp(LPSI_x_directions{2}, RPSI_x_directions{2})
                    error('LPSI_x and RPSI_x directions are different from each other')
                end
                stretch_directions_new = LPSI_x_directions;
            end
            if strcmp(variable_name, 'xcom_mpsis_y')
                [~, LPSI_y_directions] = this.getBasicVariableData('marker:LPSI_y');
                [~, RPSI_y_directions] = this.getBasicVariableData('marker:RPSI_y');
                if ~strcmp(LPSI_y_directions{1}, RPSI_y_directions{1})
                    error('LPSI_y and RPSI_y directions are different from each other')
                end
                if ~strcmp(LPSI_y_directions{2}, RPSI_y_directions{2})
                    error('LPSI_y and RPSI_y directions are different from each other')
                end
                stretch_directions_new = LPSI_y_directions;
            end   
            if strcmp(variable_name, 'mos_mpsis_x')
                [~, LPSI_x_directions] = this.getBasicVariableData('marker:LPSI_x');
                [~, RPSI_x_directions] = this.getBasicVariableData('marker:RPSI_x');
                if ~strcmp(LPSI_x_directions{1}, RPSI_x_directions{1})
                    error('LPSI_x and RPSI_x directions are different from each other')
                end
                if ~strcmp(LPSI_x_directions{2}, RPSI_x_directions{2})
                    error('LPSI_x and RPSI_x directions are different from each other')
                end
                stretch_directions_new = LPSI_x_directions;
            end
            if strcmp(variable_name, 'mos_mpsis_y')
                [~, LPSI_y_directions] = this.getBasicVariableData('marker:LPSI_y');
                [~, RPSI_y_directions] = this.getBasicVariableData('marker:RPSI_y');
                if ~strcmp(LPSI_y_directions{1}, RPSI_y_directions{1})
                    error('LPSI_y and RPSI_y directions are different from each other')
                end
                if ~strcmp(LPSI_y_directions{2}, RPSI_y_directions{2})
                    error('LPSI_y and RPSI_y directions are different from each other')
                end
                stretch_directions_new = LPSI_y_directions;
            end
            if strcmp(variable_name, 'com_rough_x')
                [~, LPSI_x_directions] = this.getBasicVariableData('marker:LPSI_x');
                [~, RPSI_x_directions] = this.getBasicVariableData('marker:RPSI_x');
                [~, LASI_x_directions] = this.getBasicVariableData('marker:LASI_x');
                [~, RASI_x_directions] = this.getBasicVariableData('marker:RASI_x');
                [~, LKNE_x_directions] = this.getBasicVariableData('marker:LKNE_x');
                [~, RKNE_x_directions] = this.getBasicVariableData('marker:RKNE_x');
                [~, LANK_x_directions] = this.getBasicVariableData('marker:LANK_x');
                [~, RANK_x_directions] = this.getBasicVariableData('marker:RANK_x');
                [~, LTOE_x_directions] = this.getBasicVariableData('marker:LTOE_x');
                [~, RTOE_x_directions] = this.getBasicVariableData('marker:RTOE_x');
                [~, LFHD_x_directions] = this.getBasicVariableData('marker:LFHD_x');
                [~, RFHD_x_directions] = this.getBasicVariableData('marker:RFHD_x');
                [~, LBHD_x_directions] = this.getBasicVariableData('marker:LBHD_x');
                [~, RBHD_x_directions] = this.getBasicVariableData('marker:RBHD_x');
                [~, C7_x_directions] = this.getBasicVariableData('marker:C7_x');

                all_directions = ...
                  [
                    LPSI_x_directions, RPSI_x_directions, LASI_x_directions, RASI_x_directions, ...
                    LKNE_x_directions, RKNE_x_directions, LANK_x_directions, RANK_x_directions, LTOE_x_directions, RTOE_x_directions, ...
                    LFHD_x_directions, RFHD_x_directions, LBHD_x_directions, RBHD_x_directions, C7_x_directions ...
                  ]';
                unique_directions = unique(cell2table(all_directions), 'rows');
                if height(unique_directions) > 1
                    error('different directions found in marker data for rough CoM estimate')
                end
                stretch_directions_new = LPSI_x_directions;
            end
            if strcmp(variable_name, 'com_rough_y')
                [~, LPSI_y_directions] = this.getBasicVariableData('marker:LPSI_y');
                [~, RPSI_y_directions] = this.getBasicVariableData('marker:RPSI_y');
                [~, LASI_y_directions] = this.getBasicVariableData('marker:LASI_y');
                [~, RASI_y_directions] = this.getBasicVariableData('marker:RASI_y');
                [~, LKNE_y_directions] = this.getBasicVariableData('marker:LKNE_y');
                [~, RKNE_y_directions] = this.getBasicVariableData('marker:RKNE_y');
                [~, LANK_y_directions] = this.getBasicVariableData('marker:LANK_y');
                [~, RANK_y_directions] = this.getBasicVariableData('marker:RANK_y');
                [~, LTOE_y_directions] = this.getBasicVariableData('marker:LTOE_y');
                [~, RTOE_y_directions] = this.getBasicVariableData('marker:RTOE_y');
                [~, LFHD_y_directions] = this.getBasicVariableData('marker:LFHD_y');
                [~, RFHD_y_directions] = this.getBasicVariableData('marker:RFHD_y');
                [~, LBHD_y_directions] = this.getBasicVariableData('marker:LBHD_y');
                [~, RBHD_y_directions] = this.getBasicVariableData('marker:RBHD_y');
                [~, C7_y_directions] = this.getBasicVariableData('marker:C7_y');

                all_directions = ...
                  [
                    LPSI_y_directions, RPSI_y_directions, LASI_y_directions, RASI_y_directions, ...
                    LKNE_y_directions, RKNE_y_directions, LANK_y_directions, RANK_y_directions, LTOE_y_directions, RTOE_y_directions, ...
                    LFHD_y_directions, RFHD_y_directions, LBHD_y_directions, RBHD_y_directions, C7_y_directions ...
                  ]';
                unique_directions = unique(cell2table(all_directions), 'rows');
                if height(unique_directions) > 1
                    error('different directions found in marker data for rough CoM estimate')
                end
                stretch_directions_new = LPSI_y_directions;
            end
            if strcmp(variable_name, 'com_rough_z')
                [~, LPSI_z_directions] = this.getBasicVariableData('marker:LPSI_z');
                [~, RPSI_z_directions] = this.getBasicVariableData('marker:RPSI_z');
                [~, LASI_z_directions] = this.getBasicVariableData('marker:LASI_z');
                [~, RASI_z_directions] = this.getBasicVariableData('marker:RASI_z');
                [~, LKNE_z_directions] = this.getBasicVariableData('marker:LKNE_z');
                [~, RKNE_z_directions] = this.getBasicVariableData('marker:RKNE_z');
                [~, LANK_z_directions] = this.getBasicVariableData('marker:LANK_z');
                [~, RANK_z_directions] = this.getBasicVariableData('marker:RANK_z');
                [~, LTOE_z_directions] = this.getBasicVariableData('marker:LTOE_z');
                [~, RTOE_z_directions] = this.getBasicVariableData('marker:RTOE_z');
                [~, LFHD_z_directions] = this.getBasicVariableData('marker:LFHD_z');
                [~, RFHD_z_directions] = this.getBasicVariableData('marker:RFHD_z');
                [~, LBHD_z_directions] = this.getBasicVariableData('marker:LBHD_z');
                [~, RBHD_z_directions] = this.getBasicVariableData('marker:RBHD_z');
                [~, C7_z_directions] = this.getBasicVariableData('marker:C7_z');

                all_directions = ...
                  [
                    LPSI_z_directions, RPSI_z_directions, LASI_z_directions, RASI_z_directions, ...
                    LKNE_z_directions, RKNE_z_directions, LANK_z_directions, RANK_z_directions, LTOE_z_directions, RTOE_z_directions, ...
                    LFHD_z_directions, RFHD_z_directions, LBHD_z_directions, RBHD_z_directions, C7_z_directions ...
                  ]';
                unique_directions = unique(cell2table(all_directions), 'rows');
                if height(unique_directions) > 1
                    error('different directions found in marker data for rough CoM estimate')
                end
                stretch_directions_new = LPSI_z_directions;
            end
            if strcmp(variable_name, 'com_rough_x_vel')
                [~, LPSI_x_directions] = this.getBasicVariableData('derivative:LPSI_x_vel');
                [~, RPSI_x_directions] = this.getBasicVariableData('derivative:RPSI_x_vel');
                [~, LASI_x_directions] = this.getBasicVariableData('derivative:LASI_x_vel');
                [~, RASI_x_directions] = this.getBasicVariableData('derivative:RASI_x_vel');
                [~, LKNE_x_directions] = this.getBasicVariableData('derivative:LKNE_x_vel');
                [~, RKNE_x_directions] = this.getBasicVariableData('derivative:RKNE_x_vel');
                [~, LANK_x_directions] = this.getBasicVariableData('derivative:LANK_x_vel');
                [~, RANK_x_directions] = this.getBasicVariableData('derivative:RANK_x_vel');
                [~, LTOE_x_directions] = this.getBasicVariableData('derivative:LTOE_x_vel');
                [~, RTOE_x_directions] = this.getBasicVariableData('derivative:RTOE_x_vel');
                [~, LFHD_x_directions] = this.getBasicVariableData('derivative:LFHD_x_vel');
                [~, RFHD_x_directions] = this.getBasicVariableData('derivative:RFHD_x_vel');
                [~, LBHD_x_directions] = this.getBasicVariableData('derivative:LBHD_x_vel');
                [~, RBHD_x_directions] = this.getBasicVariableData('derivative:RBHD_x_vel');
                [~, C7_x_directions] = this.getBasicVariableData('derivative:C7_x_vel');

                all_directions = ...
                  [
                    LPSI_x_directions, RPSI_x_directions, LASI_x_directions, RASI_x_directions, ...
                    LKNE_x_directions, RKNE_x_directions, LANK_x_directions, RANK_x_directions, LTOE_x_directions, RTOE_x_directions, ...
                    LFHD_x_directions, RFHD_x_directions, LBHD_x_directions, RBHD_x_directions, C7_x_directions ...
                  ]';
                unique_directions = unique(cell2table(all_directions), 'rows');
                if height(unique_directions) > 1
                    error('different directions found in marker data for rough CoM estimate')
                end
                stretch_directions_new = LPSI_x_directions;
            end 
            if strcmp(variable_name, 'com_rough_y_vel')
                [~, LPSI_y_directions] = this.getBasicVariableData('derivative:LPSI_y_vel');
                [~, RPSI_y_directions] = this.getBasicVariableData('derivative:RPSI_y_vel');
                [~, LASI_y_directions] = this.getBasicVariableData('derivative:LASI_y_vel');
                [~, RASI_y_directions] = this.getBasicVariableData('derivative:RASI_y_vel');
                [~, LKNE_y_directions] = this.getBasicVariableData('derivative:LKNE_y_vel');
                [~, RKNE_y_directions] = this.getBasicVariableData('derivative:RKNE_y_vel');
                [~, LANK_y_directions] = this.getBasicVariableData('derivative:LANK_y_vel');
                [~, RANK_y_directions] = this.getBasicVariableData('derivative:RANK_y_vel');
                [~, LTOE_y_directions] = this.getBasicVariableData('derivative:LTOE_y_vel');
                [~, RTOE_y_directions] = this.getBasicVariableData('derivative:RTOE_y_vel');
                [~, LFHD_y_directions] = this.getBasicVariableData('derivative:LFHD_y_vel');
                [~, RFHD_y_directions] = this.getBasicVariableData('derivative:RFHD_y_vel');
                [~, LBHD_y_directions] = this.getBasicVariableData('derivative:LBHD_y_vel');
                [~, RBHD_y_directions] = this.getBasicVariableData('derivative:RBHD_y_vel');
                [~, C7_y_directions] = this.getBasicVariableData('derivative:C7_y_vel');

                all_directions = ...
                  [
                    LPSI_y_directions, RPSI_y_directions, LASI_y_directions, RASI_y_directions, ...
                    LKNE_y_directions, RKNE_y_directions, LANK_y_directions, RANK_y_directions, LTOE_y_directions, RTOE_y_directions, ...
                    LFHD_y_directions, RFHD_y_directions, LBHD_y_directions, RBHD_y_directions, C7_y_directions ...
                  ]';
                unique_directions = unique(cell2table(all_directions), 'rows');
                if height(unique_directions) > 1
                    error('different directions found in marker data for rough CoM estimate')
                end
                stretch_directions_new = LPSI_y_directions;
            end
            if strcmp(variable_name, 'xcom_rough_x')
                com_rough_x_directions = this.stretch_variable_directions(strcmp('com_rough_x', this.stretch_variable_names), :);
                com_vel_rough_x_directions = this.stretch_variable_directions(strcmp('com_rough_x_vel', this.stretch_variable_names), :);
                
                if ~strcmp(com_rough_x_directions{1}, com_vel_rough_x_directions{1})
                    error('com_rough_x and com_vel_rough_x directions are different from each other')
                end
                if ~strcmp(com_rough_x_directions{2}, com_vel_rough_x_directions{2})
                    error('com_rough_x and com_vel_rough_x directions are different from each other')
                end
                stretch_directions_new = com_rough_x_directions;
            end
            if strcmp(variable_name, 'xcom_rough_y')
                com_rough_y_directions = this.stretch_variable_directions(strcmp('com_rough_y', this.stretch_variable_names), :);
                com_vel_rough_y_directions = this.stretch_variable_directions(strcmp('com_rough_y_vel', this.stretch_variable_names), :);
                
                if ~strcmp(com_rough_y_directions{1}, com_vel_rough_y_directions{1})
                    error('com_rough_y and com_vel_rough_y directions are different from each other')
                end
                if ~strcmp(com_rough_y_directions{2}, com_vel_rough_y_directions{2})
                    error('com_rough_y and com_vel_rough_y directions are different from each other')
                end
                stretch_directions_new = com_rough_y_directions;
            end 
            if strcmp(variable_name, 'mos_rough_x')
                stretch_directions_new = this.stretch_variable_directions(strcmp('com_rough_x', this.stretch_variable_names), :);
            end
            if strcmp(variable_name, 'mos_rough_y')
                stretch_directions_new = this.stretch_variable_directions(strcmp('com_rough_y', this.stretch_variable_names), :);
            end
            
            % gait parameters
            if strcmp(variable_name, 'step_length')
                [~, LHEE_y_directions] = this.getBasicVariableData('marker:LHEE_y');
                [~, RHEE_y_directions] = this.getBasicVariableData('marker:RHEE_y');
                if ~strcmp(LHEE_y_directions{1}, RHEE_y_directions{1})
                    error('LHEE_y and RHEE_y directions are different from each other')
                end
                if ~strcmp(LHEE_y_directions{2}, RHEE_y_directions{2})
                    error('LHEE_y and RHEE_y directions are different from each other')
                end
                stretch_directions_new = LHEE_y_directions;
            end
            if strcmp(variable_name, 'step_width')
                [~, LHEE_x_directions] = this.getBasicVariableData('marker:LHEE_x');
                [~, RHEE_x_directions] = this.getBasicVariableData('marker:RHEE_x');
                if ~strcmp(LHEE_x_directions{1}, RHEE_x_directions{1})
                    error('LHEE_x and RHEE_x directions are different from each other')
                end
                if ~strcmp(LHEE_x_directions{2}, RHEE_x_directions{2})
                    error('LHEE_x and RHEE_x directions are different from each other')
                end
                stretch_directions_new = LHEE_x_directions;
            end
            if strcmp(variable_name, 'step_time')
                stretch_directions_new = {'+'; '-'};
            end
            if strcmp(variable_name, 'pushoff_time')
                stretch_directions_new = {'+'; '-'};
            end
            if strcmp(variable_name, 'midstance_index')
                stretch_directions_new = {'+'; '-'};
            end                        
            if strcmp(variable_name, 'midstance_time')
                stretch_directions_new = {'+'; '-'};
            end                        
            if strcmp(variable_name, 'midswing_time')
                stretch_directions_new = {'+'; '-'};
            end                        
            if strcmp(variable_name, 'cadence')
                stretch_directions_new = {'+'; '-'};
            end
            if strcmp(variable_name, 'velocity')
                stretch_directions_new = this.stretch_variable_directions(strcmp('step_length', this.stretch_variable_names), :);
            end
            if strcmp(variable_name, 'band_duration')
                stretch_directions_new = this.stretch_variable_directions(strcmp('step_time', this.stretch_variable_names), :);
            end
            
            % segment angles and lengths
            if strcmp(variable_name, 'leg_length_l')
                stretch_directions_new = {'+'; '-'};
            end
            if strcmp(variable_name, 'leg_length_r')
                stretch_directions_new = {'+'; '-'};
            end
            if strcmp(variable_name, 'leg_angle_ap_l')
                stretch_directions_new = {'top leans forward'; 'top leans backward'};
                % TODO: check assumptions, similar to what I do for other angles
            end
            if strcmp(variable_name, 'leg_angle_ap_r')
                stretch_directions_new = {'top leans forward'; 'top leans backward'};
                % TODO: check assumptions, similar to what I do for other angles
            end
            if strcmp(variable_name, 'head_angle_ap')
                % determine directions
                [~, C7_y_directions] = this.getBasicVariableData('marker:C7_y');
                [~, C7_z_directions] = this.getBasicVariableData('marker:C7_z');
                [~, LBHD_y_directions] = this.getBasicVariableData('marker:LBHD_y');
                [~, RBHD_y_directions] = this.getBasicVariableData('marker:RBHD_y');
                [~, LBHD_z_directions] = this.getBasicVariableData('marker:LBHD_z');
                [~, RBHD_z_directions] = this.getBasicVariableData('marker:RBHD_z');

                unique_directions_y = unique(cell2table([C7_y_directions, LBHD_y_directions, RBHD_y_directions]'), 'rows');
                if height(unique_directions_y) > 1
                    error('different directions found in marker data for head_angle_ap estimate')
                end
                unique_directions_z = unique(cell2table([C7_z_directions, LBHD_z_directions, RBHD_z_directions]'), 'rows');
                if height(unique_directions_z) > 1
                    error('different directions found in marker data for head_angle_ap estimate')
                end

                % check assumption that y is forward-backward and z is down-up
                if ~strcmp(LBHD_y_directions{1}, 'forward')
                    error('Assuming positive y-direction for markers C7, LBHD and RBHD is "down"')
                end                  
                if ~strcmp(LBHD_y_directions{2}, 'backward')
                    error('Assuming negative y-direction for markers C7, LBHD and RBHD is "backward"')
                end                  
                if ~strcmp(LBHD_z_directions{1}, 'up')
                    error('Assuming positive z-direction for markers C7, LBHD and RBHD is "up"')
                end                  
                if ~strcmp(RBHD_z_directions{2}, 'down')
                    error('Assuming negative z-direction for markers C7, LBHD and RBHD is "down"')
                end

                % all assumptions are met, define directions
                stretch_directions_new = {'forward'; 'backward'};
            end
            if strcmp(variable_name, 'head_angle_ml')
                [~, C7_x_directions] = this.getBasicVariableData('marker:C7_x');
                [~, C7_z_directions] = this.getBasicVariableData('marker:C7_z');
                [~, LBHD_x_directions] = this.getBasicVariableData('marker:LBHD_x');
                [~, RBHD_x_directions] = this.getBasicVariableData('marker:RBHD_x');
                [~, LBHD_z_directions] = this.getBasicVariableData('marker:LBHD_z');
                [~, RBHD_z_directions] = this.getBasicVariableData('marker:RBHD_z');
                unique_directions_x = unique(cell2table([C7_x_directions, LBHD_x_directions, RBHD_x_directions]'), 'rows');
                if height(unique_directions_x) > 1
                    error('different directions found in marker data for head_angle_ap estimate')
                end
                unique_directions_z = unique(cell2table([C7_z_directions, LBHD_z_directions, RBHD_z_directions]'), 'rows');
                if height(unique_directions_z) > 1
                    error('different directions found in marker data for head_angle_ap estimate')
                end
                
                % check assumption that x is left-right and z is down-up
                if ~strcmp(LBHD_x_directions{1}, 'right')
                    error('Assuming positive x-direction for markers C7, LBHD and RBHD is "right"')
                end                  
                if ~strcmp(LBHD_x_directions{2}, 'left')
                    error('Assuming negative x-direction for markers C7, LBHD and RBHD is "left"')
                end                  
                if ~strcmp(LBHD_z_directions{1}, 'up')
                    error('Assuming positive z-direction for markers C7, LBHD and RBHD is "up"')
                end                  
                if ~strcmp(LBHD_z_directions{2}, 'down')
                    error('Assuming negative z-direction for markers C7, LBHD and RBHD is "down"')
                end

                % all assumptions are met, define directions
                stretch_directions_new = {'right'; 'left'};                    
            end
            if strcmp(variable_name, 'trunk_angle_ap')
                % determine directions
                [~, C7_y_directions] = this.getBasicVariableData('marker:C7_y');
                [~, C7_z_directions] = this.getBasicVariableData('marker:C7_z');
                [~, LPSI_y_directions] = this.getBasicVariableData('marker:LPSI_y');
                [~, RPSI_y_directions] = this.getBasicVariableData('marker:RPSI_y');
                [~, LPSI_z_directions] = this.getBasicVariableData('marker:LPSI_z');
                [~, RPSI_z_directions] = this.getBasicVariableData('marker:RPSI_z');

                unique_directions_y = unique(cell2table([C7_y_directions, LPSI_y_directions, RPSI_y_directions]'), 'rows');
                if height(unique_directions_y) > 1
                    error('different directions found in marker data for trunk_angle_ap estimate')
                end
                unique_directions_z = unique(cell2table([C7_z_directions, LPSI_z_directions, RPSI_z_directions]'), 'rows');
                if height(unique_directions_z) > 1
                    error('different directions found in marker data for trunk_angle_ap estimate')
                end

                % check assumption that y is forward-backward and z is down-up
                if ~strcmp(LPSI_y_directions{1}, 'forward')
                    error('Assuming positive y-direction for markers C7, LPSI and RPSI is "down"')
                end                  
                if ~strcmp(LPSI_y_directions{2}, 'backward')
                    error('Assuming negative y-direction for markers C7, LPSI and RPSI is "backward"')
                end                  
                if ~strcmp(LPSI_z_directions{1}, 'up')
                    error('Assuming positive z-direction for markers C7, LPSI and RPSI is "up"')
                end                  
                if ~strcmp(RPSI_z_directions{2}, 'down')
                    error('Assuming negative z-direction for markers C7, LPSI and RPSI is "down"')
                end

                % all assumptions are met, define directions
                stretch_directions_new = {'forward'; 'backward'};
            end
            if strcmp(variable_name, 'trunk_angle_ml')
                [~, C7_x_directions] = this.getBasicVariableData('marker:C7_x');
                [~, C7_z_directions] = this.getBasicVariableData('marker:C7_z');
                [~, LPSI_x_directions] = this.getBasicVariableData('marker:LPSI_x');
                [~, RPSI_x_directions] = this.getBasicVariableData('marker:RPSI_x');
                [~, LPSI_z_directions] = this.getBasicVariableData('marker:LPSI_z');
                [~, RPSI_z_directions] = this.getBasicVariableData('marker:RPSI_z');
                unique_directions_x = unique(cell2table([C7_x_directions, LPSI_x_directions, RPSI_x_directions]'), 'rows');
                if height(unique_directions_x) > 1
                    error('different directions found in marker data for trunk_angle_ap estimate')
                end
                unique_directions_z = unique(cell2table([C7_z_directions, LPSI_z_directions, RPSI_z_directions]'), 'rows');
                if height(unique_directions_z) > 1
                    error('different directions found in marker data for trunk_angle_ap estimate')
                end
                
                % check assumption that x is left-right and z is down-up
                if ~strcmp(LPSI_x_directions{1}, 'right')
                    error('Assuming positive x-direction for markers C7, LPSI and RPSI is "right"')
                end                  
                if ~strcmp(LPSI_x_directions{2}, 'left')
                    error('Assuming negative x-direction for markers C7, LPSI and RPSI is "left"')
                end                  
                if ~strcmp(LPSI_z_directions{1}, 'up')
                    error('Assuming positive z-direction for markers C7, LPSI and RPSI is "up"')
                end                  
                if ~strcmp(LPSI_z_directions{2}, 'down')
                    error('Assuming negative z-direction for markers C7, LPSI and RPSI is "down"')
                end

                % all assumptions are met, define directions
                stretch_directions_new = {'right'; 'left'};                    
                
            end
            
            if strcmp(variable_name, 'left_foot_angle_ap')
                % determine directions
                [~, LTOE_x_directions] = this.getBasicVariableData('marker:LTOE_x');
                [~, LTOE_y_directions] = this.getBasicVariableData('marker:LTOE_y');
                [~, LTOE_z_directions] = this.getBasicVariableData('marker:LTOE_z');
                [~, LTOEL_x_directions] = this.getBasicVariableData('marker:LTOEL_x');
                [~, LTOEL_y_directions] = this.getBasicVariableData('marker:LTOEL_y');
                [~, LTOEL_z_directions] = this.getBasicVariableData('marker:LTOEL_z');
                [~, LHEE_x_directions] = this.getBasicVariableData('marker:LHEE_x');
                [~, LHEE_y_directions] = this.getBasicVariableData('marker:LHEE_y');
                [~, LHEE_z_directions] = this.getBasicVariableData('marker:LHEE_z');

                unique_directions_x = unique(cell2table([LTOE_x_directions, LTOEL_x_directions, LHEE_x_directions]'), 'rows');
                if height(unique_directions_x) > 1
                    error('different directions found in marker data for left_foot_angle_ap estimate')
                end
                unique_directions_y = unique(cell2table([LTOE_y_directions, LTOEL_y_directions, LHEE_y_directions]'), 'rows');
                if height(unique_directions_y) > 1
                    error('different directions found in marker data for left_foot_angle_ap estimate')
                end
                unique_directions_z = unique(cell2table([LTOE_z_directions, LTOEL_z_directions, LHEE_z_directions]'), 'rows');
                if height(unique_directions_z) > 1
                    error('different directions found in marker data for left_foot_angle_ap estimate')
                end

                % check assumption that y is forward-backward and z is down-up
                if ~strcmp(LTOEL_x_directions{1}, 'right')
                    error('Assuming positive y-direction for markers LTOE, LTOEL and LHEE is "down"')
                end                  
                if ~strcmp(LTOEL_x_directions{2}, 'left')
                    error('Assuming negative y-direction for markers LTOE, LTOEL and LHEE is "backward"')
                end                  
                if ~strcmp(LTOEL_y_directions{1}, 'forward')
                    error('Assuming positive y-direction for markers LTOE, LTOEL and LHEE is "down"')
                end                  
                if ~strcmp(LTOEL_y_directions{2}, 'backward')
                    error('Assuming negative y-direction for markers LTOE, LTOEL and LHEE is "backward"')
                end                  
                if ~strcmp(LTOEL_z_directions{1}, 'up')
                    error('Assuming positive z-direction for markers LTOE, LTOEL and LHEE is "up"')
                end                  
                if ~strcmp(LHEE_z_directions{2}, 'down')
                    error('Assuming negative z-direction for markers LTOE, LTOEL and LHEE is "down"')
                end

                % all assumptions are met, define directions
                stretch_directions_new = {'forward'; 'backward'};

            end
            if strcmp(variable_name, 'left_foot_angle_ml')
                % determine directions
                [~, LTOE_x_directions] = this.getBasicVariableData('marker:LTOE_x');
                [~, LTOE_y_directions] = this.getBasicVariableData('marker:LTOE_y');
                [~, LTOE_z_directions] = this.getBasicVariableData('marker:LTOE_z');
                [~, LTOEL_x_directions] = this.getBasicVariableData('marker:LTOEL_x');
                [~, LTOEL_y_directions] = this.getBasicVariableData('marker:LTOEL_y');
                [~, LTOEL_z_directions] = this.getBasicVariableData('marker:LTOEL_z');
                [~, LHEE_x_directions] = this.getBasicVariableData('marker:LHEE_x');
                [~, LHEE_y_directions] = this.getBasicVariableData('marker:LHEE_y');
                [~, LHEE_z_directions] = this.getBasicVariableData('marker:LHEE_z');

                unique_directions_x = unique(cell2table([LTOE_x_directions, LTOEL_x_directions, LHEE_x_directions]'), 'rows');
                if height(unique_directions_x) > 1
                    error('different directions found in marker data for left_foot_angle_ap estimate')
                end
                unique_directions_y = unique(cell2table([LTOE_y_directions, LTOEL_y_directions, LHEE_y_directions]'), 'rows');
                if height(unique_directions_y) > 1
                    error('different directions found in marker data for left_foot_angle_ap estimate')
                end
                unique_directions_z = unique(cell2table([LTOE_z_directions, LTOEL_z_directions, LHEE_z_directions]'), 'rows');
                if height(unique_directions_z) > 1
                    error('different directions found in marker data for left_foot_angle_ap estimate')
                end

                % check assumption that y is forward-backward and z is down-up
                if ~strcmp(LTOEL_x_directions{1}, 'right')
                    error('Assuming positive y-direction for markers LTOE, LTOEL and LHEE is "down"')
                end                  
                if ~strcmp(LTOEL_x_directions{2}, 'left')
                    error('Assuming negative y-direction for markers LTOE, LTOEL and LHEE is "backward"')
                end                  
                if ~strcmp(LTOEL_y_directions{1}, 'forward')
                    error('Assuming positive y-direction for markers LTOE, LTOEL and LHEE is "down"')
                end                  
                if ~strcmp(LTOEL_y_directions{2}, 'backward')
                    error('Assuming negative y-direction for markers LTOE, LTOEL and LHEE is "backward"')
                end                  
                if ~strcmp(LTOEL_z_directions{1}, 'up')
                    error('Assuming positive z-direction for markers LTOE, LTOEL and LHEE is "up"')
                end                  
                if ~strcmp(LHEE_z_directions{2}, 'down')
                    error('Assuming negative z-direction for markers LTOE, LTOEL and LHEE is "down"')
                end

                % all assumptions are met, define directions
                stretch_directions_new = {'right'; 'left'};                    
            end
            if strcmp(variable_name, 'left_foot_angle_yaw')
                % determine directions
                [~, LTOE_x_directions] = this.getBasicVariableData('marker:LTOE_x');
                [~, LTOE_y_directions] = this.getBasicVariableData('marker:LTOE_y');
                [~, LTOEL_x_directions] = this.getBasicVariableData('marker:LTOEL_x');
                [~, LTOEL_y_directions] = this.getBasicVariableData('marker:LTOEL_y');
                [~, LHEE_x_directions] = this.getBasicVariableData('marker:LHEE_x');
                [~, LHEE_y_directions] = this.getBasicVariableData('marker:LHEE_y');

                unique_directions_x = unique(cell2table([LTOE_x_directions, LTOEL_x_directions, LHEE_x_directions]'), 'rows');
                if height(unique_directions_x) > 1
                    error('different directions found in marker data for left_foot_angle_ap estimate')
                end
                unique_directions_y = unique(cell2table([LTOE_y_directions, LTOEL_y_directions, LHEE_y_directions]'), 'rows');
                if height(unique_directions_y) > 1
                    error('different directions found in marker data for left_foot_angle_ap estimate')
                end

                % check assumption that y is forward-backward and z is down-up
                if ~strcmp(LTOEL_x_directions{1}, 'right')
                    error('Assuming positive y-direction for markers LTOE, LTOEL and LHEE is "down"')
                end                  
                if ~strcmp(LTOEL_x_directions{2}, 'left')
                    error('Assuming negative y-direction for markers LTOE, LTOEL and LHEE is "backward"')
                end                  
                if ~strcmp(LTOEL_y_directions{1}, 'forward')
                    error('Assuming positive y-direction for markers LTOE, LTOEL and LHEE is "down"')
                end                  
                if ~strcmp(LTOEL_y_directions{2}, 'backward')
                    error('Assuming negative y-direction for markers LTOE, LTOEL and LHEE is "backward"')
                end                  
                stretch_directions_new = {'clockwise'; 'counterclockwise'};
            end
            if strcmp(variable_name, 'right_foot_angle_ap')
                % determine directions
                [~, RTOE_x_directions] = this.getBasicVariableData('marker:RTOE_x');
                [~, RTOE_y_directions] = this.getBasicVariableData('marker:RTOE_y');
                [~, RTOE_z_directions] = this.getBasicVariableData('marker:RTOE_z');
                [~, RTOEL_x_directions] = this.getBasicVariableData('marker:RTOEL_x');
                [~, RTOEL_y_directions] = this.getBasicVariableData('marker:RTOEL_y');
                [~, RTOEL_z_directions] = this.getBasicVariableData('marker:RTOEL_z');
                [~, RHEE_x_directions] = this.getBasicVariableData('marker:RHEE_x');
                [~, RHEE_y_directions] = this.getBasicVariableData('marker:RHEE_y');
                [~, RHEE_z_directions] = this.getBasicVariableData('marker:RHEE_z');

                unique_directions_x = unique(cell2table([RTOE_x_directions, RTOEL_x_directions, RHEE_x_directions]'), 'rows');
                if height(unique_directions_x) > 1
                    error('different directions found in marker data for right_foot_angle_ap estimate')
                end
                unique_directions_y = unique(cell2table([RTOE_y_directions, RTOEL_y_directions, RHEE_y_directions]'), 'rows');
                if height(unique_directions_y) > 1
                    error('different directions found in marker data for right_foot_angle_ap estimate')
                end
                unique_directions_z = unique(cell2table([RTOE_z_directions, RTOEL_z_directions, RHEE_z_directions]'), 'rows');
                if height(unique_directions_z) > 1
                    error('different directions found in marker data for right_foot_angle_ap estimate')
                end

                % check assumption that y is forward-backward and z is down-up
                if ~strcmp(RTOEL_x_directions{1}, 'right')
                    error('Assuming positive y-direction for markers RTOE, RTOEL and RHEE is "down"')
                end                  
                if ~strcmp(RTOEL_x_directions{2}, 'left')
                    error('Assuming negative y-direction for markers RTOE, RTOEL and RHEE is "backward"')
                end                  
                if ~strcmp(RTOEL_y_directions{1}, 'forward')
                    error('Assuming positive y-direction for markers RTOE, RTOEL and RHEE is "down"')
                end                  
                if ~strcmp(RTOEL_y_directions{2}, 'backward')
                    error('Assuming negative y-direction for markers RTOE, RTOEL and RHEE is "backward"')
                end                  
                if ~strcmp(RTOEL_z_directions{1}, 'up')
                    error('Assuming positive z-direction for markers RTOE, RTOEL and RHEE is "up"')
                end                  
                if ~strcmp(RHEE_z_directions{2}, 'down')
                    error('Assuming negative z-direction for markers RTOE, RTOEL and RHEE is "down"')
                end

                % all assumptions are met, define directions
                stretch_directions_new = {'forward'; 'backward'};

            end
            if strcmp(variable_name, 'right_foot_angle_ml')
                % determine directions
                [~, RTOE_x_directions] = this.getBasicVariableData('marker:RTOE_x');
                [~, RTOE_y_directions] = this.getBasicVariableData('marker:RTOE_y');
                [~, RTOE_z_directions] = this.getBasicVariableData('marker:RTOE_z');
                [~, RTOEL_x_directions] = this.getBasicVariableData('marker:RTOEL_x');
                [~, RTOEL_y_directions] = this.getBasicVariableData('marker:RTOEL_y');
                [~, RTOEL_z_directions] = this.getBasicVariableData('marker:RTOEL_z');
                [~, RHEE_x_directions] = this.getBasicVariableData('marker:RHEE_x');
                [~, RHEE_y_directions] = this.getBasicVariableData('marker:RHEE_y');
                [~, RHEE_z_directions] = this.getBasicVariableData('marker:RHEE_z');

                unique_directions_x = unique(cell2table([RTOE_x_directions, RTOEL_x_directions, RHEE_x_directions]'), 'rows');
                if height(unique_directions_x) > 1
                    error('different directions found in marker data for right_foot_angle_ap estimate')
                end
                unique_directions_y = unique(cell2table([RTOE_y_directions, RTOEL_y_directions, RHEE_y_directions]'), 'rows');
                if height(unique_directions_y) > 1
                    error('different directions found in marker data for right_foot_angle_ap estimate')
                end
                unique_directions_z = unique(cell2table([RTOE_z_directions, RTOEL_z_directions, RHEE_z_directions]'), 'rows');
                if height(unique_directions_z) > 1
                    error('different directions found in marker data for right_foot_angle_ap estimate')
                end

                % check assumption that y is forward-backward and z is down-up
                if ~strcmp(RTOEL_x_directions{1}, 'right')
                    error('Assuming positive y-direction for markers RTOE, RTOEL and RHEE is "down"')
                end                  
                if ~strcmp(RTOEL_x_directions{2}, 'left')
                    error('Assuming negative y-direction for markers RTOE, RTOEL and RHEE is "backward"')
                end                  
                if ~strcmp(RTOEL_y_directions{1}, 'forward')
                    error('Assuming positive y-direction for markers RTOE, RTOEL and RHEE is "down"')
                end                  
                if ~strcmp(RTOEL_y_directions{2}, 'backward')
                    error('Assuming negative y-direction for markers RTOE, RTOEL and RHEE is "backward"')
                end                  
                if ~strcmp(RTOEL_z_directions{1}, 'up')
                    error('Assuming positive z-direction for markers RTOE, RTOEL and RHEE is "up"')
                end                  
                if ~strcmp(RHEE_z_directions{2}, 'down')
                    error('Assuming negative z-direction for markers RTOE, RTOEL and RHEE is "down"')
                end

                % all assumptions are met, define directions
                stretch_directions_new = {'right'; 'left'};                    
            end
            if strcmp(variable_name, 'right_foot_angle_yaw')
                % determine directions
                [~, RTOE_x_directions] = this.getBasicVariableData('marker:RTOE_x');
                [~, RTOE_y_directions] = this.getBasicVariableData('marker:RTOE_y');
                [~, RTOEL_x_directions] = this.getBasicVariableData('marker:RTOEL_x');
                [~, RTOEL_y_directions] = this.getBasicVariableData('marker:RTOEL_y');
                [~, RHEE_x_directions] = this.getBasicVariableData('marker:RHEE_x');
                [~, RHEE_y_directions] = this.getBasicVariableData('marker:RHEE_y');

                unique_directions_x = unique(cell2table([RTOE_x_directions, RTOEL_x_directions, RHEE_x_directions]'), 'rows');
                if height(unique_directions_x) > 1
                    error('different directions found in marker data for right_foot_angle_ap estimate')
                end
                unique_directions_y = unique(cell2table([RTOE_y_directions, RTOEL_y_directions, RHEE_y_directions]'), 'rows');
                if height(unique_directions_y) > 1
                    error('different directions found in marker data for right_foot_angle_ap estimate')
                end

                % check assumption that y is forward-backward and z is down-up
                if ~strcmp(RTOEL_x_directions{1}, 'right')
                    error('Assuming positive y-direction for markers RTOE, RTOEL and RHEE is "down"')
                end                  
                if ~strcmp(RTOEL_x_directions{2}, 'left')
                    error('Assuming negative y-direction for markers RTOE, RTOEL and RHEE is "backward"')
                end                  
                if ~strcmp(RTOEL_y_directions{1}, 'forward')
                    error('Assuming positive y-direction for markers RTOE, RTOEL and RHEE is "down"')
                end                  
                if ~strcmp(RTOEL_y_directions{2}, 'backward')
                    error('Assuming negative y-direction for markers RTOE, RTOEL and RHEE is "backward"')
                end                  
                stretch_directions_new = {'clockwise'; 'counterclockwise'};
            end
            
            % obstacle clearance
            if strcmp(variable_name, 'heel_clearance')
                [~, LHEE_z_directions] = this.getBasicVariableData('marker:LHEE_z');
                [~, RHEE_z_directions] = this.getBasicVariableData('marker:RHEE_z');
                if ~strcmp(LHEE_z_directions{1}, RHEE_z_directions{1})
                    error('LHEE_z and RHEE_z directions are different from each other')
                end
                if ~strcmp(LHEE_z_directions{2}, RHEE_z_directions{2})
                    error('LHEE_z and RHEE_z directions are different from each other')
                end
                stretch_directions_new = LHEE_z_directions;
            end
            if strcmp(variable_name, 'toes_clearance')
                [~, LTOE_z_directions] = this.getBasicVariableData('marker:LTOE_z');
                [~, RTOE_z_directions] = this.getBasicVariableData('marker:RTOE_z');
                if ~strcmp(LTOE_z_directions{1}, RTOE_z_directions{1})
                    error('LTOE_z and RTOE_z directions are different from each other')
                end
                if ~strcmp(LTOE_z_directions{2}, RTOE_z_directions{2})
                    error('LTOE_z and RTOE_z directions are different from each other')
                end
                stretch_directions_new = LTOE_z_directions;
            end
            
            if any(strcmp(this.basic_variable_load_failures, variable_name))
                stretch_directions_new = {'~'; '~'};
            end
            
            % compare against what is already on file
            stretch_directions_on_file = this.stretch_variable_directions(this_variable_index, :);
            if strcmp(stretch_directions_on_file{1}, 'TBD') && strcmp(stretch_directions_on_file{2}, 'TBD')
                % nothing is no file yet, so file the new information
                this.stretch_variable_directions(this_variable_index, :) = stretch_directions_new;
            else
                % check whether the new information matches up with what is on file
                if ~strcmp(stretch_directions_new{1}, stretch_directions_on_file{1})
                    error(['Different trials have different direction information for variable ' variable_name])
                end
                if ~strcmp(stretch_directions_new{2}, stretch_directions_on_file{2})
                    error(['Different trials have different direction information for variable ' variable_name])
                end
            end
        end
    end
end





