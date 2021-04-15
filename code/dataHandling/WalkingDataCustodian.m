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
                
            end

            % kinematics
            if this.isVariableToAnalyze('lheel_from_mpsis_initial_x')
                this.addBasicVariable('marker_trajectories')
                this.addStretchVariable('lheel_from_mpsis_initial_x')                
            end
            if this.isVariableToAnalyze('rheel_from_mpsis_initial_x')
                this.addBasicVariable('marker_trajectories')
                this.addStretchVariable('rheel_from_mpsis_initial_x')
            end
            if this.isVariableToAnalyze('mpsis_from_mpsis_initial_x')
                this.addBasicVariable('marker_trajectories')
                this.addStretchVariable('mpsis_from_mpsis_initial_x')
            end            
            if this.isVariableToAnalyze('lheel_from_com_initial_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('com_position_trajectories')
                this.addStretchVariable('lheel_from_com_initial_x')                
            end
            if this.isVariableToAnalyze('rheel_from_com_initial_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('com_position_trajectories')
                this.addStretchVariable('rheel_from_com_initial_x')
            end
            if this.isVariableToAnalyze('com_from_com_initial_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('com_position_trajectories')
                this.addStretchVariable('com_from_com_initial_x')
            end
            if this.isVariableToAnalyze('xcom_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('com_position_trajectories')
                this.addBasicVariable('com_velocity_trajectories')
                this.addStretchVariable('xcom_x')
            end
            if this.isVariableToAnalyze('xcom_y')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('com_position_trajectories')
                this.addBasicVariable('com_velocity_trajectories')
                this.addStretchVariable('xcom_y')
            end
            if this.isVariableToAnalyze('xcom_mpsis_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('derivative:LPSI_x_vel')
                this.addBasicVariable('derivative:RPSI_x_vel')
                this.addStretchVariable('xcom_mpsis_x')
            end
            if this.isVariableToAnalyze('xcom_mpsis_y')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('derivative:LPSI_y_vel')
                this.addBasicVariable('derivative:RPSI_y_vel')
                this.addStretchVariable('xcom_mpsis_y')
            end               
            if this.isVariableToAnalyze('xcom_rough_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('com_rough_x')
                this.addBasicVariable('com_rough_y')
                this.addBasicVariable('com_rough_z')
                this.addBasicVariable('com_rough_x_vel')
                this.addStretchVariable('xcom_rough_x')
            end  
            if this.isVariableToAnalyze('xcom_rough_y')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('com_rough_x')
                this.addBasicVariable('com_rough_y')
                this.addBasicVariable('com_rough_z')
                this.addBasicVariable('com_rough_y_vel')
                this.addStretchVariable('xcom_rough_y')
            end                         
            if this.isVariableToAnalyze('mos_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('com_position_trajectories')
                this.addBasicVariable('com_velocity_trajectories')
                this.addStretchVariable('xcom_x')
                this.addStretchVariable('mos_x')
            end
            if this.isVariableToAnalyze('mos_y')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('com_position_trajectories')
                this.addBasicVariable('com_velocity_trajectories')
                this.addStretchVariable('xcom_y')
                this.addStretchVariable('mos_y')
            end
            if this.isVariableToAnalyze('mos_mpsis_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('derivative:LPSI_x_vel')
                this.addBasicVariable('derivative:RPSI_x_vel')
                this.addStretchVariable('xcom_mpsis_x')
                this.addStretchVariable('mos_mpsis_x')
            end
            if this.isVariableToAnalyze('mos_mpsis_y')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('derivative:LPSI_y_vel')
                this.addBasicVariable('derivative:RPSI_y_vel')
                this.addStretchVariable('xcom_mpsis_y')
                this.addStretchVariable('mos_mpsis_y')
            end
            if this.isVariableToAnalyze('mos_rough_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('com_rough_x')
                this.addBasicVariable('com_rough_y')
                this.addBasicVariable('com_rough_z')
                this.addBasicVariable('com_rough_x_vel')
                this.addStretchVariable('xcom_rough_x')
                this.addStretchVariable('mos_rough_x')
            end
            if this.isVariableToAnalyze('mos_rough_y')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('com_rough_x')
                this.addBasicVariable('com_rough_y')
                this.addBasicVariable('com_rough_z')
                this.addBasicVariable('com_rough_y_vel')
                this.addStretchVariable('xcom_rough_y')
                this.addStretchVariable('mos_rough_y')
            end
            if this.isVariableToAnalyze('step_length')
                this.addBasicVariable('marker_trajectories')
                this.addStretchVariable('step_length')
            end
            if this.isVariableToAnalyze('step_width')
                this.addBasicVariable('marker_trajectories')
                this.addStretchVariable('step_width')
            end
            if this.isVariableToAnalyze('step_placement_x')
                this.addBasicVariable('marker_trajectories')
                this.addStretchVariable('step_placement_x')
            end
            if this.isVariableToAnalyze('stimulus_response_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('pelvis_y')
                this.addStretchVariable('step_time')
                this.addStretchVariable('pushoff_time')
                this.addStretchVariable('midstance_index')
                this.addStretchVariable('step_placement_x')
            end
            if this.isVariableToAnalyze('step_time')
                this.addStretchVariable('step_time')
            end
            if this.isVariableToAnalyze('step_duration')
                this.addStretchVariable('step_time')
                this.addStretchVariable('step_duration')
            end
            if this.isVariableToAnalyze('leg_length_l')
                this.addBasicVariable('marker_trajectories')
                this.addStretchVariable('leg_length_l')
            end
            if this.isVariableToAnalyze('leg_length_r')
                this.addBasicVariable('marker_trajectories')
                this.addStretchVariable('leg_length_r')
            end
            if this.isVariableToAnalyze('leg_angle_l')
                this.addBasicVariable('marker_trajectories')
                this.addStretchVariable('leg_angle_l')
            end
            if this.isVariableToAnalyze('leg_angle_r')
                this.addBasicVariable('marker_trajectories')
                this.addStretchVariable('leg_angle_r')
            end
            if this.isVariableToAnalyze('pushoff_time')
                this.addBasicVariable('marker_trajectories')
                this.addStretchVariable('step_time')
                this.addStretchVariable('pushoff_time')
            end
            if this.isVariableToAnalyze('midswing_event_time')
                this.addStretchVariable('midswing_event_time')
            end
            if this.isVariableToAnalyze('midstance_index')
                this.addBasicVariable('pelvis_y')
                this.addStretchVariable('step_time')
                this.addStretchVariable('pushoff_time')
                this.addStretchVariable('midstance_index')
            end
            if this.isVariableToAnalyze('midstance_time')
                this.addBasicVariable('pelvis_y')
                this.addStretchVariable('step_time')
                this.addStretchVariable('pushoff_time')
                this.addStretchVariable('midstance_time')
            end
            if this.isVariableToAnalyze('midswing_time')
                this.addStretchVariable('step_time')
                this.addStretchVariable('pushoff_time')
                this.addStretchVariable('midswing_time')
            end
            if this.isVariableToAnalyze('cadence')
                this.addStretchVariable('step_time')
                this.addStretchVariable('cadence')
            end
            if this.isVariableToAnalyze('velocity')
                this.addBasicVariable('marker_trajectories')
                this.addStretchVariable('step_length')
                this.addStretchVariable('step_time')
                this.addStretchVariable('velocity')
            end
            if this.isVariableToAnalyze('com_rough_x')
                this.addBasicVariable('com_rough_x')
                this.addStretchVariable('com_rough_x')
            end
            if this.isVariableToAnalyze('com_rough_y')
                this.addBasicVariable('com_rough_y')
                this.addStretchVariable('com_rough_y')
            end  
            if this.isVariableToAnalyze('com_rough_z')
                this.addBasicVariable('com_rough_z')
                this.addStretchVariable('com_rough_z')
            end
            if this.isVariableToAnalyze('com_rough_x_vel')
                this.addBasicVariable('com_rough_x')
                this.addBasicVariable('com_rough_x_vel')
                this.addStretchVariable('com_rough_x_vel')
            end  
            if this.isVariableToAnalyze('com_rough_y_vel')
                this.addBasicVariable('com_rough_y')
                this.addBasicVariable('com_rough_y_vel')
                this.addStretchVariable('com_rough_y_vel')
            end 
            if this.isVariableToAnalyze('head_angle_ap')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('head_angle_ap')
                this.addStretchVariable('head_angle_ap')
            end
            if this.isVariableToAnalyze('head_angle_ml')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('head_angle_ml')
                this.addStretchVariable('head_angle_ml')
            end
            if this.isVariableToAnalyze('trunk_angle_ap')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('trunk_angle_ap')
                this.addStretchVariable('trunk_angle_ap')
            end
            if this.isVariableToAnalyze('trunk_angle_ml')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('trunk_angle_ml')
                this.addStretchVariable('trunk_angle_ml')
            end
            if this.isVariableToAnalyze('pelvis_angle_ml')
                this.addBasicVariable('joint_center_trajectories')
                this.addBasicVariable('pelvis_angle_ml')
                this.addStretchVariable('pelvis_angle_ml')
            end
            if this.isVariableToAnalyze('left_leg_angle_ml')
                this.addBasicVariable('joint_center_trajectories')
                this.addBasicVariable('left_leg_angle_ml')
                this.addStretchVariable('left_leg_angle_ml')
            end
            if this.isVariableToAnalyze('right_leg_angle_ml')
                this.addBasicVariable('joint_center_trajectories')
                this.addBasicVariable('right_leg_angle_ml')
                this.addStretchVariable('right_leg_angle_ml')
            end           
            if this.isVariableToAnalyze('left_leg_angle_ap')
                this.addBasicVariable('joint_center_trajectories')
                this.addBasicVariable('left_leg_angle_ap')
                this.addStretchVariable('left_leg_angle_ap')
            end
            if this.isVariableToAnalyze('right_leg_angle_ap')
                this.addBasicVariable('joint_center_trajectories')
                this.addBasicVariable('right_leg_angle_ap')
                this.addStretchVariable('right_leg_angle_ap')
            end           
            if this.isVariableToAnalyze('left_foot_angle_ap')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('left_foot_angle_ap')
                this.addStretchVariable('left_foot_angle_ap')
            end
            if this.isVariableToAnalyze('left_foot_angle_ml')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('left_foot_angle_ml')
                this.addStretchVariable('left_foot_angle_ml')
            end
            if this.isVariableToAnalyze('left_foot_angle_yaw')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('left_foot_angle_yaw')
                this.addStretchVariable('left_foot_angle_yaw')
            end
            if this.isVariableToAnalyze('right_foot_angle_ap')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('right_foot_angle_ap')
                this.addStretchVariable('right_foot_angle_ap')
            end
            if this.isVariableToAnalyze('right_foot_angle_ml')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('right_foot_angle_ml')
                this.addStretchVariable('right_foot_angle_ml')
            end            
            if this.isVariableToAnalyze('right_foot_angle_yaw')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('right_foot_angle_yaw')
                this.addStretchVariable('right_foot_angle_yaw')
            end
            if this.isVariableToAnalyze('left_arm_angle_ap')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('left_arm_angle_ap')
                this.addStretchVariable('left_arm_angle_ap')
            end
            if this.isVariableToAnalyze('right_arm_angle_ap')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('right_arm_angle_ap')
                this.addStretchVariable('right_arm_angle_ap')
            end
            if this.isVariableToAnalyze('left_leg_angle_ap')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('left_leg_angle_ap')
                this.addStretchVariable('left_leg_angle_ap')
            end
            if this.isVariableToAnalyze('right_leg_angle_ap')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('right_leg_angle_ap')
                this.addStretchVariable('right_leg_angle_ap')
            end
            if this.isVariableToAnalyze('left_arm_phase')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('left_arm_angle_ap')
                this.addBasicVariable('left_arm_phase')
                this.addStretchVariable('left_arm_angle_ap')
                this.addStretchVariable('left_arm_phase')
            end
            if this.isVariableToAnalyze('right_arm_phase')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('right_arm_angle_ap')
                this.addBasicVariable('right_arm_phase')
                this.addStretchVariable('right_arm_angle_ap')
                this.addStretchVariable('right_arm_phase')
            end
            if this.isVariableToAnalyze('left_leg_phase')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('left_leg_angle_ap')
                this.addBasicVariable('left_leg_phase')
                this.addStretchVariable('left_leg_angle_ap')
                this.addStretchVariable('left_leg_phase')
            end
            if this.isVariableToAnalyze('right_leg_phase')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('right_leg_angle_ap')
                this.addBasicVariable('right_leg_phase')
                this.addStretchVariable('right_leg_angle_ap')
                this.addStretchVariable('right_leg_phase')
            end
            if this.isVariableToAnalyze('left_arm_right_leg_relative_phase')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('left_arm_angle_ap')
                this.addBasicVariable('left_arm_phase')
                this.addBasicVariable('right_leg_angle_ap')
                this.addBasicVariable('right_leg_phase')
                this.addBasicVariable('left_arm_right_leg_relative_phase')
                this.addStretchVariable('left_arm_angle_ap')
                this.addStretchVariable('left_arm_phase')
                this.addStretchVariable('right_leg_angle_ap')
                this.addStretchVariable('right_leg_phase')
                this.addStretchVariable('left_arm_right_leg_relative_phase')
            end
            if this.isVariableToAnalyze('right_arm_left_leg_relative_phase')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('right_arm_angle_ap')
                this.addBasicVariable('right_arm_phase')
                this.addBasicVariable('left_leg_angle_ap')
                this.addBasicVariable('left_leg_phase')
                this.addBasicVariable('right_arm_left_leg_relative_phase')
                this.addStretchVariable('right_arm_angle_ap')
                this.addStretchVariable('right_arm_phase')
                this.addStretchVariable('left_leg_angle_ap')
                this.addStretchVariable('left_leg_phase')
                this.addStretchVariable('right_arm_left_leg_relative_phase')
            end
            if this.isVariableToAnalyze('left_arm_phase_at_heelstrike')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('left_arm_angle_ap')
                this.addBasicVariable('left_arm_phase')
                this.addStretchVariable('left_arm_phase_at_heelstrike')
            end
            if this.isVariableToAnalyze('right_arm_phase_at_heelstrike')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('right_arm_angle_ap')
                this.addBasicVariable('right_arm_phase')
                this.addStretchVariable('right_arm_phase_at_heelstrike')
            end
            if this.isVariableToAnalyze('heel_clearance')
                this.addBasicVariable('marker_trajectories')
                this.addStretchVariable('heel_clearance')
            end
            if this.isVariableToAnalyze('toes_clearance')
                this.addBasicVariable('marker_trajectories')
                this.addStretchVariable('toes_clearance')
            end
            
            % protocol
            if this.isVariableToAnalyze('stimulus_state_trajectory')
                this.addBasicVariable('stimulus_state_trajectory')
                this.addStretchVariable('stimulus_state_trajectory')
            end
        end
        
        % interface
        function result = isVariableToAnalyze(this, variable_name)
            result = any(strcmp(this.variables_to_analyze, variable_name));
        end
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

            % prepare the data by loading all the basic variables from disk and calculating the required variables
            load(['analysis' filesep makeFileName(this.date, this.subject_id, trial_type, trial_number, 'availableVariables')], 'available_variables');

            % load basic variables
            for i_variable = 1 : length(variables_to_prepare)
                variable_name = variables_to_prepare{i_variable};
                
                % try loading
                [data, time, sampling_rate, labels, directions, success] = loadData(this.date, this.subject_id, trial_type, trial_number, variable_name, 'optional'); %#ok<ASGLU>
                
                % store
                if success
                    eval(['this.basic_variable_data.' variable_name ' = data;']);
                    eval(['this.time_data.' variable_name ' = time;']);
                    eval(['this.basic_variable_labels.' variable_name ' = labels;']);
                    eval(['this.basic_variable_directions.' variable_name ' = directions;']);
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
                            eval(['this.basic_variable_data.' this_variable_label ' = data;']);
                            eval(['this.time_data.' this_variable_label ' = time;']);
                            eval(['this.basic_variable_labels.' this_variable_label ' = labels;']);
                            eval(['this.basic_variable_directions.' this_variable_label ' = directions;']);
                        end
                        
                    else                    
                        % this is a true compound
                        if any(strcmp(this.basic_variable_labels.([this_variable_type '_trajectories']), this_variable_label))
                            % found this label in the data, so all is good
                            success = 1;
                        end
                    end
                end

                
                % calculate variables that can't be loaded
                
                % kinematics
                if strcmp(variable_name, 'com_rough_x')
                    component_x = 1;
                    
                    % grab required marker trajectories
                    LASI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LASI');
                    RASI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RASI');
                    LPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI');
                    RPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI');
                    LKNE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LKNE');
                    RKNE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RKNE');
                    LANK_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LANK');
                    RANK_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RANK');
                    LTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE');
                    RTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE');
                    C7_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    LFHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LFHD');
                    RFHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RFHD');
                    LBHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LBHD');
                    RBHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RBHD');
                    
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
                    this.basic_variable_data.(variable_name) = com_trajectory(:, component_x);
                    
                    % check directions
                    LPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI', 'indices');
                    RPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI', 'indices');
                    LASI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LASI', 'indices');
                    RASI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RASI', 'indices');
                    LKNE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LKNE', 'indices');
                    RKNE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RKNE', 'indices');
                    LANK_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LANK', 'indices');
                    RANK_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RANK', 'indices');
                    LTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE', 'indices');
                    RTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE', 'indices');
                    C7_indices   = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7', 'indices');
                    LFHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LFHD', 'indices');
                    RFHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RFHD', 'indices');
                    LBHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LBHD', 'indices');
                    RBHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RBHD', 'indices');
                    
                    LPSI_directions = this.basic_variable_directions.marker_trajectories(:, LPSI_indices(component_x));
                    RPSI_directions = this.basic_variable_directions.marker_trajectories(:, RPSI_indices(component_x));
                    LASI_directions = this.basic_variable_directions.marker_trajectories(:, LASI_indices(component_x));
                    RASI_directions = this.basic_variable_directions.marker_trajectories(:, RASI_indices(component_x));
                    LKNE_directions = this.basic_variable_directions.marker_trajectories(:, LKNE_indices(component_x));
                    RKNE_directions = this.basic_variable_directions.marker_trajectories(:, RKNE_indices(component_x));
                    LANK_directions = this.basic_variable_directions.marker_trajectories(:, LANK_indices(component_x));
                    RANK_directions = this.basic_variable_directions.marker_trajectories(:, RANK_indices(component_x));
                    LTOE_directions = this.basic_variable_directions.marker_trajectories(:, LTOE_indices(component_x));
                    RTOE_directions = this.basic_variable_directions.marker_trajectories(:, RTOE_indices(component_x));
                    LFHD_directions = this.basic_variable_directions.marker_trajectories(:, LFHD_indices(component_x));
                    RFHD_directions = this.basic_variable_directions.marker_trajectories(:, RFHD_indices(component_x));
                    LBHD_directions = this.basic_variable_directions.marker_trajectories(:, LBHD_indices(component_x));
                    RBHD_directions = this.basic_variable_directions.marker_trajectories(:, RBHD_indices(component_x));
                    C7_directions = this.basic_variable_directions.marker_trajectories(:, C7_indices(component_x));
                    
                    all_directions = ...
                      [
                        LPSI_directions, RPSI_directions, LASI_directions, RASI_directions, ...
                        LKNE_directions, RKNE_directions, LANK_directions, RANK_directions, LTOE_directions, RTOE_directions, ...
                        LFHD_directions, RFHD_directions, LBHD_directions, RBHD_directions, C7_directions ...
                      ]';
                    unique_directions = unique(cell2table(all_directions), 'rows');
                    if height(unique_directions) > 1
                        error('different directions found in marker data for rough CoM estimate')
                    end
                    this.basic_variable_directions.(variable_name)= LPSI_directions;
                    this.time_data.(variable_name) = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'com_rough_y')
                    component_y = 2;
                    
                    % grab required marker trajectories
                    LASI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LASI');
                    RASI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RASI');
                    LPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI');
                    RPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI');
                    LKNE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LKNE');
                    RKNE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RKNE');
                    LANK_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LANK');
                    RANK_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RANK');
                    LTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE');
                    RTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE');
                    C7_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    LFHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LFHD');
                    RFHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RFHD');
                    LBHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LBHD');
                    RBHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RBHD');
                    
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
                    this.basic_variable_data.(variable_name) = com_trajectory(:, component_y);
                    
                    % check directions
                    LPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI', 'indices');
                    RPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI', 'indices');
                    LASI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LASI', 'indices');
                    RASI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RASI', 'indices');
                    LKNE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LKNE', 'indices');
                    RKNE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RKNE', 'indices');
                    LANK_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LANK', 'indices');
                    RANK_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RANK', 'indices');
                    LTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE', 'indices');
                    RTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE', 'indices');
                    C7_indices   = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7', 'indices');
                    LFHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LFHD', 'indices');
                    RFHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RFHD', 'indices');
                    LBHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LBHD', 'indices');
                    RBHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RBHD', 'indices');
                    
                    LPSI_directions = this.basic_variable_directions.marker_trajectories(:, LPSI_indices(component_y));
                    RPSI_directions = this.basic_variable_directions.marker_trajectories(:, RPSI_indices(component_y));
                    LASI_directions = this.basic_variable_directions.marker_trajectories(:, LASI_indices(component_y));
                    RASI_directions = this.basic_variable_directions.marker_trajectories(:, RASI_indices(component_y));
                    LKNE_directions = this.basic_variable_directions.marker_trajectories(:, LKNE_indices(component_y));
                    RKNE_directions = this.basic_variable_directions.marker_trajectories(:, RKNE_indices(component_y));
                    LANK_directions = this.basic_variable_directions.marker_trajectories(:, LANK_indices(component_y));
                    RANK_directions = this.basic_variable_directions.marker_trajectories(:, RANK_indices(component_y));
                    LTOE_directions = this.basic_variable_directions.marker_trajectories(:, LTOE_indices(component_y));
                    RTOE_directions = this.basic_variable_directions.marker_trajectories(:, RTOE_indices(component_y));
                    LFHD_directions = this.basic_variable_directions.marker_trajectories(:, LFHD_indices(component_y));
                    RFHD_directions = this.basic_variable_directions.marker_trajectories(:, RFHD_indices(component_y));
                    LBHD_directions = this.basic_variable_directions.marker_trajectories(:, LBHD_indices(component_y));
                    RBHD_directions = this.basic_variable_directions.marker_trajectories(:, RBHD_indices(component_y));
                    C7_directions = this.basic_variable_directions.marker_trajectories(:, C7_indices(component_y));
                    
                    all_directions = ...
                      [
                        LPSI_directions, RPSI_directions, LASI_directions, RASI_directions, ...
                        LKNE_directions, RKNE_directions, LANK_directions, RANK_directions, LTOE_directions, RTOE_directions, ...
                        LFHD_directions, RFHD_directions, LBHD_directions, RBHD_directions, C7_directions ...
                      ]';
                    unique_directions = unique(cell2table(all_directions), 'rows');
                    if height(unique_directions) > 1
                        error('different directions found in marker data for rough CoM estimate')
                    end
                    this.basic_variable_directions.(variable_name)= LPSI_directions;
                    this.time_data.(variable_name) = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'com_rough_z')
                    component_z = 3;
                    
                    % grab required marker trajectories
                    LASI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LASI');
                    RASI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RASI');
                    LPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI');
                    RPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI');
                    LKNE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LKNE');
                    RKNE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RKNE');
                    LANK_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LANK');
                    RANK_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RANK');
                    LTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE');
                    RTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE');
                    C7_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    LFHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LFHD');
                    RFHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RFHD');
                    LBHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LBHD');
                    RBHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RBHD');
                    
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
                    this.basic_variable_data.(variable_name) = com_trajectory(:, component_z);
                    
                    % check directions
                    LPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI', 'indices');
                    RPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI', 'indices');
                    LASI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LASI', 'indices');
                    RASI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RASI', 'indices');
                    LKNE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LKNE', 'indices');
                    RKNE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RKNE', 'indices');
                    LANK_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LANK', 'indices');
                    RANK_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RANK', 'indices');
                    LTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE', 'indices');
                    RTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE', 'indices');
                    C7_indices   = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7', 'indices');
                    LFHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LFHD', 'indices');
                    RFHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RFHD', 'indices');
                    LBHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LBHD', 'indices');
                    RBHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RBHD', 'indices');
                    
                    LPSI_directions = this.basic_variable_directions.marker_trajectories(:, LPSI_indices(component_z));
                    RPSI_directions = this.basic_variable_directions.marker_trajectories(:, RPSI_indices(component_z));
                    LASI_directions = this.basic_variable_directions.marker_trajectories(:, LASI_indices(component_z));
                    RASI_directions = this.basic_variable_directions.marker_trajectories(:, RASI_indices(component_z));
                    LKNE_directions = this.basic_variable_directions.marker_trajectories(:, LKNE_indices(component_z));
                    RKNE_directions = this.basic_variable_directions.marker_trajectories(:, RKNE_indices(component_z));
                    LANK_directions = this.basic_variable_directions.marker_trajectories(:, LANK_indices(component_z));
                    RANK_directions = this.basic_variable_directions.marker_trajectories(:, RANK_indices(component_z));
                    LTOE_directions = this.basic_variable_directions.marker_trajectories(:, LTOE_indices(component_z));
                    RTOE_directions = this.basic_variable_directions.marker_trajectories(:, RTOE_indices(component_z));
                    LFHD_directions = this.basic_variable_directions.marker_trajectories(:, LFHD_indices(component_z));
                    RFHD_directions = this.basic_variable_directions.marker_trajectories(:, RFHD_indices(component_z));
                    LBHD_directions = this.basic_variable_directions.marker_trajectories(:, LBHD_indices(component_z));
                    RBHD_directions = this.basic_variable_directions.marker_trajectories(:, RBHD_indices(component_z));
                    C7_directions = this.basic_variable_directions.marker_trajectories(:, C7_indices(component_z));
                    
                    all_directions = ...
                      [
                        LPSI_directions, RPSI_directions, LASI_directions, RASI_directions, ...
                        LKNE_directions, RKNE_directions, LANK_directions, RANK_directions, LTOE_directions, RTOE_directions, ...
                        LFHD_directions, RFHD_directions, LBHD_directions, RBHD_directions, C7_directions ...
                      ]';
                    unique_directions = unique(cell2table(all_directions), 'rows');
                    if height(unique_directions) > 1
                        error('different directions found in marker data for rough CoM estimate')
                    end
                    this.basic_variable_directions.(variable_name)= LPSI_directions;
                    this.time_data.(variable_name) = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'com_rough_x_vel')
                    com_rough_x = this.getBasicVariableData('com_rough_x');
                    com_rough_x(com_rough_x==0) = NaN;
                    time = this.getTimeData('com_rough_x');
                    filter_order = this.study_settings.get('filter_order_com_vel');
                    cutoff_frequency = this.study_settings.get('filter_cutoff_com_vel');
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));
                    if any(~isnan(com_rough_x))
                        com_rough_x_vel = deriveByTime(nanfiltfilt(b, a, com_rough_x), 1/sampling_rate);
                    else
                        com_rough_x_vel = ones(size(com_rough_x)) * NaN;
                    end
                    this.basic_variable_data.com_rough_x_vel = com_rough_x_vel;
                    this.basic_variable_directions.com_rough_x_vel = this.basic_variable_directions.com_rough_x;
                    this.time_data.com_rough_x_vel = time;
                    success = 1;
                end                
                if strcmp(variable_name, 'com_rough_y_vel')
                    com_rough_y = this.getBasicVariableData('com_rough_y');
                    com_rough_y(com_rough_y==0) = NaN;
                    time = this.getTimeData('com_rough_y');
                    filter_order = this.study_settings.get('filter_order_com_vel');
                    cutoff_frequency = this.study_settings.get('filter_cutoff_com_vel');
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));
                    if any(~isnan(com_rough_y))
                        com_rough_y_vel = deriveByTime(nanfiltfilt(b, a, com_rough_y), 1/sampling_rate);
                    else
                        com_rough_y_vel = ones(size(com_rough_y)) * NaN;
                    end
                    this.basic_variable_data.com_rough_y_vel = com_rough_y_vel;
                    this.basic_variable_directions.com_rough_y_vel = this.basic_variable_directions.com_rough_y;
                    this.time_data.com_rough_y_vel = time;
                    success = 1;
                end                 
                if strcmp(variable_name, 'pelvis_y')
                    LPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI');
                    RPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI');
                    LASI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LASI');
                    RASI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RASI');
                    pelvis_trajectory = (LPSI_trajectory + RPSI_trajectory + LASI_trajectory + RASI_trajectory) * (1 / 4);
                    this.basic_variable_data.pelvis_y = pelvis_trajectory(:, 2);
                    
                    LPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI',  'indices');
                    RPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI',  'indices');
                    LASI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LASI',  'indices');
                    RASI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RASI',  'indices');
                    LPSI_directions = this.basic_variable_directions.marker_trajectories(:, LPSI_indices(2));
                    RPSI_directions = this.basic_variable_directions.marker_trajectories(:, RPSI_indices(2));
                    LASI_directions = this.basic_variable_directions.marker_trajectories(:, LASI_indices(2));
                    RASI_directions = this.basic_variable_directions.marker_trajectories(:, RASI_indices(2));
                    
                    if any(~strcmp(LPSI_directions{1}, {RPSI_directions{1}, LASI_directions{1}, RASI_directions{1}}))
                        error('LPSI, RPSI, LASI and RASI directions found in marker data are different from each other')
                    end
                    if any(~strcmp(LPSI_directions{2}, {RPSI_directions{2}, LASI_directions{2}, RASI_directions{2}}))
                        error('LPSI, RPSI, LASI and RASI directions found in marker data are different from each other')
                    end
                    this.basic_variable_directions.pelvis_y = LPSI_directions;
                    this.time_data.pelvis_y = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'c7_x')
                    C7_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    this.basic_variable_data.c7_x = C7_trajectory(:, 1);
                    C7_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7',  'indices');
                    this.basic_variable_directions.c7_x = this.basic_variable_directions.marker_trajectories(:, C7_indices(1));
                    this.time_data.c7_x = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'head_angle_ap')
                    % calculate angle trajectory
                    C7_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    LBHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LBHD');
                    RBHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RBHD');
                    MBHD_trajectory = (LBHD_trajectory + RBHD_trajectory) * 0.5;
                    head_vector_y = MBHD_trajectory(:, 2) - C7_trajectory(:, 2);
                    head_vector_z = MBHD_trajectory(:, 3) - C7_trajectory(:, 3);
                    this.basic_variable_data.head_angle_ap = atan2(head_vector_y, head_vector_z);
                    
                    % determine directions
                    C7_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7',  'indices');
                    LBHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LBHD',  'indices');
                    RBHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RBHD',  'indices');
                    C7_directions_y = this.basic_variable_directions.marker_trajectories(:, C7_indices(2));
                    LBHD_directions_y = this.basic_variable_directions.marker_trajectories(:, LBHD_indices(2));
                    RBHD_directions_y = this.basic_variable_directions.marker_trajectories(:, RBHD_indices(2));
                    C7_directions_z = this.basic_variable_directions.marker_trajectories(:, C7_indices(3));
                    LBHD_directions_z = this.basic_variable_directions.marker_trajectories(:, LBHD_indices(3));
                    RBHD_directions_z = this.basic_variable_directions.marker_trajectories(:, RBHD_indices(3));
                    
                    % check whether directions of markers are the same
                    if any(~strcmp(C7_directions_y{1}, {LBHD_directions_y{1}, RBHD_directions_y{1}}))
                        error('C7, LBHD and RBHD directions found in marker data are different from each other')
                    end
                    if any(~strcmp(C7_directions_y{2}, {LBHD_directions_y{2}, RBHD_directions_y{2}}))
                        error('C7, LBHD and RBHD directions found in marker data are different from each other')
                    end
                    if any(~strcmp(C7_directions_z{1}, {LBHD_directions_z{1}, RBHD_directions_z{1}}))
                        error('C7, LBHD and RBHD directions found in marker data are different from each other')
                    end
                    if any(~strcmp(C7_directions_z{2}, {LBHD_directions_z{2}, RBHD_directions_z{2}}))
                        error('C7, LBHD and RBHD directions found in marker data are different from each other')
                    end
                    
                    % check assumption that y is left-right and z is down-up
                    if ~strcmp(LBHD_directions_y{1}, 'forward')
                        error('Assuming positive y-direction for markers C7, LBHD and RBHD is "down"')
                    end                  
                    if ~strcmp(LBHD_directions_y{2}, 'backward')
                        error('Assuming negative y-direction for markers C7, LBHD and RBHD is "backward"')
                    end                  
                    if ~strcmp(LBHD_directions_z{1}, 'up')
                        error('Assuming positive z-direction for markers C7, LBHD and RBHD is "up"')
                    end                  
                    if ~strcmp(LBHD_directions_z{2}, 'down')
                        error('Assuming negative z-direction for markers C7, LBHD and RBHD is "down"')
                    end
                    
                    % all assumptions are met, define directions
                    this.basic_variable_directions.head_angle_ap = {'forward'; 'backward'};
                    
                    % time
                    this.time_data.head_angle_ap = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'head_angle_ml')
                    % calculate angle trajectory
                    C7_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    LBHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LBHD');
                    RBHD_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RBHD');
                    MBHD_trajectory = (LBHD_trajectory + RBHD_trajectory) * 0.5;
                    head_vector_x = MBHD_trajectory(:, 1) - C7_trajectory(:, 1);
                    head_vector_z = MBHD_trajectory(:, 3) - C7_trajectory(:, 3);
                    this.basic_variable_data.head_angle_ml = atan2(head_vector_x, head_vector_z);
                    
                    % determine directions
                    C7_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7',  'indices');
                    LBHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LBHD',  'indices');
                    RBHD_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RBHD',  'indices');
                    C7_directions_x = this.basic_variable_directions.marker_trajectories(:, C7_indices(1));
                    LBHD_directions_x = this.basic_variable_directions.marker_trajectories(:, LBHD_indices(1));
                    RBHD_directions_x = this.basic_variable_directions.marker_trajectories(:, RBHD_indices(1));
                    C7_directions_z = this.basic_variable_directions.marker_trajectories(:, C7_indices(3));
                    LBHD_directions_z = this.basic_variable_directions.marker_trajectories(:, LBHD_indices(3));
                    RBHD_directions_z = this.basic_variable_directions.marker_trajectories(:, RBHD_indices(3));
                    
                    % check whether directions of markers are the same
                    if any(~strcmp(C7_directions_x{1}, {LBHD_directions_x{1}, RBHD_directions_x{1}}))
                        error('C7, LBHD and RBHD directions found in marker data are different from each other')
                    end
                    if any(~strcmp(C7_directions_x{2}, {LBHD_directions_x{2}, RBHD_directions_x{2}}))
                        error('C7, LBHD and RBHD directions found in marker data are different from each other')
                    end
                    if any(~strcmp(C7_directions_z{1}, {LBHD_directions_z{1}, RBHD_directions_z{1}}))
                        error('C7, LBHD and RBHD directions found in marker data are different from each other')
                    end
                    if any(~strcmp(C7_directions_z{2}, {LBHD_directions_z{2}, RBHD_directions_z{2}}))
                        error('C7, LBHD and RBHD directions found in marker data are different from each other')
                    end
                    
                    % check assumption that y is left-right and z is down-up
                    if ~strcmp(LBHD_directions_x{1}, 'right')
                        error('Assuming positive x-direction for markers C7, LBHD and RBHD is "right"')
                    end                  
                    if ~strcmp(LBHD_directions_x{2}, 'left')
                        error('Assuming negative x-direction for markers C7, LBHD and RBHD is "left"')
                    end                  
                    if ~strcmp(LBHD_directions_z{1}, 'up')
                        error('Assuming positive z-direction for markers C7, LBHD and RBHD is "up"')
                    end                  
                    if ~strcmp(LBHD_directions_z{2}, 'down')
                        error('Assuming negative z-direction for markers C7, LBHD and RBHD is "down"')
                    end
                    
                    % all assumptions are met, define directions
                    this.basic_variable_directions.head_angle_ml = {'right'; 'left'};                    
                    
                    this.time_data.head_angle_ml = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'trunk_angle_ap')
                    % calculate angle trajectory
                    C7_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    LPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI');
                    RPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI');
                    MPSI_trajectory = (LPSI_trajectory + RPSI_trajectory) * 0.5;
                    trunk_vector_y = C7_trajectory(:, 2) - MPSI_trajectory(:, 2);
                    trunk_vector_z = C7_trajectory(:, 3) - MPSI_trajectory(:, 3);
                    this.basic_variable_data.trunk_angle_ap = atan2(trunk_vector_y, trunk_vector_z);
                    
                    % determine directions
                    C7_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7',  'indices');
                    LPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI',  'indices');
                    RPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI',  'indices');
                    C7_directions_y = this.basic_variable_directions.marker_trajectories(:, C7_indices(2));
                    LPSI_directions_y = this.basic_variable_directions.marker_trajectories(:, LPSI_indices(2));
                    RPSI_directions_y = this.basic_variable_directions.marker_trajectories(:, RPSI_indices(2));
                    C7_directions_z = this.basic_variable_directions.marker_trajectories(:, C7_indices(3));
                    LPSI_directions_z = this.basic_variable_directions.marker_trajectories(:, LPSI_indices(3));
                    RPSI_directions_z = this.basic_variable_directions.marker_trajectories(:, RPSI_indices(3));
                    
                    % check whether directions of markers are the same
                    if any(~strcmp(C7_directions_y{1}, {LPSI_directions_y{1}, RPSI_directions_y{1}}))
                        error('C7, LPSI and RPSI directions found in marker data are different from each other')
                    end
                    if any(~strcmp(C7_directions_y{2}, {LPSI_directions_y{2}, RPSI_directions_y{2}}))
                        error('C7, LPSI and RPSI directions found in marker data are different from each other')
                    end
                    if any(~strcmp(C7_directions_z{1}, {LPSI_directions_z{1}, RPSI_directions_z{1}}))
                        error('C7, LPSI and RPSI directions found in marker data are different from each other')
                    end
                    if any(~strcmp(C7_directions_z{2}, {LPSI_directions_z{2}, RPSI_directions_z{2}}))
                        error('C7, LPSI and RPSI directions found in marker data are different from each other')
                    end
                    
                    % check assumption that y is left-right and z is down-up
                    if ~strcmp(LPSI_directions_y{1}, 'forward')
                        error('Assuming positive y-direction for markers C7, LPSI and RPSI is "down"')
                    end                  
                    if ~strcmp(LPSI_directions_y{2}, 'backward')
                        error('Assuming negative y-direction for markers C7, LPSI and RPSI is "backward"')
                    end                  
                    if ~strcmp(LPSI_directions_z{1}, 'up')
                        error('Assuming positive z-direction for markers C7, LPSI and RPSI is "up"')
                    end                  
                    if ~strcmp(LPSI_directions_z{2}, 'down')
                        error('Assuming negative z-direction for markers C7, LPSI and RPSI is "down"')
                    end
                    
                    % all assumptions are met, define directions
                    this.basic_variable_directions.trunk_angle_ap = {'forward'; 'backward'};                    
                    
                    this.time_data.trunk_angle_ap = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'trunk_angle_ml')
                    % calculate angle trajectory
                    C7_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    LPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI');
                    RPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI');
                    MPSI_trajectory = (LPSI_trajectory + RPSI_trajectory) * 0.5;
                    trunk_vector_x = C7_trajectory(:, 1) - MPSI_trajectory(:, 1);
                    trunk_vector_z = C7_trajectory(:, 3) - MPSI_trajectory(:, 3);
                    this.basic_variable_data.trunk_angle_ml = atan2(trunk_vector_x, trunk_vector_z);
                    
                    % determine directions
                    C7_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7',  'indices');
                    LPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI',  'indices');
                    RPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI',  'indices');
                    C7_directions_x = this.basic_variable_directions.marker_trajectories(:, C7_indices(1));
                    LPSI_directions_x = this.basic_variable_directions.marker_trajectories(:, LPSI_indices(1));
                    RPSI_directions_x = this.basic_variable_directions.marker_trajectories(:, RPSI_indices(1));
                    C7_directions_z = this.basic_variable_directions.marker_trajectories(:, C7_indices(3));
                    LPSI_directions_z = this.basic_variable_directions.marker_trajectories(:, LPSI_indices(3));
                    RPSI_directions_z = this.basic_variable_directions.marker_trajectories(:, RPSI_indices(3));
                    
                    % check whether directions of markers are the same
                    if any(~strcmp(C7_directions_x{1}, {LPSI_directions_x{1}, RPSI_directions_x{1}}))
                        error('C7, LPSI and RPSI directions found in marker data are different from each other')
                    end
                    if any(~strcmp(C7_directions_x{2}, {LPSI_directions_x{2}, RPSI_directions_x{2}}))
                        error('C7, LPSI and RPSI directions found in marker data are different from each other')
                    end
                    if any(~strcmp(C7_directions_z{1}, {LPSI_directions_z{1}, RPSI_directions_z{1}}))
                        error('C7, LPSI and RPSI directions found in marker data are different from each other')
                    end
%                     if any(~strcmp(C7_directions_z{2}, {LBHD_directions_z{2}, RPSI_directions_z{2}}))
%                         error('C7, LPSI and RPSI directions found in marker data are different from each other')
%                     end
                    
                    % check assumption that y is left-right and z is down-up
                    if ~strcmp(LPSI_directions_x{1}, 'right')
                        error('Assuming positive x-direction for markers C7, LPSI and RPSI is "right"')
                    end                  
                    if ~strcmp(LPSI_directions_x{2}, 'left')
                        error('Assuming negative x-direction for markers C7, LPSI and RPSI is "left"')
                    end                  
                    if ~strcmp(LPSI_directions_z{1}, 'up')
                        error('Assuming positive z-direction for markers C7, LPSI and RPSI is "up"')
                    end                  
                    if ~strcmp(LPSI_directions_z{2}, 'down')
                        error('Assuming negative z-direction for markers C7, LPSI and RPSI is "down"')
                    end
                    
                    % all assumptions are met, define directions
                    this.basic_variable_directions.trunk_angle_ml = {'right'; 'left'};
                    
                    this.time_data.trunk_angle_ml = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'pelvis_angle_ml')
                    % calculate angle trajectory
                    left_hip_cor_trajectory = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'LHIPCOR');
                    right_hip_cor_trajectory = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'RHIPCOR');
                    pelvis_vector_x = right_hip_cor_trajectory(:, 1) - left_hip_cor_trajectory(:, 1);
                    pelvis_vector_z = right_hip_cor_trajectory(:, 3) - left_hip_cor_trajectory(:, 3);
                    this.basic_variable_data.pelvis_angle_ml = -atan2(pelvis_vector_z, pelvis_vector_x);
                    
                    % determine directions
                    LHIPCOR_indices = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'LHIPCOR',  'indices');
                    RHIPCOR_indices = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'RHIPCOR',  'indices');
                    LHIPCOR_directions_x = this.basic_variable_directions.joint_center_trajectories(:, LHIPCOR_indices(1));
                    RHIPCOR_directions_x = this.basic_variable_directions.joint_center_trajectories(:, RHIPCOR_indices(1));
                    LHIPCOR_directions_z = this.basic_variable_directions.joint_center_trajectories(:, LHIPCOR_indices(3));
                    RHIPCOR_directions_z = this.basic_variable_directions.joint_center_trajectories(:, RHIPCOR_indices(3));
                    
                    % check whether directions of markers are the same
                    if ~strcmp(LHIPCOR_directions_x{1}, RHIPCOR_directions_x{1})
                        error('LHIPCOR and RHIPCOR directions found in marker data are different from each other')
                    end
                    if ~strcmp(LHIPCOR_directions_x{2}, RHIPCOR_directions_x{2})
                        error('LHIPCOR and RHIPCOR directions found in marker data are different from each other')
                    end
                    if ~strcmp(LHIPCOR_directions_z{1}, RHIPCOR_directions_z{1})
                        error('LHIPCOR and RHIPCOR directions found in marker data are different from each other')
                    end
                    if ~strcmp(LHIPCOR_directions_z{2}, RHIPCOR_directions_z{2})
                        error('LHIPCOR and RHIPCOR directions found in marker data are different from each other')
                    end
                    
                    % check assumption that y is left-right and z is down-up
                    if ~strcmp(LHIPCOR_directions_x{1}, 'right')
                        error('Assuming positive x-direction for markers LHIPCOR and RHIPCOR is "right"')
                    end                  
                    if ~strcmp(LHIPCOR_directions_x{2}, 'left')
                        error('Assuming negative x-direction for markers LHIPCOR and RHIPCOR is "left"')
                    end                  
                    if ~strcmp(LHIPCOR_directions_z{1}, 'up')
                        error('Assuming positive z-direction for markers LHIPCOR and RHIPCOR is "up"')
                    end                  
                    if ~strcmp(LHIPCOR_directions_z{2}, 'down')
                        error('Assuming negative z-direction for markers LHIPCOR and RHIPCOR is "down"')
                    end
                    
                    % all assumptions are met, define directions
                    this.basic_variable_directions.pelvis_angle_ml = {'right'; 'left'};
                    
                    
                    this.time_data.pelvis_angle_ml = this.time_data.joint_center_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'left_leg_angle_ml')
                    % calculate angle trajectories
                    left_hip_cor_trajectory = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'LHIPCOR');
                    left_ankle_cor_trajectory = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'LANKLECOR');
                    left_leg_vector_x = left_hip_cor_trajectory(:, 1) - left_ankle_cor_trajectory(:, 1);
                    left_leg_vector_z = left_hip_cor_trajectory(:, 3) - left_ankle_cor_trajectory(:, 3);
                    this.basic_variable_data.left_leg_angle_ml = atan2(left_leg_vector_x, left_leg_vector_z);
                    
                    % determine directions
                    LHIPCOR_indices = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'LHIPCOR',  'indices');
                    LANKLECOR_indices = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'LANKLECOR',  'indices');
                    LHIPCOR_directions_x = this.basic_variable_directions.joint_center_trajectories(:, LHIPCOR_indices(1));
                    LANKLECOR_directions_x = this.basic_variable_directions.joint_center_trajectories(:, LANKLECOR_indices(1));
                    LHIPCOR_directions_z = this.basic_variable_directions.joint_center_trajectories(:, LHIPCOR_indices(3));
                    LANKLECOR_directions_z = this.basic_variable_directions.joint_center_trajectories(:, LANKLECOR_indices(3));
                    
                    % check whether directions of markers are the same
                    if ~strcmp(LHIPCOR_directions_x{1}, LANKLECOR_directions_x{1})
                        error('LHIPCOR and LANKLECOR directions found in marker data are different from each other')
                    end
                    if ~strcmp(LHIPCOR_directions_x{2}, LANKLECOR_directions_x{2})
                        error('LHIPCOR and LANKLECOR directions found in marker data are different from each other')
                    end
                    if ~strcmp(LHIPCOR_directions_z{1}, LANKLECOR_directions_z{1})
                        error('LHIPCOR and LANKLECOR directions found in marker data are different from each other')
                    end
                    if ~strcmp(LHIPCOR_directions_z{2}, LANKLECOR_directions_z{2})
                        error('LHIPCOR and LANKLECOR directions found in marker data are different from each other')
                    end
                    
                    % check assumption that y is left-right and z is down-up
                    if ~strcmp(LHIPCOR_directions_x{1}, 'right')
                        error('Assuming positive x-direction for markers LHIPCOR and LANKLECOR is "right"')
                    end                  
                    if ~strcmp(LHIPCOR_directions_x{2}, 'left')
                        error('Assuming negative x-direction for markers LHIPCOR and LANKLECOR is "left"')
                    end                  
                    if ~strcmp(LHIPCOR_directions_z{1}, 'up')
                        error('Assuming positive z-direction for markers LHIPCOR and LANKLECOR is "up"')
                    end                  
                    if ~strcmp(LHIPCOR_directions_z{2}, 'down')
                        error('Assuming negative z-direction for markers LHIPCOR and LANKLECOR is "down"')
                    end
                    
                    % all assumptions are met, define directions
                    this.basic_variable_directions.left_leg_angle_ml = {'right'; 'left'};
                    
                    this.time_data.left_leg_angle_ml = this.time_data.joint_center_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'right_leg_angle_ml')
                    % calculate angle trajectories
                    right_hip_cor_trajectory = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'RHIPCOR');
                    right_ankle_cor_trajectory = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'RANKLECOR');
                    right_leg_vector_x = right_hip_cor_trajectory(:, 1) - right_ankle_cor_trajectory(:, 1);
                    right_leg_vector_z = right_hip_cor_trajectory(:, 3) - right_ankle_cor_trajectory(:, 3);
                    this.basic_variable_data.right_leg_angle_ml = atan2(right_leg_vector_x, right_leg_vector_z);
                    
                    % determine directions
                    RHIPCOR_indices = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'RHIPCOR',  'indices');
                    RANKLECOR_indices = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'RANKLECOR',  'indices');
                    RHIPCOR_directions_x = this.basic_variable_directions.joint_center_trajectories(:, RHIPCOR_indices(1));
                    RANKLECOR_directions_x = this.basic_variable_directions.joint_center_trajectories(:, RANKLECOR_indices(1));
                    RHIPCOR_directions_z = this.basic_variable_directions.joint_center_trajectories(:, RHIPCOR_indices(3));
                    RANKLECOR_directions_z = this.basic_variable_directions.joint_center_trajectories(:, RANKLECOR_indices(3));
                    
                    % check whether directions of markers are the same
                    if ~strcmp(RHIPCOR_directions_x{1}, RANKLECOR_directions_x{1})
                        error('RHIPCOR and RANKLECOR directions found in marker data are different from each other')
                    end
                    if ~strcmp(RHIPCOR_directions_x{2}, RANKLECOR_directions_x{2})
                        error('RHIPCOR and RANKLECOR directions found in marker data are different from each other')
                    end
                    if ~strcmp(RHIPCOR_directions_z{1}, RANKLECOR_directions_z{1})
                        error('RHIPCOR and RANKLECOR directions found in marker data are different from each other')
                    end
                    if ~strcmp(RHIPCOR_directions_z{2}, RANKLECOR_directions_z{2})
                        error('RHIPCOR and RANKLECOR directions found in marker data are different from each other')
                    end
                    
                    % check assumption that y is left-right and z is down-up
                    if ~strcmp(RHIPCOR_directions_x{1}, 'right')
                        error('Assuming positive x-direction for markers RHIPCOR and RANKLECOR is "right"')
                    end                  
                    if ~strcmp(RHIPCOR_directions_x{2}, 'left')
                        error('Assuming negative x-direction for markers RHIPCOR and RANKLECOR is "left"')
                    end                  
                    if ~strcmp(RHIPCOR_directions_z{1}, 'up')
                        error('Assuming positive z-direction for markers RHIPCOR and RANKLECOR is "up"')
                    end                  
                    if ~strcmp(RHIPCOR_directions_z{2}, 'down')
                        error('Assuming negative z-direction for markers RHIPCOR and RANKLECOR is "down"')
                    end
                    
                    % all assumptions are met, define directions
                    this.basic_variable_directions.right_leg_angle_ml = {'right'; 'left'};
                    
                    this.time_data.right_leg_angle_ml = this.time_data.joint_center_trajectories;
                    success = 1;
                end               
                if strcmp(variable_name, 'left_leg_angle_ap')
                    % calculate angle trajectories
                    left_hip_cor_trajectory = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'LHIPCOR');
                    left_ankle_cor_trajectory = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'LANKLECOR');
                    left_leg_vector_y = left_hip_cor_trajectory(:, 2) - left_ankle_cor_trajectory(:, 2);
                    left_leg_vector_z = left_hip_cor_trajectory(:, 3) - left_ankle_cor_trajectory(:, 3);
                    this.basic_variable_data.left_leg_angle_ap = -atan2(left_leg_vector_y, left_leg_vector_z);
                    
                    % determine directions
                    LHIPCOR_indices = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'LHIPCOR',  'indices');
                    LANKLECOR_indices = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'LANKLECOR',  'indices');
                    LHIPCOR_directions_y = this.basic_variable_directions.joint_center_trajectories(:, LHIPCOR_indices(2));
                    LANKLECOR_directions_y = this.basic_variable_directions.joint_center_trajectories(:, LANKLECOR_indices(2));
                    LHIPCOR_directions_z = this.basic_variable_directions.joint_center_trajectories(:, LHIPCOR_indices(3));
                    LANKLECOR_directions_z = this.basic_variable_directions.joint_center_trajectories(:, LANKLECOR_indices(3));
                    
                    % check whether directions of markers are the same
                    if ~strcmp(LHIPCOR_directions_y{1}, LANKLECOR_directions_y{1})
                        error('LHIPCOR and LANKLECOR directions found in marker data are different from each other')
                    end
                    if ~strcmp(LHIPCOR_directions_y{2}, LANKLECOR_directions_y{2})
                        error('LHIPCOR and LANKLECOR directions found in marker data are different from each other')
                    end
                    if ~strcmp(LHIPCOR_directions_z{1}, LANKLECOR_directions_z{1})
                        error('LHIPCOR and LANKLECOR directions found in marker data are different from each other')
                    end
                    if ~strcmp(LHIPCOR_directions_z{2}, LANKLECOR_directions_z{2})
                        error('LHIPCOR and LANKLECOR directions found in marker data are different from each other')
                    end
                    
                    % check assumption that y is left-right and z is down-up
                    if ~strcmp(LHIPCOR_directions_y{1}, 'forward')
                        error('Assuming positive y-direction for markers LHIPCOR and LANKLECOR is "forward"')
                    end                  
                    if ~strcmp(LHIPCOR_directions_y{2}, 'backward')
                        error('Assuming negative y-direction for markers LHIPCOR and LANKLECOR is "backward"')
                    end                  
                    if ~strcmp(LHIPCOR_directions_z{1}, 'up')
                        error('Assuming positive z-direction for markers LHIPCOR and LANKLECOR is "up"')
                    end                  
                    if ~strcmp(LHIPCOR_directions_z{2}, 'down')
                        error('Assuming negative z-direction for markers LHIPCOR and LANKLECOR is "down"')
                    end
                    
                    % all assumptions are met, define directions
                    this.basic_variable_directions.left_leg_angle_ap = {'forward'; 'backward'};
                    
                    this.time_data.left_leg_angle_ap = this.time_data.joint_center_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'right_leg_angle_ap')
                    % calculate angle trajectories
                    right_hip_cor_trajectory = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'RHIPCOR');
                    right_ankle_cor_trajectory = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'RANKLECOR');
                    right_leg_vector_y = right_hip_cor_trajectory(:, 2) - right_ankle_cor_trajectory(:, 2);
                    right_leg_vector_z = right_hip_cor_trajectory(:, 3) - right_ankle_cor_trajectory(:, 3);
                    this.basic_variable_data.right_leg_angle_ap = -atan2(right_leg_vector_y, right_leg_vector_z);
                    
                    % determine directions
                    RHIPCOR_indices = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'RHIPCOR',  'indices');
                    RANKLECOR_indices = extractMarkerData(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'RANKLECOR',  'indices');
                    RHIPCOR_directions_y = this.basic_variable_directions.joint_center_trajectories(:, RHIPCOR_indices(2));
                    RANKLECOR_directions_y = this.basic_variable_directions.joint_center_trajectories(:, RANKLECOR_indices(2));
                    RHIPCOR_directions_z = this.basic_variable_directions.joint_center_trajectories(:, RHIPCOR_indices(3));
                    RANKLECOR_directions_z = this.basic_variable_directions.joint_center_trajectories(:, RANKLECOR_indices(3));
                    
                    % check whether directions of markers are the same
                    if ~strcmp(RHIPCOR_directions_y{1}, RANKLECOR_directions_y{1})
                        error('RHIPCOR and RANKLECOR directions found in marker data are different from each other')
                    end
                    if ~strcmp(RHIPCOR_directions_y{2}, RANKLECOR_directions_y{2})
                        error('RHIPCOR and RANKLECOR directions found in marker data are different from each other')
                    end
                    if ~strcmp(RHIPCOR_directions_z{1}, RANKLECOR_directions_z{1})
                        error('RHIPCOR and RANKLECOR directions found in marker data are different from each other')
                    end
                    if ~strcmp(RHIPCOR_directions_z{2}, RANKLECOR_directions_z{2})
                        error('RHIPCOR and RANKLECOR directions found in marker data are different from each other')
                    end
                    
                    % check assumption that y is left-right and z is down-up
                    if ~strcmp(RHIPCOR_directions_y{1}, 'forward')
                        error('Assuming positive x-direction for markers RHIPCOR and RANKLECOR is "forward"')
                    end                  
                    if ~strcmp(RHIPCOR_directions_y{2}, 'backward')
                        error('Assuming negative x-direction for markers RHIPCOR and RANKLECOR is "backward"')
                    end                  
                    if ~strcmp(RHIPCOR_directions_z{1}, 'up')
                        error('Assuming positive z-direction for markers RHIPCOR and RANKLECOR is "up"')
                    end                  
                    if ~strcmp(RHIPCOR_directions_z{2}, 'down')
                        error('Assuming negative z-direction for markers RHIPCOR and RANKLECOR is "down"')
                    end
                    
                    % all assumptions are met, define directions
                    this.basic_variable_directions.right_leg_angle_ap = {'forward'; 'backward'};
                    
                    this.time_data.right_leg_angle_ap = this.time_data.joint_center_trajectories;
                    success = 1;
                end               
                if strcmp(variable_name, 'left_foot_angle_ap')
                    % calculate angle trajectory
                    LTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE');
                    LTOEL_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOEL');
                    LHEE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE');
                    if isempty(LTOEL_trajectory)
                        LTOEM_trajectory = LTOE_trajectory;
                    else
                        LTOEM_trajectory = (LTOE_trajectory + LTOEL_trajectory) * 0.5;
                    end
                    foot_vector_x = LTOEM_trajectory(:, 1) - LHEE_trajectory(:, 1);
                    foot_vector_y = LTOEM_trajectory(:, 2) - LHEE_trajectory(:, 2);
                    foot_vector_z = LTOEM_trajectory(:, 3) - LHEE_trajectory(:, 3);
                    foot_vector_xy = (foot_vector_x.^2 + foot_vector_y.^2).^(0.5);
                    this.basic_variable_data.left_foot_angle_ap = atan2(foot_vector_z, foot_vector_xy);
                    
                    % determine directions
                    LTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE',  'indices');
                    LTOEL_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOEL',  'indices');
                    LHEE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE',  'indices');
                    if isempty(LTOEL_indices)
                        LTOEL_indices = LTOE_indices;
                    end
                    LTOE_directions_x = this.basic_variable_directions.marker_trajectories(:, LTOE_indices(1));
                    LTOEL_directions_x = this.basic_variable_directions.marker_trajectories(:, LTOEL_indices(1));
                    LHEE_directions_x = this.basic_variable_directions.marker_trajectories(:, LHEE_indices(1));
                    LTOE_directions_y = this.basic_variable_directions.marker_trajectories(:, LTOE_indices(2));
                    LTOEL_directions_y = this.basic_variable_directions.marker_trajectories(:, LTOEL_indices(2));
                    LHEE_directions_y = this.basic_variable_directions.marker_trajectories(:, LHEE_indices(2));
                    LTOE_directions_z = this.basic_variable_directions.marker_trajectories(:, LTOE_indices(3));
                    LTOEL_directions_z = this.basic_variable_directions.marker_trajectories(:, LTOEL_indices(3));
                    LHEE_directions_z = this.basic_variable_directions.marker_trajectories(:, LHEE_indices(3));
                    
                    % check whether directions of markers are the same
                    if any(~strcmp(LTOE_directions_x{1}, {LTOEL_directions_x{1}, LHEE_directions_x{1}}))
                        error('LTOE, LTOEL and LHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(LTOE_directions_x{2}, {LTOEL_directions_x{2}, LHEE_directions_x{2}}))
                        error('LTOE, LTOEL and LHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(LTOE_directions_y{1}, {LTOEL_directions_y{1}, LHEE_directions_y{1}}))
                        error('LTOE, LTOEL and LHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(LTOE_directions_y{2}, {LTOEL_directions_y{2}, LHEE_directions_y{2}}))
                        error('LTOE, LTOEL and LHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(LTOE_directions_z{1}, {LTOEL_directions_z{1}, LHEE_directions_z{1}}))
                        error('LTOE, LTOEL and LHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(LTOE_directions_z{2}, {LTOEL_directions_z{2}, LHEE_directions_z{2}}))
                        error('LTOE, LTOEL and LHEE directions found in marker data are different from each other')
                    end
                    
                    % check assumption that y is left-right and z is down-up
                    if ~strcmp(LTOEL_directions_x{1}, 'right')
                        error('Assuming positive y-direction for markers LTOE, LTOEL and LHEE is "right"')
                    end                  
                    if ~strcmp(LTOEL_directions_x{2}, 'left')
                        error('Assuming negative y-direction for markers LTOE, LTOEL and LHEE is "left"')
                    end                  
                    if ~strcmp(LTOEL_directions_y{1}, 'forward')
                        error('Assuming positive y-direction for markers LTOE, LTOEL and LHEE is "down"')
                    end                  
                    if ~strcmp(LTOEL_directions_y{2}, 'backward')
                        error('Assuming negative y-direction for markers LTOE, LTOEL and LHEE is "backward"')
                    end                  
                    if ~strcmp(LTOEL_directions_z{1}, 'up')
                        error('Assuming positive z-direction for markers LTOE, LTOEL and LHEE is "up"')
                    end                  
                    if ~strcmp(LTOEL_directions_z{2}, 'down')
                        error('Assuming negative z-direction for markers LTOE, LTOEL and LHEE is "down"')
                    end
                    
                    % all assumptions are met, define directions
                    this.basic_variable_directions.left_foot_angle_ap = {'up'; 'down'};
                                        
                    this.time_data.left_foot_angle_ap = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'left_foot_angle_ml')
                    % TODO: test this
                    % calculate angle trajectory
                    LTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE');
                    LTOEL_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOEL');
                    foot_vector_x = LTOEL_trajectory(:, 1) - LTOE_trajectory(:, 1);
                    foot_vector_y = LTOEL_trajectory(:, 2) - LTOE_trajectory(:, 2);
                    foot_vector_z = LTOEL_trajectory(:, 3) - LTOE_trajectory(:, 3);
                    foot_vector_xy = (foot_vector_x.^2 + foot_vector_y.^2).^(0.5);
                    this.basic_variable_data.left_foot_angle_ml = atan2(foot_vector_z, foot_vector_xy);
                    
                    % determine directions
                    LTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE',  'indices');
                    LTOEL_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOEL',  'indices');
                    LTOE_directions_x = this.basic_variable_directions.marker_trajectories(:, LTOE_indices(1));
                    LTOEL_directions_x = this.basic_variable_directions.marker_trajectories(:, LTOEL_indices(1));
                    LTOE_directions_y = this.basic_variable_directions.marker_trajectories(:, LTOE_indices(2));
                    LTOEL_directions_y = this.basic_variable_directions.marker_trajectories(:, LTOEL_indices(2));
                    LTOE_directions_z = this.basic_variable_directions.marker_trajectories(:, LTOE_indices(3));
                    LTOEL_directions_z = this.basic_variable_directions.marker_trajectories(:, LTOEL_indices(3));
                    
                    % check whether directions of markers are the same
                    if ~strcmp(LTOE_directions_x{1}, LTOEL_directions_x{1})
                        error('LTOE and LTOEL directions found in marker data are different from each other')
                    end
                    if ~strcmp(LTOE_directions_x{2}, LTOEL_directions_x{2})
                        error('LTOE and LTOEL directions found in marker data are different from each other')
                    end
                    if ~strcmp(LTOE_directions_y{1}, LTOEL_directions_y{1})
                        error('LTOE and LTOEL directions found in marker data are different from each other')
                    end
                    if ~strcmp(LTOE_directions_y{2}, LTOEL_directions_y{2})
                        error('LTOE and LTOEL directions found in marker data are different from each other')
                    end
                    if ~strcmp(LTOE_directions_z{1}, LTOEL_directions_z{1})
                        error('LTOE and LTOEL directions found in marker data are different from each other')
                    end
                    if ~strcmp(LTOE_directions_z{2}, LTOEL_directions_z{2})
                        error('LTOE and LTOEL directions found in marker data are different from each other')
                    end
                    
                    % check assumption that y is left-right and z is down-up
                    if ~strcmp(LTOEL_directions_x{1}, 'right')
                        error('Assuming positive y-direction for markers LTOE and LTOEL is "right"')
                    end                  
                    if ~strcmp(LTOEL_directions_x{2}, 'left')
                        error('Assuming negative y-direction for markers LTOE and LTOEL is "left"')
                    end                  
                    if ~strcmp(LTOEL_directions_y{1}, 'forward')
                        error('Assuming positive y-direction for markers LTOE and LTOEL is "down"')
                    end                  
                    if ~strcmp(LTOEL_directions_y{2}, 'backward')
                        error('Assuming negative y-direction for markers LTOE and LTOEL is "backward"')
                    end                  
                    if ~strcmp(LTOEL_directions_z{1}, 'up')
                        error('Assuming positive z-direction for markers LTOE and LTOEL is "up"')
                    end                  
                    if ~strcmp(LTOEL_directions_z{2}, 'down')
                        error('Assuming negative z-direction for markers LTOE and LTOEL is "down"')
                    end
                    
                    % all assumptions are met, define directions
                    this.basic_variable_directions.left_foot_angle_ml = {'right'; 'left'};
                    
                    this.time_data.left_foot_angle_ml = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'left_foot_angle_yaw')
                    % calculate angle trajectory
                    LTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE');
                    LTOEL_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOEL');
                    LHEE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE');
                    if isempty(LTOEL_trajectory)
                        LTOEM_trajectory = LTOE_trajectory;
                    else
                        LTOEM_trajectory = (LTOE_trajectory + LTOEL_trajectory) * 0.5;
                    end
                    foot_vector_x = LTOEM_trajectory(:, 1) - LHEE_trajectory(:, 1);
                    foot_vector_y = LTOEM_trajectory(:, 2) - LHEE_trajectory(:, 2);
                    this.basic_variable_data.left_foot_angle_yaw = atan2(foot_vector_x, foot_vector_y);
                    
                    % determine directions
                    LTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE',  'indices');
                    LTOEL_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOEL',  'indices');
                    LHEE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE',  'indices');
                    if isempty(LTOEL_indices)
                        LTOEL_indices = LTOE_indices;
                    end
                    LTOE_directions_x = this.basic_variable_directions.marker_trajectories(:, LTOE_indices(1));
                    LTOEL_directions_x = this.basic_variable_directions.marker_trajectories(:, LTOEL_indices(1));
                    LHEE_directions_x = this.basic_variable_directions.marker_trajectories(:, LHEE_indices(1));
                    LTOE_directions_y = this.basic_variable_directions.marker_trajectories(:, LTOE_indices(2));
                    LTOEL_directions_y = this.basic_variable_directions.marker_trajectories(:, LTOEL_indices(2));
                    LHEE_directions_y = this.basic_variable_directions.marker_trajectories(:, LHEE_indices(2));
                    LTOE_directions_z = this.basic_variable_directions.marker_trajectories(:, LTOE_indices(3));
                    LTOEL_directions_z = this.basic_variable_directions.marker_trajectories(:, LTOEL_indices(3));
                    LHEE_directions_z = this.basic_variable_directions.marker_trajectories(:, LHEE_indices(3));
                    
                    % check whether directions of markers are the same
                    if any(~strcmp(LTOE_directions_x{1}, {LTOEL_directions_x{1}, LHEE_directions_x{1}}))
                        error('LTOE, LTOEL and LHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(LTOE_directions_x{2}, {LTOEL_directions_x{2}, LHEE_directions_x{2}}))
                        error('LTOE, LTOEL and LHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(LTOE_directions_y{1}, {LTOEL_directions_y{1}, LHEE_directions_y{1}}))
                        error('LTOE, LTOEL and LHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(LTOE_directions_y{2}, {LTOEL_directions_y{2}, LHEE_directions_y{2}}))
                        error('LTOE, LTOEL and LHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(LTOE_directions_z{1}, {LTOEL_directions_z{1}, LHEE_directions_z{1}}))
                        error('LTOE, LTOEL and LHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(LTOE_directions_z{2}, {LTOEL_directions_z{2}, LHEE_directions_z{2}}))
                        error('LTOE, LTOEL and LHEE directions found in marker data are different from each other')
                    end
                    
                    % check assumption that y is left-right and z is down-up
                    if ~strcmp(LTOEL_directions_x{1}, 'right')
                        error('Assuming positive y-direction for markers LTOE, LTOEL and LHEE is "right"')
                    end                  
                    if ~strcmp(LTOEL_directions_x{2}, 'left')
                        error('Assuming negative y-direction for markers LTOE, LTOEL and LHEE is "left"')
                    end                  
                    if ~strcmp(LTOEL_directions_y{1}, 'forward')
                        error('Assuming positive y-direction for markers LTOE, LTOEL and LHEE is "down"')
                    end                  
                    if ~strcmp(LTOEL_directions_y{2}, 'backward')
                        error('Assuming negative y-direction for markers LTOE, LTOEL and LHEE is "backward"')
                    end                  
                    if ~strcmp(LTOEL_directions_z{1}, 'up')
                        error('Assuming positive z-direction for markers LTOE, LTOEL and LHEE is "up"')
                    end                  
                    if ~strcmp(LTOEL_directions_z{2}, 'down')
                        error('Assuming negative z-direction for markers LTOE, LTOEL and LHEE is "down"')
                    end
                    
                    % all assumptions are met, define directions
                    this.basic_variable_directions.left_foot_angle_yaw = {'clockwise'; 'counterclockwise'};
                                        
                    this.time_data.left_foot_angle_yaw = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'right_foot_angle_ap')
                    % calculate angle trajectory
                    RTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE');
                    RTOEL_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOEL');
                    RHEE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE');
                    if isempty(RTOEL_trajectory)
                        RTOEM_trajectory = RTOE_trajectory;
                    else
                        RTOEM_trajectory = (RTOE_trajectory + RTOEL_trajectory) * 0.5;
                    end
                    foot_vector_x = RTOEM_trajectory(:, 1) - RHEE_trajectory(:, 1);
                    foot_vector_y = RTOEM_trajectory(:, 2) - RHEE_trajectory(:, 2);
                    foot_vector_z = RTOEM_trajectory(:, 3) - RHEE_trajectory(:, 3);
                    foot_vector_xy = (foot_vector_x.^2 + foot_vector_y.^2).^(0.5);
                    this.basic_variable_data.right_foot_angle_ap = atan2(foot_vector_z, foot_vector_xy);
                    
                    % determine directions
                    RTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE',  'indices');
                    RTOEL_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOEL',  'indices');
                    RHEE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE',  'indices');
                    if isempty(RTOEL_indices)
                        RTOEL_indices = RTOE_indices;
                    end
                    RTOE_directions_x = this.basic_variable_directions.marker_trajectories(:, RTOE_indices(1));
                    RTOEL_directions_x = this.basic_variable_directions.marker_trajectories(:, RTOEL_indices(1));
                    RHEE_directions_x = this.basic_variable_directions.marker_trajectories(:, RHEE_indices(1));
                    RTOE_directions_y = this.basic_variable_directions.marker_trajectories(:, RTOE_indices(2));
                    RTOEL_directions_y = this.basic_variable_directions.marker_trajectories(:, RTOEL_indices(2));
                    RHEE_directions_y = this.basic_variable_directions.marker_trajectories(:, RHEE_indices(2));
                    RTOE_directions_z = this.basic_variable_directions.marker_trajectories(:, RTOE_indices(3));
                    RTOEL_directions_z = this.basic_variable_directions.marker_trajectories(:, RTOEL_indices(3));
                    RHEE_directions_z = this.basic_variable_directions.marker_trajectories(:, RHEE_indices(3));
                    
                    % check whether directions of markers are the same
                    if any(~strcmp(RTOE_directions_x{1}, {RTOEL_directions_x{1}, RHEE_directions_x{1}}))
                        error('RTOE, RTOEL and RHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(RTOE_directions_x{2}, {RTOEL_directions_x{2}, RHEE_directions_x{2}}))
                        error('RTOE, RTOEL and RHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(RTOE_directions_y{1}, {RTOEL_directions_y{1}, RHEE_directions_y{1}}))
                        error('RTOE, RTOEL and RHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(RTOE_directions_y{2}, {RTOEL_directions_y{2}, RHEE_directions_y{2}}))
                        error('RTOE, RTOEL and RHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(RTOE_directions_z{1}, {RTOEL_directions_z{1}, RHEE_directions_z{1}}))
                        error('RTOE, RTOEL and RHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(RTOE_directions_z{2}, {RTOEL_directions_z{2}, RHEE_directions_z{2}}))
                        error('RTOE, RTOEL and RHEE directions found in marker data are different from each other')
                    end
                    
                    % check assumption that y is left-right and z is down-up
                    if ~strcmp(RTOEL_directions_x{1}, 'right')
                        error('Assuming positive y-direction for markers RTOE, RTOEL and RHEE is "right"')
                    end                  
                    if ~strcmp(RTOEL_directions_x{2}, 'left')
                        error('Assuming negative y-direction for markers RTOE, RTOEL and RHEE is "left"')
                    end                  
                    if ~strcmp(RTOEL_directions_y{1}, 'forward')
                        error('Assuming positive y-direction for markers RTOE, RTOEL and RHEE is "down"')
                    end                  
                    if ~strcmp(RTOEL_directions_y{2}, 'backward')
                        error('Assuming negative y-direction for markers RTOE, RTOEL and RHEE is "backward"')
                    end                  
                    if ~strcmp(RTOEL_directions_z{1}, 'up')
                        error('Assuming positive z-direction for markers RTOE, RTOEL and RHEE is "up"')
                    end                  
                    if ~strcmp(RTOEL_directions_z{2}, 'down')
                        error('Assuming negative z-direction for markers RTOE, RTOEL and RHEE is "down"')
                    end
                    
                    % all assumptions are met, define directions
                    this.basic_variable_directions.right_foot_angle_ap = {'up'; 'down'};                    
                    
                    this.time_data.right_foot_angle_ap = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'right_foot_angle_ml')
                    % TODO: test this
                    % calculate angle trajectory
                    RTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE');
                    RTOEL_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOEL');
                    foot_vector_x = RTOEL_trajectory(:, 1) - RTOE_trajectory(:, 1);
                    foot_vector_y = RTOEL_trajectory(:, 2) - RTOE_trajectory(:, 2);
                    foot_vector_z = RTOEL_trajectory(:, 3) - RTOE_trajectory(:, 3);
                    foot_vector_xy = (foot_vector_x.^2 + foot_vector_y.^2).^(0.5);
                    this.basic_variable_data.right_foot_angle_ml = -atan2(foot_vector_z, foot_vector_xy);
                    
                    % determine directions
                    RTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE',  'indices');
                    RTOEL_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOEL',  'indices');
                    RTOE_directions_x = this.basic_variable_directions.marker_trajectories(:, RTOE_indices(1));
                    RTOEL_directions_x = this.basic_variable_directions.marker_trajectories(:, RTOEL_indices(1));
                    RTOE_directions_y = this.basic_variable_directions.marker_trajectories(:, RTOE_indices(2));
                    RTOEL_directions_y = this.basic_variable_directions.marker_trajectories(:, RTOEL_indices(2));
                    RTOE_directions_z = this.basic_variable_directions.marker_trajectories(:, RTOE_indices(3));
                    RTOEL_directions_z = this.basic_variable_directions.marker_trajectories(:, RTOEL_indices(3));
                    
                    % check whether directions of markers are the same
                    if ~strcmp(RTOE_directions_x{1}, RTOEL_directions_x{1})
                        error('RTOE and RTOEL directions found in marker data are different from each other')
                    end
                    if ~strcmp(RTOE_directions_x{2}, RTOEL_directions_x{2})
                        error('RTOE and RTOEL directions found in marker data are different from each other')
                    end
                    if ~strcmp(RTOE_directions_y{1}, RTOEL_directions_y{1})
                        error('RTOE and RTOEL directions found in marker data are different from each other')
                    end
                    if ~strcmp(RTOE_directions_y{2}, RTOEL_directions_y{2})
                        error('RTOE and RTOEL directions found in marker data are different from each other')
                    end
                    if ~strcmp(RTOE_directions_z{1}, RTOEL_directions_z{1})
                        error('RTOE and RTOEL directions found in marker data are different from each other')
                    end
                    if ~strcmp(RTOE_directions_z{2}, RTOEL_directions_z{2})
                        error('RTOE and RTOEL directions found in marker data are different from each other')
                    end
                    
                    % check assumption that y is left-right and z is down-up
                    if ~strcmp(RTOEL_directions_x{1}, 'right')
                        error('Assuming positive y-direction for markers RTOE and RTOEL is "right"')
                    end                  
                    if ~strcmp(RTOEL_directions_x{2}, 'left')
                        error('Assuming negative y-direction for markers RTOE and RTOEL is "left"')
                    end                  
                    if ~strcmp(RTOEL_directions_y{1}, 'forward')
                        error('Assuming positive y-direction for markers RTOE and RTOEL is "down"')
                    end                  
                    if ~strcmp(RTOEL_directions_y{2}, 'backward')
                        error('Assuming negative y-direction for markers RTOE and RTOEL is "backward"')
                    end                  
                    if ~strcmp(RTOEL_directions_z{1}, 'up')
                        error('Assuming positive z-direction for markers RTOE and RTOEL is "up"')
                    end                  
                    if ~strcmp(RTOEL_directions_z{2}, 'down')
                        error('Assuming negative z-direction for markers RTOE and RTOEL is "down"')
                    end
                    
                    % all assumptions are met, define directions
                    this.basic_variable_directions.right_foot_angle_ml = {'right'; 'left'};
                    
                    this.time_data.right_foot_angle_ml = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'right_foot_angle_yaw')
                    % calculate angle trajectory
                    RTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE');
                    RTOEL_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOEL');
                    RHEE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE');
                    if isempty(RTOEL_trajectory)
                        RTOEM_trajectory = RTOE_trajectory;
                    else
                        RTOEM_trajectory = (RTOE_trajectory + RTOEL_trajectory) * 0.5;
                    end
                    foot_vector_x = RTOEM_trajectory(:, 1) - RHEE_trajectory(:, 1);
                    foot_vector_y = RTOEM_trajectory(:, 2) - RHEE_trajectory(:, 2);
                    this.basic_variable_data.right_foot_angle_yaw = atan2(foot_vector_x, foot_vector_y);
                    
                    % determine directions
                    RTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE',  'indices');
                    RTOEL_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOEL',  'indices');
                    RHEE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE',  'indices');
                    if isempty(RTOEL_indices)
                        RTOEL_indices = RTOE_indices;
                    end
                    RTOE_directions_x = this.basic_variable_directions.marker_trajectories(:, RTOE_indices(1));
                    RTOEL_directions_x = this.basic_variable_directions.marker_trajectories(:, RTOEL_indices(1));
                    RHEE_directions_x = this.basic_variable_directions.marker_trajectories(:, RHEE_indices(1));
                    RTOE_directions_y = this.basic_variable_directions.marker_trajectories(:, RTOE_indices(2));
                    RTOEL_directions_y = this.basic_variable_directions.marker_trajectories(:, RTOEL_indices(2));
                    RHEE_directions_y = this.basic_variable_directions.marker_trajectories(:, RHEE_indices(2));
                    RTOE_directions_z = this.basic_variable_directions.marker_trajectories(:, RTOE_indices(3));
                    RTOEL_directions_z = this.basic_variable_directions.marker_trajectories(:, RTOEL_indices(3));
                    RHEE_directions_z = this.basic_variable_directions.marker_trajectories(:, RHEE_indices(3));
                    
                    % check whether directions of markers are the same
                    if any(~strcmp(RTOE_directions_x{1}, {RTOEL_directions_x{1}, RHEE_directions_x{1}}))
                        error('RTOE, RTOEL and RHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(RTOE_directions_x{2}, {RTOEL_directions_x{2}, RHEE_directions_x{2}}))
                        error('RTOE, RTOEL and RHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(RTOE_directions_y{1}, {RTOEL_directions_y{1}, RHEE_directions_y{1}}))
                        error('RTOE, RTOEL and RHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(RTOE_directions_y{2}, {RTOEL_directions_y{2}, RHEE_directions_y{2}}))
                        error('RTOE, RTOEL and RHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(RTOE_directions_z{1}, {RTOEL_directions_z{1}, RHEE_directions_z{1}}))
                        error('RTOE, RTOEL and RHEE directions found in marker data are different from each other')
                    end
                    if any(~strcmp(RTOE_directions_z{2}, {RTOEL_directions_z{2}, RHEE_directions_z{2}}))
                        error('RTOE, RTOEL and RHEE directions found in marker data are different from each other')
                    end
                    
                    % check assumption that y is left-right and z is down-up
                    if ~strcmp(RTOEL_directions_x{1}, 'right')
                        error('Assuming positive y-direction for markers RTOE, RTOEL and RHEE is "right"')
                    end                  
                    if ~strcmp(RTOEL_directions_x{2}, 'left')
                        error('Assuming negative y-direction for markers RTOE, RTOEL and RHEE is "left"')
                    end                  
                    if ~strcmp(RTOEL_directions_y{1}, 'forward')
                        error('Assuming positive y-direction for markers RTOE, RTOEL and RHEE is "down"')
                    end                  
                    if ~strcmp(RTOEL_directions_y{2}, 'backward')
                        error('Assuming negative y-direction for markers RTOE, RTOEL and RHEE is "backward"')
                    end                  
                    if ~strcmp(RTOEL_directions_z{1}, 'up')
                        error('Assuming positive z-direction for markers RTOE, RTOEL and RHEE is "up"')
                    end                  
                    if ~strcmp(RTOEL_directions_z{2}, 'down')
                        error('Assuming negative z-direction for markers RTOE, RTOEL and RHEE is "down"')
                    end
                    
                    % all assumptions are met, define directions
                    this.basic_variable_directions.right_foot_angle_yaw = {'clockwise'; 'counterclockwise'};
                                        
                    this.time_data.right_foot_angle_yaw = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'left_arm_right_leg_relative_phase')
                    left_arm_phase = this.getBasicVariableData('left_arm_phase');
                    right_leg_phase = this.getBasicVariableData('right_leg_phase');
                    
                    left_arm_right_leg_relative_phase = normalizeAngle(left_arm_phase - right_leg_phase);
                    
                    % store
                    this.basic_variable_data.left_arm_right_leg_relative_phase = left_arm_right_leg_relative_phase;
                    this.basic_variable_directions.left_arm_right_leg_relative_phase = {'lead'; 'lag'};
                    this.time_data.left_arm_right_leg_relative_phase = this.time_data.marker_trajectories;
                    success = 1;
                end             
                if strcmp(variable_name, 'right_arm_left_leg_relative_phase')
                    right_arm_phase = this.getBasicVariableData('right_arm_phase');
                    left_leg_phase = this.getBasicVariableData('left_leg_phase');
                    
                    right_arm_left_leg_relative_phase = normalizeAngle(right_arm_phase - left_leg_phase);
                    
                    % store
                    this.basic_variable_data.right_arm_left_leg_relative_phase = right_arm_left_leg_relative_phase;
                    this.basic_variable_directions.right_arm_left_leg_relative_phase = {'lead'; 'lag'};
                    this.time_data.right_arm_left_leg_relative_phase = this.time_data.marker_trajectories;
                    success = 1;
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
                if ~success & ~any(strcmp(variable_name, this.basic_variable_load_failures))
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
                    stretch_directions = [];
                    
                    % extract time
                    this_stretch_times = stretch_times(i_stretch, :);
                    this_stretch_start_time = this_stretch_times(1);
                    this_stretch_end_time = this_stretch_times(end);
                    
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
                        for i_band = number_of_bands : -1 : 1
                            % going backward makes a difference for the
                            % junction points between two steps. We want to
                            % use data from the earlier step for BoS, so we
                            % go backward
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
                        
                        % now calculate XCoM - BoS
                        xcom_x = stretch_variables{strcmp(this.stretch_variable_names, 'xcom_x')}(:, i_stretch);
                        stretch_data = xcom_x - bos_x_data;
                    end
                    if strcmp(variable_name, 'mos_y')
                        % first calculate base of support
                        LTOEL_y = this.getTimeNormalizedData('marker:LTOEL_y', this_stretch_times);
                        RTOEL_y = this.getTimeNormalizedData('marker:RTOEL_y', this_stretch_times);
                        bos_y_data = zeros(size(LTOEL_y));
                        for i_band = number_of_bands : -1 : 1
                            % going backward makes a difference for the
                            % junction points between two steps. We want to
                            % use data from the earlier step for BoS, so we
                            % go backward
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
                        for i_band = number_of_bands : -1 : 1
                            % going backward makes a difference for the
                            % junction points between two steps. We want to
                            % use data from the earlier step for BoS, so we
                            % go backward
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
                        
                        % now calculate XCoM - BoS
                        xcom_mpsis_x = stretch_variables{strcmp(this.stretch_variable_names, 'xcom_mpsis_x')}(:, i_stretch);
                        stretch_data = xcom_mpsis_x - bos_x_data;
                    end
                    if strcmp(variable_name, 'mos_mpsis_y')
                        % first calculate base of support
                        LTOEL_y = this.getTimeNormalizedData('marker:LTOEL_y', this_stretch_times);
                        RTOEL_y = this.getTimeNormalizedData('marker:RTOEL_y', this_stretch_times);
                        bos_y_data = zeros(size(LTOEL_y));
                        for i_band = number_of_bands : -1 : 1
                            % going backward makes a difference for the
                            % junction points between two steps. We want to
                            % use data from the earlier step for BoS, so we
                            % go backward
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
                        
                        % now calculate XCoM - BoS
                        xcom_mpsis_y = stretch_variables{strcmp(this.stretch_variable_names, 'xcom_mpsis_y')}(:, i_stretch);
                        stretch_data = xcom_mpsis_y - bos_y_data;
                    end        
                    if strcmp(variable_name, 'xcom_rough_x')
                        % get mpsis position
                        com_rough_x =  this.getTimeNormalizedData('com_rough_x', this_stretch_times);
                        com_rough_y =  this.getTimeNormalizedData('com_rough_y', this_stretch_times);
                        com_rough_z =  this.getTimeNormalizedData('com_rough_z', this_stretch_times);
                        
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
                        com_rough_x_vel =  this.getTimeNormalizedData('com_rough_x_vel', this_stretch_times);
                        omega_0 = (9.81 * leg_length_data.^(-1)).^(0.5);
%                         omega_0 = sqrt(9.81/leg_length);
                       
                        stretch_data = com_rough_x + omega_0.^(-1) .* com_rough_x_vel;
                    end
                    if strcmp(variable_name, 'xcom_rough_y')
                        % get mpsis position
                        com_rough_x =  this.getTimeNormalizedData('com_rough_x', this_stretch_times);
                        com_rough_y =  this.getTimeNormalizedData('com_rough_y', this_stretch_times);
                        com_rough_z =  this.getTimeNormalizedData('com_rough_z', this_stretch_times);
                        
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
                        com_rough_y_vel =  this.getTimeNormalizedData('com_rough_y_vel', this_stretch_times);
                        omega_0 = (9.81 * leg_length_data.^(-1)).^(0.5);
%                         omega_0 = sqrt(9.81/leg_length);
                       
                        stretch_data = com_rough_y + omega_0.^(-1) .* com_rough_y_vel;
                    end
                    if strcmp(variable_name, 'mos_rough_x')
                        % first calculate base of support
                        LTOEL_x = this.getTimeNormalizedData('marker:LTOEL_x', this_stretch_times);
                        RTOEL_x = this.getTimeNormalizedData('marker:RTOEL_x', this_stretch_times);
                        bos_x_data = zeros(size(LTOEL_x));
                        for i_band = number_of_bands : -1 : 1
                            % going backward makes a difference for the
                            % junction points between two steps. We want to
                            % use data from the earlier step for BoS, so we
                            % go backward
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
                        
                        % now calculate XCoM - BoS
                        xcom_rough_x = stretch_variables{strcmp(this.stretch_variable_names, 'xcom_rough_x')}(:, i_stretch);
                        stretch_data = xcom_rough_x - bos_x_data;
                    end  
                    if strcmp(variable_name, 'mos_rough_y')
                        % first calculate base of support
                        LTOEL_y = this.getTimeNormalizedData('marker:LTOEL_y', this_stretch_times);
                        RTOEL_y = this.getTimeNormalizedData('marker:RTOEL_y', this_stretch_times);
                        bos_y_data = zeros(size(LTOEL_y));
                        for i_band = number_of_bands : -1 : 1
                            % going backward makes a difference for the
                            % junction points between two steps. We want to
                            % use data from the earlier step for BoS, so we
                            % go backward
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
                        
                        % now calculate XCoM - BoS
                        xcom_rough_y = stretch_variables{strcmp(this.stretch_variable_names, 'xcom_rough_y')}(:, i_stretch);
                        stretch_data = xcom_rough_y - bos_y_data;
                    end
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
                    if strcmp(variable_name, 'step_time')
                        stretch_data = diff(this_stretch_times)';
                    end
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
                    if strcmp(variable_name, 'leg_angle_l')
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
                    if strcmp(variable_name, 'leg_angle_r')
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
                    
                    if strcmp(variable_name, 'band_duration')
                        stretch_data = diff(this_stretch_times)';
                    end
                    if strcmp(variable_name, 'pushoff_time')
                        % load events
                        event_data = load(['analysis' filesep makeFileName(this.date, this.subject_id, this.trial_type, this.trial_number, 'events.mat')]);
                        left_pushoff_times = event_data.event_data{strcmp(event_data.event_labels, 'left_pushoff')};
                        right_pushoff_times = event_data.event_data{strcmp(event_data.event_labels, 'right_pushoff')};
                        
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
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
                    if strcmp(variable_name, 'midswing_event_time')
                        % load events
                        event_data = load(['analysis' filesep makeFileName(this.date, this.subject_id, this.trial_type, this.trial_number, 'events.mat')]);
                        left_midswing_times = event_data.event_data{strcmp(event_data.event_labels, 'left_midswing')};
                        right_midswing_times = event_data.event_data{strcmp(event_data.event_labels, 'right_midswing')};
                        
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                % find first left push-off after band start
                                band_start_time = this_stretch_times(i_band);
                                band_end_time = this_stretch_times(i_band+1);
                                this_midswing_time = min(left_midswing_times(left_midswing_times >= band_start_time));
                                if this_midswing_time >= band_end_time
                                    this_midswing_time = band_start_time;
                                end
                                stretch_data(i_band) = this_midswing_time - band_start_time;
                                
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                band_start_time = this_stretch_times(i_band);
                                band_end_time = this_stretch_times(i_band+1);
                                this_midswing_time = min(right_midswing_times(right_midswing_times >= band_start_time));
                                if this_midswing_time >= band_end_time
                                    this_midswing_time = band_start_time;
                                end
                                stretch_data(i_band) = this_midswing_time - band_start_time;
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stretch_data(i_band) = NaN;
                            end
                        end
                    end
                    if strcmp(variable_name, 'midstance_index')
                        LANK_y = this.getTimeNormalizedData('marker:LANK_y', this_stretch_times);
                        RANK_y = this.getTimeNormalizedData('marker:RANK_y', this_stretch_times);
                        pelvis_y = this.getTimeNormalizedData('pelvis_y', this_stretch_times);
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
                        pelvis_y = this.getTimeNormalizedData('pelvis_y', this_stretch_times);
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
                    if strcmp(variable_name, 'left_arm_phase') || strcmp(variable_name, 'right_arm_phase') || strcmp(variable_name, 'left_leg_phase') || strcmp(variable_name, 'right_leg_phase')
     
                        variable_time = this.getTimeData(variable_name);
                        variable_data = this.getBasicVariableData(variable_name);

                        band_time_indices = zeros(size(this_stretch_times));
                        for i_band_time = 1 : length(this_stretch_times)
                            [~, time_index] = min(abs(variable_time - this_stretch_times(i_band_time)));
                            band_time_indices(i_band_time) = time_index;
                        end
                        time_extracted = variable_time(band_time_indices(1) : band_time_indices(end));
                        data_extracted = variable_data(band_time_indices(1) : band_time_indices(end));
                        band_time_indices_local = band_time_indices - band_time_indices(1) + 1;
                                  
                        % make sure there's no leaf change within this data stretch
                        data_extracted_groomed = data_extracted;
                        for i_time = 2 : length(time_extracted)
                            if abs(data_extracted_groomed(i_time) - data_extracted_groomed(i_time-1)) > 5
                                % normalize
                                while data_extracted_groomed(i_time) - data_extracted_groomed(i_time-1) <= -pi
                                    data_extracted_groomed(i_time) = data_extracted_groomed(i_time) + 2*pi;
                                end
                                while data_extracted_groomed(i_time) - data_extracted_groomed(i_time-1) > pi
                                    data_extracted_groomed(i_time) = data_extracted_groomed(i_time) - 2*pi;
                                end
                            end
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
                            data_normalized = spline(time_extracted, data_extracted_groomed, time_normalized);
                        else
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
                            data_normalized = time_normalized * NaN;
                        end
                        
                        % re-normalize angle
                        stretch_data = normalizeAngle(data_normalized);
                    end
                    
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
            if strcmp(variable_name, 'com_rough_x_vel')
                com_rough_x_directions = this.basic_variable_directions.com_rough_x;
                stretch_directions_new = com_rough_x_directions;
            end 
            if strcmp(variable_name, 'com_rough_y_vel')
                com_rough_y_directions = this.basic_variable_directions.com_rough_y;
                stretch_directions_new = com_rough_y_directions;
            end
            if strcmp(variable_name, 'xcom_rough_x')
                com_rough_x_directions = this.basic_variable_directions.com_rough_x;
                stretch_directions_new = com_rough_x_directions;
            end
            if strcmp(variable_name, 'xcom_rough_y')
                com_rough_y_directions = this.basic_variable_directions.com_rough_y;
                stretch_directions_new = com_rough_y_directions;
            end 
            if strcmp(variable_name, 'mos_rough_x')
                com_rough_x_directions = this.basic_variable_directions.com_rough_x;
                stretch_directions_new = com_rough_x_directions;
            end
            if strcmp(variable_name, 'mos_rough_y')
                com_rough_y_directions = this.basic_variable_directions.com_rough_y;
                stretch_directions_new = com_rough_y_directions;
            end
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
            if strcmp(variable_name, 'step_time')
                stretch_directions_new = {'+'; '-'};
            end
            if strcmp(variable_name, 'leg_length_l')
                stretch_directions_new = {'+'; '-'};
            end
            if strcmp(variable_name, 'leg_length_r')
                stretch_directions_new = {'+'; '-'};
            end
            if strcmp(variable_name, 'leg_angle_l')
                stretch_directions_new = {'top leans forward'; 'top leans backward'};
            end
            if strcmp(variable_name, 'leg_angle_r')
                stretch_directions_new = {'top leans forward'; 'top leans backward'};
            end
            if strcmp(variable_name, 'pushoff_time')
                stretch_directions_new = {'+'; '-'};
            end
            if strcmp(variable_name, 'midswing_event_time')
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
            if strcmp(variable_name, 'left_arm_phase') || strcmp(variable_name, 'right_arm_phase') || strcmp(variable_name, 'left_leg_phase') || strcmp(variable_name, 'right_leg_phase')
                % TODO: not tested yet
                stretch_directions_new = this.basic_variable_directions.(variable_name);
            end
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





