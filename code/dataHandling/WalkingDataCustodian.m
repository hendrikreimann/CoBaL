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
            
            % load this information from the subjects.mat and studySettings.txt files
            this.data_directory = data_directory;
            load([data_directory filesep 'subjectInfo.mat'], 'date', 'subject_id');
            this.date = date;
            this.subject_id = subject_id;
            
            this.subject_info = load([data_directory filesep 'subjectInfo.mat']);
            this.subject_settings = SettingsCustodian([data_directory filesep 'subjectSettings.txt']);
            % load settings
            study_settings_file = '';
            if exist([data_directory filesep '..' filesep 'studySettings.txt'], 'file')
                study_settings_file = [data_directory filesep '..' filesep 'studySettings.txt'];
            end    
            if exist([data_directory filesep '..' filesep '..' filesep 'studySettings.txt'], 'file')
                study_settings_file = [data_directory filesep '..' filesep '..' filesep 'studySettings.txt'];
            end
            this.study_settings = SettingsCustodian(study_settings_file);
            emg_normalization_file_name = [data_directory filesep 'analysis' filesep makeFileName(date, subject_id, 'emgNormalization.mat')];
            if exist(emg_normalization_file_name, 'file')
                emg_normalization_data = load(emg_normalization_file_name);
                this.emg_normalization_values = emg_normalization_data.emg_normalization_values;
                this.emg_normalization_labels = emg_normalization_data.emg_variable_names;
            end
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
                    
                    this.addBasicVariable([this_variable_type '_trajectories'])
                    this.addStretchVariable(this_variable_name)
                end
            end
            
            
            % for each possible variable to analyze, list the basic and required variables required to calculate it
            if this.isVariableToAnalyze('marker_trajectories')
                this.addBasicVariable('marker_trajectories')
            end
            if this.isVariableToAnalyze('joint_center_trajectories')
                this.addBasicVariable('joint_center_trajectories')
            end
            if this.isVariableToAnalyze('com_trajectories')
                this.addBasicVariable('com_trajectories')
            end
            
            % kinematics
            if this.isVariableToAnalyze('lheel_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lheel_x')
                this.addStretchVariable('lheel_x')
            end
            if this.isVariableToAnalyze('rheel_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('rheel_x')
                this.addStretchVariable('rheel_x')
            end
%             heel position relative to mpsis
            if this.isVariableToAnalyze('lheel_from_mpsis_initial_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lheel_x')
                this.addBasicVariable('mpsis_x')
                this.addStretchVariable('lheel_from_mpsis_initial_x')                
            end
            if this.isVariableToAnalyze('lheel_from_com_initial_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lheel_x')
                this.addBasicVariable('com')
                this.addBasicVariable('com_x')
                this.addStretchVariable('lheel_from_com_initial_x')                
            end
            if this.isVariableToAnalyze('rheel_from_mpsis_initial_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('rheel_x')
                this.addBasicVariable('mpsis_x')
                this.addStretchVariable('rheel_from_mpsis_initial_x')
            end
            if this.isVariableToAnalyze('rheel_from_com_initial_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('rheel_x')
                this.addBasicVariable('com')
                this.addBasicVariable('com_x')
                this.addStretchVariable('rheel_from_com_initial_x')
            end
            if this.isVariableToAnalyze('mpsis_from_mpsis_initial_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('mpsis_x')
                this.addStretchVariable('mpsis_from_mpsis_initial_x')
            end            
            if this.isVariableToAnalyze('com_from_com_initial_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('com')
                this.addBasicVariable('com_x')
                this.addStretchVariable('com_from_com_initial_x')
            end    
            if this.isVariableToAnalyze('lheel_y')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lheel_y')
                this.addStretchVariable('lheel_y')
            end
            if this.isVariableToAnalyze('rheel_y')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('rheel_y')
                this.addStretchVariable('rheel_y')
            end
            if this.isVariableToAnalyze('lankle_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lankle_x')
                this.addStretchVariable('lankle_x')
            end
            if this.isVariableToAnalyze('rankle_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('rankle_x')
                this.addStretchVariable('rankle_x')
            end
            if this.isVariableToAnalyze('lankle_y')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lankle_y')
                this.addStretchVariable('lankle_y')
            end
            if this.isVariableToAnalyze('rankle_y')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('rankle_y')
                this.addStretchVariable('rankle_y')
            end    
            if this.isVariableToAnalyze('ltoel_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('ltoel_x')
                this.addStretchVariable('ltoel_x')               
            end
            if this.isVariableToAnalyze('step_length')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lheel_y')
                this.addBasicVariable('rheel_y')
                this.addStretchVariable('lheel_y')
                this.addStretchVariable('rheel_y')
                this.addStretchVariable('step_length')
            end
            if this.isVariableToAnalyze('step_width')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lheel_x')
                this.addBasicVariable('rheel_x')
                this.addStretchVariable('lheel_x')
                this.addStretchVariable('rheel_x')
                this.addStretchVariable('step_width')
            end
            if this.isVariableToAnalyze('step_placement_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lheel_x')
                this.addBasicVariable('rheel_x')
                this.addStretchVariable('lheel_x')
                this.addStretchVariable('rheel_x')
                this.addStretchVariable('step_placement_x')
            end
            if this.isVariableToAnalyze('stimulus_response_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lheel_x')
                this.addBasicVariable('rheel_x')
                this.addBasicVariable('lankle_y')
                this.addBasicVariable('rankle_y')
                this.addBasicVariable('pelvis_y')
                this.addStretchVariable('step_time')
                this.addStretchVariable('pushoff_time')
                this.addStretchVariable('midstance_index')
                this.addStretchVariable('lheel_x')
                this.addStretchVariable('rheel_x')
                this.addStretchVariable('step_placement_x')
            end
            if this.isVariableToAnalyze('step_time')
                this.addStretchVariable('step_time')
            end
            if this.isVariableToAnalyze('step_duration')
                this.addStretchVariable('step_time')
                this.addStretchVariable('step_duration')
            end
            if this.isVariableToAnalyze('pushoff_time')
                this.addBasicVariable('marker_trajectories')
                this.addStretchVariable('step_time')
                this.addStretchVariable('pushoff_time')
            end
            if this.isVariableToAnalyze('midstance_index')
                this.addBasicVariable('lankle_y')
                this.addBasicVariable('rankle_y')
                this.addBasicVariable('pelvis_y')
                this.addStretchVariable('step_time')
                this.addStretchVariable('pushoff_time')
                this.addStretchVariable('midstance_index')
            end
            if this.isVariableToAnalyze('cadence')
                this.addStretchVariable('step_time')
                this.addStretchVariable('cadence')
            end
            if this.isVariableToAnalyze('velocity')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lheel_y')
                this.addBasicVariable('rheel_y')
                this.addStretchVariable('lheel_y')
                this.addStretchVariable('rheel_y')
                this.addStretchVariable('step_length')
                this.addStretchVariable('step_time')
                this.addStretchVariable('velocity')
            end
            if this.isVariableToAnalyze('cop_from_mpsis_x')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_x')
                this.addStretchVariable('cop_x')
                this.addBasicVariable('mpsis_x')
                this.addStretchVariable('mpsis_x')
                this.addStretchVariable('cop_from_mpsis_x')
            end
            if this.isVariableToAnalyze('cop_from_mpsis_y')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_y')
                this.addStretchVariable('cop_y')
                this.addBasicVariable('mpsis_y')
                this.addStretchVariable('mpsis_y')
                this.addStretchVariable('cop_from_mpsis_y')
            end
            if this.isVariableToAnalyze('mpsis_x')
                this.addBasicVariable('mpsis_x')
                this.addStretchVariable('mpsis_x')
            end
            if this.isVariableToAnalyze('mpsis_y')
                this.addBasicVariable('mpsis_y')
                this.addStretchVariable('mpsis_y')
            end
            if this.isVariableToAnalyze('mpsis_x_vel')
                this.addBasicVariable('mpsis_x')
                this.addBasicVariable('mpsis_x_vel')
                this.addStretchVariable('mpsis_x_vel')
            end
            if this.isVariableToAnalyze('mpsis_y_vel')
                this.addBasicVariable('mpsis_y')
                this.addBasicVariable('mpsis_y_vel')
                this.addStretchVariable('mpsis_y_vel')
            end
            if this.isVariableToAnalyze('cop_from_com_x')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_x')
                this.addStretchVariable('cop_x')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_x')
                this.addStretchVariable('com_x')
                this.addStretchVariable('cop_from_com_x')
            end
            if this.isVariableToAnalyze('cop_from_com_y')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_y')
                this.addStretchVariable('cop_y')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_y')
                this.addStretchVariable('com_y')
                this.addStretchVariable('cop_from_com_y')
            end
            if this.isVariableToAnalyze('cop_to_com_x')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_x')
                this.addStretchVariable('cop_x')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_x')
                this.addStretchVariable('com_x')
                this.addStretchVariable('cop_to_com_x')
            end
            if this.isVariableToAnalyze('cop_to_com_y')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_y')
                this.addStretchVariable('cop_y')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_y')
                this.addStretchVariable('com_y')
                this.addStretchVariable('cop_to_com_y')
            end
            if this.isVariableToAnalyze('init_heel_to_com_x')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lheel_x')
                this.addStretchVariable('lheel_x')
                this.addBasicVariable('rheel_x')
                this.addStretchVariable('rheel_x')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_x')
                this.addStretchVariable('com_x')
                this.addStretchVariable('init_heel_to_com_x')
            end
            if this.isVariableToAnalyze('cop_to_com_vel_scaled_x')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_x')
                this.addBasicVariable('com_x_vel')
                this.addStretchVariable('com_x_vel')
                this.addBasicVariable('com_z')
                this.addStretchVariable('com_z')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_x')
                this.addBasicVariable('cop_x_vel')
                this.addStretchVariable('cop_x')
                this.addStretchVariable('com_x')
                this.addStretchVariable('cop_x_vel')
                this.addStretchVariable('cop_to_com_vel_scaled_x')
            end
            if this.isVariableToAnalyze('cop_to_com_vel_scaled_y')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_y')
                this.addBasicVariable('com_y_vel')
                this.addStretchVariable('com_y_vel')
                this.addBasicVariable('com_z')
                this.addStretchVariable('com_z')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_y')
                this.addBasicVariable('cop_y_vel')
                this.addStretchVariable('cop_y_vel')
                this.addStretchVariable('cop_to_com_vel_scaled_y')
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
            if this.isVariableToAnalyze('com_x')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_x')
                this.addStretchVariable('com_x')
            end
            if this.isVariableToAnalyze('com_y')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_y')
                this.addStretchVariable('com_y')
            end
            if this.isVariableToAnalyze('com_z')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com_z')
                this.addStretchVariable('com_z')
            end
            if this.isVariableToAnalyze('com_x_vel')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_x')
                this.addBasicVariable('com_x_vel')
                this.addStretchVariable('com_x_vel')
            end
            if this.isVariableToAnalyze('com_x_vel_scaled')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_x')
                this.addBasicVariable('com_x_vel')
                this.addStretchVariable('com_x_vel')
                this.addBasicVariable('com_z')
                this.addStretchVariable('com_z')
                this.addStretchVariable('com_x_vel_scaled')
            end
            if this.isVariableToAnalyze('com_y_vel')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_y')
                this.addBasicVariable('com_y_vel')
                this.addStretchVariable('com_y_vel')
            end
            if this.isVariableToAnalyze('com_z_vel')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_z')
                this.addBasicVariable('com_z_vel')
                this.addStretchVariable('com_z_vel')
            end
            if this.isVariableToAnalyze('com_x_acc')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_x')
                this.addBasicVariable('com_x_vel')
                this.addBasicVariable('com_x_acc')
                this.addStretchVariable('com_x_acc')
            end
            if this.isVariableToAnalyze('com_y_acc')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_y')
                this.addBasicVariable('com_y_vel')
                this.addBasicVariable('com_y_acc')
                this.addStretchVariable('com_y_acc')
            end
            if this.isVariableToAnalyze('com_z_acc')
%                 this.addBasicVariable('com_trajectories') % in-house kinematics
                this.addBasicVariable('com_position_trajectories') % opensim kinematics
                this.addBasicVariable('com')
                this.addBasicVariable('com_z')
                this.addBasicVariable('com_z_vel')
                this.addBasicVariable('com_z_acc')
                this.addStretchVariable('com_z_acc')
            end
            if this.isVariableToAnalyze('heel_clearance')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('lheel_y')
                this.addBasicVariable('lheel_z')
                this.addBasicVariable('rheel_y')
                this.addBasicVariable('rheel_z')
                this.addStretchVariable('heel_clearance')
            end
            if this.isVariableToAnalyze('toes_clearance')
                this.addBasicVariable('marker_trajectories')
                this.addBasicVariable('ltoes_y')
                this.addBasicVariable('ltoes_z')
                this.addBasicVariable('rtoes_y')
                this.addBasicVariable('rtoes_z')
                this.addStretchVariable('toes_clearance')
            end
            % force plate
            if this.isVariableToAnalyze('copl_y')
                this.addBasicVariable('left_foot_cop_world')
                this.addBasicVariable('copl_y')
                this.addStretchVariable('copl_y')
            end
            if this.isVariableToAnalyze('copl_x')
                this.addBasicVariable('left_foot_cop_world')
                this.addBasicVariable('copl_x')
                this.addStretchVariable('copl_x')
            end
            if this.isVariableToAnalyze('copr_y')
                this.addBasicVariable('right_foot_cop_world')
                this.addBasicVariable('copr_y')
                this.addStretchVariable('copr_y')
            end
            if this.isVariableToAnalyze('copr_x')
                this.addBasicVariable('right_foot_cop_world')
                this.addBasicVariable('copr_x')
                this.addStretchVariable('copr_x')
            end
            if this.isVariableToAnalyze('cop_y')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_y')
                this.addStretchVariable('cop_y')
            end
            if this.isVariableToAnalyze('cop_x')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_x')
                this.addStretchVariable('cop_x')
            end
            if this.isVariableToAnalyze('cop_y_vel')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_y')
                this.addBasicVariable('cop_y_vel')
                this.addStretchVariable('cop_y_vel')
            end
            if this.isVariableToAnalyze('cop_y_acc')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_y')
                this.addBasicVariable('cop_y_vel')
                this.addBasicVariable('cop_y_acc')
                this.addStretchVariable('cop_y_acc')
            end
            if this.isVariableToAnalyze('cop_x_vel')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_x')
                this.addBasicVariable('cop_x_vel')
                this.addStretchVariable('cop_x_vel')
            end
            if this.isVariableToAnalyze('cop_x_acc')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_x')
                this.addBasicVariable('cop_x_vel')
                this.addBasicVariable('cop_x_acc')
                this.addStretchVariable('cop_x_acc')
            end
            if this.isVariableToAnalyze('fxl')
                this.addBasicVariable('left_foot_wrench_world')
                this.addBasicVariable('fxl')
                this.addStretchVariable('fxl')
            end
            if this.isVariableToAnalyze('fyl')
                this.addBasicVariable('left_foot_wrench_world')
                this.addBasicVariable('fyl')
                this.addStretchVariable('fyl')
            end
            if this.isVariableToAnalyze('fzl')
                this.addBasicVariable('left_foot_wrench_world')
                this.addBasicVariable('fzl')
                this.addStretchVariable('fzl')
            end
            if this.isVariableToAnalyze('mxl')
                this.addBasicVariable('left_foot_wrench_world')
                this.addBasicVariable('mxl')
                this.addStretchVariable('mxl')
            end
            if this.isVariableToAnalyze('myl')
                this.addBasicVariable('left_foot_wrench_world')
                this.addBasicVariable('myl')
                this.addStretchVariable('myl')
            end
            if this.isVariableToAnalyze('mzl')
                this.addBasicVariable('left_foot_wrench_world')
                this.addBasicVariable('mzl')
                this.addStretchVariable('mzl')
            end
            if this.isVariableToAnalyze('fxr')
                this.addBasicVariable('right_foot_wrench_world')
                this.addBasicVariable('fxr')
                this.addStretchVariable('fxr')
            end
            if this.isVariableToAnalyze('fyr')
                this.addBasicVariable('right_foot_wrench_world')
                this.addBasicVariable('fyr')
                this.addStretchVariable('fyr')
            end
            if this.isVariableToAnalyze('fzr')
                this.addBasicVariable('right_foot_wrench_world')
                this.addBasicVariable('fzr')
                this.addStretchVariable('fzr')
            end
            if this.isVariableToAnalyze('mxr')
                this.addBasicVariable('right_foot_wrench_world')
                this.addBasicVariable('mxr')
                this.addStretchVariable('mxr')
            end
            if this.isVariableToAnalyze('myr')
                this.addBasicVariable('right_foot_wrench_world')
                this.addBasicVariable('myr')
                this.addStretchVariable('myr')
            end
            if this.isVariableToAnalyze('mzr')
                this.addBasicVariable('right_foot_wrench_world')
                this.addBasicVariable('mzr')
                this.addStretchVariable('mzr')
            end
            if this.isVariableToAnalyze('fx')
                this.addBasicVariable('total_forceplate_wrench_world')
                this.addBasicVariable('fx')
                this.addStretchVariable('fx')
            end
            if this.isVariableToAnalyze('fy')
                this.addBasicVariable('total_forceplate_wrench_world')
                this.addBasicVariable('fy')
                this.addStretchVariable('fy')
            end
            if this.isVariableToAnalyze('fz')
                this.addBasicVariable('total_forceplate_wrench_world')
                this.addBasicVariable('fz')
                this.addStretchVariable('fz')
            end
            if this.isVariableToAnalyze('mx')
                this.addBasicVariable('total_forceplate_wrench_world')
                this.addBasicVariable('mx')
                this.addStretchVariable('mx')
            end
            if this.isVariableToAnalyze('my')
                this.addBasicVariable('total_forceplate_wrench_world')
                this.addBasicVariable('my')
                this.addStretchVariable('my')
            end
            if this.isVariableToAnalyze('mz')
                this.addBasicVariable('total_forceplate_wrench_world')
                this.addBasicVariable('mz')
                this.addStretchVariable('mz')
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
                name_to_use = [this_variable_split{1} '_trajectories'];
            else
                name_to_use = variable_name;
            end
            
%             eval(['time_data = this.time_data.' name_to_use ';']);
            time_data = this.time_data.(name_to_use);
        end
        function [variable_data, variable_directions] = getBasicVariableData(this, variable_name)
            % unpack name
            if any(variable_name==':')
                this_variable_split = strsplit(variable_name, ':');
                name_to_use = [this_variable_split{1} '_trajectories'];
                label_to_use = this_variable_split{2};
                trajectory_data = this.basic_variable_data.(name_to_use);
                trajectory_labels = this.basic_variable_labels.(name_to_use);
                variable_data = trajectory_data(:, strcmp(trajectory_labels, label_to_use));
                
            else
                name_to_use = variable_name;
%                 eval(['variable_data = this.basic_variable_data.' name_to_use ';']);
                variable_data = this.basic_variable_data.(variable_name);
            end
            
            if nargout > 1
%                 eval(['variable_directions = this.basic_variable_directions.' name_to_use ';']);
                variable_directions = this.basic_variable_directions.(name_to_use);
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
            load(['analysis' filesep makeFileName(this.subject_info.date, this.subject_info.subject_id, trial_type, trial_number, 'availableVariables')], 'available_variables');
            
            % load basic variables
            for i_variable = 1 : length(variables_to_prepare)
                variable_name = variables_to_prepare{i_variable};
                
                % try loading
                [data, time, sampling_rate, labels, directions, success] = loadData(this.subject_info.date, this.subject_info.subject_id, trial_type, trial_number, variable_name, 'optional'); %#ok<ASGLU>
                
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
                    if any(strcmp(this.basic_variable_labels.([this_variable_type '_trajectories']), this_variable_label))
                        % found this label in the data, so all is good
                        success = 1;
                    end
                end

                
                % calculate variables that can't be loaded
                
                % kinematics
                if strcmp(variable_name, 'lheel_x')
                    LHEE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE');
                    this.basic_variable_data.lheel_x = LHEE_trajectory(:, 1);
                    LHEE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE',  'indices');
                    this.basic_variable_directions.lheel_x = this.basic_variable_directions.marker_trajectories(:, LHEE_indices(1));
                    this.time_data.lheel_x = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'rheel_x')
                    RHEE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE');
                    this.basic_variable_data.rheel_x = RHEE_trajectory(:, 1);
                    RHEE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE',  'indices');
                    this.basic_variable_directions.rheel_x = this.basic_variable_directions.marker_trajectories(:, RHEE_indices(1));
                    this.time_data.rheel_x = this.time_data.marker_trajectories;
                    success = 1;
                end 
                if strcmp(variable_name, 'lheel_y')
                    LHEE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE');
                    this.basic_variable_data.lheel_y = LHEE_trajectory(:, 2);
                    LHEE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE',  'indices');
                    this.basic_variable_directions.lheel_y = this.basic_variable_directions.marker_trajectories(:, LHEE_indices(2));
                    this.time_data.lheel_y = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'rheel_y')
                    RHEE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE');
                    this.basic_variable_data.rheel_y = RHEE_trajectory(:, 2);
                    RHEE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE',  'indices');
                    this.basic_variable_directions.rheel_y = this.basic_variable_directions.marker_trajectories(:, RHEE_indices(2));
                    this.time_data.rheel_y = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'lheel_z')
                    LHEE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE');
                    this.basic_variable_data.lheel_z = LHEE_trajectory(:, 3);
                    LHEE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE',  'indices');
                    this.basic_variable_directions.lheel_z = this.basic_variable_directions.marker_trajectories(:, LHEE_indices(3));
                    this.time_data.lheel_z = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'rheel_z')
                    RHEE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE');
                    this.basic_variable_data.rheel_z = RHEE_trajectory(:, 3);
                    RHEE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE',  'indices');
                    this.basic_variable_directions.rheel_z = this.basic_variable_directions.marker_trajectories(:, RHEE_indices(3));
                    this.time_data.rheel_z = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'ltoes_y')
                    LTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE');
                    this.basic_variable_data.ltoes_y = LTOE_trajectory(:, 2);
                    LTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE',  'indices');
                    this.basic_variable_directions.ltoes_y = this.basic_variable_directions.marker_trajectories(:, LTOE_indices(2));
                    this.time_data.ltoes_y = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'rtoes_y')
                    RTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE');
                    this.basic_variable_data.rtoes_y = RTOE_trajectory(:, 2);
                    RTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE',  'indices');
                    this.basic_variable_directions.rtoes_y = this.basic_variable_directions.marker_trajectories(:, RTOE_indices(2));
                    this.time_data.rtoes_y = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'ltoes_z')
                    LTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE');
                    this.basic_variable_data.ltoes_z = LTOE_trajectory(:, 3);
                    LTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOE',  'indices');
                    this.basic_variable_directions.ltoes_z = this.basic_variable_directions.marker_trajectories(:, LTOE_indices(3));
                    this.time_data.ltoes_z = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'rtoes_z')
                    RTOE_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE');
                    this.basic_variable_data.rtoes_z = RTOE_trajectory(:, 3);
                    RTOE_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RTOE',  'indices');
                    this.basic_variable_directions.rtoes_z = this.basic_variable_directions.marker_trajectories(:, RTOE_indices(3));
                    this.time_data.rtoes_z = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'lankle_x')
                    LANK_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LANK');
                    this.basic_variable_data.lankle_x = LANK_trajectory(:, 1);
                    LANK_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LANK',  'indices');
                    this.basic_variable_directions.lankle_x = this.basic_variable_directions.marker_trajectories(:, LANK_indices(1));
                    this.time_data.lankle_x = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'rankle_x')
                    RANK_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RANK');
                    this.basic_variable_data.rankle_x = RANK_trajectory(:, 1);
                    RANK_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RANK',  'indices');
                    this.basic_variable_directions.rankle_x = this.basic_variable_directions.marker_trajectories(:, RANK_indices(1));
                    this.time_data.rankle_x = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'lankle_y')
                    LANK_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LANK');
                    this.basic_variable_data.lankle_y = LANK_trajectory(:, 2);
                    LANK_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LANK',  'indices');
                    this.basic_variable_directions.lankle_y = this.basic_variable_directions.marker_trajectories(:, LANK_indices(2));
                    this.time_data.lankle_y = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'rankle_y')
                    RANK_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RANK');
                    this.basic_variable_data.rankle_y = RANK_trajectory(:, 2);
                    RANK_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RANK',  'indices');
                    this.basic_variable_directions.rankle_y = this.basic_variable_directions.marker_trajectories(:, RANK_indices(2));
                    this.time_data.rankle_y = this.time_data.marker_trajectories;
                    success = 1;
                end
                 if strcmp(variable_name, 'ltoel_x')
                    LTOEL_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOEL');
                    this.basic_variable_data.ltoe_x = LTOEL_trajectory(:, 1);
                    LTOEL_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LTOEL',  'indices');
                    this.basic_variable_directions.ltoel_x = this.basic_variable_directions.marker_trajectories(:, LTOEL_indices(1));
                    this.time_data.ltoel_x = this.time_data.marker_trajectories;
                    success = 1;             
                 end
                if strcmp(variable_name, 'lasis_y')
                    LASI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LASI');
                    this.basic_variable_data.lasis_y = LASI_trajectory(:, 2);
                    LASI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LASI',  'indices');
                    this.basic_variable_directions.lasis_y = this.basic_variable_directions.marker_trajectories(:, LASI_indices(2));
                    this.time_data.lasis_y = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'rasis_y')
                    RASI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RASI');
                    this.basic_variable_data.rasis_y = RASI_trajectory(:, 2);
                    RASI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RASI',  'indices');
                    this.basic_variable_directions.rasis_y = this.basic_variable_directions.marker_trajectories(:, RASI_indices(2));
                    this.time_data.rasis_y = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'mpsis_x')
                    LPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI');
                    RPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI');
                    MPSI_trajectory = (LPSI_trajectory + RPSI_trajectory) * 0.5;
                    this.basic_variable_data.mpsis_x = MPSI_trajectory(:, 1);
                    LPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI',  'indices');
                    RPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI',  'indices');
                    LPSI_directions = this.basic_variable_directions.marker_trajectories(:, LPSI_indices(1));
                    RPSI_directions = this.basic_variable_directions.marker_trajectories(:, RPSI_indices(1));
                    if ~strcmp(LPSI_directions{1}, RPSI_directions{1})
                        error('LPSI and RPSI directions found in marker data are different from each other')
                    end
                    if ~strcmp(LPSI_directions{2}, RPSI_directions{2})
                        error('LPSI and RPSI directions found in marker data are different from each other')
                    end
                    this.basic_variable_directions.mpsis_x = LPSI_directions;
                    this.time_data.mpsis_x = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'mpsis_y')
                    LPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI');
                    RPSI_trajectory = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI');
                    MPSI_trajectory = (LPSI_trajectory + RPSI_trajectory) * 0.5;
                    this.basic_variable_data.mpsis_y = MPSI_trajectory(:, 2);
                    LPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI',  'indices');
                    RPSI_indices = extractMarkerData(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI',  'indices');
                    LPSI_directions = this.basic_variable_directions.marker_trajectories(:, LPSI_indices(2));
                    RPSI_directions = this.basic_variable_directions.marker_trajectories(:, RPSI_indices(2));
                    if ~strcmp(LPSI_directions{1}, RPSI_directions{1})
                        error('LPSI and RPSI directions found in marker data are different from each other')
                    end
                    if ~strcmp(LPSI_directions{2}, RPSI_directions{2})
                        error('LPSI and RPSI directions found in marker data are different from each other')
                    end
                    this.basic_variable_directions.mpsis_y = LPSI_directions;
                    this.time_data.mpsis_y = this.time_data.marker_trajectories;
                    success = 1;
                end
                if strcmp(variable_name, 'mpsis_x_vel')
                    mpsis_x = this.getBasicVariableData('mpsis_x');
                    mpsis_x(mpsis_x==0) = NaN;
                    time = this.getTimeData('mpsis_x');
                    filter_order = this.study_settings.get('filter_order_mpsis_vel');
                    cutoff_frequency = this.study_settings.get('filter_cutoff_mpsis_vel');
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));
                    if any(~isnan(mpsis_x))
                        mpsis_x_vel = deriveByTime(nanfiltfilt(b, a, mpsis_x), 1/sampling_rate);
                    else
                        mpsis_x_vel = ones(size(mpsis_x)) * NaN;
                    end
                    this.basic_variable_data.mpsis_x_vel = mpsis_x_vel;
                    this.basic_variable_directions.mpsis_x_vel = this.basic_variable_directions.mpsis_x;
                    this.time_data.mpsis_x_vel = time;
                    success = 1;
                end
                if strcmp(variable_name, 'mpsis_y_vel')
                    mpsis_y = this.getBasicVariableData('mpsis_y');
                    mpsis_y(mpsis_y==0) = NaN;
                    time = this.getTimeData('mpsis_y');
                    filter_order = this.study_settings.get('filter_order_mpsis_vel');
                    cutoff_frequency = this.study_settings.get('filter_cutoff_mpsis_vel');
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;
                    if any(~isnan(mpsis_y))
                        mpsis_y_vel = deriveByTime(nanfiltfilt(b, a, mpsis_y), 1/sampling_rate);
                    else
                        mpsis_y_vel = ones(size(mpsis_y)) * NaN;
                    end
                    this.basic_variable_data.mpsis_y_vel = mpsis_y_vel;
                    this.basic_variable_directions.mpsis_y_vel = this.basic_variable_directions.mpsis_y;
                    this.time_data.mpsis_y_vel = time;
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
                    this.basic_variable_data.head_angle_ap = rad2deg(atan2(head_vector_y, head_vector_z));
                    
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
                    this.basic_variable_data.head_angle_ml = rad2deg(atan2(head_vector_x, head_vector_z));
                    
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
                    this.basic_variable_data.trunk_angle_ap = rad2deg(atan2(trunk_vector_y, trunk_vector_z));
                    
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
                    this.basic_variable_data.trunk_angle_ml = rad2deg(atan2(trunk_vector_x, trunk_vector_z));
                    
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
                    this.basic_variable_data.pelvis_angle_ml = -rad2deg(atan2(pelvis_vector_z, pelvis_vector_x));
                    
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
                    this.basic_variable_data.left_leg_angle_ml = rad2deg(atan2(left_leg_vector_x, left_leg_vector_z));
                    
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
                    this.basic_variable_data.right_leg_angle_ml = rad2deg(atan2(right_leg_vector_x, right_leg_vector_z));
                    
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
                    this.basic_variable_data.left_leg_angle_ap = -rad2deg(atan2(left_leg_vector_y, left_leg_vector_z));
                    
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
                    this.basic_variable_data.right_leg_angle_ap = -rad2deg(atan2(right_leg_vector_y, right_leg_vector_z));
                    
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
                    this.basic_variable_data.left_foot_angle_ap = rad2deg(atan2(foot_vector_z, foot_vector_xy));
                    
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
                    this.basic_variable_data.left_foot_angle_ml = rad2deg(atan2(foot_vector_z, foot_vector_xy));
                    
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
                    this.basic_variable_data.left_foot_angle_yaw = rad2deg(atan2(foot_vector_x, foot_vector_y));
                    
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
                    this.basic_variable_data.right_foot_angle_ap = rad2deg(atan2(foot_vector_z, foot_vector_xy));
                    
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
                    this.basic_variable_data.right_foot_angle_ml = -rad2deg(atan2(foot_vector_z, foot_vector_xy));
                    
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
                    this.basic_variable_data.right_foot_angle_yaw = rad2deg(atan2(foot_vector_x, foot_vector_y));
                    
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
                
                if strcmp(variable_name, 'com')
                    if dataIsAvailable(this.subject_info.date, this.subject_info.subject_id, trial_type, trial_number, 'com_trajectories')
                        % in-house kinematics - TODO: not tested yet
                        [data, time, ~, labels, directions] = loadData(this.subject_info.date, this.subject_info.subject_id, trial_type, trial_number, 'com_trajectories');
                        this.basic_variable_data.com = extractMarkerData(data, labels, 'BODYCOM');
                        com_indices = extractMarkerData(data, labels, 'BODYCOM',  'indices');
                        this.basic_variable_directions.com = directions(:, com_indices);
                        this.time_data.com = time;
                        success = 1;
                    elseif dataIsAvailable(this.subject_info.date, this.subject_info.subject_id, trial_type, trial_number, 'com_position_trajectories')
                        [data, time, ~, labels, directions] = loadData(this.subject_info.date, this.subject_info.subject_id, trial_type, trial_number, 'com_position_trajectories');
                        this.basic_variable_data.com = extractMarkerData(data, labels, 'center_of_mass');
                        com_indices = extractMarkerData(data, labels, 'center_of_mass',  'indices');
                        this.basic_variable_directions.com = directions(:, com_indices);
                        this.time_data.com = time;
                        success = 1;
                    else
                        error('CoM data needed, but found neither in-house nor OpenSim data')
                    end
                end
                if strcmp(variable_name, 'com_x')
                    this.basic_variable_data.com_x = this.basic_variable_data.com(:, 1);
                    this.basic_variable_directions.com_x = this.basic_variable_directions.com(:, 1);
                    this.time_data.com_x = this.time_data.com;
                    success = 1;
                end
                if strcmp(variable_name, 'com_y')
                    this.basic_variable_data.com_y = this.basic_variable_data.com(:, 2);
                    this.basic_variable_directions.com_y = this.basic_variable_directions.com(:, 2);
                    this.time_data.com_y = this.time_data.com;
                    success = 1;
                end
                if strcmp(variable_name, 'com_z')
                    this.basic_variable_data.com_z = this.basic_variable_data.com(:, 3);
                    this.basic_variable_directions.com_z = this.basic_variable_directions.com(:, 3);
                    this.time_data.com_z = this.time_data.com;
                    success = 1;
                end
                if strcmp(variable_name, 'com_x_vel')
                    com_x = this.getBasicVariableData('com_x');
                    com_x(com_x==0) = NaN;
                    time = this.getTimeData('com_x');
                    filter_order = this.study_settings.get('filter_order_com_vel');
                    cutoff_frequency = this.study_settings.get('filter_cutoff_com_vel');
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));
                    if any(~isnan(com_x))
                        com_x_vel = deriveByTime(nanfiltfilt(b, a, com_x), 1/sampling_rate);
                    else
                        com_x_vel = ones(size(com_x)) * NaN;
                    end
                    this.basic_variable_data.com_x_vel = com_x_vel;
                    this.basic_variable_directions.com_x_vel = this.basic_variable_directions.com_x;
                    this.time_data.com_x_vel = time;
                    success = 1;
                end
                if strcmp(variable_name, 'com_y_vel')
                    com_y = this.getBasicVariableData('com_y');
                    com_y(com_y==0) = NaN;
                    time = this.getTimeData('com_y');
                    filter_order = this.study_settings.get('filter_order_com_vel');
                    cutoff_frequency = this.study_settings.get('filter_cutoff_com_vel');
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;
                    if any(~isnan(com_y))
                        com_y_vel = deriveByTime(nanfiltfilt(b, a, com_y), 1/sampling_rate);
                    else
                        com_y_vel = ones(size(com_y)) * NaN;
                    end
                    this.basic_variable_data.com_y_vel = com_y_vel;
                    this.basic_variable_directions.com_y_vel = this.basic_variable_directions.com_y;
                    this.time_data.com_y_vel = time;
                    success = 1;
                end
                if strcmp(variable_name, 'com_z_vel')
                    com_z = this.getBasicVariableData('com_z');
                    com_z(com_z==0) = NaN;
                    time = this.getTimeData('com_z');
                    filter_order = this.study_settings.get('filter_order_com_vel');
                    cutoff_frequency = this.study_settings.get('filter_cutoff_com_vel');
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;
                    if any(~isnan(com_z))
                        com_z_vel = deriveByTime(nanfiltfilt(b, a, com_z), 1/sampling_rate);
                    else
                        com_z_vel = ones(size(com_z)) * NaN;
                    end
                    this.basic_variable_data.com_z_vel = com_z_vel;
                    this.basic_variable_directions.com_z_vel = this.basic_variable_directions.com_z;
                    this.time_data.com_z_vel = time;
                    success = 1;
                end
                if strcmp(variable_name, 'com_x_acc')
                    com_x_vel = this.getBasicVariableData('com_x_vel');
                    time = this.getTimeData('com_x_vel');
                    filter_order = this.study_settings.get('filter_order_com_acc');
                    cutoff_frequency = this.study_settings.get('filter_cutoff_com_acc');
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));
                    if any(~isnan(com_x_vel))
                        com_x_acc = deriveByTime(nanfiltfilt(b, a, com_x_vel), 1/sampling_rate);
                    else
                        com_x_acc = ones(size(com_x_vel)) * NaN;
                    end
                    this.basic_variable_data.com_x_acc = com_x_acc;
                    this.basic_variable_directions.com_x_acc = this.basic_variable_directions.com_x_vel;
                    this.time_data.com_x_acc = time;
                    success = 1;
                end
                if strcmp(variable_name, 'com_y_acc')
                    com_y_vel = this.getBasicVariableData('com_y_vel');
                    time = this.getTimeData('com_y_vel');
                    filter_order = this.study_settings.get('filter_order_com_acc');
                    cutoff_frequency = this.study_settings.get('filter_cutoff_com_acc');
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;
                    if any(~isnan(com_y_vel))
                        com_y_acc = deriveByTime(nanfiltfilt(b, a, com_y_vel), 1/sampling_rate);
                    else
                        com_y_acc = ones(size(com_y_vel)) * NaN;
                    end
                    this.basic_variable_data.com_y_acc = com_y_acc;
                    this.basic_variable_directions.com_y_acc = this.basic_variable_directions.com_y_vel;
                    this.time_data.com_y_acc = time;
                    success = 1;
                end
                if strcmp(variable_name, 'com_z_acc')
                    com_z_vel = this.getBasicVariableData('com_z_vel');
                    time = this.getTimeData('com_z_vel');
                    filter_order = this.study_settings.get('filter_order_com_acc');
                    cutoff_frequency = this.study_settings.get('filter_cutoff_com_acc');
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;
                    if any(~isnan(com_z_vel))
                        com_z_acc = deriveByTime(nanfiltfilt(b, a, com_z_vel), 1/sampling_rate);
                    else
                        com_z_acc = ones(size(com_z_vel)) * NaN;
                    end
                    this.basic_variable_data.com_z_acc = com_z_acc;
                    this.basic_variable_directions.com_z_acc = this.basic_variable_directions.com_z_vel;
                    this.time_data.com_z_acc = time;
                    success = 1;
                end
            
                % force plate
                if strcmp(variable_name, 'copl_y')
                    [left_foot_cop_world_data, left_foot_cop_world_directions] = this.getBasicVariableData('left_foot_cop_world');
                    this.basic_variable_data.copl_y = left_foot_cop_world_data(:, 2);
                    this.basic_variable_directions.copl_y = left_foot_cop_world_directions(:, 2);
                    this.time_data.copl_y = this.time_data.left_foot_cop_world;
                    success = 1;
                end
                if strcmp(variable_name, 'copl_x')
                    [left_foot_cop_world_data, left_foot_cop_world_directions] = this.getBasicVariableData('left_foot_cop_world');
                    this.basic_variable_data.copl_x = left_foot_cop_world_data(:, 1);
                    this.basic_variable_directions.copl_x = left_foot_cop_world_directions(:, 1);
                    this.time_data.copl_x = this.time_data.left_foot_cop_world;
                    success = 1;
                end
                if strcmp(variable_name, 'copr_y')
                    [right_foot_cop_world_data, right_foot_cop_world_directions] = this.getBasicVariableData('right_foot_cop_world');
                    this.basic_variable_data.copr_y = right_foot_cop_world_data(:, 2);
                    this.basic_variable_directions.copr_y = right_foot_cop_world_directions(:, 2);
                    this.time_data.copr_y = this.time_data.right_foot_cop_world;
                    success = 1;
                end
                if strcmp(variable_name, 'copr_x')
                    [right_foot_cop_world_data, right_foot_cop_world_directions] = this.getBasicVariableData('right_foot_cop_world');
                    this.basic_variable_data.copr_x = right_foot_cop_world_data(:, 1);
                    this.basic_variable_directions.copr_x = right_foot_cop_world_directions(:, 1);
                    this.time_data.copr_x = this.time_data.right_foot_cop_world;
                    success = 1;
                end
                if strcmp(variable_name, 'cop_y')
                    [total_forceplate_cop_world_data, total_forceplate_cop_world_directions] = this.getBasicVariableData('total_forceplate_cop_world');
                    this.basic_variable_data.cop_y = total_forceplate_cop_world_data(:, 2);
                    this.basic_variable_directions.cop_y = total_forceplate_cop_world_directions(:, 2);
                    this.time_data.cop_y = this.time_data.total_forceplate_cop_world;
                    success = 1;
                end
                if strcmp(variable_name, 'cop_x')
                    [total_forceplate_cop_world_data, total_forceplate_cop_world_directions] = this.getBasicVariableData('total_forceplate_cop_world');
                    this.basic_variable_data.cop_x = total_forceplate_cop_world_data(:, 1);
                    this.basic_variable_directions.cop_x = total_forceplate_cop_world_directions(:, 1);
                    this.time_data.cop_x = this.time_data.total_forceplate_cop_world;
                    success = 1;
                end
                if strcmp(variable_name, 'cop_y_vel')
                    cop_y = this.getBasicVariableData('cop_y');
                    cop_y(cop_y==0) = NaN;
                    time = this.getTimeData('cop_y');
                    filter_order = this.study_settings.get('force_plate_derivative_filter_order');
                    cutoff_frequency = this.study_settings.get('force_plate_derivative_filter_cutoff');
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;
                    cop_y_vel = deriveByTime(nanfiltfilt(b, a, cop_y), 1/sampling_rate);
                    cop_y_vel(isnan(cop_y_vel)) = 0;
                    this.basic_variable_data.cop_y_vel = cop_y_vel;
                    this.basic_variable_directions.cop_y_vel = this.basic_variable_directions.cop_y;
                    this.time_data.cop_y_vel = time;
                    success = 1;
                end
                if strcmp(variable_name, 'cop_y_acc')
                    cop_y_vel = this.getBasicVariableData('cop_y_vel');
                    cop_y_vel(cop_y_vel==0) = NaN;
                    time = this.getTimeData('cop_y_vel');
                    filter_order = this.study_settings.get('force_plate_derivative_filter_order');
                    cutoff_frequency = this.study_settings.get('force_plate_derivative_filter_cutoff');
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;
                    cop_y_acc = deriveByTime(nanfiltfilt(b, a, cop_y_vel), 1/sampling_rate);
                    cop_y_acc(isnan(cop_y_acc)) = 0;
                    this.basic_variable_data.cop_y_acc = cop_y_acc;
                    this.basic_variable_directions.cop_y_acc = this.basic_variable_directions.cop_y;
                    this.time_data.cop_y_acc = time;
                    success = 1;
                end
                if strcmp(variable_name, 'cop_x_vel')
                    cop_x = this.getBasicVariableData('cop_x');
                    cop_x(cop_x==0) = NaN;
                    time = this.getTimeData('cop_x');
                    filter_order = this.study_settings.get('force_plate_derivative_filter_order');
                    cutoff_frequency = this.study_settings.get('force_plate_derivative_filter_cutoff');
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;
                    cop_x_vel = deriveByTime(nanfiltfilt(b, a, cop_x), 1/sampling_rate);
                    cop_x_vel(isnan(cop_x_vel)) = 0;
                    this.basic_variable_data.cop_x_vel = cop_x_vel;
                    this.basic_variable_directions.cop_x_vel = this.basic_variable_directions.cop_x;
                    this.time_data.cop_x_vel = time;
                    success = 1;
                end
                if strcmp(variable_name, 'cop_x_acc')
                    cop_x_vel = this.getBasicVariableData('cop_x_vel');
                    cop_x_vel(cop_x_vel==0) = NaN;
                    time = this.getTimeData('cop_x_vel');
                    filter_order = this.study_settings.get('force_plate_derivative_filter_order');
                    cutoff_frequency = this.study_settings.get('force_plate_derivative_filter_cutoff');
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;
                    cop_x_acc = deriveByTime(nanfiltfilt(b, a, cop_x_vel), 1/sampling_rate);
                    cop_x_acc(isnan(cop_x_acc)) = 0;
                    this.basic_variable_data.cop_x_acc = cop_x_acc;
                    this.basic_variable_directions.cop_x_acc = this.basic_variable_directions.cop_x;
                    this.time_data.cop_x_acc = time;
                    success = 1;
                end
                if strcmp(variable_name, 'fxl')
                    [left_foot_wrench_world_data, left_foot_wrench_world_directions] = this.getBasicVariableData('left_foot_wrench_world');
                    this.basic_variable_data.fxl = left_foot_wrench_world_data(:, 1);
                    this.basic_variable_directions.fxl = left_foot_wrench_world_directions(:, 1);
                    this.time_data.fxl = this.time_data.left_foot_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'fyl')
                    [left_foot_wrench_world_data, left_foot_wrench_world_directions] = this.getBasicVariableData('left_foot_wrench_world');
                    this.basic_variable_data.fyl = left_foot_wrench_world_data(:, 2);
                    this.basic_variable_directions.fyl = left_foot_wrench_world_directions(:, 2);
                    this.time_data.fyl = this.time_data.left_foot_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'fzl')
                    [left_foot_wrench_world_data, left_foot_wrench_world_directions] = this.getBasicVariableData('left_foot_wrench_world');
                    this.basic_variable_data.fzl = left_foot_wrench_world_data(:, 3);
                    this.basic_variable_directions.fzl = left_foot_wrench_world_directions(:, 3);
                    this.time_data.fzl = this.time_data.left_foot_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'mxl')
                    [left_foot_wrench_world_data, left_foot_wrench_world_directions] = this.getBasicVariableData('left_foot_wrench_world');
                    this.basic_variable_data.mxl = left_foot_wrench_world_data(:, 4);
                    this.basic_variable_directions.mxl = left_foot_wrench_world_directions(:, 4);
                    this.time_data.mxl = this.time_data.left_foot_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'myl')
                    [left_foot_wrench_world_data, left_foot_wrench_world_directions] = this.getBasicVariableData('left_foot_wrench_world');
                    this.basic_variable_data.myl = left_foot_wrench_world_data(:, 5);
                    this.basic_variable_directions.myl = left_foot_wrench_world_directions(:, 5);
                    this.time_data.myl = this.time_data.left_foot_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'mzl')
                    [left_foot_wrench_world_data, left_foot_wrench_world_directions] = this.getBasicVariableData('left_foot_wrench_world');
                    this.basic_variable_data.mzl = left_foot_wrench_world_data(:, 6);
                    this.basic_variable_directions.mzl = left_foot_wrench_world_directions(:, 6);
                    this.time_data.mzl = this.time_data.left_foot_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'fxr')
                    [right_foot_wrench_world_data, right_foot_wrench_world_directions] = this.getBasicVariableData('right_foot_wrench_world');
                    this.basic_variable_data.fxr = right_foot_wrench_world_data(:, 1);
                    this.basic_variable_directions.fxr = right_foot_wrench_world_directions(:, 1);
                    this.time_data.fxr = this.time_data.right_foot_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'fyr')
                    [right_foot_wrench_world_data, right_foot_wrench_world_directions] = this.getBasicVariableData('right_foot_wrench_world');
                    this.basic_variable_data.fyr = right_foot_wrench_world_data(:, 2);
                    this.basic_variable_directions.fyr = right_foot_wrench_world_directions(:, 2);
                    this.time_data.fyr = this.time_data.right_foot_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'fzr')
                    [right_foot_wrench_world_data, right_foot_wrench_world_directions] = this.getBasicVariableData('right_foot_wrench_world');
                    this.basic_variable_data.fzr = right_foot_wrench_world_data(:, 3);
                    this.basic_variable_directions.fzr = right_foot_wrench_world_directions(:, 3);
                    this.time_data.fzr = this.time_data.right_foot_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'mxr')
                    [right_foot_wrench_world_data, right_foot_wrench_world_directions] = this.getBasicVariableData('right_foot_wrench_world');
                    this.basic_variable_data.mxr = right_foot_wrench_world_data(:, 4);
                    this.basic_variable_directions.mxr = right_foot_wrench_world_directions(:, 4);
                    this.time_data.mxr = this.time_data.right_foot_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'myr')
                    [right_foot_wrench_world_data, right_foot_wrench_world_directions] = this.getBasicVariableData('right_foot_wrench_world');
                    this.basic_variable_data.myr = right_foot_wrench_world_data(:, 5);
                    this.basic_variable_directions.myr = right_foot_wrench_world_directions(:, 5);
                    this.time_data.myr = this.time_data.right_foot_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'mzr')
                    [right_foot_wrench_world_data, right_foot_wrench_world_directions] = this.getBasicVariableData('right_foot_wrench_world');
                    this.basic_variable_data.mzr = right_foot_wrench_world_data(:, 6);
                    this.basic_variable_directions.mzr = right_foot_wrench_world_directions(:, 6);
                    this.time_data.mzr = this.time_data.right_foot_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'fx')
                    [total_forceplate_wrench_world_data, total_forceplate_wrench_world_directions] = this.getBasicVariableData('total_forceplate_wrench_world');
                    this.basic_variable_data.fx = total_forceplate_wrench_world_data(:, 1);
                    this.basic_variable_directions.fx = total_forceplate_wrench_world_directions(:, 1);
                    this.time_data.fx = this.time_data.total_forceplate_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'fy')
                    [total_forceplate_wrench_world_data, total_forceplate_wrench_world_directions] = this.getBasicVariableData('total_forceplate_wrench_world');
                    this.basic_variable_data.fy = total_forceplate_wrench_world_data(:, 2);
                    this.basic_variable_directions.fy = total_forceplate_wrench_world_directions(:, 2);
                    this.time_data.fy = this.time_data.total_forceplate_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'fz')
                    [total_forceplate_wrench_world_data, total_forceplate_wrench_world_directions] = this.getBasicVariableData('total_forceplate_wrench_world');
                    this.basic_variable_data.fz = total_forceplate_wrench_world_data(:, 3);
                    this.basic_variable_directions.fz = total_forceplate_wrench_world_directions(:, 3);
                    this.time_data.fz = this.time_data.total_forceplate_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'mx')
                    [total_forceplate_wrench_world_data, total_forceplate_wrench_world_directions] = this.getBasicVariableData('total_forceplate_wrench_world');
                    this.basic_variable_data.mx = total_forceplate_wrench_world_data(:, 4);
                    this.basic_variable_directions.mx = total_forceplate_wrench_world_directions(:, 4);
                    this.time_data.mx = this.time_data.total_forceplate_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'my')
                    [total_forceplate_wrench_world_data, total_forceplate_wrench_world_directions] = this.getBasicVariableData('total_forceplate_wrench_world');
                    this.basic_variable_data.my = total_forceplate_wrench_world_data(:, 5);
                    this.basic_variable_directions.my = total_forceplate_wrench_world_directions(:, 5);
                    this.time_data.my = this.time_data.total_forceplate_wrench_world;
                    success = 1;
                end
                if strcmp(variable_name, 'mz')
                    [total_forceplate_wrench_world_data, total_forceplate_wrench_world_directions] = this.getBasicVariableData('total_forceplate_wrench_world');
                    this.basic_variable_data.mz = total_forceplate_wrench_world_data(:, 6);
                    this.basic_variable_directions.mz = total_forceplate_wrench_world_directions(:, 6);
                    this.time_data.mz = this.time_data.total_forceplate_wrench_world;
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
%                     
                    if strcmp(variable_name, 'lheel_from_mpsis_initial_x')
                        lheel_x = this.getTimeNormalizedData('lheel_x', this_stretch_times);
                        mpsis_x = this.getTimeNormalizedData('mpsis_x', this_stretch_times);
                        stretch_data = lheel_x - mpsis_x(1);
                    end
                    if strcmp(variable_name, 'lheel_from_com_initial_x')
                        lheel_x = this.getTimeNormalizedData('lheel_x', this_stretch_times);
                        com_x = this.getTimeNormalizedData('com_x', this_stretch_times);
                        stretch_data = lheel_x - com_x(1);
                    end
                    if strcmp(variable_name, 'rheel_from_mpsis_initial_x')
                        rheel_x = this.getTimeNormalizedData('rheel_x', this_stretch_times);
                        mpsis_x = this.getTimeNormalizedData('mpsis_x', this_stretch_times);
                        stretch_data = rheel_x - mpsis_x(1);
                    end
                    if strcmp(variable_name, 'rheel_from_com_initial_x')
                        rheel_x = this.getTimeNormalizedData('rheel_x', this_stretch_times);
                        com_x = this.getTimeNormalizedData('com_x', this_stretch_times);
                        stretch_data = rheel_x - com_x(1);
                    end
                    if strcmp(variable_name, 'mpsis_from_mpsis_initial_x')
                        mpsis_x = this.getTimeNormalizedData('mpsis_x', this_stretch_times);
                        stretch_data = mpsis_x - mpsis_x(1);
                    end
                    if strcmp(variable_name, 'com_from_com_initial_x')
                       com_x =  this.getTimeNormalizedData('com_x', this_stretch_times);
                       stretch_data = com_x - com_x(1);
                    end
                    if strcmp(variable_name, 'step_length')
                        lheel_y = this.getTimeNormalizedData('lheel_y', this_stretch_times);
                        rheel_y = this.getTimeNormalizedData('rheel_y', this_stretch_times);
                        
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            [~, band_end_indices] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stretch_data(i_band) = lheel_y(band_end_indices) - rheel_y(band_end_indices);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stretch_data(i_band) = rheel_y(band_end_indices) - lheel_y(band_end_indices);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stretch_data(i_band) = NaN;
                            end
                        end
                    end
                    if strcmp(variable_name, 'step_width')
                        lheel_x = this.getTimeNormalizedData('lheel_x', this_stretch_times);
                        rheel_x = this.getTimeNormalizedData('rheel_x', this_stretch_times);
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            [~, band_end_indices] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stretch_data(i_band) = abs(lheel_x(band_end_indices) - rheel_x(band_end_indices));
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stretch_data(i_band) = abs(rheel_x(band_end_indices) - lheel_x(band_end_indices));
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stretch_data(i_band) = NaN;
                            end
                        end
                    end
                    if strcmp(variable_name, 'step_placement_x')
                        lheel_x = this.getTimeNormalizedData('lheel_x', this_stretch_times);
                        rheel_x = this.getTimeNormalizedData('rheel_x', this_stretch_times);
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            [~, band_end_indices] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stretch_data(i_band) = lheel_x(band_end_indices) - rheel_x(band_end_indices);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stretch_data(i_band) = rheel_x(band_end_indices) - lheel_x(band_end_indices);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stretch_data(i_band) = NaN;
                            end
                        end
                    end
                    if strcmp(variable_name, 'step_time')
                        stretch_data = diff(this_stretch_times)';
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
                    if strcmp(variable_name, 'midstance_index')
                        lankle_x = this.getTimeNormalizedData('lankle_y', this_stretch_times);
                        rankle_x = this.getTimeNormalizedData('rankle_y', this_stretch_times);
                        pelvis_y = this.getTimeNormalizedData('pelvis_y', this_stretch_times);
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stance_ankle_this_band = rankle_x(band_start_index : band_end_index);
                                pelvis_this_band = pelvis_y(band_start_index : band_end_index);
                                [~, zero_crossing_index] = min(abs(stance_ankle_this_band - pelvis_this_band));
                                stretch_data(i_band) = band_start_index - 1 + zero_crossing_index;
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stance_ankle_this_band = lankle_x(band_start_index : band_end_index);
                                pelvis_this_band = pelvis_y(band_start_index : band_end_index);
                                [~, zero_crossing_index] = min(abs(stance_ankle_this_band - pelvis_this_band));
                                stretch_data(i_band) = band_start_index - 1 + zero_crossing_index;
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stretch_data(i_band) = NaN;
                            end
                        end
                        
%                         if strcmp(stance_foot_data{i_stretch}, 'STANCE_RIGHT')
%                             stance_ankle_y = this.getTimeNormalizedData('rankle_y', this_stretch_times);
%                         end
%                         if strcmp(stance_foot_data{i_stretch}, 'STANCE_LEFT')
%                             stance_ankle_y = this.getTimeNormalizedData('lankle_y', this_stretch_times);
%                         end
%                         if strcmp(stance_foot_data{i_stretch}, 'STANCE_BOTH')
%                             stance_ankle_y = NaN;
%                         end
%                         stretch_data = find(stance_ankle_y < pelvis_y, 1, 'first');
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
                    if strcmp(variable_name, 'cop_from_mpsis_x')
                        mpsis_x = stretch_variables{strcmp(this.stretch_variable_names, 'mpsis_x')}(:, i_stretch);
                        cop_x = stretch_variables{strcmp(this.stretch_variable_names, 'cop_x')}(:, i_stretch);
                        stretch_data = cop_x - mpsis_x;
                    end
                    if strcmp(variable_name, 'cop_from_mpsis_y')
                        mpsis_y = stretch_variables{strcmp(this.stretch_variable_names, 'mpsis_y')}(:, i_stretch);
                        cop_y = stretch_variables{strcmp(this.stretch_variable_names, 'cop_y')}(:, i_stretch);
                        stretch_data = cop_y - mpsis_y;
                    end
                    if strcmp(variable_name, 'cop_from_com_x')
                        com_x = stretch_variables{strcmp(this.stretch_variable_names, 'com_x')}(:, i_stretch);
                        cop_x = stretch_variables{strcmp(this.stretch_variable_names, 'cop_x')}(:, i_stretch);
                        stretch_data = cop_x - com_x;
                    end
                    if strcmp(variable_name, 'cop_from_com_y')
                        com_y = stretch_variables{strcmp(this.stretch_variable_names, 'com_y')}(:, i_stretch);
                        cop_y = stretch_variables{strcmp(this.stretch_variable_names, 'cop_y')}(:, i_stretch);
                        stretch_data = cop_y - com_y;
                    end
                    if strcmp(variable_name, 'cop_to_com_x')
                        com_x = stretch_variables{strcmp(this.stretch_variable_names, 'com_x')}(:, i_stretch);
                        cop_x = stretch_variables{strcmp(this.stretch_variable_names, 'cop_x')}(:, i_stretch);
                        stretch_data = com_x - cop_x;
                    end
                    if strcmp(variable_name, 'cop_to_com_y')
                        com_y = stretch_variables{strcmp(this.stretch_variable_names, 'com_y')}(:, i_stretch);
                        cop_y = stretch_variables{strcmp(this.stretch_variable_names, 'cop_y')}(:, i_stretch);
                        stretch_data = com_y - cop_y;
                    end
                    if strcmp(variable_name, 'init_heel_to_com_x')
                        com_x = stretch_variables{strcmp(this.stretch_variable_names, 'com_x')}(:, i_stretch);
                        
                        lheel_x = this.getTimeNormalizedData('lheel_x', this_stretch_times);
                        rheel_x = this.getTimeNormalizedData('rheel_x', this_stretch_times);
                        stretch_data = zeros(size(rheel_x)) * NaN;
                        
                        
                        for i_band = 1 : number_of_bands
                            [band_start_index, band_end_index] = getBandIndices(i_band, this.number_of_time_steps_normalized);
                            
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                stretch_data(band_start_index : band_end_index) = com_x(band_start_index : band_end_index) - rheel_x(band_start_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                stretch_data(band_start_index : band_end_index) = com_x(band_start_index : band_end_index) - lheel_x(band_start_index);
                            end
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                % do nothing, leave NaNs
                            end
                        end
                    end
                    
                    
                    if strcmp(variable_name, 'cop_to_com_vel_scaled_x')
                        com_z = stretch_variables{strcmp(this.stretch_variable_names, 'com_z')}(:, i_stretch);
                        com_x_vel = stretch_variables{strcmp(this.stretch_variable_names, 'com_x_vel')}(:, i_stretch);
                        cop_x_vel = stretch_variables{strcmp(this.stretch_variable_names, 'cop_x_vel')}(:, i_stretch);
                        
                        g = 9.81;
                        omega = (g./com_z).^(0.5);
                        
                        stretch_data = (com_x_vel - cop_x_vel) ./ omega;
                        
                        
                        com_x = stretch_variables{strcmp(this.stretch_variable_names, 'com_x')}(:, i_stretch);
                        cop_x = stretch_variables{strcmp(this.stretch_variable_names, 'cop_x')}(:, i_stretch);
                        pos = com_x - cop_x;
                        
                        
                    end
                    if strcmp(variable_name, 'cop_to_com_vel_scaled_y')
                        com_y_vel = stretch_variables{strcmp(this.stretch_variable_names, 'com_y_vel')}(:, i_stretch);
                        cop_y_vel = stretch_variables{strcmp(this.stretch_variable_names, 'cop_y_vel')}(:, i_stretch);
                        stretch_data = com_y_vel - cop_y_vel;
                    end
                    if strcmp(variable_name, 'com_x_vel_scaled')
                        com_z = stretch_variables{strcmp(this.stretch_variable_names, 'com_z')}(:, i_stretch);
                        com_x_vel = stretch_variables{strcmp(this.stretch_variable_names, 'com_x_vel')}(:, i_stretch);
                        
                        g = 9.81;
                        omega = (g./com_z).^(0.5);
                        
                        stretch_data = com_x_vel ./ omega;
                    end
                    
                    if strcmp(variable_name, 'heel_clearance')
                        stretch_data = zeros(number_of_bands, 1);
                        for i_band = 1 : number_of_bands
                            if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_BOTH')
                                stretch_data(i_band) = NaN;
                            else
                                % get relevant data
                                if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_RIGHT')
                                    swing_marker_y_complete = this.getBasicVariableData('lheel_y');
                                    swing_marker_z_complete = this.getBasicVariableData('lheel_z');
                                    variable_time = this.getTimeData('lheel_y');
                                end
                                if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                    swing_marker_y_complete = this.getBasicVariableData('rheel_y');
                                    swing_marker_z_complete = this.getBasicVariableData('rheel_z');
                                    variable_time = this.getTimeData('rheel_y');
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
                                    obstacle_pos_y = this.subject_info.near_distance;
                                end
                                if strcmp(relevant_condition_data{i_stretch}, 'OBS_FAR')
                                    obstacle_pos_y = this.subject_info.far_distance;
                                end
                                [~, crossing_time_step] = min(abs(swing_heel_marker_y - obstacle_pos_y));

                                % extract swing foot heel marker
                                if strcmp(relevant_condition_data{i_stretch}, 'OBS_NO')
                                    stretch_data(i_band) = NaN;
                                else
                                    stretch_data(i_band) = swing_heel_marker_z(crossing_time_step) - this.subject_info.obstacle_height;
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
                                    swing_marker_y_complete = this.getBasicVariableData('ltoes_y');
                                    swing_marker_z_complete = this.getBasicVariableData('ltoes_z');
                                    variable_time = this.getTimeData('ltoes_y');
                                end
                                if strcmp(stance_foot_data{i_stretch, i_band}, 'STANCE_LEFT')
                                    swing_marker_y_complete = this.getBasicVariableData('rtoes_y');
                                    swing_marker_z_complete = this.getBasicVariableData('rtoes_z');
                                    variable_time = this.getTimeData('rtoes_y');
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
                                    obstacle_pos_y = this.subject_info.toe_marker + this.subject_info.near_distance;
                                end
                                if strcmp(relevant_condition_data{i_stretch}, 'OBS_FAR')
                                    obstacle_pos_y = this.subject_info.toe_marker + this.subject_info.far_distance;
                                end
                                [~, crossing_time_step] = min(abs(swing_toes_marker_y - obstacle_pos_y));

                                % extract swing foot toes marker
                                if strcmp(relevant_condition_data{i_stretch}, 'OBS_NO')
                                    stretch_data(i_band) = NaN;
                                else
                                    stretch_data(i_band) = swing_toes_marker_z(crossing_time_step) - this.subject_info.obstacle_height;
                                end
                            end
                        end
                    end
                    if strcmp(variable_name, 'copl_x') || strcmp(variable_name, 'copl_y') || strcmp(variable_name, 'copr_x') || strcmp(variable_name, 'copr_y')
                        % whole stretches might be zero, deal with this in a special way.
                        band_times = this_stretch_times;

                        % extract data
                        try
                            variable_time = this.getTimeData(variable_name);
                            variable_data = this.getBasicVariableData(variable_name);

                            band_time_indices = zeros(size(band_times));
                            for i_band_time = 1 : length(band_times)
                                [~, time_index] = min(abs(variable_time - band_times(i_band_time)));
                                band_time_indices(i_band_time) = time_index;
                            end

                            data_extracted = variable_data(band_time_indices(1) : band_time_indices(end));
                            band_time_indices_local = band_time_indices - band_time_indices(1) + 1;
                            if any(isnan(data_extracted))
                                exception = MException('CoBaL:NaN', 'Data contains NaN values.');
                                throw(exception)
                            end
                        catch error
                            disp(['Error while processing variable ''' variable_name ''''])
                            throw(error)
                        end
                        
                        % check whether last band might be zeros only
                        last_band_zeros_only = 0;
                        data_last_band = data_extracted(band_time_indices_local(i_band) : band_time_indices_local(i_band+1));
                        if ~any(data_last_band(2:end)~=0)
                            last_band_zeros_only = 1;
                        end
                        
                        if last_band_zeros_only
                            data_normalized_first_bands = this.getTimeNormalizedData(variable_name, band_times(1:end-1));
                            data_normalized_last_band = zeros(this.number_of_time_steps_normalized-1, 1);
                            stretch_data = [data_normalized_first_bands; data_normalized_last_band];
                        else
                            stretch_data = this.getTimeNormalizedData(variable_name, band_times);
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
                lheel_x_directions = this.basic_variable_directions.lheel_x;
                mpsis_x_directions = this.basic_variable_directions.mpsis_x;
                if ~strcmp(lheel_x_directions{1}, mpsis_x_directions{1})
                    error('lheel_x and mpsis_x directions are different from each other')
                end
                if ~strcmp(lheel_x_directions{2}, mpsis_x_directions{2})
                    error('lheel_x and mpsis_x directions are different from each other')
                end
                stretch_directions_new = lheel_x_directions;
            end
            if strcmp(variable_name, 'lheel_from_com_initial_x')
                lheel_x_directions = this.basic_variable_directions.lheel_x;
                com_x_directions = this.basic_variable_directions.com_x;
                if ~strcmp(lheel_x_directions{1}, com_x_directions{1})
                    error('lheel_x and com_x directions are different from each other')
                end
                if ~strcmp(lheel_x_directions{2}, com_x_directions{2})
                    error('lheel_x and com_x directions are different from each other')
                end
                stretch_directions_new = lheel_x_directions;
            end
            if strcmp(variable_name, 'rheel_from_mpsis_initial_x')
                rheel_x_directions = this.basic_variable_directions.rheel_x;
                mpsis_x_directions = this.basic_variable_directions.mpsis_x;
                if ~strcmp(rheel_x_directions{1}, mpsis_x_directions{1})
                    error('rheel_x and mpsis_x directions are different from each other')
                end
                if ~strcmp(rheel_x_directions{2}, mpsis_x_directions{2})
                    error('rheel_x and mpsis_x directions are different from each other')
                end
                stretch_directions_new = rheel_x_directions;
            end
            if strcmp(variable_name, 'rheel_from_com_initial_x')
                rheel_x_directions = this.basic_variable_directions.rheel_x;
                com_x_directions = this.basic_variable_directions.com_x;
                if ~strcmp(rheel_x_directions{1}, com_x_directions{1})
                    error('rheel_x and com_x directions are different from each other')
                end
                if ~strcmp(rheel_x_directions{2}, com_x_directions{2})
                    error('rheel_x and com_x directions are different from each other')
                end
                stretch_directions_new = rheel_x_directions;
            end
            if strcmp(variable_name, 'mpsis_from_mpsis_initial_x')
                mpsis_x_directions = this.basic_variable_directions.mpsis_x;
                stretch_directions_new = mpsis_x_directions;
            end
            if strcmp(variable_name, 'com_from_com_initial_x')
                com_x_directions = this.basic_variable_directions.com_x;
                stretch_directions_new = com_x_directions;
            end
            if strcmp(variable_name, 'step_length')
                lheel_y_directions = this.basic_variable_directions.lheel_y;
                rheel_y_directions = this.basic_variable_directions.rheel_y;
                if ~strcmp(lheel_y_directions{1}, rheel_y_directions{1})
                    error('lheel_y and rheel_y directions are different from each other')
                end
                if ~strcmp(lheel_y_directions{2}, rheel_y_directions{2})
                    error('lheel_y and rheel_y directions are different from each other')
                end
                stretch_directions_new = lheel_y_directions;
            end
            if strcmp(variable_name, 'step_width')
                lheel_x_directions = this.basic_variable_directions.lheel_x;
                rheel_x_directions = this.basic_variable_directions.rheel_x;
                if ~strcmp(lheel_x_directions{1}, rheel_x_directions{1})
                    error('lheel_x and rheel_x directions are different from each other')
                end
                if ~strcmp(lheel_x_directions{2}, rheel_x_directions{2})
                    error('lheel_x and rheel_x directions are different from each other')
                end
                stretch_directions_new = lheel_x_directions;
            end
            if strcmp(variable_name, 'step_placement_x')
                lheel_x_directions = this.basic_variable_directions.lheel_x;
                rheel_x_directions = this.basic_variable_directions.rheel_x;
                if ~strcmp(lheel_x_directions{1}, rheel_x_directions{1})
                    error('lheel_x and rheel_x directions are different from each other')
                end
                if ~strcmp(lheel_x_directions{2}, rheel_x_directions{2})
                    error('lheel_x and rheel_x directions are different from each other')
                end
                stretch_directions_new = lheel_x_directions;
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
            if strcmp(variable_name, 'cop_from_mpsis_x')
                mpsis_x_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'mpsis_x'), :);
                cop_x_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'cop_x'), :);
                if ~strcmp(mpsis_x_directions{1}, cop_x_directions{1})
                    error('mpsis_x and cop_x directions are different from each other')
                end
                if ~strcmp(mpsis_x_directions{2}, cop_x_directions{2})
                    error('mpsis_x and cop_x directions are different from each other')
                end
                stretch_directions_new = mpsis_x_directions;
            end
            if strcmp(variable_name, 'cop_from_mpsis_y')
                mpsis_y_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'mpsis_y'), :);
                cop_y_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'cop_y'), :);
                if ~strcmp(mpsis_y_directions{1}, cop_y_directions{1})
                    error('mpsis_y and cop_y directions are different from each other')
                end
                if ~strcmp(mpsis_y_directions{2}, cop_y_directions{2})
                    error('mpsis_y and cop_y directions are different from each other')
                end
                stretch_directions_new = mpsis_y_directions;
            end
            if strcmp(variable_name, 'cop_from_com_x')
                com_x_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'com_x'), :);
                cop_x_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'cop_x'), :);
                if ~strcmp(com_x_directions{1}, cop_x_directions{1})
                    error('com_x and cop_x directions are different from each other')
                end
                if ~strcmp(com_x_directions{2}, cop_x_directions{2})
                    error('com_x and cop_x directions are different from each other')
                end
                stretch_directions_new = com_x_directions;
            end
            if strcmp(variable_name, 'cop_from_com_y')
                com_y_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'com_y'), :);
                cop_y_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'cop_y'), :);
                if ~strcmp(com_y_directions{1}, cop_y_directions{1})
                    error('com_y and cop_y directions are different from each other')
                end
                if ~strcmp(com_y_directions{2}, cop_y_directions{2})
                    error('com_y and cop_y directions are different from each other')
                end
                stretch_directions_new = com_y_directions;
            end
            if strcmp(variable_name, 'cop_to_com_x')
                com_x_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'com_x'), :);
                cop_x_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'cop_x'), :);
                if ~strcmp(com_x_directions{1}, cop_x_directions{1})
                    error('com_x and cop_x directions are different from each other')
                end
                if ~strcmp(com_x_directions{2}, cop_x_directions{2})
                    error('com_x and cop_x directions are different from each other')
                end
                stretch_directions_new = com_x_directions;
            end
            if strcmp(variable_name, 'cop_to_com_y')
                com_y_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'com_y'), :);
                cop_y_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'cop_y'), :);
                if ~strcmp(com_y_directions{1}, cop_y_directions{1})
                    error('com_y and cop_y directions are different from each other')
                end
                if ~strcmp(com_y_directions{2}, cop_y_directions{2})
                    error('com_y and cop_y directions are different from each other')
                end
                stretch_directions_new = com_y_directions;
            end
            if strcmp(variable_name, 'init_heel_to_com_x')
                com_x_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'com_x'), :);
                stretch_directions_new = com_x_directions;
            end
            if strcmp(variable_name, 'com_x_vel_scaled')
                com_x_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'com_x'), :);
                stretch_directions_new = com_x_directions;
            end
            if strcmp(variable_name, 'cop_to_com_vel_scaled_x')
                com_x_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'com_x'), :);
                cop_x_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'cop_x'), :);
                if ~strcmp(com_x_directions{1}, cop_x_directions{1})
                    error('com_x and cop_x directions are different from each other')
                end
                if ~strcmp(com_x_directions{2}, cop_x_directions{2})
                    error('com_x and cop_x directions are different from each other')
                end
                stretch_directions_new = com_x_directions;
            end
            if strcmp(variable_name, 'cop_to_com_vel_scaled_y')
                com_y_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'com_y'), :);
                cop_y_directions = this.stretch_variable_directions(strcmp(this.stretch_variable_names, 'cop_y'), :);
                if ~strcmp(com_y_directions{1}, cop_y_directions{1})
                    error('com_y and cop_y directions are different from each other')
                end
                if ~strcmp(com_y_directions{2}, cop_y_directions{2})
                    error('com_y and cop_y directions are different from each other')
                end
                stretch_directions_new = com_y_directions;
            end
            if strcmp(variable_name, 'heel_clearance')
                lheel_z_directions = this.basic_variable_directions.lheel_z;
                rheel_z_directions = this.basic_variable_directions.rheel_z;
                if ~strcmp(lheel_z_directions{1}, rheel_z_directions{1})
                    error('lheel_z and rheel_z directions are different from each other')
                end
                if ~strcmp(lheel_z_directions{2}, rheel_z_directions{2})
                    error('lheel_z and rheel_z directions are different from each other')
                end
                stretch_directions_new = lheel_z_directions;
            end
            if strcmp(variable_name, 'toes_clearance')
                ltoes_z_directions = this.basic_variable_directions.ltoes_z;
                rtoes_z_directions = this.basic_variable_directions.rtoes_z;
                if ~strcmp(ltoes_z_directions{1}, rtoes_z_directions{1})
                    error('ltoes_z and rtoes_z directions are different from each other')
                end
                if ~strcmp(ltoes_z_directions{2}, rtoes_z_directions{2})
                    error('ltoes_z and rtoes_z directions are different from each other')
                end
                stretch_directions_new = ltoes_z_directions;
            end
            if strcmp(variable_name, 'copl_x') || strcmp(variable_name, 'copl_y') || strcmp(variable_name, 'copr_x') || strcmp(variable_name, 'copr_y')
                stretch_directions_new = this.basic_variable_directions.(variable_name);
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





