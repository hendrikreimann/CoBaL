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
        date = [];
        subject_id = [];
        variables_to_analyze = {};
        basic_variable_names = {};
        stretch_variable_names = {};
        
        % these are private and should only be accessed using get functions
        number_of_time_steps_normalized;
        basic_variable_data;
        basic_variable_labels;
        stretch_variable_data;
        time_data;
        
        subject_settings;
        study_settings;
        emg_normalization_values;
        emg_normalization_labels;
    end
    
    methods
        % constructor
%         function this = WalkingDataCustodian(date, subject_id, variables_to_analyze)
        function this = WalkingDataCustodian(variables_to_analyze)
            % load this information from the subjects.mat and studySettings.txt files
            load('subjectInfo.mat', 'date', 'subject_id');
            this.subject_settings = loadSettingsFile('subjectSettings.txt');
            % load settings
            study_settings_file = '';
            if exist(['..' filesep 'studySettings.txt'], 'file')
                study_settings_file = ['..' filesep 'studySettings.txt'];
            end    
            if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
                study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
            end
            this.study_settings = loadSettingsFile(study_settings_file);
            emg_normalization_file_name = ['analysis' filesep makeFileName(date, subject_id, 'emgNormalization.mat')];
            if exist(emg_normalization_file_name, 'file')
                emg_normalization_data = load(emg_normalization_file_name);
                this.emg_normalization_values = emg_normalization_data.emg_normalization_values;
                this.emg_normalization_labels = emg_normalization_data.emg_variable_names;
            end
            if nargin < 1
                variables_to_analyze = this.study_settings.variables_to_analyze;
            end
            
            this.date = date;
            this.subject_id = subject_id;
            this.variables_to_analyze = variables_to_analyze;
            this.number_of_time_steps_normalized = this.study_settings.number_of_time_steps_normalized;
            
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
            end
        end
        function determineVariables(this)
            % for each possible variable to analyze, list the basic and required variables required to calculate it
            
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
            if this.isVariableToAnalyze('step_time')
                this.addStretchVariable('step_time')
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
            if this.isVariableToAnalyze('cop_from_com_x')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_x')
                this.addStretchVariable('cop_x')
                this.addBasicVariable('com_trajectories')
                this.addBasicVariable('com_x')
                this.addStretchVariable('com_x')
                this.addStretchVariable('cop_from_com_x')
            end
            if this.isVariableToAnalyze('cop_from_com_y')
                this.addBasicVariable('total_forceplate_cop_world')
                this.addBasicVariable('cop_y')
                this.addStretchVariable('cop_y')
                this.addBasicVariable('com_trajectories')
                this.addBasicVariable('com_y')
                this.addStretchVariable('com_y')
                this.addStretchVariable('cop_from_com_y')
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
                this.addBasicVariable('com_trajectories')
                this.addBasicVariable('com_x')
                this.addStretchVariable('com_x')
            end
            if this.isVariableToAnalyze('com_y')
                this.addBasicVariable('com_trajectories')
                this.addBasicVariable('com_y')
                this.addStretchVariable('com_y')
            end
            if this.isVariableToAnalyze('com_z')
                this.addBasicVariable('com_trajectories')
                this.addBasicVariable('com_z')
                this.addStretchVariable('com_z')
            end
            if this.isVariableToAnalyze('com_x_vel')
                this.addBasicVariable('com_trajectories')
                this.addBasicVariable('com_x')
                this.addBasicVariable('com_x_vel')
                this.addStretchVariable('com_x_vel')
            end
            if this.isVariableToAnalyze('com_y_vel')
                this.addBasicVariable('com_trajectories')
                this.addBasicVariable('com_y')
                this.addBasicVariable('com_y_vel')
                this.addStretchVariable('com_y_vel')
            end
            if this.isVariableToAnalyze('com_z_vel')
                this.addBasicVariable('com_trajectories')
                this.addBasicVariable('com_z')
                this.addBasicVariable('com_z_vel')
                this.addStretchVariable('com_z_vel')
            end
            % joint angles
            if this.isVariableToAnalyze('lumbar_roll_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('lumbar_roll_angle')
                this.addStretchVariable('lumbar_roll_angle')
            end
            if this.isVariableToAnalyze('lumbar_pitch_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('lumbar_pitch_angle')
                this.addStretchVariable('lumbar_pitch_angle')
            end
            if this.isVariableToAnalyze('lumbar_yaw_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('lumbar_yaw_angle')
                this.addStretchVariable('lumbar_yaw_angle')
            end
            if this.isVariableToAnalyze('cevical_roll_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('cevical_roll_angle')
                this.addStretchVariable('cevical_roll_angle')
            end
            if this.isVariableToAnalyze('cevical_pitch_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('cevical_pitch_angle')
                this.addStretchVariable('cevical_pitch_angle')
            end
            if this.isVariableToAnalyze('cevical_yaw_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('cevical_yaw_angle')
                this.addStretchVariable('cevical_yaw_angle')
            end
            if this.isVariableToAnalyze('left_hip_abduction_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('left_hip_abduction_angle')
                this.addStretchVariable('left_hip_abduction_angle')
            end
            if this.isVariableToAnalyze('left_hip_flexion_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('left_hip_flexion_angle')
                this.addStretchVariable('left_hip_flexion_angle')
            end
            if this.isVariableToAnalyze('left_hip_introtation_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('left_hip_introtation_angle')
                this.addStretchVariable('left_hip_introtation_angle')
            end
            if this.isVariableToAnalyze('left_knee_flexion_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('left_knee_flexion_angle')
                this.addStretchVariable('left_knee_flexion_angle')
            end
            if this.isVariableToAnalyze('left_knee_extrotation_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('left_knee_extrotation_angle')
                this.addStretchVariable('left_knee_extrotation_angle')
            end
            if this.isVariableToAnalyze('left_ankle_inversion_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('left_ankle_inversion_angle')
                this.addStretchVariable('left_ankle_inversion_angle')
            end
            if this.isVariableToAnalyze('left_ankle_dorsiflexion_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('left_ankle_dorsiflexion_angle')
                this.addStretchVariable('left_ankle_dorsiflexion_angle')
            end
            if this.isVariableToAnalyze('right_hip_abduction_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('right_hip_abduction_angle')
                this.addStretchVariable('right_hip_abduction_angle')
            end
            if this.isVariableToAnalyze('right_hip_flexion_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('right_hip_flexion_angle')
                this.addStretchVariable('right_hip_flexion_angle')
            end
            if this.isVariableToAnalyze('right_hip_introtation_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('right_hip_introtation_angle')
                this.addStretchVariable('right_hip_introtation_angle')
            end
            if this.isVariableToAnalyze('right_knee_flexion_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('right_knee_flexion_angle')
                this.addStretchVariable('right_knee_flexion_angle')
            end
            if this.isVariableToAnalyze('right_knee_extrotation_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('right_knee_extrotation_angle')
                this.addStretchVariable('right_knee_extrotation_angle')
            end
            if this.isVariableToAnalyze('right_ankle_inversion_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('right_ankle_inversion_angle')
                this.addStretchVariable('right_ankle_inversion_angle')
            end
            if this.isVariableToAnalyze('right_ankle_dorsiflexion_angle')
                this.addBasicVariable('joint_angle_trajectories')
                this.addBasicVariable('right_ankle_dorsiflexion_angle')
                this.addStretchVariable('right_ankle_dorsiflexion_angle')
            end
            % joint torques
            if this.isVariableToAnalyze('lumbar_roll_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('lumbar_roll_torque')
                this.addStretchVariable('lumbar_roll_torque')
            end
            if this.isVariableToAnalyze('lumbar_pitch_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('lumbar_pitch_torque')
                this.addStretchVariable('lumbar_pitch_torque')
            end
            if this.isVariableToAnalyze('lumbar_yaw_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('lumbar_yaw_torque')
                this.addStretchVariable('lumbar_yaw_torque')
            end
            if this.isVariableToAnalyze('cevical_roll_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('cevical_roll_torque')
                this.addStretchVariable('cevical_roll_torque')
            end
            if this.isVariableToAnalyze('cevical_pitch_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('cevical_pitch_torque')
                this.addStretchVariable('cevical_pitch_torque')
            end
            if this.isVariableToAnalyze('cevical_yaw_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('cevical_yaw_torque')
                this.addStretchVariable('cevical_yaw_torque')
            end
            if this.isVariableToAnalyze('left_hip_abduction_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('left_hip_abduction_torque')
                this.addStretchVariable('left_hip_abduction_torque')
            end
            if this.isVariableToAnalyze('left_hip_flexion_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('left_hip_flexion_torque')
                this.addStretchVariable('left_hip_flexion_torque')
            end
            if this.isVariableToAnalyze('left_hip_introtation_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('left_hip_introtation_torque')
                this.addStretchVariable('left_hip_introtation_torque')
            end
            if this.isVariableToAnalyze('left_knee_flexion_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('left_knee_flexion_torque')
                this.addStretchVariable('left_knee_flexion_torque')
            end
            if this.isVariableToAnalyze('left_knee_extrotation_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('left_knee_extrotation_torque')
                this.addStretchVariable('left_knee_extrotation_torque')
            end
            if this.isVariableToAnalyze('left_ankle_inversion_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('left_ankle_inversion_torque')
                this.addStretchVariable('left_ankle_inversion_torque')
            end
            if this.isVariableToAnalyze('left_ankle_dorsiflexion_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('left_ankle_dorsiflexion_torque')
                this.addStretchVariable('left_ankle_dorsiflexion_torque')
            end
            if this.isVariableToAnalyze('right_hip_abduction_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('right_hip_abduction_torque')
                this.addStretchVariable('right_hip_abduction_torque')
            end
            if this.isVariableToAnalyze('right_hip_flexion_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('right_hip_flexion_torque')
                this.addStretchVariable('right_hip_flexion_torque')
            end
            if this.isVariableToAnalyze('right_hip_introtation_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('right_hip_introtation_torque')
                this.addStretchVariable('right_hip_introtation_torque')
            end
            if this.isVariableToAnalyze('right_knee_flexion_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('right_knee_flexion_torque')
                this.addStretchVariable('right_knee_flexion_torque')
            end
            if this.isVariableToAnalyze('right_knee_extrotation_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('right_knee_extrotation_torque')
                this.addStretchVariable('right_knee_extrotation_torque')
            end
            if this.isVariableToAnalyze('right_ankle_inversion_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('right_ankle_inversion_torque')
                this.addStretchVariable('right_ankle_inversion_torque')
            end
            if this.isVariableToAnalyze('right_ankle_dorsiflexion_torque')
                this.addBasicVariable('joint_torque_trajectories_belt')
                this.addBasicVariable('right_ankle_dorsiflexion_torque')
                this.addStretchVariable('right_ankle_dorsiflexion_torque')
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
            % emg
            if this.isVariableToAnalyze('left_glut_med')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_glut_med')
                this.addStretchVariable('left_glut_med')
            end
            if this.isVariableToAnalyze('left_delt_ant')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_delt_ant')
                this.addStretchVariable('left_delt_ant')
            end
            if this.isVariableToAnalyze('left_tibi_ant')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_tibi_ant')
                this.addStretchVariable('left_tibi_ant')
            end
            if this.isVariableToAnalyze('left_gastroc_med')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_gastroc_med')
                this.addStretchVariable('left_gastroc_med')
            end
            if this.isVariableToAnalyze('left_pero_lng')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_pero_lng')
                this.addStretchVariable('left_pero_lng')
            end
            if this.isVariableToAnalyze('left_tfl')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_tfl')
                this.addStretchVariable('left_tfl')
            end
            if this.isVariableToAnalyze('right_glut_med')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_glut_med')
                this.addStretchVariable('right_glut_med')
            end
            if this.isVariableToAnalyze('right_delt_ant')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_delt_ant')
                this.addStretchVariable('right_delt_ant')
            end
            if this.isVariableToAnalyze('right_tibi_ant')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_tibi_ant')
                this.addStretchVariable('right_tibi_ant')
            end
            if this.isVariableToAnalyze('right_gastroc_med')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_gastroc_med')
                this.addStretchVariable('right_gastroc_med')
            end
            if this.isVariableToAnalyze('right_pero_lng')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_pero_lng')
                this.addStretchVariable('right_pero_lng')
            end
            if this.isVariableToAnalyze('right_tfl')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_tfl')
                this.addStretchVariable('right_tfl')
            end
            if this.isVariableToAnalyze('left_glut_med_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_glut_med')
                this.addStretchVariable('left_glut_med')
                this.addStretchVariable('left_glut_med_rescaled')
            end
            if this.isVariableToAnalyze('left_delt_ant_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_delt_ant')
                this.addStretchVariable('left_delt_ant')
                this.addStretchVariable('left_delt_ant_rescaled')
            end
            if this.isVariableToAnalyze('left_tibi_ant_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_tibi_ant')
                this.addStretchVariable('left_tibi_ant')
                this.addStretchVariable('left_tibi_ant_rescaled')
            end
            if this.isVariableToAnalyze('left_gastroc_med_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_gastroc_med')
                this.addStretchVariable('left_gastroc_med')
                this.addStretchVariable('left_gastroc_med_rescaled')
            end
            if this.isVariableToAnalyze('left_pero_lng_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_pero_lng')
                this.addStretchVariable('left_pero_lng')
                this.addStretchVariable('left_pero_lng_rescaled')
            end
            if this.isVariableToAnalyze('left_tfl_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('left_tfl')
                this.addStretchVariable('left_tfl')
                this.addStretchVariable('left_tfl_rescaled')
            end
            if this.isVariableToAnalyze('right_glut_med_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_glut_med')
                this.addStretchVariable('right_glut_med')
                this.addStretchVariable('right_glut_med_rescaled')
            end
            if this.isVariableToAnalyze('right_delt_ant_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_delt_ant')
                this.addStretchVariable('right_delt_ant')
                this.addStretchVariable('right_delt_ant_rescaled')
            end
            if this.isVariableToAnalyze('right_tibi_ant_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_tibi_ant')
                this.addStretchVariable('right_tibi_ant')
                this.addStretchVariable('right_tibi_ant_rescaled')
            end
            if this.isVariableToAnalyze('right_gastroc_med_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_gastroc_med')
                this.addStretchVariable('right_gastroc_med')
                this.addStretchVariable('right_gastroc_med_rescaled')
            end
            if this.isVariableToAnalyze('right_pero_lng_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_pero_lng')
                this.addStretchVariable('right_pero_lng')
                this.addStretchVariable('right_pero_lng_rescaled')
            end
            if this.isVariableToAnalyze('right_tfl_rescaled')
                this.addBasicVariable('emg_trajectories')
                this.addBasicVariable('right_tfl')
                this.addStretchVariable('right_tfl')
                this.addStretchVariable('right_tfl_rescaled')
            end
            
            
  
        end
        
        % interface
        function result = isVariableToAnalyze(this, variable_name)
            result = any(strcmp(this.variables_to_analyze, variable_name));
        end
        function result = isBasicVariable(this, variable_name)
            result = any(strcmp(this.basic_variable_names, variable_name));
        end
        function time_data = getTimeData(this, variable_name) %#ok<STOUT,INUSL>
            eval(['time_data = this.time_data.' variable_name ';']);
        end
        function variable_data = getBasicVariableData(this, variable_name) %#ok<STOUT,INUSL>
            eval(['variable_data = this.basic_variable_data.' variable_name ';']);
        end
        
        function prepareBasicVariables(this, condition, trial, variables_to_prepare)
            if nargin < 4
                variables_to_prepare = this.basic_variable_names;
            end
            
            % clear out old data
            this.basic_variable_data = struct;
            this.basic_variable_labels = struct;
            this.stretch_variable_data = struct;
            this.time_data = struct;
            
            % prepare the data by loading all the basic variables from disk and calculating the required variables
            load(['analysis' filesep makeFileName(this.date, this.subject_id, condition, trial, 'availableVariables')]);
            
            % load basic variables
            for i_variable = 1 : length(variables_to_prepare)
                variable_name = variables_to_prepare{i_variable};
                
                % try loading
                [data, time, sampling_rate, labels, success] = loadData(this.date, this.subject_id, condition, trial, variable_name, 'optional'); %#ok<ASGLU>
                
                % store
                if success
                    eval(['this.basic_variable_data.' variable_name ' = data;']);
                    eval(['this.time_data.' variable_name ' = time;']);
                    eval(['this.basic_variable_labels.' variable_name ' = labels;']);
                end
                
                % calculate variables that can't be loaded
                
                % kinematics
                if strcmp(variable_name, 'lheel_x')
                    LHEE_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE');
                    this.basic_variable_data.lheel_x = LHEE_trajectory(:, 1);
                    this.time_data.lheel_x = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'rheel_x')
                    RHEE_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE');
                    this.basic_variable_data.rheel_x = RHEE_trajectory(:, 1);
                    this.time_data.rheel_x = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'lheel_y')
                    LHEE_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE');
                    this.basic_variable_data.lheel_y = LHEE_trajectory(:, 2);
                    this.time_data.lheel_y = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'rheel_y')
                    RHEE_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE');
                    this.basic_variable_data.rheel_y = RHEE_trajectory(:, 2);
                    this.time_data.rheel_y = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'mpsis_x')
                    LPSI_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI');
                    RPSI_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI');
                    MPSI_trajectory = (LPSI_trajectory + RPSI_trajectory) * 0.5;
                    this.basic_variable_data.mpsis_x = MPSI_trajectory(:, 1);
                    this.time_data.mpsis_x = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'mpsis_y')
                    LPSI_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI');
                    RPSI_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI');
                    MPSI_trajectory = (LPSI_trajectory + RPSI_trajectory) * 0.5;
                    this.basic_variable_data.mpsis_y = MPSI_trajectory(:, 2);
                    this.time_data.mpsis_y = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'c7_x')
                    c7_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    this.basic_variable_data.c7_x = c7_trajectory(:, 2);
                    this.time_data.c7_x = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'head_angle_ap')
                    c7_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    LBHD_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LBHD');
                    RBHD_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RBHD');
                    MBHD_trajectory = (LBHD_trajectory + RBHD_trajectory) * 0.5;
                    head_vector_y = MBHD_trajectory(:, 2) - c7_trajectory(:, 2);
                    head_vector_z = MBHD_trajectory(:, 3) - c7_trajectory(:, 3);
                    this.basic_variable_data.head_angle_ap = rad2deg(atan2(head_vector_y, head_vector_z));
                    this.time_data.head_angle_ap = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'head_angle_ml')
                    c7_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    LBHD_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LBHD');
                    RBHD_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RBHD');
                    MBHD_trajectory = (LBHD_trajectory + RBHD_trajectory) * 0.5;
                    head_vector_x = MBHD_trajectory(:, 1) - c7_trajectory(:, 1);
                    head_vector_z = MBHD_trajectory(:, 3) - c7_trajectory(:, 3);
                    this.basic_variable_data.head_angle_ml = rad2deg(atan2(head_vector_x, head_vector_z));
                    this.time_data.head_angle_ml = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'trunk_angle_ap')
                    c7_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    LPSI_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI');
                    RPSI_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI');
                    MPSI_trajectory = (LPSI_trajectory + RPSI_trajectory) * 0.5;
                    trunk_vector_y = c7_trajectory(:, 2) - MPSI_trajectory(:, 2);
                    trunk_vector_z = c7_trajectory(:, 3) - MPSI_trajectory(:, 3);
                    this.basic_variable_data.trunk_angle_ap = rad2deg(atan2(trunk_vector_y, trunk_vector_z));
                    this.time_data.trunk_angle_ap = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'trunk_angle_ml')
                    c7_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'C7');
                    LPSI_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LPSI');
                    RPSI_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RPSI');
                    MPSI_trajectory = (LPSI_trajectory + RPSI_trajectory) * 0.5;
                    trunk_vector_x = c7_trajectory(:, 1) - MPSI_trajectory(:, 1);
                    trunk_vector_z = c7_trajectory(:, 3) - MPSI_trajectory(:, 3);
                    this.basic_variable_data.trunk_angle_ml = rad2deg(atan2(trunk_vector_x, trunk_vector_z));
                    this.time_data.trunk_angle_ml = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'pelvis_angle_ml')
                    left_hip_cor_trajectory = extractMarkerTrajectories(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'LHIPCOR');
                    right_hip_cor_trajectory = extractMarkerTrajectories(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'RHIPCOR');
                    pelvis_vector_x = right_hip_cor_trajectory(:, 1) - left_hip_cor_trajectory(:, 1);
                    pelvis_vector_z = right_hip_cor_trajectory(:, 3) - left_hip_cor_trajectory(:, 3);
                    this.basic_variable_data.pelvis_angle_ml = -rad2deg(atan2(pelvis_vector_z, pelvis_vector_x));
                    this.time_data.pelvis_angle_ml = this.time_data.joint_center_trajectories;
                end
                if strcmp(variable_name, 'left_leg_angle_ml')
                    left_hip_cor_trajectory = extractMarkerTrajectories(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'LHIPCOR');
                    left_ankle_cor_trajectory = extractMarkerTrajectories(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'LANKLECOR');
                    left_leg_vector_x = left_hip_cor_trajectory(:, 1) - left_ankle_cor_trajectory(:, 1);
                    left_leg_vector_z = left_hip_cor_trajectory(:, 3) - left_ankle_cor_trajectory(:, 3);
                    this.basic_variable_data.left_leg_angle_ml = rad2deg(atan2(left_leg_vector_x, left_leg_vector_z));
                    this.time_data.left_leg_angle_ml = this.time_data.joint_center_trajectories;
                end
                if strcmp(variable_name, 'right_leg_angle_ml')
                    right_hip_cor_trajectory = extractMarkerTrajectories(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'RHIPCOR');
                    right_ankle_cor_trajectory = extractMarkerTrajectories(this.basic_variable_data.joint_center_trajectories, this.basic_variable_labels.joint_center_trajectories, 'RANKLECOR');
                    right_leg_vector_x = right_hip_cor_trajectory(:, 1) - right_ankle_cor_trajectory(:, 1);
                    right_leg_vector_z = right_hip_cor_trajectory(:, 3) - right_ankle_cor_trajectory(:, 3);
                    this.basic_variable_data.right_leg_angle_ml = rad2deg(atan2(right_leg_vector_x, right_leg_vector_z));
                    this.time_data.right_leg_angle_ml = this.time_data.joint_center_trajectories;
                end
                if strcmp(variable_name, 'left_arm_right_leg_relative_phase')
                    left_arm_phase = this.getBasicVariableData('left_arm_phase');
                    right_leg_phase = this.getBasicVariableData('right_leg_phase');
                    
                    left_arm_right_leg_relative_phase = normalizeAngle(left_arm_phase - right_leg_phase);
                    
                    % store
                    this.basic_variable_data.left_arm_right_leg_relative_phase = left_arm_right_leg_relative_phase;
                    this.time_data.left_arm_right_leg_relative_phase = this.time_data.marker_trajectories;
                end             
                if strcmp(variable_name, 'right_arm_left_leg_relative_phase')
                    right_arm_phase = this.getBasicVariableData('right_arm_phase');
                    left_leg_phase = this.getBasicVariableData('left_leg_phase');
                    
                    right_arm_left_leg_relative_phase = normalizeAngle(right_arm_phase - left_leg_phase);
                    
                    % store
                    this.basic_variable_data.right_arm_left_leg_relative_phase = right_arm_left_leg_relative_phase;
                    this.time_data.right_arm_left_leg_relative_phase = this.time_data.marker_trajectories;
                end
                if strcmp(variable_name, 'com_x')
                    com_trajectory = extractMarkerTrajectories(this.basic_variable_data.com_trajectories, this.basic_variable_labels.com_trajectories, 'BODYCOM');
                    this.basic_variable_data.com_x = com_trajectory(:, 1);
                    this.time_data.com_x = this.time_data.com_trajectories;
                end
                if strcmp(variable_name, 'com_y')
                    com_trajectory = extractMarkerTrajectories(this.basic_variable_data.com_trajectories, this.basic_variable_labels.com_trajectories, 'BODYCOM');
                    this.basic_variable_data.com_y = com_trajectory(:, 2);
                    this.time_data.com_y = this.time_data.com_trajectories;
                end
                if strcmp(variable_name, 'com_z')
                    com_trajectory = extractMarkerTrajectories(this.basic_variable_data.com_trajectories, this.basic_variable_labels.com_trajectories, 'BODYCOM');
                    this.basic_variable_data.com_z = com_trajectory(:, 3);
                    this.time_data.com_z = this.time_data.com_trajectories;
                end
                if strcmp(variable_name, 'com_x_vel')
                    com_x = this.getBasicVariableData('com_x');
                    time = this.getTimeData('com_x');
                    sampling_rate = 1/median(diff(time));
                    com_x_vel = deriveByTime(nanfiltfilt(b, a, com_x), 1/sampling_rate);
                    this.basic_variable_data.com_x_vel = com_x_vel;
                    this.time_data.com_x_vel = time;
                end
                if strcmp(variable_name, 'com_y_vel')
                    com_y = this.getBasicVariableData('com_y');
                    time = this.getTimeData('com_y');
                    sampling_rate = 1/median(diff(time));
                    com_y_vel = deriveByTime(nanfiltfilt(b, a, com_y), 1/sampling_rate);
                    this.basic_variable_data.com_y_vel = com_y_vel;
                    this.time_data.com_y_vel = time;
                end
                if strcmp(variable_name, 'com_z_vel')
                    com_z = this.getBasicVariableData('com_z');
                    time = this.getTimeData('com_z');
                    sampling_rate = 1/median(diff(time));
                    com_z_vel = deriveByTime(nanfiltfilt(b, a, com_z), 1/sampling_rate);
                    this.basic_variable_data.com_z_vel = com_z_vel;
                    this.time_data.com_x_vel = time;
                end
                
                % joint angles
                if strcmp(variable_name, 'lumbar_roll_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'lumbar joint - sideways bending (right/left)');
                    this.basic_variable_data.lumbar_roll_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.lumbar_roll_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'lumbar_pitch_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'lumbar joint - forward/backward bending');
                    this.basic_variable_data.lumbar_pitch_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.lumbar_pitch_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'lumbar_yaw_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'lumbar joint - internal rotation (right/left)');
                    this.basic_variable_data.lumbar_yaw_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.lumbar_yaw_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'cevical_roll_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'cervical joint - sideways bending (right/left)');
                    this.basic_variable_data.cevical_roll_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.cevical_roll_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'cevical_pitch_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'cervical joint - forward/backward bending');
                    this.basic_variable_data.cevical_pitch_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.cevical_pitch_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'cevical_yaw_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'cervical joint - internal rotation (right/left)');
                    this.basic_variable_data.cevical_yaw_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.cevical_yaw_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'left_hip_abduction_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'left hip ab/adduction');
                    this.basic_variable_data.left_hip_abduction_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.left_hip_abduction_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'left_hip_flexion_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'left hip flexion/extension');
                    this.basic_variable_data.left_hip_flexion_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.left_hip_flexion_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'left_hip_introtation_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'left hip internal/external rotation');
                    this.basic_variable_data.left_hip_introtation_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.left_hip_introtation_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'left_knee_flexion_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'left knee flexion/extension');
                    this.basic_variable_data.left_knee_flexion_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.left_knee_flexion_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'left_knee_extrotation_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'left knee external/internal rotation');
                    this.basic_variable_data.left_knee_extrotation_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.left_knee_extrotation_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'left_ankle_inversion_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'left ankle inversion/eversion');
                    this.basic_variable_data.left_ankle_inversion_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.left_ankle_inversion_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'left_ankle_dorsiflexion_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'left ankle dorsi/plantarflexion');
                    this.basic_variable_data.left_ankle_dorsiflexion_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.left_ankle_dorsiflexion_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'right_hip_abduction_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'right hip ab/adduction');
                    this.basic_variable_data.right_hip_abduction_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.right_hip_abduction_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'right_hip_flexion_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'right hip flexion/extension');
                    this.basic_variable_data.right_hip_flexion_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.right_hip_flexion_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'right_hip_introtation_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'right hip internal/external rotation');
                    this.basic_variable_data.right_hip_introtation_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.right_hip_introtation_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'right_knee_flexion_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'right knee flexion/extension');
                    this.basic_variable_data.right_knee_flexion_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.right_knee_flexion_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'right_knee_extrotation_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'right knee external/internal rotation');
                    this.basic_variable_data.right_knee_extrotation_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.right_knee_extrotation_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'right_ankle_inversion_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'right ankle inversion/eversion');
                    this.basic_variable_data.right_ankle_inversion_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.right_ankle_inversion_angle = this.time_data.joint_angle_trajectories;
                end
                if strcmp(variable_name, 'right_ankle_dorsiflexion_angle')
                    joint_angle_indicator = strcmp(this.basic_variable_labels.joint_angle_trajectories, 'right ankle dorsi/plantarflexion');
                    this.basic_variable_data.right_ankle_dorsiflexion_angle = rad2deg(this.basic_variable_data.joint_angle_trajectories(:, joint_angle_indicator));
                    this.time_data.right_ankle_dorsiflexion_angle = this.time_data.joint_angle_trajectories;
                end
                
                % joint torques
                if strcmp(variable_name, 'lumbar_roll_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'lumbar joint - sideways bending (right/left)');
                    this.basic_variable_data.lumbar_roll_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.lumbar_roll_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'lumbar_pitch_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'lumbar joint - forward/backward bending');
                    this.basic_variable_data.lumbar_pitch_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.lumbar_pitch_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'lumbar_yaw_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'lumbar joint - internal rotation (right/left)');
                    this.basic_variable_data.lumbar_yaw_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.lumbar_yaw_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'cevical_roll_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'cervical joint - sideways bending (right/left)');
                    this.basic_variable_data.cevical_roll_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.cevical_roll_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'cevical_pitch_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'cervical joint - forward/backward bending');
                    this.basic_variable_data.cevical_pitch_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.cevical_pitch_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'cevical_yaw_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'cervical joint - internal rotation (right/left)');
                    this.basic_variable_data.cevical_yaw_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.cevical_yaw_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'left_hip_abduction_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'left hip ab/adduction');
                    this.basic_variable_data.left_hip_abduction_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.left_hip_abduction_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'left_hip_flexion_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'left hip flexion/extension');
                    this.basic_variable_data.left_hip_flexion_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.left_hip_flexion_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'left_hip_introtation_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'left hip internal/external rotation');
                    this.basic_variable_data.left_hip_introtation_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.left_hip_introtation_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'left_knee_flexion_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'left knee flexion/extension');
                    this.basic_variable_data.left_knee_flexion_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.left_knee_flexion_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'left_knee_extrotation_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'left knee external/internal rotation');
                    this.basic_variable_data.left_knee_extrotation_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.left_knee_extrotation_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'left_ankle_inversion_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'left ankle inversion/eversion');
                    this.basic_variable_data.left_ankle_inversion_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.left_ankle_inversion_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'left_ankle_dorsiflexion_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'left ankle dorsi/plantarflexion');
                    this.basic_variable_data.left_ankle_dorsiflexion_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.left_ankle_dorsiflexion_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'right_hip_abduction_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'right hip ab/adduction');
                    this.basic_variable_data.right_hip_abduction_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.right_hip_abduction_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'right_hip_flexion_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'right hip flexion/extension');
                    this.basic_variable_data.right_hip_flexion_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.right_hip_flexion_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'right_hip_introtation_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'right hip internal/external rotation');
                    this.basic_variable_data.right_hip_introtation_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.right_hip_introtation_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'right_knee_flexion_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'right knee flexion/extension');
                    this.basic_variable_data.right_knee_flexion_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.right_knee_flexion_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'right_knee_extrotation_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'right knee external/internal rotation');
                    this.basic_variable_data.right_knee_extrotation_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.right_knee_extrotation_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'right_ankle_inversion_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'right ankle inversion/eversion');
                    this.basic_variable_data.right_ankle_inversion_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.right_ankle_inversion_torque = this.time_data.joint_torque_trajectories_belt;
                end
                if strcmp(variable_name, 'right_ankle_dorsiflexion_torque')
                    joint_torque_indicator = strcmp(this.basic_variable_labels.joint_torque_trajectories_belt, 'right ankle dorsi/plantarflexion');
                    this.basic_variable_data.right_ankle_dorsiflexion_torque = rad2deg(this.basic_variable_data.joint_torque_trajectories_belt(:, joint_torque_indicator));
                    this.time_data.right_ankle_dorsiflexion_torque = this.time_data.joint_torque_trajectories_belt;
                end
                
                % force plate
                if strcmp(variable_name, 'copl_y')
                    left_foot_cop_world = this.getBasicVariableData('left_foot_cop_world');
                    this.basic_variable_data.copl_y = left_foot_cop_world(:, 2);
                    this.time_data.copl_y = this.time_data.left_foot_cop_world;
                end
                if strcmp(variable_name, 'copl_x')
                    left_foot_cop_world = this.getBasicVariableData('left_foot_cop_world');
                    this.basic_variable_data.copl_x = left_foot_cop_world(:, 1);
                    this.time_data.copl_x = this.time_data.left_foot_cop_world;
                end
                if strcmp(variable_name, 'copr_y')
                    right_foot_cop_world = this.getBasicVariableData('right_foot_cop_world');
                    this.basic_variable_data.copr_y = right_foot_cop_world(:, 2);
                    this.time_data.copr_y = this.time_data.right_foot_cop_world;
                end
                if strcmp(variable_name, 'copr_x')
                    right_foot_cop_world = this.getBasicVariableData('right_foot_cop_world');
                    this.basic_variable_data.copr_x = right_foot_cop_world(:, 1);
                    this.time_data.copr_x = this.time_data.right_foot_cop_world;
                end
                if strcmp(variable_name, 'cop_y')
                    total_forceplate_cop_world = this.getBasicVariableData('total_forceplate_cop_world');
                    this.basic_variable_data.cop_y = total_forceplate_cop_world(:, 2);
                    this.time_data.cop_y = this.time_data.total_forceplate_cop_world;
                end
                if strcmp(variable_name, 'cop_x')
                    total_forceplate_cop_world = this.getBasicVariableData('total_forceplate_cop_world');
                    this.basic_variable_data.cop_x = total_forceplate_cop_world(:, 1);
                    this.time_data.cop_x = this.time_data.total_forceplate_cop_world;
                end
                if strcmp(variable_name, 'cop_y_vel')
                    cop_y = this.getBasicVariableData('cop_y');
                    cop_y(cop_y==0) = NaN;
                    time = this.getTimeData('cop_y');
                    filter_order = this.study_settings.force_plate_derivative_filter_order;
                    cutoff_frequency = this.study_settings.force_plate_derivative_filter_cutoff;
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;
                    cop_y_vel = deriveByTime(nanfiltfilt(b, a, cop_y), 1/sampling_rate);
                    cop_y_vel(isnan(cop_y_vel)) = 0;
                    this.basic_variable_data.cop_y_vel = cop_y_vel;
                    this.time_data.cop_y_vel = time;
                end
                if strcmp(variable_name, 'cop_y_acc')
                    cop_y_vel = this.getBasicVariableData('cop_y_vel');
                    cop_y_vel(cop_y_vel==0) = NaN;
                    time = this.getTimeData('cop_y_vel');
                    filter_order = this.study_settings.force_plate_derivative_filter_order;
                    cutoff_frequency = this.study_settings.force_plate_derivative_filter_cutoff;
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;
                    cop_y_acc = deriveByTime(nanfiltfilt(b, a, cop_y_vel), 1/sampling_rate);
                    cop_y_acc(isnan(cop_y_acc)) = 0;
                    this.basic_variable_data.cop_y_acc = cop_y_acc;
                    this.time_data.cop_y_acc = time;
                end
                if strcmp(variable_name, 'cop_x_vel')
                    cop_x = this.getBasicVariableData('cop_x');
                    cop_x(cop_x==0) = NaN;
                    time = this.getTimeData('cop_x');
                    filter_order = this.study_settings.force_plate_derivative_filter_order;
                    cutoff_frequency = this.study_settings.force_plate_derivative_filter_cutoff;
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;
                    cop_x_vel = deriveByTime(nanfiltfilt(b, a, cop_x), 1/sampling_rate);
                    cop_x_vel(isnan(cop_x_vel)) = 0;
                    this.basic_variable_data.cop_x_vel = cop_x_vel;
                    this.time_data.cop_x_vel = time;
                end
                if strcmp(variable_name, 'cop_x_acc')
                    cop_x_vel = this.getBasicVariableData('cop_x_vel');
                    cop_x_vel(cop_x_vel==0) = NaN;
                    time = this.getTimeData('cop_x_vel');
                    filter_order = this.study_settings.force_plate_derivative_filter_order;
                    cutoff_frequency = this.study_settings.force_plate_derivative_filter_cutoff;
                    sampling_rate = 1/median(diff(time));
                    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;
                    cop_x_acc = deriveByTime(nanfiltfilt(b, a, cop_x_vel), 1/sampling_rate);
                    cop_x_acc(isnan(cop_x_acc)) = 0;
                    this.basic_variable_data.cop_x_acc = cop_x_acc;
                    this.time_data.cop_x_acc = time;
                end
                if strcmp(variable_name, 'fxl')
                    left_foot_wrench_world = this.getBasicVariableData('left_foot_wrench_world');
                    this.basic_variable_data.fxl = left_foot_wrench_world(:, 1);
                    this.time_data.fxl = this.time_data.left_foot_wrench_world;
                end
                if strcmp(variable_name, 'fyl')
                    left_foot_wrench_world = this.getBasicVariableData('left_foot_wrench_world');
                    this.basic_variable_data.fyl = left_foot_wrench_world(:, 2);
                    this.time_data.fyl = this.time_data.left_foot_wrench_world;
                end
                if strcmp(variable_name, 'fzl')
                    left_foot_wrench_world = this.getBasicVariableData('left_foot_wrench_world');
                    this.basic_variable_data.fzl = left_foot_wrench_world(:, 3);
                    this.time_data.fzl = this.time_data.left_foot_wrench_world;
                end
                if strcmp(variable_name, 'mxl')
                    left_foot_wrench_world = this.getBasicVariableData('left_foot_wrench_world');
                    this.basic_variable_data.mxl = left_foot_wrench_world(:, 4);
                    this.time_data.mxl = this.time_data.left_foot_wrench_world;
                end
                if strcmp(variable_name, 'myl')
                    left_foot_wrench_world = this.getBasicVariableData('left_foot_wrench_world');
                    this.basic_variable_data.myl = left_foot_wrench_world(:, 5);
                    this.time_data.myl = this.time_data.left_foot_wrench_world;
                end
                if strcmp(variable_name, 'mzl')
                    left_foot_wrench_world = this.getBasicVariableData('left_foot_wrench_world');
                    this.basic_variable_data.mzl = left_foot_wrench_world(:, 6);
                    this.time_data.mzl = this.time_data.left_foot_wrench_world;
                end
                if strcmp(variable_name, 'fxr')
                    right_foot_wrench_world = this.getBasicVariableData('right_foot_wrench_world');
                    this.basic_variable_data.fxr = right_foot_wrench_world(:, 1);
                    this.time_data.fxr = this.time_data.right_foot_wrench_world;
                end
                if strcmp(variable_name, 'fyr')
                    right_foot_wrench_world = this.getBasicVariableData('right_foot_wrench_world');
                    this.basic_variable_data.fyr = right_foot_wrench_world(:, 2);
                    this.time_data.fyr = this.time_data.right_foot_wrench_world;
                end
                if strcmp(variable_name, 'fzr')
                    right_foot_wrench_world = this.getBasicVariableData('right_foot_wrench_world');
                    this.basic_variable_data.fzr = right_foot_wrench_world(:, 3);
                    this.time_data.fzr = this.time_data.right_foot_wrench_world;
                end
                if strcmp(variable_name, 'mxr')
                    right_foot_wrench_world = this.getBasicVariableData('right_foot_wrench_world');
                    this.basic_variable_data.mxr = right_foot_wrench_world(:, 4);
                    this.time_data.mxr = this.time_data.right_foot_wrench_world;
                end
                if strcmp(variable_name, 'myr')
                    right_foot_wrench_world = this.getBasicVariableData('right_foot_wrench_world');
                    this.basic_variable_data.myr = right_foot_wrench_world(:, 5);
                    this.time_data.myr = this.time_data.right_foot_wrench_world;
                end
                if strcmp(variable_name, 'mzr')
                    right_foot_wrench_world = this.getBasicVariableData('right_foot_wrench_world');
                    this.basic_variable_data.mzr = right_foot_wrench_world(:, 6);
                    this.time_data.mzr = this.time_data.right_foot_wrench_world;
                end
                if strcmp(variable_name, 'fx')
                    total_forceplate_wrench_world = this.getBasicVariableData('total_forceplate_wrench_world');
                    this.basic_variable_data.fx = total_forceplate_wrench_world(:, 1);
                    this.time_data.fx = this.time_data.total_forceplate_wrench_world;
                end
                if strcmp(variable_name, 'fy')
                    total_forceplate_wrench_world = this.getBasicVariableData('total_forceplate_wrench_world');
                    this.basic_variable_data.fy = total_forceplate_wrench_world(:, 2);
                    this.time_data.fy = this.time_data.total_forceplate_wrench_world;
                end
                if strcmp(variable_name, 'fz')
                    total_forceplate_wrench_world = this.getBasicVariableData('total_forceplate_wrench_world');
                    this.basic_variable_data.fz = total_forceplate_wrench_world(:, 3);
                    this.time_data.fz = this.time_data.total_forceplate_wrench_world;
                end
                if strcmp(variable_name, 'mx')
                    total_forceplate_wrench_world = this.getBasicVariableData('total_forceplate_wrench_world');
                    this.basic_variable_data.mx = total_forceplate_wrench_world(:, 4);
                    this.time_data.mx = this.time_data.total_forceplate_wrench_world;
                end
                if strcmp(variable_name, 'my')
                    total_forceplate_wrench_world = this.getBasicVariableData('total_forceplate_wrench_world');
                    this.basic_variable_data.my = total_forceplate_wrench_world(:, 5);
                    this.time_data.my = this.time_data.total_forceplate_wrench_world;
                end
                if strcmp(variable_name, 'mz')
                    total_forceplate_wrench_world = this.getBasicVariableData('total_forceplate_wrench_world');
                    this.basic_variable_data.mz = total_forceplate_wrench_world(:, 6);
                    this.time_data.mz = this.time_data.total_forceplate_wrench_world;
                end
                % emg
                if strcmp(variable_name, 'left_glut_med')
                    this.basic_variable_data.left_glut_med = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'left_glut_med'));
                    this.time_data.left_glut_med = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'left_delt_ant')
                    this.basic_variable_data.left_delt_ant = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'left_delt_ant'));
                    this.time_data.left_delt_ant = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'left_tibi_ant')
                    this.basic_variable_data.left_tibi_ant = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'left_tibi_ant'));
                    this.time_data.left_tibi_ant = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'left_gastroc_med')
                    this.basic_variable_data.left_gastroc_med = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'left_gastroc_med'));
                    this.time_data.left_gastroc_med = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'left_pero_lng')
                    this.basic_variable_data.left_pero_lng = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'left_pero_lng'));
                    this.time_data.left_pero_lng = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'left_tfl')
                    this.basic_variable_data.left_tfl = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'left_tfl'));
                    this.time_data.left_tfl = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'right_glut_med')
                    this.basic_variable_data.right_glut_med = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'right_glut_med'));
                    this.time_data.right_glut_med = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'right_delt_ant')
                    this.basic_variable_data.right_delt_ant = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'right_delt_ant'));
                    this.time_data.right_delt_ant = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'right_tibi_ant')
                    this.basic_variable_data.right_tibi_ant = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'right_tibi_ant'));
                    this.time_data.right_tibi_ant = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'right_gastroc_med')
                    this.basic_variable_data.right_gastroc_med = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'right_gastroc_med'));
                    this.time_data.right_gastroc_med = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'right_pero_lng')
                    this.basic_variable_data.right_pero_lng = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'right_pero_lng'));
                    this.time_data.right_pero_lng = this.time_data.emg_trajectories;
                end
                if strcmp(variable_name, 'right_tfl')
                    this.basic_variable_data.right_tfl = this.basic_variable_data.emg_trajectories(:, strcmp(this.basic_variable_labels.emg_trajectories, 'right_tfl'));
                    this.time_data.right_tfl = this.time_data.emg_trajectories;
                end
            end
        end
        function stretch_variables = calculateStretchVariables(this, stretch_start_times, stretch_end_times, condition_stance_foot_list, variables_to_calculate)
            if nargin < 5
                variables_to_calculate = this.stretch_variable_names;
            end
            
            number_of_stretch_variables = length(variables_to_calculate);
            number_of_stretches = length(stretch_start_times);
            stretch_variables = cell(number_of_stretch_variables, 1);
            
            for i_variable = 1 : number_of_stretch_variables
                variable_name = variables_to_calculate{i_variable};
                
                % extract and normalize data from stretches
                for i_stretch = 1 : number_of_stretches
                    stretch_data = [];
                    
                    % time
                    this_stretch_start_time = stretch_start_times(i_stretch);
                    this_stretch_end_time = stretch_end_times(i_stretch);
                    
                    % calculate normalized stretch data for the basic variables
                    if this.isBasicVariable(variable_name)
                        stretch_data = this.getTimeNormalizedData(variable_name, this_stretch_start_time, this_stretch_end_time);
                    end
                
                    % calculate stretch variables that are not basic variables
                    if strcmp(variable_name, 'step_length')
                        lheel_y = this.getTimeNormalizedData('lheel_y', this_stretch_start_time, this_stretch_end_time);
                        rheel_y = this.getTimeNormalizedData('rheel_y', this_stretch_start_time, this_stretch_end_time);
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')
                            stretch_data = lheel_y(end) - rheel_y(end);
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')
                            stretch_data = rheel_y(end) - lheel_y(end);
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_BOTH')
                            stretch_data = NaN;
                        end
                    end
                    if strcmp(variable_name, 'step_width')
                        lheel_x = this.getTimeNormalizedData('lheel_x', this_stretch_start_time, this_stretch_end_time);
                        rheel_x = this.getTimeNormalizedData('rheel_x', this_stretch_start_time, this_stretch_end_time);
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')
                            stretch_data = abs(lheel_x(end) - rheel_x(1));
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')
                            stretch_data = abs(rheel_x(end) - lheel_x(1));
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_BOTH')
                            stretch_data = NaN;
                        end
                    end
                    if strcmp(variable_name, 'step_placement_x')
                        lheel_x = this.getTimeNormalizedData('lheel_x', this_stretch_start_time, this_stretch_end_time);
                        rheel_x = this.getTimeNormalizedData('rheel_x', this_stretch_start_time, this_stretch_end_time);
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')
                            stretch_data = lheel_x(end) - rheel_x(1);
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')
                            stretch_data = rheel_x(end) - lheel_x(1);
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_BOTH')
                            stretch_data = NaN;
                        end
                    end
                    if strcmp(variable_name, 'step_time')
                        stretch_data = this_stretch_end_time - this_stretch_start_time;
                    end
                    if strcmp(variable_name, 'cadence')
                        second_to_minute = 1/60;
                        step_time = stretch_variables{strcmp(this.stretch_variable_names, 'step_time')}(i_stretch);
                        stretch_data = (step_time * second_to_minute)^(-1);
                    end
                    if strcmp(variable_name, 'velocity')
                        step_time = stretch_variables{strcmp(this.stretch_variable_names, 'step_time')}(i_stretch);
                        step_length = stretch_variables{strcmp(this.stretch_variable_names, 'step_length')}(i_stretch);
                        stretch_data = step_length / step_time;
                    end
                    if strcmp(variable_name, 'left_arm_phase') || strcmp(variable_name, 'right_arm_phase') || strcmp(variable_name, 'left_leg_phase') || strcmp(variable_name, 'right_leg_phase')
                        % extract data
                        variable_time = this.getTimeData(variable_name);
                        variable_data = this.getBasicVariableData(variable_name);
                        [~, start_index] = min(abs(variable_time - this_stretch_start_time));
                        [~, end_index] = min(abs(variable_time - this_stretch_end_time));
                        time_extracted = variable_time(start_index : end_index);
                        data_extracted = variable_data(start_index : end_index);
                        
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
                        
                        % normalize time
                        time_normalized = linspace(time_extracted(1), time_extracted(end), this.number_of_time_steps_normalized)';
                        data_normalized = spline(time_extracted, data_extracted_groomed, time_normalized);
                        
                        % re-normalize angle
                        stretch_data = normalizeAngle(data_normalized);
                    end
                    if strcmp(variable_name, 'left_arm_phase_at_heelstrike')
                        left_arm_phase = stretch_variables{strcmp(variables_to_calculate, 'left_arm_phase')}(:, i_stretch);
                        stretch_data = left_arm_phase(1);
                        
                        % normalize angle, but take care that the leaf change is not too close to this value
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')
                            stretch_data = normalizeAngle(stretch_data, 0);
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')
                            stretch_data = normalizeAngle(stretch_data, pi);
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_BOTH')
                            stretch_data = NaN;
                        end
                    end
                    if strcmp(variable_name, 'right_arm_phase_at_heelstrike')
                        right_arm_phase = stretch_variables{strcmp(variables_to_calculate, 'right_arm_phase')}(:, i_stretch);
                        stretch_data = right_arm_phase(1);
                        
                        % normalize angle, but take care that the leaf change is not too close to this value
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_RIGHT')
                            stretch_data = normalizeAngle(stretch_data, pi);
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_LEFT')
                            stretch_data = normalizeAngle(stretch_data, 0);
                        end
                        if strcmp(condition_stance_foot_list{i_stretch}, 'STANCE_BOTH')
                            stretch_data = NaN;
                        end
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
                    if strcmp(variable_name, 'left_glut_med_rescaled')
                        left_glut_med = this.getTimeNormalizedData('left_glut_med', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'left_glut_med'));
                        stretch_data = left_glut_med * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'left_delt_ant_rescaled')
                        left_delt_ant = this.getTimeNormalizedData('left_delt_ant', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'left_delt_ant'));
                        stretch_data = left_delt_ant * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'left_tibi_ant_rescaled')
                        left_tibi_ant = this.getTimeNormalizedData('left_tibi_ant', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'left_tibi_ant'));
                        stretch_data = left_tibi_ant * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'left_gastroc_med_rescaled')
                        left_gastroc_med = this.getTimeNormalizedData('left_gastroc_med', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'left_gastroc_med'));
                        stretch_data = left_gastroc_med * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'left_pero_lng_rescaled')
                        left_pero_lng = this.getTimeNormalizedData('left_pero_lng', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'left_pero_lng'));
                        stretch_data = left_pero_lng * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'left_tfl_rescaled')
                        left_tfl = this.getTimeNormalizedData('left_tfl', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'left_tfl'));
                        stretch_data = left_tfl * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'right_glut_med_rescaled')
                        right_glut_med = this.getTimeNormalizedData('right_glut_med', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'right_glut_med'));
                        stretch_data = right_glut_med * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'right_delt_ant_rescaled')
                        right_delt_ant = this.getTimeNormalizedData('right_delt_ant', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'right_delt_ant'));
                        stretch_data = right_delt_ant * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'right_tibi_ant_rescaled')
                        right_tibi_ant = this.getTimeNormalizedData('right_tibi_ant', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'right_tibi_ant'));
                        stretch_data = right_tibi_ant * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'right_gastroc_med_rescaled')
                        right_gastroc_med = this.getTimeNormalizedData('right_gastroc_med', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'right_gastroc_med'));
                        stretch_data = right_gastroc_med * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'right_pero_lng_rescaled')
                        right_pero_lng = this.getTimeNormalizedData('right_pero_lng', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'right_pero_lng'));
                        stretch_data = right_pero_lng * 1 / normalization_value;
                    end
                    if strcmp(variable_name, 'right_tfl_rescaled')
                        right_tfl = this.getTimeNormalizedData('right_tfl', this_stretch_start_time, this_stretch_end_time);
                        normalization_value = this.emg_normalization_values(strcmp(this.emg_normalization_labels, 'right_tfl'));
                        stretch_data = right_tfl * 1 / normalization_value;
                    end
                    
                    % store in cell
                    stretch_variables{i_variable} = [stretch_variables{i_variable} stretch_data];
                end
                
            end
            
            
        end
        function data_normalized = getTimeNormalizedData(this, variable_name, start_time, end_time)
            % extract data
            variable_time = this.getTimeData(variable_name);
            variable_data = this.getBasicVariableData(variable_name);
            [~, start_index] = min(abs(variable_time - start_time));
            [~, end_index] = min(abs(variable_time - end_time));
            time_extracted = variable_time(start_index : end_index);
            data_extracted = variable_data(start_index : end_index);
                
            % normalize data in time
            if ~isempty(time_extracted)
                time_normalized = linspace(time_extracted(1), time_extracted(end), this.number_of_time_steps_normalized)';
                data_normalized = spline(time_extracted, data_extracted, time_normalized);
            else
                data_normalized = zeros(this.number_of_time_steps_normalized, 1) * NaN;
            end
        end
    end
end







