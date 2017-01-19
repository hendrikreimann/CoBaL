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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% input

% description

classdef WalkingDataCustodian < handle
    properties
        date = [];
        subject_id = [];
        variables_to_analyze = {};
        basic_variable_names = {};
        required_variable_names = {};
        
        % these are private and should only be accessed using get functions
        number_of_time_steps_normalized;
        basic_variable_data;
        basic_variable_labels;
        required_variable_data;
        time_data;
    end
    
    methods
        % constructor
%         function this = WalkingDataCustodian(date, subject_id, variables_to_analyze)
        function this = WalkingDataCustodian()
            % load this information from the subjects.mat and studySettings.txt files
            load('subjectInfo.mat', 'date', 'subject_id');
            study_settings = loadSettingsFile(['..' filesep 'studySettings.txt']);
            
            this.date = date;
            this.subject_id = subject_id;
            this.variables_to_analyze = study_settings.variables_to_analyze;
            this.number_of_time_steps_normalized = study_settings.number_of_time_steps_normalized;
            
            this.determineVariables();
        end
        
        % initialization
        function addBasicVariable(this, variable_name)
            if ~any(strcmp(this.basic_variable_names, variable_name))
                this.basic_variable_names = [this.basic_variable_names; variable_name];
            end
        end
        function addRequiredVariable(this, variable_name)
            if ~any(strcmp(this.required_variable_names, variable_name))
                this.required_variable_names = [this.required_variable_names; variable_name];
            end
            
            % TODO: make sure that after I added this one, all the variables are in the right order, so that
            % variables that are needed later are calculated first, e.g. c7 pos and MPSIS pos are calculated before
            % trunk angle
            
        end
        function determineVariables(this)
            % for each possible variable to analyze, list the basic and required variables required to calculate it
            if this.variablesToAnalyzeContains('lheel_x_pos')
                this.addBasicVariable('marker_trajectories')
                this.addRequiredVariable('lheel_x_pos')
            end
            if this.variablesToAnalyzeContains('rheel_x_pos')
                this.addBasicVariable('marker_trajectories')
                this.addRequiredVariable('rheel_x_pos')
            end
        end
        
        % interface
        function contains = variablesToAnalyzeContains(this, variable_name)
            contains = any(strcmp(this.variables_to_analyze, variable_name));
        end
        function time_data = getTimeData(this, variable_name)
            eval(['time_data = this.time_data.' variable_name ';']);
        end
        function variable_data = getRequiredVariableData(this, variable_name)
            eval(['variable_data = this.required_variable_data.' variable_name ';']);
        end
        
        function prepareData(this, condition, trial)
            % clear out old data
            this.basic_variable_data = struct;
            this.basic_variable_labels = struct;
            this.required_variable_data = struct;
            this.time_data = struct;
            
            % prepare the data by loading all the basic variables from disk and calculating the required variables
            load(['analysis' filesep makeFileName(this.date, this.subject_id, condition, trial, 'availableVariables')]);
            
            % load basic variables
            for i_variable = 1 : length(this.basic_variable_names)
                variable_name = this.basic_variable_names{i_variable};
                
                % load data
                [data, time, sampling_rate, labels] = loadData(this.date, this.subject_id, condition, trial, variable_name);
                
                % store
                eval(['this.basic_variable_data.' variable_name ' = data;']);
                eval(['this.time_data.' variable_name ' = time;']);
                eval(['this.basic_variable_labels.' variable_name ' = labels;']);
            end
            
            % calculate required variables
            for i_variable = 1 : length(this.required_variable_names)
                variable_name = this.required_variable_names{i_variable};
                this.calculateRequiredVariable(variable_name);
            end
        end
        function calculateRequiredVariable(this, variable_name)
            if strcmp(variable_name, 'lheel_x_pos')
                LHEE_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'LHEE');
                variable_data = LHEE_trajectory(:, 1);
                variable_time = this.time_data.marker_trajectories;
                eval(['this.required_variable_data.' variable_name ' = variable_data;']);
                eval(['this.time_data.' variable_name ' = variable_time;']);
            end
            if strcmp(variable_name, 'rheel_x_pos')
                LHEE_trajectory = extractMarkerTrajectories(this.basic_variable_data.marker_trajectories, this.basic_variable_labels.marker_trajectories, 'RHEE');
                variable_data = LHEE_trajectory(:, 1);
                variable_time = this.time_data.marker_trajectories;
                eval(['this.required_variable_data.' variable_name ' = variable_data;']);
                eval(['this.time_data.' variable_name ' = variable_time;']);
            end
        end
        
        function data_normalized = getTimeNormalizedData(this, variable_name, start_time, end_time)
            % extract data
            variable_time = this.getTimeData(variable_name);
            variable_data = this.getRequiredVariableData(variable_name);
            [~, start_index] = min(abs(variable_time - start_time));
            [~, end_index] = min(abs(variable_time - end_time));
            time_extracted = variable_time(start_index : end_index);
            data_extracted = variable_data(start_index : end_index);
                
            % normalize data in time
            time_normalized = linspace(time_extracted(1), time_extracted(end), this.number_of_time_steps_normalized)';
            data_normalized = spline(time_extracted, data_extracted, time_normalized);
        end
    end
end







