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


function processInverseKinematicsInOpensim(varargin)
    [trial_type_list, trial_number_list] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    parse(parser, varargin{:})
    
    % load settings
    study_settings_file = '';
    if exist(['..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep 'studySettings.txt'];
    end    
    if exist(['..' filesep '..' filesep 'studySettings.txt'], 'file')
        study_settings_file = ['..' filesep '..' filesep 'studySettings.txt'];
    end
    study_settings = SettingsCustodian(study_settings_file);
    subject_settings = SettingsCustodian('subjectSettings.txt');
    subject_info = load('subjectInfo.mat');
    
    %% set up
    import org.opensim.modeling.*
    generic_setup_file_ik = [getCobalPath filesep 'resources' filesep 'opensim' filesep 'CoBaLWalker50_setupInverseKinematics.xml'];
    data_root = [pwd filesep 'opensim'];
    ikTool = InverseKinematicsTool(generic_setup_file_ik);
    model_file = [data_root filesep subject_info.date '_' subject_info.subject_id '.osim'];

    % Load the model and initialize
    model = Model(model_file);
    model.initSystem();

    % Tell Tool to use the loaded model
    ikTool.setModel(model);

    %% process
    for i_condition = 1 : length(trial_type_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            %% load data
            trial_type = trial_type_list{i_condition};
            marker_file = [data_root filesep 'marker' filesep makeFileName(subject_info.date, subject_info.subject_id, trial_type, i_trial, 'marker.trc')];
            marker_data = MarkerData(marker_file);
            output_file = [data_root filesep 'inverseKinematics' filesep makeFileName(subject_info.date, subject_info.subject_id, trial_type, i_trial, 'inverseKinematics.mot')];
            
            % Get initial and final time 
            time_start = marker_data.getStartFrameTime();
            time_end = marker_data.getLastFrameTime();

            % Setup the ikTool for this trial
            ikTool.setName(makeFileName(subject_info.date, subject_info.subject_id, trial_type, i_trial, 'marker.trc'));
            ikTool.setMarkerDataFileName(marker_file);
            ikTool.setStartTime(time_start);
            ikTool.setEndTime(time_end);
            ikTool.setOutputMotionFileName(output_file);

            % Save the settings in a setup file
            setup_file = [data_root filesep 'setupFiles' filesep makeFileName(subject_info.date, subject_info.subject_id, trial_type, i_trial, 'setupIK.xml')];
            ikTool.print(setup_file);

            % Run IK
            ikTool.run();
            disp(['Performing inverse kinematics: condition ' trial_type_list{i_condition} ', Trial ' num2str(i_trial) ' completed']);                
            
            % move log file
            error_file_source = [pwd filesep makeFileName(subject_info.date, subject_info.subject_id, trial_type, i_trial, 'marker.trc') '_ik_marker_errors.sto'];
            error_file_destination = [data_root filesep 'logs' filesep makeFileName(subject_info.date, subject_info.subject_id, trial_type, i_trial, 'ikErrors.sto')];
            movefile(error_file_source, error_file_destination);
            
            
        end
    end



end