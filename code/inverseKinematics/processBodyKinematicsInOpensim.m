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


function processBodyKinematicsInOpensim(varargin)
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
    subject_settings = loadSettingsFromFile('subject');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');

    if ~directoryExists(['opensim' filesep 'bodyKinematics'])
        mkdir(['opensim' filesep 'bodyKinematics'])
    end
    
    %% set up
    import org.opensim.modeling.*
    generic_setup_file_analysis = [getCobalPath filesep 'resources' filesep 'opensim' filesep 'CoBaLWalker50_setupAnalysis.xml'];
    DOMnode = xmlread(generic_setup_file_analysis);
    document_node = DOMnode.getFirstChild;
    
    % set name of resulting model
    analyze_tool_node = document_node.getElementsByTagName("AnalyzeTool").item(0);
    analyze_tool_node.setAttribute("name", makeFileName(collection_date, subject_id, 'setupAnalysis.xml'))
    
    % set results directory
    analyze_tool_node.getElementsByTagName("model_file").item(0).getFirstChild.setData(makeFileName(collection_date, subject_id, '.osim'));
    analyze_tool_node.getElementsByTagName("force_set_files").item(0).getFirstChild.setData('');
    analyze_tool_node.getElementsByTagName("results_directory").item(0).getFirstChild.setData(['opensim' filesep 'bodyKinematics']);
    
    % save subject-specific scale file
    setup_file_analysis = [pwd filesep 'opensim' filesep makeFileName(collection_date, subject_id, 'setupAnalysis.xml')];
    xmlwrite(setup_file_analysis, DOMnode);       
    
    %% process
    analysis_tool = AnalyzeTool(setup_file_analysis);
    
    data_root = [pwd filesep 'opensim'];
    for i_condition = 1 : length(trial_type_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            % load data to get initial and final time
            trial_type = trial_type_list{i_condition};
            marker_file = [data_root filesep '..' filesep 'processed' filesep makeFileName(collection_date, subject_id, trial_type, i_trial, 'markerTrajectories.mat')];
            load(marker_file, 'time_mocap');
            time_start = time_mocap(1);
            time_end = time_mocap(end);
            analysis_tool.setInitialTime(time_start);
            analysis_tool.setFinalTime(time_end);
            
            % run the tool for this trial
%             inverse_kinematics_file_name = ['opensim' filesep 'inverseKinematics' filesep makeFileName(collection_date, subject_id, trial_type, i_trial, 'inverseKinematics.mot')];
            inverse_kinematics_file_name = [data_root filesep 'inverseKinematics' filesep makeFileName(collection_date, subject_id, trial_type, i_trial, 'inverseKinematics.mot')];
            analysis_tool.setCoordinatesFileName(inverse_kinematics_file_name);
            analysis_tool.run();
            
            % print and rename results
            output_file = ['opensim' filesep 'bodyKinematics' filesep makeFileName(collection_date, subject_id, trial_type, i_trial, 'bodyKinematics')];
            analysis_tool.printResults(output_file);
            movefile([output_file '_BodyKinematics_pos_global.sto'], [output_file 'Pos.sto']);
            movefile([output_file '_BodyKinematics_vel_global.sto'], [output_file 'Vel.sto']);
            movefile([output_file '_BodyKinematics_acc_global.sto'], [output_file 'Acc.sto']);
        end
    end



end