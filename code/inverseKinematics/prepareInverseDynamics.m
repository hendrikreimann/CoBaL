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

function prepareInverseDynamics(varargin)
    % parse arguments
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', false)
    parse(parser, varargin{:})
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    trials_to_process_table = makeTrialsToProcessTable(condition_list, trial_number_list);

    %% prepare
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');

    %% process
    
    for i_trial = 1 : height(trials_to_process_table)
        % create container for data
        trial_data = struct;
        
        % fill in basics
        trial_data.trial_type = trials_to_process_table{i_trial, 'trial type'};
        trial_data.trial_number = trials_to_process_table{i_trial, 'trial number'};

        % load belt speed data and integrate to belt translation
        collection_date = subject_settings.get('collection_date');
        subject_id = subject_settings.get('subject_id');
        plc_data = load(['processed' filesep makeFileName(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'PLCData.mat')]);
        belt_velocity = (plc_data.belt_speed_left_trajectory + plc_data.belt_speed_right_trajectory) * 0.5;
        belt_translation = cumtrapz(plc_data.time, belt_velocity);
        
        % load inverse kinematics
        load_file_name = ['opensim' filesep 'inverseKinematics' filesep makeFileName(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'inverseKinematics.mot')];
        ik_data = readOpensimMot(load_file_name, 11);
        
        % resample belt translation to inverse kinematics time
        belt_translation_resampled = spline(plc_data.time, belt_translation, ik_data.time);
        
        % add belt translation to ap-component of pelvis
        pelvis_ap_column = strcmp(ik_data.labels, 'pelvis_tx');
        ik_data.trajectories(:, pelvis_ap_column) = ik_data.trajectories(:, pelvis_ap_column) + belt_translation_resampled;
        
        % save inverse kinematics
        save_file_name = makeFileName(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'inverseKinematics.mot');
        save_file_path = [pwd filesep 'opensim' filesep 'inverseDynamics' filesep];
        writeOpensimMot(save_file_name, save_file_path, ik_data.time, ik_data.trajectories, ik_data.labels, ik_data.header);
        
        % load forceplate data
        load_file_name = ['opensim' filesep 'forceplate' filesep makeFileName(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'grfIndividual.mot')];
        force_data = readOpensimMot(load_file_name, 6);
        
        % add belt translation to ap-component of CoP
        left_plate_cop_ap_column = strcmp(force_data.labels, 'left_plate_force_px');
        right_plate_cop_ap_column = strcmp(force_data.labels, 'right_plate_force_px');
        force_data.trajectories(:, left_plate_cop_ap_column) = force_data.trajectories(:, left_plate_cop_ap_column) + belt_translation_resampled;
        force_data.trajectories(:, right_plate_cop_ap_column) = force_data.trajectories(:, right_plate_cop_ap_column) + belt_translation_resampled;
        
        % save updated forceplate data in 
        save_file_name = makeFileName(collection_date, subject_id, trial_data.trial_type, trial_data.trial_number, 'grfIndividual.mot');
        writeOpensimMot(save_file_name, save_file_path, force_data.time, force_data.trajectories, force_data.labels, force_data.header);        
        
        disp(['processed ' load_file_name ' and saved as inverseDynamics' filesep save_file_name])
    end
    return
    %% transform to belt space
    for i_condition = 1 : length(types_to_analyze)
        trials_to_process = trials_to_analyze{i_condition};
        for i_trial = trials_to_process
            % load data
            this_condition = types_to_analyze{i_condition};
            if any(strcmp(study_settings.get('conditions_to_transform_to_belt_space', 1), this_condition))
                % extract data for new structure
                if exist(['processed' filesep makeFileName(date, subject_id, trial_type, i_trial, 'plcData')], 'file')
                    load(['processed' filesep makeFileName(date, subject_id, trial_type, i_trial, 'plcData')])
                else
                    error(['Failed to load PLC data file for condition ' trial_type ', trial ' num2str(i_trial)]);
                end
                time_belts = time_plcData - time_plcData(1);

                % calculate shift
                belt_speed_trajectory = mean([belt_speed_left_trajectory belt_speed_right_trajectory], 2);
                delta_t = diff(time_belts);
                belt_position_trajectory_plcData = zeros(size(belt_speed_trajectory));
                for i_time = 2 : length(belt_speed_trajectory)
                    belt_position_trajectory_plcData(i_time) = belt_position_trajectory_plcData(i_time-1) + delta_t(i_time-1) * belt_speed_trajectory(i_time-1);
                end

                % apply shift to marker trajectories
                file_name_raw = ['raw' filesep makeFileName(date, subject_id, this_condition, i_trial, 'markerTrajectoriesRaw.mat')];
                load(file_name_raw);
                marker_trajectories = marker_trajectories_raw;
                belt_position_trajectory_mocap = spline(time_belts, belt_position_trajectory_plcData, time_mocap)';
                for i_marker = 1 : size(marker_headers, 2)
                    marker_trajectories(:, (i_marker-1)*3+2) = marker_trajectories(:, (i_marker-1)*3+2) + belt_position_trajectory_mocap;
                end

                save_folder = 'processed';
                save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'markerTrajectories.mat');
                save ...
                  ( ...
                    [save_folder filesep save_file_name], ...
                    'marker_trajectories', ...
                    'time_mocap', ...
                    'sampling_rate_mocap', ...
                    'marker_labels' ...
                  );
                addAvailableData('marker_trajectories', 'time_mocap', 'marker_labels', save_folder, save_file_name);
                disp(['Transformed marker data in ' file_name_raw ' to belt space and saved to ' file_name_shifted])                    
            end






%             % apply shift to forceplate trajectories
%             load(makeFileName(date, subject_id, condition, i_trial, 'forceplateTrajectories'));
%             belt_position_trajectory_forceplate = spline(time_belts, belt_position_trajectory_plcData, time_forceplate)';
%             
%             for i_time = 1 : length(time_forceplate)
%                 % define forceplate rotation and translation
%                 world_to_Acb_rotation = [1 0 0; 0 1 0; 0 0 1];
%                 world_to_Acb_translation = [0.5588; 0; 0];
%                 world_to_Acb_trafo = [world_to_Acb_rotation world_to_Acb_translation; 0 0 0 1];
%                 world_to_Acb_adjoint = rigidToAdjointTransformation(world_to_Acb_trafo);
% 
%                 % transform
%                 left_forceplate_wrench_Acb = (world_to_Acb_adjoint' * left_forceplate_wrench_world')';
%                 right_forceplate_wrench_Acb = (world_to_Acb_adjoint' * right_forceplate_wrench_world')';
% 
%             end
% 
%             % calculate wrenches and CoP for complete plate
%             total_forceplate_wrench_Acb = left_forceplate_wrench_Acb + right_forceplate_wrench_Acb;
%             copx_trajectory = - total_forceplate_wrench_Acb(:, 5) ./ total_forceplate_wrench_Acb(:, 3);
%             copy_trajectory = total_forceplate_wrench_Acb(:, 4) ./ total_forceplate_wrench_Acb(:, 3);
%             total_forceplate_cop_Acb = [copx_trajectory copy_trajectory];
%                 
%             % re-zero CoP for low loads
%             left_forceplate_low_load_indicator = copxl_trajectory == 0;
%             left_forceplate_cop_world(left_forceplate_low_load_indicator, :) = 0;
%             right_forceplate_low_load_indicator = copxr_trajectory == 0;
%             right_forceplate_cop_world(right_forceplate_low_load_indicator, :) = 0;
%             total_forceplate_low_load_indicator = copx_trajectory == 0;
%             total_forceplate_cop_world(total_forceplate_low_load_indicator, :) = 0;            

        end
    end
end

function data = readOpensimMot(filename, number_of_header_lines)
    % import data
    imported_data = importdata(filename, '\t', number_of_header_lines);
    data.header = imported_data.textdata(1:end-1, 1);
    data.labels = imported_data.textdata(end, 2:end);
    data.time = imported_data.data(:, 1);
    data.trajectories = imported_data.data(:, 2:end);
end

function data = writeOpensimMot(filename, path, time, data, labels, header)
    % prepare
    data_to_save = [time, data];
    number_of_columns = size(data_to_save, 2);
    
    % open file and write header
    fid = fopen([path filename], 'w');
    for i_line = 1 : length(header)
        fprintf(fid, header{i_line});
        fprintf(fid, '\n');
    end

    % write column labels
    column_headings = 'time\t';
    for i_column = 1 : length(labels)
        column_headings = [column_headings labels{i_column} '\t'];
    end
    column_headings = [column_headings '\n'];
    fprintf(fid, column_headings);

    % write data
    data_format = [];
    for i = 1:number_of_columns
        data_format = [data_format '%10.8f\t'];
    end
    data_format = [data_format '\n'];
    fprintf(fid,data_format, data_to_save');

    fclose(fid);
end
