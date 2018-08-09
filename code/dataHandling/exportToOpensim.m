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


function exportToOpensim(varargin)
    % parse arguments
    [trial_type_list, trial_number_list, excluded_trial_type_list, excluded_trial_number_list] = parseTrialArguments(varargin{:});
    parser = inputParser;
    addParameter(parser, 'type', 'all')
    parser.KeepUnmatched = true;
    parse(parser, varargin{:})
    type = parser.Results.type;
    
    if ~directoryExists('opensim')
        mkdir('opensim')
    end
    if ~directoryExists(['opensim' filesep 'forceplate'])
        mkdir(['opensim' filesep 'forceplate'])
    end
    if ~directoryExists(['opensim' filesep 'marker'])
        mkdir(['opensim' filesep 'marker'])
    end
    if ~directoryExists(['opensim' filesep 'setupFiles'])
        mkdir(['opensim' filesep 'setupFiles'])
    end
    if ~directoryExists(['opensim' filesep 'inverseKinematics'])
        mkdir(['opensim' filesep 'inverseKinematics'])
    end
    if ~directoryExists(['opensim' filesep 'logs'])
        mkdir(['opensim' filesep 'logs'])
    end
    meter_to_millimeter = 1e3;
    
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
    
    % add static trial to list
    static_trial_type = subject_settings.get('static_reference_trial_type');
    static_trial_number = subject_settings.get('static_reference_trial_number');
    static_trial_index_in_type_list = strcmp(excluded_trial_type_list, static_trial_type);
    if ~any(static_trial_index_in_type_list)
        error(['Trial type "' static_trial_type '" specified as static reference in subjectSettings.txt, but no data found for this type.'])
    end
    if ~any(static_trial_number == excluded_trial_number_list{static_trial_index_in_type_list})
        error(['Trial number ' num2str(static_trial_number) ' of type "' static_trial_type '" specified as static reference in subjectSettings.txt, but no data found for this combination.'])
    end
    trial_type_list = [trial_type_list; static_trial_type];
    trial_number_list = [trial_number_list; static_trial_number];
    
    % define transformations
    transformation_file = [getCobalPath filesep 'resources' filesep 'opensim' filesep 'CobalOpensimTransformations.txt'];
    transformation_settings = SettingsCustodian(transformation_file);
    cobal_to_opensim_rotation = [transformation_settings.get('opensim_x_in_cobal')', transformation_settings.get('opensim_y_in_cobal')', transformation_settings.get('opensim_z_in_cobal')'];
    cobal_to_opensim_translation = transformation_settings.get('cobal_to_opensim_translation')';
%     cobal_to_opensim_rotation = [0 0 1; 1 0 0; 0 1 0];
%     cobal_to_opensim_translation = [0; 0; 0];
    cobal_to_opensim_trafo = [cobal_to_opensim_rotation cobal_to_opensim_translation; 0 0 0 1];
    cobal_to_opensim_adjoint = rigidToAdjointTransformation(cobal_to_opensim_trafo);
    
    %% process
    for i_condition = 1 : length(trial_type_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            this_trial_type = trial_type_list{i_condition};
            % forceplate data
            if strcmp(type, 'forceplate') || strcmp(type, 'all')
                % load data
                forceplate_file_name = makeFileName(subject_info.date, subject_info.subject_id, this_trial_type, i_trial, 'forceplateTrajectories.mat');
                if exist(['processed' filesep forceplate_file_name], 'file')
                    load(['processed' filesep forceplate_file_name]); %#ok<LOAD>

                    % load time information from associated marker file
                    marker_file_name = makeFileName(subject_info.date, subject_info.subject_id, this_trial_type, i_trial, 'markerTrajectories.mat');
                    load(['processed' filesep marker_file_name], 'time_mocap');

                    % resample force trajectories to mocap time
                    wrench_left_world = spline(time_forceplate, left_forceplate_wrench_world', time_mocap)';
                    wrench_right_world = spline(time_forceplate, right_forceplate_wrench_world', time_mocap)';
                    cop_left_world = [spline(time_forceplate, left_forceplate_cop_world', time_mocap)', zeros(size(time_mocap))];
                    cop_right_world = [spline(time_forceplate, right_forceplate_cop_world', time_mocap)', zeros(size(time_mocap))];

                    % transform wrenches and CoP trajectories from world to opensim coordinate frame
                    wrench_left_opensim = (cobal_to_opensim_adjoint' * wrench_left_world')';
                    wrench_right_opensim = (cobal_to_opensim_adjoint' * wrench_right_world')';

                    cop_left_opensim = (cobal_to_opensim_rotation' * cop_left_world')';
                    cop_right_opensim = (cobal_to_opensim_rotation' * cop_right_world')';

                    % invert to change from human-generated-forces to ground-reaction-forces
                    wrench_left_opensim = - wrench_left_opensim;
                    wrench_right_opensim = - wrench_right_opensim;

                    %% save .mot file
                    time_mocap = time_mocap - time_mocap(1); % first entry has to be zero or OpenSim will throw an exception

                    % write the new header
                    data_to_save = ...
                      [ ...
                        time_mocap, ...
                        wrench_left_opensim(:, 1:3), ...
                        cop_left_opensim, ...
                        wrench_right_opensim(:, 1:3), ...
                        cop_right_opensim, ...
                        wrench_left_opensim(:, 4:6), ...
                        wrench_right_opensim(:, 4:6) ...
                      ];

                    save_folder = ['opensim' filesep 'forceplate'];
                    save_file_name =  makeFileName(subject_info.date, subject_info.subject_id, this_trial_type, i_trial, 'grf.mot');
                    fid = fopen([save_folder filesep save_file_name],'w');

                    fprintf(fid, 'name %s\n', save_file_name);
                    fprintf(fid, 'datacolumns %d\n', 19);  % total number of datacolumns
                    fprintf(fid, 'datarows %d\n', length(time_mocap)); % number of datarows
                    fprintf(fid, 'range %f %f\n', time_mocap(1), time_mocap(end)); % range of time data
                    fprintf(fid, 'endheader\n');

                    column_headings = ...
                      [ ...
                        'time\t', ...
                        'ground_force_vx\t', ...
                        'ground_force_vy\t', ...
                        'ground_force_vz\t', ...
                        'ground_force_px\t', ...
                        'ground_force_py\t', ...
                        'ground_force_pz\t', ...
                        '1_ground_force_vx\t', ...
                        '1_ground_force_vy\t', ...
                        '1_ground_force_vz\t', ...
                        '1_ground_force_px\t', ...
                        '1_ground_force_py\t', ...
                        '1_ground_force_pz\t', ...
                        'ground_torque_x\t', ...
                        'ground_torque_y\t', ...
                        'ground_torque_z\t', ...
                        '1_ground_torque_x\t', ...
                        '1_ground_torque_y\t', ...
                        '1_ground_torque_z\t', ...
                        '\n'
                      ];
                    fprintf(fid, column_headings);

                    data_format = [];
                    for i = 1:19
                        data_format = [data_format '%10.6f\t'];
                    end
                    data_format = [data_format '\n'];

                    fprintf(fid,data_format, data_to_save');

                    fclose(fid);
                    disp(['processed ' forceplate_file_name ' and saved as ' save_file_name])
                else
                    disp(['failed to load ' forceplate_file_name ' - skipped'])
                end
            end    
            
            % marker data
            if strcmp(type, 'marker') || strcmp(type, 'all')
                % load time information from associated marker file
                marker_file_name = makeFileName(subject_info.date, subject_info.subject_id, this_trial_type, i_trial, 'markerTrajectories.mat');
                if exist(['processed' filesep marker_file_name], 'file')
                    load(['processed' filesep marker_file_name]);

                    number_of_data_points = size(marker_trajectories, 1);
                    number_of_markers = size(marker_trajectories, 2) / 3;

                    % transform to mm
                    marker_trajectories = marker_trajectories * meter_to_millimeter;

                    % transform to OpenSim coordinate frame
                    marker_trajectories_opensim = zeros(size(marker_trajectories));
                    for i_marker = 1 : number_of_markers
                        this_marker_trajectory_world = marker_trajectories(:, (i_marker-1)*3 + [1 2 3]);
                        this_marker_trajectory_opensim = (cobal_to_opensim_rotation' * this_marker_trajectory_world')'; 
                        % TODO: this implicitly assumes that the translation is zero, which is explicitly coded in the
                        % forceplate transformation above. If that is ever changed, it has to be changed here as well
                        marker_trajectories_opensim(:, (i_marker-1)*3 + [1 2 3]) = this_marker_trajectory_opensim;
                    end

                    % extract marker names
                    marker_names = marker_labels(1 : 3 : end);
                    for i_marker = 1 : number_of_markers
                        marker_names{i_marker} = marker_names{i_marker}(1:end-2);
                    end

                    % save
                    time_mocap = time_mocap - time_mocap(1); % first entry has to be zero or OpenSim will throw an exception
                    save_folder = ['opensim' filesep 'marker'];
                    save_file_name = [save_folder filesep makeFileName(subject_info.date, subject_info.subject_id, this_trial_type, i_trial, 'marker.trc')];
                    fid = fopen(save_file_name,'w');
                    fid_1 = fopen(save_file_name, 'w');

                    % first write the header data
                    fprintf(fid_1,'PathFileType\t4\t(X/Y/Z)\t %s\n', save_file_name);
                    fprintf(fid_1,'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
                    fprintf(fid_1,'%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n', sampling_rate_mocap, sampling_rate_mocap, number_of_data_points, number_of_markers, 'mm', sampling_rate_mocap, 1, number_of_data_points); 

                    % marker names
                    fprintf(fid_1, 'Frame#\tTime\t');
                    for i_marker = 1 : number_of_markers
                        fprintf(fid_1, [marker_names{i_marker} '\t\t\t']);
                    end
                    fprintf(fid_1, '\n');

                    % sub-header
                    fprintf(fid_1, '\t\t');
                    for i_marker = 1 : number_of_markers
                        fprintf(fid_1, ['X' num2str(i_marker) '\t' 'Y' num2str(i_marker) '\t' 'Z' num2str(i_marker) '\t']);
                    end
                    fprintf(fid_1, '\n');
                    fprintf(fid_1, '\n');

                    % then write the output marker data
                    format_text = '%i\t%2.4f\t';
                    data_out = [(1:number_of_data_points)', time_mocap, marker_trajectories_opensim];
                    for i_time = 1 : number_of_data_points
                        fprintf(fid_1, format_text, data_out(i_time, :));
                        fprintf(fid_1, '\n');
                    end

                    % close the file
                    fclose(fid_1);
                    disp(['processed ' marker_file_name ' and saved as ' save_file_name])
                else
                    disp(['failed to load ' marker_file_name ' - skipped'])
                end
            end
        end
    end
end
