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
    
    % load settings
    subject_settings = loadSettingsFromFile('subject');
    study_settings = loadSettingsFromFile('study');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');
    
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
    % TODO: the naming is appropriate for moving points, but not for changing coordinate frames (opposite). Leave this
    % for now, but be aware that this might cause confusion
    cobal_to_opensim_scale = transformation_settings.get('cobal_to_opensim_scale');
    
    %% process
    for i_condition = 1 : length(trial_type_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            this_trial_type = trial_type_list{i_condition};
            % forceplate data
            if strcmp(type, 'forceplate') || strcmp(type, 'all')
                % load data
                forceplate_file_name = makeFileName(collection_date, subject_id, this_trial_type, i_trial, 'forceplateTrajectories.mat');
                if exist(['processed' filesep forceplate_file_name], 'file')
                    % load data
                    [forceplate_data, time_forceplate, sampling_rate_forceplate, labels_forceplate, directions_forceplate, success] = loadData(collection_date, subject_id, this_trial_type, i_trial, 'forceplate_trajectories');
                    [~, time_mocap] = loadData(collection_date, subject_id, this_trial_type, i_trial, 'marker_trajectories');
                    
                    % get indices for total forceplate
                    forceplate_table = study_settings.getTable('forceplate_table', 1);
                    total_forceplate_wrench_world_indices = ...
                      [ ...
                        find(strcmp(labels_forceplate, 'fx')) ...
                        find(strcmp(labels_forceplate, 'fy')) ...
                        find(strcmp(labels_forceplate, 'fz')) ...
                        find(strcmp(labels_forceplate, 'mx')) ...
                        find(strcmp(labels_forceplate, 'my')) ...
                        find(strcmp(labels_forceplate, 'mz')) ...
                      ];
                    total_forceplate_cop_world_indices = ...
                      [ ...
                        find(strcmp(labels_forceplate, 'copx')) ...
                        find(strcmp(labels_forceplate, 'copy')) ...
                      ];
                  
                    % extract wrench and cop
                    total_forceplate_wrench_world = forceplate_data(:, total_forceplate_wrench_world_indices);
                    total_forceplate_cop_world = forceplate_data(:, total_forceplate_cop_world_indices);
                    
                    % resample to motion capture time
                    total_wrench_world_resampled = spline(time_forceplate, total_forceplate_wrench_world', time_mocap)';
                    total_cop_world_resampled = [spline(time_forceplate, total_forceplate_cop_world', time_mocap)', zeros(size(time_mocap))];
                    
                    % transform wrenches and CoP trajectories from world to opensim coordinate frame
                    total_wrench_opensim = (cobal_to_opensim_adjoint' * total_wrench_world_resampled')';
                    total_cop_opensim = (cobal_to_opensim_rotation' * total_cop_world_resampled')';
                    
                    % invert to change from human-generated-forces to ground-reaction-forces
                    total_wrench_opensim = - total_wrench_opensim;
                  
                    % store result in output data
                    % not needed for now
                    
                    % do the same for individual forceplates
                    number_of_forceplates = size(forceplate_table, 1);
                    
                    % prepare output data
                    data_to_save = zeros(length(time_mocap), number_of_forceplates*9 + 1);
                    data_to_save = time_mocap - time_mocap(1);
                    
                    % loop through forceplates
                    for i_forceplate = 1 : number_of_forceplates
                        % get indices
                        this_forceplate_label = forceplate_table.label{i_forceplate};
                        this_forceplate_wrench_world_indices = ...
                          [ ...
                            find(strcmp(labels_forceplate, ['fx_' this_forceplate_label])) ...
                            find(strcmp(labels_forceplate, ['fy_' this_forceplate_label])) ...
                            find(strcmp(labels_forceplate, ['fz_' this_forceplate_label])) ...
                            find(strcmp(labels_forceplate, ['mx_' this_forceplate_label])) ...
                            find(strcmp(labels_forceplate, ['my_' this_forceplate_label])) ...
                            find(strcmp(labels_forceplate, ['mz_' this_forceplate_label])) ...
                          ];
                        this_forceplate_cop_world_indices = ...
                          [ ...
                            find(strcmp(labels_forceplate, ['copx_' this_forceplate_label])) ...
                            find(strcmp(labels_forceplate, ['copy_' this_forceplate_label])) ...
                          ];
                        
                        % extract wrench and cop
                        this_forceplate_wrench_world = forceplate_data(:, this_forceplate_wrench_world_indices);
                        this_forceplate_cop_world = forceplate_data(:, this_forceplate_cop_world_indices);
                        
                        % resample to motion capture time
                        this_wrench_world_resampled = spline(time_forceplate, this_forceplate_wrench_world', time_mocap)';
                        this_cop_world_resampled = [spline(time_forceplate, this_forceplate_cop_world', time_mocap)', zeros(size(time_mocap))];
                        
                        % transform wrenches and CoP trajectories from world to opensim coordinate frame
                        this_wrench_opensim = (cobal_to_opensim_adjoint' * this_wrench_world_resampled')';
                        this_cop_opensim = (cobal_to_opensim_rotation' * this_cop_world_resampled')';
                        
                        % invert to change from human-generated-forces to ground-reaction-forces
                        this_wrench_opensim = - this_wrench_opensim;
                        
                        % store result in output data
                        this_force_vector = this_wrench_opensim(:, 1:3);
                        this_force_position = this_cop_opensim;
                        this_torque_vector = this_wrench_opensim(:, 4:6);
                        
                        force_vector_indices = 1 + (i_forceplate-1)*6 + (1:3);
                        force_position_indices = 1 + (i_forceplate-1)*6 + (4:6);
                        torque_vector_indices = 1 + number_of_forceplates*6 + (i_forceplate-1)*3 + (1:3);
                        
                        data_to_save(:, force_vector_indices) = this_force_vector;
                        data_to_save(:, force_position_indices) = this_force_position;
                        data_to_save(:, torque_vector_indices) = this_torque_vector;
                        
                    end

                    % assemble column headings
                    column_headings = ['time\t'];
                    for i_forceplate = 1 : number_of_forceplates
                        this_forceplate_label = forceplate_table.label{i_forceplate};
                        
                        % add force vector and position
                        this_column_headings = ...
                          [ ...
                            [this_forceplate_label '_force_vx\t'], ...
                            [this_forceplate_label '_force_vy\t'], ...
                            [this_forceplate_label '_force_vz\t'], ...
                            [this_forceplate_label '_force_px\t'], ...
                            [this_forceplate_label '_force_py\t'], ...
                            [this_forceplate_label '_force_pz\t'], ...
                          ];
                        
                        column_headings = [column_headings this_column_headings]; %#ok<AGROW>
                    end
                    for i_forceplate = 1 : number_of_forceplates
                        this_forceplate_label = forceplate_table.label{i_forceplate};
                        
                        % add torque vector
                        this_column_headings = ...
                          [ ...
                            [this_forceplate_label '_torque_x\t'], ...
                            [this_forceplate_label '_torque_y\t'], ...
                            [this_forceplate_label '_torque_z\t'], ...
                          ];
                        
                        column_headings = [column_headings this_column_headings]; %#ok<AGROW>
                    end
                    
                    % finish with line break
                    column_headings = [column_headings '\n']; %#ok<AGROW>
                    
                    %% save .mot file -- individual force plates
%                     time_mocap = time_mocap - time_mocap(1); % first entry has to be zero or OpenSim will throw an exception

%                     data_to_save = ...
%                       [ ...
%                         time_mocap, ...
%                         wrench_left_opensim(:, 1:3), ...
%                         cop_left_opensim, ...
%                         wrench_right_opensim(:, 1:3), ...
%                         cop_right_opensim, ...
%                         wrench_left_opensim(:, 4:6), ...
%                         wrench_right_opensim(:, 4:6), ...
%                       ];

                    save_folder = ['opensim' filesep 'forceplate'];
                    save_file_name =  makeFileName(collection_date, subject_id, this_trial_type, i_trial, 'grfIndividual.mot');
                    fid = fopen([save_folder filesep save_file_name],'w');

                    number_of_data_columns = 1 + number_of_forceplates * 9;
                    fprintf(fid, 'name %s\n', save_file_name);
                    fprintf(fid, 'datacolumns %d\n', number_of_data_columns);  % total number of datacolumns
                    fprintf(fid, 'datarows %d\n', length(time_mocap)); % number of datarows
                    fprintf(fid, 'range %f %f\n', 0, time_mocap(end)-time_mocap(1)); % range of time data
                    fprintf(fid, 'endheader\n');

                    % write the new header
%                     column_headings = ...
%                       [ ...
%                         'time\t', ...
%                         'left_plate_force_vx\t', ...
%                         'left_plate_force_vy\t', ...
%                         'left_plate_force_vz\t', ...
%                         'left_plate_force_px\t', ...
%                         'left_plate_force_py\t', ...
%                         'left_plate_force_pz\t', ...
%                         'right_plate_force_vx\t', ...
%                         'right_plate_force_vy\t', ...
%                         'right_plate_force_vz\t', ...
%                         'right_plate_force_px\t', ...
%                         'right_plate_force_py\t', ...
%                         'right_plate_force_pz\t', ...
%                         'left_plate_torque_x\t', ...
%                         'left_plate_torque_y\t', ...
%                         'left_plate_torque_z\t', ...
%                         'right_plate_torque_x\t', ...
%                         'right_plate_torque_y\t', ...
%                         'right_plate_torque_z\t', ...
%                         '\n'
%                       ];
                    fprintf(fid, column_headings);

                    data_format = [];
                    for i = 1:number_of_data_columns
                        data_format = [data_format '%10.6f\t']; %#ok<AGROW>
                    end
                    data_format = [data_format '\n'];

                    fprintf(fid, data_format, data_to_save'); %#ok<AGROW>

                    fclose(fid);
                    disp(['processed ' forceplate_file_name ' and saved as ' save_file_name])
                    
                    %% save .mot file -- total forces
%                     time_mocap = time_mocap - time_mocap(1); % first entry has to be zero or OpenSim will throw an exception
% 
%                     data_to_save = ...
%                       [ ...
%                         time_mocap, ...
%                         total_wrench_opensim(:, 1:3), ...
%                         total_cop_opensim, ...
%                         total_wrench_opensim(:, 4:6) ...
%                       ];
% 
%                     save_folder = ['opensim' filesep 'forceplate'];
%                     save_file_name =  makeFileName(collection_date, subject_id, this_trial_type, i_trial, 'grfTotal.mot');
%                     fid = fopen([save_folder filesep save_file_name],'w');
% 
%                     fprintf(fid, 'name %s\n', save_file_name);
%                     fprintf(fid, 'datacolumns %d\n', 10);  % total number of datacolumns
%                     fprintf(fid, 'datarows %d\n', length(time_mocap)); % number of datarows
%                     fprintf(fid, 'range %f %f\n', time_mocap(1), time_mocap(end)); % range of time data
%                     fprintf(fid, 'endheader\n');
% 
%                     % write the new header
%                     column_headings = ...
%                       [ ...
%                         'time\t', ...
%                         'total_force_vx\t', ...
%                         'total_force_vy\t', ...
%                         'total_force_vz\t', ...
%                         'total_force_px\t', ...
%                         'total_force_py\t', ...
%                         'total_force_pz\t', ...
%                         'total_torque_x\t', ...
%                         'total_torque_y\t', ...
%                         'total_torque_z\t', ...
%                         '\n'
%                       ];
%                     fprintf(fid, column_headings);
% 
%                     data_format = [];
%                     for i = 1:28
%                         data_format = [data_format '%10.6f\t'];
%                     end
%                     data_format = [data_format '\n'];
% 
%                     fprintf(fid,data_format, data_to_save');
% 
%                     fclose(fid);
                    disp(['processed ' forceplate_file_name ' and saved as ' save_file_name])                    
                else
                    disp(['failed to load ' forceplate_file_name ' - skipped'])
                end
            end    
            
            % marker data
            if strcmp(type, 'marker') || strcmp(type, 'all')
                % load time information from associated marker file
                marker_file_name = makeFileName(collection_date, subject_id, this_trial_type, i_trial, 'markerTrajectories.mat');
                if exist(['processed' filesep marker_file_name], 'file')
                    load(['processed' filesep marker_file_name]);

                    number_of_data_points = size(marker_trajectories, 1);
                    number_of_markers = size(marker_trajectories, 2) / 3;

                    % transform to mm
                    marker_trajectories = marker_trajectories * cobal_to_opensim_scale;

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
                    save_file_name = [save_folder filesep makeFileName(collection_date, subject_id, this_trial_type, i_trial, 'marker.trc')];
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
