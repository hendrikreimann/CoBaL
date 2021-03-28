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

% this script transform raw data from ascii into matlab format

function importFromOpensim(varargin)
    % parse arguments
    parser = inputParser;
    addParameter(parser, 'type', 'all')
    sources_default = {['opensim' filesep 'inverseKinematics'], ['opensim' filesep 'bodyKinematics']};
%     sources_default = {['opensim' filesep 'bodyKinematics']};
    addParameter(parser, 'sources', sources_default)
    parser.KeepUnmatched = true;
    parse(parser, varargin{:})
    sources = parser.Results.sources;
    
    
    direction_file = [getCobalPath filesep 'resources' filesep 'opensim' filesep 'OpensimDirections.txt'];
    direction_settings = SettingsCustodian(direction_file);
    
    % define transformations
    transformation_file = [getCobalPath filesep 'resources' filesep 'opensim' filesep 'CobalOpensimTransformations.txt'];
    transformation_settings = SettingsCustodian(transformation_file);
    cobal_to_opensim_rotation = [transformation_settings.get('opensim_x_in_cobal')', transformation_settings.get('opensim_y_in_cobal')', transformation_settings.get('opensim_z_in_cobal')'];
    cobal_to_opensim_translation = transformation_settings.get('cobal_to_opensim_translation')';
    opensim_to_cobal_rotation = cobal_to_opensim_rotation';
    opensim_to_cobal_translation = - cobal_to_opensim_translation;
    
    cobal_to_opensim_trafo = [cobal_to_opensim_rotation cobal_to_opensim_translation; 0 0 0 1];
    opensim_to_cobal_trafo = [opensim_to_cobal_rotation opensim_to_cobal_translation; 0 0 0 1];
    opensim_to_cobal_adjoint = rigidToAdjointTransformation(opensim_to_cobal_trafo);
    % TODO: the naming is appropriate for moving points, but not for changing coordinate frames (opposite). Leave this
    % for now, but be aware that this might cause confusion
    
    for i_source = 1 : length(sources)
        source_dir = sources{i_source};
        if exist(source_dir, 'dir')
            % get list of files to import from this directory
            clear file_name_list_mot;
            clear file_name_list_sto;
            data_dir_mot = dir([source_dir filesep '*.mot']);
            [file_name_list_mot{1:length(data_dir_mot)}] = deal(data_dir_mot.name);
            data_dir_sto = dir([source_dir filesep '*.sto']);
            [file_name_list_sto{1:length(data_dir_sto)}] = deal(data_dir_sto.name);
            file_name_list = [file_name_list_mot; file_name_list_sto];

            % go through files and import
            number_of_files = length(file_name_list);
            for i_file = 1 : number_of_files
                % file name stuff
                data_file_name = file_name_list{i_file};
                
                if strcmp(data_file_name(end-3:end), '.mot')
                    [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(data_file_name);
                    number_of_header_lines = 11;

                    % import data
                    [imported_data, delimiter, number_of_header_lines_returned] = importdata([source_dir filesep data_file_name], '\t', number_of_header_lines);
                    data_headers = imported_data.textdata(9, :);
                    data_type = reformatDataType(imported_data.textdata{1, 1});

                    time = imported_data.data(:, 1);
                    sampling_rate = 1/median(diff(time));

                    data_trajectories_opensim = imported_data.data(:, 2:end);
                    labels_opensim = data_headers(2:end);
                    directions_opensim = createDirectionsForOpensimData(labels_opensim, direction_settings);
                    [data_trajectories, labels, directions] = transformSpatialData(data_trajectories_opensim, labels_opensim, directions_opensim, cobal_to_opensim_trafo);

                    % save
                    variables_to_save = struct;
                    variables_to_save.(data_type) = data_trajectories;
                    variables_to_save.time = time;
                    variables_to_save.sampling_rate = sampling_rate;
                    variables_to_save.labels = labels;
                    variables_to_save.directions = directions; %#ok<STRNU>
                    
                    save_folder = 'processed';
                    save_file_name = makeFileName(date, subject_id, trial_type, trial_number, [file_type '.mat']);
                    save([save_folder filesep save_file_name], '-struct', 'variables_to_save');
                    addAvailableData ...
                      ( ...
                        data_type, ...
                        'time', ...
                        'sampling_rate', ...
                        '_labels', ...
                        '_directions', ...
                        save_folder, ...
                        save_file_name ...
                      );


                    disp(['imported ' source_dir filesep data_file_name])
                end
                
                if strcmp(data_file_name(end-3:end), '.sto')
                        [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(data_file_name);
                    number_of_header_lines = 19;

                    % import data
                    [imported_data, delimiter, number_of_header_lines_returned] = importdata([source_dir filesep data_file_name], '\t', number_of_header_lines);
                    data_headers = imported_data.textdata(14, :);
                    data_type = reformatDataType(imported_data.textdata{1, 1});

                    time = imported_data.data(:, 1);
                    sampling_rate = 1/median(diff(time));

                    data_trajectories_opensim = imported_data.data(:, 2:end);
                    labels_opensim = data_headers(2:end);
                    directions_opensim = createDirectionsForOpensimData(labels_opensim, direction_settings);
                    [data_trajectories, labels, directions] = transformSpatialData(data_trajectories_opensim, labels_opensim, directions_opensim, cobal_to_opensim_trafo);
                    
                    % save
                    variables_to_save = struct;
                    variables_to_save.(data_type) = data_trajectories;
                    variables_to_save.time = time;
                    variables_to_save.sampling_rate = sampling_rate;
                    variables_to_save.labels = labels;
                    variables_to_save.directions = directions;
                    
                    save_folder = 'processed';
                    save_file_name = makeFileName(date, subject_id, trial_type, trial_number, [file_type '.mat']);
                    save([save_folder filesep save_file_name], '-struct', 'variables_to_save');
                    addAvailableData ...
                      ( ...
                        data_type, ...
                        'time', ...
                        'sampling_rate', ...
                        '_labels', ...
                        '_directions', ...
                        save_folder, ...
                        save_file_name ...
                      );
                  
                    disp(['imported ' source_dir filesep data_file_name])
                end
                
                
            end

            disp(['imported ' num2str(number_of_files) ' files'])        

        else
            error(['Failed to locate directory ' source_dir ' to process .mot files'])



        end
    end

end


function variable_name_out = reformatDataType(data_type_in)
    variable_name_out = 'DataNotSpecified';
    if strcmp(data_type_in, 'Coordinates')
        variable_name_out = 'joint_angle_trajectories';
    end
    if strcmp(data_type_in, 'Positions')
        variable_name_out = 'com_position_trajectories';
    end
    if strcmp(data_type_in, 'Velocities')
        variable_name_out = 'com_velocity_trajectories';
    end
    if strcmp(data_type_in, 'Accelerations')
        variable_name_out = 'com_acceleration_trajectories';
    end

end












