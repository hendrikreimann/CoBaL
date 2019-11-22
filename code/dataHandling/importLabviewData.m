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

% input:
% Experimental data files generated by e.g. Nexus, QTM, labview etc.
% Files are expected to be in subfolders named <subject code>_<source type>, e.g. XYZ_labview
% any source type is viable, but if it is not in the default list, it must be specified as a name-value pair,
% e.g. importAscii('sources', 'someSource') would look for data in the subfolder XYZ_someSource
%
% output:
% Files containing the same data in .mat format, with some additional information about where they came from.
% Output files will be saved to folders "raw" and "processed".

function importLabviewData()


    %% prepare
    milliseconds_to_seconds = 1e-3;

    % create folders if necessary
    if ~directoryExists('raw')
        mkdir('raw')
    end
    if ~directoryExists('processed')
        mkdir('processed')
    end
    if ~directoryExists('analysis')
        mkdir('analysis')
    end

    labview_source_dir = 'labview';

    %% import data
    % get list of files to import from this directory
    clear file_name_list;
    data_dir_csv = dir([labview_source_dir filesep '*csv']);
    [file_name_list{1:length(data_dir_csv)}] = deal(data_dir_csv.name);

    % import labview saved data
    number_of_files = length(file_name_list);
    for i_file = 1 : number_of_files
        % file name stuff
        data_file_name = file_name_list{i_file};
        [date, subject_id, trial_type, trial_number, file_type, success] = getFileParameters(data_file_name);
        if success
            imported_data = importdata([labview_source_dir filesep data_file_name], ',', 2);
            labview_trajectories = imported_data.data; %#ok<NASGU>

            % extract headers
            column_name_string = imported_data.textdata{1, 1};
            labview_header = strsplit(column_name_string, ',');
            number_of_data_columns = size(imported_data.textdata, 2);

            % extract data into properly named variables
            variables_to_save = struct();
            variables_to_save_list = {};
            for i_column = 1 : number_of_data_columns
                variable_name = [strrep(labview_header{i_column}, ' ', '_'), '_trajectory'];
                extract_string = ['variables_to_save.' variable_name ' = labview_trajectories(:, i_column);'];
                eval(extract_string);
                variables_to_save_list = [variables_to_save_list; variable_name]; %#ok<AGROW>
            end

            % take special care of time, transform to seconds and rename according to file type
            variables_to_save.time = variables_to_save.time_trajectory * milliseconds_to_seconds; % transform to seconds
            variables_to_save.time = variables_to_save.time - variables_to_save.time(1); % zero
            variables_to_save = rmfield(variables_to_save, 'time_trajectory'); % this is not a variable, so remove from list
            variables_to_save_list(strcmp(variables_to_save_list, 'time_trajectory')) = []; %#ok<AGROW>
            if isfield(variables_to_save, 'QTM_time_stamp_trajectory')
                variables_to_save.qtmTimeStamp = variables_to_save.QTM_time_stamp_trajectory * milliseconds_to_seconds; % transform to seconds
                variables_to_save = rmfield(variables_to_save, 'QTM_time_stamp_trajectory'); % this is not a variable, so remove from list
                variables_to_save_list{strcmp(variables_to_save_list, 'QTM_time_stamp_trajectory')} = 'qtmTimeStamp'; %#ok<AGROW>
                first_qtm_time_stamp_first = variables_to_save.qtmTimeStamp(find(~isnan(variables_to_save.qtmTimeStamp), 1, 'first')); %#ok<NASGU>
                first_qtm_time_stamp_last = variables_to_save.qtmTimeStamp(find(~isnan(variables_to_save.qtmTimeStamp), 1, 'last')); %#ok<NASGU>
            end

            % add data source
            variables_to_save.data_source = 'labview';

            % add sampling rate
            variables_to_save.sampling_rate = NaN;

            % save
            save_folder = 'processed';
            save_file_name = makeFileName(date, subject_id, trial_type, trial_number, file_type);
            save([save_folder filesep save_file_name], '-struct', 'variables_to_save');

            for i_variable = 1 : length(variables_to_save_list)
                if ~checkDataAvailability(date, subject_id, trial_type, trial_number, variables_to_save_list{i_variable})
                    addAvailableData ...
                      ( ...
                        variables_to_save_list{i_variable}, ...
                        'time', ...
                        'sampling_rate', ...
                        variables_to_save_list{i_variable}, ...
                        {'~', '~'}, ... % placeholder for direction
                        save_folder, ...
                        save_file_name ...
                      );
                end
            end

            disp(['imported ' labview_source_dir filesep data_file_name ' and saved as ' save_folder filesep save_file_name])                
        end

    end


end





















  
