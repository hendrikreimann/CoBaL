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

function [data_folder_list, subject_list] = determineDataStructure(subjects)
    % determine subject list and current folder type
    
    current_file_path = pwd;
    
    if exist([current_file_path filesep 'subjectSettings.txt'], 'file')
        % data folders contain subjectSettings.txt files
        current_folder_type = 'data';
    else
        % we're not in a data folder
        if exist([current_file_path filesep 'studySettings.txt'], 'file')
            % study folders contain studySettings.txt files
            current_folder_type = 'study';
        else
            % we're neither in a data folder nor in a study folder
            if exist([current_file_path filesep '..' filesep 'studySettings.txt'], 'file')
                % but one level above is a study folder, so we must be in a subject folder containing data folders
                current_folder_type = 'subject';
            else
                error('Something is wrong with the folder structure. Current folder does not appear to be a study, subject or data folder.')
            end
        end
    end
    
    % determine data folders
    if strcmp(current_folder_type, 'data')
        path_split = strsplit(pwd, filesep);
        subject_list = path_split(end-1);
        data_folder_list = {pwd};
    end
    if strcmp(current_folder_type, 'subject')
        % get list of everything in the current directory
        things_in_current_folder = dir;
        % get a logical vector that tells which is a directory
        dir_flags = [things_in_current_folder.isdir];
        % extract only those that are directories.
        dir_name = {things_in_current_folder.name};
        folder_list = dir_name(dir_flags);
        % remove folders called 'figures'
        folder_list(strcmp(folder_list, 'figures')) = [];
        
        % store data folders in lists
        number_of_data_folders = length(folder_list);
        path_split = strsplit(pwd, filesep);
        subject = path_split{end};
        subject_list = cell(number_of_data_folders, 1);
        data_folder_list = cell(number_of_data_folders, 1);
        for i_folder = 1 : number_of_data_folders
            subject_list{i_folder} = subject;
            data_folder_list{i_folder} = [pwd filesep folder_list{i_folder}];
        end
    end
    if strcmp(current_folder_type, 'study')
        % no list was passed, so go through folders
        
        % get list of everything in the current directory
        things_in_current_folder = dir;
        % get a logical vector that tells which is a directory
        dir_flags = [things_in_current_folder.isdir];
        % extract only those that are directories.
        dir_name = {things_in_current_folder.name};
        folder_list = dir_name(dir_flags);
        
        % include only folders that were specified in subject argument
        if ~isempty(subjects)
            similar_folders = find(ismember(folder_list,subjects));
            folder_list = folder_list(similar_folders);
        end
        
        number_of_data_folders = length(folder_list);
        subjects = {};
        for i_folder = 1 : number_of_data_folders
            this_folder_name = folder_list{i_folder};
            if exist([current_file_path filesep this_folder_name filesep 'subjectSettings.txt'], 'file')
                subjects = [subjects; this_folder_name];
            end
        end
        
        % check each subject folder for data folders
        subject_list = {};
        data_folder_list = {};
        for i_subject = 1 : length(subjects)
            subject_path = [pwd filesep subjects{i_subject}];
            if exist([subject_path filesep 'subjectSettings.txt'], 'file')
                % subject folder is also a data folder
                subject_list = [subject_list; subjects{i_subject}];
                data_folder_list = [data_folder_list; subject_path];
            else
                % subject folder contains data folders
                things_in_subject_folder = dir(subject_path);
                % get a logical vector that tells which is a directory
                dir_flags = [things_in_subject_folder.isdir];
                % extract only those that are directories.
                dir_name = {things_in_subject_folder.name};
                folder_list = dir_name(dir_flags);
                % remove pointers to upper level directories
                folder_list(1:2) = [];
                
                % remove figure folder within subject folder
                if any(strcmp(folder_list, 'figures'))
                    figure_folder_index = strcmp(folder_list, 'figures');
                    folder_list(figure_folder_index) = [];
                end

                % store data folders in lists
                number_of_data_folders = length(folder_list);
                path_split = strsplit(subject_path, filesep);
                subject = path_split{end};
%                 subject_list = cell(number_of_data_folders, 1);
%                 data_folder_list = cell(number_of_data_folders, 1);
                for i_folder = 1 : number_of_data_folders
                    this_subject_list{i_folder} = subject;
                    this_data_folder_list{i_folder} = [subject_path filesep folder_list{i_folder}];
                end              
                subject_list = [subject_list; this_subject_list'];
                data_folder_list = [data_folder_list; this_data_folder_list'];
            end
        end
    end


end
