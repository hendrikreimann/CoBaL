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

% plot results

% input
% ... .mat for each subject

function calculateCohortStatistics(varargin)
    %% parse input
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'subjects', [])
    parse(parser, varargin{:})
    subjects = parser.Results.subjects;

    %% determine subjects and data folders
    [data_folder_list, subject_list] = determineDataStructure(subjects);
    if isempty(subjects)
        subjects = subject_list;
    end
    
    %% collect data from all data folders
    height_list_loaded = zeros(length(data_folder_list), 1);
    weight_list_loaded = zeros(length(data_folder_list), 1);
    age_list_loaded = zeros(length(data_folder_list), 1);
    gender_list_loaded = cell(length(data_folder_list), 1);
    subject_id_from_data_folder_list = cell(length(data_folder_list), 1);
    
    for i_folder = 1 : length(data_folder_list)
        % load data
        data_path = data_folder_list{i_folder};
        load([data_path filesep 'subjectInfo.mat'], 'height', 'weight', 'age', 'gender', 'subject_id');
        height_list_loaded(i_folder) = height;
        weight_list_loaded(i_folder) = weight;
        age_list_loaded(i_folder) = age;
        gender_list_loaded{i_folder} = gender;
        subject_id_from_data_folder_list{i_folder} = subject_id;
    end
    
    % extract data for requested subjects
    height_list = zeros(length(subjects), 1);
    weight_list = zeros(length(subjects), 1);
    age_list = zeros(length(subjects), 1);
    gender_list = cell(length(subjects), 1);
    for i_subject = 1 : length(subjects)
        data_entry_index = find(strcmp(subjects(i_subject), subject_id_from_data_folder_list), 1, 'first');
        height_list(i_subject) = height_list_loaded(data_entry_index);
        weight_list(i_subject) = weight_list_loaded(data_entry_index);
        age_list(i_subject) = age_list_loaded(data_entry_index);
        gender_list{i_subject} = gender_list_loaded{data_entry_index};
    end
    
    %% calculate stats
    mean_height = mean(height_list);
    std_height = std(height_list);
    mean_weight = mean(weight_list);
    std_weight = std(weight_list);
    mean_age = mean(age_list);
    std_age = std(age_list);
    
    % count gender
    unique_gender_entries = unique(gender_list);
    gender_count = zeros(length(unique_gender_entries), 1);
    for i_gender = 1 : length(unique_gender_entries)
        gender_count(i_gender) = sum(strcmp(unique_gender_entries{i_gender}, gender_list));
    end
    
    %% report
    disp(['height = ' num2str(mean_height) ' ' char(177) ' ' num2str(std_height)])
    disp(['weight = ' num2str(mean_weight) ' ' char(177) ' ' num2str(std_weight)])
    disp(['age    = ' num2str(mean_age) ' ' char(177) ' ' num2str(std_age)])
    gender_report_string = [num2str(gender_count(1)) ' ' unique_gender_entries{1}];
    for i_gender = 2 : length(unique_gender_entries)
        gender_report_string = [gender_report_string ', ' num2str(gender_count(i_gender)) ' ' unique_gender_entries{i_gender}];
    end
    disp(gender_report_string)
end





