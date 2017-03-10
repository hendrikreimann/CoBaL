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

% figure out folder structure and create folders on remote
subject_settings = SettingsCustodian('subjectSettings.txt');
path_split = strsplit(pwd, filesep);
subject_label = path_split{end};
study_label = path_split{end-1};

% folder_creation_string = ...
%   [ ...
%     '''' ...
%     'mkdir data;' ...
%     'mkdir data/' study_label ';' ...
%     'mkdir data/' study_label '/' subject_label ';' ...
%     'mkdir data/' study_label '/' subject_label '/processed;' ...
%     'mkdir data/' study_label '/' subject_label '/analysis;' ...
%     '''' ...
%   ];
% 
% system(['ssh tuf79669@compute.temple.edu ' folder_creation_string]);
% 
% % copy data
% target_dir = ['data/' study_label];
% system(['scp ../studySettings.txt tuf79669@compute.temple.edu:' target_dir]);
%  
% target_dir = ['data/' study_label '/' subject_label];
% system(['scp subjectSettings.txt tuf79669@compute.temple.edu:' target_dir]);
% system(['scp subjectInfo.mat tuf79669@compute.temple.edu:' target_dir]);
% system(['scp subjectModel.mat tuf79669@compute.temple.edu:' target_dir]);

target_dir = ['data/' study_label '/' subject_label '/processed'];
system(['scp processed/*kinematicTrajectories.mat tuf79669@compute.temple.edu:' target_dir]);
system(['scp processed/*markerTrajectories.mat tuf79669@compute.temple.edu:' target_dir]);

target_dir = ['data/' study_label '/' subject_label '/analysis'];
system(['scp analysis/*availableVariables.mat tuf79669@compute.temple.edu:' target_dir]);





