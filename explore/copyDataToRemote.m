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

% set remote
% remote = 'tuf79669@compute.temple.edu';
remote = 'jhrei@128.175.94.240';


% figure out folder structure and create folders on remote
subject_settings = SettingsCustodian('subjectSettings.txt');
path_split = strsplit(pwd, filesep);
subject_label = path_split{end};
study_label = path_split{end-1};

% create folders
folder_creation_string = ...
  [ ...
    '''' ...
    'mkdir data;' ...
    'mkdir data/' study_label ';' ...
    'mkdir data/' study_label '/' subject_label ';' ...
    'mkdir data/' study_label '/' subject_label '/processed;' ...
    'mkdir data/' study_label '/' subject_label '/analysis;' ...
    '''' ...
  ];
system(['ssh ' remote ' ' folder_creation_string]);

% copy study settings
target_dir = ['data/' study_label];
system(['scp ../studySettings.txt ' remote ':' target_dir]);


% copy subject settings, info and model files
target_dir = ['data/' study_label '/' subject_label];
system(['scp subjectSettings.txt ' remote ':' target_dir]);
system(['scp *.mat ' remote ':' target_dir]);

% copy data
target_dir = ['data/' study_label '/' subject_label '/processed'];
system(['scp processed/*kinematicTrajectories.mat ' remote ':' target_dir]);
system(['scp processed/*markerTrajectories.mat ' remote ':' target_dir]);
% system(['scp processed/*forceplateTrajectories.mat ' remote ':' target_dir]);
% system(['scp processed/*labviewData.mat ' remote ':' target_dir]);
% system(['scp processed/*.mat ' remote ':' target_dir]);

target_dir = ['data/' study_label '/' subject_label '/analysis'];
system(['scp analysis/*.mat ' remote ':' target_dir]);


processOnRemote


