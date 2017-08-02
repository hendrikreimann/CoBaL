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
remote = 'jhrei@128.175.143.182';

% figure out folder structure
subject_settings = SettingsCustodian('subjectSettings.txt');
path_split = strsplit(pwd, filesep);
subject_label = path_split{end};
study_label = path_split{end-1};

% copy data
source_dir = ['data/' study_label '/' subject_label '/processed/'];
% system(['scp ' remote ':' source_dir '*kinematicTrajectories.mat processed/']);
system(['scp ' remote ':' source_dir '*dynamicTrajectories.mat processed/']);
source_dir = ['data/' study_label '/' subject_label '/analysis/'];
system(['scp ' remote ':' source_dir '*.mat analysis/']);





