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

trials_to_process_string = '12 : 20';

% figure out folder structure
subject_settings = SettingsCustodian('subjectSettings.txt');
path_split = strsplit(pwd, filesep);
subject_label = path_split{end};
study_label = path_split{end-1};

% create temporary file with commands
file_id = fopen('tmp_compute.m', 'w');
fprintf(file_id, 'addpath ~/CoBaL/helper\n');
fprintf(file_id, 'addpath ~/CoBaL/processing\n');
fprintf(file_id, 'addpath ~/miscMatlabStuff\n');
fprintf(file_id, 'addpath ~/KinematicChain\n');
fprintf(file_id, 'addpath ~/ScrewGeometry\n');
fprintf(file_id, 'parpool(32)\n');
fprintf(file_id, 'calculateKinematicTrajectories(''use_parallel'', true)\n');
fprintf(file_id, 'optimizeKinematicTrajectories(''use_parallel'', true)\n');
% fprintf(file_id, ['calculateKinematicTrajectories(''use_parallel'', true, ''trials'', ' trials_to_process_string ')\n']);
% fprintf(file_id, ['optimizeKinematicTrajectories(''use_parallel'', true, ''trials'', ' trials_to_process_string ')\n']);
fclose(file_id);

% copy command file to compute
target_dir = ['data/' study_label '/' subject_label];
system(['scp tmp_compute.m tuf79669@compute.temple.edu:' target_dir]);
delete tmp_compute.m

% % execute command file
% execution_string = ...
%   [ ...
%     '''' ...
%     'module switch matlab/R2012a matlab/R2017a' ...
%     'cd data/' study_label '/' subject_label ';' ...
%     'nohup matlab -r < tmp_compute.m -logfile mat_out.log & ' ...
%     '''' ...
%   ];
% system(['ssh tuf79669@compute.temple.edu ' execution_string]);




