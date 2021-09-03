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

function importVU(varargin)
    % parse arguments
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'visualize', true)
    addParameter(parser, 'source', 'orig')
    addParameter(parser, 'file', '')
    parse(parser, varargin{:})
    settings.visualize = parser.Results.visualize;
    settings.file = parser.Results.file;

    %% prepare
    % set some parameters
    settings.millimeter_to_meter = 1e-3;

%     study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');

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
    settings.source_dir = 'VU';
    settings.collection_date = subject_settings.get('collection_date');
    settings.subject_id = subject_settings.get('subject_id');
    
    % load protocol file if available
    if exist('protocolInfo.mat', 'file')
        settings.protocol_file_available = 1;
        settings.protocol_info = load('protocolInfo.mat');
    else
        settings.protocol_file_available = 0;
    end

    %% import data
    if isempty(settings.file)
        data_dir_mat = dir([settings.source_dir filesep '*mat']);
        [file_name_list{1:length(data_dir_mat)}] = deal(data_dir_mat.name);
    else
        file_name_list = {settings.file};
    end

    % go through raw files and import
    number_of_files = length(file_name_list);
    for i_file = 1 : number_of_files
        data_file_name = file_name_list{i_file};
        importSingleFile(data_file_name, settings);
    end

end


function importSingleFile(data_file_name, settings)
    % file info
    file_info = struct;
    file_info.collection_date = settings.collection_date;
    file_info.subject_id = settings.subject_id;
    
    file_name_trunk = data_file_name(1:end-4);
    file_name_split = strsplit(file_name_trunk, '_');
    file_info.trial_number = str2double(file_name_split{2});
    
    trial_type = settings.protocol_info.trial_type(settings.protocol_info.trial_number == file_info.trial_number);
    
    if isempty(trial_type)
        disp(['did not find ' settings.source_dir, filesep, data_file_name ' in protocol table, skipping import'])
    else
        file_info.trial_type = trial_type{1};

        % load data from file
        vu_data = load([settings.source_dir, filesep, data_file_name]);

        % import data
        importTrialDataForceplate(vu_data, file_info);
        importTrialDataMarker(vu_data, file_info);

        % report to command window
        disp(['imported ' settings.source_dir, filesep, data_file_name])
    end    
end



function importTrialDataForceplate(vu_data, file_info)
    % extract data
    force_vu = vu_data.DBforces;
    total_forceplate_wrench_world = [force_vu.combined.glob.F, force_vu.combined.glob.M];
    total_forceplate_cop_world = force_vu.combined.glob.CoP(:, [1 2]);
    left_forceplate_wrench_world = [force_vu.left.glob.F, force_vu.left.glob.M];
    left_forceplate_cop_world = force_vu.left.glob.CoP(:, [1 2]);
    right_forceplate_wrench_world = [force_vu.right.glob.F, force_vu.right.glob.M];
    right_forceplate_cop_world = force_vu.right.glob.CoP(:, [1 2]);

    % labels
    left_foot_wrench_labels = {'fxl', 'fyl', 'fzl', 'mxl', 'myl', 'mzl'};
    right_foot_wrench_labels = {'fxr', 'fyr', 'fzr', 'mxr', 'myr', 'mzr'};
    total_wrench_labels = {'fx', 'fy', 'fz', 'mx', 'my', 'mz'};
    left_foot_cop_labels = {'copxl', 'copyl'};
    right_foot_cop_labels = {'copxr', 'copyr'};
    total_cop_labels = {'copx', 'copy'};

    % directions for wrench
    wrench_directions = cell(2, 6);
    [wrench_directions{1, [1 4]}] = deal('right');
    [wrench_directions{2, [1 4]}] = deal('left');
    [wrench_directions{1, [2 5]}] = deal('forward');
    [wrench_directions{2, [2 5]}] = deal('backward');
    [wrench_directions{1, [3 6]}] = deal('up');
    [wrench_directions{2, [3 6]}] = deal('down');

    % directions for cop
    cop_directions = cell(2, 2);
    [cop_directions{1, 1}] = deal('right');
    [cop_directions{2, 1}] = deal('left');
    [cop_directions{1, 2}] = deal('forward');
    [cop_directions{2, 2}] = deal('backward');

    % consolidate into single variable
    forceplate_raw_trajectories = ...
      [ ...
        total_forceplate_wrench_world, ...
        total_forceplate_cop_world, ...
        left_forceplate_wrench_world, ...
        left_forceplate_cop_world, ...
        right_forceplate_wrench_world, ...
        right_forceplate_cop_world ...
      ];
    forceplate_labels = ...
      [ ...
        total_wrench_labels, ...
        total_cop_labels, ...
        left_foot_wrench_labels, ...
        left_foot_cop_labels, ...
        right_foot_wrench_labels, ...
        right_foot_cop_labels ...
      ];
    forceplate_directions = ...
      [ ...
        wrench_directions, ...
        cop_directions, ...
        wrench_directions, ...
        cop_directions, ...
        wrench_directions, ...
        cop_directions ...
      ];

    % add auxiliary data
    sampling_rate = vu_data.DBfs;
    time = (1 : size(forceplate_raw_trajectories, 1))' / sampling_rate;

    % save forceplate data
    save_folder = 'raw';
    data_source = 'VU';
    save_file_name = makeFileName(file_info.collection_date, file_info.subject_id, file_info.trial_type, file_info.trial_number, 'forceplateTrajectoriesRaw.mat');
    save ...
        ( ...
        [save_folder filesep save_file_name], ...
        'forceplate_raw_trajectories', ...
        'forceplate_labels', ...
        'forceplate_directions', ...
        'time', ...
        'sampling_rate', ...
        'data_source' ...
        );
    addAvailableData('forceplate_raw_trajectories', 'time', 'sampling_rate', '_forceplate_labels', '_forceplate_directions', save_folder, save_file_name);
end

function importTrialDataMarker(vu_data, file_info)
    markers_vu = vu_data.traj;
    number_of_frames = size(markers_vu.segment(1).blm, 1);
    
    % prep containers and auxiliary data
    marker_raw_trajectories = [];
    marker_labels = {};
    marker_directions = {};
    sampling_rate = vu_data.traj_fs;
    time = (1 : number_of_frames)' * 1/sampling_rate;
    
    % head
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'head', 'LFHD');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'head', 'RFHD');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'head', 'LBHD');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'head', 'RBHD');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    % trunk
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'thorax', 'C7');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'thorax', 'CLAV');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    % left upper arm    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'uparm_left', 'LSHO');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'uparm_left', 'LUPA');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'uparm_left', 'LELB');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    % left lower arm    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'lowarm_left', 'LFRA');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'lowarm_left', 'LWRA');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'lowarm_left', 'LWRB');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    % right upper arm    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'uparm_right', 'RSHO');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'uparm_right', 'RUPA');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'uparm_right', 'RELB');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    % right lower arm    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'lowarm_right', 'RFRA');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'lowarm_right', 'RWRA');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'lowarm_right', 'RWRB');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    % pelvis
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'pelvis', 'LASI');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'pelvis', 'RASI');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'pelvis', 'MPSIS');
    this_marker_labels = {'LPSI_x', 'LPSI_y', 'LPSI_z'}; % use MPSI marker to get LPSI
    this_marker_trajectories(:, 1) = this_marker_trajectories(:, 1) - 0.03; % but shift 3 cm to the left
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'pelvis', 'MPSIS');
    this_marker_labels = {'RPSI_x', 'RPSI_y', 'RPSI_z'}; % use MPSI marker to get LPSI
    this_marker_trajectories(:, 1) = this_marker_trajectories(:, 1) + 0.03; % but shift 3 cm to the right
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];

    % left thigh
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'thigh_left', 'LTHIP');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'thigh_left', 'LTHI');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'thigh_left', 'LTHID');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'thigh_left', 'LKNE');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'thigh_left', 'LKNEM');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    % left shank
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'calf_left', 'LTIBP');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'calf_left', 'LTIB');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'calf_left', 'LTIBD');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'calf_left', 'LANK');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'calf_left', 'LANKM');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    % left foot
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'foot_left', 'LHEE');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'foot_left', 'LTOE');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'foot_left', 'LTOEL');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'foot_left', 'LTOETIP');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];

    % right thigh
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'thigh_right', 'RTHIP');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'thigh_right', 'RTHI');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'thigh_right', 'RTHID');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'thigh_right', 'RKNE');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'thigh_right', 'RKNEM');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    % right shank
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'calf_right', 'RTIBP');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'calf_right', 'RTIB');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'calf_right', 'RTIBD');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'calf_right', 'RANK');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'calf_right', 'RANKM');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    % right foot
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'foot_right', 'RHEE');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'foot_right', 'RTOE');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'foot_right', 'RTOEL');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    [this_marker_trajectories, this_marker_labels, this_marker_directions] = extractSingleMarkerTrajectory(markers_vu, 'foot_right', 'RTOETIP');
    marker_raw_trajectories = [marker_raw_trajectories this_marker_trajectories];
    marker_labels = [marker_labels this_marker_labels];
    marker_directions = [marker_directions this_marker_directions];
    
    % save
    save_folder = 'raw';
    save_file_name = makeFileName(file_info.collection_date, file_info.subject_id, file_info.trial_type, file_info.trial_number, 'markerTrajectoriesRaw.mat');
    save ...
        ( ...
        [save_folder filesep save_file_name], ...
        'marker_raw_trajectories', ...
        'time', ...
        'sampling_rate', ...
        'marker_labels', ...
        'marker_directions' ...
        );
    addAvailableData('marker_raw_trajectories', 'time', 'sampling_rate', '_marker_labels', '_marker_directions', save_folder, save_file_name);
end

function [marker_trajectories, marker_labels, marker_directions] = extractSingleMarkerTrajectory(data, segment_label, marker_label)
    % extract segment names
    number_of_segments = size(data.segment, 2);
    segment_names = cell(number_of_segments, 1);
    for i_segment = 1 : number_of_segments
        segment_names(i_segment) = data.segment(i_segment).name;
    end

    % extract marker data
    segment_data = data.segment(strcmp(segment_names, segment_label));
    marker_index = find(strcmp(segment_data.blm_name, marker_label));
    component_indices = (marker_index-1)*3 + (1 : 3);
    marker_trajectories = segment_data.blm(:, component_indices);
    
    % define marker labels
    marker_labels = {[marker_label '_x'], [marker_label '_y'], [marker_label '_z']};
    
    % define marker directions
    marker_directions = [{'right'; 'left'}, {'forward'; 'backward'}, {'up'; 'down'}];
    
end


















  
