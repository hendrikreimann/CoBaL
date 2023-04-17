%     This file is part of the CoBaL code base
%     Copyright (C) 2022 Hendrik Reimann <hendrikreimann@gmail.com>
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

function createVirtualPSISMarkersVU()

    % load settings
    subject_settings = SettingsCustodian('subjectSettings.txt');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id');

    % load static reference file
    static_reference_file_name = ['processed' filesep makeFileName(collection_date, subject_id, subject_settings.get('static_reference_trial_type'), subject_settings.get('static_reference_trial_number'), 'markerTrajectories')];
    static_reference_data = load(static_reference_file_name);

    % find first time step where all markers are available
    i_time = 1;
    while any(isnan(static_reference_data.marker_trajectories(i_time, :)))
        i_time = i_time + 1;
    end
    marker_reference = static_reference_data.marker_trajectories(i_time, :);
    
    % extract relevant data and indices
    LASI_reference = extractMarkerData(marker_reference, static_reference_data.marker_labels, 'LASI');
    LASI_indices = extractMarkerData(marker_reference, static_reference_data.marker_labels, 'LASI', 'indices');
    RASI_reference = extractMarkerData(marker_reference, static_reference_data.marker_labels, 'RASI');
    RASI_indices = extractMarkerData(marker_reference, static_reference_data.marker_labels, 'RASI', 'indices');
    LPSI_indices = extractMarkerData(marker_reference, static_reference_data.marker_labels, 'LPSI', 'indices');
    RPSI_indices = extractMarkerData(marker_reference, static_reference_data.marker_labels, 'RPSI', 'indices');
    MPSIS_reference = extractMarkerData(marker_reference, static_reference_data.marker_labels, 'MPSIS');
    
    % replace RPSI and LPSI in reference with shifted MPSI
    shift = 0.03;
    LPSI_reference_new = MPSIS_reference + [-shift 0 0];
    RPSI_reference_new = MPSIS_reference + [shift 0 0];
    marker_reference(LPSI_indices) = LPSI_reference_new;
    marker_reference(RPSI_indices) = RPSI_reference_new;
    
    % rigid body fill the RPSI and LPSI trajectories for the static trial
    visualize = 0;
    static_reference_trajectories_new_tmp = rigidBodyFillMarkerFromReference ...
      ( ...
        static_reference_data.marker_trajectories, ...
        static_reference_data.marker_labels, ...
        static_reference_data.marker_directions, ...
        marker_reference, ...
        static_reference_data.marker_labels, ...
        static_reference_data.marker_directions, ...
        'LPSI', ...
        'MPSIS', ...
        'RASI', ...
        'LASI', ...
        visualize ...
      );
    static_reference_trajectories_new = rigidBodyFillMarkerFromReference ...
      ( ...
        static_reference_trajectories_new_tmp, ...
        static_reference_data.marker_labels, ...
        static_reference_data.marker_directions, ...
        marker_reference, ...
        static_reference_data.marker_labels, ...
        static_reference_data.marker_directions, ...
        'RPSI', ...
        'MPSIS', ...
        'RASI', ...
        'LASI', ...
        visualize ...
      );
    
    % save static reference data
    static_reference_data.marker_trajectories = static_reference_trajectories_new;
    save(static_reference_file_name, '-struct', 'static_reference_data');
    
    % fill everything else
    rigidBodyFillTrialFromReference('LPSI', 'MPSIS', 'RASI', 'LASI');
    rigidBodyFillTrialFromReference('RPSI', 'MPSIS', 'RASI', 'LASI');
end