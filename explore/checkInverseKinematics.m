%     This file is part of the CoBaL code base
%     Copyright (C) 2020 Hendrik Reimann <hendrikreimann@gmail.com>
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

% specify trial to run
trial_type = 'quietEC';
trial_number = 4;
marker_to_compare_label = 'R_ear';

% trial_type = 'ECFoamCog';
% trial_number = 1;
% marker_to_compare_label = 'LEYE';

components_to_compare = [1];

% load data
subject_info = load('subjectInfo.mat');
model_data = load('subjectModel.mat');
model = model_data.kinematic_tree;
marker_data = load(['processed' filesep makeFileName(subject_info.date, subject_info.subject_id, trial_type, trial_number, 'markerTrajectories')]);
angle_data = load(['processed' filesep makeFileName(subject_info.date, subject_info.subject_id, trial_type, trial_number, 'kinematicTrajectories')]);
marker_to_compare_trajectory_measured = extractMarkerData(marker_data.marker_trajectories, marker_data.marker_labels, marker_to_compare_label);

% run forward kinematics
number_of_time_steps = length(marker_data.time_mocap);
marker_indices = extractMarkerData(model.markerReferencePositions, model.markerLabels, marker_to_compare_label, 'indices');
marker_to_compare_trajectory_from_model = zeros(number_of_time_steps, 3);
for i_time = 1 : number_of_time_steps
    joint_angles_here = angle_data.joint_angle_trajectories(i_time, :)';
    model.jointAngles = joint_angles_here;
    model.updateConfiguration();
    marker_positions_here = model.exportMarkerPositions();
    marker_to_compare_position_here = marker_positions_here(marker_indices);
    marker_to_compare_trajectory_from_model(i_time, :) = marker_to_compare_position_here;
end

% visualize
figure;
number_of_components = length(components_to_compare);
for i_component = 1 : length(components_to_compare)
    subplot(number_of_components, 1, i_component); 
    title([subject_info.subject_id ', ' trial_type zeroPrefixedIntegerString(trial_number, 3) ', ' marker_to_compare_label ' marker, component ' num2str(components_to_compare(i_component)) ' (mean free)']); hold on;
    plot(marker_data.time_mocap, marker_to_compare_trajectory_from_model(:, components_to_compare(i_component)) - mean(marker_to_compare_trajectory_from_model(:, components_to_compare(i_component))));
    plot(marker_data.time_mocap, marker_to_compare_trajectory_measured(:, components_to_compare(i_component)) - mean(marker_to_compare_trajectory_measured(:, components_to_compare(i_component))));
    legend('reconstructed from model', 'measured')
    xlabel('time (s)'); ylabel('position (m)');
end








