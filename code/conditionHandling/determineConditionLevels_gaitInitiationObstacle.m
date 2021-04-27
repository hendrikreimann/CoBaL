%     This file is part of the CoBaL code base
%     Copyright (C) 2019 Hendrik Reimann <hendrikreimann@gmail.com>
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

function [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_gaitInitiationObstacle(study_settings, subject_settings, trial_data)
    
    % allocate output variables
    conditions_trial = struct;
    event_variables_to_save = struct;
    
    % determine start and end
    stance_foot_data = {'STANCE_BOTH', 'STANCE_BOTH', 'STANCE_LEFT'};
    removal_flags = 0;

    init_time = trial_data.right_pushoff_times(1) - 1; % assume heel-off happened at least one second before toes-off, so start looking at that point
    end_time = trial_data.right_touchdown_times(1);
    unload_time = trial_data.right_pushoff_times(1);

    % determine unload time as maximal backward-right shift of the CoP, following Halliday et al, Gait and Posture 8 (1998) 8?14
    [~, start_time_index_forceplate] = min(abs(trial_data.time_forceplate - init_time));
    [~, unload_time_index_forceplate] = min(abs(trial_data.time_forceplate - unload_time));
    cop_data_relevant = trial_data.cop_world_trajectory(start_time_index_forceplate : unload_time_index_forceplate, :);
    time_forceplate_relevant = trial_data.time_forceplate(start_time_index_forceplate : unload_time_index_forceplate);
    [~, release_time_index_forceplate] = max(cop_data_relevant(:, 1));
    release_time = time_forceplate_relevant(release_time_index_forceplate);
    start_time = release_time - 0.5;


    stretch_start_times = trial_data.right_pushoff_times(1) - 1; % HR: this is probably not right anymore
    stretch_end_times = trial_data.right_touchdown_times(1);
    stretch_pushoff_times = 0;
    condition_experimental_list = {trial_data.condition_experimental};
    stretch_times = [start_time release_time unload_time end_time];

%     if visualize
%         for i_trigger = 1 : length(stretch_start_times)
%             if strcmp(stance_foot_data(i_trigger), 'STANCE_LEFT')
%                 stretch_indicator_height = 0.01;
%             else
%                 stretch_indicator_height = -0.01;
%             end
% 
%             plot([stretch_start_times(i_trigger) stretch_end_times(i_trigger)], [1 1]*stretch_indicator_height, 'linewidth', 3);
%         end
%     end

    % add new variables to be saved
    conditions_trial.condition_experimental_list = condition_experimental_list;
    event_variables_to_save.stretch_pushoff_times = stretch_pushoff_times;
    event_variables_to_save.stretch_times = stretch_times;
    event_variables_to_save.stance_foot_data = stance_foot_data;




    
    
end

