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

function [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_posturalTransitions(trial_data)
    % allocate output variables
    conditions_trial = struct;
    event_variables_to_save = struct;
    removal_flags = 0;

    if strcmp(trial_data.trial_type, "GI")
        onset_times = trial_data.loaded_events_data.event_data{strcmp(trial_data.loaded_events_data.event_labels, 'onset')};
        heel_off_lead_times = trial_data.loaded_events_data.event_data{strcmp(trial_data.loaded_events_data.event_labels, 'heel_off_lead')};
        heel_contact_lead_times = trial_data.loaded_events_data.event_data{strcmp(trial_data.loaded_events_data.event_labels, 'heel_contact_lead')};
        heel_off_trail_times = trial_data.loaded_events_data.event_data{strcmp(trial_data.loaded_events_data.event_labels, 'heel_off_trail')};
        stretch_times = [onset_times(1) heel_off_lead_times(1) heel_contact_lead_times(1) heel_off_trail_times(1)];
        % TODO: hard-coded for now, this should depend on the preferred leg or on actual data
        stance_foot_data = {'STANCE_BOTH', 'STANCE_LEFT', 'STANCE_LEFT', 'STANCE_LEFT'};
    end

    % add new variables to be saved
    conditions_trial.condition_experimental_list = trial_data.trial_type;
    event_variables_to_save.stretch_times = stretch_times;
    event_variables_to_save.stance_foot_data = stance_foot_data;
end

