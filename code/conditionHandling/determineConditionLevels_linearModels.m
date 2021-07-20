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

function [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_linearModels(trial_data)
    bands_per_stretch = 2;

    stretch_start_times = trial_data.trigger_times;
    left_touchdown_times = trial_data.loaded_events_data.event_data{strcmp(trial_data.loaded_events_data.event_labels, 'left_touchdown')};
    right_touchdown_times = trial_data.loaded_events_data.event_data{strcmp(trial_data.loaded_events_data.event_labels, 'right_touchdown')};
    number_of_stretches = length(stretch_start_times);
    stretch_times = zeros(number_of_stretches, bands_per_stretch+1);
    first_stance_leg_list = cell(number_of_stretches, 1);
    stance_foot_data = cell(number_of_stretches, bands_per_stretch);
    removal_flags = false(number_of_stretches, 1);
    for i_stretch = 1 : number_of_stretches
        % pull out relevant events
        this_stretch_start = stretch_start_times(i_stretch);
        this_right_touchdown = min(right_touchdown_times(right_touchdown_times > this_stretch_start));
        this_left_touchdown = min(left_touchdown_times(left_touchdown_times > this_stretch_start));
        
        % determine whether this stretch is started by a left or a right heel-strike
        if ismember(this_stretch_start, left_touchdown_times)
            stance_foot_data(i_stretch, :) = {'STANCE_LEFT', 'STANCE_RIGHT'};
            this_stretch_times = [this_stretch_start this_right_touchdown this_left_touchdown];
            first_stance_leg_list{i_stretch} = 'LEFT';
        elseif ismember(this_stretch_start, right_touchdown_times)
            stance_foot_data(i_stretch, :) = {'STANCE_RIGHT', 'STANCE_LEFT'};
            this_stretch_times = [this_stretch_start this_left_touchdown this_right_touchdown];
            first_stance_leg_list{i_stretch} = 'RIGHT';
        else
            error(['stretch start time ' num2str(this_stretch_start) ' is neither a left nor a right touchdown.'])
        end
        
        % check for problems
        if ~issorted(this_stretch_times)
            removal_flags(i_stretch) = 1;
        end
        stretch_times(i_stretch, :) = this_stretch_times;
    end

    event_variables_to_save.stretch_times = stretch_times;
    event_variables_to_save.stance_foot_data = stance_foot_data;

    % conditions
    conditions_trial = struct;
    conditions_trial.first_stance_leg_list = first_stance_leg_list;
end

