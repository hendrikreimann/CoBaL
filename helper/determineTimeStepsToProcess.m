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

function time_steps_to_optimize = determineTimeStepsToProcess(date, subject_id, condition, trial_number, padding)
    load(['processed' filesep makeFileName(date, subject_id, condition, trial_number, 'markerTrajectories')]);
    load(['analysis' filesep makeFileName(date, subject_id, condition, trial_number, 'relevantDataStretches')]);
    
    time_steps_to_optimize_indicator = false(length(time_mocap), 1);
    for i_stretch = 1 : length(stretch_start_times)
        % get unpadded stretches as indices
        this_stretch_start_time = stretch_start_times(i_stretch);
        this_stretch_end_time = stretch_end_times(i_stretch);
        [~, this_stretch_start_index_mocap] = min(abs(this_stretch_start_time - time_mocap));
        [~, this_stretch_end_index_mocap] = min(abs(this_stretch_end_time - time_mocap));
        
        % add padding 
        this_stretch_start_index_mocap_padded = max([this_stretch_start_index_mocap - padding, 1]);
        this_stretch_end_index_mocap_padded = min([this_stretch_end_index_mocap + padding, length(time_mocap)]);
        
        % flip indicators
        time_steps_to_optimize_indicator(this_stretch_start_index_mocap_padded : this_stretch_end_index_mocap_padded) = true;
    end
    
    time_steps_to_optimize = find(time_steps_to_optimize_indicator);
    
    if iscolumn(time_steps_to_optimize)
        time_steps_to_optimize = time_steps_to_optimize';
    end
end