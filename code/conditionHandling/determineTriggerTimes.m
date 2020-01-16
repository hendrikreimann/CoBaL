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

function trigger_times = determineTriggerTimes(study_settings, trial_data)
    trigger_times = [];

    experimental_paradigm = study_settings.get('experimental_paradigm');
    perturbation_paradigms = ...
      { ...
        'Vision', 'CadenceVision', 'CognitiveLoadVision', 'SR_VisualStim', ...
        'GVS', 'CadenceGVS', 'FatigueGVS', 'CognitiveLoadGvs', 'OculusLaneRestriction' ...
      };
    if any(strcmp(experimental_paradigm, perturbation_paradigms))
        % find the time steps where the stimulus state crosses a threshold
        stimulus_threshold = 1.5;
        trigger_indices_labview = find(diff(sign(trial_data.stimulus_state_trajectory - stimulus_threshold)) > 0) + 2;
        trigger_times = trial_data.time_stimulus(trigger_indices_labview);
    end
    
    if strcmp(experimental_paradigm, 'Stochastic Resonance')
        trigger_times = trial_data.left_touchdown_times(1:end-1);
    end

end