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

function trial_data = determineTriggerTimes(study_settings, trial_data)
    trial_data.trigger_times = [];

    experimental_paradigm = study_settings.get('experimental_paradigm');
    paradigms_with_perturbation = ...
      { ...
        'Vision', 'CadenceVision', 'CognitiveLoadVision', 'SR_VisualStim', ...
        'GVS', 'CadenceGVS', 'FatigueGVS', 'CognitiveLoadGvs', 'OculusLaneRestriction' ...
      };
    if any(strcmp(experimental_paradigm, paradigms_with_perturbation))
        % find the time steps where the stimulus state crosses a threshold
        stimulus_threshold = 1.5;
        trigger_indices_labview = find(diff(sign(trial_data.stimulus_state_trajectory - stimulus_threshold)) > 0) + 2;
        trial_data.trigger_times = trial_data.time_stimulus(trigger_indices_labview);
    end
    
    
    if strcmp(experimental_paradigm, 'Vision_old') || strcmp(experimental_paradigm, 'GVS_old')
        % find the time steps where the stimulus state crosses a threshold
        stimulus_threshold = 0.5;
        trigger_indices_stimulus = find(diff(sign(trial_data.stimulus_state_trajectory - stimulus_threshold)) > 0) + 1;

        %
        epsilon = 1e-5;
        % remove weird noise in illusion trajectory (check labview
        % for odd behavior i.e. wait time and illusion_trajectory)
        stim_start_indices_stimulus = find(diff(sign(abs(trial_data.illusion_trajectory) - epsilon)) > 0) + 1;
        i_stim = 1;
        while i_stim ~= length(stim_start_indices_stimulus)
            if trial_data.illusion_trajectory(stim_start_indices_stimulus(i_stim) + 5) == 0
                stim_start_indices_stimulus(i_stim) = [];
                i_stim = i_stim - 1;
            end
            i_stim = i_stim + 1;
        end
        trigger_indices_stimulus = trigger_indices_stimulus(1 : length(stim_start_indices_stimulus)); % in case a stim is triggered, but not recorded

        trial_data.trigger_times = trial_data.time_stimulus(trigger_indices_stimulus);
        trial_data.stim_start_times = trial_data.time_stimulus(stim_start_indices_stimulus);
        trial_data.stim_start_indices_stimulus = stim_start_indices_stimulus;
    end
    if strcmp(experimental_paradigm, 'Stochastic Resonance')
        trial_data.trigger_times = trial_data.left_touchdown_times(1:end-1);
    end
    if strcmp(experimental_paradigm, 'GvsOverground')
        % find the time steps where the first forceplate vertical force crosses a threshold
        stimulus_threshold = 20;
%         [first_forceplate_wrench_trajectory, time_forceplate] = loadData(collection_date, subject_id, condition_list{i_condition}, i_trial, 'left_foot_wrench_world');
%         vertical_force_trajectory = first_forceplate_wrench_trajectory(:, 3);

        trigger_indices_forceplate = find(diff(sign(-trial_data.vertical_force_trajectory - stimulus_threshold)) > 0) + 2;
        trial_data.trigger_times = trial_data.time_forceplate(trigger_indices_forceplate);
    end

end