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

function [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels(study_settings, subject_settings, trial_data)
    experimental_paradigm = study_settings.get('experimental_paradigm');
    
    % allocate
    conditions_trial = struct;
    event_variables_to_save = struct;
    removal_flags = zeros(size(trial_data.trigger_times));
    
    intermittent_perturbation_paradigms = {'Vision', 'CadenceVision', 'GVS', 'CadenceGVS', 'FatigueGVS', 'OculusLaneRestriction', 'CognitiveLoadVision', 'CognitiveLoadGvs', 'SR_VisualStim'};
    if any(strcmp(experimental_paradigm, intermittent_perturbation_paradigms))
        [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_intermittentPerturbations(study_settings, subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'Stochastic Resonance')
        [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_stochasticResonance(subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'Vision_old')
        [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_visionOld(study_settings, subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'GVS_old')
        [conditions_trial, event_variables_to_save, removal_flags] = determineConditionLevels_gvsOld(study_settings, subject_settings, trial_data);
    end





end