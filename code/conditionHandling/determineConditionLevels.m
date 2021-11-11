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

function [conditions_trial, event_variables_to_save, removal_flags] ...
    = determineConditionLevels(study_settings, subject_settings, trial_data)
    experimental_paradigm = study_settings.get('experimental_paradigm');
    subject_id = subject_settings.get('subject_id');
    gender = subject_settings.get('gender');
    conditions_table = study_settings.get('conditions');
    
    % allocate
    conditions_trial = struct;
    event_variables_to_save = struct;
    event_variables_to_save.stretch_times = [];
    removal_flags = zeros(size(trial_data.trigger_times));
    
    % determine levels particular to the experimental paradigm
    paradigms_with_intermittent_perturbation = ...
      { ...
        'Vision', 'CadenceVision', 'GVS', 'CadenceGVS', 'FatigueGVS', ...
        'OculusLaneRestriction', 'CognitiveLoadVision', 'CognitiveLoadGvs', 'SR_VisualStim' ...
      };
    if any(strcmp(experimental_paradigm, paradigms_with_intermittent_perturbation))
        [conditions_trial, event_variables_to_save, removal_flags] ...
            = determineConditionLevels_intermittentPerturbations(study_settings, subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'Normal Walking')
        [conditions_trial, event_variables_to_save, removal_flags] ...
            = determineConditionLevels_normalWalking(trial_data);
    end
    if strcmp(experimental_paradigm, 'Linear Models')
        [conditions_trial, event_variables_to_save, removal_flags] ...
            = determineConditionLevels_linearModels(trial_data, study_settings);
    end
    if strcmp(experimental_paradigm, 'Stochastic Resonance')
        [conditions_trial, event_variables_to_save, removal_flags] ...
            = determineConditionLevels_stochasticResonance(subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'Normal Walking nGVS')
        [conditions_trial_normal, event_variables_to_save_normal, removal_flags_normal] ...
            = determineConditionLevels_normalWalking(trial_data);
        [conditions_trial_ngvs, event_variables_to_save_ngvs, removal_flags_ngvs] ...
            = determineConditionLevels_Ngvs(subject_settings, trial_data);
        
        conditions_trial = mergeConditionStruct(conditions_trial_normal, conditions_trial_ngvs);
        event_variables_to_save = mergeConditionStruct(event_variables_to_save_normal, event_variables_to_save_ngvs);
        removal_flags = removal_flags_normal & removal_flags_ngvs;
    end
    if strcmp(experimental_paradigm, 'Self Pacing Comparison')
        [conditions_trial_normal, event_variables_to_save_normal, removal_flags_normal] ...
            = determineConditionLevels_normalWalking(trial_data);
        [conditions_trial_selfpacing, event_variables_to_save_selfpacing, removal_flags_selfpacing] ...
            = determineConditionLevels_selfpacing(subject_settings, trial_data);
        
        conditions_trial = mergeConditionStruct(conditions_trial_normal, conditions_trial_selfpacing);
        event_variables_to_save = mergeConditionStruct(event_variables_to_save_normal, event_variables_to_save_selfpacing);
        removal_flags = removal_flags_normal & removal_flags_selfpacing;
    end
    if strcmp(experimental_paradigm, 'nGVS_Vision')
        [conditions_trial_vision, event_variables_to_save_vision, removal_flags_vision] ...
            = determineConditionLevels_intermittentPerturbations(study_settings, subject_settings, trial_data);
        [conditions_trial_ngvs, event_variables_to_save_ngvs, removal_flags_ngvs] ...
            = determineConditionLevels_Ngvs(subject_settings, trial_data);
        
        conditions_trial = mergeConditionStruct(conditions_trial_vision, conditions_trial_ngvs);
        event_variables_to_save = mergeConditionStruct(event_variables_to_save_vision, event_variables_to_save_ngvs);
        removal_flags = removal_flags_vision & removal_flags_ngvs;
    end
    if strcmp(experimental_paradigm, 'Vision_old')
        [conditions_trial, event_variables_to_save, removal_flags] ...
            = determineConditionLevels_visionOld(study_settings, subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'GVS_old')
        [conditions_trial, event_variables_to_save, removal_flags] ...
            = determineConditionLevels_gvsOld(study_settings, subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'GvsOverground')
        [conditions_trial, event_variables_to_save, removal_flags] ...
            = determineConditionLevels_gvsOverground(study_settings, subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'GaitInitiationObstacle')
        [conditions_trial, event_variables_to_save, removal_flags] ...
            = determineConditionLevels_gaitInitiationObstacle(study_settings, subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'Vision Stochastic')
        [conditions_trial, event_variables_to_save, removal_flags] ...
            = determineConditionLevels_visionStochastic(study_settings, subject_settings, trial_data);
    end
    if strcmp(experimental_paradigm, 'APDM')
        [conditions_trial, event_variables_to_save, removal_flags] ...
            = determineConditionLevels_normalWalking(trial_data);
    end
    if strcmp(experimental_paradigm, 'postural transitions')
        [conditions_trial, event_variables_to_save, removal_flags] ...
            = determineConditionLevels_posturalTransitions(trial_data);
    end
    
    % add affected_side if required
    if any(strcmp(conditions_table(:, 1), 'affected_side'))
        conditions_trial ...
            = determineConditionLevels_affectedSide(subject_settings, trial_data, conditions_trial);
    end
    
    % add group if required
    if any(strcmp(conditions_table(:, 1), 'group'))
        conditions_trial ...
            = determineConditionLevels_group(subject_settings, trial_data, conditions_trial);
    end
    
    % add block if required
    if any(strcmp(conditions_table(:, 1), 'block'))
        conditions_trial ...
            = determineConditionLevels_block(subject_settings, trial_data, conditions_trial);
    end
    
    % add type from filename if required
    if any(strcmp(conditions_table(:, 1), 'type'))
        conditions_trial ...
            = determineConditionLevels_type(trial_data, conditions_trial);
    end
    
    % add metronome from metronome filename if required
    if any(strcmp(conditions_table(:, 1), 'metronome'))
        conditions_trial ...
            = determineConditionLevels_metronome(trial_data, conditions_trial);
    end
    
    % add subject
    if any(strcmp(conditions_table(:, 1), 'subject'))
        condition_subject_list = cell(size(event_variables_to_save.stretch_times, 1), 1);
        for i_stretch = 1 : length(condition_subject_list)
            condition_subject_list{i_stretch} = subject_id;
        end
        conditions_trial.subject_list = condition_subject_list;
    end
    
    % add gender
    if any(strcmp(conditions_table(:, 1), 'gender'))
        condition_gender_list = cell(size(event_variables_to_save.stretch_times, 1), 1);
        for i_stretch = 1 : length(condition_gender_list)
            condition_gender_list{i_stretch} = gender;
        end
        conditions_trial.gender_list = condition_gender_list;
    end
end

function merged_struct = mergeConditionStruct(struct_one, struct_two)
    merged_struct = struct_one;
    names = fieldnames(struct_two);
    for i_field = 1 : length(names)
        this_field_name = names{i_field};
        this_field_data = struct_two.(this_field_name);
        
        if isfield(struct_one, this_field_name)
            % field is already present in 
            if ~isequal(struct_one.(this_field_name), struct_two.(this_field_name))
                error(['Field "' this_field_name '" already present with different data'])
            end
        end
        
        merged_struct.(this_field_name) = this_field_data;
    end

end
