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

function conditions_trial = determineConditionLevels_affectedSide(subject_settings, trial_data, conditions_trial)
    number_of_triggers = length(trial_data.trigger_times);
    affected_side = subject_settings.get('affected_side');
    if isempty(affected_side)
        warning('"affected_side" not specified in subject settings')
    end
    if ~(strcmp(affected_side, 'left') || strcmp(affected_side, 'right')) || strcmp(affected_side, 'neither')
        warning ...
          ( ...
            [ ...
              '"affected_side: ' affected_side '" not recognized in subject settings. ' ...
              'Possible entries are "left" and "right".' ...
            ])
    end
    condition_affected_side_list = cell(number_of_triggers, 1);
    for i_stretch = 1 : number_of_triggers
        if isfield(conditions_trial, 'trigger_foot_list')
            % trigger foot is a factor, so check whether the trigger foot was the more-affected foot or not
            if strcmp(affected_side, 'left')
                if strcmp(conditions_trial.trigger_foot_list(i_stretch), 'TRIGGER_LEFT')
                    condition_affected_side_list{i_stretch} = 'TRIGGER_AFFECTED';
                elseif strcmp(conditions_trial.trigger_foot_list(i_stretch), 'TRIGGER_RIGHT')
                    condition_affected_side_list{i_stretch} = 'TRIGGER_UNAFFECTED';
                end
            elseif strcmp(affected_side, 'right')
                if strcmp(conditions_trial.trigger_foot_list(i_stretch), 'TRIGGER_RIGHT')
                    condition_affected_side_list{i_stretch} = 'TRIGGER_AFFECTED';
                elseif strcmp(conditions_trial.trigger_foot_list(i_stretch), 'TRIGGER_LEFT')
                    condition_affected_side_list{i_stretch} = 'TRIGGER_UNAFFECTED';
                end
            end
        else
            % trigger foot is not a factor, so we simply store the label for the affected side
            condition_affected_side_list{i_stretch} = affected_side;
        end
        
        
    end
    conditions_trial.affected_side_list = condition_affected_side_list;
    
end

