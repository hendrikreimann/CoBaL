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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.s

function [A, ADot] = createConstraintMatrix_pointConstraints ...
  ( ...
    plant, ... 
    left_foot_constraint, ...
    right_foot_constraint ...
  )
    if any(isnan([left_foot_constraint, right_foot_constraint]))
        A = NaN;
        ADot = NaN;
        return
    end

    % left foot
    left_heel_body_jacobian = plant.bodyJacobians{plant.getEndEffectorIndex('left heel')};
    left_heel_body_jacobian_dot = plant.bodyJacobianTemporalDerivatives{plant.getEndEffectorIndex('left heel')};
    left_toes_body_jacobian = plant.bodyJacobians{plant.getEndEffectorIndex('LTOESEEF')};
    left_toes_body_jacobian_dot = plant.bodyJacobianTemporalDerivatives{plant.getEndEffectorIndex('LTOESEEF')};
    if left_foot_constraint == 0
        % left foot is in swing
        left_foot_constraint = [];
        left_foot_constraint_dot = [];
    elseif left_foot_constraint == 1
        % left foot is in heelstrike
        left_foot_constraint = left_heel_body_jacobian(1:3, :);
        left_foot_constraint_dot = left_heel_body_jacobian_dot(1:3, :);
    elseif left_foot_constraint == 2
        % left foot is in pushoff
        left_foot_constraint = left_toes_body_jacobian(1:3, :);
        left_foot_constraint_dot = left_toes_body_jacobian_dot(1:3, :);
    else
        error('Left foot constraint number must be an integer between 0 and 2')
    end
    
    % right foot
    right_heel_body_jacobian = plant.bodyJacobians{plant.getEndEffectorIndex('right heel')};
    right_heel_body_jacobian_dot = plant.bodyJacobianTemporalDerivatives{plant.getEndEffectorIndex('right heel')};
    right_toes_body_jacobian = plant.bodyJacobians{plant.getEndEffectorIndex('RTOESEEF')};
    right_toes_body_jacobian_dot = plant.bodyJacobianTemporalDerivatives{plant.getEndEffectorIndex('RTOESEEF')};
    
    if right_foot_constraint == 0
        % right foot is in swing
        right_foot_constraint = [];
        right_foot_constraint_dot = [];
    elseif right_foot_constraint == 1
        % right foot is in heelstrike
        right_foot_constraint = right_heel_body_jacobian(1:3, :);
        right_foot_constraint_dot = right_heel_body_jacobian_dot(1:3, :);
    elseif right_foot_constraint == 2
        % right foot is in pushoff
        right_foot_constraint = right_toes_body_jacobian(1:3, :);
        right_foot_constraint_dot = right_toes_body_jacobian_dot(1:3, :);
    else
        error('Right foot constraint number must be an integer between 0 and 2')
    end
    

%     right_heel_contact = constraint_binary(1);
%     right_toes_contact = constraint_binary(2);
%     left_heel_contact = constraint_binary(3);
%     left_toes_contact = constraint_binary(4);
% 
    
%     if ~right_heel_contact && ~right_toes_contact
%         % no constraint at all
%         right_foot_constraint = [];
%         right_foot_constraint_dot = [];
%     elseif right_heel_contact && ~right_toes_contact
%         % constrain heel sliding, leave heel rotation free
%         right_foot_constraint = right_heel_body_jacobian(1:3, :);
%         right_foot_constraint_dot = right_heel_body_jacobian_dot(1:3, :);
%     elseif ~right_heel_contact && right_toes_contact
%         % no toe sliding, leave toe rotation free
%         right_foot_constraint = right_toes_body_jacobian(1:3, :);
%         right_foot_constraint_dot = right_toes_body_jacobian_dot(1:3, :);
%     elseif right_heel_contact && right_toes_contact
%         % full contact, fully constrain foot segment
%         right_foot_constraint = right_heel_body_jacobian;
%         right_foot_constraint_dot = right_heel_body_jacobian_dot;
%     end
% 
%     if ~left_heel_contact && ~left_toes_contact
%         % no constraint at all
%         left_foot_constraint = [];
%         left_foot_constraint_dot = [];
%     elseif left_heel_contact && ~left_toes_contact
%         % constrain heel sliding, leave heel rotation free
%         left_foot_constraint = left_heel_body_jacobian(1:3, :);
%         left_foot_constraint_dot = left_heel_body_jacobian_dot(1:3, :);
%     elseif ~left_heel_contact && left_toes_contact
%         % no toe sliding, leave toe rotation free
%         left_foot_constraint = left_toes_body_jacobian(1:3, :);
%         left_foot_constraint_dot = left_toes_body_jacobian_dot(1:3, :);
%     elseif left_heel_contact && left_toes_contact
%         % full contact, fully constrain foot segment
%         left_foot_constraint = left_heel_body_jacobian;
%         left_foot_constraint_dot = left_heel_body_jacobian_dot;
%     end

    % concatenate to a single constraint matrix
    if numel(right_foot_constraint) + numel(left_foot_constraint) == 0
        A = zeros(1, plant.numberOfJoints);
        ADot = zeros(1, plant.numberOfJoints);
    else
        A = [right_foot_constraint; left_foot_constraint];
        ADot = [right_foot_constraint_dot; left_foot_constraint_dot];
    end

end