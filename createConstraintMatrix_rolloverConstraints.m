function [A, ADot] = createConstraintMatrix_rolloverConstraints ...
  ( ...
    plant, ... 
    rollover_shapes, ...
    left_foot_constraint_index, ...
    right_foot_constraint_index, ...
    theta_left, ...
    theta_right, ...
    T_foot_to_world_left, ...
    T_foot_to_world_right ...
  )
    left_ankle_joint_index = 12;
    right_ankle_joint_index = 18;

    % get contact point in foot coordinates
    contact_point_foot = [0; rollover_shapes.left_fit_theta_x(theta_left); rollover_shapes.left_fit_theta_y(theta_left); 1];
    contact_point_world = T_foot_to_world_left * contact_point_foot;
    
    foot_frame_body_Jacobian_left = plant.calculateArbitraryFrameBodyJacobian(T_foot_to_world_left, left_ankle_joint_index);
    ankle_point_Jacobian_left = plant.calculateArbitraryPointJacobian(T_foot_to_world_left(:, 4), left_ankle_joint_index, 'world');
    ankle_point_lateral_Jacobian_left = foot_frame_body_Jacobian_left(1, :);
    foot_vector_rotation_Jacobian_left = foot_frame_body_Jacobian_left(5, :);
    contact_point_Jacobian_left = plant.calculateArbitraryPointJacobian(contact_point_world, left_ankle_joint_index, 'world');

    left_foot_constraint_matrix = ...
      [ ...
        ankle_point_lateral_Jacobian_left; ...
        foot_vector_rotation_Jacobian_left; ...
        contact_point_Jacobian_left ...
      ];
  
% test by projecting an arbitrary joint velocity onto this constraint
A = left_foot_constraint_matrix;
P_A = A' * (A * A')^(-1) * A; % projection onto the space spanned by A
P_A_orth = eye(plant.numberOfJoints) - P_A;

joint_velocity_unconstrained = ones(plant.numberOfJoints, 1);
joint_velocity_constrained = (P_A_orth * joint_velocity_unconstrained);
    

    
    
    
    
    
    
    
    
    
    % left foot
    left_heel_body_jacobian = plant.bodyJacobians{1};
    left_heel_body_jacobian_dot = plant.bodyJacobianTemporalDerivatives{1};
    left_toes_body_jacobian = plant.bodyJacobians{2};
    left_toes_body_jacobian_dot = plant.bodyJacobianTemporalDerivatives{2};
    if left_foot_constraint_index == 0
        % left foot is in swing
        left_foot_constraint_matrix = [];
        left_foot_constraint_dot = [];
    elseif left_foot_constraint_index == 1
        % left foot is in heelstrike
        left_foot_constraint_matrix = left_heel_body_jacobian(1:5, :);
        left_foot_constraint_dot = left_heel_body_jacobian_dot(1:5, :);
    elseif left_foot_constraint_index == 2
        % left foot is in pushoff
        left_foot_constraint_matrix = left_toes_body_jacobian(1:5, :);
        left_foot_constraint_dot = left_toes_body_jacobian_dot(1:5, :);
    else
        error('Left foot constraint number must be an integer between 0 and 2')
    end
    
    % right foot
    right_heel_body_jacobian = plant.bodyJacobians{4};
    right_heel_body_jacobian_dot = plant.bodyJacobianTemporalDerivatives{4};
    right_toes_body_jacobian = plant.bodyJacobians{5};
    right_toes_body_jacobian_dot = plant.bodyJacobianTemporalDerivatives{5};
    
    if right_foot_constraint_index == 0
        % right foot is in swing
        right_foot_constraint_matrix = [];
        right_foot_constraint_dot = [];
    elseif right_foot_constraint_index == 1
        % right foot is in heelstrike
        right_foot_constraint_matrix = right_heel_body_jacobian(1:5, :);
        right_foot_constraint_dot = right_heel_body_jacobian_dot(1:5, :);
    elseif right_foot_constraint_index == 2
        % right foot is in pushoff
        right_foot_constraint_matrix = right_toes_body_jacobian(1:5, :);
        right_foot_constraint_dot = right_toes_body_jacobian_dot(1:5, :);
    else
        error('Right foot constraint number must be an integer between 0 and 2')
    end

    % concatenate to a single constraint matrix
    if numel(right_foot_constraint_matrix) + numel(left_foot_constraint_matrix) == 0
        A = zeros(1, plant.numberOfJoints);
        ADot = zeros(1, plant.numberOfJoints);
    else
        A = [right_foot_constraint_matrix; left_foot_constraint_matrix];
        ADot = [right_foot_constraint_dot; left_foot_constraint_dot];
    end








%     constraint_binary = de2bi(constraintDecimal, 4);
%     right_heel_contact = constraint_binary(1);
%     right_toes_contact = constraint_binary(2);
%     left_heel_contact = constraint_binary(3);
%     left_toes_contact = constraint_binary(4);
% 
%     right_heel_body_jacobian = plant.bodyJacobians{1};
%     right_heel_body_jacobian_dot = plant.bodyJacobianTemporalDerivatives{1};
%     right_toes_body_jacobian = plant.bodyJacobians{2};
%     right_toes_body_jacobian_dot = plant.bodyJacobianTemporalDerivatives{2};
%     left_heel_body_jacobian = plant.bodyJacobians{4};
%     left_heel_body_jacobian_dot = plant.bodyJacobianTemporalDerivatives{4};
%     left_toes_body_jacobian = plant.bodyJacobians{5};
%     left_toes_body_jacobian_dot = plant.bodyJacobianTemporalDerivatives{5};
%     
%     if ~right_heel_contact && ~right_toes_contact
%         % no constraint at all
%         right_foot_constraint = [];
%         right_foot_constraint_dot = [];
%     elseif right_heel_contact && ~right_toes_contact
%         % constrain heel sliding, leave heel rotation free around hinge
%         right_foot_constraint = right_heel_body_jacobian(1:5, :);
%         right_foot_constraint_dot = right_heel_body_jacobian_dot(1:5, :);
%     elseif ~right_heel_contact && right_toes_contact
%         % no toe sliding, leave toe rotation free around hinge
%         right_foot_constraint = right_toes_body_jacobian(1:5, :);
%         right_foot_constraint_dot = right_toes_body_jacobian_dot(1:5, :);
%     elseif right_heel_contact && right_toes_contact
%         % full contact, fully constrain foot segment but leave toe segment free
%         right_foot_constraint = right_heel_body_jacobian;
%         right_foot_constraint_dot = right_heel_body_jacobian_dot;
%     end
% 
%     if ~left_heel_contact && ~left_toes_contact
%         % no constraint at all
%         left_foot_constraint = [];
%         left_foot_constraint_dot = [];
%     elseif left_heel_contact && ~left_toes_contact
%         % constrain heel sliding, leave heel rotation free around hinge
%         left_foot_constraint = left_heel_body_jacobian(1:5, :);
%         left_foot_constraint_dot = left_heel_body_jacobian_dot(1:5, :);
%     elseif ~left_heel_contact && left_toes_contact
%         % no toe sliding, leave toe rotation free around hinge
%         left_foot_constraint = left_toes_body_jacobian(1:5, :);
%         left_foot_constraint_dot = left_toes_body_jacobian_dot(1:5, :);
%     elseif left_heel_contact && left_toes_contact
%         % full contact, fully constrain foot segment but leave toe segment free
%         left_foot_constraint = left_heel_body_jacobian;
%         left_foot_constraint_dot = left_heel_body_jacobian_dot;
%     end
% 
%     % concatenate to a single constraint matrix
%     if numel(right_foot_constraint) + numel(left_foot_constraint) == 0
%         A = zeros(1, plant.numberOfJoints);
%         ADot = zeros(1, plant.numberOfJoints);
%     else
%         A = [right_foot_constraint; left_foot_constraint];
%         ADot = [right_foot_constraint_dot; left_foot_constraint_dot];
%     end

end