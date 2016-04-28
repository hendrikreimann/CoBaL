function [A, ADot] = createConstraintMatrix_bodyVelocityConstraints ...
  ( ...
    plant, ... 
    left_foot_constraint, ...
    right_foot_constraint, ...
    phi_left, ...
    rho_left, ...
    phi_right, ...
    rho_right, ...
    V_body_left_fits_heelstrike, ...
    V_body_left_fits_pushoff, ...
    V_body_right_fits_heelstrike, ...
    V_body_right_fits_pushoff ...
  )

    % left foot
    if left_foot_constraint == 0
        % left foot is in swing
        left_foot_constraint = [];
        left_foot_constraint_dot = [];
    else 
        left_ankle_body_jacobian = plant.bodyJacobians{3};
        left_ankle_body_jacobian_dot = plant.bodyJacobianTemporalDerivatives{3};
        
        % calculate body velocity for lateral roll
        left_ankle_rho_rotation_joint_index = 12;
        rho_rotation_center_local = [0; 0; 0; 1]; % assume the ankle is the center of a sphere resting on the ground
        left_sphere_transformation = reshape(plant.jointTransformations{left_ankle_rho_rotation_joint_index}, 4, 4);
        rho_rotation_center_world = left_sphere_transformation * rho_rotation_center_local;
        p_virtual_contact = [rho_rotation_center_world(1:2); 0; 1];
        J_body_virtual_contact = plant.calculateArbitraryFrameBodyJacobian([eye(4, 3) p_virtual_contact], left_ankle_rho_rotation_joint_index);
        A_rho = J_body_virtual_contact([1:3 6], :);
        P_A_rho = A_rho' * (A_rho * A_rho')^(-1) * A_rho; % projection onto the space spanned by the columns of A
        P_A_rho_orth = eye(plant.numberOfJoints) - P_A_rho; % projection onto the null space of A
        joint_velocity_rho_change_unconstrained = zeros(plant.numberOfJoints, 1);
        joint_velocity_rho_change_unconstrained(left_ankle_rho_rotation_joint_index) = 1;
        joint_velocity_rho_change_constrained = (P_A_rho_orth * joint_velocity_rho_change_unconstrained);
        body_velocity_rho_change_constrained = left_ankle_body_jacobian * joint_velocity_rho_change_constrained;
        V_rho_body = body_velocity_rho_change_constrained * 1 / body_velocity_rho_change_constrained(4);
        
        % calculate body velocity direction for anterior roll
        V_phi_body = zeros(6, 1);
        if left_foot_constraint == 1
            V_phi_body(1) = feval(V_body_left_fits_heelstrike{1},[phi_left, rho_left]);
            V_phi_body(2) = feval(V_body_left_fits_heelstrike{2},[phi_left, rho_left]);
            V_phi_body(3) = feval(V_body_left_fits_heelstrike{3},[phi_left, rho_left]);
            V_phi_body(4) = feval(V_body_left_fits_heelstrike{4},[phi_left, rho_left]);
            V_phi_body(5) = feval(V_body_left_fits_heelstrike{5},[phi_left, rho_left]);
            V_phi_body(6) = feval(V_body_left_fits_heelstrike{6},[phi_left, rho_left]);
        elseif left_foot_constraint==2
            V_phi_body(1) = feval(V_body_left_fits_pushoff{1},[phi_left, rho_left]);
            V_phi_body(2) = feval(V_body_left_fits_pushoff{2},[phi_left, rho_left]);
            V_phi_body(3) = feval(V_body_left_fits_pushoff{3},[phi_left, rho_left]);
            V_phi_body(4) = feval(V_body_left_fits_pushoff{4},[phi_left, rho_left]);
            V_phi_body(5) = feval(V_body_left_fits_pushoff{5},[phi_left, rho_left]);
            V_phi_body(6) = feval(V_body_left_fits_pushoff{6},[phi_left, rho_left]);
        else            
            error('Left foot constraint number must be an integer between 0 and 2')
        end
            
            
        % transform to joint constraints
        V_body_allowed = [V_phi_body V_rho_body];
        [~, ~, V] = svd(V_body_allowed');
        C = V(:, 3:6)'; % orthogonal complement of the allowed body velocity directions
        left_foot_constraint = C * left_ankle_body_jacobian;
        left_foot_constraint_dot = C * left_ankle_body_jacobian_dot; % TODO: add C_dot * left_ankle_body_jacobian
        
%         W_rho = pinv(left_ankle_body_jacobian) * V_rho_body;
%         W_phi = pinv(left_ankle_body_jacobian) * V_phi_body;
%         [~, ~, V_body] = svd(left_ankle_body_jacobian);
%         B = V_body(:, 7 : end); % this is the null space of the body Jacobian, i.e. the space of joint changes that do not move the foot
%         C = [B W_phi W_rho]; % this is the space of unconstrained joint changes
%         [~, ~, V_constraint] = svd(C');
%         left_foot_constraint = V_constraint(:, plant.numberOfJoints-3 : plant.numberOfJoints)';        
%         left_foot_constraint_dot = zeros(size(left_foot_constraint)); % ACHTUNG: I just put zeros here for now, this has to be calculated appropriately later!
    end
    
    % right foot
    if right_foot_constraint == 0
        % right foot is in swing
        right_foot_constraint = [];
        right_foot_constraint_dot = [];
    else 
        right_ankle_body_jacobian = plant.bodyJacobians{6};
        right_ankle_body_jacobian_dot = plant.bodyJacobianTemporalDerivatives{6};
        
        % calculate body velocity for lateral roll
        right_ankle_rho_rotation_joint_index = 18;
        rho_rotation_center_local = [0; 0; 0; 1]; % assume the ankle is the center of a sphere resting on the ground
        right_sphere_transformation = reshape(plant.jointTransformations{right_ankle_rho_rotation_joint_index}, 4, 4);
        rho_rotation_center_world = right_sphere_transformation * rho_rotation_center_local;
        p_virtual_contact = [rho_rotation_center_world(1:2); 0; 1];
        J_body_virtual_contact = plant.calculateArbitraryFrameBodyJacobian([eye(4, 3) p_virtual_contact], right_ankle_rho_rotation_joint_index);
        A_rho = J_body_virtual_contact([1:3 6], :);
        P_A_rho = A_rho' * (A_rho * A_rho')^(-1) * A_rho; % projection onto the space spanned by the columns of A
        P_A_rho_orth = eye(plant.numberOfJoints) - P_A_rho; % projection onto the null space of A
        joint_velocity_rho_change_unconstrained = zeros(plant.numberOfJoints, 1);
        joint_velocity_rho_change_unconstrained(right_ankle_rho_rotation_joint_index) = 1;
        joint_velocity_rho_change_constrained = (P_A_rho_orth * joint_velocity_rho_change_unconstrained);
        body_velocity_rho_change_constrained = right_ankle_body_jacobian * joint_velocity_rho_change_constrained;
        V_rho_body = body_velocity_rho_change_constrained * 1 / body_velocity_rho_change_constrained(4);
        
        % calculate body velocity direction for anterior roll
        V_phi_body = zeros(6, 1);
        if right_foot_constraint == 1
            V_phi_body(1) = feval(V_body_right_fits_heelstrike{1},[phi_right, rho_right]);
            V_phi_body(2) = feval(V_body_right_fits_heelstrike{2},[phi_right, rho_right]);
            V_phi_body(3) = feval(V_body_right_fits_heelstrike{3},[phi_right, rho_right]);
            V_phi_body(4) = feval(V_body_right_fits_heelstrike{4},[phi_right, rho_right]);
            V_phi_body(5) = feval(V_body_right_fits_heelstrike{5},[phi_right, rho_right]);
            V_phi_body(6) = feval(V_body_right_fits_heelstrike{6},[phi_right, rho_right]);
        elseif right_foot_constraint==2
            V_phi_body(1) = feval(V_body_right_fits_pushoff{1},[phi_right, rho_right]);
            V_phi_body(2) = feval(V_body_right_fits_pushoff{2},[phi_right, rho_right]);
            V_phi_body(3) = feval(V_body_right_fits_pushoff{3},[phi_right, rho_right]);
            V_phi_body(4) = feval(V_body_right_fits_pushoff{4},[phi_right, rho_right]);
            V_phi_body(5) = feval(V_body_right_fits_pushoff{5},[phi_right, rho_right]);
            V_phi_body(6) = feval(V_body_right_fits_pushoff{6},[phi_right, rho_right]);
        else            
            error('Left foot constraint number must be an integer between 0 and 2')
        end
            
            
        % transform to joint constraints
        V_body_allowed = [V_phi_body V_rho_body];
        [~, ~, V] = svd(V_body_allowed');
        C = V(:, 3:6)'; % orthogonal complement of the allowed body velocity directions
        right_foot_constraint = C * right_ankle_body_jacobian;
        right_foot_constraint_dot = C * right_ankle_body_jacobian_dot; % TODO: add C_dot * right_ankle_body_jacobian
        
%         W_rho = pinv(right_ankle_body_jacobian) * V_rho_body;
%         W_phi = pinv(right_ankle_body_jacobian) * V_phi_body;
%         [~, ~, V_body] = svd(right_ankle_body_jacobian);
%         B = V_body(:, 7 : end); % this is the null space of the body Jacobian, i.e. the space of joint changes that do not move the foot
%         C = [B W_phi W_rho]; % this is the space of unconstrained joint changes
%         [~, ~, V_constraint] = svd(C');
%         right_foot_constraint = V_constraint(:, plant.numberOfJoints-3 : plant.numberOfJoints)';        
%         right_foot_constraint_dot = zeros(size(right_foot_constraint)); % ACHTUNG: I just put zeros here for now, this has to be calculated appropriately later!
    end

    % concatenate to a single constraint matrix
    if numel(right_foot_constraint) + numel(left_foot_constraint) == 0
        A = zeros(1, plant.numberOfJoints);
        ADot = zeros(1, plant.numberOfJoints);
    else
        A = [right_foot_constraint; left_foot_constraint];
        ADot = [right_foot_constraint_dot; left_foot_constraint_dot];
    end
end