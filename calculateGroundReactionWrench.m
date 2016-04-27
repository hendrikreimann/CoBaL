function groundReactionWrench = calculateGroundReactionWrench(plant, constraintTorque, applicationFrame)

    % transform into ground reaction wrench automatically - go backwards, use only previous joints
    ground_reaction_wrenches = zeros(6, plant.numberOfJoints);
    resistable_torques = zeros(plant.numberOfJoints, plant.numberOfJoints);
    resisted_torques = zeros(plant.numberOfJoints, plant.numberOfJoints);
    remaining_torques = zeros(plant.numberOfJoints, plant.numberOfJoints);
    remaining_torque = constraintTorque;
    for i_joint = plant.numberOfJoints : - 1 : 1
        % calculate which part of the remaining torque can be resisted at this segment
        application_frame_fixed_to_joint_i_Jacobian = plant.calculateArbitraryFrameBodyJacobian(applicationFrame, i_joint);
        [~, ~, V_J] = svd(application_frame_fixed_to_joint_i_Jacobian);
        k_J = rank(application_frame_fixed_to_joint_i_Jacobian);
        B_J = V_J(:, 1:k_J); % basis of the range space of J
        C_J = V_J(:, k_J+1:end); % basis of the null space of J
        C_J(i_joint, :) = 0; % basis of the null space of the desired projection
        [~, ~, V_N] = svd(C_J');

        k_N = rank(C_J');
        B_N = V_N(:, 1:k_N); % basis of the null space of the desired projection
        C_N = V_N(:, k_N+1:end); % basis of the orthogonal complement of the null space of the desired projection


    %     i_joint = i_joint
    %     rank(B_J)
    %     rank(C_J)

        P = B_J * (C_N'*B_J)^(-1) * C_N'; % projection onto the range space of J, and the i-th component lying completely NOT in the null space of the projection (i.e. it will be left invariant)
        resistable_torque = P * remaining_torque;
        remaining_torque_next = remaining_torque - resistable_torque;


        % calculate the wrench at the i-th segment that resists this part of the torque
        resisting_wrench = pinv(application_frame_fixed_to_joint_i_Jacobian') * resistable_torque;
        ground_reaction_wrenches(:, i_joint) = resisting_wrench;

    %     check = (origin_frame_fixed_to_joint_i_Jacobian' * resisting_wrench)'
    %     resistable_torque'
    %     
        remaining_torque = remaining_torque_next;
        resistable_torques(:, i_joint) = resistable_torque;
        resisted_torques(:, i_joint) = (application_frame_fixed_to_joint_i_Jacobian' * resisting_wrench)';
        remaining_torques(:, i_joint) = remaining_torque;
    end
    groundReactionWrench = sum(ground_reaction_wrenches, 2);


end





















