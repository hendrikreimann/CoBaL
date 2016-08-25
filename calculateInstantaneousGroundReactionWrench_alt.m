function groundReactionWrench = calculateInstantaneousGroundReactionWrench_alt(plant, constraintTorque, applicationFrame, applicationJoint)
    if isempty(constraintTorque)
        groundReactionWrench = [];
    else

        J_contact = plant.calculateArbitraryFrameBodyJacobian(applicationFrame, applicationJoint);
        groundReactionWrench = pinv(J_contact') * constraintTorque;
    end



end





















