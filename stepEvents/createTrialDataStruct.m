function trialData = createTrialDataStruct
    trialData = struct;

    trialData.time_mocap = [];
    trialData.sampling_rate_mocap = NaN;

    trialData.left_heel_z_pos_trajectory = [];
    trialData.left_heel_z_vel_trajectory = [];
    trialData.left_heel_z_acc_trajectory = [];

    trialData.right_heel_z_pos_trajectory = [];
    trialData.right_heel_z_vel_trajectory = [];
    trialData.right_heel_z_acc_trajectory = [];

    trialData.left_toes_z_pos_trajectory = [];
    trialData.left_toes_z_vel_trajectory = [];
    trialData.left_toes_z_acc_trajectory = [];

    trialData.right_toes_z_pos_trajectory = [];
    trialData.right_toes_z_vel_trajectory = [];
    trialData.right_toes_z_acc_trajectory = [];
    

end