% figure; plot(left_ankle_acceleration_strategy(:, 1:6));
figure; plot(left_ankle_acceleration_strategy(:, 7:12));
% figure; plot(left_hip_acceleration_strategy(:, 1:6));
figure; plot(left_hip_acceleration_strategy(:, 7:12));
% figure; plot(right_ankle_acceleration_strategy(:, 1:6));
figure; plot(right_ankle_acceleration_strategy(:, 13:18));
% figure; plot(right_hip_acceleration_strategy(:, 1:6));
figure; plot(right_hip_acceleration_strategy(:, 13:18));
return
% A_left_trajectory_unfolded = zeros(256, 5, 38);
% J_c_left_trajectory_unfolded = zeros(256, 2, 38);
% for i_time = 1 : 256
%     A_left_trajectory_unfolded(i_time, :, :) = A_left_trajectory{i_time};
%     J_c_left_trajectory_unfolded(i_time, :, :) = J_c_left_trajectory{i_time};
% end

% for i_row = 1 : 5
%     figure; hold on
%     for i_joint = 1 : 38;
%         plot3(1:256, A_left_trajectory_unfolded(:, i_row, i_joint), i_joint*ones(256, 1));
%     end
% end

% for i_row = 1 : 2
%     figure; hold on
%     for i_joint = 1 : 38;
%         plot3(1:256, J_c_left_trajectory_unfolded(:, i_row, i_joint), i_joint*ones(256, 1));
%     end
% end

scene_limits = [-0.8 0.8; 56 58; -0.1 2];
stick_figure = showMarkerComparisonStickFigure(plant, joint_angles_mean_left_control, [], scene_limits);
stick_figure.showLinkMassEllipsoids = false;
stick_figure.update;

scene_limits = [-0.8 0.8; 50 52; -0.1 2];
stick_figure = showMarkerComparisonStickFigure(plant, joint_angles_mean_right_control, [], scene_limits);
stick_figure.showLinkMassEllipsoids = false;
stick_figure.update;




