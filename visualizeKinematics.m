%

% flags
plot_reconstruction_errors          = 0;
show_marker_comparison_stick_figure = 1;
plot_marker_paths                   = 0;
plot_foot_point_paths               = 0;
plot_foot_point_trajectories        = 0;
plot_angle_trajectories             = 0;

close_figures = 0;

trial_to_show = 3;

segments_to_plot = [6 9 10 12 15 16 18 21 24];
% segments_to_plot = [21 24];
segment_labels = ...
  { ... 
    'pelvis', ...
    'left thigh', 'left shank', 'left foot', ...
    'right thigh', 'right shank', 'right foot', ...
    'torso and arms', 'head and neck' ...
  };

marker_colors = {[1 0 0], [0 1 0], [0 0 1], [1 0.9 0], [1 0.5 0], [1 0 1]};
shade_factor = 0.6;

% load data
load subjectInfo.mat;
model_file_name = makeFileName(date, subject_id, 'model');
load(model_file_name);
marker_trajectories_file_name = makeFileName(date, subject_id, 'walking', trial_to_show, 'markerTrajectories');
load(marker_trajectories_file_name);
angle_trajectories_file_name = makeFileName(date, subject_id, 'walking', trial_to_show, 'angleTrajectories');
load(angle_trajectories_file_name);
kinematic_trajectories_file_name = makeFileName(date, subject_id, 'walking', trial_to_show, 'kinematicTrajectories');
load(kinematic_trajectories_file_name);


% figure save directory
% figure_save_dir = [data_root directorySeparator '..' directorySeparator 'figures' directorySeparator 'figures_trial_' num2str(trial_to_show)];
% [status, message, message_id] = mkdir(figure_save_dir);
%% plot_reconstruction_errors
if plot_reconstruction_errors
    position_map = ...
      [ ...
        600 900 500 400;
        600 450 500 400;
        600   0 500 400;
        1100 900 500 400;
        1100 450 500 400;
        1100   0 500 400;
        1600 900 500 400;
        1600 450 500 400;
        1600   0 500 400;
        2100 900 500 400;
        2100 450 500 400;
        2100   0 500 400;
      ];
    for i_joint = segments_to_plot
        figure('position', position_map(segments_to_plot==i_joint, :));
        number_of_markers = plant.getNumberOfMarkers(i_joint);

        for i_marker = 1 : number_of_markers
            marker_indices = plant.getMarkerExportIndices(i_joint, i_marker);
            marker_index = (marker_indices(1)-1) / 3 + 1;
            plot(time_mocap, marker_reconstruction_error(:, marker_index), '-', 'color', plant.getMarkerVisualizationColor(i_joint, i_marker), 'linewidth', 2)
%             semilogy(time, marker_reconstruction_error{trial_to_show}(:, marker_index), '-', 'color', marker_colors{i_marker}, 'linewidth', 2)
            hold on;
        end
%         title(['marker reconstruction errors, segment ' num2str(i_joint)]); xlabel('time (s)'); ylabel('error (m)'); 
        title(['marker reconstruction errors' segment_labels(segments_to_plot==i_joint)]); xlabel('time (s)'); ylabel('error (m)'); 
        

        save_file_name = ['reconstructionErrors_trial' num2str(trial_to_show) '_segment' num2str(i_joint)];
%         savefig([figure_save_dir directorySeparator save_file_name]);
        if close_figures; delete(gcf); end
    end
end

%% show_marker_comparison_stick_figure
if show_marker_comparison_stick_figure

%     scene_limits = 1*[-2.5 1.5; 0 1; -0.1 2];
%     scene_limits = 1*[-2.5 1.5; -1 2; -0.1 2];
    scene_limits = 1*[-0.8 0.8; 0 2; -0.1 2];
    stick_figure = showMarkerComparisonStickFigure(plant, angle_trajectories, marker_trajectories, scene_limits);
    
end
return
%% plot marker paths
if plot_marker_paths
    for i_joint = segments_to_plot
        figure; axes; hold on; axis equal;
        title(['marker trajectories, segment ' num2str(i_joint)]); 
        xlabel('x'); ylabel('y'); zlabel('z'); 
        number_of_markers = plant.getNumberOfMarkers(i_joint);
        for i_marker = 1 : number_of_markers
            marker_indices = plant.getMarkerExportIndices(i_joint, i_marker);
            plot3( ...
                   marker_trajectories(:, marker_indices(1)), ...
                   marker_trajectories(:, marker_indices(2)), ...
                   marker_trajectories(:, marker_indices(3)), ...
                   'color', marker_colors{i_marker}, ...
                   'linestyle', '-', ...
                   'linewidth', 1 ...
                 );
            plot3( ...
                   marker_trajectories_reconstructed(:, marker_indices(1)), ...
                   marker_trajectories_reconstructed(:, marker_indices(2)), ...
                   marker_trajectories_reconstructed(:, marker_indices(3)), ...
                   'color', marker_colors{i_marker}, ...
                   'linestyle', '--', ...
                   'linewidth', 2 ...
                 );
        end
        save_file_name = ['markerPaths_trial' num2str(trial_to_show) '_segment' num2str(i_joint)];
%         savefig([figure_save_dir directorySeparator save_file_name]);
        if close_figures; delete(gcf); end
    end
end

%% plot foot point paths
if plot_foot_point_paths
    % right foot
    figure; axes; hold on; axis equal;
    title('right foot ankle, heel and toes trajectories'); 
    xlabel('x'); ylabel('y'); zlabel('z');
    % shank markers
    shank_markers_joint_index = 10;
    number_of_markers = plant.getNumberOfMarkers(shank_markers_joint_index);
    for i_marker = 1 : number_of_markers
        marker_indices = plant.getMarkerExportIndices(shank_markers_joint_index, i_marker);
        plot3( ...
               marker_trajectories_processed(:, marker_indices(1)), ...
               marker_trajectories_processed(:, marker_indices(2)), ...
               marker_trajectories_processed(:, marker_indices(3)), ...
               'color', marker_colors{i_marker} * shade_factor, ...
               'linestyle', '-', ...
               'linewidth', 1 ...
             );
    end
    % foot markers
    left_foot_markers_joint_index = 12;
    number_of_markers = plant.getNumberOfMarkers(left_foot_markers_joint_index);
    for i_marker = 1 : number_of_markers
        marker_indices = plant.getMarkerExportIndices(left_foot_markers_joint_index, i_marker);
        plot3( ...
               marker_trajectories_processed(:, marker_indices(1)), ...
               marker_trajectories_processed(:, marker_indices(2)), ...
               marker_trajectories_processed(:, marker_indices(3)), ...
               'color', marker_colors{i_marker} * shade_factor, ...
               'linestyle', '-', ...
               'linewidth', 1 ...
             );
    end
    if show_optimized
        plot3( ...
               right_ankle_joint_trajectory_reconstructed_from_optimized(:, 1), ...
               right_ankle_joint_trajectory_reconstructed_from_optimized(:, 2), ...
               right_ankle_joint_trajectory_reconstructed_from_optimized(:, 3), ...
               'color', 'c', ...
               'linestyle', '-', ...
               'linewidth', 2 ...
             );
        plot3( ...
               right_heel_trajectory_reconstructed_from_optimized(:, 1), ...
               right_heel_trajectory_reconstructed_from_optimized(:, 2), ...
               right_heel_trajectory_reconstructed_from_optimized(:, 3), ...
               'color', 'm', ...
               'linestyle', '-', ...
               'linewidth', 2 ...
             );
        plot3( ...
               right_toes_trajectory_reconstructed_from_optimized(:, 1), ...
               right_toes_trajectory_reconstructed_from_optimized(:, 2), ...
               right_toes_trajectory_reconstructed_from_optimized(:, 3), ...
               'color', 'm', ...
               'linestyle', '-', ...
               'linewidth', 2 ...
             );
    end
    save_file_name = ['rightFootPointPaths_trial' num2str(trial_to_show)];
%     savefig([figure_save_dir directorySeparator save_file_name]);
    if close_figures
        delete(gcf);
    end
end

%% plot foot point trajectories
if plot_foot_point_trajectories
    
    % right foot
    figure; axes; hold on;
    title('right foot ankle, heel and toes trajectories'); 
    xlabel('time (s)'); ylabel('z-position (m)');
    plot(time, right_ankle_joint_trajectory_reconstructed_from_optimized(:, 3), 'r', 'linewidth', 2, 'displayname', 'right ankle joint');
    plot(time, right_heel_trajectory_reconstructed_from_optimized(:, 3), 'g', 'linewidth', 2, 'displayname', 'right heel');
    plot(time, right_toes_trajectory_reconstructed_from_optimized(:, 3), 'b', 'linewidth', 2, 'displayname', 'right toes');
    legend('toggle')
    right_foot_markers_joint_index = 12;
    number_of_markers = plant.getNumberOfMarkers(right_foot_markers_joint_index);
    for i_marker = 1 : number_of_markers
        marker_indices = plant.getMarkerExportIndices(right_foot_markers_joint_index, i_marker);
        plot(time, marker_trajectories_processed(:, marker_indices(3)), '-', 'color', marker_colors{i_marker} * shade_factor, 'linewidth', 1);
    end
    
    % left foot
    figure; axes; hold on;
    title('left foot ankle, heel and toes trajectories'); 
    xlabel('time (s)'); ylabel('z-position (m)');
    plot(time, left_ankle_joint_trajectory_reconstructed_from_optimized(:, 3), 'r-', 'linewidth', 2, 'displayname', 'left ankle joint');
    plot(time, left_heel_trajectory_reconstructed_from_optimized(:, 3), 'g-', 'linewidth', 2, 'displayname', 'left heel');
    plot(time, left_toes_trajectory_reconstructed_from_optimized(:, 3), 'b-', 'linewidth', 2, 'displayname', 'left toes');
    legend('toggle')
    left_foot_markers_joint_index = 18;
    number_of_markers = plant.getNumberOfMarkers(left_foot_markers_joint_index);
    for i_marker = 1 : number_of_markers
        marker_indices = plant.getMarkerExportIndices(left_foot_markers_joint_index, i_marker);
        plot(time, marker_trajectories_processed(:, marker_indices(3)), '-', 'color', marker_colors{i_marker} * shade_factor, 'linewidth', 1);
    end
    
end

%% plot angle trajectories
if plot_angle_trajectories
    
    figure('position', [0 1050 600 300]); axes; hold on; title(['angle trajectories, virtual joints, trial ' num2str(trial_to_show)]); 
    nan_indicator = ~any(isnan(joint_angle_trajectory_optimized), 2);
    notnan_indicator = any(isnan(joint_angle_trajectory_optimized), 2);
    z_nan = zeros(size(time)); z_nan(nan_indicator) = NaN;
    z_notnan = zeros(size(time)); z_notnan(notnan_indicator) = NaN;
    plot(joint_angle_trajectory_optimized(:, 1), 'linewidth', 2, 'displayname', 'optimized - trunk configuration, x-translation')
    plot(joint_angle_trajectory_optimized(:, 2), 'linewidth', 2, 'displayname', 'optimized - trunk configuration, y-translation')
    plot(joint_angle_trajectory_optimized(:, 3), 'linewidth', 2, 'displayname', 'optimized - trunk configuration, z-translation')
    plot(joint_angle_trajectory_optimized(:, 4), 'linewidth', 2, 'displayname', 'optimized - trunk configuration, z-rotation') 
    plot(joint_angle_trajectory_optimized(:, 5), 'linewidth', 2, 'displayname', 'optimized - trunk configuration, x-rotation')
    plot(joint_angle_trajectory_optimized(:, 6), 'linewidth', 2, 'displayname', 'optimized - trunk configuration, y-rotation')
    plot(z_nan, 'color', [1 0 0], 'linewidth', 1, 'marker', '+', 'displayname', 'non-existing data indicator')
    plot(z_notnan, 'color', [0 1 0], 'linewidth', 1, 'displayname', 'existing data indicator')
    legend('show', 'location', 'SE')
    save_file_name = ['jointAngles_trial' num2str(trial_to_show) '_virtual'];
%     savefig([figure_save_dir directorySeparator save_file_name]);
    if close_figures
        delete(gcf);
    end
    
    figure('position', [0 700 600 300]); axes; hold on;
    title(['angle trajectories, right leg joints, trial ' num2str(trial_to_show)])
    plot(time, joint_angle_trajectory_optimized(:, 7), 'linewidth', 2, 'displayname', 'optimized - hip internal rotation')
    plot(time, joint_angle_trajectory_optimized(:, 8), 'linewidth', 2, 'displayname', 'optimized - hip ab/adduction')
    plot(time, joint_angle_trajectory_optimized(:, 9), 'linewidth', 2, 'displayname', 'optimized - hip flexion/extension')
    plot(time, joint_angle_trajectory_optimized(:, 10), 'linewidth', 2, 'displayname', 'optimized - knee flexion/extension')
    plot(time, joint_angle_trajectory_optimized(:, 11), 'linewidth', 2, 'displayname', 'optimized - ankle plantar/dorsiflexion')
    plot(time, joint_angle_trajectory_optimized(:, 12), 'linewidth', 2, 'displayname', 'optimized - ankle inversion/eversion')
    legend('toggle')
    save_file_name = ['jointAngles_trial' num2str(trial_to_show) '_rightLeg'];
%     savefig([figure_save_dir directorySeparator save_file_name]);
    if close_figures
        delete(gcf);
    end

    figure('position', [0 350 600 300]); axes; hold on;
    title(['angle trajectories, left leg joints, trial ' num2str(trial_to_show)])
    plot(time, joint_angle_trajectory_optimized(:, 13), 'linewidth', 2, 'displayname', 'optimized - hip internal rotation')
    plot(time, joint_angle_trajectory_optimized(:, 14), 'linewidth', 2, 'displayname', 'optimized - hip ab/adduction')
    plot(time, joint_angle_trajectory_optimized(:, 15), 'linewidth', 2, 'displayname', 'optimized - hip flexion/extension')
    plot(time, joint_angle_trajectory_optimized(:, 16), 'linewidth', 2, 'displayname', 'optimized - knee flexion/extension')
    plot(time, joint_angle_trajectory_optimized(:, 17), 'linewidth', 2, 'displayname', 'optimized - ankle plantar/dorsiflexion')
    plot(time, joint_angle_trajectory_optimized(:, 18), 'linewidth', 2, 'displayname', 'optimized - ankle inversion/eversion')
    legend('toggle')
    save_file_name = ['jointAngles_trial' num2str(trial_to_show) '_leftLeg'];
%     savefig([figure_save_dir directorySeparator save_file_name]);
    if close_figures
        delete(gcf);
    end

%     figure('position', [0 0 600 300]); axes; hold on;
%     title(['angle trajectories, lumbar and cervical joints, trial ' num2str(trial_to_show)])
%     plot(time, joint_angle_trajectory_optimized(:, 19), 'linewidth', 2, 'displayname', 'optimized - l5 - z')
%     plot(time, joint_angle_trajectory_optimized(:, 20), 'linewidth', 2, 'displayname', 'optimized - l5 - x')
%     plot(time, joint_angle_trajectory_optimized(:, 21), 'linewidth', 2, 'displayname', 'optimized - l5 - y')
%     plot(time, joint_angle_trajectory_optimized(:, 22), 'linewidth', 2, 'displayname', 'optimized - neck - z')
%     plot(time, joint_angle_trajectory_optimized(:, 23), 'linewidth', 2, 'displayname', 'optimized - neck - x')
%     plot(time, joint_angle_trajectory_optimized(:, 24), 'linewidth', 2, 'displayname', 'optimized - neck - y')
%     legend('toggle')
%     save_file_name = ['jointAngles_trial' num2str(trial_to_show) '_leftLeg'];
%     savefig([figure_save_dir directorySeparator save_file_name]);
%     if close_figures
%         delete(gcf);
%     end
end



