%

% flags
plot_reconstruction_errors          = 0;
show_marker_comparison_stick_figure = 1;
plot_marker_paths                   = 0;
plot_angle_trajectories             = 0;

close_figures = 0;

trial_to_show = 2;

segments_to_plot = ...
  [ ...
    6 ...
    9 10 12 ...
    15 16 18 ...
    21 24 ...
    27 29 31 ...
    34 36 38 ...
  ];
% segments_to_plot = [21 24];
segment_labels = ...
  { ... 
    'pelvis', ...
    'left thigh', 'left shank', 'left foot', ...
    'right thigh', 'right shank', 'right foot', ...
    'torso and arms', 'head and neck', ...
    'left upper arm', 'left lower arm', 'left hand', ...
    'right upper arm', 'right lower arm', 'right hand', ...
  };

marker_colors = {[1 0 0], [0 1 0], [0 0 1], [1 0.9 0], [1 0.5 0], [1 0 1], [1 0 0], [0 1 0], [0 0 1], [1 0.9 0], [1 0.5 0], [1 0 1], [1 0 0], [0 1 0], [0 0 1], [1 0.9 0], [1 0.5 0], [1 0 1], [1 0 0], [0 1 0], [0 0 1], [1 0.9 0], [1 0.5 0], [1 0 1]};
shade_factor = 0.6;

% load data
load subjectInfo.mat;
load(makeFileName(date, subject_id, 'model'));
load(makeFileName(date, subject_id, 'walking', trial_to_show, 'markerTrajectories'));
load(makeFileName(date, subject_id, 'walking', trial_to_show, 'angleTrajectories'));
% load(makeFileName(date, subject_id, 'walking', trial_to_show, 'kinematicTrajectories'));


angle_plot_groups = ...
  { ...
     1 : 6; ...
     7 : 12; ...
     13 : 18; ...
     19 : 24 ...
   };

%% plot_reconstruction_errors
if plot_reconstruction_errors
    reconstruction_plot_axes = [];
    for i_joint = segments_to_plot
        figure; reconstruction_plot_axes_segment = axes; hold on;
        number_of_markers = plant.getNumberOfMarkers(i_joint);

        for i_marker = 1 : number_of_markers
            marker_indices = plant.getMarkerExportIndices(i_joint, i_marker);
            marker_index = (marker_indices(1)-1) / 3 + 1;
            plot(time_mocap, marker_reconstruction_error(:, marker_index), '-', 'color', plant.getMarkerVisualizationColor(i_joint, i_marker), 'linewidth', 2)
%             semilogy(time, marker_reconstruction_error{trial_to_show}(:, marker_index), '-', 'color', marker_colors{i_marker}, 'linewidth', 2)
        end
        title(['marker reconstruction errors' segment_labels(segments_to_plot==i_joint)]); xlabel('time (s)'); ylabel('error (m)'); 
        reconstruction_plot_axes = [reconstruction_plot_axes reconstruction_plot_axes_segment];
    end
end

%% show_marker_comparison_stick_figure
if show_marker_comparison_stick_figure

    scene_limits = [-0.8 0.8; 0 20; -0.1 2];
    scene_limits = [-0.8 0.8; 0 2; -0.1 2];
    stick_figure = showMarkerComparisonStickFigure(plant, joint_angle_trajectories, marker_trajectories, scene_limits);
    stick_figure.showLinkMassEllipsoids = false;
    stick_figure.update;
end

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
    end
end


%% plot angle trajectories
if plot_angle_trajectories
    number_of_groups = size(angle_plot_groups, 1);
    group_axes = zeros(1, number_of_groups);
    for i_group = 1 : number_of_groups
        figure; group_axes(i_group) = axes; hold on
        for i_joint = angle_plot_groups{i_group}
            plot(time_mocap, joint_angle_trajectories(:, i_joint), 'linewidth', 2, 'displayname', plant.jointLabels{i_joint})
        end
        legend('show', 'location', 'SE')
    end
    linkaxes(group_axes, 'x');
    distFig('rows', number_of_groups);
    
    

end



