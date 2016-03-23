% estimates the center of rotation (CoR) of a joint using the SCoRE algorithm from Ehrig et al., "A
% survey of formal methods for determining the centre of rotation of ball joints. J. Biomech 39:2798-2809." 
% 
% Hendrik Reimann, 2014
% Department of Kinesiology
% Temple University
% hendrikreimann@gmail.com

% in:
% markerClusterOne - N x 3M matrix with 3d coordinates of M markers at N points in time
% markerClusterTwo - N x 3M matrix with 3d coordinates of M markers at N points in time

% out:
% CoR = center of rotation, 3x1 matrix
% meanError = mean difference between mean radius and current radius (should be constant) for all markers
%

function [CoR, meanError, AoR] = estimateJointKinematics(markerClusterOneReference, markerClusterTwoReference, markerClusterOneTrajectory, markerClusterTwoTrajectory, degreesOfFreedom)

    % remove the NaN entries in the trajectories
    nan_entries = any(isnan([markerClusterOneTrajectory markerClusterTwoTrajectory]), 2);
    markerClusterOneTrajectory(nan_entries, :) = [];
    markerClusterTwoTrajectory(nan_entries, :) = [];
    
    
    if nargin < 5
        degreesOfFreedom = 3;
    end

    AoR = zeros(3, degreesOfFreedom);

    number_of_time_frames = size(markerClusterOneTrajectory, 1);
    number_of_markers_one = length(markerClusterOneReference) / 3;
    number_of_markers_two = length(markerClusterTwoReference) / 3;
    
    cluster_one_reference_world = reshape(markerClusterOneReference, 3, number_of_markers_one);
    cluster_two_reference_world = reshape(markerClusterTwoReference, 3, number_of_markers_two);

    cluster_one_reference_center_world = mean(cluster_one_reference_world, 2);
    cluster_two_reference_center_world = mean(cluster_two_reference_world, 2);

    % calculate rotations and translations and save them to a combined matrix
    A = zeros(3*number_of_time_frames, 6);
    b = zeros(3*number_of_time_frames, 1);
    for i_time = 1:number_of_time_frames
        % calculate cluster one rotation and translation
        cluster_one_current_world = reshape(markerClusterOneTrajectory(i_time, :), 3, number_of_markers_one);
        p_cluster_one_to_world = mean(cluster_one_current_world, 2);
        R_cluster_one_to_world = rotationMatrixFromMarkers(cluster_one_reference_world, cluster_one_current_world);

        % calculate cluster two rotation and translation
        cluster_two_current_world = reshape(markerClusterTwoTrajectory(i_time, :), 3, number_of_markers_two);
        p_cluster_two_to_world = mean(cluster_two_current_world, 2);
        R_cluster_two_to_world = rotationMatrixFromMarkers(cluster_two_reference_world, cluster_two_current_world);

        % put it into the appropriate place in the A and b matrices
        first = 3*(i_time-1)+1;
        last = 3*(i_time-1)+3;
        A(first:last, 1:3) = R_cluster_one_to_world;
        A(first:last, 4:6) = -R_cluster_two_to_world;
        b(first:last) = p_cluster_two_to_world - p_cluster_one_to_world; 
    end

    c = A\b;
    
    if degreesOfFreedom == 1
    % play around with the null space
        [~, ~, V] = svd(A);
        nullspace = V(:, end);
        aor_cluster_one = nullspace(1:3); aor_cluster_one = aor_cluster_one * 1 / norm(aor_cluster_one);
        aor_cluster_two = nullspace(4:6); aor_cluster_two = aor_cluster_two * 1 / norm(aor_cluster_two);
        AoR = (aor_cluster_one + aor_cluster_two) * 0.5; AoR = AoR * 1 / norm(AoR);        
    end

    cor_cluster_one = c(1:3);
    cor_cluster_two = c(4:6);
    CoR_one_world = cor_cluster_one + cluster_one_reference_center_world;
    CoR_two_world = cor_cluster_two + cluster_two_reference_center_world;
    
    CoR = mean([CoR_one_world CoR_two_world], 2);

    check = A*c - b;
    check_reorg = reshape(check, 3, number_of_time_frames);
    check_rms = (mean(sum(check_reorg.^2, 1))).^(0.5);

    
    
    
    meanError = check_rms / 3;

    % todo: this is just what I think is probably the mean error. Check
    % that by actually calculating the mean distance to the moving CoR


%     figure; axes; hold on
%     for i_marker = 1 : number_of_markers_one
%         plot3(markerClusterOneTrajectory(:, (i_marker-1)*3+1), markerClusterOneTrajectory(:, (i_marker-1)*3+2), markerClusterOneTrajectory(:, (i_marker-1)*3+3));
%     end
%     for i_marker = 1 : number_of_markers_two
%         plot3(markerClusterTwoTrajectory(:, (i_marker-1)*3+1), markerClusterTwoTrajectory(:, (i_marker-1)*3+2), markerClusterTwoTrajectory(:, (i_marker-1)*3+3));
%     end
%     xlabel('x'); ylabel('y'), zlabel('z');



end