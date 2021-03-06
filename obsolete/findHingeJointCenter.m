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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% compare the kinematic tree against the kinematic chain



% function hinge_joint_center = findHingeJointCenter(point_on_link_one, point_on_link_two, point_on_axis, point_on_axis_to_center_distance, direction)
function hinge_joint_center = findHingeJointCenter(point_on_link_one, point_on_link_two, point_on_axis, point_on_axis_to_center_distance, r_A, r_B, direction)
    % finds the center of a hinge joint, given a point on each adjacent link, a point on the axis of rotation and the
    % distance from that point to the center

    % let the link joint angle at which the two points are closest to each other be 0. Then there is a unique intersection
    % of the line between the two points and the joint axis (proof?). The center is defined as that point
    
    % direction specifies an assumption about which of two candidate points to pick by giving the direction of rotation 
    % around the AB-axis to get from M to the joint center. Options are "positive" or "negative".
    
    % for the right arm, given by shoulder, wrist, lateral epicondyle, this should be "positive"
    % for the left arm, given by shoulder, wrist, lateral epicondyle, this should be "negative"

    % notes from 12.2.2017
    % the following should be possible in theory
    % 1. define axes of rotation for the hip where the z-axis lies along the hipcor --> kneemarker direction
    % 2. find angles for xy-rotation in this frame so that the knee marker fits the data
    % 3. find angles for z-rotation in this frame and knee flexion to fit the ankle joint center
    %    problem: I don't know the knee flexion axis - does this matter?
    % 4. this should allow me to calculate the knee joint center position
    
    try_this_method = false;
    if try_this_method
        A = point_on_link_one;
        B = point_on_link_two;
        M = point_on_axis;
        
        % define axes of rotation
        omega_3 = normVector(point_on_link_one - point_on_axis);
        omega_1_prime = [1; 0; 0];
        omega_2_prime = [0; 1; 0];
        hip_rotation_basis_prime = orthogonalizeBasis([omega_3 omega_1_prime omega_2_prime]);
        omega_1 = hip_rotation_basis_prime(:, 2);
        omega_2 = hip_rotation_basis_prime(:, 3);
        omega_3 = hip_rotation_basis_prime(:, 1);
        hip_rotation_basis = [omega_1 omega_2 omega_3];
        
        % find angles for rotation between these two frames
        hip_angles = eulerAnglesFromRotationMatrixZXY(hip_rotation_basis);
        
        
        % calculate left hip flexion-extension
        left_hip_cor_to_left_knee_cor_current_world = left_knee_cor_current - left_hip_cor_current;
        R_world_to_7 = left_hip_direction_matrix_inverse * pelvis_joint_rotation^(-1);
        left_hip_cor_to_left_knee_cor_current_7 = R_world_to_7 * left_hip_cor_to_left_knee_cor_current_world;
        theta_7 = atan2(left_hip_cor_to_left_knee_cor_current_7(2), -left_hip_cor_to_left_knee_cor_current_7(3))
        joint_angles(7) = theta_7;

        % left hip ab-adduction
        R_7 = expAxis(left_hip_direction_matrix(:, 1), joint_angles(7));
        R_world_to_8 = left_hip_direction_matrix_inverse * R_7^(-1) * pelvis_joint_rotation^(-1);
        left_hip_cor_to_left_knee_cor_current_8 = R_world_to_8 * left_hip_cor_to_left_knee_cor_current_world;
        theta_8 = atan2(-left_hip_cor_to_left_knee_cor_current_8(1), -left_hip_cor_to_left_knee_cor_current_8(3));
        joint_angles(8) = theta_8;
        
        
    end
    
    % but the following should also work
    % 1. the position of hip and ankle joint centers should already be sufficient to give me the hip and knee flexion angles
    % not really, because I don't have a reference. All I can get is the knee flexion angle. Does that help me already?
    % and no, I can't even get the knee angle, just 'some' angle, so forget this
    
    
    A = point_on_link_one;
    B = point_on_link_two;
    M = point_on_axis;
    
    if nargin < 5
        direction = 'positive';
    end
    
    
    if false


        % find the point P on the line segment AB that lies on a plane normal to AB intersecting M
        p = A;
        u1 = B - p;
        u2 = M - p;
        u3 = cross(u1, u2);
        A_qr = [u1 u2 u3];

        Q = zeros(3, 3);
        R_qr = zeros(size(3, 3));
        for i_col = 1 : 3
            v = A_qr(:, i_col);
            for i_row = 1 : i_col-1
                R_qr(i_row, i_col) = Q(:, i_row)' * A_qr(:, i_col);
                v = v - R_qr(i_row, i_col) * Q(:, i_row);
            end
            R_qr(i_col, i_col) = norm(v);
            Q(:, i_col) = v/R_qr(i_col, i_col);
        end
        R = Q;
        T_world_to_new = inv([R p; [0 0 0 1]]); % T transforms from world coordinates into my new frame

        A_prime = T_world_to_new * [A; 1];
        B_prime = T_world_to_new * [B; 1];
        M_prime = eye(3, 4) * T_world_to_new * [M; 1];
        P_prime = [M_prime(1); 0; 0];
        P = eye(3, 4) * T_world_to_new^(-1) * [P_prime; 1];

        % find the angle between the plane through ABM and ABE
        E_to_M = point_on_axis_to_center_distance;
        P_to_M = norm(P-M);

        sin_alpha = E_to_M/P_to_M;
        if abs(sin_alpha) < 1
            alpha = asin(sin_alpha);
        else
            if sin_alpha > 1
                alpha = pi/2;
            else
                alpha = -pi/2;
            end
        end

            P_to_E = cos(alpha) * P_to_M;

        alpha_deg = rad2deg(alpha);

        % rotate M back around AB by that angle to find the joint center
        xi = generateTwistCoordinates(P, normVector(B-A));
        if strcmp(direction, 'positive')
            M_rotated = eye(3, 4) * expTwist(xi, alpha) * [M,; 1];
        elseif strcmp(direction, 'negative')
            M_rotated = eye(3, 4) * expTwist(xi, -alpha) * [M,; 1];
        else
            error('direction must be "positive" or "negative"');
        end

        % M_rotated lies on the vector from P to E
        E = P + P_to_E * normVector(M_rotated - P);
        hinge_joint_center = E;
    
    end
   
    
    
    
    
    
    % alternative: just make the simplified assumption that the joint axis is perpendicular to the plane through
    % the three given points
    if false
        joint_axis = normVector(cross(A-M, B-M));
        if strcmp(direction, 'positive')
            hinge_joint_center = point_on_axis + point_on_axis_to_center_distance * joint_axis;
        elseif strcmp(direction, 'negative')
            hinge_joint_center = point_on_axis - point_on_axis_to_center_distance * joint_axis;
        else
            error('direction must be "positive" or "negative"');
        end
    end
    
    
    % alternative - this actually works. See pictures taken on Nov. 26 for illustration.
    
    % we know that E lies on a sphere of radius r_A around A, and on a sphere of radius r_B around B
    % these spheres intersect in a circle 
    % the center of that circle is on the line between A and B
    % the radius of that circle is r_C
    
    % find r_C, given the sides of the triangle a, b, c
    r_C = triangleHeight(norm(B-A), r_A, r_B);
    
    alpha = asin(r_C / r_A)
    beta = asin(r_C / r_B)
    gamma = pi - alpha - beta
    
    % find C
    AC_norm = sqrt(r_A^2 - r_C^2);
    if r_A^2 - r_C^2 < 0
        warning('distance between provided points is larger than possible')
    end
    C_from_A = A + AC_norm * normVector(B-A);
    BC_norm = sqrt(r_B^2 - r_C^2);
    if r_B^2 - r_C^2 < 0
        warning('distance between provided points is larger than possible')
    end
    C_from_B = B + BC_norm * normVector(A-B);
    C = mean([C_from_A C_from_B], 2);
    
    % switch to coordinates with C as the origin
    p = C;
    u3 = B - p;
    u1 = M - p;
    u2 = cross(u3, u1);
    A_qr = [u1 u2 u3];

    Q = zeros(3, 3);
    R_qr = zeros(size(3, 3));
    for i_col = 1 : 3
        v = A_qr(:, i_col);
        for i_row = 1 : i_col-1
            R_qr(i_row, i_col) = Q(:, i_row)' * A_qr(:, i_col);
            v = v - R_qr(i_row, i_col) * Q(:, i_row);
        end
        R_qr(i_col, i_col) = norm(v);
        Q(:, i_col) = v/R_qr(i_col, i_col);
    end
    R = Q;
    T_new_to_world = ([R p; [0 0 0 1]]); % T transforms from world coordinates into my new frame
    if any(any(isnan(T_new_to_world)))
        hinge_joint_center = zeros(3, 1) * NaN;
    else
        
        T_world_to_new = T_new_to_world^(-1); % T transforms from world coordinates into my new frame



        C_prime = eye(3, 4) * T_world_to_new * [C; 1];
        M_prime = eye(3, 4) * T_world_to_new * [M; 1];
        A_prime = eye(3, 4) * T_world_to_new * [A; 1];

        r_M = point_on_axis_to_center_distance;

        h = triangleHeight(norm(M-C), r_C, r_M);

        % find H
        CH_norm = sqrt(r_C^2 - h^2);
        if r_C^2 - h^2 < 0
            warning('distance between joint center candidates and marker is larger than possible')
        end
        H_from_C = C_prime + CH_norm * normVector(M_prime-C_prime);
        MH_norm = sqrt(r_M^2 - h^2);
        if r_M^2 - h^2 < 0
            warning('distance between joint center candidates and marker is larger than possible')
        end
        H_from_M = M_prime + MH_norm * normVector(C_prime-M_prime);
        H_prime = mean([H_from_C H_from_M], 2);

        E_1_prime = [H_prime(1); h; 0];
        E_2_prime = [H_prime(1); -h; 0];

        E_1 = eye(3, 4) * T_new_to_world * [E_1_prime; 1]
        E_2 = eye(3, 4) * T_new_to_world * [E_2_prime; 1]


        % do E_1 and E_2 actually lie on the appropriate spheres?
        r_A;
        r_B;
        r_M;

        r_A_check_1 = norm(A - E_1);
        r_A_check_2 = norm(A - E_2);
        r_B_check_1 = norm(B - E_1);
        r_B_check_2 = norm(B - E_2);
        r_M_check_1 = norm(M - E_1);
        r_M_check_2 = norm(M - E_2);

        % test relationships
        BC = C - B;
        CE_1 = E_1 - C;
        CE_2 = E_2 - C;
        E_1M = M - E_1;
        E_2M = M - E_2;
        check_1 = dot(cross(BC, CE_1), E_1M); % this should always be > 0
        check_2 = dot(cross(BC, CE_2), E_2M); % this should always be < 0

        % now figure out which point to choose
        % we have two options, depending on where the marker is located relative to the plane spanned by A, B and E
        % C is the point on the line segment AB closest to E
        % if the cross product BC x CE points in the general direction of M, choose the upper triangle, i.e. E_1
        % if the cross product BC x CE points away from of M, choose the lower triangle, i.e. E_1

        % for the left arm with A = shoulder, B = wrist and E = elbow, M = lateral humeral epicondyle, this should be negative
        % for the right arm with A = shoulder, B = wrist and E = elbow, M = lateral humeral epicondyle, this should be positive

        % for the left leg with A = hip, B = ankle and E = knee, M = lateral femoral epicondyle, this should be positive
        % for the right leg with A = hip, B = ankle and E = knee, M = lateral femoral epicondyle, this should be negative

        % note: this does not work satisfactorily for small angles in the hinge joint when the joint centers are not estimated accurately

        if strcmp(direction, 'positive')
            hinge_joint_center = E_1;
        elseif strcmp(direction, 'negative')
            hinge_joint_center = E_2;
        else
            error('direction must be "positive" or "negative"');
        end
    end
end






