
% script to explore the step changes

% load an exemplary configuration 
% check out what happens to the foot when rotating the knee/hip 
% maybe calculate the required torques / work to do that

visualize = 1;

% define data to load
data_root = '/Users/reimajbi/Drive_UD/Vision_HY';
subject = 'DJB';
condition = 'walking';

% define condition to check
stance_foot_to_check = 'STANCE_RIGHT';
index_to_check = 'CONTROL';
example_to_use = 1;

% load model and data
% load([data_root filesep subject filesep 'subjectInfo.mat'], 'date', 'subject_id')
% load([data_root filesep subject filesep 'subjectModel.mat'])
% load([data_root filesep subject filesep 'analysis' filesep makeFileName(date, subject_id, 'results')])

% find desired condition
stance_foot_indicator = strcmp(condition_stance_foot_list_session, stance_foot_to_check);
index_indicator = strcmp(condition_index_list_session, index_to_check);
check_condition_indicator = stance_foot_indicator & index_indicator;
check_condition_indices = find(check_condition_indicator);
check_condition_example = check_condition_indices(example_to_use);
origin_trial = origin_trial_list_session(check_condition_example);
origin_stretch_start = origin_start_time_list_session(check_condition_example);
origin_stretch_end = origin_end_time_list_session(check_condition_example);

% load data
% load([data_root filesep subject filesep 'processed' filesep makeFileName(date, subject_id, condition, origin_trial, 'kinematicTrajectories')]);

% extract stretch end and surrounding time
[~, index_now] = min(abs(time_mocap - origin_stretch_end));
theta = joint_angle_trajectories_optimized(index_now, :)';
kinematic_tree.jointAngles = theta;
kinematic_tree.updateConfiguration();

left_heel_position_now = kinematic_tree.endEffectorPositions{strcmp(kinematic_tree.endEffectorLabels, 'left heel')};

% rotate left hip and look at effect on lateral left heel position
observed_left_hip_rotation_delta = deg2rad(1);
theta_left_hip_rotation = theta;
theta_left_hip_rotation(strcmp(kinematic_tree.jointLabels, 'left hip internal/external rotation')) = theta_left_hip_rotation(strcmp(kinematic_tree.jointLabels, 'left hip internal/external rotation')) + observed_left_hip_rotation_delta;
kinematic_tree.jointAngles = theta_left_hip_rotation;
kinematic_tree.updateConfiguration();
left_heel_position_left_hip_rotation = kinematic_tree.endEffectorPositions{strcmp(kinematic_tree.endEffectorLabels, 'left heel')};
delta_left_heel_by_delta_left_hip_rotation = (left_heel_position_left_hip_rotation - left_heel_position_now);

% rotate right knee and look at effect on lateral left heel position
observed_right_knee_rotation_delta = -deg2rad(0.5);
J_right_foot = kinematic_tree.bodyJacobians{strcmp(kinematic_tree.endEffectorLabels, 'right heel')};
J_right_leg = [zeros(7, 13) eye(7) zeros(7, 18)];
theta_change = pinv([J_right_foot; J_right_leg]) * [zeros(6, 1); [0 0 0 0 observed_right_knee_rotation_delta 0 0]'];
theta_right_knee_rotation = theta + theta_change;
kinematic_tree.jointAngles = theta_right_knee_rotation;
kinematic_tree.updateConfiguration();
left_heel_position_right_knee_rotation = kinematic_tree.endEffectorPositions{strcmp(kinematic_tree.endEffectorLabels, 'left heel')};
delta_left_heel_by_delta_right_knee_rotation = (left_heel_position_right_knee_rotation - left_heel_position_now);

delta_left_heel_by_delta_both_rotations = delta_left_heel_by_delta_left_hip_rotation + delta_left_heel_by_delta_right_knee_rotation;
all_rotation_deltas_mm = [delta_left_heel_by_delta_left_hip_rotation delta_left_heel_by_delta_right_knee_rotation delta_left_heel_by_delta_both_rotations] * 1000;




% estimate differentials
delta_theta = 0.001;
theta_left_hip_abduction = theta;
theta_left_hip_abduction(strcmp(kinematic_tree.jointLabels, 'left hip ab/adduction')) = theta_left_hip_abduction(strcmp(kinematic_tree.jointLabels, 'left hip ab/adduction')) + delta_theta;
kinematic_tree.jointAngles = theta_left_hip_abduction;
kinematic_tree.updateConfiguration();
left_heel_position_left_hip_abduction = kinematic_tree.endEffectorPositions{strcmp(kinematic_tree.endEffectorLabels, 'left heel')};
d_left_heel_by_d_left_hip_abduction = (left_heel_position_left_hip_abduction - left_heel_position_now) * 1 / delta_theta;

theta_left_hip_rotation = theta;
theta_left_hip_rotation(strcmp(kinematic_tree.jointLabels, 'left hip internal/external rotation')) = theta_left_hip_rotation(strcmp(kinematic_tree.jointLabels, 'left hip internal/external rotation')) + delta_theta;
kinematic_tree.jointAngles = theta_left_hip_rotation;
kinematic_tree.updateConfiguration();
left_heel_position_left_hip_rotation = kinematic_tree.endEffectorPositions{strcmp(kinematic_tree.endEffectorLabels, 'left heel')};
d_left_heel_by_d_left_hip_rotation = (left_heel_position_left_hip_rotation - left_heel_position_now) * 1 / delta_theta;

theta_right_knee_rotation = theta;
theta_right_knee_rotation(strcmp(kinematic_tree.jointLabels, 'right knee external/internal rotation')) = theta_right_knee_rotation(strcmp(kinematic_tree.jointLabels, 'right knee external/internal rotation')) + delta_theta;
kinematic_tree.jointAngles = theta_right_knee_rotation;
kinematic_tree.updateConfiguration();
left_heel_position_right_knee_rotation = kinematic_tree.endEffectorPositions{strcmp(kinematic_tree.endEffectorLabels, 'left heel')};
d_left_heel_by_d_right_knee_rotation = (left_heel_position_right_knee_rotation - left_heel_position_now) * 1 / delta_theta;








if visualize
    % show stick figure of geometric model
    pelvis_center_pos = kinematic_tree.jointPositions{6};
    scene_bound = repmat(pelvis_center_pos, 1, 2) + 2*[-0.5 0.5; -0.5 0.5; -1 1];
%     stick_figure = KinematicTreeController(kinematic_tree, scene_bound, 'ellipsoid');
        stick_figure = KinematicTreeController(kinematic_tree, scene_bound, 'none');


end






return

























