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

% show population results

plot_detailed       = 1;
plot_overview       = 0;
plot_episodes       = 0;
include_control     = 0;

show_legend         = 1;
dictate_axes        = 0;
mark_pushoff        = 0;

save_figures        = 0;

error_shades = 'cinv';
% error_shades = 'std';

% define subjects
subjects = {'DXT', 'EFU', 'GHJ', 'RON', 'RRB', 'YMU'};
% subjects = {'DXT', 'EFU', 'RON', 'RRB', 'YMU'};
subjects = {'DXT'};
% subjects = {'RON'};
% subjects = {'BRC', 'RTZ', 'XDQ', 'XEA'};
% subjects = {'BRC'};
% subjects = {'RTZ'};
% subjects = {'XDQ'};
% subjects = {'XEA'};

% subjects = {'CVX'};
% subjects = {'WAU'};
% subjects = {'YPQ'};
% subjects = {'CVX', 'WAU', 'YPQ'};
% subjects = {'JXG'};

% subjects = {'GGU'};
% subjects = {'XYC'};
% subjects = {'LDZ'};
% subjects = {'STD'};

subjects = {'GGU', 'XYC', 'LDZ', 'STD'};

subjects = {'C'};
subjects = {'DDPilot'};

%% choose variables to plot
% variable info contains the following columns
% variable name | display name | y-label with unit | save label
continuous_variable_info = {};
discrete_variable_info = {};

% markers
% the cell array should have the following entries in each line: variable name, variable label, unit, file label for saving, forced axis scale, positive direction label, negative direction label
continuous_variable_info = [continuous_variable_info; {'lheel_x_pos_normalized_all', 'left heel pos, ml', 'heel pos (m)', 'lheelpos', [-0.2 0.2], 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'rheel_x_pos_normalized_all', 'right heel pos, ml', 'heel pos (m)', 'rheelpos', 0, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'lheel_x_pos_com_normalized_all', 'left heel pos, ml', 'heel pos (m)', 'lheelpos', 0, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'rheel_x_pos_com_normalized_all', 'right heel pos, ml', 'heel pos (m)', 'rheelpos', 0, 'right', 'left'}];

% continuous_variable_info = [continuous_variable_info; {'trunk_angle_ml_normalized_all', 'trunk angle, ml', 'angle (deg)', 'trunkangleml', [-10 16], 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'lleg_angle_ml_normalized_all', 'left leg angle, ml', 'angle (deg)', 'llegangleml', 0, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'rleg_angle_ml_normalized_all', 'right leg angle, ml', 'angle (deg)', 'rlegangleml', 0, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'trunk_angle_ap_normalized_all', 'trunk angle, ap', 'angle (deg)', 'trunkangleml', [-5 25], 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'com_x_pos_normalized_all', 'CoM, ml', 'com (m)', 'comml', [-5 5], 'left', 'right'}];

% continuous_variable_info = [continuous_variable_info; {'lheel_x_pos_response', 'left heel pos response, ml', 'heel pos (m)', 'lheelposRsp', [-0.05 0.05], 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'rheel_x_pos_response', 'right heel pos response, ml', 'heel pos (m)', 'rheelposRsp', [-0.05 0.05], 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'lheel_x_pos_mpsis_response', 'left heel pos response, ml, rel. to Mpsis', 'heel pos (m)', 'lheelposRsp', [-0.03 0.03], 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'rheel_x_pos_mpsis_response', 'right heel pos response, ml, rel. to Mpsis', 'heel pos (m)', 'rheelposRsp', [-0.03 0.03], 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'lheel_x_pos_com_response', 'left heel pos response, ml, rel. to CoM', 'heel pos (m)', 'lheelposRsp', [-0.02 0.02], 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'rheel_x_pos_com_response', 'right heel pos response, ml, rel. to CoM', 'heel pos (m)', 'rheelposRsp', [-0.02 0.02], 'right', 'left'}];

% continuous_variable_info = [continuous_variable_info; {'trunk_angle_ml_response', 'trunk angle response, ml', 'angle (deg)', 'trunkanglemlRsp', [-5 5], 'cw', 'c-cw'}];
% continuous_variable_info = [continuous_variable_info; {'lleg_angle_ml_response', 'left leg angle response, ml', 'angle (deg)', 'lleganglemlRsp', [-5 5], 'cw', 'c-cw'}];
% continuous_variable_info = [continuous_variable_info; {'rleg_angle_ml_response', 'right leg angle response, ml', 'angle (deg)', 'rleganglemlRsp', [-5 5], 'cw', 'c-cw'}];
% continuous_variable_info = [continuous_variable_info; {'com_x_pos_response', 'CoM response, ml', 'com (m)', 'commlRsp', [-0.1 0.1], 'left', 'right'}];

% forceplate
% continuous_variable_info = [continuous_variable_info; {'cop_x_normalized_all', 'total CoP, ml', 'CoP (m)', 'copx', 0, 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'cop_x_response', 'total CoP response, ml', 'CoP (m)', 'copxRsp', [-0.02 0.02], 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'cop_x_stancefoot_response', 'total CoP response, ml, rel. to stance foot', 'CoP (m)', 'copxRspStancefoot', [-0.02 0.02], 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'cop_x_mpsis_response', 'total CoP response, ml, rel. to MPSIS', 'CoP (m)', 'copxRspMpsis', [-0.004 0.008], 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'cop_x_com_response', 'total CoP response, ml, rel. to CoM', 'CoP (m)', 'copxRspMpsis', [-0.09 0.09], 'right', 'left'}];
% continuous_variable_info = [continuous_variable_info; {'f_x_normalized_all', 'total force, ml', 'f (N)', 'fx', [-120 120], '?', '?'}];
% continuous_variable_info = [continuous_variable_info; {'f_x_response', 'total force response, ml', 'f (N)', 'fxRsp', [-30 30], '?', '?'}];
% continuous_variable_info = [continuous_variable_info; {'f_z_normalized_all', 'total vertical force', 'f (N)', 'fz', 0, '?', '?'}];
% continuous_variable_info = [continuous_variable_info; {'f_z_response', 'total force vertical response', 'f (N)', 'fzRsp', 0.025, '?', '?'}];
% continuous_variable_info = [continuous_variable_info; {'m_y_normalized_all', 'total moment, ml', 'm (Nm)', 'my', 0, '?', '?'}];
% continuous_variable_info = [continuous_variable_info; {'m_y_response', 'total moment response, ml', 'm (Nm)', 'myRsp', 100, '?', '?'}];

% EMG
% continuous_variable_info = [continuous_variable_info; {'lglutmed_normalized_all', 'left Gluteus Medius', 'EMG', 'lglutmed', [0 0], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'ltibiant_normalized_all', 'left Tibialis Anterior', 'EMG', 'ltibiant', [0 0], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'lgastroc_normalized_all', 'left Gastrocnemius Medialis', 'EMG', 'lgastroc', [0 0], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'lperolng_normalized_all', 'left Peroneus Longus', 'EMG', 'lperolng', [0 0], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'ltnsrflt_normalized_all', 'left TFL', 'EMG', 'ltnsrflt', [0 0], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'rglutmed_normalized_all', 'right Gluteus Medius', 'EMG', 'rglutmed', [0 0], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'rtibiant_normalized_all', 'right Tibialis Anterior', 'EMG', 'rtibiant', [0 0], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'rgastroc_normalized_all', 'right Gastrocnemius Medialis', 'EMG', 'rgastroc', [0 0], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'rperolng_normalized_all', 'right Peroneus Longus', 'EMG', 'rperolng', [0 0], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'rtnsrflt_normalized_all', 'right TFL', 'EMG', 'rtnsrflt', [0 0], '+', '-'}];

% continuous_variable_info = [continuous_variable_info; {'lglutmed_rescaled_response', 'left Gluteus Medius', 'EMG', 'lglutmedRescaledRsp', [-2.5 2.5], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'ltibiant_rescaled_response', 'left Tibialis Anterior', 'EMG', 'ltibiantRescaledRsp', [-2.5 2.5], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'lgastroc_rescaled_response', 'left Gastrocnemius Medialis', 'EMG', 'lgastrocRescaledRsp', [-2.5 2.5], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'lperolng_rescaled_response', 'left Peroneus Longus', 'EMG', 'lperolngRescaledRsp', [-2.5 2.5], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'rglutmed_rescaled_response', 'right Gluteus Medius', 'EMG', 'rglutmedRescaledRsp', [-2.5 2.5], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'rtibiant_rescaled_response', 'right Tibialis Anterior', 'EMG', 'rtibiantRescaledRsp', [-2.5 2.5], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'rgastroc_rescaled_response', 'right Gastrocnemius Medialis', 'EMG', 'rgastrocRescaledRsp',[-2.5 2.5], '+', '-'}];
% continuous_variable_info = [continuous_variable_info; {'rperolng_rescaled_response', 'right Peroneus Longus', 'EMG', 'rperolngRescaledRsp', [-2.5 2.5], '+', '-'}];



% the cell array should have the following entries in each line: variable name, variable label, unit, file label for saving, forced axis scale, positive direction label, negative direction label
% discrete_variable_info = [discrete_variable_info; {'step_times_all', 'step times', '(s)', 'steptimes', [0.41 0.69], '', ''}];
% discrete_variable_info = [discrete_variable_info; {'step_length_all', 'step length', '(m)', 'steplength', [0.5 0.85], 'front', 'back'}];
% discrete_variable_info = [discrete_variable_info; {'step_width_all', 'step width', '(m)', 'stepwidth', [0.0 0.25], 'right', 'left'}];
% discrete_variable_info = [discrete_variable_info; {'step_speed_all', 'step speed', '(m/s)', 'stepspeed', [0.8 2], '', ''}];
% discrete_variable_info = [discrete_variable_info; {'cadence_all', 'cadence', '(steps/min)', 'cadence', [80 150], '', ''}];
% discrete_variable_info = [discrete_variable_info; {'foot_placement_mpsis_all', 'foot placement rel. to MPSIS', '(m)', 'stepmpsis', [0.0 0.25], 'right', 'left'}];
% discrete_variable_info = [discrete_variable_info; {'step_width_response', 'step width response', '(m)', 'stepwidthrsp', [0.0 0.25], 'right', 'left'}];
% discrete_variable_info = [discrete_variable_info; {'foot_placement_world_response', 'foot placement response', '(m)', 'stepworldrsp', [-0.14 0.14], 'right', 'left'}];
% discrete_variable_info = [discrete_variable_info; {'foot_placement_stancefoot_response', 'foot placement rel. to stancefoot response', '(m)', 'stepstancefootrsp', [-0.14 0.14], 'right', 'left'}];
% discrete_variable_info = [discrete_variable_info; {'foot_placement_mpsis_response', 'foot placement rel. to MPSIS response', '(m)', 'stepmpsisrsp', [-0.14 0.14], 'right', 'left'}];


%% choose conditions to plot
condition_labels = {'stance foot', 'perturbation', 'delay', 'index', 'experimental'};
condition_column_index = find(strcmp(condition_labels, 'index'));
condition_column_stancefoot = find(strcmp(condition_labels, 'stance foot'));

% conditions_control = ...
%   {
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'walking'; ...
%     'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'CONTROL', 'walking'; ...
%   };
% 
% % for vision
% conditions_to_plot = ...
%   {
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'ONE', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'ONE', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'ONE', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'ONE', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'TWO', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'TWO', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'TWO', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'TWO', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'THREE', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'THREE', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'THREE', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'THREE', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'FOUR', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'FOUR', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'FOUR', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'FOUR', 'walking'; ...
%   };
% comparison_to_make = 2; % perturbation only
condition_column_index = find(strcmp(condition_labels, 'index'));
condition_column_stancefoot = find(strcmp(condition_labels, 'stance foot'));

% for phase-dependent GVS
conditions_control = ...
  {
    'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'walking'; ...
    'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'CONTROL', 'walking'; ...
  };

conditions_to_plot = ...
  {
    'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'ONE', 'walking'; ...
    'STANCE_RIGHT', 'ILLUSION_RIGHT', '150ms', 'ONE', 'walking'; ...
    'STANCE_RIGHT', 'ILLUSION_RIGHT', '450ms', 'ONE', 'walking'; ...
    'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'ONE', 'walking'; ...
    'STANCE_RIGHT', 'ILLUSION_LEFT', '150ms', 'ONE', 'walking'; ...
    'STANCE_RIGHT', 'ILLUSION_LEFT', '450ms', 'ONE', 'walking'; ...
    'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'TWO', 'walking'; ...
    'STANCE_LEFT', 'ILLUSION_RIGHT', '150ms', 'TWO', 'walking'; ...
    'STANCE_LEFT', 'ILLUSION_RIGHT', '450ms', 'TWO', 'walking'; ...
    'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'TWO', 'walking'; ...
    'STANCE_LEFT', 'ILLUSION_LEFT', '150ms', 'TWO', 'walking'; ...
    'STANCE_LEFT', 'ILLUSION_LEFT', '450ms', 'TWO', 'walking'; ...
    'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'THREE', 'walking'; ...
    'STANCE_RIGHT', 'ILLUSION_RIGHT', '150ms', 'THREE', 'walking'; ...
    'STANCE_RIGHT', 'ILLUSION_RIGHT', '450ms', 'THREE', 'walking'; ...
    'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'THREE', 'walking'; ...
    'STANCE_RIGHT', 'ILLUSION_LEFT', '150ms', 'THREE', 'walking'; ...
    'STANCE_RIGHT', 'ILLUSION_LEFT', '450ms', 'THREE', 'walking'; ...
    'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'FOUR', 'walking'; ...
    'STANCE_LEFT', 'ILLUSION_RIGHT', '150ms', 'FOUR', 'walking'; ...
    'STANCE_LEFT', 'ILLUSION_RIGHT', '450ms', 'FOUR', 'walking'; ...
    'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'FOUR', 'walking'; ...
    'STANCE_LEFT', 'ILLUSION_LEFT', '150ms', 'FOUR', 'walking'; ...
    'STANCE_LEFT', 'ILLUSION_LEFT', '450ms', 'FOUR', 'walking'; ...
  };
% comparison_to_make = 2; % perturbation only
comparison_to_make = 3; % delay only



% for phase-dependent GVS
% conditions_control = ...
%   {
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'walking'; ...
%     'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'CONTROL', 'walking'; ...
%   };

% % first step left stance

% % first step left stance
% conditions_to_plot = ...
%   {
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'ONE', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '150ms', 'ONE', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '450ms', 'ONE', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'ONE', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '150ms', 'ONE', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '450ms', 'ONE', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'TWO', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '150ms', 'TWO', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '450ms', 'TWO', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'TWO', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '150ms', 'TWO', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '450ms', 'TWO', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'THREE', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '150ms', 'THREE', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '450ms', 'THREE', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'THREE', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '150ms', 'THREE', 'walking'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '450ms', 'THREE', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'FOUR', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '150ms', 'FOUR', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '450ms', 'FOUR', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'FOUR', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '150ms', 'FOUR', 'walking'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '450ms', 'FOUR', 'walking'; ...
%   };
% comparison_to_make = 2; % perturbation only% conditions_to_plot = ...
%   {
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'ONE'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'ONE'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'TWO'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'TWO'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'THREE'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'THREE'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'FOUR'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'FOUR'; ...
%   };






% first step left stance

% first step right stance
% conditions_to_plot = ...
%   {
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'ONE'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'ONE'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'TWO'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'TWO'; ...
%     'STANCE_RIGHT', 'ILLUSION_RIGHT', '0ms', 'THREE'; ...
%     'STANCE_RIGHT', 'ILLUSION_LEFT', '0ms', 'THREE'; ...
%     'STANCE_LEFT', 'ILLUSION_RIGHT', '0ms', 'FOUR'; ...
%     'STANCE_LEFT', 'ILLUSION_LEFT', '0ms', 'FOUR'; ...
%   };

% % for ArmSense
% conditions_control = {};
% conditions_to_plot = ...
%   {
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'baselineOG'; ...
%     'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'CONTROL', 'baselineOG'; ...
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'baselineTM'; ...
%     'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'CONTROL', 'baselineTM'; ...
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'feedback'; ...
%     'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'CONTROL', 'feedback'; ...
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'postTM'; ...
%     'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'CONTROL', 'postTM'; ...
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'postOG'; ...
%     'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'CONTROL', 'postOG'; ...
%   };
% conditions_to_plot = ...
%   {
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'baselineOG'; ...
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'baselineTM'; ...
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'feedback'; ...
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'postTM'; ...
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'postOG'; ...
%   };
% comparison_to_make = 5; % experimental

% for Obstacle
% conditions_control = {};
% conditions_to_plot = ...
%   {
%     'STANCE_BOTH', 'CONTROL', 'CONTROL', 'ONE', 'NEAR'; ...
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'TWO', 'NEAR'; ...
%     'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'THREE', 'NEAR'; ...
%     'STANCE_BOTH', 'CONTROL', 'CONTROL', 'ONE', 'FAR'; ...
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'TWO', 'FAR'; ...
%     'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'THREE', 'FAR'; ...
%     'STANCE_BOTH', 'CONTROL', 'CONTROL', 'ONE', 'NO'; ...
%     'STANCE_LEFT', 'CONTROL', 'CONTROL', 'TWO', 'NO'; ...
%     'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'THREE', 'NO'; ...
%   };
% comparison_to_make = 5; % experimental

% for subconcussion pilot
conditions_control = {};
conditions_to_plot = ...
  {
    'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'MP'; ...
    'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'CONTROL', 'MP'; ...
    'STANCE_LEFT', 'CONTROL', 'CONTROL', 'CONTROL', 'None'; ...
    'STANCE_RIGHT', 'CONTROL', 'CONTROL', 'CONTROL', 'None'; ...
  };
comparison_to_make = 5; % experimental






number_of_conditions_control = size(conditions_control, 1);
number_of_conditions_to_plot = size(conditions_to_plot, 1);

%% define comparisons and episodes

% here we parse the list of conditions to plot and extract groups of conditions that should be compared against each
% other, according to the specified comparison to make

% comparison_to_make = 3; % delay only
use_control = ~isempty(conditions_control);
comparison_indices = {};
control_conditions_for_comparisons = [];
conditions_already_compared = [];
while length(conditions_already_compared) < number_of_conditions_to_plot
    % start with the first available condition
    i_condition = 1;
    while ismember(i_condition, conditions_already_compared)
        i_condition = i_condition + 1;
    end
    
    this_comparison = i_condition; % this is the first condition in this episode, more will be added
    % search for conditions that differ from this one only in the one we're comparing
    for j_condition = 1 : number_of_conditions_to_plot
        if i_condition ~= j_condition
            % check which conditions labels agree between these two conditions
            comparison_table = zeros(1, length(condition_labels)); % this is a table indicating equality between the two conditions in questions
            for i_label = 1 : length(condition_labels)
                comparison_table(i_label) = strcmp(conditions_to_plot{i_condition, i_label}, conditions_to_plot{j_condition, i_label});
            end

            % look at the relevant entries of the comparison table
            comparison_table_relevant = comparison_table;
            comparison_table_relevant(comparison_to_make) = [];
            if all(comparison_table_relevant)
                this_comparison = [this_comparison, j_condition];
            end
        end
    end
    comparison_indices = [comparison_indices; this_comparison];
    conditions_already_compared = [conditions_already_compared this_comparison];
end

% now we go through the groups of conditions and make list of indices that form episodes

% define episodes
episode_first_stretch_indices = find(strcmp(conditions_to_plot(:, 4), 'ONE'));
episode_indices = {};
comparisons_already_used = [];
number_of_comparisons = length(comparison_indices);
while length(comparisons_already_used) < number_of_comparisons
    % start with the first available comparison
    i_comparison = 1;
    while ismember(i_comparison, comparisons_already_used)
        i_comparison = i_comparison + 1;
    end
    
    this_episode = i_comparison; % this is the first comparison in this episode, more will be added
    
    % search for comparisons that differ from this one in only the step number
    base_comparison = comparison_indices{i_comparison};
    example_condition_in_base_comparison = base_comparison(1);
    example_condition_in_base_comparison_labels = conditions_to_plot(example_condition_in_base_comparison, :);
    for j_comparison = 1 : number_of_comparisons
        if i_comparison ~= j_comparison
            this_comparison = comparison_indices{j_comparison};
            example_condition_in_this_comparison = this_comparison(1);
            example_condition_in_this_comparison_labels = conditions_to_plot(example_condition_in_this_comparison, :);
            % check which conditions labels agree between these two conditions
            comparison_table = zeros(1, length(condition_labels)); % this is a table indicating equality between the two conditions in questions
            for i_label = 1 : length(condition_labels)
                comparison_table(i_label) = strcmp(example_condition_in_base_comparison_labels{i_label}, example_condition_in_this_comparison_labels{i_label});
            end
            
            % look at the relevant entries of the comparison table
            comparison_table_relevant = comparison_table;
            comparison_table_relevant([condition_column_stancefoot condition_column_index comparison_to_make]) = [];
            if all(comparison_table_relevant)
                % check if the stance foot is alternating
                if strcmp(example_condition_in_base_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                    if strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'TWO') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                        this_episode = [this_episode, j_comparison];
                    elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'THREE') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                        this_episode = [this_episode, j_comparison];
                    elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'FOUR') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                        this_episode = [this_episode, j_comparison];
                    end
                elseif strcmp(example_condition_in_base_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                    if strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'TWO') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                        this_episode = [this_episode, j_comparison];
                    elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'THREE') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_LEFT')
                        this_episode = [this_episode, j_comparison];
                    elseif strcmp(example_condition_in_this_comparison_labels(condition_column_index), 'FOUR') && strcmp(example_condition_in_this_comparison_labels(condition_column_stancefoot), 'STANCE_RIGHT')
                        this_episode = [this_episode, j_comparison];
                    end
                end
            end
        end
    end
    episode_indices = [episode_indices; this_episode];
    comparisons_already_used = [comparisons_already_used this_episode];
    
end


%% collect data from all subjects
step_times_data = [];
pushoff_times_data = [];
condition_stance_foot_data = {};
condition_perturbation_data = {};
condition_delay_data = {};
condition_index_data = {};
condition_experimental_data = {};
continuous_variable_data = cell(size(continuous_variable_info, 1), 1);
discrete_variable_data = cell(size(discrete_variable_info, 1), 1);

for i_subject = 1 : length(subjects)
    % load subject data
    path = strsplit(pwd, filesep);
    if strcmp(path(end), subjects{i_subject})
        data_path = '';
    else
        data_path = [subjects{i_subject} filesep];
    end
    load([data_path 'subjectInfo.mat']);
    load([data_path 'analysis' filesep date '_' subject_id '_resultsConditions.mat']);
    load([data_path 'analysis' filesep date '_' subject_id '_resultsBalance.mat']);
%     load([data_path 'analysis' filesep date '_' subject_id '_resultsArmswing.mat']);
    load([data_path 'analysis' filesep date '_' subject_id '_resultsForceplate.mat']);
%     load([data_path 'analysis' filesep date '_' subject_id '_resultsEmg.mat']);
    
    % time
    step_times_data = [step_times_data step_times_all];
    pushoff_times_data = [pushoff_times_data pushoff_times_all];
    
    % conditions
    condition_stance_foot_data = [condition_stance_foot_data; condition_stance_foot_list_all];
    condition_perturbation_data = [condition_perturbation_data; condition_perturbation_list_all];
    condition_delay_data = [condition_delay_data; condition_delay_list_all];
    condition_index_data = [condition_index_data; condition_index_list_all];
    condition_experimental_data = [condition_experimental_data; condition_experimental_list_all];
    for i_variable = 1 : size(discrete_variable_info, 1)
        evalstring = ['variable_data = ' discrete_variable_info{i_variable, 1} ';'];
        eval(evalstring);
        discrete_variable_data{i_variable} = [discrete_variable_data{i_variable} variable_data];
    end
    for i_variable = 1 : size(continuous_variable_info, 1)
        evalstring = ['variable_data = ' continuous_variable_info{i_variable, 1} ';'];
        eval(evalstring);
        continuous_variable_data{i_variable} = [continuous_variable_data{i_variable} variable_data];
    end


end
step_time_mean = mean(step_times_all);
pushoff_time_mean = mean(pushoff_times_all);
time_normalized = linspace(0, step_time_mean, number_of_time_steps_normalized);

%% do plots
color_control = [0.3 0.1 1];
colors_comparison = ...
  [ ...
    [1 0.3 0.1] * 0.7; ...
    [0.3 1 0.1] * 0.7; ...
    [0.3 0.1 1] * 0.7; ...
  ]; % should have one row per condition in the comparison

colors_comparison = ...
  [ ...
    [1 0.3 0.1] * 0.7; ...
    [0.3 1 0.1] * 0.7; ...
    lightenColor([0.3 0.1 1], 0.1); ...
  ]; % should have one row per condition in the comparison

% colors_comparison = ...
%   [ ...
%     [241 90 34] * 1/255; ...
%     [0 166 81] * 1/255; ...
%     [0 173 220] * 1/255; ...
%   ]; % should have one row per condition in the comparison

% colors_comparison = ...
%   [ ...
%     [191 0 0] * 1/255; ...
%     [64 0 146] * 1/255; ...
%     [255 178 0] * 1/255; ...
%     [1 168 5] * 1/255; ...
%     [0 202 229] * 1/255; ...
%   ]; % should have one row per condition in the comparison

double_stance_color = [0 0 0];
double_stance_alpha = 0.05;

%% plot detailed
if plot_detailed
    % discrete variables
    for i_variable = 1 : size(discrete_variable_info, 1)
        points_to_plot = discrete_variable_data{i_variable, 1};
        for i_comparison = 1 : length(comparison_indices);
            % find correct condition indicator for control
            this_comparison = comparison_indices{i_comparison};
            representant_condition_index = this_comparison(1);
            if strcmp(conditions_to_plot(representant_condition_index, 1), 'STANCE_LEFT')
                applicable_control_condition_indices = 1;
            elseif strcmp(conditions_to_plot(representant_condition_index, 1), 'STANCE_RIGHT')
                applicable_control_condition_indices = 2;
            end
            stance_foot_indicator = strcmp(condition_stance_foot_data, conditions_control(applicable_control_condition_indices, 1));
            perturbation_indicator = strcmp(condition_perturbation_data, conditions_control(applicable_control_condition_indices, 2));
            delay_indicator = strcmp(condition_delay_data, conditions_control(applicable_control_condition_indices, 3));
            index_indicator = strcmp(condition_index_data, conditions_control(applicable_control_condition_indices, 4));
            this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator;
            
            figure; axes; hold on;
%             % plot control
%             histogram ...
%               ( ...
%                 points_to_plot(this_condition_indicator), ...
%                 'edgecolor', color_control, ...
%                 'facecolor', lightenColor(color_control, 0.5), ...
%                 'DisplayName', 'control' ...
%               )
            
            
            
            
            this_comparison = comparison_indices{i_comparison};
            for i_condition = 1 : length(this_comparison)
                % find correct condition indicator
                condition_identifier = conditions_to_plot(this_comparison(i_condition), :);
                stance_foot_indicator = strcmp(condition_stance_foot_data, condition_identifier{1});
                perturbation_indicator = strcmp(condition_perturbation_data, condition_identifier{2});
                delay_indicator = strcmp(condition_delay_data, condition_identifier{3});
                index_indicator = strcmp(condition_index_data, condition_identifier{4});
                this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator;
                
                label_string = conditions_to_plot{comparison_indices{i_comparison}(i_condition), comparison_to_make};
                histogram ...
                  ( ...
                    points_to_plot(this_condition_indicator), ...
                    'edgecolor', colors_comparison(i_condition, :), ...
                    'facecolor', lightenColor(colors_comparison(i_condition, :), 0.5), ...
                    'DisplayName', label_string ...
                  )
            end

            % annotate
            if show_legend
                legend('toggle')
            end
            title_string = discrete_variable_info{i_variable, 2};
            for i_label = 1 : length(condition_labels);
                if i_label ~= comparison_to_make
                    title_string = [title_string ' - ' conditions_to_plot{comparison_indices{i_comparison}(1), i_label}];
                end
            end
            title(title_string); set(gca, 'Fontsize', 12)
        end
    end    
    
    
    
    % continuous variables
    for i_variable = 1 : size(continuous_variable_info, 1)
        trajectories_to_plot = continuous_variable_data{i_variable, 1};
        for i_comparison = 1 : length(comparison_indices);
            % find correct condition indicator for control
            this_comparison = comparison_indices{i_comparison};
            
            figure; axes; hold on;
            if ~isempty(conditions_control)
                % find correct condition indicator for control
                this_comparison = comparison_indices{i_comparison};
                representant_condition_index = this_comparison(1);
                if strcmp(conditions_to_plot(representant_condition_index, 1), 'STANCE_LEFT')
                    applicable_control_condition_indices = 1;
                elseif strcmp(conditions_to_plot(representant_condition_index, 1), 'STANCE_RIGHT')
                    applicable_control_condition_indices = 2;
                end
                stance_foot_indicator = strcmp(condition_stance_foot_data, conditions_control(applicable_control_condition_indices, 1));
                perturbation_indicator = strcmp(condition_perturbation_data, conditions_control(applicable_control_condition_indices, 2));
                delay_indicator = strcmp(condition_delay_data, conditions_control(applicable_control_condition_indices, 3));
                index_indicator = strcmp(condition_index_data, conditions_control(applicable_control_condition_indices, 4));
                experimental_indicator = strcmp(condition_experimental_data, conditions_control(applicable_control_condition_indices, 5));
                this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator;

                % plot control
                plot ...
                  ( ...
                    time_normalized, trajectories_to_plot(:, this_condition_indicator), ...
                    'HandleVisibility', 'off', ...
                    'color', lightenColor(color_control, 0.5) ...
                  );
                control_mean_plot = plot ...
                  ( ...
                    time_normalized, mean(trajectories_to_plot(:, this_condition_indicator), 2), ...
                    'DisplayName', 'control', ...
                    'linewidth', 5, ...
                    'color', color_control ...
                  );
            end
            
            % plot stimulus
            condition_mean_plots = zeros(1, length(this_comparison));
            for i_condition = 1 : length(this_comparison)
                % find correct condition indicator
                condition_identifier = conditions_to_plot(this_comparison(i_condition), :);
                stance_foot_indicator = strcmp(condition_stance_foot_data, condition_identifier{1});
                perturbation_indicator = strcmp(condition_perturbation_data, condition_identifier{2});
                delay_indicator = strcmp(condition_delay_data, condition_identifier{3});
                index_indicator = strcmp(condition_index_data, condition_identifier{4});
                experimental_indicator = strcmp(condition_experimental_data, condition_identifier{5});
                this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator;
                plot ...
                  ( ...
                    time_normalized, trajectories_to_plot(:, this_condition_indicator), ...
                    'HandleVisibility', 'off', ...
                    'color', lightenColor(colors_comparison(i_condition, :), 0.5) ...
                  );
                label_string = conditions_to_plot{comparison_indices{i_comparison}(i_condition), comparison_to_make};
                condition_mean_plots(i_condition) = plot ...
                  ( ...
                    time_normalized, mean(trajectories_to_plot(:, this_condition_indicator), 2), ...
                    'DisplayName', label_string, ...
                    'linewidth', 5, ...
                    'color', colors_comparison(i_condition, :) ...
                  );

            end

            % reorder
            if ~isempty(conditions_control)
                uistack(control_mean_plot, 'top');
            end
            for i_condition = 1 : length(this_comparison)
                uistack(condition_mean_plots(i_condition), 'top');
            end

            if dictate_axes
                set(gca, 'xlim', [time_normalized(1), time_normalized(end)]);
                set(gca, 'ylim', [continuous_variable_info{i_variable, 5}(1), continuous_variable_info{i_variable, 5}(2)]);
            end
            
            % annotate
            if show_legend
                legend('toggle')
            end
            title_string = continuous_variable_info{i_variable, 2};
            for i_label = 1 : length(condition_labels);
                if i_label ~= comparison_to_make
                    title_string = [title_string ' - ' conditions_to_plot{comparison_indices{i_comparison}(1), i_label}];
                end
            end
            title(title_string); set(gca, 'Fontsize', 12)
        end
    end
end



%% plot overview
if plot_overview

    for i_variable = 1 : size(discrete_variable_info, 1)
        points_to_plot = discrete_variable_data{i_variable, 1};
        for i_comparison = 1 : length(comparison_indices);
            % find correct condition indicator for control
            this_comparison = comparison_indices{i_comparison};
            representant_condition_index = this_comparison(1);
            
            % make condition labels for box plot
            condition_labels_for_boxplot = cell(length(points_to_plot), 1);
            condition_indicator_comparison = false(length(points_to_plot), length(this_comparison));
            
            if include_control
                % control
                this_comparison = comparison_indices{i_comparison};
                representant_condition_index = this_comparison(1);
                if strcmp(conditions_to_plot(representant_condition_index, 1), 'STANCE_LEFT')
                    applicable_control_condition_indices = 1;
                elseif strcmp(conditions_to_plot(representant_condition_index, 1), 'STANCE_RIGHT')
                    applicable_control_condition_indices = 2;
                end
                stance_foot_indicator = strcmp(condition_stance_foot_data, conditions_control(applicable_control_condition_indices, 1));
                perturbation_indicator = strcmp(condition_perturbation_data, conditions_control(applicable_control_condition_indices, 2));
                delay_indicator = strcmp(condition_delay_data, conditions_control(applicable_control_condition_indices, 3));
                index_indicator = strcmp(condition_index_data, conditions_control(applicable_control_condition_indices, 4));
                experimental_indicator = strcmp(condition_experimental_data, conditions_control(applicable_control_condition_indices, 5));
                condition_indicator_control = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator;
                % figure out condition label
                label_string = 'CONTROL';
                % place condition label into label cell array
                for i_point = 1 : length(condition_labels_for_boxplot)
                    if condition_indicator_control(i_point)
                        condition_labels_for_boxplot{i_point} = label_string;
                    end
                end
            end
            
            % stimulus
            for i_condition = 1 : length(this_comparison)
                % find correct condition indicator
                condition_identifier = conditions_to_plot(this_comparison(i_condition), :);
                stance_foot_indicator = strcmp(condition_stance_foot_data, condition_identifier{1});
                perturbation_indicator = strcmp(condition_perturbation_data, condition_identifier{2});
                delay_indicator = strcmp(condition_delay_data, condition_identifier{3});
                index_indicator = strcmp(condition_index_data, condition_identifier{4});
                experimental_indicator = strcmp(condition_experimental_data, condition_identifier{5});
                this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator;
                condition_indicator_comparison(:, i_condition) = this_condition_indicator;
                
                % figure out condition label
                label_string = conditions_to_plot{comparison_indices{i_comparison}(i_condition), comparison_to_make};
                % place condition label into label cell array
                for i_point = 1 : length(condition_labels_for_boxplot)
                    if this_condition_indicator(i_point)
                        condition_labels_for_boxplot{i_point} = label_string;
                    end
                end
            end            
            
            % prune data points that we don't want to look at in this plot
            if include_control
                condition_indicator_stimulus_or_control = [condition_indicator_control condition_indicator_comparison];
                indices_for_this_comparison = any(condition_indicator_stimulus_or_control, 2);
                conditions_pruned = condition_indicator_stimulus_or_control(indices_for_this_comparison, :);
                points_to_plot_pruned = points_to_plot(indices_for_this_comparison);
                condition_labels_for_boxplot_pruned = condition_labels_for_boxplot(indices_for_this_comparison);
                group_order = [conditions_to_plot(comparison_indices{i_comparison}, comparison_to_make); 'CONTROL'];
            else
                indices_for_this_comparison = any(condition_indicator_comparison, 2);
                conditions_pruned = condition_indicator_comparison(indices_for_this_comparison, :);
                points_to_plot_pruned = points_to_plot(indices_for_this_comparison);
                condition_labels_for_boxplot_pruned = condition_labels_for_boxplot(indices_for_this_comparison);
                group_order = conditions_to_plot(comparison_indices{i_comparison}, comparison_to_make);
            end
            
            % make plot
            figure; axes; hold on;
            box_plot_data = boxplot(points_to_plot_pruned, condition_labels_for_boxplot_pruned, 'grouporder', group_order, 'widths', 0.8);
            
            % color the boxes
            experimental_conditions = group_order;
            setBoxPlotColors(gca, box_plot_data, group_order, experimental_conditions, colors_comparison);

            % annotate
            title_string = discrete_variable_info{i_variable, 2};
            filename_string = discrete_variable_info{i_variable, 4};
            for i_label = 1 : length(condition_labels);
                if i_label ~= comparison_to_make
                    title_string = [title_string ' - ' strrep(conditions_to_analyze{comparison_indices{i_comparison}(1), i_label}, '_', ' ')];
                    filename_string = [filename_string '_' conditions_to_analyze{comparison_indices{i_comparison}(1), i_label}];
                end
            end
            title(title_string, 'interpreter', 'LaTeX'); set(gca, 'Fontsize', 12)            

            if dictate_axes
                xlimits = get(gca, 'xlim')
                set(gca, 'ylim', [discrete_variable_info{i_variable, 5}(1), discrete_variable_info{i_variable, 5}(2)]);
            end
            
            xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
            postext = text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), [discrete_variable_info{i_variable, 6} ' $\rightarrow$'] , 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right', 'interpreter', 'LaTeX');
            negtext = text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), ['$\leftarrow$ ' discrete_variable_info{i_variable, 7}], 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left', 'interpreter', 'LaTeX');
            
            % save
            if save_figures
                % figure out folders
                if ~exist('figures', 'dir')
                    mkdir('figures')
                end
                filename = ['figures/' filename_string '.eps'];
                saveas(gcf, filename, 'epsc2')
                
                % make labels invisible and save again
                zero_plot = plot(xlimits, [0 0], 'color', [0.7 0.7 0.7]);
                uistack(zero_plot, 'bottom')
                set(postext, 'visible', 'off');
                set(negtext, 'visible', 'off');
                set(get(gca, 'xaxis'), 'visible', 'off');
                set(get(gca, 'yaxis'), 'visible', 'off');
                set(get(gca, 'xlabel'), 'visible', 'off');
                set(get(gca, 'ylabel'), 'visible', 'off');
                set(get(gca, 'title'), 'visible', 'off');
                set(gca, 'xticklabel', '');
                set(gca, 'yticklabel', '');
                set(gca, 'position', [0 0 1 1]);
                filename = ['figures/' filename_string '_naked.eps'];
                saveas(gcf, filename, 'epsc2')
                
                close(gcf)
            end
        end
    end    
    

    for i_variable = 1 : size(continuous_variable_info, 1)
        trajectories_to_plot = continuous_variable_data{i_variable, 1};
        for i_comparison = 1 : length(comparison_indices);
            % find correct condition indicator for control
            this_comparison = comparison_indices{i_comparison};
            representant_condition_index = this_comparison(1);
            
            figure; axes; hold on;
            legend_handles = [];
            legend_data = {};
            if use_control
                stance_condition = conditions_to_plot(representant_condition_index, 1);
                applicable_control_condition_indices = find(strcmp(conditions_control(:, 1), stance_condition));

                stance_foot_indicator = strcmp(condition_stance_foot_data, conditions_control(applicable_control_condition_indices, 1));
                perturbation_indicator = strcmp(condition_perturbation_data, conditions_control(applicable_control_condition_indices, 2));
                delay_indicator = strcmp(condition_delay_data, conditions_control(applicable_control_condition_indices, 3));
                index_indicator = strcmp(condition_index_data, conditions_control(applicable_control_condition_indices, 4));
                experimental_indicator = strcmp(condition_experimental_data, conditions_control(applicable_control_condition_indices, 5));
                this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator;
                
                % plot control
                current_plots = shadedErrorBar ...
                  ( ...
                    time_normalized, ...
                    mean(trajectories_to_plot(:, this_condition_indicator), 2), ...
                    cinv(trajectories_to_plot(:, this_condition_indicator), 2), ...
                    { ...
                      'color', color_control, ...
                      'linewidth', 6 ...
                    }, ...
                    1 ...
                  );
                legend_handles = [legend_handles, current_plots.mainLine];
                legend_data = [legend_data, 'control'];
            end
            
            

            
            % plot stimulus
            this_comparison = comparison_indices{i_comparison};
            condition_mean_plots = zeros(1, length(this_comparison));
            for i_condition = 1 : length(this_comparison)
                % find correct condition indicator
                condition_identifier = conditions_to_plot(this_comparison(i_condition), :);
                stance_foot_indicator = strcmp(condition_stance_foot_data, condition_identifier{1});
                perturbation_indicator = strcmp(condition_perturbation_data, condition_identifier{2});
                delay_indicator = strcmp(condition_delay_data, condition_identifier{3});
                index_indicator = strcmp(condition_index_data, condition_identifier{4});
                experimental_indicator = strcmp(condition_experimental_data, condition_identifier{5});
                this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator;
                
                if strcmp(error_shades, 'cinv')
                    error_curves = cinv(trajectories_to_plot(:, this_condition_indicator), 2);
                elseif strcmp(error_shades, 'std')
                    error_curves = std(trajectories_to_plot(:, this_condition_indicator), 1, 2);
                end
                current_plots = shadedErrorBar ...
                  ( ...
                    time_normalized, ...
                    mean(trajectories_to_plot(:, this_condition_indicator), 2), ...
                    error_curves, ...
                    { ...
                      'color', colors_comparison(i_condition, :), ...
                      'linewidth', 6 ...
                    }, ...
                    1 ...
                  );
                legend_handles = [legend_handles, current_plots.mainLine];
                label_string = strrep(conditions_to_plot{comparison_indices{i_comparison}(i_condition), comparison_to_make}, '_', ' ');
                legend_data = [legend_data, label_string];
            end

            % annotate
            if show_legend
                this_legend = legend(legend_handles, legend_data);
            end
            title_string = continuous_variable_info{i_variable, 2};
            filename_string = continuous_variable_info{i_variable, 4};
            for i_label = 1 : length(condition_labels);
                if i_label ~= comparison_to_make
                    title_string = [title_string ' - ' strrep(conditions_to_plot{comparison_indices{i_comparison}(1), i_label}, '_', ' ')];
                    filename_string = [filename_string '_' conditionStringToFilename(conditions_to_plot{comparison_indices{i_comparison}(1), i_label})];
                end
            end
            title(title_string, 'interpreter', 'LaTeX'); set(gca, 'Fontsize', 12)
            
            set(gca, 'xlim', [time_normalized(1), time_normalized(end)]);
            xlabel('normalized time (s)');
            ylabel(continuous_variable_info{i_variable, 3});
            
            if dictate_axes
                set(gca, 'xlim', [time_normalized(1), time_normalized(end)]);
                set(gca, 'ylim', [continuous_variable_info{i_variable, 5}(1), continuous_variable_info{i_variable, 5}(2)]);
            end
            
            xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
            postext = text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), [continuous_variable_info{i_variable, 6} ' $\rightarrow$'] , 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right', 'interpreter', 'LaTeX');
            negtext = text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), ['$\leftarrow$ ' continuous_variable_info{i_variable, 7}], 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left', 'interpreter', 'LaTeX');
            
            if mark_pushoff
                pushoff_patch_x = [0 pushoff_time_mean pushoff_time_mean 0];
                pushoff_patch_y = [ylimits(1) ylimits(1) ylimits(2) ylimits(2)];
                double_stance_patch_vertices = [pushoff_patch_x'; pushoff_patch_y'];
                patch_handle = ...
                    patch ...
                      ( ...
                        pushoff_patch_x, ...
                        pushoff_patch_y, ...
                        double_stance_color, ...
                        'EdgeColor', 'none', ...
                        'FaceAlpha', double_stance_alpha ...
                      ); 
                uistack(patch_handle, 'bottom')
            end
            
            % save
            if save_figures
                % figure out folders
                if ~exist('figures', 'dir')
                    mkdir('figures')
                end
                filename = ['figures/' filename_string '.eps'];
                saveas(gcf, filename, 'epsc2')
                
                % make labels invisible and save again
                set(postext, 'visible', 'off');
                set(negtext, 'visible', 'off');
                set(get(gca, 'xaxis'), 'visible', 'off');
                set(get(gca, 'yaxis'), 'visible', 'off');
                set(get(gca, 'xlabel'), 'visible', 'off');
                set(get(gca, 'ylabel'), 'visible', 'off');
                set(get(gca, 'title'), 'visible', 'off');
                set(gca, 'xticklabel', '');
                set(gca, 'yticklabel', '');
                set(gca, 'position', [0 0 1 1]);
                filename = ['figures/' filename_string '_naked.eps'];
                saveas(gcf, filename, 'epsc2');
                
                close(gcf)
            end
        end
    end    
    
    
    
    
    
end

%% plot episodes
if plot_episodes

    for i_variable = 1 : size(continuous_variable_info, 1)
        trajectories_to_plot = continuous_variable_data{i_variable, 1};
        for i_episode = 1 : length(episode_indices);
            % extract data for steps
            this_episode = episode_indices{i_episode};
            figure; axes; hold on;
            legend_handles = [];
            legend_data = {};
            
            for i_comparison = 1 : length(this_episode)
                
                
                this_comparison = comparison_indices{this_episode(i_comparison)};
                condition_mean_plots = zeros(1, length(this_comparison));
                for i_condition = 1 : length(this_comparison)
                    % find correct condition indicator
                    condition_identifier = conditions_to_plot(this_comparison(i_condition), :);
                    stance_foot_indicator = strcmp(condition_stance_foot_data, condition_identifier{1});
                    perturbation_indicator = strcmp(condition_perturbation_data, condition_identifier{2});
                    delay_indicator = strcmp(condition_delay_data, condition_identifier{3});
                    index_indicator = strcmp(condition_index_data, condition_identifier{4});
                    experimental_indicator = strcmp(condition_experimental_data, condition_identifier{5});
                    this_condition_indicator = stance_foot_indicator & perturbation_indicator & delay_indicator & index_indicator & experimental_indicator;
                    
                    if strcmp(condition_identifier{4}, 'ONE')
                        time_episode = time_normalized;
                    elseif strcmp(condition_identifier{4}, 'TWO')
                        time_episode = time_normalized + time_normalized(end);
                    elseif strcmp(condition_identifier{4}, 'THREE')
                        time_episode = time_normalized + time_normalized(end) * 2;
                    elseif strcmp(condition_identifier{4}, 'FOUR')
                        time_episode = time_normalized + time_normalized(end) * 3;
                    end
                    
                    if strcmp(error_shades, 'cinv')
                        error_curves = cinv(trajectories_to_plot(:, this_condition_indicator), 2);
                    elseif strcmp(error_shades, 'std')
                        error_curves = std(trajectories_to_plot(:, this_condition_indicator), 1, 2);
                    end
%                     current_plots = shadedErrorBar ...
%                       ( ...
%                         time_episode, ...
%                         mean(trajectories_to_plot(:, this_condition_indicator), 2), ...
%                         error_curves, ...
%                         { ...
%                           'color', colors_comparison(i_condition, :), ...
%                           'linewidth', 6 ...
%                         }, ...
%                         1 ...
%                       );
                    current_plots = plot ...
                      ( ...
                        time_episode, ...
                        mean(trajectories_to_plot(:, this_condition_indicator), 2), ...
                        'color', colors_comparison(i_condition, :), ...
                        'linewidth', 6 ...
                      );
                end
            end            

            % plot zero
            zero_line = plot([0 time_normalized(end) * 4], [0 0], 'color', [1 1 1]*0.7)
            uistack(zero_line, 'bottom')

            
            % annotate % TODO: fix the title
            if show_legend
                this_legend = legend(legend_handles, legend_data);
            end
            title_string = continuous_variable_info{i_variable, 2};
            filename_string = continuous_variable_info{i_variable, 4};
            for i_label = 1 : length(condition_labels);
                if i_label ~= comparison_to_make
                    title_string = [title_string ' - ' strrep(conditions_to_plot{comparison_indices{this_episode(i_comparison)}(1), i_label}, '_', ' ')];
                    if i_label == condition_column_index
                        filename_string = [filename_string '_all'];
                    else
                        filename_string = [filename_string '_' conditionStringToFilename(conditions_to_plot{comparison_indices{this_episode(i_comparison)}(1), i_label})];
                    end
                end
            end
            title(title_string, 'interpreter', 'LaTeX'); set(gca, 'Fontsize', 12)
            
            set(gca, 'xlim', [0 time_normalized(end) * 4]);
            xlabel('normalized time');
            ylabel(continuous_variable_info{i_variable, 3});
            
            if dictate_axes
%                 set(gca, 'xlim', [0 time_normalized(end) * 4]);
                set(gca, 'xlim', [0 1.2]); % for grant figure
                set(gca, 'ylim', [continuous_variable_info{i_variable, 5}(1), continuous_variable_info{i_variable, 5}(2)]);
            end
            
            xlimits = get(gca, 'xlim'); ylimits = get(gca, 'ylim');
            postext = text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(2), [continuous_variable_info{i_variable, 6} ' $\rightarrow$'] , 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'right', 'interpreter', 'LaTeX');
            negtext = text(xlimits(1) - (xlimits(2)-xlimits(1))*0.12, ylimits(1), ['$\leftarrow$ ' continuous_variable_info{i_variable, 7}], 'rotation', 90, 'Fontsize', 24, 'horizontalalignment', 'left', 'interpreter', 'LaTeX');
            
            if mark_pushoff
                for i_step = 1 : 4
%                     pushoff_time_percentage = pushoff_time_mean / time_normalized(end) * 100;
%                     pushoff_patch_x = [0 pushoff_time_percentage pushoff_time_percentage 0] + (i_step-1)*100;
%                     pushoff_patch_y = [ylimits(1) ylimits(1) ylimits(2) ylimits(2)];
                    heelstrike_time = time_normalized(end) * (i_step-1);
                    pushoff_time = heelstrike_time + pushoff_time_mean;
                    pushoff_patch_x = [heelstrike_time pushoff_time pushoff_time heelstrike_time];
                    pushoff_patch_y = [ylimits(1) ylimits(1) ylimits(2) ylimits(2)];

                    double_stance_patch_vertices = [pushoff_patch_x'; pushoff_patch_y'];
                    patch_handle = ...
                        patch ...
                          ( ...
                            pushoff_patch_x, ...
                            pushoff_patch_y, ...
                            double_stance_color, ...
                            'EdgeColor', 'none', ...
                            'FaceAlpha', double_stance_alpha ...
                          ); 
                    uistack(patch_handle, 'bottom')
                end
            end
            
            % save
            if save_figures
                % figure out folders
                if ~exist('figures', 'dir')
                    mkdir('figures')
                end
                filename = ['figures/' filename_string '.eps'];
                saveas(gcf, filename, 'epsc2')
                
                % make labels invisible and save again
                set(postext, 'visible', 'off');
                set(negtext, 'visible', 'off');
                set(get(gca, 'xaxis'), 'visible', 'off');
                set(get(gca, 'yaxis'), 'visible', 'off');
                set(get(gca, 'xlabel'), 'visible', 'off');
                set(get(gca, 'ylabel'), 'visible', 'off');
                set(get(gca, 'title'), 'visible', 'off');
                set(gca, 'xticklabel', '');
                set(gca, 'yticklabel', '');
                set(gca, 'position', [0 0 1 1]);
                filename = ['figures/' filename_string '_naked.eps'];
                saveas(gcf, filename, 'epsc2');
                
                close(gcf)
            end            
            
        end
    end    
    
    
end













