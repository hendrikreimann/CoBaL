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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% this function applies several basic processing steps to experimental data, e.g. filtering

% input: 
% Experimental data files generated by importAscii.m, in the subfolder "raw"
%
% output: 
% multiple files with processed data for each trial, in the subfolder "processed"


function preprocessForceplateData(varargin)
    % parse arguments
    [types_to_analyze, trials_to_analyze, types_to_exclude, trials_to_exclude] = parseTrialArguments(varargin{:});
    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, 'type', 'all')
    parse(parser, varargin{:})

    % add excluded trials back in, because while we don't want to analyze them, we still want to pre-process them
    types_to_analyze = [types_to_analyze; types_to_exclude];
    trials_to_analyze = [trials_to_analyze; trials_to_exclude];
    
    % load settings
    study_settings = loadSettingsFromFile('study');
    subject_settings = loadSettingsFromFile('subject');
    
    % extract forceplate settings
    forceplate_table = study_settings.getTable('forceplate_table');
    number_of_forceplates = size(forceplate_table, 1);
    
    data_dir = dir(['raw' filesep '*_forceplateTrajectoriesRaw.mat']);
    clear file_name_list;
    [file_name_list{1:length(data_dir)}] = deal(data_dir.name);
    number_of_files = length(file_name_list);
    for i_trial = 1 : number_of_files
        raw_forceplate_file_name = file_name_list{i_trial};
        [date, subject_id, trial_type, trial_number] = getFileParameters(raw_forceplate_file_name);

        % does the user want to process this file?
        if any(strcmp(trial_type, types_to_analyze))
            % condition is set to be processed, now check trial number
            trial_number_list_this_condition = trials_to_analyze{strcmp(trial_type, types_to_analyze)};
            if ismember(trial_number, trial_number_list_this_condition)
                % load raw data
                loaded_data = load(['raw' filesep raw_forceplate_file_name]);
                
                % figure out variable
                if isfield(loaded_data, 'forceplate_trajectories_raw')
                    raw_forceplate_trajectories = loaded_data.forceplate_trajectories_raw;
                elseif isfield(loaded_data, 'forceplate_raw_trajectories')
                    raw_forceplate_trajectories = loaded_data.forceplate_raw_trajectories;
                end
                
                % remove data points that are either unreasonably large or 0
                bad_data_points = false(size(raw_forceplate_trajectories, 1), 1);
                if subject_settings.get('remove_forceplate_zero_data_points', 1)
                    zero_data_points = any(raw_forceplate_trajectories==0, 2);
                    bad_data_points = bad_data_points | zero_data_points;
                    
                    % store problematic data information
                    saveProblemInformation(date, subject_id, trial_type, trial_number, 'zero in forceplate data', loaded_data.time_forceplate, zero_data_points)
                end
                if subject_settings.get('remove_forceplate_large_data_points', 1)
                    large_data_points = any(abs(raw_forceplate_trajectories)>1e13, 2);
                    bad_data_points = bad_data_points | large_data_points;
                    
                    % store problematic data information
                    saveProblemInformation(date, subject_id, trial_type, trial_number, 'large value in forceplate data', loaded_data.time_forceplate, large_data_points)
                end
                raw_forceplate_trajectories_without_bad = raw_forceplate_trajectories;
                raw_forceplate_trajectories_without_bad(bad_data_points, :) = NaN;
                
                % filter
                filter_order_low = study_settings.get('force_plate_filter_order');
                cutoff_frequency_low = study_settings.get('force_plate_filter_cutoff');
                if ~isempty(filter_order_low) && ~isempty(cutoff_frequency_low)
                    [b_lowpass, a_lowpass] = butter(filter_order_low, cutoff_frequency_low/(loaded_data.sampling_rate_forceplate/2), 'low');
                    forceplate_trajectories_groomed = nanfiltfilt(b_lowpass, a_lowpass, raw_forceplate_trajectories_without_bad);
                else
                    forceplate_trajectories_groomed = raw_forceplate_trajectories;
                end
                
                % calculate CoP and total wrench
                fz_threshold = subject_settings.get('forceplate_load_threshold');
                wrenches = {};
                forceplate_trajectories = [];
                forceplate_labels = {};
                forceplate_directions = {};
                for i_forceplate = 1 : number_of_forceplates
                    % extract data for this forceplate
                    this_plate_label = forceplate_table.label{i_forceplate};
                    this_plate_fx = forceplate_trajectories_groomed(:, strcmp(loaded_data.forceplate_labels, ['fx_' this_plate_label]));
                    this_plate_fy = forceplate_trajectories_groomed(:, strcmp(loaded_data.forceplate_labels, ['fy_' this_plate_label]));
                    this_plate_fz = forceplate_trajectories_groomed(:, strcmp(loaded_data.forceplate_labels, ['fz_' this_plate_label]));
                    this_plate_mx = forceplate_trajectories_groomed(:, strcmp(loaded_data.forceplate_labels, ['mx_' this_plate_label]));
                    this_plate_my = forceplate_trajectories_groomed(:, strcmp(loaded_data.forceplate_labels, ['my_' this_plate_label]));
                    this_plate_mz = forceplate_trajectories_groomed(:, strcmp(loaded_data.forceplate_labels, ['mz_' this_plate_label]));
                    this_plate_wrench = [this_plate_fx this_plate_fy this_plate_fz this_plate_mx this_plate_my this_plate_mz];
                    wrenches = [wrenches; this_plate_wrench]; %#ok<AGROW>
                    
                    % calculate CoP
                    this_place_copx = - this_plate_my ./ this_plate_fz;
                    this_place_copy =   this_plate_mx ./ this_plate_fz;
                    this_plate_low_load_indicator = (this_plate_fz < fz_threshold);
                    this_place_copx(this_plate_low_load_indicator, :) = 0;
                    this_place_copy(this_plate_low_load_indicator, :) = 0;
                    
                    % define labels
                    this_plate_labels = ...
                      { ...
                        ['fx_' this_plate_label], ...
                        ['fy_' this_plate_label], ...
                        ['fz_' this_plate_label], ...
                        ['mx_' this_plate_label], ...
                        ['my_' this_plate_label], ...
                        ['mz_' this_plate_label], ...
                        ['copx_' this_plate_label], ...
                        ['copy_' this_plate_label], ...
                      };
                    
                    % define directions
                    this_wrench_directions = ...
                      { ...
                        study_settings.get('direction_x_pos'), ...
                        study_settings.get('direction_y_pos'), ...
                        study_settings.get('direction_z_pos'), ....
                        study_settings.get('direction_x_pos'), ...
                        study_settings.get('direction_y_pos'), ...
                        study_settings.get('direction_z_pos'); ....
                        study_settings.get('direction_x_neg'), ...
                        study_settings.get('direction_y_neg'), ...
                        study_settings.get('direction_z_neg'), ....
                        study_settings.get('direction_x_neg'), ...
                        study_settings.get('direction_y_neg'), ...
                        study_settings.get('direction_z_neg'); ....
                      };
                    this_cop_directions = ...
                      { ...
                        study_settings.get('direction_x_pos'), ...
                        study_settings.get('direction_y_pos'); ...
                        study_settings.get('direction_x_neg'), ...
                        study_settings.get('direction_y_neg'), ...
                      };
                    
                    
                    % store
                    forceplate_trajectories = [forceplate_trajectories this_plate_wrench this_place_copx this_place_copy]; %#ok<AGROW>
                    forceplate_labels = [forceplate_labels this_plate_labels]; %#ok<AGROW>
                    forceplate_directions = [forceplate_directions this_wrench_directions this_cop_directions]; %#ok<AGROW>
                end
                
                % calculate total wrench and CoP
                total_forceplate_wrench = wrenches{1};
                for i_forceplate = 2 : number_of_forceplates
                    total_forceplate_wrench = total_forceplate_wrench + wrenches{i_forceplate};
                end
                copx_total = - total_forceplate_wrench(:, 5) ./ total_forceplate_wrench(:, 3);
                copy_total = total_forceplate_wrench(:, 4) ./ total_forceplate_wrench(:, 3);
                low_load_indicator = (total_forceplate_wrench(:, 3) < fz_threshold);
                copx_total(low_load_indicator, :) = 0;
                copy_total(low_load_indicator, :) = 0;
                
                % define labels and directions for total wrench
               total_labels = ...
                  { ...
                    'fx', ...
                    'fy', ...
                    'fz', ...
                    'mx', ...
                    'my', ...
                    'mz', ...
                    'copx', ...
                    'copy', ...
                  };
                total_wrench_directions = ...
                  { ...
                    study_settings.get('direction_x_pos'), ...
                    study_settings.get('direction_y_pos'), ...
                    study_settings.get('direction_z_pos'), ....
                    study_settings.get('direction_x_pos'), ...
                    study_settings.get('direction_y_pos'), ...
                    study_settings.get('direction_z_pos'); ....
                    study_settings.get('direction_x_neg'), ...
                    study_settings.get('direction_y_neg'), ...
                    study_settings.get('direction_z_neg'), ....
                    study_settings.get('direction_x_neg'), ...
                    study_settings.get('direction_y_neg'), ...
                    study_settings.get('direction_z_neg'); ....
                  };
                total_cop_directions = ...
                  { ...
                    study_settings.get('direction_x_pos'), ...
                    study_settings.get('direction_y_pos'); ...
                    study_settings.get('direction_x_neg'), ...
                    study_settings.get('direction_y_neg'), ...
                  };

                % store total wrench and CoP
                forceplate_trajectories = [total_forceplate_wrench copx_total copy_total forceplate_trajectories]; %#ok<AGROW>
                forceplate_labels = [total_labels forceplate_labels]; %#ok<AGROW>
                forceplate_directions = [total_wrench_directions total_cop_directions forceplate_directions]; %#ok<AGROW>

%                 % extract
%                 fxl_trajectory = forceplate_trajectories(:, 1);
%                 fyl_trajectory = forceplate_trajectories(:, 2);
%                 fzl_trajectory = forceplate_trajectories(:, 3);
%                 mxl_trajectory = forceplate_trajectories(:, 4);
%                 myl_trajectory = forceplate_trajectories(:, 5);
%                 mzl_trajectory = forceplate_trajectories(:, 6);
%                 fxr_trajectory = forceplate_trajectories(:, 7);
%                 fyr_trajectory = forceplate_trajectories(:, 8);
%                 fzr_trajectory = forceplate_trajectories(:, 9);
%                 mxr_trajectory = forceplate_trajectories(:, 10);
%                 myr_trajectory = forceplate_trajectories(:, 11);
%                 mzr_trajectory = forceplate_trajectories(:, 12);
% 

%                 if strcmp(loaded_data.data_source, 'nexus')
%                     % transform forceplate data to CoBaL world frame A_cw
%                     left_forceplate_wrench_Acl = [fxl_trajectory fyl_trajectory fzl_trajectory mxl_trajectory myl_trajectory mzl_trajectory];
% %                     left_forceplate_cop_Acl = [copxl_trajectory copyl_trajectory zeros(size(copxl_trajectory))];
%                     right_forceplate_wrench_Acr = [fxr_trajectory fyr_trajectory fzr_trajectory mxr_trajectory myr_trajectory mzr_trajectory];
% %                     right_forceplate_cop_Acr = [copxr_trajectory copyr_trajectory zeros(size(copxr_trajectory))];
% 
%                     % define forceplate rotation and translation
%                     Acr_to_world_rotation = [-1 0 0; 0 1 0; 0 0 -1];
%                     Acr_to_world_translation = [0.5588; 0; 0]; % origin of Acw in Acr frame
%                     Acr_to_world_trafo = [Acr_to_world_rotation Acr_to_world_translation; 0 0 0 1];
%                     Acl_to_world_rotation = [-1 0 0; 0 1 0; 0 0 -1];
%                     Acl_to_world_translation = [-0.5588; 0; 0]; % origin of Acw in Acl frame
%                     Acl_to_world_trafo = [Acl_to_world_rotation Acl_to_world_translation; 0 0 0 1];
%                     Acr_to_world_adjoint = rigidToAdjointTransformation(Acr_to_world_trafo);
%                     Acl_to_world_adjoint = rigidToAdjointTransformation(Acl_to_world_trafo);
% 
%                     % transform
%                     left_forceplate_wrench_world = (Acl_to_world_adjoint' * left_forceplate_wrench_Acl')';
% %                     left_forceplate_cop_world = (eye(2, 4) * Acl_to_world_trafo * [left_forceplate_cop_Acl ones(size(left_forceplate_cop_Acl, 1), 1)]')';
%                     right_forceplate_wrench_world = (Acr_to_world_adjoint' * right_forceplate_wrench_Acr')';
% %                     right_forceplate_cop_world = (eye(2, 4) * Acr_to_world_trafo * [right_forceplate_cop_Acr ones(size(right_forceplate_cop_Acr, 1), 1)]')';
%                 elseif strcmp(loaded_data.data_source, 'qtm')
%                     % transform forceplate data to CoBaL world frame A_cw
%                     left_forceplate_wrench_Acl = [fxl_trajectory fyl_trajectory fzl_trajectory mxl_trajectory myl_trajectory mzl_trajectory];
%                     right_forceplate_wrench_Acr = [fxr_trajectory fyr_trajectory fzr_trajectory mxr_trajectory myr_trajectory mzr_trajectory];
%                     
%                     % define forceplate rotation and translation
%                     lab_orientation = study_settings.get('lab_orientation', 1);
%                     if strcmp(lab_orientation, 'normal')
%                         Acr_to_world_rotation = [-1 0 0; 0 1 0; 0 0 -1];
%                         Acl_to_world_rotation = [-1 0 0; 0 1 0; 0 0 -1];
%                     elseif strcmp(lab_orientation, 'inverse')
%                         Acl_to_world_rotation = [1 0 0; 0 -1 0; 0 0 -1];
%                         Acr_to_world_rotation = [1 0 0; 0 -1 0; 0 0 -1];
%                     else
%                         error('Setting "lab_orientation" in studySettings.txt must either be "normal" or "inverse"')                        
%                     end
%                     
%                     
%                     Acr_to_world_translation = loaded_data.forceplate_location_right';
%                     Acr_to_world_trafo = [Acr_to_world_rotation Acr_to_world_translation; 0 0 0 1];
%                     
% %                     Acl_to_world_rotation = [-.999991 -.001155 .004082; -.001156 .999999 -.000193; -.004082 -.000197 -.999992];
%                     Acl_to_world_translation = loaded_data.forceplate_location_left';
%                     Acl_to_world_trafo = [Acl_to_world_rotation Acl_to_world_translation; 0 0 0 1];
%                     
%                     Acr_to_world_adjoint = rigidToAdjointTransformation(Acr_to_world_trafo);
%                     Acl_to_world_adjoint = rigidToAdjointTransformation(Acl_to_world_trafo);
% 
%                     % transform
%                     left_forceplate_wrench_world = (Acl_to_world_adjoint' * left_forceplate_wrench_Acl')';
%                     right_forceplate_wrench_world = (Acr_to_world_adjoint' * right_forceplate_wrench_Acr')';
% 
%                 elseif strcmp(loaded_data.data_source, 'neurocom')
%                     % transform forceplate data to VEPO world frame A_vw
%                     left_forceplate_wrench_Anl = [fxl_trajectory fyl_trajectory fzl_trajectory mxl_trajectory myl_trajectory mzl_trajectory];
% %                     left_forceplate_cop_Anl = [copxl_trajectory copyl_trajectory zeros(size(copxl_trajectory))];
%                     right_forceplate_wrench_Anr = [fxr_trajectory fyr_trajectory fzr_trajectory mxr_trajectory myr_trajectory mzr_trajectory];
% %                     right_forceplate_cop_Anr = [copxr_trajectory copyr_trajectory zeros(size(copxr_trajectory))];
% 
%                     % define forceplate rotation and translation
%                     inchToMeter = 0.0254;
% 
%                     Anl_to_world_rotation = [-1 0 0; 0 1 0; 0 0 -1];
%                     Anl_to_world_translation = [5; -10; 0] * inchToMeter;
% %                     Anl_to_world_translation = [5; 10; 0] * 0;
%                     Anl_to_world_trafo = [Anl_to_world_rotation Anl_to_world_translation; 0 0 0 1];
%                     Anl_to_world_adjoint = rigidToAdjointTransformation(Anl_to_world_trafo);
% 
%                     Anr_to_world_rotation = [-1 0 0; 0 1 0; 0 0 -1];
%                     Anr_to_world_translation = [-5; -10; 0] * inchToMeter;
%                     Anr_to_world_trafo = [Anr_to_world_rotation Anr_to_world_translation; 0 0 0 1];
%                     Anr_to_world_adjoint = rigidToAdjointTransformation(Anr_to_world_trafo);
% 
%                     % transform
%                     left_forceplate_wrench_world = (Anl_to_world_adjoint' * left_forceplate_wrench_Anl')';
% %                     left_forceplate_cop_world = (eye(2, 4) * Anl_to_world_trafo * [left_forceplate_cop_Anl ones(size(left_forceplate_cop_Anl, 1), 1)]')';
%                     right_forceplate_wrench_world = (Anr_to_world_adjoint' * right_forceplate_wrench_Anr')';
% %                     right_forceplate_cop_world = (eye(2, 4) * Anr_to_world_trafo * [right_forceplate_cop_Anr ones(size(right_forceplate_cop_Anr, 1), 1)]')';
% 
%                 else
%                     error(['data source "' loaded_data.data_source '" not recognized'])
%                 end

                % calculate wrench for complete plate
%                 total_forceplate_wrench_world = left_forceplate_wrench_world + right_forceplate_wrench_world;

                % calculate CoP
%                 copxl_raw = - myl_trajectory ./ fzl_trajectory;
%                 copyl_raw = mxl_trajectory ./ fzl_trajectory;
%                 copxr_raw = - myr_trajectory ./ fzr_trajectory;
%                 copyr_raw = mxr_trajectory ./ fzr_trajectory;


%                 copxl_world = - left_forceplate_wrench_world(:, 5) ./ left_forceplate_wrench_world(:, 3);
%                 copyl_world = left_forceplate_wrench_world(:, 4) ./ left_forceplate_wrench_world(:, 3);
%                 copxr_world = - right_forceplate_wrench_world(:, 5) ./ right_forceplate_wrench_world(:, 3);
%                 copyr_world = right_forceplate_wrench_world(:, 4) ./ right_forceplate_wrench_world(:, 3);
%                 copx_world = - total_forceplate_wrench_world(:, 5) ./ total_forceplate_wrench_world(:, 3);
%                 copy_world = total_forceplate_wrench_world(:, 4) ./ total_forceplate_wrench_world(:, 3);
% %                 copxr_world = - left_forceplate_wrench_world(:, 5) ./ left_forceplate_wrench_world(:, 3);
% %                 copyr_world = left_forceplate_wrench_world(:, 4) ./ left_forceplate_wrench_world(:, 3);
% %                 copx_world = - left_forceplate_wrench_world(:, 5) ./ left_forceplate_wrench_world(:, 3);
% %                 copy_world = left_forceplate_wrench_world(:, 4) ./ left_forceplate_wrench_world(:, 3);
%                 % re-zero CoP for low loads
%                 left_forceplate_low_load_indicator = (fzl_trajectory < fz_threshold);
% %                 copxl_raw(left_forceplate_low_load_indicator, :) = 0;
% %                 copyl_raw(left_forceplate_low_load_indicator, :) = 0;
%                 copxl_world(left_forceplate_low_load_indicator, :) = 0;
%                 copyl_world(left_forceplate_low_load_indicator, :) = 0;
% 
%                 right_forceplate_low_load_indicator = (fzr_trajectory < fz_threshold);
% %                 copxr_raw(right_forceplate_low_load_indicator, :) = 0;
% %                 copyr_raw(right_forceplate_low_load_indicator, :) = 0;
%                 copxr_world(right_forceplate_low_load_indicator, :) = 0;
%                 copyr_world(right_forceplate_low_load_indicator, :) = 0;
% 
%                 total_forceplate_low_load_indicator = (left_forceplate_low_load_indicator & right_forceplate_low_load_indicator);
%                 copx_world(total_forceplate_low_load_indicator, :) = 0;
%                 copy_world(total_forceplate_low_load_indicator, :) = 0;
% 
%                 left_forceplate_cop_world = [copxl_world copyl_world];
%                 right_forceplate_cop_world = [copxr_world copyr_world];
%                 total_forceplate_cop_world = [copx_world copy_world];
% 
%                 % directions for wrench
%                 wrench_directions = cell(2, 6);
%                 [wrench_directions{1, [1 4]}] = deal('right');
%                 [wrench_directions{2, [1 4]}] = deal('left');
%                 [wrench_directions{1, [2 5]}] = deal('forward');
%                 [wrench_directions{2, [2 5]}] = deal('backward');
%                 [wrench_directions{1, [3 6]}] = deal('up');
%                 [wrench_directions{2, [3 6]}] = deal('down');
% 
%                 % directions for cop
%                 cop_directions = cell(2, 2);
%                 [cop_directions{1, 1}] = deal('right');
%                 [cop_directions{2, 1}] = deal('left');
%                 [cop_directions{1, 2}] = deal('forward');
%                 [cop_directions{2, 2}] = deal('backward');
% 
%                 % labels
%                 left_foot_wrench_labels = {'fxl', 'fyl', 'fzl', 'mxl', 'myl', 'mzl'};
%                 right_foot_wrench_labels = {'fxr', 'fyr', 'fzr', 'mxr', 'myr', 'mzr'};
%                 total_wrench_labels = {'fx', 'fy', 'fz', 'mx', 'my', 'mz'};
%                 left_foot_cop_labels = {'copxl', 'copyl'};
%                 right_foot_cop_labels = {'copxr', 'copyr'};
%                 total_cop_labels = {'copx', 'copy'};
                
                % extract
                time = loaded_data.time_forceplate;
                sampling_rate = loaded_data.sampling_rate_forceplate;
                
                % consolidate into single variable
%                 forceplate_trajectories = ...
%                   [ ...
%                     total_forceplate_wrench_world, ...
%                     total_forceplate_cop_world, ...
%                     left_forceplate_wrench_world, ...
%                     left_forceplate_cop_world, ...
%                     right_forceplate_wrench_world, ...
%                     right_forceplate_cop_world ...
%                   ];
%                 forceplate_labels = ...
%                   [ ...
%                     total_wrench_labels, ...
%                     total_cop_labels, ...
%                     left_foot_wrench_labels, ...
%                     left_foot_cop_labels, ...
%                     right_foot_wrench_labels, ...
%                     right_foot_cop_labels ...
%                   ];
%                 forceplate_directions = ...
%                   [ ...
%                     wrench_directions, ...
%                     cop_directions, ...
%                     wrench_directions, ...
%                     cop_directions, ...
%                     wrench_directions, ...
%                     cop_directions ...
%                   ];

                % save
                save_folder = 'processed';
%                 save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'forceplateTrajectoriesOld.mat');
%                 save ...
%                   ( ...
%                     [save_folder filesep save_file_name], ...
%                     'left_forceplate_wrench_world', ...
%                     'left_forceplate_cop_world', ...
%                     'right_forceplate_wrench_world', ...
%                     'right_forceplate_cop_world', ...
%                     'left_foot_wrench_world', ...
%                     'left_foot_cop_world', ...
%                     'right_foot_wrench_world', ...
%                     'right_foot_cop_world', ...
%                     'total_forceplate_wrench_world', ...
%                     'total_forceplate_cop_world', ...
%                     'fxl_trajectory', ...
%                     'fyl_trajectory', ...
%                     'fzl_trajectory', ...
%                     'mxl_trajectory', ...
%                     'myl_trajectory', ...
%                     'mzl_trajectory', ...
%                     'fxr_trajectory', ...
%                     'fyr_trajectory', ...
%                     'fzr_trajectory', ...
%                     'mxr_trajectory', ...
%                     'myr_trajectory', ...
%                     'mzr_trajectory', ...
%                     'time_forceplate', ...
%                     'sampling_rate_forceplate', ...
%                     'wrench_directions', ...
%                     'cop_directions', ...
%                     'left_foot_wrench_labels', ...
%                     'right_foot_wrench_labels', ...
%                     'total_wrench_labels', ...
%                     'left_foot_cop_labels', ...
%                     'right_foot_cop_labels', ...
%                     'total_cop_labels', ...
%                     'cop_trajectories', ...
%                     'cop_labels', ...
%                     'cop_all_directions' ...
%                   );
% 
%                 addAvailableData('left_foot_wrench_world', 'time_forceplate', 'sampling_rate_forceplate', '_left_foot_wrench_labels', '_wrench_directions', save_folder, save_file_name);
%                 addAvailableData('left_foot_cop_world', 'time_forceplate', 'sampling_rate_forceplate', '_left_foot_cop_labels', '_cop_directions', save_folder, save_file_name);
%                 addAvailableData('right_foot_wrench_world', 'time_forceplate', 'sampling_rate_forceplate', '_right_foot_wrench_labels', '_wrench_directions', save_folder, save_file_name);
%                 addAvailableData('right_foot_cop_world', 'time_forceplate', 'sampling_rate_forceplate', '_right_foot_cop_labels', '_cop_directions', save_folder, save_file_name);
% 
%                 addAvailableData('total_forceplate_wrench_world', 'time_forceplate', 'sampling_rate_forceplate', '_total_wrench_labels', '_wrench_directions', save_folder, save_file_name);
%                 addAvailableData('total_forceplate_cop_world', 'time_forceplate', 'sampling_rate_forceplate', '_total_cop_labels', '_cop_directions', save_folder, save_file_name);
% 
%                 addAvailableData('fxl_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'fxl_trajectory', {'right', 'left'}, save_folder, save_file_name);
%                 addAvailableData('fyl_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'fyl_trajectory', {'forward', 'backward'}, save_folder, save_file_name);
%                 addAvailableData('fzl_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'fzl_trajectory', {'up', 'down'}, save_folder, save_file_name);
%                 addAvailableData('mxl_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'mxl_trajectory', {'right', 'left'}, save_folder, save_file_name);
%                 addAvailableData('myl_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'myl_trajectory', {'forward', 'backward'}, save_folder, save_file_name);
%                 addAvailableData('mzl_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'mzl_trajectory', {'up', 'down'}, save_folder, save_file_name);
% 
%                 addAvailableData('fxr_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'fxr_trajectory', {'right', 'left'}, save_folder, save_file_name);
%                 addAvailableData('fyr_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'fyr_trajectory', {'forward', 'backward'}, save_folder, save_file_name);
%                 addAvailableData('fzr_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'fzr_trajectory', {'up', 'down'}, save_folder, save_file_name);
%                 addAvailableData('mxr_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'mxr_trajectory', {'right', 'left'}, save_folder, save_file_name);
%                 addAvailableData('myr_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'myr_trajectory', {'forward', 'backward'}, save_folder, save_file_name);
%                 addAvailableData('mzr_trajectory', 'time_forceplate', 'sampling_rate_forceplate', 'mzr_trajectory', {'up', 'down'}, save_folder, save_file_name);
% 
%                 addAvailableData('cop_trajectories', 'time_forceplate', 'sampling_rate_forceplate', '_cop_labels', cop_all_directions, save_folder, save_file_name);
                
                % save in new format
                save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'forceplateTrajectories.mat');
                save ...
                  ( ...
                    [save_folder filesep save_file_name], ...
                    'forceplate_trajectories', ...
                    'forceplate_labels', ...
                    'time', ...
                    'forceplate_directions', ...
                    'sampling_rate' ...
                  );

                addAvailableData('forceplate_trajectories', 'time', 'sampling_rate', '_forceplate_labels', '_forceplate_directions', save_folder, save_file_name);
                
                disp(['processed ' raw_forceplate_file_name ' and saved as ' [save_folder filesep save_file_name]])        
            end
        end
    end
end

