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
    forceplate_table = study_settings.getTable('forceplate_table', 1);
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
                filter_order_low = study_settings.get('force_plate_filter_order', 1);
                cutoff_frequency_low = study_settings.get('force_plate_filter_cutoff', 1);
                if ~isempty(filter_order_low) && ~isempty(cutoff_frequency_low)
                    [b_lowpass, a_lowpass] = butter(filter_order_low, cutoff_frequency_low/(loaded_data.sampling_rate_forceplate/2), 'low');
                    forceplate_trajectories_groomed = nanfiltfilt(b_lowpass, a_lowpass, raw_forceplate_trajectories_without_bad);
                else
                    forceplate_trajectories_groomed = raw_forceplate_trajectories;
                end
                
                % collect and sort
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
                    this_plate_copx = - this_plate_my ./ this_plate_fz;
                    this_plate_copy =   this_plate_mx ./ this_plate_fz;
                    this_plate_low_load_indicator = (abs(this_plate_fz) < fz_threshold);
                    this_plate_copx(this_plate_low_load_indicator, :) = 0;
                    this_plate_copy(this_plate_low_load_indicator, :) = 0;

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
                        study_settings.get('direction_x_pos', 1), ...
                        study_settings.get('direction_y_pos', 1), ...
                        study_settings.get('direction_z_pos', 1), ....
                        study_settings.get('direction_x_pos', 1), ...
                        study_settings.get('direction_y_pos', 1), ...
                        study_settings.get('direction_z_pos', 1); ....
                        study_settings.get('direction_x_neg', 1), ...
                        study_settings.get('direction_y_neg', 1), ...
                        study_settings.get('direction_z_neg', 1), ....
                        study_settings.get('direction_x_neg', 1), ...
                        study_settings.get('direction_y_neg', 1), ...
                        study_settings.get('direction_z_neg', 1); ....
                      };
                    this_cop_directions = ...
                      { ...
                        study_settings.get('direction_x_pos', 1), ...
                        study_settings.get('direction_y_pos', 1); ...
                        study_settings.get('direction_x_neg', 1), ...
                        study_settings.get('direction_y_neg', 1), ...
                      };


                    % store
                    forceplate_trajectories = [forceplate_trajectories this_plate_wrench this_plate_copx this_plate_copy]; %#ok<AGROW>
                    forceplate_labels = [forceplate_labels this_plate_labels]; %#ok<AGROW>
                    forceplate_directions = [forceplate_directions this_wrench_directions this_cop_directions]; %#ok<AGROW>
                end

                if strcmp(study_settings.get('combined_forceplate_data_source', 1), 'none')
                    
                elseif strcmp(study_settings.get('combined_forceplate_data_source', 1), 'from_data')
                    combined_plate_labels = study_settings.get('combined_forceplate_labels');
                    this_plate_fx = forceplate_trajectories_groomed(:, strcmp(loaded_data.forceplate_labels, combined_plate_labels{1}));
                    this_plate_fy = forceplate_trajectories_groomed(:, strcmp(loaded_data.forceplate_labels, combined_plate_labels{2}));
                    this_plate_fz = forceplate_trajectories_groomed(:, strcmp(loaded_data.forceplate_labels, combined_plate_labels{3}));
                    this_plate_mx = forceplate_trajectories_groomed(:, strcmp(loaded_data.forceplate_labels, combined_plate_labels{4}));
                    this_plate_my = forceplate_trajectories_groomed(:, strcmp(loaded_data.forceplate_labels, combined_plate_labels{5}));
                    this_plate_mz = forceplate_trajectories_groomed(:, strcmp(loaded_data.forceplate_labels, combined_plate_labels{6}));
                    total_forceplate_wrench = [this_plate_fx this_plate_fy this_plate_fz this_plate_mx this_plate_my this_plate_mz];

                    % calculate CoP
                    copx_total = forceplate_trajectories_groomed(:, strcmp(loaded_data.forceplate_labels, combined_plate_labels{7}));
                    copy_total = forceplate_trajectories_groomed(:, strcmp(loaded_data.forceplate_labels, combined_plate_labels{8}));
                    
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
                        study_settings.get('direction_x_pos', 1), ...
                        study_settings.get('direction_y_pos', 1), ...
                        study_settings.get('direction_z_pos', 1), ....
                        study_settings.get('direction_x_pos', 1), ...
                        study_settings.get('direction_y_pos', 1), ...
                        study_settings.get('direction_z_pos', 1); ....
                        study_settings.get('direction_x_neg', 1), ...
                        study_settings.get('direction_y_neg', 1), ...
                        study_settings.get('direction_z_neg', 1), ....
                        study_settings.get('direction_x_neg', 1), ...
                        study_settings.get('direction_y_neg', 1), ...
                        study_settings.get('direction_z_neg', 1); ....
                      };
                    total_cop_directions = ...
                      { ...
                        study_settings.get('direction_x_pos', 1), ...
                        study_settings.get('direction_y_pos', 1); ...
                        study_settings.get('direction_x_neg', 1), ...
                        study_settings.get('direction_y_neg', 1), ...
                      };

                    % store total wrench and CoP
                    forceplate_trajectories = [total_forceplate_wrench copx_total copy_total forceplate_trajectories]; %#ok<AGROW>
                    forceplate_labels = [total_labels forceplate_labels]; %#ok<AGROW>
                    forceplate_directions = [total_wrench_directions total_cop_directions forceplate_directions]; %#ok<AGROW>                    
                elseif strcmp(study_settings.get('combined_forceplate_data_source', 1), 'calculate')
                    % calculate total wrench and CoP
                    total_forceplate_wrench = wrenches{1};
                    for i_forceplate = 2 : number_of_forceplates
                        total_forceplate_wrench = total_forceplate_wrench + wrenches{i_forceplate};
                    end
                    copx_total = - total_forceplate_wrench(:, 5) ./ total_forceplate_wrench(:, 3);
                    copy_total = total_forceplate_wrench(:, 4) ./ total_forceplate_wrench(:, 3);
                    low_load_indicator = (abs(total_forceplate_wrench(:, 3)) < fz_threshold);
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
                        study_settings.get('direction_x_pos', 1), ...
                        study_settings.get('direction_y_pos', 1), ...
                        study_settings.get('direction_z_pos', 1), ....
                        study_settings.get('direction_x_pos', 1), ...
                        study_settings.get('direction_y_pos', 1), ...
                        study_settings.get('direction_z_pos', 1); ....
                        study_settings.get('direction_x_neg', 1), ...
                        study_settings.get('direction_y_neg', 1), ...
                        study_settings.get('direction_z_neg', 1), ....
                        study_settings.get('direction_x_neg', 1), ...
                        study_settings.get('direction_y_neg', 1), ...
                        study_settings.get('direction_z_neg', 1); ....
                      };
                    total_cop_directions = ...
                      { ...
                        study_settings.get('direction_x_pos', 1), ...
                        study_settings.get('direction_y_pos', 1); ...
                        study_settings.get('direction_x_neg', 1), ...
                        study_settings.get('direction_y_neg', 1), ...
                      };

                    % store total wrench and CoP
                    forceplate_trajectories = [total_forceplate_wrench copx_total copy_total forceplate_trajectories]; %#ok<AGROW>
                    forceplate_labels = [total_labels forceplate_labels]; %#ok<AGROW>
                    forceplate_directions = [total_wrench_directions total_cop_directions forceplate_directions]; %#ok<AGROW>
                else
                    warning(['Cannot interpret study settings entry "combined_forceplate_data_source: ' study_settings.get('combined_forceplate_data_source', 1) '", combined forcplate data not processed' ])
                end
                
                % extract
                time = loaded_data.time_forceplate;
                sampling_rate = loaded_data.sampling_rate_forceplate;
                
                % save
                save_folder = 'processed';
                
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

