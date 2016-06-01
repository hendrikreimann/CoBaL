

current_directory = pwd;
load('subjectInfo.mat');

%% process data from Labview
data_dir = dir('*_labviewTrajectories.mat');
clear file_name_list;
[file_name_list{1:length(data_dir)}] = deal(data_dir.name);

number_of_files = length(file_name_list);
for i_trial = 1 : number_of_files
    % load data
    raw_forceplate_file_name = file_name_list{i_trial};
    [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(raw_forceplate_file_name);
    load(raw_forceplate_file_name);

    % do not process if there is Nexus data available for this trial
    if ~exist(makeFileName(date, subject_id, trial_type, trial_number, 'forceplateTrajectoriesRaw.mat'), 'file')
        time_forceplate = time_labview;
        sampling_rate_forceplate = NaN;
        
        % extract and process data
        left_forceplate_wrench_Acl = [fxl_trajectory fyl_trajectory fzl_trajectory mxl_trajectory myl_trajectory mzl_trajectory];
        left_forceplate_cop_Acl = [copxl_trajectory copyl_trajectory zeros(size(copxl_trajectory))];
        right_forceplate_wrench_Acr = [fxr_trajectory fyr_trajectory fzr_trajectory mxr_trajectory myr_trajectory mzr_trajectory];
        right_forceplate_cop_Acr = [copxr_trajectory copyr_trajectory zeros(size(copxr_trajectory))];

        % define forceplate rotation and translation
        Acr_to_Acw_rotation = [-1 0 0; 0 1 0; 0 0 -1];
        Acr_to_Acw_translation = [0.5588; 0; 0];
        Acr_to_Acw_trafo = [Acr_to_Acw_rotation Acr_to_Acw_translation; 0 0 0 1];
        Acl_to_Acw_rotation = [-1 0 0; 0 1 0; 0 0 -1];
        Acl_to_Acw_translation = [-0.5588; 0; 0];
        Acl_to_Acw_trafo = [Acl_to_Acw_rotation Acl_to_Acw_translation; 0 0 0 1];
        Acr_to_Acw_adjoint = rigidToAdjointTransformation(Acr_to_Acw_trafo);
        Acl_to_Acw_adjoint = rigidToAdjointTransformation(Acl_to_Acw_trafo);

        % transform
        left_forceplate_wrench_Acw = (Acl_to_Acw_adjoint' * left_forceplate_wrench_Acl')';
        left_forceplate_cop_Acw = (eye(2, 4) * Acl_to_Acw_trafo * [left_forceplate_cop_Acl ones(size(left_forceplate_cop_Acl, 1), 1)]')';
        right_forceplate_wrench_Acw = (Acr_to_Acw_adjoint' * right_forceplate_wrench_Acr')';
        right_forceplate_cop_Acw = (eye(2, 4) * Acr_to_Acw_trafo * [right_forceplate_cop_Acr ones(size(right_forceplate_cop_Acr, 1), 1)]')';

        % calculate wrenches and CoP for complete plate
        total_forceplate_wrench_Acw = left_forceplate_wrench_Acw + right_forceplate_wrench_Acw;
        copx_trajectory = - total_forceplate_wrench_Acw(:, 5) ./ total_forceplate_wrench_Acw(:, 3);
        copy_trajectory = total_forceplate_wrench_Acw(:, 4) ./ total_forceplate_wrench_Acw(:, 3);
        total_forceplate_cop_Acw = [copx_trajectory copy_trajectory];
        
        % re-zero CoP for low loads
        left_forceplate_low_load_indicator = copxl_trajectory == 0;
        left_forceplate_cop_Acw(left_forceplate_low_load_indicator, :) = 0;
        right_forceplate_low_load_indicator = copxr_trajectory == 0;
        right_forceplate_cop_Acw(right_forceplate_low_load_indicator, :) = 0;
        total_forceplate_low_load_indicator = copx_trajectory == 0;
        total_forceplate_cop_Acw(total_forceplate_low_load_indicator, :) = 0;
        
        
%         % visualize
%         figure; axes; hold on
%         plot(time_forceplate, total_forceplate_cop_Acw(:, 1));
%         plot(time_forceplate, left_forceplate_cop_Acw(:, 1));
%         plot(time_forceplate, right_forceplate_cop_Acw(:, 1));
%         legend('total', 'left', 'right')
%         
%         figure; axes; hold on
%         plot(time_forceplate, total_forceplate_cop_Acw(:, 2));
%         plot(time_forceplate, left_forceplate_cop_Acw(:, 2));
%         plot(time_forceplate, right_forceplate_cop_Acw(:, 2));
%         legend('total', 'left', 'right')
        
        
        % save
        save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'forceplateTrajectories');
        save ...
          ( ...
            save_file_name, ...
            'left_forceplate_wrench_Acw', ...
            'left_forceplate_cop_Acw', ...
            'right_forceplate_wrench_Acw', ...
            'right_forceplate_cop_Acw', ...
            'total_forceplate_wrench_Acw', ...
            'total_forceplate_cop_Acw', ...
            'fxl_trajectory', ...
            'fyl_trajectory', ...
            'fzl_trajectory', ...
            'mxl_trajectory', ...
            'myl_trajectory', ...
            'mzl_trajectory', ...
            'copxl_trajectory', ...
            'copyl_trajectory', ...
            'fxr_trajectory', ...
            'fyr_trajectory', ...
            'fzr_trajectory', ...
            'mxr_trajectory', ...
            'myr_trajectory', ...
            'mzr_trajectory', ...
            'copxr_trajectory', ...
            'copyr_trajectory', ...
            'time_forceplate', ...
            'sampling_rate_forceplate' ...
          );
        disp(['processed ' raw_forceplate_file_name ' and saved as ' save_file_name])        
    end
end

%% process data from Nexus
data_dir = dir('*_forceplateTrajectoriesRaw.mat');
clear file_name_list;
[file_name_list{1:length(data_dir)}] = deal(data_dir.name);


number_of_files = length(file_name_list);
for i_trial = 1 : number_of_files

    % load data
    raw_forceplate_file_name = file_name_list{i_trial};
    [date, subject_id, trial_type, trial_number, file_type] = getFileParameters(raw_forceplate_file_name);
    load(raw_forceplate_file_name);
    
    % filter
    filter_order_low = 4;
    cutoff_frequency_low = 50; % in Hz
    [b_lowpass, a_lowpass] = butter(filter_order_low, cutoff_frequency_low/(sampling_rate_forceplate/2), 'low');
    forceplate_trajectories_filtered = filtfilt(b_lowpass, a_lowpass, forceplate_trajectories_raw);
    
    
    % extract
    fxl_trajectory = forceplate_trajectories_filtered(:, 1);
    fyl_trajectory = forceplate_trajectories_filtered(:, 2);
    fzl_trajectory = forceplate_trajectories_filtered(:, 3);
    mxl_trajectory = forceplate_trajectories_filtered(:, 4);
    myl_trajectory = forceplate_trajectories_filtered(:, 5);
    mzl_trajectory = forceplate_trajectories_filtered(:, 6);
    fxr_trajectory = forceplate_trajectories_filtered(:, 7);
    fyr_trajectory = forceplate_trajectories_filtered(:, 8);
    fzr_trajectory = forceplate_trajectories_filtered(:, 9);
    mxr_trajectory = forceplate_trajectories_filtered(:, 10);
    myr_trajectory = forceplate_trajectories_filtered(:, 11);
    mzr_trajectory = forceplate_trajectories_filtered(:, 12);
    
    % calculate CoP
    copxl_trajectory = - myl_trajectory ./ fzl_trajectory;
    copyl_trajectory = mxl_trajectory ./ fzl_trajectory;
    copxr_trajectory = - myr_trajectory ./ fzr_trajectory;
    copyr_trajectory = mxr_trajectory ./ fzr_trajectory;
    
    % zero CoP for no contact times
    fz_threshold = 60;
    copxl_trajectory(fzl_trajectory < fz_threshold) = 0;
    copyl_trajectory(fzl_trajectory < fz_threshold) = 0;
    copxr_trajectory(fzr_trajectory < fz_threshold) = 0;
    copyr_trajectory(fzr_trajectory < fz_threshold) = 0;
    
    % transform forceplate data to CoBaL world frame A_cw
    left_forceplate_wrench_Acl = [fxl_trajectory fyl_trajectory fzl_trajectory mxl_trajectory myl_trajectory mzl_trajectory];
    left_forceplate_cop_Acl = [copxl_trajectory copyl_trajectory zeros(size(copxl_trajectory))];
    right_forceplate_wrench_Acr = [fxr_trajectory fyr_trajectory fzr_trajectory mxr_trajectory myr_trajectory mzr_trajectory];
    right_forceplate_cop_Acr = [copxr_trajectory copyr_trajectory zeros(size(copxr_trajectory))];

    % define forceplate rotation and translation
    Acr_to_Acw_rotation = [-1 0 0; 0 1 0; 0 0 -1];
    Acr_to_Acw_translation = [0.5588; 0; 0];
    Acr_to_Acw_trafo = [Acr_to_Acw_rotation Acr_to_Acw_translation; 0 0 0 1];
    Acl_to_Acw_rotation = [-1 0 0; 0 1 0; 0 0 -1];
    Acl_to_Acw_translation = [-0.5588; 0; 0];
    Acl_to_Acw_trafo = [Acl_to_Acw_rotation Acl_to_Acw_translation; 0 0 0 1];
    Acr_to_Acw_adjoint = rigidToAdjointTransformation(Acr_to_Acw_trafo);
    Acl_to_Acw_adjoint = rigidToAdjointTransformation(Acl_to_Acw_trafo);

    % transform
    left_forceplate_wrench_Acw = (Acl_to_Acw_adjoint' * left_forceplate_wrench_Acl')';
    left_forceplate_cop_Acw = (eye(2, 4) * Acl_to_Acw_trafo * [left_forceplate_cop_Acl ones(size(left_forceplate_cop_Acl, 1), 1)]')';
    right_forceplate_wrench_Acw = (Acr_to_Acw_adjoint' * right_forceplate_wrench_Acr')';
    right_forceplate_cop_Acw = (eye(2, 4) * Acr_to_Acw_trafo * [right_forceplate_cop_Acr ones(size(right_forceplate_cop_Acr, 1), 1)]')';

    % calculate wrenches and CoP for complete plate
    total_forceplate_wrench_Acw = left_forceplate_wrench_Acw + right_forceplate_wrench_Acw;
    copx_trajectory = - total_forceplate_wrench_Acw(:, 5) ./ total_forceplate_wrench_Acw(:, 3);
    copy_trajectory = total_forceplate_wrench_Acw(:, 4) ./ total_forceplate_wrench_Acw(:, 3);
    total_forceplate_cop_Acw = [copx_trajectory copy_trajectory];

    % re-zero CoP for low loads
    left_forceplate_low_load_indicator = copxl_trajectory == 0;
    left_forceplate_cop_Acw(left_forceplate_low_load_indicator, :) = 0;
    right_forceplate_low_load_indicator = copxr_trajectory == 0;
    right_forceplate_cop_Acw(right_forceplate_low_load_indicator, :) = 0;
    total_forceplate_low_load_indicator = copx_trajectory == 0;
    total_forceplate_cop_Acw(total_forceplate_low_load_indicator, :) = 0;


%     % visualize
%     figure; axes; hold on
%     plot(time_forceplate, total_forceplate_cop_Acw(:, 1));
%     plot(time_forceplate, left_forceplate_cop_Acw(:, 1));
%     plot(time_forceplate, right_forceplate_cop_Acw(:, 1));
%     legend('total', 'left', 'right')
% 
%     figure; axes; hold on
%     plot(time_forceplate, total_forceplate_cop_Acw(:, 2));
%     plot(time_forceplate, left_forceplate_cop_Acw(:, 2));
%     plot(time_forceplate, right_forceplate_cop_Acw(:, 2));
%     legend('total', 'left', 'right')    
    
    % re-zero CoP for low loads
    left_forceplate_low_load_indicator = copxl_trajectory == 0;
    left_forceplate_cop_Acw(left_forceplate_low_load_indicator, :) = 0;
    right_forceplate_low_load_indicator = copxr_trajectory == 0;
    right_forceplate_cop_Acw(right_forceplate_low_load_indicator, :) = 0;
    
    % save
    save_file_name = makeFileName(date, subject_id, trial_type, trial_number, 'forceplateTrajectories');
    save ...
      ( ...
        save_file_name, ...
        'left_forceplate_wrench_Acw', ...
        'left_forceplate_cop_Acw', ...
        'right_forceplate_wrench_Acw', ...
        'right_forceplate_cop_Acw', ...
        'total_forceplate_wrench_Acw', ...
        'total_forceplate_cop_Acw', ...
        'fxl_trajectory', ...
        'fyl_trajectory', ...
        'fzl_trajectory', ...
        'mxl_trajectory', ...
        'myl_trajectory', ...
        'mzl_trajectory', ...
        'copxl_trajectory', ...
        'copyl_trajectory', ...
        'fxr_trajectory', ...
        'fyr_trajectory', ...
        'fzr_trajectory', ...
        'mxr_trajectory', ...
        'myr_trajectory', ...
        'mzr_trajectory', ...
        'copxr_trajectory', ...
        'copyr_trajectory', ...
        'time_forceplate', ...
        'sampling_rate_forceplate' ...
      );
    disp(['processed ' raw_forceplate_file_name ' and saved as ' save_file_name])        
    
    % delete raw data file
    delete(raw_forceplate_file_name);

end

