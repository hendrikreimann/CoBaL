function calculateComTrajectory(varargin)
    [condition_list, trial_number_list] = parseTrialArguments(varargin{:});
    
    load('subjectInfo.mat', 'date', 'subject_id');
    load('subjectModel.mat');
    
    for i_condition = 1 : length(condition_list)
        trials_to_process = trial_number_list{i_condition};
        for i_trial = trials_to_process
            % load data
            tic
            
            condition = condition_list{i_condition};
            load(['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'markerTrajectories')]);

            % calculate
            joint_center_trajectories = zeros(size(marker_trajectories, 1), length(joint_center_headers)*3);
            for i_time = 1 : size(marker_trajectories, 1)
                % calculate virtual marker positions
                joint_center_trajectories(i_time, :) = ...
                    calculateVirtualMarkerPositions ...
                      ( ...
                        marker_reference, ...
                        marker_trajectories(i_time, :), ...
                        marker_headers, ...
                        joint_center_reference, ...
                        joint_center_headers ...
                      );
                
                
            end
            
            com_file_name = ['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'comTrajectories')];
            save ...
              ( ...
                com_file_name, ...
                'joint_center_trajectories' ...
              );

          toc
          return
        end
    end
end