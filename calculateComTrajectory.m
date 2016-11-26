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
            virtual_marker_trajectories = zeros(size(marker_trajectories, 1), length(virtual_marker_headers)*3);
            for i_time = 1 : size(marker_trajectories, 1)
                % calculate virtual marker positions
                virtual_marker_trajectories(i_time, :) = ...
                    calculateVirtualMarkerPositions ...
                      ( ...
                        marker_reference, ...
                        marker_trajectories(i_time, :), ...
                        marker_headers, ...
                        virtual_marker_reference, ...
                        virtual_marker_headers ...
                      );
                
                
            end
            
            com_file_name = ['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'comTrajectories')];
            save ...
              ( ...
                com_file_name, ...
                'virtual_marker_trajectories' ...
              );

          toc
          return
        end
    end
end