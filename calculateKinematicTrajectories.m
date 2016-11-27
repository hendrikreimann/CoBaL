function calculateKinematicTrajectories(varargin)
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

            com_labels = [segment_labels 'BODY'];
            for i_label = 1 : length(com_labels)
                com_labels{i_label} = [com_labels{i_label} 'COM'];
            end
            
            % calculate
            joint_center_trajectories = zeros(size(marker_trajectories, 1), length(joint_center_headers)*3);
            com_trajectories = zeros(size(marker_trajectories, 1), length(com_labels)*3);
            for i_time = 1 : size(marker_trajectories, 1)
                % calculate joint center positions
                marker_current = marker_trajectories(i_time, :);
                joint_center_current = ...
                    calculateJointCenterPositions ...
                      ( ...
                        marker_reference, ...
                        marker_current, ...
                        marker_headers, ...
                        joint_center_reference, ...
                        joint_center_headers ...
                      );
                joint_center_trajectories(i_time, :) = joint_center_current;
                
                % calculate segment centers of mass
                number_of_segments = length(segment_coms_mcs);
                mcs_to_wcs_transformations = calculateMcsToWcsTransformations([marker_current joint_center_current], [marker_headers joint_center_headers], markers_by_segment);
                segment_coms_wcs = cell(number_of_segments, 1);
                for i_segment = 1 : number_of_segments
                    segment_com_mcs = [segment_coms_mcs{i_segment}; 1];
                    T_mcs_to_wcs = mcs_to_wcs_transformations{i_segment};
                    segment_coms_wcs{i_segment} = eye(3, 4) * T_mcs_to_wcs * segment_com_mcs;
                end
                
                % calculate whole body center of mass
                body_com = [0; 0; 0];
                for i_segment = 1 : number_of_segments
                    body_com = body_com + segment_masses(i_segment) * segment_coms_wcs{i_segment};
                end
                body_com = body_com * 1 / sum(segment_masses);
                
                % export centers of mass
                for i_segment = 1 : number_of_segments
                    com_trajectories(i_time, (i_segment-1)*3+1 : (i_segment-1)*3+3) = segment_coms_wcs{i_segment};
                end
                com_trajectories(i_time, number_of_segments*3+1 : number_of_segments*3+3) = body_com;
            end
            
            com_file_name = ['processed' filesep makeFileName(date, subject_id, condition, i_trial, 'kinematicTrajectories')];
            save ...
              ( ...
                com_file_name, ...
                'joint_center_trajectories', ...
                'com_trajectories', ...
                'com_labels' ...
              );

          toc
          return
        end
    end
end