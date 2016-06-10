function ...
  [ ...
    Jacobian, ...
    correlation_c, ...
    correlation_p, ...
    step_response, ...
    stimulus_response ...
  ] = ...
  calculateStepResponse ...
  ( ...
    footPositionTrajectoriesStimulus, ...
    footPositionTrajectoriesReference, ...
    comPositionTrajectoriesStimulus, ...
    comPositionTrajectoriesReference, ...
    comVelocityTrajectoriesStimulus, ...
    comVelocityTrajectoriesReference, ...
    index_mid, ...
    index_end ...
  )
    if nargin < 8
        index_end = size(footPositionTrajectoriesStimulus, 1);
    end
    if nargin < 7
        index_mid = round(index_end/2);
    end
    
    % extract relevant information
    com_positions_mid_ref = comPositionTrajectoriesReference(index_mid, :);
    com_velocities_mid_ref = comVelocityTrajectoriesReference(index_mid, :);
    foot_positions_end_ref = footPositionTrajectoriesReference(index_end, :);
    com_positions_mid_stm = comPositionTrajectoriesStimulus(index_mid, :);
    com_velocities_mid_stm = comVelocityTrajectoriesStimulus(index_mid, :);
    foot_positions_end_stm = footPositionTrajectoriesStimulus(index_end, :);
    
    % calculate difference from reference mean
    com_positions_mid_delta_ref = com_positions_mid_ref - mean(com_positions_mid_ref);
    com_velocities_mid_delta_ref = com_velocities_mid_ref - mean(com_velocities_mid_ref);
    foot_positions_end_delta_ref = foot_positions_end_ref - mean(foot_positions_end_ref);
    com_positions_mid_delta_stm = com_positions_mid_stm - mean(com_positions_mid_ref);
    com_velocities_mid_delta_stm = com_velocities_mid_stm - mean(com_velocities_mid_ref);
    foot_positions_end_delta_stm = foot_positions_end_stm - mean(foot_positions_end_ref);
    
    % calculate regression coefficients
    input_matrix = ...
      [ ...
        com_positions_mid_delta_ref; ...
        com_velocities_mid_delta_ref; ...
      ];
    output_matrix = ...
      [ ...
        foot_positions_end_delta_ref; ...
      ];
    Jacobian = output_matrix * pinv(input_matrix);
    
    % evaluate regression
    [correlation_c, correlation_p] = corr(input_matrix', output_matrix');

    % remove expected part
    step_response = foot_positions_end_delta_stm;
    if correlation_p(2) > 0.05
        expected_response = Jacobian(1) * com_positions_mid_delta_stm;
    else
        expected_response = Jacobian(1) * com_positions_mid_delta_stm + Jacobian(2) * com_velocities_mid_delta_stm;
    end
    stimulus_response = step_response - expected_response;
    
%     figure; hold on
%     plot(step_response, 'x-')
%     plot(expected_response, 'x-')
%     plot(stimulus_response, 'o-')
%     legend('absolute response', 'expected response', 'stimulus response')
    
    
    
    
%     figure;
%     plot(com_positions_mid_delta_ref, foot_positions_end_delta_ref, 'x')
%     axis equal
% 
%     figure;
%     plot(com_velocities_mid_delta_ref, foot_positions_end_delta_ref, 'x')
%     axis equal




end