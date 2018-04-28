% function ...
%   [  ...
%     x_trajectory, ...
%     p_trajectory, ...
%     u_trajectory ...
%   ] = ...
function ...
  gvsModel_xiDotControlled ...
    ( ...
      alphaXiDot, ...
      gamma_p, ...
      gamma_v, ...
      gamma_a, ...
      sigma_motor, ...
      sigma_head_pos, ...
      sigma_head_vel, ...
      sigma_head_acc, ...
      time_step, ...
      T_total, ...
      rngSeed ...
    )

%     if nargin == 0
%         alphaXiDot = 5;
%         gamma_p = 1e-1;
%         gamma_v = 2e-2;
%         gamma_a = 5e-3;
%         sigma_motor = 1e-6;
%         sigma_head_pos = 4.9e-5;
%         sigma_head_vel = 3.7e-5;
%         sigma_head_acc = 3.7e-5;
%         time_step = 0.001;
%         T_total = 2;
%         rngSeed = 0;
%     end
    
    if nargin == 0
        alphaXiDot = 5;
        gamma_p = 1e-1;
        gamma_v = 1e-2;
        gamma_a = 1e-2;
        gamma_p = 1e-0;
        gamma_v = 1e-1;
        gamma_a = 1e-1;
        sigma_motor = 1e-5;
        sigma_head_pos = 4.9e-5;
        sigma_head_vel = 3.7e-5;
        sigma_head_acc = 3.7e-5;
        time_step = 0.001;
        T_total = 5;
        rngSeed = 0;
    end


    % set parameters
    number_of_dofs = 3;
    time = time_step : time_step : T_total;
    time_delay = 0.150; % in seconds
    number_of_delay_time_steps = round(time_delay / time_step);
    number_of_time_steps = length(time);
    gvs_bias_start = 0.5053;
    gvs_bias_ramp_time = 0.0;
    
    z_c = 0.814;
%     z_c = 0.6;
    g = 9.81;
    omega = sqrt(g/z_c);
    
    x_reference = 0.016;
    k_proportional = 0.5;
%     k_proportional = 0;

    % initialize system
    stream = RandStream('mt19937ar');
    reset(stream, rngSeed);
    RandStream.setGlobalStream(stream)

    A_system = [0 1 0; 0 0 1; 0 omega^2 0];
    B_system = [0; 0; -omega^2];
    C_system = [1 0 -1/omega^2];

    % control
    control_style = 'constant';
%     control_style = 'proportional';
    
    sigma_step = 12e-3;
    sigma_init = 7.3e-3;
    
%     sigma_step = 0;
%     sigma_init = 0;
    
%     sigma_step = 0;
    b_offset = 0.0055; % results in step width of ~5cm
%     T_step = 0.6247;
    T_step = 0.5053;
    T_n = T_step;
    t_n = T_step;
    n = 1;
    step_number_trajectory = zeros(1, number_of_time_steps);
    step_number_trajectory(1) = n;

    % trajectories
    gvs_bias_ramp_number_of_time_steps = gvs_bias_ramp_time * time_step^(-1);
    gvs_bias_ramp = sinusSigmoid(1:gvs_bias_ramp_number_of_time_steps, 1, gvs_bias_ramp_number_of_time_steps);
    gvs_bias_ramp_time_steps = (1:gvs_bias_ramp_number_of_time_steps) + round(gvs_bias_start*time_step^(-1)) - round(gvs_bias_ramp_number_of_time_steps/2);
    gvs_bias_p_trajectory = zeros(1, number_of_time_steps);
    gvs_bias_p_trajectory(time > gvs_bias_start) = gamma_p;
    gvs_bias_p_trajectory(gvs_bias_ramp_time_steps) = gvs_bias_ramp*gamma_p;
    gvs_bias_v_trajectory = zeros(1, number_of_time_steps);
    gvs_bias_v_trajectory(time > gvs_bias_start) = gamma_v;
    gvs_bias_v_trajectory(gvs_bias_ramp_time_steps) = gvs_bias_ramp*gamma_v;
    gvs_bias_a_trajectory = zeros(1, number_of_time_steps);
    gvs_bias_a_trajectory(time > gvs_bias_start) = gamma_a;
    gvs_bias_a_trajectory(gvs_bias_ramp_time_steps) = gvs_bias_ramp*gamma_a;
    
    % % TF: what exactly is happening here?
    % Kalman filter
    F_kalman = eye(3) + time_step * A_system;
    B_kalman = time_step * B_system;
    H = eye(3); % assume full observability

    R_kalman = diag([sigma_head_pos^2, sigma_head_vel^2, sigma_head_acc^2]);                                % covariance of the observation noise
    Q_kalman = B_kalman*B_kalman'*sigma_motor;                                                              % covariance of the process noise

    P_aposteriori = zeros(3);             % initial error covariance matrix
    
    % initialize
    x_trajectory = zeros(number_of_dofs, number_of_time_steps);
    p_trajectory = zeros(1, number_of_time_steps);
    u_trajectory = zeros(1, number_of_time_steps);
    r_trajectory = zeros(number_of_dofs, number_of_time_steps);
    xi_trajectory = zeros(1, number_of_time_steps);
    x_trajectory(1, 1) = sigma_init * randn(1);
    p_trajectory(1) = x_trajectory(1, 1);

    % initialize sensor variables
    z_trajectory = zeros(size(H, 1), number_of_time_steps);
    x_hat_delayed_trajectory = zeros(number_of_dofs, number_of_time_steps);
    z_delayed_trajectory = zeros(size(H, 1), number_of_time_steps);
    x_hat_predict_linear_trajectory = zeros(number_of_dofs, number_of_time_steps);
    x_hat_predict_nonlinear_trajectory = zeros(number_of_dofs, number_of_time_steps);
    xi_hat_predict_nonlinear_trajectory = zeros(1, number_of_time_steps);
    p_hat_predict_nonlinear_trajectory = zeros(1, number_of_time_steps);

    z_trajectory(:, 1) = H*x_trajectory(:, 1);
    x_hat_delayed_trajectory(:, 1) = x_trajectory(:, 1);
    z_delayed_trajectory(:, 1) = H*x_trajectory(:, 1);
    x_hat_predict_linear_trajectory(:, 1) = x_trajectory(:, 1);
    x_hat_predict_nonlinear_trajectory(:, 1) = x_trajectory(:, 1);

    P_kalman = P_aposteriori;
    P_kalman_delayed = P_aposteriori;
    P_kalman_delayed_trajectory = zeros(number_of_time_steps, 9);
    P_kalman_delayed_trajectory(1, :) = reshape(P_aposteriori, 1, 9);
    
    for i_time = 2 : number_of_time_steps
        % time
        t = time(i_time);
        i_time_delayed = max([i_time - number_of_delay_time_steps, 1]);

        % estimate
    %     x_estimate = x_trajectory(:, i_time-1);
        x_estimate = x_hat_predict_nonlinear_trajectory(:, i_time-1);

        % control
        c_rel_expected = r_trajectory(1, i_time-1) - p_hat_predict_nonlinear_trajectory(i_time-1);
        v_expected = r_trajectory(2, i_time-1);
        a_expected = r_trajectory(3, i_time-1);

        xi_expected = c_rel_expected + v_expected/omega;
        xi_dot_expected = v_expected + a_expected/omega;

        xi_hat_predict_nonlinear = x_hat_predict_nonlinear_trajectory(1, i_time-1) + x_hat_predict_nonlinear_trajectory(2, i_time-1)/omega;
        xi_dot_hat_predict_nonlinear = x_hat_predict_nonlinear_trajectory(2, i_time-1) + x_hat_predict_nonlinear_trajectory(3, i_time-1)/omega;

        u = alphaXiDot * (xi_dot_hat_predict_nonlinear - xi_dot_expected);
%         u = alphaXiDot * (xi_hat_predict_nonlinear - xi_expected);

        w = mvnrnd(B_kalman * u, Q_kalman, 1)';
        w_noiseless = B_kalman * u;

        % iterate real system
        x_trajectory(:, i_time) = F_kalman * x_trajectory(:, i_time-1) + w; % discrete version
        u_trajectory(:, i_time) = u;
        if i_time == number_of_delay_time_steps*2
    %             x_trajectory(1, i_time) = 0.01;
        end
        xi_trajectory(i_time) = x_trajectory(1, i_time) + x_trajectory(2, i_time)/omega;


        % iterate internal model
        r_trajectory(:, i_time) = F_kalman * r_trajectory(:, i_time-1);

        % output
        p_trajectory(i_time) = C_system * x_trajectory(:, i_time);
        p_hat_predict_nonlinear_trajectory(i_time) = C_system * x_hat_predict_nonlinear_trajectory(:, i_time);

        % step control
        if t > t_n(n)
            % step finished, calculate new step
            n = n+1;

            xi_hat = x_estimate(1) + 1/omega * x_estimate(2);
            
            % constant offset control
            if strcmp(control_style, 'constant')
                p_new = xi_hat - b_offset * (-1)^n;
            elseif strcmp(control_style, 'proportional')
                % offset plus proportional control
                p_new = xi_hat - b_offset * (-1)^n + k_proportional * (xi_hat - x_reference);
            else
                error('Control style not valid.')
            end
            if n > 2
                % apply step control noise
                p_new = p_new + sigma_step * randn(1);
            end
            
            a_new = omega^2*(x_trajectory(1, i_time) - p_new);
            a_new_desired = omega^2*(x_estimate(1) - p_new);

            % reset ankle torque modulation
    %         q = zeros(size(q));

            T_n(n) = T_step;
            t_n(n) = t + T_step;

            % update system
            x_trajectory(3, i_time) = a_new;
            r_trajectory(1, i_time) = x_estimate(1);
            r_trajectory(2, i_time) = x_estimate(2);
            r_trajectory(3, i_time) = a_new_desired;
        end
        step_number_trajectory(i_time) = n;





        % kalman observation and update
        v = randn(stream, 1, 1) * sqrt(R_kalman);
        v = mvnrnd(zeros(size(H, 1), 1), R_kalman)';

        z_trajectory(:, i_time) = H * x_trajectory(:, i_time) + v;
        z_trajectory(1, i_time) = z_trajectory(1, i_time) + gvs_bias_p_trajectory(i_time);
        z_trajectory(2, i_time) = z_trajectory(2, i_time) + gvs_bias_v_trajectory(i_time);
        z_trajectory(3, i_time) = z_trajectory(3, i_time) + gvs_bias_a_trajectory(i_time);

        % kalman observation and update - delayed
        z_delayed_trajectory(:, i_time) = z_trajectory(:, i_time_delayed);
        x_hat_delayed_trajectory(:, i_time) = F_kalman * x_hat_delayed_trajectory(:, i_time-1) + B_kalman * u_trajectory(:, i_time_delayed);            % predicted state estimate



        P_kalman_delayed = F_kalman * P_kalman_delayed * F_kalman' + Q_kalman;
        y_tilde = z_delayed_trajectory(:, i_time) - H * x_hat_delayed_trajectory(:, i_time);                    % measurement residual
        S_lqr = H * P_kalman_delayed * H' + R_kalman;                                                                 % innovation covariance
        K_kalman = P_kalman_delayed * H' * pinv(S_lqr);                                                                % Kalman gain
        x_hat_delayed_trajectory(:, i_time) = x_hat_delayed_trajectory(:, i_time) + K_kalman*y_tilde;    % aposteriori estimate
        P_kalman_delayed = (eye(number_of_dofs) - K_kalman * H) * P_kalman_delayed * (eye(number_of_dofs) - K_kalman * H)' + K_kalman*R_kalman*K_kalman';
        if i_time <= number_of_delay_time_steps + 1
            z_delayed_trajectory(:, i_time) = z_trajectory(:, 1);
            x_hat_delayed_trajectory(:, i_time) = x_hat_delayed_trajectory(:, 1);
        end

        % x_hat_delayed_trajectory is not making the step
        if step_number_trajectory(i_time_delayed+1) > step_number_trajectory(i_time_delayed)
            xi_hat = x_hat_delayed_trajectory(1, i_time-1) + 1/omega * x_hat_delayed_trajectory(2, i_time-1);
%             p_new = xi_hat - b_offset * (-1)^step_number_trajectory(i_time);
            
            if strcmp(control_style, 'constant')
                p_new = xi_hat - b_offset * (-1)^step_number_trajectory(i_time);
            elseif strcmp(control_style, 'proportional')
                % offset plus proportional control
                p_new = xi_hat - b_offset * (-1)^step_number_trajectory(i_time) + k_proportional * (xi_hat - x_reference);
            else
                error('Control style not valid.')
            end
            
            
            a_new = omega^2*(x_hat_delayed_trajectory(1, i_time-1) - p_new);
            x_hat_delayed_trajectory(3, i_time) = a_new;

    %         P_kalman_delayed = zeros(3);

        end

        % linear prediction
        x_dot_predict = A_system * x_hat_delayed_trajectory(:, i_time) + B_system * u_trajectory(:, i_time);
        x_hat_predict_linear_trajectory(:, i_time) = x_hat_delayed_trajectory(:, i_time) + number_of_delay_time_steps * time_step * x_dot_predict;

        % non-linear prediction
        x_predict_nonlinear = x_hat_delayed_trajectory(:, i_time);
        for j_time = i_time_delayed+1 : i_time - 1
    %         t_predicted = time(j_time);
            u_current = u_trajectory(j_time);
            x_dot_predict_nonlinear = A_system * x_predict_nonlinear + B_system * u_current;
            x_predict_nonlinear = x_predict_nonlinear + time_step*x_dot_predict_nonlinear;

            % was a step made here?
            if step_number_trajectory(j_time) > step_number_trajectory(j_time-1)
                xi_hat = x_predict_nonlinear(1) + 1/omega * x_predict_nonlinear(2);
                if strcmp(control_style, 'constant')
                    p_new = xi_hat - b_offset * (-1)^step_number_trajectory(j_time);
                elseif strcmp(control_style, 'proportional')
                    % offset plus proportional control
                    p_new = xi_hat - b_offset * (-1)^step_number_trajectory(j_time) + k_proportional * (xi_hat - x_reference);
                else
                    error('Control style not valid.')
                end     
                a_new = omega^2*(x_predict_nonlinear(1) - p_new);
                x_predict_nonlinear(3) = a_new;
            end

        end
        x_hat_predict_nonlinear_trajectory(:, i_time) = x_predict_nonlinear;    
        xi_hat_predict_nonlinear_trajectory(i_time) = x_hat_predict_nonlinear_trajectory(1, i_time) + x_hat_predict_nonlinear_trajectory(2, i_time)/omega;

        P_kalman_delayed_trajectory(i_time, :) = reshape(P_kalman_delayed, 1, 9);
    
    end

    visualize_trajectories = 1;
    if visualize_trajectories
        % lateral position
        set(0,'defaulttextinterpreter','latex')
        figure; axes; hold on; title('lateral position trajectories', 'fontsize', 24)
        plot(time, x_trajectory(1, :), 'linewidth', 3, 'displayname', '$c$');
    %     plot(time, xi_trajectory(1, :), 'linewidth', 3, 'displayname', '$\xi_1$');
    %     plot(time, p_trajectory_reference(1, :), 'linewidth', 1, 'displayname', '$p_1 - reference$');
%         plot(time, z_delayed_trajectory(1, :), 'displayname', '$z_d$');
        plot(time, x_hat_delayed_trajectory(1, :), 'displayname', '$\widehat c_d$');
        plot(time, x_hat_predict_nonlinear_trajectory(1, :), 'displayname', '$\widehat c_{nlp}$');
    %     plot(time, x_sensed_trajectory(1, :), 'displayname', '$\widehat x_1$');
    %     plot(time, z_delayed_trajectory(1, :), 'displayname', '$z_1$');
    %     plot(time, xi_trajectory(1, :), 'displayname', '$\xi$');
    %     plot(time, xi_hat_predict_nonlinear_trajectory(1, :), 'displayname', '$\widehat \xi$');
    %     plot(time, xi_dot_trajectory(1, :), 'displayname', '$\dot \xi_1$');
    %     plot(time, r_trajectory(1, :), 'linewidth', 3, 'displayname', '$r$');
        plot(time, p_trajectory, 'linewidth', 1, 'displayname', '$p$');
        my_legend = legend('show');
        set(my_legend, 'Interpreter', 'Latex', 'fontsize', 12);
        xlabel('time (s)', 'fontsize', 24)
        ylabel('lateral position', 'fontsize', 24)
    end




end