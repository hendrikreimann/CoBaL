function timeDerivative = deriveByTime(x, time)
    % disable warnings - we deal with NaNs explicitly here
    warning_id = 'MATLAB:chckxy:IgnoreNaN';
    warning_state_old = warning('query', warning_id);
    warning('off', warning_id)

    % make time a column if necessary
    if isrow(time)
        time = time';
    end
    
    if min(size(x)) == 1
        % only one variable
        if ~iscolumn(x)
            x = transpose(x);
        end
        
        if length(time) == 1
            % dt was provided instead of time, so generate a time vector
            time = (1 : length(x))' * time;
        end
        
    else
        % multiple variables
        if length(time) == 1
            % assume that time is the first data dimension
            time = (1 : size(x, 1))' * time;
        elseif length(time) == size(x, 2)
            % time vector was given, but time is the secodn data dimension, correct that
            x = transpose(x);
        elseif length(time) == size(x, 1)
            % all good
        else
            error('lenght of time vector does not correspond to size of the data');
        end
    end
        
    % derive
    timeDerivative = zeros(size(x));
    for i_dim = 1 : size(x, 2)

        dt = diff(time);
        dx = diff(x(:, i_dim));
        derivative = dx./dt;

        % resample
        time_resampling = time(1:end-1) + dt*0.5;
        derivative_resampled = spline(time_resampling, derivative, time);

        % remove NaN entries
        derivative_resampled(isnan(x(:, i_dim))) = NaN;

        timeDerivative(:, i_dim) = derivative_resampled;

    end
        

    % re-enable warnings if they were enabled before
    if strcmp(warning_state_old.state, 'on')
        warning('on', warning_id);
    end

    
    
    
    
    
    
    
%     if size(x, 1) == 1 & size(x, 2) > 1
%         x = x';
%     end
%     
%     if size(x, 2) > 1
%         timeDerivetive = zeros(size(x));
%         for i_dim = 1 : size(x, 2)
%             timeDerivative(:, i_dim) = deriveByTime(x(:, i_dim), time);
%         end
%         return
%     end
% 
%     % check whether time is a vector or a scalar
%     if length(time) == 1
%         time = (1 : length(x))' * time;
%     end
%     dt = diff(time);
%     dx = diff(x);
%     derivative = dx./dt;
% 
%     % resample
%     time_resampling = time(1:end-1) + dt*0.5;
%     derivative_resampled = spline(time_resampling, derivative, time);
%     
%     % remove NaN entries
%     derivative_resampled(isnan(x)) = NaN;
%     
% %     figure; hold on;
% %     plot(time_resampling, derivative, 'linewidth', 2);
% %     plot(time, derivative_resampled)
%     
%     
%     timeDerivative = derivative_resampled;
    
end

