function [fit_data, fit_period, fit_offset_x, fit_offset_y, fit_amplitude] = sinusoidFit(time, data, period_init, offset_x_init, offset_y_init, amplitude_init)
% function sinusoidFit(time, data, period_init, offset_x_init, offset_y_init, amplitude_init)

    % generate initial fit to check this out
    fit_data_init = generateSinusoid(time, period_init, offset_x_init, offset_y_init, amplitude_init);
    parameters_init = [offset_x_init, offset_y_init, amplitude_init];
    
%     figure; hold on
%     plot(time, data); 
    

    options = optimset ...
        ( ...
            'GradObj', 'off', ...
            'Display','off', ...
            'LargeScale', 'off', ...
            'DerivativeCheck', 'on' ...
        );
%     parameters_fit = fminunc(@sawtoothDifference, parameters_init, options);
    % time_offset, amplitude_offset, amplitude_gain
    lower_bound = [-inf, -inf, 0];
    upper_bound = [inf, inf, inf]; 
    parameters_fit = fmincon(@sinusoidDifference, parameters_init, [], [], [], [], lower_bound, upper_bound, [], options);
    
    fit_period = 70;
    fit_offset_x = parameters_fit(1);
    fit_offset_y = parameters_fit(2);
    fit_amplitude = parameters_fit(3);
    
    fit_data = generateSinusoid(time, fit_period, fit_offset_x, fit_offset_y, fit_amplitude);
%     
% 
%     figure; hold on
%     plot(time, data);
%     plot(time, fit_data, 'linewidth', 2);


    function f = sinusoidDifference(parameters)
        period = 70;
        offset_x = parameters(1);
        offset_y = parameters(2);
        amplitude = parameters(3);

        sinusoid = generateSinusoid(time, period, offset_x, offset_y, amplitude);
        
%         plot(time, sinusoid, 'linewidth', 2)

        f = mean((sinusoid - data).^2);
    end


end





