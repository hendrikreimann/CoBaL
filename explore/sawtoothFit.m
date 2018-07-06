

function [fit_data, fit_period, fit_offset_x, fit_offset_y, fit_amplitude] = sawtoothFit(time, data, period_init, offset_x_init, offset_y_init, amplitude_init)
% function sawtoothFit(time, data, period_init, offset_x_init, offset_y_init, amplitude_init)

    % generate initial fit to check this out
    fit_data_init = generateSawtooth(time, period_init, offset_x_init, offset_y_init, amplitude_init);
    parameters_init = [period_init, offset_x_init, offset_y_init, amplitude_init];
    
%     figure; hold on
%     plot(time, data);

    options = optimset ...
        ( ...
            'GradObj', 'off', ...
            'Display','off', ...
            'LargeScale', 'off', ...
            'DerivativeCheck', 'on', ...
            'UseParallel', 'always' ...
        );
%     parameters_fit = fminunc(@sawtoothDifference, parameters_init, options);
    lower_bound = [0, -inf, -inf, 0];
    upper_bound = [inf, inf, inf, inf]; 
    parameters_fit = fmincon(@sawtoothDifference, parameters_init, [], [], [], [], lower_bound, upper_bound, [], options);
    
    fit_period = parameters_fit(1);
    fit_offset_x = parameters_fit(2);
    fit_offset_y = parameters_fit(3);
    fit_amplitude = parameters_fit(4);
    
    fit_data = generateSawtooth(time, fit_period, fit_offset_x, fit_offset_y, fit_amplitude);
%     
% 
%     figure; hold on
%     plot(time, data);
% %     plot(time, fit_data_init, 'linewidth', 2);
%     plot(time, fit_data, 'linewidth', 2);


    function f = sawtoothDifference(parameters)
        period = parameters(1);
        offset_x = parameters(2);
        offset_y = parameters(3);
        amplitude = parameters(4);

        sawtooth = generateSawtooth(time, period, offset_x, offset_y, amplitude);
        
%         plot(time, sawtooth, 'linewidth', 2)

        f = mean((sawtooth - data).^2);
    end


end





