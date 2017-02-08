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

function [data_normalized, time_derivative_normalized] = normalizePeriodicVariable(data, time, peaks)
    data_normalized = zeros(size(data)) * NaN;
    time_derivative_normalized = zeros(size(data)) * NaN;
    
    % calculate and normalize velocity
    sampling_rate = median(diff(time))^(-1);
    filter_order = 4;
    cutoff_frequency = 10; % cutoff frequency, in Hz
    [b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;
    velocity = deriveByTime(nanfiltfilt(b, a, data), 1/sampling_rate);
    
    for i_peak = 1 : length(peaks)-1
        % normalize data
        data_stretch = data(peaks(i_peak) : peaks(i_peak+1));
        [~, peak_index] = max(data_stretch);
        data_upswing = data_stretch(1 : peak_index);
        data_downswing = data_stretch(peak_index : end);
        data_upswing_normalized = (data_upswing - min(data_upswing)) * 2 / (max(data_upswing) - min(data_upswing)) - 1;
        data_downswing_normalized = (data_downswing - min(data_downswing)) * 2 / (max(data_downswing) - min(data_downswing)) - 1;
        data_stretch_normalized = zeros(size(data_stretch));
        data_stretch_normalized(1 : peak_index) = data_upswing_normalized;
        data_stretch_normalized(peak_index : length(data_stretch)) = data_downswing_normalized;
        
        % normalize velocity
        velocity_stretch = velocity(peaks(i_peak) : peaks(i_peak+1));
        velocity_upswing = velocity_stretch(1 : peak_index);
        velocity_downswing = velocity_stretch(peak_index : end);
        [~, peak_index_upswing] = max(velocity_upswing);
        [~, peak_index_downswing] = min(velocity_downswing);
        velocity_upswing_accel = velocity_upswing(1 : peak_index_upswing);
        velocity_upswing_decel = velocity_upswing(peak_index_upswing : end);
        velocity_downswing_accel = velocity_downswing(1 : peak_index_downswing);
        velocity_downswing_decel = velocity_downswing(peak_index_downswing : end);
        
        velocity_upswing_accel_normalized = (velocity_upswing_accel - min(velocity_upswing_accel)) / (max(velocity_upswing_accel) - min(velocity_upswing_accel));
        velocity_upswing_decel_normalized = (velocity_upswing_decel - min(velocity_upswing_decel)) / (max(velocity_upswing_decel) - min(velocity_upswing_decel));
        velocity_downswing_accel_normalized = (velocity_downswing_accel - min(velocity_downswing_accel)) / (max(velocity_downswing_accel) - min(velocity_downswing_accel)) - 1;
        velocity_downswing_decel_normalized = (velocity_downswing_decel - min(velocity_downswing_decel)) / (max(velocity_downswing_decel) - min(velocity_downswing_decel)) - 1;
        
%         velocity_upswing_normalized = (velocity_upswing - min(velocity_upswing)) / (max(velocity_upswing) - min(velocity_upswing));
%         velocity_downswing_normalized = (velocity_downswing - min(velocity_downswing)) / (max(velocity_downswing) - min(velocity_downswing)) - 1;
%         velocity_stretch_normalized = zeros(size(velocity_stretch));
%         velocity_stretch_normalized(1 : peak_index) = velocity_upswing_normalized;
%         velocity_stretch_normalized(peak_index : length(velocity_stretch)) = velocity_downswing_normalized;
        
        velocity_stretch_normalized = zeros(size(data_stretch));
        velocity_stretch_normalized(1 : peak_index_upswing) = velocity_upswing_accel_normalized;
        velocity_stretch_normalized(peak_index_upswing : length(velocity_upswing)) = velocity_upswing_decel_normalized;

        velocity_stretch_normalized(length(velocity_upswing) + (1 : peak_index_downswing) - 1) = velocity_downswing_accel_normalized;
        velocity_stretch_normalized(length(velocity_upswing) + (peak_index_downswing : length(velocity_downswing)) - 1) = velocity_downswing_decel_normalized;
        
        % store
        data_normalized(peaks(i_peak) : peaks(i_peak+1)) = data_stretch_normalized;
        time_derivative_normalized(peaks(i_peak) : peaks(i_peak+1)) = velocity_stretch_normalized;
        
%         figure; hold on
% %         plot(time_stretch, data_stretch_normalized);
%         plot(1 : peak_index, velocity_upswing);
%         plot(peak_index : length(velocity_stretch), velocity_downswing);
%         plot(1 : peak_index_upswing, velocity_upswing_accel_normalized);
%         plot(peak_index_upswing : length(velocity_upswing), velocity_upswing_decel_normalized);
%         plot(length(velocity_upswing) + (1 : peak_index_downswing) - 1, velocity_downswing_accel_normalized);
%         plot(length(velocity_upswing) + (peak_index_downswing : length(velocity_downswing)) - 1, velocity_downswing_decel_normalized);
        
%         plot(time_stretch, velocity_stretch_normalized); hold on
    end
    
   
%     figure; hold on
%     plot(time, data_normalized);
%     plot(time, time_derivative_normalized);
end
