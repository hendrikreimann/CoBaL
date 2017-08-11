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

% generate data
sampling_rate = 100;
T = 10;
time = 1/sampling_rate : 1/sampling_rate : T;
x_clean = sin(time*2*pi)';
x_nan = x_clean;
x_nan(161 : 212) = NaN;
x_nan(530 : 556) = NaN;
x_nan(711 : 735) = NaN;

filter_order = 10;
cutoff_frequency = 5;
[b, a] = butter(filter_order, cutoff_frequency/(sampling_rate/2));	% set filter parameters for butterworth filter: 2=order of filter;

d1 = designfilt ...
       ( ...
         'lowpassiir', ...
         'FilterOrder', 6, ...
         'HalfPowerFrequency', 1, ...
         'DesignMethod', 'butter' ...
       );
%          'HalfPowerFrequency', 0.15, ...


y_clean = nanDeriveByTime(x_clean, 1/sampling_rate);
y_nan = nanDeriveByTime(x_nan, 1/sampling_rate);
% y_filt = nanDeriveByTime(nanfiltfilt(b, a, x_nan), 1/sampling_rate);
y_filt = nanDeriveByTime(nanfiltfilt(d1, x_nan), 1/sampling_rate);
z_clean = nanDeriveByTime(y_clean, 1/sampling_rate);
z_nan = nanDeriveByTime(y_nan, 1/sampling_rate);
% z_filt = nanDeriveByTime(nanfiltfilt(b, a, y_nan), 1/sampling_rate);
z_filt = nanDeriveByTime(nanfiltfilt(d1, y_nan), 1/sampling_rate);



figure; hold on
plot(time, x_clean)
plot(time, x_nan)
plot(time, y_clean)
plot(time, y_nan)
plot(time, z_filt)
plot(time, z_clean)
plot(time, z_nan)
plot(time, z_filt)









