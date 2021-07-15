function [ x ] = frequencySmoothing_multisine( y )
%SiPsmooth_more1 - Spectrum smoothing by averaging adjacent frequency points with
%an increasing number of points averaged with increasing frequency. Input is
%a vector of 10 vlaues. Return is a column vector of 10 values.
%
%   
% x(1:16)
% x=[y(1) mean(y(1:2)) y(2) mean(y(2:3)) mean(y(3:4)) mean(y(4:5)) mean(y(5:7)) mean(y(6:9)) mean(y(8:11)) mean(y(11:14)),...
% 			mean(y(15:18)) mean(y(18:23)) mean(y(23:29)) mean(y(29:35)) mean(y(35:45)),...
% 			mean(y(45:60))]; 
x = ...
  [ ...
    y(1) ...                1
    mean(y(1:2)) ...        2
    mean(y(1:3)) ...        3
    mean(y(1:4)) ...        4
    mean(y(2:4)) ...        5
    mean(y(3:5)) ...        6
    mean(y(4:6)) ...        7
    mean(y(5:8)) ...        8
    mean(y(6:9)) ...        9
    mean(y(8:10)) ...       10
  ]; 

end

