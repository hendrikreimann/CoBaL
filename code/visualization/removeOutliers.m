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


function [data_pruned, outlier_indices] = removeOutliers(data, percentile)

    if nargin<2
        percentile = 25;
    end

    % remove outliers
    data_quartile_1 = prctile(data, percentile);
    data_quartile_3 = prctile(data, 100-percentile);
    data_iqr = iqr(data);
    data_upper_inner_fence_no = data_quartile_3 + 1.5*data_iqr;
    data_lower_inner_fence_no = data_quartile_1 - 1.5*data_iqr;
    outlier_indices = data>data_upper_inner_fence_no | data<data_lower_inner_fence_no;
    data_pruned = data(~outlier_indices);

end