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

function [V_para, V_perp, S, k_para, k_perp, E_para, E_perp, SigmaHat] = calculateUcmVariance(data, taskJacobian)
%
% Input
% data = (numberOfDofs x numberOfTrials) matrix containing numberOfTrials rows of observations
%
% Output
% V_para = variance parallel to the UCM (GEV, normalized by DoF)
% V_perp = variance perpendicular to the UCM (NGEV, normalized by DoF)
% S = UCM signature, ratio of normalized variance parallel and perpendicular to the UCM
% 
% ---------------------------------------------------------------------

numberOfDofs = size(data, 1);
numberOfTrials = size(data, 2);

% estimate covariance matrix
dataMeanFree = data - repmat(mean(data, 2), 1, size(data, 2));
SigmaHat = 1/numberOfTrials * (dataMeanFree*dataMeanFree');

% calculate bases of UCM and orthogonal space
k_perp = size(taskJacobian, 1);                             % number of controlled DoFs
k_para = numberOfDofs - k_perp;                             % number of uncontrolled DoFs
[~, ~, V] = svd(taskJacobian);                              % use singular value decomposition to get the basis
E_perp = V(:, 1:k_perp);                                    % E_perp contains the base vectors of the range space
E_para = V(:, k_perp+1:end);                                % E_para contains the base vectors of the null space
V_perp = 1/k_perp * trace(E_perp' * SigmaHat * E_perp);
V_para = 1/k_para * trace(E_para' * SigmaHat * E_para);
S = V_para / V_perp;


