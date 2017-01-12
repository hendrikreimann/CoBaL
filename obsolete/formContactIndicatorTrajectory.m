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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.% compare the kinematic tree against the kinematic chain

function contactIndicator = formContactIndicatorTrajectory(liftoffIndices, touchdownIndices, numberOfTimeSteps)

    
    if isempty(liftoffIndices) && isempty(touchdownIndices)
        contactIndicator = ones(numberOfTimeSteps, 1);
    else
        contactIndicator = zeros(numberOfTimeSteps, 1);
        
        % check for emptyness
        if isempty(liftoffIndices) && isempty(touchdownIndices)
            return
        elseif isempty(liftoffIndices) && ~isempty(touchdownIndices)
            contactIndicator(touchdownIndices(1):end) = 1;
            return
        elseif ~isempty(liftoffIndices) && isempty(touchdownIndices)
            contactIndicator(1:liftoffIndices(1)) = 1;
            return
        end
        
        % check first one
        if liftoffIndices(1) < touchdownIndices(1)
            contactIndicator(1:liftoffIndices(1)) = 1;
        end
        % mark stretch between touchdowns and liftoff as contact
        for i_touchdown = 1 : length(touchdownIndices)
            liftoff_indices_relative = liftoffIndices-touchdownIndices(i_touchdown);
            liftoff_indices_relative(liftoff_indices_relative<0) = inf;
            next_liftoff = min(liftoff_indices_relative) + touchdownIndices(i_touchdown);
            if isinf(next_liftoff)
                next_liftoff = numberOfTimeSteps;
            end
            contactIndicator(touchdownIndices(i_touchdown) : next_liftoff) = 1;
        end
    end



end

%#ok<*AGROW>
