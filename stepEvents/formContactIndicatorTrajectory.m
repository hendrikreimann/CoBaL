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
