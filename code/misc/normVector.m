% norm a vector

function normedVector = normVector(vector)
    normedVector = zeros(size(vector));
    for i_col = 1 : size(vector, 2)
        normedVector(:, i_col) = vector(:, i_col) * norm(vector(:, i_col))^(-1);
    end
end

