function B = orthogonalizeBasis(A)

    [n_dim, n_vec] = size(A);
     
    Q = zeros(n_dim, n_vec);
    R_qr = zeros(n_dim, n_vec);
    for i_col = 1 : n_vec
        v = A(:, i_col);
        for i_row = 1 : i_col-1
            R_qr(i_row, i_col) = Q(:, i_row)' * A(:, i_col);
            v = v - R_qr(i_row, i_col) * Q(:, i_row);
        end
        R_qr(i_col, i_col) = norm(v);
        Q(:, i_col) = v/R_qr(i_col, i_col);
    end
    B = Q;
end