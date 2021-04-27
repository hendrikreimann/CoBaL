function confidence_interval_radius = cinv(data, dim)
    if nargin < 2
        dim = 1;
    end
    N = size(data, dim);
    confidence_interval_radius = tinv(0.975, N-1) * nanstd(data, 1, dim)/sqrt(N);
end