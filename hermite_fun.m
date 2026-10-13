function [hfunc_final, hfunc_all] = hermite_fun(N, q, t, T)
% HERMITE_FUN Computes an Hermite function of order q.
%   [hfunc_final, hfunc_all] = hermite_fun(N, q, t, T)
%
%   Args:
%       N (int): Length of the function in samples.
%       q (int): Order of the Hermite function.
%       t (double, optional): Values on which to compute the function.
%           Defaults to a centered vector if not provided.
%       T (double, optional): Scale of the Gaussian function.
%           Defaults to sqrt(N) if not provided.
%
%   Returns:
%       hfunc_final (double): The Hermite function of order q.
%       hfunc_all (double): All Hermite functions up to order q.

    % Handle optional arguments
    if nargin < 3 || isempty(t)
        t = (0:(N-1)) - floor(N/2);
    end
    if nargin < 4 || isempty(T)
        T = sqrt(N);
    end

    % Call the Hermite polynomial function
    [~, all_hp] = hermite_poly(sqrt(2 * pi) * t / T, q - 1);

    % Compute the Gaussian component
    gaussian_basic = exp(-pi * (t / T).^2) / sqrt(T);

    % Initialize the output matrix for all functions
    hfunc_all = zeros(q, length(t));

    % Loop to calculate each Hermite function
    for i = 1:q
        % Normalization constant
        Cnorm = sqrt(factorial(i - 1) * (2^(i - 1 - 0.5)));
        
        % Check if Cnorm is not zero to avoid division by zero
        if Cnorm > 1e-10
            hfunc_all(i, :) = gaussian_basic .* all_hp(i, :) / Cnorm;
        else
            hfunc_all(i, :) = 0; % Or handle as needed
        end
    end
    
    % Return the final function (order q-1 in the polynomial)
    hfunc_final = hfunc_all(q, :);
end