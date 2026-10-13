function [hp, all_hp] = hermite_poly(t, n)
% HERMITE_POLY Generates a Hermite polynomial of order n on the vector t.
%   [hp, all_hp] = hermite_poly(t, n)
%
%   Args:
%       t (double): A real-valued vector on which to compute the Hermite 
%                   polynomials.
%       n (int): Order of the Hermite polynomial.
%
%   Returns:
%       hp (double): Returns an array with the Hermite polynomial computed on t.
%       all_hp (double): Returns all Hermite polynomials up to order n.

    num_t = length(t);
    all_hp = zeros(n + 1, num_t);

    if n == 0
        hp = ones(1, num_t);
        all_hp(1, :) = hp;
    else
        hp = 2 * t;
        all_hp(1, :) = ones(1, num_t);
        all_hp(2, :) = hp;

        for i = 3:(n + 1)
            % Recurrence relation: H_i(x) = 2x*H_{i-1}(x) - 2(i-1)*H_{i-2}(x)
            % Note: indices in MATLAB are 1-based, so i-1 corresponds to H_{i-2}
            % and i-2 corresponds to H_{i-3} from the Python code logic.
            hp = 2 * t .* all_hp(i - 1, :) - 2 * (i - 2) * all_hp(i - 2, :);
            all_hp(i, :) = hp;
        end
    end
    
    % The last calculated hp is the one of order n.
    hp = all_hp(n+1, :);
end