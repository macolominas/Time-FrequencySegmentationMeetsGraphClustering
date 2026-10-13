function X = sampleSTFT(S, M, t, f, tau)
% sampleSTFT - Generate a TF point cloud from a thresholded STFT
%
% INPUTS:
%   S   : STFT matrix (frequency x time)
%   M   : number of points to sample
%   t   : time vector (1 x Nt)
%   f   : frequency vector (1 x Nf)
%   tau : STFT magnitude threshold
%
% OUTPUT:
%   X   : M x 2 matrix
%         Each row is [time, frequency]
%
% The probability of selecting a TF point is proportional
% to the magnitude of the STFT, considering only points
% whose magnitude is larger than tau.
%
% Sampling is performed with replacement.

    % Dimensions
    [Nf, Nt] = size(S);

    if length(t) ~= Nt
        error('Length of t must match the number of columns of S.');
    end

    if length(f) ~= Nf
        error('Length of f must match the number of rows of S.');
    end

    if tau < 0
        error('tau must be non-negative.');
    end

    % STFT magnitude
    A = abs(S);

    % Keep only coefficients above the threshold
    mask = A > tau;

    if ~any(mask(:))
        error('No STFT coefficients are above the threshold tau.');
    end

    % Linear indices of the surviving TF points
    idx_valid = find(mask);

    % Corresponding amplitudes
    w = A(idx_valid);

    % Normalize weights to obtain a probability distribution
    p = w / sum(w);

    % Cumulative distribution
    cdf = cumsum(p);
    cdf(end) = 1;   % Avoid numerical problems

    % Generate M uniform random numbers
    r = rand(M,1);

    % Inverse transform sampling
    % histc returns the bin containing each random number.
    [~, idx_sample] = histc(r, [0; cdf]);

    % Protect against the case r == 1
    idx_sample(idx_sample == 0) = 1;
    idx_sample(idx_sample > length(idx_valid)) = length(idx_valid);

    % Map back to the original STFT matrix
    idx = idx_valid(idx_sample);

    % Convert linear indices to frequency/time indices
    [ifreq, itime] = ind2sub([Nf, Nt], idx);

    % Construct the point cloud
    X = [t(itime).' f(ifreq).'];

end
