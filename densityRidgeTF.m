function [ridge,lambda2,V2,gradS,H] = densityRidgeTF(S,t,f,threshold)
%DENSITYRIDGETF Density ridge extraction directly from a TF representation.
%
%   [ridge,lambda2,V2,gradS,H] = densityRidgeTF(S,t,f,threshold)
%
%   INPUTS
%       S         : Nf x Nt non-negative TF representation
%       t         : 1 x Nt time vector
%       f         : 1 x Nf frequency vector
%       threshold : relative threshold (e.g. 0.01)
%
%   OUTPUTS
%       ridge     : Nf x Nt logical matrix containing ridge points
%       lambda2   : Nf x Nt smallest Hessian eigenvalue
%       V2        : Nf x Nt x 2 eigenvector associated with lambda2
%       gradS     : Nf x Nt x 2 gradient [dS/dt,dS/df]
%       H         : Nf x Nt x 2 x 2 Hessian
%
%   The density ridge of dimension 1 is defined by
%
%       v2' * grad(p) = 0
%
%   together with
%
%       lambda2 < 0,
%
%   where v2 is the eigenvector associated with the smallest
%   eigenvalue of the Hessian.
%
%   The input TF representation is normalized and interpreted
%   as a density.
%
%   -------------------------------------------------------------------

%% Check input

S = double(S);

if any(S(:) < 0)
    error('S must be non-negative.');
end

[Nf,Nt] = size(S);

if length(t) ~= Nt
    error('Length of t must match the number of columns of S.');
end

if length(f) ~= Nf
    error('Length of f must match the number of rows of S.');
end

%% Normalize TFR

p = S ./ sum(S(:));

%% Derivatives

% First derivatives
[dp_df,dp_dt] = gradient(p);%,f,t);

gradS = zeros(Nf,Nt,2);
gradS(:,:,1) = dp_dt;
gradS(:,:,2) = dp_df;

% Second derivatives
[d2p_df2, d2p_dfdt] = gradient(dp_df);%,f,t);
[~,d2p_dt2]          = gradient(dp_dt);%,f,t);

%% Hessian

H = zeros(Nf,Nt,2,2);

H(:,:,1,1) = d2p_dt2;
H(:,:,1,2) = d2p_dfdt;
H(:,:,2,1) = d2p_dfdt;
H(:,:,2,2) = d2p_df2;

%% Eigenvalues/eigenvectors

lambda2 = zeros(Nf,Nt);
V2      = zeros(Nf,Nt,2);

for k = 1:Nt
    for l = 1:Nf

        Hlocal = squeeze(H(l,k,:,:));

        [V,D] = eig(Hlocal);

        eigvals = diag(D);

        % Smallest eigenvalue
        [lambda_min,idx] = min(eigvals);

        lambda2(l,k) = lambda_min;

        v = V(:,idx);

        % Store eigenvector
        V2(l,k,1) = v(1);
        V2(l,k,2) = v(2);
    end
end

%% Ridge condition

% v2' * grad(p)
ridge_function = ...
      V2(:,:,1).*gradS(:,:,1) ...
    + V2(:,:,2).*gradS(:,:,2);

%% Threshold

if nargin < 4 || isempty(threshold)
    threshold = 0;
end

threshold_value = threshold * max(p(:));

active = p >= threshold_value;

%% Ridge points

% Numerical tolerance for the ridge condition.
% We use a relative tolerance with respect to the local gradient.

grad_norm = sqrt( ...
      gradS(:,:,1).^2 ...
    + gradS(:,:,2).^2 );

tol = 0.05 * grad_norm;

ridge = ...
    active ...
    & (lambda2 < 0) ...
    & (abs(ridge_function) <= tol);

end