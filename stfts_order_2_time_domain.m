function [Fg,Ftg,Ft2g] = stfts_order_2_time_domain(s,sigma,ft,bt)
% stfts : computes 3 STFTs of a signal with different windows in time domain.
%
% INPUTS:
%   s: real or complex signal
%   sigma: window parameter
%   ft: frequency bins
%   bt: time bins
%
% OUTPUTS:
%   Fg: the short-time Fourier transform (STFT) with a Gaussian window g(t)
%   = exp(-pi*t^2/sigma^2)/sigma
%   Ftng: STFT with t^ng.



n = length(s);
%s = s(:);

% Optional parameters
% if nargin<5
%     ft = 1:n/2;
%     bt = 1:n;
% end
nb = length(bt);
neta = length(ft);

% S = fftshift(fft(s));
% f = -n/2:n/2-1;

t = 0:n-1;
f = 0:1/n:1-1/n;


% Initialization
Fg = zeros(neta,nb);
Ftg = zeros(neta,nb);
Ft2g = zeros(neta,nb);


for b = 1:length(bt) % For each column
    
    % Computes STFTs
    % para t^n g
    g = exp(-pi*(t-(b-1)).^2/sigma^2)/sigma;
    tmp = exp(1i*2*pi*f*(b-1)).*fft(s.*g);
    Fg(:,b) = tmp(ft+1);
    tmp = exp(1i*2*pi*f*(b-1)).*fft(s.*((t-b)).*g);
    Ftg(:,b) = tmp(ft+1)/n;
    tmp = exp(1i*2*pi*f*(b-1)).*fft(s.*((t-b).^2).*g);
    Ft2g(:,b) = tmp(ft+1)/n^2;
    
end

end






