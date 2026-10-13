clear all;
close all;

%% Create a signal
rng(0)
N = 2^8;
f0 = 0.25; % Location in normalized frequency
t0 = N*0.5; % Location in time
t = (0:N-1)-t0;
order = 20;
x = hermite_fun(N, order, t, sqrt(2*N)).*cos(2*pi*f0*t);         

%% Contaminate signal with white noise
SNRin = 120;
noise = randn(size(x));
s = sigmerge(x.',noise.', SNRin);

%% Compute Spectrogram
N = length(s); % Compute L for Gaussian window
Nfft = 2*N ;   % number of FFT bins
fmax = 0.5;
prec = 10^(-6);
[w, T] = roundgauss(Nfft, prec);
F = tfrstft(s,1:N,Nfft,w,0);
F = F(1:round(fmax*Nfft),:);
S = abs(F).^2;
S = flipud(S);

%%
figure()
imshow(S)

