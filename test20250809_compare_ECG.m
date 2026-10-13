% Generación de la señal
% close all
clear all
% clc


load('ECG.mat');

x = ecg1;

x = x(:)';

N = length(x);

% Cálculo de la STFT y SST
sigma = 10/sqrt(N);%0.05;
fmax = 0.5; %maximum normalized frequency for STFT evaluation
redun = 1; %redundancy for the STFT computation (frequency bins equals N*redun*fmax)
bt = 1:N; %time axis
Nfft = N;
ft = 1/redun:1/redun:round(fmax*Nfft); %frequency axis

% x = x + 0.316*std(x)*randn(1,N);

[Fg, Ftg, Ft2g] = stfts_order_2_time_domain(x,N*sigma,ft,bt-1);%% OJO ACA que es en dominio temporal
% [Fg, Ftg, Ft2g] = stfts_order_2(x,sigma,ft,bt);
gamma = 3*sqrt(2)*median(abs(real(Fg(:))))/.6745;
[SST1,SST2] = my_synchro_order_2(Fg,ft,sigma,Ftg,Ft2g,gamma);

filter = ones(3,3)/9;
absFg_s = conv2(abs(Fg),filter,'same');
delta = sqrt(max(absFg_s(:))*median(absFg_s(:))/2);

S_a = absFg_s>delta;
S_b = abs(Fg)>gamma;
% S = abs(SST2);
% figure; imagesc(abs(S_a));% colormap(gray); colorbar;
% title('MatrixPlot[S]');
% axis square


% Building the graphs
side = 3;
G_a = buildGraph(S_a,side,1);
G_b = buildGraph(S_b,side,1);
G_c = buildGraph(abs(Fg),side,gamma^2);
G_d = buildGraph(abs(SST2)>redun*gamma/sigma,side,1);
G_e = buildGraph(abs(SST2),side,redun*gamma^2/sigma);

% Componentes conexas
comp_a = conncomp(G_a); % devuelve vector con la clase de cada nodo
numComp_a = max(comp_a);

comp_sizes_a = histcounts(comp_a,1:(numComp_a+1));
[~,valid_comps_a] = sort(comp_sizes_a,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_a1 = ComponentMask(S_a, comp_a, valid_comps_a(1));
S_a2 = ComponentMask(S_a, comp_a, valid_comps_a(2));
S_a3 = ComponentMask(S_a, comp_a, valid_comps_a(3));
S_a4 = ComponentMask(S_a, comp_a, valid_comps_a(4));


% Componentes conexas
comp_b = conncomp(G_b); % devuelve vector con la clase de cada nodo
numComp_b = max(comp_b);

comp_sizes_b = histcounts(comp_b,1:(numComp_b+1));
[~,valid_comps_b] = sort(comp_sizes_b,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_b1 = ComponentMask(S_b, comp_b, valid_comps_b(1));
S_b2 = ComponentMask(S_b, comp_b, valid_comps_b(2));
S_b3 = ComponentMask(S_b, comp_b, valid_comps_b(3));
S_b4 = ComponentMask(S_b, comp_b, valid_comps_b(4));


% Componentes conexas
comp_c = conncomp(G_c); % devuelve vector con la clase de cada nodo
numComp_c = max(comp_c);

comp_sizes_c = histcounts(comp_c,1:(numComp_c+1));
[~,valid_comps_c] = sort(comp_sizes_c,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_c1 = ComponentMask(abs(Fg), comp_c, valid_comps_c(1));
S_c2 = ComponentMask(abs(Fg), comp_c, valid_comps_c(2));
S_c3 = ComponentMask(abs(Fg), comp_c, valid_comps_c(3));
S_c4 = ComponentMask(abs(Fg), comp_c, valid_comps_c(4));

% Componentes conexas
comp_d = conncomp(G_d); % devuelve vector con la clase de cada nodo
numComp_d = max(comp_d);

comp_sizes_d = histcounts(comp_d,1:(numComp_d+1));
[~,valid_comps_d] = sort(comp_sizes_d,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_d1 = ComponentMask(abs(Fg), comp_d, valid_comps_d(1));
S_d2 = ComponentMask(abs(Fg), comp_d, valid_comps_d(2));
S_d3 = ComponentMask(abs(Fg), comp_d, valid_comps_d(3));
S_d4 = ComponentMask(abs(Fg), comp_d, valid_comps_d(4));


% Componentes conexas
comp_e = conncomp(G_e); % devuelve vector con la clase de cada nodo
numComp_e = max(comp_e);

comp_sizes_e = histcounts(comp_e,1:(numComp_e+1));
[~,valid_comps_e] = sort(comp_sizes_e,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_e1 = ComponentMask(abs(SST2), comp_e, valid_comps_e(1));
S_e2 = ComponentMask(abs(SST2), comp_e, valid_comps_e(2));
S_e3 = ComponentMask(abs(SST2), comp_e, valid_comps_e(3));
S_e4 = ComponentMask(abs(SST2), comp_e, valid_comps_e(4));




figure;
subplot(241)
imagesc(abs(Fg))
title('STFT')

subplot(242)
imagesc(S_a1 + 2*S_a2 + 3*S_a3 + 4*S_a4)
title('method A')

subplot(243)
imagesc(S_b1 + 2*S_b2 + 3*S_b3 + 4*S_b4)
title('method B')

subplot(244)
imagesc(S_c1 + 2*S_c2 + 3*S_c3 + 4*S_c4)
title('method C')

subplot(245)
imagesc(S_d1 + 2*S_d2 + 3*S_d3 + 4*S_d4)
title('method D')

subplot(246)
imagesc(S_e1 + 2*S_e2 + 3*S_e3 + 4*S_e4)
title('method E')





