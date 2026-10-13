% Generación de la señal
% close all
clear all
% clc

N = 512;
t = 0:1/N:1-1/N;


%x = cos(2*pi*(75*t + (10/pi)*sin(2*pi*t))) + [zeros(1,128)  cos(2*pi*(128*(t(129:768/2)) + 64*t(129:768/2).^2))   zeros(1,128)];
% x = cos(2*pi*100*t);

% % x = cos(2*pi*(100*t + (25/pi)*sin(2*pi*t))) + [zeros(1,180)  cos(2*pi*(75*(t(181:404)) + 64*t(181:404).^2))   zeros(1,512-404)];
% w = [zeros(1,160)  0.5-.5*cos(2*pi*5*(t(1:round(N/10)))) ones(1,140) 0.5+.5*cos(2*pi*5*(t(1:round(N/10)))) zeros(1,514-404)];
% % x = cos(2*pi*(100*t + (25/pi)*sin(2*pi*t))) + w.*(cos(2*pi*(75*t + 64*t.^2)));
% x1 = cos(2*pi*(100*t + (25/pi)*sin(2*pi*t)));
% x2 = w.*(cos(2*pi*(75*t + 64*t.^2)));
% x = x1 + x2;


w = [0.5-.5*cos(2*pi*5*(t(1:round(N/10)))) ones(1,410) 0.5+.5*cos(2*pi*5*(t(1:round(N/10))))];
x1 = w.*[cos(2*pi*(100*t + 75*t.^2))];
% para funcion de hermite---------------------------
f0 = 0.175; % Location in normalized frequency
t0 = N*0.75; % Location in time
t = (0:N-1)-t0;
order = 20;%tenía 20
x2 = hermite_fun(N, order, t, sqrt(N)).*cos(2*pi*f0*t);
x2 = std(x1)*x2/std(x2);
x = x1 + x2;
%------------------------------------

SNR = 40;
amp = 10^(-SNR/20);
xn = x + amp*std(x)*randn(1,N);


% Cálculo de la STFT y SST
sigma = 1/sqrt(N);%0.05;
fmax = 0.5; %maximum normalized frequency for STFT evaluation
redun = 1; %redundancy for the STFT computation (frequency bins equals N*redun*fmax)
bt = 1:N; %time axis
Nfft = N;
ft = 1/redun:1/redun:round(fmax*Nfft); %frequency axis

% x = x + 0.3162*std(x)*randn(1,N);

[Fg, Ftg, Ft2g] = stfts_order_2(xn,sigma,ft,bt);
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
G_c = buildGraph(abs(Fg),side,gamma/2);
G_d = buildGraph(abs(Fg),side,gamma^2);
G_e = buildGraph(abs(SST2),side,0.5*gamma/sigma);
% G_e = buildGraph(abs(SST2)>0.5*gamma/sigma,side,1);

% Componentes conexas
comp_a = conncomp(G_a); % devuelve vector con la clase de cada nodo
numComp_a = max(comp_a);

comp_sizes_a = histcounts(comp_a,1:(numComp_a+1));
[~,valid_comps_a] = sort(comp_sizes_a,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_a1 = ComponentMask(S_a, comp_a, valid_comps_a(1));
S_a2 = ComponentMask(S_a, comp_a, valid_comps_a(2));

xr_1a = 2*sigma*real(sum(Fg.*S_a1));
xr_2a = 2*sigma*real(sum(Fg.*S_a2));

% Componentes conexas
comp_b = conncomp(G_b); % devuelve vector con la clase de cada nodo
numComp_b = max(comp_b);

comp_sizes_b = histcounts(comp_b,1:(numComp_b+1));
[~,valid_comps_b] = sort(comp_sizes_b,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_b1 = ComponentMask(S_b, comp_b, valid_comps_b(1));
S_b2 = ComponentMask(S_b, comp_b, valid_comps_b(2));

xr_1b = 2*sigma*real(sum(Fg.*S_b1));
xr_2b = 2*sigma*real(sum(Fg.*S_b2));


% Componentes conexas
comp_c = conncomp(G_c); % devuelve vector con la clase de cada nodo
numComp_c = max(comp_c);

comp_sizes_c = histcounts(comp_c,1:(numComp_c+1));
[~,valid_comps_c] = sort(comp_sizes_c,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_c1 = ComponentMask(abs(Fg), comp_c, valid_comps_c(1));
S_c2 = ComponentMask(abs(Fg), comp_c, valid_comps_c(2));

xr_1c = 2*sigma*real(sum(Fg.*S_c1));
xr_2c = 2*sigma*real(sum(Fg.*S_c2));


% Componentes conexas
comp_d = conncomp(G_d); % devuelve vector con la clase de cada nodo
numComp_d = max(comp_d);

comp_sizes_d = histcounts(comp_d,1:(numComp_d+1));
[~,valid_comps_d] = sort(comp_sizes_d,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_d1 = ComponentMask(abs(Fg), comp_d, valid_comps_d(1));
S_d2 = ComponentMask(abs(Fg), comp_d, valid_comps_d(2));

xr_1d = 2*sigma*real(sum(Fg.*S_d1));
xr_2d = 2*sigma*real(sum(Fg.*S_d2));

% Componentes conexas
comp_e = conncomp(G_e); % devuelve vector con la clase de cada nodo
numComp_e = max(comp_e);

comp_sizes_e = histcounts(comp_e,1:(numComp_e+1));
[~,valid_comps_e] = sort(comp_sizes_e,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_e1 = ComponentMask(abs(SST2), comp_e, valid_comps_e(1));
S_e2 = ComponentMask(abs(SST2), comp_e, valid_comps_e(2));

xr_1e = 2*sigma*real(sum(SST2.*S_e1));
xr_2e = 2*sigma*real(sum(SST2.*S_e2));

figure;
subplot(251)
imagesc(S_a1)
text(25,25,num2str(20*log10(norm(x1)/norm(x1 - xr_1a))),'color',[1 1 1],'fontsize',20);
title('mode 1, method A')
text(25,240,['input SNR = ' num2str(SNR)],'color',[1 1 1],'fontsize',20);
subplot(256)
imagesc(S_a2)
text(25,25,num2str(20*log10(norm(x2)/norm(x2 - xr_2a))),'color',[1 1 1],'fontsize',20);
title('mode 2, method A')

subplot(252)
imagesc(S_b1)
text(25,25,num2str(20*log10(norm(x1)/norm(x1 - xr_1b))),'color',[1 1 1],'fontsize',20);
title('mode 1, method B')
subplot(257)
imagesc(S_b2)
text(25,25,num2str(20*log10(norm(x2)/norm(x2 - xr_2b))),'color',[1 1 1],'fontsize',20);
title('mode 2, method B')

subplot(253)
imagesc(S_c1)
text(25,25,num2str(20*log10(norm(x1)/norm(x1 - xr_1c))),'color',[1 1 1],'fontsize',20);
title('mode 1, method C')
subplot(258)
imagesc(S_c2)
text(25,25,num2str(20*log10(norm(x2)/norm(x2 - xr_2c))),'color',[1 1 1],'fontsize',20);
title('mode 2, method C')

subplot(254)
imagesc(S_d1)
text(25,25,num2str(20*log10(norm(x1)/norm(x1 - xr_1d))),'color',[1 1 1],'fontsize',20);
title('mode 1, method D')
subplot(259)
imagesc(S_d2)
text(25,25,num2str(20*log10(norm(x2)/norm(x2 - xr_2d))),'color',[1 1 1],'fontsize',20);
title('mode 2, method D')

subplot(255)
imagesc(S_e1)
text(25,25,num2str(20*log10(norm(x1)/norm(x1 - xr_1e))),'color',[1 1 1],'fontsize',20);
title('mode 1, method E')
subplot(2,5,10)
imagesc(S_e2)
text(25,25,num2str(20*log10(norm(x2)/norm(x2 - xr_2e))),'color',[1 1 1],'fontsize',20);
title('mode 2, method E')


%error y correlación...








