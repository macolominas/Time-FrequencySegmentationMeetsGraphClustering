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

% phip = 225 - 500*(t-0.5).^2;
% figure; plot(t,phip)
% 
x1 = zeros(size(t));

for n = 50:150
    x1 = x1 + sin(2*pi*n*(t-0.5))/25;
end;

% 
% x1 = zeros(size(t));
% x1(75) = 10;
% 
% x2 = [zeros(1,150) cos(2*pi*(125*t(151:end) + 75*exp(-5*t(151:end))))];
% 
% x3 = [zeros(1,150) cos(2*pi*(225*t(151:end) + 100*exp(-5*t(151:end))))];
% 

x2 = cos(2*pi*(225*t - (500/3)*(t-0.5).^3));

x = x1 + x2;% + x3;


SNR = 20;
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
G_c = buildGraph(abs(Fg),side,gamma^2);
G_d = buildGraph(abs(SST2)>gamma^2/sigma,side,1);
G_e = buildGraph(abs(SST2),side,gamma^2/sigma);


% Componentes conexas
comp_a = conncomp(G_a); % devuelve vector con la clase de cada nodo
numComp_a = max(comp_a);

comp_sizes_a = histcounts(comp_a,1:(numComp_a+1));
[~,valid_comps_a] = sort(comp_sizes_a,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_a1 = ComponentMask(S_a, comp_a, valid_comps_a(1));
S_a2 = ComponentMask(S_a, comp_a, valid_comps_a(2));
%S_a3 = ComponentMask(S_a, comp_a, valid_comps_a(3));

% Componentes conexas
comp_b = conncomp(G_b); % devuelve vector con la clase de cada nodo
numComp_b = max(comp_b);

comp_sizes_b = histcounts(comp_b,1:(numComp_b+1));
[~,valid_comps_b] = sort(comp_sizes_b,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_b1 = ComponentMask(S_b, comp_b, valid_comps_b(1));
S_b2 = ComponentMask(S_b, comp_b, valid_comps_b(2));
%S_b3 = ComponentMask(S_b, comp_b, valid_comps_b(3));

% Componentes conexas
comp_c = conncomp(G_c); % devuelve vector con la clase de cada nodo
numComp_c = max(comp_c);

comp_sizes_c = histcounts(comp_c,1:(numComp_c+1));
[~,valid_comps_c] = sort(comp_sizes_c,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_c1 = ComponentMask(abs(Fg), comp_c, valid_comps_c(1));
S_c2 = ComponentMask(abs(Fg), comp_c, valid_comps_c(2));
%S_c3 = ComponentMask(abs(Fg), comp_c, valid_comps_c(3));

% Componentes conexas
comp_d = conncomp(G_d); % devuelve vector con la clase de cada nodo
numComp_d = max(comp_d);

comp_sizes_d = histcounts(comp_d,1:(numComp_d+1));
[~,valid_comps_d] = sort(comp_sizes_d,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_d1 = ComponentMask(abs(Fg), comp_d, valid_comps_d(1));
S_d2 = ComponentMask(abs(Fg), comp_d, valid_comps_d(2));
%S_d3 = ComponentMask(abs(Fg), comp_d, valid_comps_d(3));

% Componentes conexas
comp_e = conncomp(G_e); % devuelve vector con la clase de cada nodo
numComp_e = max(comp_e);

comp_sizes_e = histcounts(comp_e,1:(numComp_e+1));
[~,valid_comps_e] = sort(comp_sizes_e,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_e1 = ComponentMask(abs(SST2), comp_e, valid_comps_e(1));
S_e2 = ComponentMask(abs(SST2), comp_e, valid_comps_e(2));
%S_e3 = ComponentMask(abs(SST2), comp_e, valid_comps_e(3));



figure;

subplot(441)
imagesc(abs(Fg))
title('STFT')
set(gca,'ydir','normal')

subplot(442)
imagesc(S_c1 + 2*S_c2);% + 3*S_c3)
title('method A')
% set(gca,'colormap',parula(4))
set(gca,'colormap',jet)
cmap = get(gca,'colormap');
set(gca,'ydir','normal')

subplot(443)
plot(real(sum(Fg.*S_c1))-25,'color',cmap(21,:));
hold on
plot(real(sum(Fg.*S_c2))-50,'color',cmap(42,:));
% plot(real(sum(Fg.*S_c3)),'color',cmap(64,:));
xlim([1 512])

figure;

subplot(241)
imagesc(abs(Fg))
title('STFT')
set(gca,'ydir','normal')
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
xlabel('time')
ylabel('frequency')
set(gca,'fontsize',15)

subplot(242)
imagesc(S_c1 + 2*S_c2);% + 3*S_c3)
title('Method A')
set(gca,'colormap',jet)
cmap = get(gca,'colormap');
% cmap = colormap(jet(5));
set(gca,'ydir','normal')
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
xlabel('time')
ylabel('frequency')
set(gca,'fontsize',15)

subplot(243)
plot(real(sum(Fg.*S_c1))-25,'color',cmap(21,:));
hold on
plot(real(sum(Fg.*S_c2))-50,'color',cmap(42,:));
% plot(real(sum(Fg.*S_c3)),'color',cmap(64,:));
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
xlabel('time')
title('Extracted Modes')
ylim([-65 115])
xlim([1 N])
set(gca,'fontsize',15)









