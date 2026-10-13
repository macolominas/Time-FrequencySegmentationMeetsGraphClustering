% Generación de la señal
% close all
clear all
% clc

load('EEG.mat');

x = EEG(:)';

N = length(x);

fs = 100;
t = 0:1/fs:N/fs-1/fs;

tam_fuente = 16;

% Cálculo de la STFT y SST
sigma = 1/sqrt(N);%0.05;
fmax = 0.25; %maximum normalized frequency for STFT evaluation
redun = 1; %redundancy for the STFT computation (frequency bins equals N*redun*fmax)
bt = 1:N; %time axis
Nfft = N;
ft = 1/redun:1/redun:round(fmax*Nfft); %frequency axis

% x = x + 0.5623*std(x)*randn(1,N);

[Fg, Ftg, Ft2g] = stfts_order_2_time_domain(x,N*sigma,ft,bt-1);%% OJO ACA que es en dominio temporal
% [Fg, Ftg, Ft2g] = stfts_order_2(x,sigma,ft,bt);
% gamma = 3*sqrt(2)*median(abs(real(Fg(:))))/.6745;
gamma = 3*median(abs(real(Fg(:))))/.6745;
[SST1,SST2] = my_synchro_order_2(Fg,ft,sigma,Ftg,Ft2g,gamma);

% filter = ones(3,3)/9;
% absFg_s = conv2(abs(Fg),filter,'same');
% delta = sqrt(max(absFg_s(:))*median(absFg_s(:))/2);
% 
% S_a = absFg_s>delta;
% S_b = abs(Fg)>gamma;
% S = abs(SST2);
% figure; imagesc(abs(S_a));% colormap(gray); colorbar;
% title('MatrixPlot[S]');
% axis square


% Building the graphs
side = 3;
Fg_aux = zeros(size(Fg));
Fg_aux(12*N/fs:16*N/fs,:) = Fg(12*N/fs:16*N/fs,:);
G_c = buildGraph(abs(Fg_aux),side,gamma^2);


% Componentes conexas
comp_c = conncomp(G_c); % devuelve vector con la clase de cada nodo
numComp_c = max(comp_c);

comp_sizes_c = histcounts(comp_c,1:(numComp_c+1));
[~,valid_comps_c] = sort(comp_sizes_c,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_c1 = ComponentMask(abs(Fg), comp_c, valid_comps_c(1));
S_c2 = ComponentMask(abs(Fg), comp_c, valid_comps_c(2));
S_c3 = ComponentMask(abs(Fg), comp_c, valid_comps_c(3));
S_c4 = ComponentMask(abs(Fg), comp_c, valid_comps_c(4));
S_c5 = ComponentMask(abs(Fg), comp_c, valid_comps_c(5));
S_c6 = ComponentMask(abs(Fg), comp_c, valid_comps_c(6));
S_c7 = ComponentMask(abs(Fg), comp_c, valid_comps_c(7));
S_c8 = ComponentMask(abs(Fg), comp_c, valid_comps_c(8));
S_c9 = ComponentMask(abs(Fg), comp_c, valid_comps_c(9));


load('experto1_comienzo.mat');
load('experto1_duracion.mat');

experto1_comienzo = experto1_comienzo - 520;

figure;
subplot(811)
plot(t,x,'k'); hold on
for n = 1:5
    plot(experto1_comienzo(n):1/fs:experto1_comienzo(n)+experto1_duracion(n),x(fs*experto1_comienzo(n):fs*(experto1_comienzo(n)+experto1_duracion(n))),'r')
end;
title('Input EEG Signal')
set(gca,'xticklabel',[])
set(gca,'fontsize',tam_fuente)

subplot(8,1,[2 4])
imagesc(t,ft*fs/N,abs(Fg))
set(gca,'ydir','normal')
hold on
plot(t,11*ones(size(t)),'m--')
plot(t,16*ones(size(t)),'m--')
title('STFT')
ylabel('frequency (Hz)')
set(gca,'xticklabel',[])
set(gca,'fontsize',tam_fuente)

subplot(8,1,[5 7])
imagesc(t,ft*fs/N,S_c1 + 2*S_c2 + 4*S_c4 + 5*S_c5 + 6*S_c6 + 9*S_c9)
set(gca,'ydir','normal')
ylabel('frequency (Hz)')
title('Detected Masks')
set(gca,'xticklabel',[])
set(gca,'fontsize',tam_fuente)

subplot(8,1,8)
x1 = real(sum(Fg.*S_c1)*0.5/N);
plot(t(x1~=0),x1(x1~=0),'r')
hold on

x2 = real(sum(Fg.*S_c2)*0.5/N);
plot(t(x2~=0),x2(x2~=0),'r')

x4 = real(sum(Fg.*S_c4)*0.5/N);
plot(t(x4~=0),x4(x4~=0),'r')

x5 = real(sum(Fg.*S_c5)*0.5/N);
plot(t(x5~=0),x5(x5~=0),'r')

x6 = real(sum(Fg.*S_c6)*0.5/N);
plot(t(x6~=0),x6(x6~=0),'r')

x9 = real(sum(Fg.*S_c9)*0.5/N);
plot(t(x9~=0),x9(x9~=0),'r')

xlim([0 30])
xlabel('time (s)')
title('Isolated Spindles')
set(gca,'fontsize',tam_fuente)


