% Generación de la señal
% close all
clear all
% clc

N = 512;
t = 0:1/N:1-1/N;

SNR = 5:10:35;
amp = 10.^(-SNR/20);

I = 5;

SNR_a = zeros(length(SNR),I);
SNR_b = zeros(length(SNR),I);
SNR_c = zeros(length(SNR),I);
SNR_d = zeros(length(SNR),I);
SNR_e = zeros(length(SNR),I);
% SNR_f = zeros(length(SNR),I);
% SNR_g = zeros(length(SNR),I);

error_1_a = zeros(length(SNR),I);
error_2_a = zeros(length(SNR),I);
error_1_b = zeros(length(SNR),I);
error_2_b = zeros(length(SNR),I);
error_1_c = zeros(length(SNR),I);
error_2_c = zeros(length(SNR),I);
error_1_d = zeros(length(SNR),I);
error_2_d = zeros(length(SNR),I);
error_1_e = zeros(length(SNR),I);
error_2_e = zeros(length(SNR),I);
% error_1_f = zeros(length(SNR),I);
% error_2_f = zeros(length(SNR),I);
% error_1_g = zeros(length(SNR),I);
% error_2_g = zeros(length(SNR),I);


%x = cos(2*pi*(75*t + (10/pi)*sin(2*pi*t))) + [zeros(1,128)  cos(2*pi*(128*(t(129:768/2)) + 64*t(129:768/2).^2))   zeros(1,128)];
% x = cos(2*pi*100*t);

% para funcion de hermite---------------------------
% f0 = 0.25; % Location in normalized frequency
% t0 = N*0.5; % Location in time
% t = (0:N-1)-t0;
% order = 45;%tenía 20
% x = hermite_fun(N, order, t, sqrt(N)).*cos(2*pi*f0*t);  
%------------------------------------

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

% % x = cos(2*pi*(100*t + (25/pi)*sin(2*pi*t))) + [zeros(1,180)  cos(2*pi*(75*(t(181:404)) + 64*t(181:404).^2))   zeros(1,512-404)];
% w = [zeros(1,160)  0.5-.5*cos(2*pi*5*(t(1:round(N/10)))) ones(1,140) 0.5+.5*cos(2*pi*5*(t(1:round(N/10)))) zeros(1,514-404)];
% % x = cos(2*pi*(100*t + (25/pi)*sin(2*pi*t))) + w.*(cos(2*pi*(75*t + 64*t.^2)));
% x1 = cos(2*pi*(100*t + (25/pi)*sin(2*pi*t)));
% x2 = w.*(cos(2*pi*(75*t + 64*t.^2))); %75 t es la base de frecuencia
% x = x1 + x2;


% Para el cálculo de la STFT y SST
sigma = 1/sqrt(N);%0.05;
fmax = 0.5; %maximum normalized frequency for STFT evaluation
redun = 1; %redundancy for the STFT computation (frequency bins equals N*redun*fmax)
bt = 1:N; %time axis
Nfft = N;
ft = 1/redun:1/redun:round(fmax*Nfft); %frequency axis
%----------------------------------

side = 3; %Para el graph
filter = ones(3,3)/9;% para filtrar al modulo de la STFT

for i = 1:I
ruido = randn(1,N);
for j = 1:length(SNR)
xn = x + amp(j)*std(x)*ruido;

[Fg, Ftg, Ft2g] = stfts_order_2(xn,sigma,ft,bt);
gamma = 3*sqrt(2)*median(abs(real(Fg(:))))/.6745;
[SST1,SST2] = my_synchro_order_2(Fg,ft,sigma,Ftg,Ft2g,gamma);



absFg_s = conv2(abs(Fg),filter,'same');
delta = sqrt(max(absFg_s(:))*median(absFg_s(:))/2);
S_e = absFg_s>delta;
S_b = abs(Fg)>gamma;

% Building the graphs
G_e = buildGraph(S_e,side,1);
G_b = buildGraph(S_b,side,1);
G_a = buildGraph(abs(Fg),side,gamma^2);
% G_d = buildGraph(abs(Fg),side,gamma^2);
G_d = buildGraph(abs(SST2)>gamma/sigma,side,1);
% G_f = buildGraph(abs(SST2),side,0.5*gamma/sigma);
G_c = buildGraph(abs(SST2),side,gamma^2/sigma^2);


% Componentes conexas
comp_a = conncomp(G_a); % devuelve vector con la clase de cada nodo
numComp_a = max(comp_a);

comp_sizes_a = histcounts(comp_a,1:(numComp_a+1));
[~,valid_comps_a] = sort(comp_sizes_a,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_a1 = ComponentMask(abs(Fg), comp_a, valid_comps_a(1));
S_a2 = ComponentMask(abs(Fg), comp_a, valid_comps_a(2));

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

S_c1 = ComponentMask(abs(SST2), comp_c, valid_comps_c(1));
S_c2 = ComponentMask(abs(SST2), comp_c, valid_comps_c(2));

xr_1c = 2*sigma*real(sum(SST2.*S_c1));
xr_2c = 2*sigma*real(sum(SST2.*S_c2));

% Componentes conexas
comp_d = conncomp(G_d); % devuelve vector con la clase de cada nodo
numComp_d = max(comp_d);

comp_sizes_d = histcounts(comp_d,1:(numComp_d+1));
[~,valid_comps_d] = sort(comp_sizes_d,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_d1 = ComponentMask(abs(SST2), comp_d, valid_comps_d(1));
S_d2 = ComponentMask(abs(SST2), comp_d, valid_comps_d(2));

xr_1d = 2*sigma*real(sum(SST2.*S_d1));
xr_2d = 2*sigma*real(sum(SST2.*S_d2));

% Componentes conexas
comp_e = conncomp(G_e); % devuelve vector con la clase de cada nodo
numComp_e = max(comp_e);

comp_sizes_e = histcounts(comp_e,1:(numComp_e+1));
[~,valid_comps_e] = sort(comp_sizes_e,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_e1 = ComponentMask(S_e, comp_e, valid_comps_e(1));
S_e2 = ComponentMask(S_e, comp_e, valid_comps_e(2));

xr_1e = 2*sigma*real(sum(Fg.*S_e1));
xr_2e = 2*sigma*real(sum(Fg.*S_e2));

% % Componentes conexas
% comp_f = conncomp(G_f); % devuelve vector con la clase de cada nodo
% numComp_f = max(comp_f);
% 
% comp_sizes_f = histcounts(comp_f,1:(numComp_f+1));
% [~,valid_comps_f] = sort(comp_sizes_f,'descend');%nos quedamos con los componentes con más nodos en orden decreciente
% 
% S_f1 = ComponentMask(abs(SST2), comp_f, valid_comps_f(1));
% S_f2 = ComponentMask(abs(SST2), comp_f, valid_comps_f(2));
% 
% xr_1f = 2*sigma*real(sum(SST2.*S_f1));
% xr_2f = 2*sigma*real(sum(SST2.*S_f2));
% 
% % Componentes conexas
% comp_g = conncomp(G_g); % devuelve vector con la clase de cada nodo
% numComp_g = max(comp_g);
% 
% comp_sizes_g = histcounts(comp_g,1:(numComp_g+1));
% [~,valid_comps_g] = sort(comp_sizes_g,'descend');%nos quedamos con los componentes con más nodos en orden decreciente
% 
% S_g1 = ComponentMask(abs(SST2), comp_g, valid_comps_g(1));
% S_g2 = ComponentMask(abs(SST2), comp_g, valid_comps_g(2));
% 
% xr_1g = 2*sigma*real(sum(SST2.*S_g1));
% xr_2g = 2*sigma*real(sum(SST2.*S_g2));

aux1 = norm(x1 - xr_1a);
aux2 = norm(x1 - xr_2a);
error_1_a(j,i) = min(aux1,aux2);

aux1 = norm(x2 - xr_1a);
aux2 = norm(x2 - xr_2a);
error_2_a(j,i) = min(aux1,aux2);

aux1 = norm(x1 - xr_1b);
aux2 = norm(x1 - xr_2b);
error_1_b(j,i) = min(aux1,aux2);

aux1 = norm(x2 - xr_1b);
aux2 = norm(x2 - xr_2b);
error_2_b(j,i) = min(aux1,aux2);

aux1 = norm(x1 - xr_1c);
aux2 = norm(x1 - xr_2c);
error_1_c(j,i) = min(aux1,aux2);

aux1 = norm(x2 - xr_1c);
aux2 = norm(x2 - xr_2c);
error_2_c(j,i) = min(aux1,aux2);

aux1 = norm(x1 - xr_1d);
aux2 = norm(x1 - xr_2d);
error_1_d(j,i) = min(aux1,aux2);

aux1 = norm(x2 - xr_1d);
aux2 = norm(x2 - xr_2d);
error_2_d(j,i) = min(aux1,aux2);

aux1 = norm(x1 - xr_1e);
aux2 = norm(x1 - xr_2e);
error_1_e(j,i) = min(aux1,aux2);

aux1 = norm(x2 - xr_1e);
aux2 = norm(x2 - xr_2e);
error_2_e(j,i) = min(aux1,aux2);

% aux1 = norm(x1 - xr_1f);
% aux2 = norm(x1 - xr_2f);
% error_1_f(j,i) = min(aux1,aux2);
% 
% aux1 = norm(x2 - xr_1f);
% aux2 = norm(x2 - xr_2f);
% error_2_f(j,i) = min(aux1,aux2);
% 
% aux1 = norm(x1 - xr_1g);
% aux2 = norm(x1 - xr_2g);
% error_1_g(j,i) = min(aux1,aux2);
% 
% aux1 = norm(x2 - xr_1g);
% aux2 = norm(x2 - xr_2g);
% error_2_g(j,i) = min(aux1,aux2);

SNR_a(j,i) = 20*log10(norm(x)/norm(x - xr_1a - xr_2a));
SNR_b(j,i) = 20*log10(norm(x)/norm(x - xr_1b - xr_2b));
SNR_c(j,i) = 20*log10(norm(x)/norm(x - xr_1c - xr_2c));
SNR_d(j,i) = 20*log10(norm(x)/norm(x - xr_1d - xr_2d));
SNR_e(j,i) = 20*log10(norm(x)/norm(x - xr_1e - xr_2e));
% SNR_f(j,i) = 20*log10(norm(x)/norm(x - xr_1f - xr_2f));
% SNR_g(j,i) = 20*log10(norm(x)/norm(x - xr_1g - xr_2g));

end;%i
end;%j

SNR_1_a = 20*log10(norm(x1)./error_1_a);
SNR_2_a = 20*log10(norm(x2)./error_2_a);
SNR_1_b = 20*log10(norm(x1)./error_1_b);
SNR_2_b = 20*log10(norm(x2)./error_2_b);
SNR_1_c = 20*log10(norm(x1)./error_1_c);
SNR_2_c = 20*log10(norm(x2)./error_2_c);
SNR_1_d = 20*log10(norm(x1)./error_1_d);
SNR_2_d = 20*log10(norm(x2)./error_2_d);
SNR_1_e = 20*log10(norm(x1)./error_1_e);
SNR_2_e = 20*log10(norm(x2)./error_2_e);
% SNR_1_f = 20*log10(norm(x1)./error_1_f);
% SNR_2_f = 20*log10(norm(x2)./error_2_f);
% SNR_1_g = 20*log10(norm(x1)./error_1_g);
% SNR_2_g = 20*log10(norm(x2)./error_2_g);


figure;
subplot(441)
imagesc(abs(Fg))
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('STFT')
xlabel('time')
ylabel('frequency')
subplot(442)
plot(SNR,median(error_1_a,2),'-o')
hold on;
plot(SNR,median(error_1_b,2),'-d')
plot(SNR,median(error_1_c,2),'-s')
plot(SNR,median(error_1_d,2),'-^')
plot(SNR,median(error_1_e,2),'-<')
% plot(SNR,median(error_1_f,2),'->')
% plot(SNR,median(error_1_g,2),'-v')
title('Median Error for sinusoidal chirp')
legend('A','B','C','D','E','location','northeast')
xlabel('SNR in')
subplot(443)
plot(SNR,median(error_2_a,2),'-o')
hold on;
plot(SNR,median(error_2_b,2),'-d')
plot(SNR,median(error_2_c,2),'-s')
plot(SNR,median(error_2_d,2),'-^')
plot(SNR,median(error_2_e,2),'-<')
% plot(SNR,median(error_2_f,2),'->')
% plot(SNR,median(error_2_g,2),'-V')
title('Median Error for linear chirp')
%legend('A','B','C','D','E','F','G')
xlabel('SNR in')
subplot(444)
plot(SNR,median(SNR_a,2),'-o')
hold on;
plot(SNR,median(SNR_b,2),'-d')
plot(SNR,median(SNR_c,2),'-s')
plot(SNR,median(SNR_d,2),'-^')
plot(SNR,median(SNR_e,2),'-<')
% plot(SNR,median(SNR_f,2),'->')
% plot(SNR,median(SNR_g,2),'-v')
title('Global SNR out')
%legend('A','B','C','D','E','F','G')
xlabel('SNR in')

% subplot(4,4,6)
% boxplot([error_1_a(1,:)' error_1_b(1,:)' error_1_c(1,:)' error_1_d(1,:)' error_1_f(1,:)' error_1_g(1,:)'])
% set(gca,'xticklabel',{'A','B','C','D','F','G'})
% 
% subplot(4,4,10)
% boxplot([error_1_a(4,:)' error_1_b(4,:)' error_1_c(4,:)' error_1_d(4,:)' error_1_f(4,:)' error_1_g(4,:)'])
% set(gca,'xticklabel',{'A','B','C','D','F','G'})
% 
% subplot(4,4,14)
% boxplot([error_1_a(5,:)' error_1_b(8,:)' error_1_c(8,:)' error_1_f(8,:)' error_1_g(8,:)'])
% set(gca,'xticklabel',{'A','B','C','F','G'})
% 
% 
%error y correlación...


figure;
subplot(4,12,1:3)
% subplot(441)
imagesc(abs(Fg))
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('STFT')
xlabel('time')
ylabel('frequency')
set(gca,'fontsize',15)

subplot(4,12,4:6)
% subplot(442)
plot(SNR,median(error_1_a,2),'-o')
hold on;
plot(SNR,median(error_1_b,2),'-d')
plot(SNR,median(error_1_c,2),'-s')
plot(SNR,median(error_1_d,2),'-^')
plot(SNR,median(error_1_e,2),'-<')
% plot(SNR,median(error_1_f,2),'->')
% plot(SNR,median(error_1_g,2),'-v')
title('Median Error for sinusoidal chirp')
legend('A','B','C','D','E','location','northeast')
xlim([5 40])
xlabel('SNR in')
set(gca,'fontsize',15)

subplot(4,12,7:9)
% subplot(443)
plot(SNR,median(error_2_a,2),'-o')
hold on;
plot(SNR,median(error_2_b,2),'-d')
plot(SNR,median(error_2_c,2),'-s')
plot(SNR,median(error_2_d,2),'-^')
plot(SNR,median(error_2_e,2),'-<')
% plot(SNR,median(error_2_f,2),'->')
% plot(SNR,median(error_2_g,2),'-V')
title('Median Error for linear chirp')
xlim([5 40])
%legend('A','B','C','D','E','F','G')
xlabel('SNR in')
set(gca,'fontsize',15)

subplot(4,12,13:15)
% subplot(441)
imagesc(abs(SST2))
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('SST2')
xlabel('time')
ylabel('frequency')
set(gca,'fontsize',15)

subplot(4,12,16)
boxplot([error_1_a(1,:)' error_1_b(1,:)' error_1_e(1,:)'])
set(gca,'xticklabel',{'A','B','E'})
set(gca,'yticklabel',[])
title('5 dB')
set(gca,'fontsize',15)

subplot(4,12,17)
boxplot([error_1_a(4,:)' error_1_b(4,:)' error_1_e(4,:)'])
set(gca,'xticklabel',{'A','B','E'})
set(gca,'yticklabel',[])
title('20 dB')
set(gca,'fontsize',15)

subplot(4,12,18)
boxplot([error_1_a(8,:)' error_1_b(8,:)' error_1_c(8,:)' error_1_d(8,:)' error_1_e(8,:)'])
set(gca,'xticklabel',{'A','B','C','D','E'})
set(gca,'yticklabel',[])
title('40 dB')
set(gca,'fontsize',15)

subplot(4,12,19)
boxplot([error_2_a(1,:)' error_2_b(1,:)' error_2_e(1,:)'])
set(gca,'xticklabel',{'A','B','E'})
set(gca,'yticklabel',[])
title('5 dB')
set(gca,'fontsize',15)

subplot(4,12,20)
boxplot([error_2_a(4,:)' error_2_b(4,:)' error_2_e(4,:)'])
set(gca,'xticklabel',{'A','B','E'})
set(gca,'yticklabel',[])
title('20 dB')
set(gca,'fontsize',15)

subplot(4,12,21)
boxplot([error_2_a(8,:)' error_2_b(8,:)' error_2_c(8,:)' error_2_d(8,:)' error_2_e(8,:)'])
set(gca,'xticklabel',{'A','B','C','D','E'})
set(gca,'yticklabel',[])
title('40 dB')
set(gca,'fontsize',15)






figure;
subplot(4,12,1:3)
% subplot(441)
imagesc(abs(Fg))
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('STFT')
xlabel('time')
ylabel('frequency')
set(gca,'fontsize',15)

subplot(4,12,4:6)
% subplot(442)
plot(SNR,median(error_1_a,2),'-o')
hold on;
plot(SNR,median(error_1_b,2),'-d')
plot(SNR,median(error_1_c,2),'-s')
plot(SNR,median(error_1_d,2),'-^')
plot(SNR,median(error_1_e,2),'-<')
% plot(SNR,median(error_1_f,2),'->')
% plot(SNR,median(error_1_g,2),'-v')
title('Median Error for sinusoidal chirp')
legend('A','B','C','D','E','location','northeast','orientation','horizontal')
xlim([5 40])
xlabel('SNR in')
set(gca,'fontsize',15)

subplot(4,12,7:9)
% subplot(443)
plot(SNR,median(error_2_a,2),'-o')
hold on;
plot(SNR,median(error_2_b,2),'-d')
plot(SNR,median(error_2_c,2),'-s')
plot(SNR,median(error_2_d,2),'-^')
plot(SNR,median(error_2_e,2),'-<')
% plot(SNR,median(error_2_f,2),'->')
% plot(SNR,median(error_2_g,2),'-V')
title('Median Error for linear chirp')
xlim([5 40])
%legend('A','B','C','D','E','F','G')
xlabel('SNR in')
set(gca,'fontsize',15)

subplot(4,12,13:15)
% subplot(441)
imagesc(abs(SST2))
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('SST2')
xlabel('time')
ylabel('frequency')
set(gca,'fontsize',15)

subplot(4,12,16)
boxplot([error_1_a(1,:)' error_1_b(1,:)' error_1_e(1,:)'])
set(gca,'xticklabel',{'A','B','E'})
set(gca,'yticklabel',[])
title('5 dB')
set(gca,'fontsize',15)

subplot(4,12,17)
boxplot([error_1_a(4,:)' error_1_b(4,:)' error_1_e(4,:)'])
set(gca,'xticklabel',{'A','B','E'})
set(gca,'yticklabel',[])
title('20 dB')
set(gca,'fontsize',15)

subplot(4,12,18)
boxplot([error_1_c(8,:)' error_1_d(8,:)'])
set(gca,'xticklabel',{'C','D'})
set(gca,'yticklabel',[])
title('40 dB')
set(gca,'fontsize',15)

subplot(4,12,19)
boxplot([error_2_a(1,:)' error_2_b(1,:)' error_2_e(1,:)'])
set(gca,'xticklabel',{'A','B','E'})
set(gca,'yticklabel',[])
title('5 dB')
set(gca,'fontsize',15)

subplot(4,12,20)
boxplot([error_2_a(4,:)' error_2_b(4,:)' error_2_e(4,:)'])
set(gca,'xticklabel',{'A','B','E'})
set(gca,'yticklabel',[])
title('20 dB')
set(gca,'fontsize',15)

subplot(4,12,21)
boxplot([error_2_c(8,:)' error_2_d(8,:)'])
set(gca,'xticklabel',{'C','D'})
set(gca,'yticklabel',[])
title('40 dB')
set(gca,'fontsize',15)







