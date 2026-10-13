% Generación de la señal
% close all
clear all
clc

N = 512;
t = 0:1/N:1-1/N;

SNR = 5:5:35;
amp = 10.^(-SNR/20);

I = 100;

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
error_1_rd = zeros(length(SNR),I);
error_2_rd = zeros(length(SNR),I);
error_1_z = zeros(length(SNR),I);
error_2_z = zeros(length(SNR),I);
% error_1_f = zeros(length(SNR),I);
% error_2_f = zeros(length(SNR),I);
% error_1_g = zeros(length(SNR),I);
% error_2_g = zeros(length(SNR),I);


x1 = zeros(size(t));

for n = 50:150
    x1 = x1 + sin(2*pi*n*(t-0.5))/25;
end;


x2 = cos(2*pi*(225*t - (500/3)*(t-0.5).^3));

x = x1 + x2;





% para ridge detection y ribbon---------------
jump = 4;
d = 27;
xr_1rd = zeros(1,N);
xr_2rd = zeros(1,N);
% %--------------------------------------------

% para zeros----------------------------------
LB = 2;



%---------------------------------------------


% Para el cálculo de la STFT y SST
sigma = 1/sqrt(N);%0.05;
fmax = 0.5; %maximum normalized frequency for STFT evaluation
redun = 1; %redundancy for the STFT computation (frequency bins equals N*redun*fmax)
bt = 1:N; %time axis
Nfft = N;
ft = 1/redun:1/redun:round(fmax*Nfft); %frequency axis
%----------------------------------

side = 1; %Para el graph
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

%-----ZEROS---------------------------
cerosFg = find_spectrogram_zeros(abs(Fg));
TRI = delaunay(cerosFg);
[~,~,TRI_LADOS] = describe_triangles(TRI,cerosFg,N,1/sigma);
longTriangulos=zeros(size(TRI,1),1);
for k =1:size(TRI,1)
    if any(TRI_LADOS(k,:)>1.7)        
        longTriangulos(k)=1;
    end
end
TRIselected=TRI(logical(longTriangulos),:);
mask = mask_from_triangles(Fg,TRIselected,cerosFg);
G_z = buildGraph(abs(mask),side,1);

comp_z = conncomp(G_z); % devuelve vector con la clase de cada nodo
numComp_z = max(comp_z);

comp_sizes_z = histcounts(comp_z,1:(numComp_z+1));
[~,valid_comps_z] = sort(comp_sizes_z,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

S_z1 = ComponentMask(mask, comp_z, valid_comps_z(1));
S_z2 = ComponentMask(mask, comp_z, valid_comps_z(2));

xr_1z = 2*sigma*real(sum(Fg.*S_z1));
xr_2z = 2*sigma*real(sum(Fg.*S_z2));

%-----------------------------------------


%------ classical RD --------------
Faux = Fg;
c1 = exridge(abs(Faux).*(abs(Faux)>gamma),0,0,jump);
for n = 1:N
    xr_1rd(n) = sum(Faux(max(1,c1(n)-d):min(c1(n)+d,fmax*Nfft),n));
    Faux(max(1,c1(n)-d):min(c1(n)+d,fmax*Nfft),n) = 0;
end;
xr_1rd = 2*sigma*real(xr_1rd);
c2 = exridge(abs(Faux).*(abs(Faux)>gamma),0,0,jump);
for n = 1:N
    xr_2rd(n) = sum(Faux(max(1,c2(n)-d):min(c2(n)+d,fmax*Nfft),n));
%     Faux(max(1,c2(n)-d):min(c2(n)+d,fmax*Nfft),n) = 0;
end;
xr_2rd = 2*sigma*real(xr_2rd);
%-----------------------------------


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

aux1 = norm(x1 - xr_1rd);
aux2 = norm(x1 - xr_2rd);
error_1_rd(j,i) = min(aux1,aux2);

aux1 = norm(x2 - xr_1rd);
aux2 = norm(x2 - xr_2rd);
error_2_rd(j,i) = min(aux1,aux2);

aux1 = norm(x1 - xr_1z);
aux2 = norm(x1 - xr_2z);
error_1_z(j,i) = min(aux1,aux2);

aux1 = norm(x2 - xr_1z);
aux2 = norm(x2 - xr_2z);
error_2_z(j,i) = min(aux1,aux2);

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
plot(SNR,mean(error_1_a,2),'-o')
hold on;
plot(SNR,mean(error_1_b,2),'-d')
plot(SNR,mean(error_1_c,2),'-s')
plot(SNR,mean(error_1_d,2),'-^')
plot(SNR,mean(error_1_e,2),'-<')
plot(SNR,mean(error_1_rd,2),'->')
plot(SNR,mean(error_1_z,2),'-x')
title('Mean Error for sinusoidal chirp')
legend('A','B','C','D','E','location','northeast')
xlabel('SNR in')
subplot(443)
plot(SNR,mean(error_2_a,2),'-o')
hold on;
plot(SNR,mean(error_2_b,2),'-d')
plot(SNR,mean(error_2_c,2),'-s')
plot(SNR,mean(error_2_d,2),'-^')
plot(SNR,mean(error_2_e,2),'-<')
plot(SNR,mean(error_2_rd,2),'->')
plot(SNR,mean(error_2_z,2),'-x')
title('Mean Error for linear chirp')
%legend('A','B','C','D','E','F','G')
xlabel('SNR in')
subplot(444)
plot(SNR,median(SNR_a,2),'-o')
hold on;
plot(SNR,median(SNR_b,2),'-d')
plot(SNR,median(SNR_c,2),'-s')
plot(SNR,median(SNR_d,2),'-^')
plot(SNR,median(SNR_e,2),'-<')
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





% figure;
% subplot(4,12,1:3)
% imagesc(abs(Fg))
% set(gca,'ydir','normal')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% title('STFT')
% xlabel('time')
% ylabel('frequency')
% set(gca,'fontsize',15)
% 
% subplot(4,12,4:6)
% plot(SNR,median(error_1_a,2),'-o')
% hold on;
% plot(SNR,median(error_1_b,2),'-d')
% plot(SNR,median(error_1_c,2),'-s')
% plot(SNR,median(error_1_d,2),'-^')
% plot(SNR,median(error_1_e,2),'-<')
% plot(SNR,median(error_1_rd,2),'->')
% plot(SNR,median(error_1_z,2),'-x')
% title('Median Error for sinusoidal chirp')
% legend('A','B','C','D','E','R','Z','location','northeast','orientation','horizontal','NumColumns',2)
% xlim([5 35])
% xlabel('SNR in')
% set(gca,'fontsize',15)
% 
% subplot(4,12,7:9)
% plot(SNR,median(error_2_a,2),'-o')
% hold on;
% plot(SNR,median(error_2_b,2),'-d')
% plot(SNR,median(error_2_c,2),'-s')
% plot(SNR,median(error_2_d,2),'-^')
% plot(SNR,median(error_2_e,2),'-<')
% plot(SNR,median(error_2_rd,2),'->')
% plot(SNR,median(error_2_z,2),'-x')
% title('Median Error for linear chirp')
% xlim([5 35])
% %legend('A','B','C','D','E','F','G')
% xlabel('SNR in')
% set(gca,'fontsize',15)
% 
% subplot(4,12,13:15)
% imagesc(abs(SST2))
% set(gca,'ydir','normal')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% title('SST2')
% xlabel('time')
% ylabel('frequency')
% set(gca,'fontsize',15)
% 
% subplot(4,12,16)
% boxplot([error_1_a(1,:)' error_1_b(1,:)' error_1_e(1,:)' error_1_rd(1,:)'])
% set(gca,'xticklabel',{'A','B','E','R'})
% set(gca,'yticklabel',[])
% title('5 dB')
% set(gca,'fontsize',15)
% 
% subplot(4,12,17)
% boxplot([error_1_a(4,:)' error_1_b(4,:)' error_1_e(4,:)' error_1_rd(4,:)' error_1_z(4,:)'])
% set(gca,'xticklabel',{'A','B','E','R','Z'})
% set(gca,'yticklabel',[])
% title('20 dB')
% set(gca,'fontsize',15)
% 
% subplot(4,12,18)
% boxplot([error_1_a(7,:)' error_1_b(7,:)' error_1_e(7,:)' error_1_rd(7,:)'])
% set(gca,'xticklabel',{'A','B','E','R'})
% set(gca,'yticklabel',[])
% title('35 dB')
% set(gca,'fontsize',15)
% 
% subplot(4,12,19)
% boxplot([error_2_a(1,:)' error_2_b(1,:)' error_2_e(1,:)' error_2_rd(1,:)'])
% set(gca,'xticklabel',{'A','B','E','R'})
% set(gca,'yticklabel',[])
% title('5 dB')
% set(gca,'fontsize',15)
% 
% subplot(4,12,20)
% boxplot([error_2_a(4,:)' error_2_b(4,:)' error_2_e(4,:)' error_2_rd(4,:)' error_2_z(4,:)'])
% set(gca,'xticklabel',{'A','B','E','R','Z'})
% set(gca,'yticklabel',[])
% title('20 dB')
% set(gca,'fontsize',15)
% 
% subplot(4,12,21)
% boxplot([error_2_a(7,:)' error_2_b(7,:)' error_2_e(7,:)' error_2_rd(7,:)'])
% set(gca,'xticklabel',{'A','B','E','R'})
% set(gca,'yticklabel',[])
% title('35 dB')
% set(gca,'fontsize',15)
% 




figure;
subplot(4,12,1:3)
imagesc(abs(Fg))
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('STFT')
xlabel('time')
ylabel('frequency')
set(gca,'fontsize',15)

subplot(4,12,4:6)
plot(SNR,median(error_1_a,2),'-o')
hold on;
plot(SNR,median(error_1_b,2),'-d')
plot(SNR,median(error_1_c,2),'-s')
plot(SNR,median(error_1_d,2),'-^')
plot(SNR,median(error_1_e,2),'-<')
plot(SNR,median(error_1_rd,2),'->')
plot(SNR,median(error_1_z,2),'-x')
title('Median Error for linear chirp')
legend('A','B','C','D','E','R','Z','location','northeast','orientation','horizontal','NumColumns',2)
xlim([5 35])
xlabel('SNR in')
set(gca,'fontsize',15)

subplot(4,12,7:9)
plot(SNR,median(error_2_a,2),'-o')
hold on;
plot(SNR,median(error_2_b,2),'-d')
plot(SNR,median(error_2_c,2),'-s')
plot(SNR,median(error_2_d,2),'-^')
plot(SNR,median(error_2_e,2),'-<')
plot(SNR,median(error_2_rd,2),'->')
plot(SNR,median(error_2_z,2),'-x')
title('Median Error for Hermite function')
xlim([5 35])
%legend('A','B','C','D','E','F','G')
xlabel('SNR in')
set(gca,'fontsize',15)

subplot(4,12,13:15)
imagesc(abs(SST2))
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('SST2')
xlabel('time')
ylabel('frequency')
set(gca,'fontsize',15)

subplot(4,12,16)
boxplot([error_1_a(1,:)' error_1_b(1,:)' error_1_e(1,:)' error_1_rd(1,:)'])
set(gca,'xticklabel',{'A','B','E','R'})
set(gca,'yticklabel',[])
title('5 dB')
set(gca,'fontsize',15)

subplot(4,12,17)
boxplot([error_1_a(4,:)' error_1_b(4,:)' error_1_e(4,:)' error_1_rd(4,:)' error_1_z(4,:)'])
set(gca,'xticklabel',{'A','B','E','R','Z'})
set(gca,'yticklabel',[])
title('20 dB')
set(gca,'fontsize',15)

subplot(4,12,18)
boxplot([error_1_a(7,:)' error_1_b(7,:)' error_1_e(7,:)' error_1_rd(7,:)'])
set(gca,'xticklabel',{'A','B','E','R'})
set(gca,'yticklabel',[])
title('35 dB')
set(gca,'fontsize',15)

subplot(4,12,19)
boxplot([error_2_a(1,:)' error_2_b(1,:)' error_2_e(1,:)' error_2_rd(1,:)' error_2_z(1,:)'])
set(gca,'xticklabel',{'A','B','E','R','Z'})
set(gca,'yticklabel',[])
title('5 dB')
set(gca,'fontsize',15)

subplot(4,12,20)
boxplot([error_2_a(4,:)' error_2_b(4,:)' error_2_e(4,:)' error_2_z(4,:)'])
set(gca,'xticklabel',{'A','B','E','Z'})
set(gca,'yticklabel',[])
title('20 dB')
set(gca,'fontsize',15)

subplot(4,12,21)
boxplot([error_2_a(7,:)' error_2_b(7,:)' error_2_e(7,:)'])
set(gca,'xticklabel',{'A','B','E'})
set(gca,'yticklabel',[])
title('40 dB')
set(gca,'fontsize',15)

























figure;
subplot(3,9,1:3)
imagesc(abs(Fg))
axis equal
axis tight
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('STFT')
xlabel('time')
ylabel('frequency')
set(gca,'fontsize',15)

subplot(3,9,4:6)
plot(SNR,median(error_1_a,2),'-o')
hold on;
plot(SNR,mean(error_1_b,2),'-d')
plot(SNR,mean(error_1_c,2),'-s')
plot(SNR,mean(error_1_d,2),'-^')
plot(SNR,mean(error_1_e,2),'-<')
plot(SNR,mean(error_1_rd,2),'->')
plot(SNR,mean(error_1_z,2),'-x')
title('Median Error for linear chirp')
legend('A','B','C','D','E','R','Z','location','northeast','orientation','horizontal','NumColumns',2)
xlim([5 35])
xlabel('SNR in')
set(gca,'fontsize',15)

subplot(3,9,7:9)
plot(SNR,mean(error_2_a,2),'-o')
hold on;
plot(SNR,mean(error_2_b,2),'-d')
plot(SNR,mean(error_2_c,2),'-s')
plot(SNR,mean(error_2_d,2),'-^')
plot(SNR,mean(error_2_e,2),'-<')
plot(SNR,mean(error_2_rd,2),'->')
plot(SNR,mean(error_2_z,2),'-x')
title('Median Error for Hermite function')
xlim([5 35])
%legend('A','B','C','D','E','F','G')
xlabel('SNR in')
set(gca,'fontsize',15)

subplot(3,9,10:12)
imagesc(abs(SST2))
axis equal
axis tight
set(gca,'ydir','normal')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('SST2')
xlabel('time')
ylabel('frequency')
set(gca,'fontsize',15)

subplot(3,9,13)
boxplot([error_1_a(1,:)' error_1_b(1,:)' error_1_e(1,:)' error_1_rd(1,:)'])
set(gca,'xticklabel',{'A','B','E','R'})
set(gca,'yticklabel',[])
title('5 dB')
set(gca,'fontsize',15)

subplot(3,9,14)
boxplot([error_1_a(4,:)' error_1_b(4,:)' error_1_e(4,:)' error_1_rd(4,:)' error_1_z(4,:)'])
set(gca,'xticklabel',{'A','B','E','R','Z'})
set(gca,'yticklabel',[])
title('20 dB')
set(gca,'fontsize',15)

subplot(3,9,15)
boxplot([error_1_a(7,:)' error_1_b(7,:)' error_1_e(7,:)' error_1_rd(7,:)'])
set(gca,'xticklabel',{'A','B','E','R'})
set(gca,'yticklabel',[])
title('35 dB')
set(gca,'fontsize',15)

subplot(3,9,16)
boxplot([error_2_a(1,:)' error_2_b(1,:)' error_2_e(1,:)' error_2_rd(1,:)' error_2_z(1,:)'])
set(gca,'xticklabel',{'A','B','E','R','Z'})
set(gca,'yticklabel',[])
title('5 dB')
set(gca,'fontsize',15)

subplot(3,9,17)
boxplot([error_2_a(4,:)' error_2_b(4,:)' error_2_e(4,:)' error_2_z(4,:)'])
set(gca,'xticklabel',{'A','B','E','Z'})
set(gca,'yticklabel',[])
title('20 dB')
set(gca,'fontsize',15)

subplot(3,9,18)
boxplot([error_2_a(7,:)' error_2_b(7,:)' error_2_e(7,:)'])
set(gca,'xticklabel',{'A','B','E'})
set(gca,'yticklabel',[])
title('35 dB')
set(gca,'fontsize',15)


subplot(3,9,19:21)
plot(t,x)
hold on
plot(t,x1-5)
plot(t,x2-10)
set(gca,'xtick',[])
set(gca,'ytick',[])
xlabel('time')



