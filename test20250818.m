%% Generación de la matriz S
close all
clear all
% clc

N = 512;
t = 0:1/N:1-1/N;


%x = cos(2*pi*(75*t + (10/pi)*sin(2*pi*t))) + [zeros(1,128)  cos(2*pi*(128*(t(129:768/2)) + 64*t(129:768/2).^2))   zeros(1,128)];
% x = cos(2*pi*100*t);

% x = cos(2*pi*(100*t + (25/pi)*sin(2*pi*t))) + [zeros(1,180)  cos(2*pi*(75*(t(181:404)) + 64*t(181:404).^2))   zeros(1,512-404)];
w = [zeros(1,160)  0.5-.5*cos(2*pi*5*(t(1:round(N/10)))) ones(1,140) 0.5+.5*cos(2*pi*5*(t(1:round(N/10)))) zeros(1,514-404)];
x = cos(2*pi*(100*t + (25/pi)*sin(2*pi*t))) + w.*(cos(2*pi*(75*t + 64*t.^2)));
xn = x + 0.316*std(x)*randn(1,N);

% para funcion de hermite---------------------------
% f0 = 0.25; % Location in normalized frequency
% t0 = N*0.5; % Location in time
% t = (0:N-1)-t0;
% order = 45;%tenía 20
% x = hermite_fun(N, order, t, sqrt(N)).*cos(2*pi*f0*t);  
%------------------------------------

sigma = 1/sqrt(N);%0.05;
fmax = 0.5; %maximum normalized frequency for STFT evaluation
redun = 1; %redundancy for the STFT computation (frequency bins equals N*redun*fmax)
bt = 1:N; %time axis
Nfft = N;
ft = 1/redun:1/redun:round(fmax*Nfft); %frequency axis

% x = x + 0.3162*std(x)*randn(1,N);

[Fg, Ftg, Ft2g] = stfts_order_2(xn,sigma,ft,bt);
gamma = 3*sqrt(2)*median(abs(real(Fg(:))))/.6745
delta = sqrt(max(abs(Fg(:))*median(abs(Fg(:))))/2)
[SST1,SST2] = my_synchro_order_2(Fg,ft,sigma,Ftg,Ft2g,gamma);

figure; imagesc(abs(SST1))

figure; imagesc(abs(SST2))


% S = 2*abs(SST2)*sigma;
% S = (abs(SST2)>0);
% S = abs(Fg)>gamma;
S = abs(SST2);
figure; imagesc(abs(S));% colormap(gray); colorbar;
title('MatrixPlot[S]');
axis square

%% Funciones de indexación
ind = @(a,b) size(S,2)*(a-1) + b;

revind = @(k) ...
    (mod(k,size(S,2))==0) .* [k/size(S,2), size(S,2)] + ...
    (mod(k,size(S,2))~=0) .* [(k - mod(k,size(S,2)))/size(S,2) + 1, mod(k,size(S,2))];

%% Dimensiones
box1 = size(S,1); % filas
box2 = size(S,2); % columnas

%% Construcción de matriz dispersa M (grafo de adyacencias)
I = [];
J = [];
lado = 3;
for i = 1:box1
    for j = 1:box2
        for k1 = -lado:lado
            for k2 = -lado:lado
                i2 = i + k1;
                j2 = j + k2;
                if i2 >= 1 && i2 <= box1 && j2 >= 1 && j2 <= box2
                    if S(i,j) * S(i2,j2) >= gamma^2
                        I(end+1) = ind(i,j); %#ok<SAGROW>
                        J(end+1) = ind(i2,j2); %#ok<SAGROW>
                    end
                end
            end
        end
    end
end

M = sparse(I,J,1, box1*box2, box1*box2);

%% Grafo
G = graph(M);

figure; plot(G);
title('AdjacencyGraph[M]');

%% Componentes conexas
comp = conncomp(G); % devuelve vector con la clase de cada nodo
numComp = max(comp);


%% Dibujar máscaras de componentes
% Seleccionar solo componentes con más de 1 nodo
comp_sizes = histcounts(comp,1:(numComp+1));
% valid_comps = find(comp_sizes > 1);
[~,valid_comps] = sort(comp_sizes,'descend');%nos quedamos con los componentes con más nodos en orden decreciente

for idx = 1:min(3,length(valid_comps)) % como en el código, 3 primeras
    comp_id = valid_comps(idx);
    nodes = find(comp==comp_id);
    
    % Convertir índices lineales a (fila,columna)
    coords = arrayfun(@(x) revind(x), nodes, 'UniformOutput', false);
    coords = cell2mat(coords'); % cada fila: [fila,columna]
    
    S2 = zeros(box1,box2);
    for c = 1:size(coords,1)
        S2(coords(c,1), coords(c,2)) = 1;
    end
    
    nodesComp = find(comp==comp_id);

    figure; plot(subgraph(G,nodesComp));
    title(['Subgraph: Componente conexa ' num2str(comp_id)]);

    
    figure; imagesc(S2); colormap(gray); colorbar;
    title(['Componente conexa ' num2str(comp_id)]);
end


absFg = abs(Fg);
delta = sqrt(max(absFg(:)*median(absFg(:))) /2);

figure; subplot(121);
imagesc(abs(Fg)>gamma)
xr_1 = 2*sigma*real(sum(Fg.*(absFg>gamma)));
text(25,25,num2str(20*log10(norm(x)/norm(x-xr_1))),'color',[1 1 1],'fontsize',20)
subplot(122);
imagesc(abs(Fg)>delta)
xr_2 = 2*sigma*real(sum(Fg.*(absFg>delta)));
text(25,25,num2str(20*log10(norm(x)/norm(x-xr_2))),'color',[1 1 1],'fontsize',20)
