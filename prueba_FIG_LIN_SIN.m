% --- Inicialización de la figura con proporciones verticales ideales ---
figure('Units', 'normalized', 'Position', [0.05, 0.05, 0.70, 0.88]);

% =========================================================================
% CONFIGURACIÓN DE LA GEOMETRÍA (Simula una grilla 4x4 usando solo 3 columnas)
% =========================================================================
x_min = 0.05;
x_max = 0.72; 
ancho_util = x_max - x_min; % Espacio neto real ocupado (deja libre el último cuarto)

% Límites de la columna del tercio izquierdo dentro de la zona útil
w_izq = ancho_util * 0.33;  

% Límites del bloque derecho (2/3 restantes de la zona útil)
x_bloque_der = x_min + w_izq + 0.04; % Espacio de separación entre izquierda y derecha
ancho_bloque_der = x_max - x_bloque_der;
w_grafica_error = ancho_bloque_der * 0.47; % Ancho individual de cada gráfico de error
gap_bloques_der = ancho_bloque_der * 0.06;  % Separación entre bloque Chirp y Hermite

% Alturas fijas compartidas
y_sup = 0.54;
h_sup = 0.39;
y_inf = 0.08;
h_inf = 0.38;

tam_fuente = 16; % Fuente grande y legible

% =========================================================================
% 1. COLUMNA DEL TERCIO IZQUIERDO (Tres filas independientes - NADA SE BORRA)
% =========================================================================
% Fila Superior Izquierda: STFT
axes('Position', [x_min, 0.68, w_izq, 0.25]);
imagesc(abs(Fg))
axis equal; axis tight;
set(gca, 'ydir', 'normal', 'xtick', [], 'ytick', [])
title('STFT'); xlabel('time'); ylabel('frequency'); set(gca, 'fontsize', tam_fuente)

% Fila Media Izquierda: SST2
axes('Position', [x_min, 0.38, w_izq, 0.25]);
imagesc(abs(SST2))
axis equal; axis tight;
set(gca, 'ydir', 'normal', 'xtick', [], 'ytick', [])
title('SST2'); xlabel('time'); ylabel('frequency'); set(gca, 'fontsize', tam_fuente)

% Fila Inferior Izquierda: Tres plots desplazados verticalmente
axes('Position', [x_min, 0.08, w_izq, 0.25]);
plot(t, x)
hold on
plot(t, x1-5)
plot(t, x2-10)
set(gca, 'xtick', [], 'ytick', [])
xlabel('time'); set(gca, 'fontsize', tam_fuente)
ylabel('amplitude (a. u.)'); set(gca, 'fontsize', tam_fuente)
text(0.05,1.9,'mixed signal','fontsize', tam_fuente)
text(0.05,1.9-5,'sinusoidale chirp','fontsize', tam_fuente)
text(0.05,1.9-11,'linear chirp','fontsize', tam_fuente)
ylim([-12 3])



% =========================================================================
% DATOS Y CONFIGURACIÓN DE LOS 6 BOXPLOTS (Ajustado a 4, 5, 4, 4, 5, 4 columnas)
% =========================================================================
box_data = { ...
    [error_1_a(1, :)' error_1_b(1, :)' error_1_e(1, :)' error_1_rd(1, :)'], ...          % 4 columnas
    [error_1_a(4, :)' error_1_b(4, :)' error_1_e(4, :)' error_1_rd(4, :)' error_1_z(4, :)'], ... % 5 columnas
    [error_1_a(7, :)' error_1_b(7, :)' error_1_e(7, :)' error_1_rd(7, :)'], ...          % 4 columnas
    [error_2_a(1, :)' error_2_b(1, :)' error_2_e(1, :)' error_2_rd(1, :)'], ...          % 4 columnas (Cambió a 4)
    [error_2_a(4, :)' error_2_b(4, :)' error_2_e(4, :)' error_2_rd(4, :)' error_2_z(4, :)'], ... % 5 columnas (Cambió a 5)
    [error_2_a(7, :)' error_2_b(7, :)' error_2_e(7, :)' error_2_rd(7, :)'] ...           % 4 columnas (Cambió a 4)
};

box_labels = { ...
    {'A', 'B', 'E', 'R'}, ...
    {'A', 'B', 'E', 'R', 'Z'}, ...
    {'A', 'B', 'E', 'R'}, ...
    {'A', 'B', 'E', 'R'}, ...       % 4 etiquetas
    {'A', 'B', 'E', 'R', 'Z'}, ...  % 5 etiquetas
    {'A', 'B', 'E', 'R'} ...        % 4 etiquetas
};

box_titles = {'5 dB', '20 dB', '35 dB', '5 dB', '20 dB', '35 dB'};
columnas_por_box = cellfun(@(x) size(x, 2), box_data);
gap_entre_boxplots = 0.006; 


% =========================================================================
% 2. BLOQUE DERECHO - CHIRP (Gráfica de Error Superior + 3 Boxplots Inferiores)
% =========================================================================
x_start_chirp = x_bloque_der;

% Error Linear Chirp (Arriba)
axes('Position', [x_start_chirp, y_sup, w_grafica_error, h_sup]);
plot(SNR, median(error_1_a, 2), '-o'); hold on;
plot(SNR, median(error_1_b, 2), '-d'); plot(SNR, median(error_1_c, 2), '-s');
plot(SNR, median(error_1_d, 2), '-^'); plot(SNR, median(error_1_e, 2), '-<');
plot(SNR, median(error_1_rd, 2), '->'); plot(SNR, median(error_1_z, 2), '-x');
title('Median Error for sinusoidal chirp')
legend('STFT-MulPol (A)', 'STFT-BinPol (B)', 'SST2-MulPol (C)', 'SST2-BinPol (D)', 'Method E', 'Classic RD (R)', 'Spec. zeros (Z)', 'location', 'northeast', 'orientation', 'vertical', 'NumColumns', 1)
legend boxoff
xlim([5 35]); xlabel('SNR in'); set(gca, 'fontsize', tam_fuente)

% 3 Boxplots de Chirp (Abajo alineados con su gráfica superior)
columnas_totales_chirp = sum(columnas_por_box(1:3));
ancho_neto_chirp = w_grafica_error - (2 * gap_entre_boxplots);
x_box = x_start_chirp;

for i = 1:3
    w_celda = (columnas_por_box(i) / columnas_totales_chirp) * ancho_neto_chirp;
    
    ax_box = axes('Position', [x_box, y_inf, w_celda, h_inf]);
    boxplot(box_data{i});
    
    set(ax_box, 'xticklabel', box_labels{i})
    set(ax_box, 'yticklabel', [])
    title(box_titles{i})
    set(ax_box, 'fontsize', tam_fuente)
    
    x_box = x_box + w_celda + gap_entre_boxplots;
end


% =========================================================================
% 3. BLOQUE DERECHO - HERMITE (Gráfica de Error Superior + 3 Boxplots Inferiores)
% =========================================================================
x_start_hermite = x_start_chirp + w_grafica_error + gap_bloques_der;

% Error Hermite Function (Arriba)
axes('Position', [x_start_hermite, y_sup, w_grafica_error, h_sup]);
plot(SNR, median(error_2_a, 2), '-o'); hold on;
plot(SNR, median(error_2_b, 2), '-d'); plot(SNR, median(error_2_c, 2), '-s');
plot(SNR, median(error_2_d, 2), '-^'); plot(SNR, median(error_2_e, 2), '-<');
plot(SNR, median(error_2_rd, 2), '->'); plot(SNR, median(error_2_z, 2), '-x');
title('Median Error for linear chirp')
xlim([5 35]); xlabel('SNR in'); set(gca, 'fontsize', tam_fuente)

% 3 Boxplots de Hermite (Abajo alineados con su gráfica superior)
columnas_totales_hermite = sum(columnas_por_box(4:6));
ancho_neto_hermite = w_grafica_error - (2 * gap_entre_boxplots);
x_box = x_start_hermite;

for i = 4:6
    w_celda = (columnas_por_box(i) / columnas_totales_hermite) * ancho_neto_hermite;
    
    ax_box = axes('Position', [x_box, y_inf, w_celda, h_inf]);
    boxplot(box_data{i});
    
    set(ax_box, 'xticklabel', box_labels{i})
    set(ax_box, 'yticklabel', [])
    title(box_titles{i})
    set(ax_box, 'fontsize', tam_fuente)
    
    x_box = x_box + w_celda + gap_entre_boxplots;
end
