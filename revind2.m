function coords = revind2(k,S)
% REVIND Convierte índices lineales a coordenadas (fila,columna).
%    coords = revind(k,S)
%
% Entrada:
%   k : escalar, vector o matriz de índices lineales (>=1)
%   S : matriz de referencia para dimensiones
%
% Salida:
%   coords : matriz Nx2, cada fila es [fila,columna]

nCols = size(S,2);

row = floor((k-1)./nCols) + 1;
col = mod(k-1,nCols) + 1;

coords = [row(:), col(:)];

end

% Ejemplo
% S = zeros(50,100);
% coords = revind([104,205,306],S);
% coords =
%     1     4
%     2     5
%     3     6
