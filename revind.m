function coords = revind(k,S)
% Convierte índice lineal k a coordenadas (fila,columna)
nCols = size(S,2);

if mod(k,nCols) == 0
    coords = [k/nCols, nCols];
else
    coords = [(k - mod(k,nCols))/nCols + 1, mod(k,nCols)];
end

end
