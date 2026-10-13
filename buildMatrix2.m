function S = buildMatrix2(nRows,nCols)
% Genera la matriz S de forma vectorizada

% Crear grillas de índices
[jGrid,iGrid] = ndgrid(1:nRows, 1:nCols);

% Calcular los términos de forma vectorizada
term1 = exp(-(5*cos(jGrid/3) - 0.1*iGrid).^2);
term2 = exp(-(2*jGrid - iGrid - 65).^2);
term3 = 0.6*abs(rand(nRows,nCols));

% Combinar
S = term1 + term2 + term3;

end
