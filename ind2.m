function k = ind2(a,b,S)
% IND Convierte coordenadas (a,b) en índices lineales.
%    k = ind(a,b,S)
%
% Soporta escalares, vectores o matrices de a y b.

nCols = size(S,2);
k = nCols.*(a-1) + b;

end

% ejemplo
% S = zeros(50,100);
% k = ind([1 2 3],[4 5 6],S); % devuelve [104,205,306]
