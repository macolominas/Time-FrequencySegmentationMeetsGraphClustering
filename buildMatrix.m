function S = buildMatrix(nRows,nCols)
% Genera la matriz S del código original en Mathematica

S = zeros(nRows,nCols);
for j = 1:nRows
    for i = 1:nCols
        term1 = exp(-(5*cos(j/3) - 0.1*i)^2);
        term2 = exp(-(2*j - i - 65)^2);
        term3 = 0.6*abs(rand); 
        S(j,i) = term1 + term2 + term3;
    end
end

end
