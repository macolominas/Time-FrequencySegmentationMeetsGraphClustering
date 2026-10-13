function G = buildGraph_parallel(S)
% Construye el grafo a partir de la matriz S (paralelizado con parfor)

box1 = size(S,1);
box2 = size(S,2);

% Prealocar celdas (cada iteración de parfor genera sus aristas)
edgesCell = cell(box1*box2,1);

parfor idx = 1:(box1*box2)
    [i,j] = ind2sub([box1,box2], idx); 
    localI = [];
    localJ = [];
    for k1 = -3:3
        for k2 = -3:3
            i2 = i + k1;
            j2 = j + k2;
            if i2 >= 1 && i2 <= box1 && j2 >= 1 && j2 <= box2
                if S(i,j) * S(i2,j2) >= 0.5
                    localI(end+1) = ind(i,j,S); %#ok<AGROW>
                    localJ(end+1) = ind(i2,j2,S); %#ok<AGROW>
                end
            end
        end
    end
    edgesCell{idx} = [localI(:), localJ(:)];
end

% Concatenar resultados de cada worker
edges = vertcat(edgesCell{:});
I = edges(:,1);
J = edges(:,2);

% Matriz dispersa y grafo
M = sparse(I,J,1, box1*box2, box1*box2);
G = graph(M);

end
