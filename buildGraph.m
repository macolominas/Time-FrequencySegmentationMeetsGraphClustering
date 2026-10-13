function G = buildGraph(S,side,thr)
% Construye el grafo a partir de la matriz S

box1 = size(S,1);
box2 = size(S,2);

I = [];
J = [];

for i = 1:box1
    for j = 1:box2
        for k1 = -side:side
            for k2 = -side:side
                i2 = i + k1;
                j2 = j + k2;
                if i2 >= 1 && i2 <= box1 && j2 >= 1 && j2 <= box2
                    if S(i,j) * S(i2,j2) >= thr
                        I(end+1) = ind(i,j,S); %#ok<AGROW>
                        J(end+1) = ind(i2,j2,S); %#ok<AGROW>
                    end
                end
            end
        end
    end
end

M = sparse(I,J,1, box1*box2, box1*box2);
G = graph(M);

end
