function plotComponentMask(S, comp, comp_id)
% Dibuja la máscara binaria de una componente conexa

box1 = size(S,1);
box2 = size(S,2);

nodes = find(comp==comp_id);

coords = arrayfun(@(x) revind(x,S), nodes, 'UniformOutput', false);
coords = cell2mat(coords'); 

S2 = zeros(box1,box2);
for c = 1:size(coords,1)
    S2(coords(c,1), coords(c,2)) = 1;
end

figure; imagesc(S2); colormap(gray); colorbar;
title(['Componente conexa ' num2str(comp_id)]);

end
