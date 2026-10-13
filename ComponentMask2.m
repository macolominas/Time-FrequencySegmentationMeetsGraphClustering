function S2 = ComponentMask2(S, comp, comp_id)
% PLOTCOMPONENTMASK Dibuja la máscara binaria de una componente conexa.
%
%   plotComponentMask(S, comp, comp_id)
%
%   S       : matriz original
%   comp    : vector de etiquetas de componentes (salida de conncomp)
%   comp_id : identificador de la componente a visualizar

box1 = size(S,1);
box2 = size(S,2);

% Nodos de la componente
nodes = find(comp==comp_id);

% Convertir a coordenadas de forma vectorizada
coords = revind(nodes,S);

% Construir máscara binaria
S2 = zeros(box1,box2);
idxLinear = sub2ind([box1,box2], coords(:,1), coords(:,2));
S2(idxLinear) = 1;

% Mostrar
% figure; imagesc(S2); colormap(gray); colorbar;
% title(['Componente conexa ' num2str(comp_id)]);

end
