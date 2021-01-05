% 寻找点的邻域点，没有定向
function [oneRingP, vertex_valence] = NearP(faces)

% faces -  nf * 3
% oneRingP - 邻域顶点
% vertex_valence - 顶点的度

%半边-面id
nf = size(faces,1);
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1];
f_name = [1:nf, 1:nf, 1:nf]';
half_face = sparse(X, Y, f_name);

np = max(max(faces));
oneRingP{np} = [];
vertex_valence(np) = 0;
for P = 1:np
    neighbor_P = find(half_face(:,P));%返回非零元素位置，点P的邻域点（没有顺序）
    vertex_valence(P) = length(neighbor_P); 
    oneRingP{P} = neighbor_P; 
end
end
