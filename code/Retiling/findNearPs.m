function [oneRingP, vertex_valence] = findNearPs(faces)

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
    %初始化，第一个点
    neighbor_P = find(half_face(:,P));%返回非零元素位置，点P的邻域点（没有顺序）
    vertex_valence(P) = length(neighbor_P);

    a = neighbor_P(1);
    oneRingP{P}(1) = a;
    t = 2;
    neighbor_res = neighbor_P(2:end);
    %找 b 的下一个点
    while  ~isempty(neighbor_res)        
        a_nextP = find(half_face(a,neighbor_res));
        
        k = 1;
        b = neighbor_res(a_nextP(k)); %顺时针还是逆时针不知道
        while half_face(a,b) ~= half_face(P, a) % 当网格时封闭网格是需要考虑
            k = k+1;
            b = neighbor_res(a_nextP(k));    % 如果顶点没有oneRing，不能直接用!!!
        end
        oneRingP{P}(t) = b; t = t+1;
        a = b; 
        neighbor_res = neighbor_res(neighbor_res ~= b);
    end
end
end
