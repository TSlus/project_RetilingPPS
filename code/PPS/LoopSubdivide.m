function [vertices, faces, hedge_face] = LoopSubdivide(vertices, faces, hedge_face_old)
nf = size(faces, 1); np = size(vertices, 1);
nhe = 3 * nf; % 半边数
ne = nhe / 2;
if mod(nhe, 2) ~= 0
    warning('网格半边数数目不是偶数，网格不是封闭网格。')
    return;
end

% 半边结构
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1];
% 半边
% he = zeros(nhe, 2);
% he(:,1) = X; he(:,2) = Y;
% 半边--面
F = [1:nf, 1:nf, 1:nf]';
hedge_face = sparse(X, Y, F);
% 半边--对点
Z = [x3; x1; x2];
hedge_third_p = sparse(X, Y, Z);
% 半边--反向半边的索引
hedge_oppo_name = sparse(Y, X, (1:3*nf)');

% 每个原始顶点的度
v_valence(np) = 0;
nearPs{np} = [];
for P = 1:np
    neighbor_P = find(hedge_face(:,P));
    nearPs{P} = neighbor_P;           % P的邻域点
    v_valence(P) = length(neighbor_P);  % P的度
end

% 思路：半边循环，将每个面分裂成四个
flag_he = zeros(nhe, 1); % 记录半边上是否插入了新点，同时用它记录插入的点
% hedge_midpoint = hedge_face; % 在每条半边上插入一个点

% 面循环
idxp = np + 1; % 插入的新点所以
faces_new = zeros(4*nf, 3);
vertices_new = zeros(np + ne, 3);
flag_v = vertices_new(:, 1); % 记录顶点是否已经更新
face_name = zeros(12*nf, 3); % 记录新产生半边，所在原始面
for i = 1:nf
    a = faces(i, 1); b = faces(i, 2); c = faces(i, 3);
    %% 加入新点k1, k2, k3
    % add_point;
    % 第一条半边情况
    if flag_he(i) % 如果第一条半边上已经插入了点
        k1 = flag_he(i);
    else          % 如果第一条半边上还没插入新点
        k1 = idxp;
        flag_he(i) = idxp;
        % 对面半边也要插入相应的点
        name = hedge_oppo_name(a, b);
        flag_he(name) = idxp;
        idxp = idxp + 1;
    end
    
    % 第二条半边情况
    if flag_he(i + nf) % 如果第二条半边上已经插入了点
        k2 = flag_he(i + nf);
    else          % 如果第二条半边上还没插入新点
        k2 = idxp;
        flag_he(i + nf) = idxp;
        % 对面半边也要插入相应的点
        name = hedge_oppo_name(b, c);
        flag_he(name) = idxp;
        idxp = idxp + 1;
    end
    
    % 第三条半边情况
    if flag_he(i + 2*nf) % 如果第一条半边上已经插入了点
        k3 = flag_he(i + 2*nf);
    else          % 如果第一条半边上还没插入新点
        k3 = idxp;
        flag_he(i + 2*nf) = idxp;
        % 对面半边也要插入相应的点
        name = hedge_oppo_name(c, a);
        flag_he(name) = idxp;
        idxp = idxp + 1;
    end
    
    %%
    % 1.更新面
    face_add  = [a, k1, k3;
        k1, b, k2;
        k3, k2, c;
        k1, k2, k3];
    faces_new(4*i - 3:4*i,:) = face_add;
    
    y1 = [face_add(:,1),face_add(:,2)]; y2 = [face_add(:,2),face_add(:,3)];
    y3 = [face_add(:,3),face_add(:,1)];
    face_name((12*i-11):12*i, [1,2]) = [y1; y2; y3];
    face_name((12*i-11):12*i, 3) = hedge_face_old(a, b);
    
    %%
    % 2.更新点
    % update_vertex;
    % a点
    if flag_v(a) == 0
        flag_v(a) = 1;
        va = vertices(a,:);
        u = v_valence(a);
        beta = 1/u*(5/8-(3/8+1/4*cos(2*pi/u))^2);
        vertices_new(a,:) = (1-u*beta)*va + ones(1,u) * beta * vertices(nearPs{a},:);
    end
    
    % b点
    if flag_v(b) == 0
        flag_v(b) = 1;
        va = vertices(b,:);
        u = v_valence(b);
        beta = 1/u*(5/8-(3/8+1/4*cos(2*pi/u))^2);
        vertices_new(b,:) = (1-u*beta)*va + ones(1,u) * beta * vertices(nearPs{b},:);
    end
    
    % c点
    if flag_v(c) == 0
        flag_v(c) = 1;
        va = vertices(c,:);
        u = v_valence(c);
        beta = 1/u*(5/8-(3/8+1/4*cos(2*pi/u))^2);
        vertices_new(c,:) = (1-u*beta)*va + ones(1,u) * beta * vertices(nearPs{c},:);
    end
    
    %
    % k1
    if flag_v(k1) == 0
        flag_v(k1) = 1;
        v1 = vertices(a,:); v2 = vertices(b, :); v3 = vertices(c, :);
        d = hedge_third_p(b, a); v4 = vertices(d,:);
        vertices_new(k1, :) = (v1 + v2)*3/8 + (v3 + v4)/8;
    end
    
    % k2
    if flag_v(k2) == 0
        flag_v(k2) = 1;
        v1 = vertices(b,:); v2 = vertices(c, :); v3 = vertices(a, :);
        d = hedge_third_p(c, b); v4 = vertices(d,:);
        vertices_new(k2, :) = (v1 + v2)*3/8 + (v3 + v4)/8;
    end
    
    % k3
    if flag_v(k3) == 0
        flag_v(k3) = 1;
        v1 = vertices(c,:); v2 = vertices(a, :); v3 = vertices(b, :);
        d = hedge_third_p(a, c); v4 = vertices(d,:);
        vertices_new(k3, :) = (v1 + v2)*3/8 + (v3 + v4)/8;
    end
    
    
end

vertices = vertices_new;
faces = faces_new;
hedge_face = sparse(face_name(:, 1), face_name(:, 2), face_name(:, 3));
end


