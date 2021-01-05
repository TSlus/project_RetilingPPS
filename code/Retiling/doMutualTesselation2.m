function faces_Mutual = doMutualTesselation2(mymesh, nameF_cand, vertices_cand) %#ok<INUSL>
mesh_to_VFdetail;
vertices; faces;

nCand_NZ = length(nameF_cand);

%% 1.先考察每一个面上有哪些 candidate points
% t: 每一列是指示矩阵，用来给矩阵分配内存
t = ones(1, nf);
tem = zeros(1, 10); % 每个面，不超过10个candidate点
cell_p(1:nf) = {tem}; % 每一个面上有哪些 candidate points
for i = 1:nCand_NZ
    fi = nameF_cand(i);
    ti = t(fi); % 第fi个面的指标
    cell_p{fi}(ti) = i;
    t(fi) = ti + 1;
end
% 取出非零点索引
for i = 1:nf
    ti = t(i);
    cell_p{i} = cell_p{i}(1:(ti-1));
end

%% 2.在每个面做三角剖分
N = 10;
faces_Mutual = zeros((2*N + 1)*nf, 3);
t = 1;
for i = 1:nf
    n_cand_every_face = length(cell_p{i});
    if n_cand_every_face == 0
        faces_Mutual(t, :) = faces(i, :);
        t = t+1;
        
    elseif n_cand_every_face == 1
        a = faces(i, 1); b = faces(i, 2); c = faces(i, 3);
        id_i = cell_p{i} + np; 
        faces_one_cand = [a, b, id_i; b, c, id_i; c, a, id_i]; 
        faces_Mutual(t:t+2, :) = faces_one_cand;
        t = t+3;
        
    else % 三角形面上多余1个candidate point
        idx1 = faces(i,:);
        idx2 = cell_p{i};
        vert1 = vertices(idx1,:);
        vert2 = vertices_cand(idx2,:);
        idx = [idx1, idx2 + np];
        vert = [vert1; vert2]; nv = size(vert, 1);
        
        e1 = vert(1,:)/norm(vert(1,:)); n1 = norm_face(i,:)/norm(norm_face(i,:));
        e2 = cross(n1, e1);
        x_data = dot(vert, repmat(e1, nv, 1),2); %size(x_data)
        y_data = dot(vert, repmat(e2, nv, 1),2); %size(y_data)
        T = delaunayTriangulation(x_data, y_data);
        f_add = idx(T.ConnectivityList);
        % 保证生成的 mesh 是封闭定向一致的
%         if norm_face(i,3) < 0
%             f_add = f_add(:,[2,1,3]);
%         end
        n_add = size(f_add,1);
        faces_Mutual(t:(t+ n_add - 1),:) = f_add;
        t = t + n_add;  
    end
end
faces_Mutual = faces_Mutual(1:t-1, :);

end