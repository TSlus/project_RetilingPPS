function forces = ComputeRepulsiveForce_adj(mymesh, nameF_cand, vertices_cand, FSD, VSD)
mesh_to_VFdetail;
radius = mymesh.radius;

nCand = length(nameF_cand);
k_0 = find([nameF_cand, 0] == 0, 1);
nCand_NZ = k_0 - 1;

t = ones(nf, 1);
NEvery = 8; % 每个面，不超过6个candidate点
arry_p = zeros(nf, NEvery);
for i = 1:nCand_NZ
    fi = nameF_cand(i);
    ti = t(fi); % 第fi个面的指标
    arry_p(fi, ti) = i;
    t(fi) = ti + 1;
end

% arry_p 与 FSD结合
% arry_p:记录每个面上candidate points
% FSD:   记录与面 f 接近的面（不包含临近三个面）
FSD_arry_p{nf} = [];
for i = 1:nf
    fis = FSD{i};
    temp = arry_p(fis, :);
    FSD_arry_p{i} = temp(temp > 0)';
end

[CanF_Rotate, CDP_rotate] = CandidateAfterRotate3(vertices_cand(1:nCand_NZ,:), ...
    vertices, faces, nameF_cand(1:nCand_NZ), norm_face, hedge_face);

forces = zeros(nCand, 3);
for i = 1:nf
    if arry_p(i, 1) == 0; continue; end
    norm_fi = norm_face(i, :);
    
    % 计算forces2的准备工作
    % 将面 i 周围三个面上的点旋转过来
    a = faces(i,1); b = faces(i,2); c = faces(i,3);
%     f1 = full(hedge_face(b, a)); f2 = full(hedge_face(c, b)); f3 = full(hedge_face(a, c));
%     cand_id1 = arry_p(f1,:); cand_id2 = arry_p(f2,:); cand_id3 = arry_p(f3,:);
    cand_id1 = arry_p(hedge_face(b, a),:); cand_id2 = arry_p(hedge_face(c, b),:); 
    cand_id3 = arry_p(hedge_face(a, c),:);
    
    cand_id1 = cand_id1(cand_id1 > 0); n1 = length(cand_id1);
    cand_id2 = cand_id2(cand_id2 > 0); n2 = length(cand_id2);
    cand_id3 = cand_id3(cand_id3 > 0); n3 = length(cand_id3);
    ns = n1+n2+n3;
    
    cand_id = [cand_id1(:); cand_id2(:); cand_id3(:)];
    candf_id = [b*ones(n1,1); c*ones(n2,1); a*ones(n3,1)];
    rotated_idx = CanF_Rotate(sub2ind(size(CanF_Rotate), cand_id, candf_id));
    rotate_ps = CDP_rotate(rotated_idx, :);
    
    % 计算faces34的准备工作
    % forces3
    idx_others3 = FSD_arry_p{i}; % 与平面 i 较近的 candidate points
    % forces4
    idx_others4 = VSD{i} + nCand_NZ; % 与平面 i 较近的 CSP
    idx_others = [idx_others3(:); idx_others4(:)];
    n_others = length(idx_others);
    
    proj_ps = project_point_to_triangle3(vertices_cand(idx_others,:), ...
        vertices(faces(i,1),:), repmat(norm_fi, n_others, 1));
    
    % 根据面 i 上点的个数不同方式计算
    if arry_p(i, 2) == 0 % 第i个面上只有一个点
        cand_repul = arry_p(i, 1);
        p_to_repul = vertices_cand(cand_repul, :);
        
        forces2 = p_to_repul - rotate_ps; % forces2
        forces34 = p_to_repul - proj_ps; % forces34
        % 没有 forces1
        
        mforces = zeros(ns + n_others, 3);
        mforces(1:ns,:) = forces2; mforces(ns + 1:ns + n_others,:) = forces34;
        
        mforces_norm = sqrt(sum(mforces.^2, 2));
        idx_r = ((mforces_norm > 0) & (mforces_norm < radius));
        mforces = mforces(idx_r, :);
        mforces_norm = mforces_norm(idx_r,:);
        
        mforces = mforces ./ mforces_norm;
        
        forces(cand_repul, :) = sum((radius - mforces_norm) .* mforces, 1);
        %     elseif arry_p(i, 3) == 0
        %         cand_repul = arry_p(i, 1:2);
        %         p_to_repul = vertices_cand(cand_repul, :);
        %
        %
        %     elseif arry_p(i, 4) == 0
        %         cand_repul = arry_p(i, 1:3);
        %         p_to_repul = vertices_cand(cand_repul, :);
    else
        cand_repul = arry_p(i, :); cand_repul = cand_repul(cand_repul > 0);
        n_cani = length(cand_repul);
        p_to_repul = vertices_cand(cand_repul, :);
        
        for j = 1:n_cani
            forces1 = p_to_repul(j,:) - p_to_repul; % forces1
            forces2 = p_to_repul(j,:) - rotate_ps; % forces2
            forces34 = p_to_repul(j,:) - proj_ps; % forces34
            
            mforces = zeros(n_cani + ns + n_others, 3);
            mforces(1:n_cani, :) = forces1;
            mforces(n_cani+1:n_cani + ns,:) = forces2; 
            mforces(n_cani + ns + 1:n_cani + ns + n_others,:) = forces34;
            
            mforces_norm = sqrt(sum(mforces.^2, 2));
            idx_r = ((mforces_norm > 0) & (mforces_norm < radius));
            mforces = mforces(idx_r, :);
            mforces_norm = mforces_norm(idx_r,:);
            
            mforces = mforces ./ mforces_norm;
            
            forces(cand_repul(j), :) = sum((radius - mforces_norm) .* mforces, 1);
        end
        
    end
    
end

end