% 如果继续改进，就是减少在寻找邻域点时，find函数使用
% 只在第一次将每个点，邻域点全找出来
% 然后在去除顶点过程中，进行更新
function [vertices_Mutual, faces_Mutual, n_rem] = RemovingOldVertices_cpp(...
    mymesh, vertices_cand, faces_Mutual, CSP_idx) %#ok<INUSL>
mesh_to_VFdetail;

nCand = size(vertices_cand, 1);
vertices_Mutual = [vertices; vertices_cand];

remove_idx = zeros(1, np); % 被移除就赋值为1
idx_change = 1:np;

% faces_old = faces_faces_Mutual;
for i = 1:np
    if ismember(i, CSP_idx); continue; end
    
    vi = vertices(i, :);
    [r_i, ~] = find(faces_Mutual == i);
    delete_fM = zeros(1, size(faces_Mutual, 1));
    delete_fM(r_i) = 1;
    tri_1 = faces_Mutual(r_i,:);
    tri_turn = tri_1(:, [2,3,1]);
    % nearPsp 用来判断选中的三角形哪条边为critical edge
    nf_sp = size(tri_1, 1);
    nearPsp = sparse(tri_1(:), tri_turn(:), [1:nf_sp, 1:nf_sp, 1:nf_sp]');
    nearP = findNearP_sp(nearPsp, i); n_near = length(nearP);
    
    % 计算顶点 i 的法向量
    v_nearp = vertices_Mutual(nearP, :);
    normals_pi = cross(v_nearp - vi, v_nearp([2:n_near, 1], :) - vi, 2);
    normals_pi_norm = sqrt(sum(normals_pi.^2, 2));
    normal = sum(normals_pi ./ normals_pi_norm, 1); 
    normal = normal./norm(normal);
    normal0 = normal;
    
    % 检查投影
    n_test = 0; success = 0;
    while n_test < 15
        switch n_test
            case 0
            case 1; normal = [.8507, .4472, .2764];
            case 2; normal = [-.8507, .4472, .2764];
            case 3; normal = [.8507, -.4472, .2764];
            case 4; normal = [.8507, .4472, -.2764];
            case 5; normal = [.5257, .4472, .7236];
            case 6; normal = [-.5257, .4472, .7236];
            case 7; normal = [.5257, -.4472, .7236];
            case 8; normal = [.5257, .4472, -.7236];
            case 9; normal = [.0, .4472, .8944];
            case 10; normal = [.0, -.4472, .8944];
            case 11; normal = [.0, .4472, -.8944];
            case 12; normal = [0, 1, 0];
            case 13; normal = [1, 0,0];
            case 14; normal = [0, 0, 1];
        end
        if sum(normal.*normal0) < 0; normal = -normal; end
        
        base1 = project_point_to_triangle3(v_nearp(1, :), vi, normal);
        base1 = base1/norm(base1);
        base2 = cross(normal, base1);
        
        % 先投影，再计算二维坐标对应
        proj_nears = project_point_to_triangle3(v_nearp, vi, normal);
        xy = zeros(n_near, 2);
        xy(:, 1) = sum(proj_nears.*base1, 2);
        xy(:, 2) = sum(proj_nears.*base2, 2);
        
        % 判断投影方向是否可用
        is_proj_flag = ProjectFormNormal(xy);
        if is_proj_flag
            success = 1;
            break;
        end
        n_test = n_test+1;
    end
    
    if ~success; continue; end % 如果投影不成功，点 i 就不能被去除
    
    % 检查 critical edges
    re_flag = 1; loopf1 = 0;
    for ip = nearP
        [rip, ~] = find(faces_Mutual == ip);
        triRing = faces_Mutual(rip,:);
        
        for j = 1:size(triRing, 1)
            a = triRing(j, 1);
            idx1 = find(nearP == a, 1); if isempty(idx1); continue; end % 362.708s 553591729次
            b = triRing(j, 2);
            idx2 = find(nearP == b, 1); if isempty(idx2); continue; end
            c = triRing(j, 3);
            idx3 = find(nearP == c, 1); if isempty(idx3); continue; end
            %n_ce = n_ce + 1; if n_ce == n_near+1; loopf1 = 1; break; end % critical edge 最多有度nsp个
            
            aif = nearPsp(c, b); bif = nearPsp(a, c); cif = nearPsp(b, a);
            ceif = find(~[aif, bif, cif]);
            if length(ceif) ~= 1 % 如果三条边都不是critical edge
                continue;
            end
            
            idxs = [idx3, idx2, idx1; idx1, idx3, idx2; idx2, idx1, idx3];
            idxs = idxs(ceif, :); % c, b, a 在邻域点钟的位置
            e1 = xy(idxs(2), :) - xy(idxs(3), :);
            e2 = xy(idxs(1), :) - xy(idxs(3), :);
            mcross = cross([e1, 0], [e2, 0]); 
            if mcross(3) < 0 % 点1不可以去除，不做三角化
                re_flag = 0; loopf1 = 1; break;
            end
        end
        if loopf1 ; break; end
    end
    
    
    if ~re_flag  % 如果不能点移除
        continue;
    end
    
    polyin = polyshape(xy(:,1), xy(:,2));
    T = triangulation(polyin);
    
    % 加入的面
    adj_near = [1, n_near:-1:2];
    adj_idxs = nearP(adj_near);
    faces_add2 = adj_idxs(T.ConnectivityList);
    % 新增的面
    faces_add1 = faces_Mutual(delete_fM == 0, :);
    faces_Mutual = [faces_add1; faces_add2];
    
    % 调整索引
    idx_change(i+1:np) = idx_change(i+1:np) - 1;
    % 记录成功去除的点
    remove_idx(i) = 1;
    
%     faces_old = faces_faces_Mutual;
end
remain_p = remove_idx == 0;
vertices_Mutual = [vertices(remain_p, :); vertices_cand];
n_rem = length(find(remain_p));

idx_change2 = (np+1):(np + nCand); 
my_idx_change = [idx_change, idx_change2 - (np - n_rem)];

faces_Mutual = my_idx_change(faces_Mutual);
end

