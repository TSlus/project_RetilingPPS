function [faces_Mutual, faces_delete, delete_arry, faces_add, add_arry, iC_near] = ...
    removeNonNearV(iC, vertices_Mutual, faces_Mutual)

nf_M = size(faces_Mutual, 1);
delete_fM = zeros(1, nf_M);

nc = length(iC);
faces_delete = zeros(nf_M, 4); delete_arry = ones(nc+1, 1); td = 1;
faces_add = zeros(nf_M, 3); add_arry = ones(nc+1, 1); ta = 1;
iC_near{nc} = [];
for k = 1:nc
    % 邻域点
    [r_i, ~] = find(faces_Mutual == iC(k)); n_tri = length(r_i);    
    faces_i = faces_Mutual(r_i,:);
    try 
    near_p = findNearP(faces_i, iC(k)); iC_near{k} = near_p;
    catch
        ik = iC(k)
    end
    num_nearp = length(near_p);

    % 删除
    delete_fM(r_i) = 1; 
    faces_delete(td:(td+n_tri-1),:) = [faces_i, r_i]; % 第四列放面索引
    td = td + n_tri;
    delete_arry(k+1) = delete_arry(k) + n_tri;
    
    % “投影”
    t = (1:num_nearp) * 2 * pi / num_nearp;
    x = cos(t); y = sin(t);
    xy = [x; y]; xy = xy';
    di_vec = vertices_Mutual(near_p,:) - vertices_Mutual(iC(k),:);
    di = sqrt(sum(di_vec.^2, 2));
    
    % 三角化
    
    xy = xy .* di;
    polyin = polyshape(xy(:,1), xy(:,2));
    T = triangulation(polyin);
    adj_face = near_p([1, num_nearp:-1:2]);
    faces_add_i = adj_face(T.ConnectivityList); 
    
    n_addF = size(faces_add_i, 1);
    faces_add(ta:(ta + n_addF - 1), :) = faces_add_i; ta = ta + n_addF;
    add_arry(k+1) = add_arry(k) + n_addF;
    % 利用之前的半边信息，判断顶点是否可以移除，如果不可移除，
    % 整批需要单独移除, 二分法?
end

faces_delete = faces_delete(1:(td-1),:);
faces_add = faces_add(1:(ta-1), :);

faces_rem = faces_Mutual(delete_fM == 0, :);
% aa = length(find(delete_fM == 0))
faces_Mutual = [faces_rem; faces_add];
end 