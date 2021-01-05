function [faces_Mutual, stateP] = unRemovedV(faces_Mutual_old, hedge_old, iC,...
    faces_delete, delete_arry, faces_add, add_arry, iC_near, stateP)

% 在移除的这一批顶点中，有一些点不能被移除
% 就一个一个顶点检查

% 把旧的邻域半边也取出来，判断半边重复数是否大于1
n_fM_old = size(faces_Mutual_old, 1);
delete_fM = zeros(1,n_fM_old);
faces_add2 = faces_add; t2 = 1;
n = length(iC); nc = 1;
for i = 1:n
    near_i = iC_near{i};
    n_near = length(near_i);
    
    test_he = zeros(n_near, n_near); % 测试有无半边重复
    temp = repmat(1:n_near, n_near, 1);
    xhe = temp(:); temp = temp'; yhe = temp(:);
    
    % 1.原始半边（含邻域点之间）
    x1 = near_i(xhe); y1 = near_i(yhe); 
    test_he(sub2ind(size(test_he), xhe, yhe))=hedge_old(sub2ind(size(hedge_old), x1, y1));
    % 不对称就报错
    test1 = test_he - test_he';
    if ~isempty(find(test1, 1)); disp('矩阵不对称'); end
    nhe_old = sum(test_he(:)==1);
    
    adj_idx = 1:max(near_i); adj_idx(near_i) = 1:n_near;
    % 2.删除的半边
    idx_de = zeros(n_near, 2);
    idx_de(:,1) = near_i; idx_de(:,2) = near_i([2:n_near,1]);
    test_he(sub2ind(size(test_he), adj_idx(idx_de(:,1)), adj_idx(idx_de(:,2)))) = 0;
    
    d1 = delete_arry(i); d2 = delete_arry(i+1);

    % 3.新增的半边，
    a1 = add_arry(i); a2 = add_arry(i+1);
    n_add = a2 - a1;
    faces_add_i = faces_add(a1:(a2-1), :);     
    faces_add_turn = faces_add_i(:,[2,3,1]);
    test_he(sub2ind(size(test_he), adj_idx(faces_add_i(:)), adj_idx(faces_add_turn(:)))) = 1;
    
    n_now = sum(test_he(:)==1);
    
    % 判断条件：移除前后半边数关系，半边矩阵是否对称
    if n_now ~= (nhe_old - n_near + 3*n_add) || ~isempty(find(test_he-test_he',1))
        stateP(iC(i)) = -1;
        continue;
    end
    
    mri = faces_delete(d1:(d2-1),4);
    delete_fM(mri) = 1; % 移除的面
    faces_add2(t2:(t2+n_add-1),:) = faces_add_i; t2 = t2 + n_add;
    nc = nc+1;
    
end
% if nc == n+1;disp('这一批全部通过。'); end

faces_add2 = faces_add2(1:(t2-1),:);
faces_rem = faces_Mutual_old(delete_fM == 0, :);

% 删除可以去除的面过后，还要加入新的面！
faces_Mutual = [faces_rem; faces_add2];
end


% my test
% function [faces_Mutual, stateP] = unRemovedV(faces_Mutual_old, hedge_old, iC,...
%     faces_delete, delete_arry, faces_add, add_arry, iC_near, stateP)
% 
% % 在移除的这一批顶点中，有一些点不能被移除
% % 就一个一个顶点检查
% 
% % 把旧的邻域半边也取出来，判断半边重复数是否大于1
% n_fM_old = size(faces_Mutual_old, 1);
% delete_fM = zeros(1,n_fM_old);
% faces_add2 = faces_add; t2 = 1;
% n = length(iC); nc = 1;
% for i = 1:n
%     % near_i = find(hedge_old(:,iC(i))); % near_i排序！！！！！
%     near_i = iC_near{i};
%     n_near = length(near_i);
%     % 原始邻域半边（包含邻域点之间），移除原始半边，新三角面（增加）
%     test_he = zeros(n_near, n_near); % 测试有无半边重复
%     temp = repmat(1:n_near, n_near, 1);
%     xhe = temp(:); temp = temp'; yhe = temp(:);
%     
%     % 1.原始半边（含邻域点之间），是n_near+1和点
%     x1 = near_i(xhe); y1 = near_i(yhe); 
%     test_he(sub2ind(size(test_he), xhe, yhe))=hedge_old(sub2ind(size(hedge_old), x1, y1));
%     % 不对称就报错
%     test1 = test_he - test_he';
%     if ~isempty(find(test1, 1)); disp('矩阵不对称'); end
%     nhe_old = sum(test_he(:)==1);
%     
%     adj_idx = 1:max(near_i); adj_idx(near_i) = 1:n_near;
%     % 2.删除的半边
%     idx_de = zeros(n_near, 2);
%     idx_de(:,1) = near_i; idx_de(:,2) = near_i([2:n_near,1]);
%     test_he(sub2ind(size(test_he), adj_idx(idx_de(:,1)), adj_idx(idx_de(:,2)))) = 0;
% %     adj_idx(idx_de(:,1))
% %     adj_idx(idx_de(:,2))
%     
%     d1 = delete_arry(i); d2 = delete_arry(i+1);
%     % 借用了已知删除三角形
% %     faces_delete_i = faces_delete(d1:(d2-1), [1,2,3]); 
% %     faces_delete_turn = faces_delete_i(:,[2,3,1]); 
% %     idxs = zeros(3*(d2 - d1), 2);
% %     idxs(:,1) = faces_delete_i(:); idxs(:,2) = faces_delete_turn(:);
% %     % 除去和顶点iC(i)有关的半边
% %     [ri, ~] = find(idxs == iC(i));
% %     idxs = removerows(idxs, ri);
% %     if size(idxs,1) ~= n_near; disp('移除的半边不符');end
% %     
% %     test_he(sub2ind(size(test_he), adj_idx(idxs(:,1)), adj_idx(idxs(:,2)))) = 0;
% 
%     test2 = zeros(n_near, n_near);
%     % 3.新增的半边，
%     a1 = add_arry(i); a2 = add_arry(i+1);
%     n_add = a2 - a1;
%     faces_add_i = faces_add(a1:(a2-1), :); 
% %     sort(unique(faces_add_i(:)))
% %     sort(near_i)
%     
%     faces_add_turn = faces_add_i(:,[2,3,1]);
% %     test2(sub2ind(size(test2), adj_idx(faces_add_i(:)), adj_idx(faces_add_turn(:)))) = 1;
% %     test_he = test_he + test2;
%     test_he(sub2ind(size(test_he), adj_idx(faces_add_i(:)), adj_idx(faces_add_turn(:)))) = 1;
%     
%     n_now = sum(test_he(:)==1);
%     % 3.2新增的半边里不能包含反向邻域半边
% %     temp = test2(sub2ind(size(test2), adj_idx(idx_de(:,2)), adj_idx(idx_de(:,1))));
%     test2(sub2ind(size(test2), adj_idx(idx_de(:,2)), adj_idx(idx_de(:,1)))) = 1;
%     
%     % 正向半边全部有
%     % temp2 = test2(sub2ind(size(test2), adj_idx(idx_de(:,1)), adj_idx(idx_de(:,2))));
% %     
% %     if max(max(test_he)) > 1 || ~isempty(find(temp,1)) || sum(~temp2) ~= 0  ...
% %             ||~isempty(find(test_he-test_he',1))
%     
% %     if max(max(test_he)) > 1 || ~isempty(find(test_he-test_he',1)) ...
% %             || max(max(test2)) > 1 || ~isempty(find(test2-test2',1)) || (nhe_old-n_near+3*n_add)~= n_now
%     if n_now ~= (nhe_old - n_near + 3*n_add) || ~isempty(find(test_he-test_he',1))
%         stateP(iC(i)) = -1;
%         continue;
%     end
%     % 三角化时如果有警告，就不移
%     
%     % 确定新生成的三角形定向
%     
%     
%     mri = faces_delete(d1:(d2-1),4);
%     delete_fM(mri) = 1; % 移除的面
%     faces_add2(t2:(t2+n_add-1),:) = faces_add_i; t2 = t2 + n_add;
%     nc = nc+1;
%     
% end
% if nc == n+1;disp('这一批全部通过。'); end
% 
% faces_add2 = faces_add2(1:(t2-1),:);
% faces_rem = faces_Mutual_old(delete_fM == 0, :);
% 
% % 删除可以去除的面过后，还要加入新的面！
% faces_Mutual = [faces_rem; faces_add2];
% end
% 
