function [faces_add2, re_flag, mcontinue, test_odd] = addFace(...
    xy, nearP, crirical_tri, crirical_tri3, test_odd)

faces_add2 = [];
n_near = length(nearP);
re_flag = 1; mcontinue = 0;
lastwarn('');
warning ('off');
polyin = polyshape(xy(:,1), xy(:,2));
T = triangulation(polyin);

[msg, ~] = lastwarn;
if ~isempty(msg)
    lastwarn('');
    warning ('on');
    mcontinue = 1;
    return;
end
warning ('on');

% 加入的面
adj_near = [1, n_near:-1:2];
adj_idxs = nearP(adj_near);
faces_add2 = adj_idxs(T.ConnectivityList);

% 再次小心判断，如果新加入的面中，只有一个面与ce对点相关
% 由于本身ce的存在，就形成了2度点
for j = 1:size(crirical_tri, 1)
    valence_ce = find(faces_add2 == crirical_tri(j,4));
    if length(valence_ce) == 1; re_flag = 0; break; end
end
if ~re_flag; return; end

% 不仅是二度点，对于crirical_tri和新加入的三角面，
% 查看是否有半边被用两次以上
crirical_tri = crirical_tri(:, 1:3);

nh_new1 = size(crirical_tri,1)*3;
nh_new2 = size(faces_add2,1)*3; % 考察的半边
nh_new3 = size(crirical_tri3,1);

nh_new = nh_new1 + nh_new2 + nh_new3;
halfedge_2idx = zeros(nh_new, 2);
crirical_tri_turn = crirical_tri(:, [2,3,1]);
halfedge_2idx(1:nh_new1, :) = [crirical_tri(:), crirical_tri_turn(:)];

faces_add2_turn = faces_add2(:, [2,3,1]);
halfedge_2idx((nh_new1+1):(nh_new1 + nh_new2), :) = [faces_add2(:), faces_add2_turn(:)];

halfedge_2idx((nh_new1 + nh_new2+1):nh_new, :) = crirical_tri3;

% halfedge_2idx = halfedge_2idx(1:(nh_new1 + nh_new2),:);
% 构造半边稀疏矩阵
halfedges = sparse(halfedge_2idx(:,1), halfedge_2idx(:,2), 1);
if full(max(max(halfedges))) > 1
    re_flag = 0;
    % disp('半边重复判定成功一次。')
    test_odd = test_odd + 1;
    % 之前错误，主要是因为新增的面和空3半边重复
    % 新增的面，与老式半边重复

end
% 空3反向半边没有判定，重复半边，就需要定义得更加细致。

% % 将邻域反向半边也加入进来（其实不必，这样最多生成2度点，前面2度点已经判断过了）
% crirical_tri = crirical_tri(:, 1:3);
% nh_new1 = n_near;
% nh_new2 = size(faces_add2,1)*3; % 考察的半边
% nh_new3 = size(crirical_tri3,1)*2;
% nh_new4 = size(crirical_tri,1)*2;

end