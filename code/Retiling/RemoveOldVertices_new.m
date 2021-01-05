% remove old vertices
function [vertices_Mutual, faces_Mutual, n_remain, P2, CSP_idx] = RemoveOldVertices_new(...
    mymesh, vertices_cand, faces_Mutual, CSP_idx)

vertices = mymesh.vertices; np = size(vertices, 1);
nCand = size(vertices_cand, 1);
vertices_Mutual = [vertices; vertices_cand];

edges =[faces_Mutual(:,[1,2]); faces_Mutual(:,[2,3]); faces_Mutual(:,[3,1])];
% 半边
hedge = sparse(edges(:,1), edges(:,2), 1);

iA = 1:np; iA(CSP_idx) = 0; iA = iA(iA>0);
stateP = -1*ones(1,np); stateP(iA) = 1; % 取值：-1，未被移除；0，已经移除；1，待移除
P2 = [];
%%
while ~isempty(iA)
    faces_Mutual_old = faces_Mutual;
    hedge_old = hedge; % max(max(hedge))
    
    %% 分批移除
    iC = findNonNearV(iA, hedge); % 将iC中的点移除掉
    stateP(iC) = 0; % 如果不能被移除，还会进行修改
    
    try
    %% 将 iC 中的点同时移除
    [faces_Mutual, faces_delete, delete_arry, faces_add, add_arry, iC_near] = ...
        removeNonNearV(iC, vertices_Mutual, faces_Mutual);
    
    catch
        %figure(22)
        save myodddata vertices_Mutual faces_Mutual;
        warning('出错。')
    end
    
    % 重新构建稀疏矩阵，并判断这一批顶点是否可以移除
    hedge = sparse([faces_Mutual(:,1);faces_Mutual(:,2);faces_Mutual(:,3)],...
        [faces_Mutual(:,2);faces_Mutual(:,3);faces_Mutual(:,1)], 1);
    
    %% 如果这一批不能移除，就逐一判断，移除
    if max(max(hedge)) >1
        %% 逐一判断，移除
        [faces_Mutual, stateP] = unRemovedV(faces_Mutual_old, hedge_old, iC,...
            faces_delete, delete_arry, faces_add, add_arry, iC_near, stateP);
        
        % 记得更新半边信息
        hedge = sparse([faces_Mutual(:,1);faces_Mutual(:,2);faces_Mutual(:,3)],...
            [faces_Mutual(:,2);faces_Mutual(:,3);faces_Mutual(:,1)], 1);
        
        %% 解决特殊情况，二度点在函数 unRemovedV 中没有全部移除，应该是有三点共线？？？
        if max(max(hedge)) > 1
            if max(max(hedge)) > 2; disp('半边数大于2，不可能'); end
%             % disp('半边仍然重复。');
%             
%             nC = max(max(faces_Mutual));
%             rx = zeros(1, size(faces_Mutual,1));
%             % 一般只会调用一次，所以耗时不多？
%             for j = 1:nC
%                 [fi,~] = find(faces_Mutual == j);
%                 if length(fi) == 2
%                     rx(fi) = 1;
%                 end
%             end
%             
%             faces_Mutual = faces_Mutual(rx == 0,:);
%             hedge = sparse([faces_Mutual(:,1);faces_Mutual(:,2);faces_Mutual(:,3)],...
%                 [faces_Mutual(:,2);faces_Mutual(:,3);faces_Mutual(:,1)], 1);
        [faces_Mutual,hedge, P2] = removeValence2(faces_Mutual, hedge);
        end
        
        %
        if max(max(hedge)) == 2; disp('2度顶点没有去除干净。'); end
    end
    
    % iA中的iC点移除，更新
    iA = find(stateP == 1);
end

%% 修改面索引，去除孤立顶点
idx_change = 1:np; % 用于修改最后的顶点索引
removeP = find(stateP == 0);
for j = removeP
    idx_change(j+1:np) = idx_change(j+1:np) - 1;% adjust indexes
end

n_remove = length(removeP);
n_remain = np - length(removeP);

%% adjust points
vid = 1:np; vid(removeP) = 0;
vertices_Mutual = [vertices(vid > 0, :); vertices_cand];

%% adjust faces
idx_change2 = (np+1):(np + nCand);
my_idx_change = [idx_change, idx_change2 - n_remove];
faces_Mutual = my_idx_change(faces_Mutual);

%% adjust CSP_idx
CSP_idx = my_idx_change(CSP_idx);


end


%% my test
% function [vertices_Mutual, faces_Mutual, n_remain] = RemoveOldVertices_new(...
%     mymesh, vertices_cand, faces_Mutual, CSP_idx)
%
% vertices = mymesh.vertices; np = size(vertices, 1);
% nCand = size(vertices_cand, 1);
% vertices_Mutual = [vertices; vertices_cand];
%
% edges =[faces_Mutual(:,[1,2]); faces_Mutual(:,[2,3]); faces_Mutual(:,[3,1])];
% % 半边
% x = edges(:,1); y = edges(:,2);
% hedge = sparse(x, y, 1);
%
% iA = 1:np; iA(CSP_idx) = 0; iA = iA(iA>0);
% stateP = -1*ones(1,np); stateP(iA) = 1; % 取值：-1，未被移除；0，已经移除；1，待移除
% while ~isempty(iA)
%     faces_Mutual_old = faces_Mutual;
%     hedge_old = hedge; % max(max(hedge))
%
%     % 分批移除
%     iC = findNonNearV(iA, hedge); % 将iC中的点移除掉
%     stateP(iC) = 0; % 如果不能被移除，还会进行修改
%
%     try
%         % 将 iC 中的点同时移除
%         [faces_Mutual, faces_delete, delete_arry, faces_add, add_arry, iC_near] = ...
%             removeNonNearV(iC, vertices_Mutual, faces_Mutual);
%     catch
%         trimesh(faces_Mutual, vertices_Mutual(:,1),vertices_Mutual(:,2),vertices_Mutual(:,3));
%         [x, ~] = find(hedge == 2);
%         axis equal; hold on
%         plot3(vertices_Mutual(x,1),vertices_Mutual(x,2),vertices_Mutual(x,3),'r*');
%
%         n_remain = -1;
%         return;
%     end
%
%     % 重新构建稀疏矩阵，判断这一批顶点是否可以移除
%     hedge = sparse([faces_Mutual(:,1);faces_Mutual(:,2);faces_Mutual(:,3)],...
%         [faces_Mutual(:,2);faces_Mutual(:,3);faces_Mutual(:,1)], 1);
%
%     % 如果这一批不能移除，就逐一判断，移除
%     if max(max(hedge)) >1
%         [faces_Mutual, stateP] = unRemovedV(faces_Mutual_old, hedge_old, iC,...
%             faces_delete, delete_arry, faces_add, add_arry, iC_near, stateP);
%
%         % 记得更新半边信息
%         hedge = sparse([faces_Mutual(:,1);faces_Mutual(:,2);faces_Mutual(:,3)],...
%             [faces_Mutual(:,2);faces_Mutual(:,3);faces_Mutual(:,1)], 1);
%         if max(max(hedge)) > 1
%             if max(max(hedge)) > 2; disp('半边数大于3，不可能'); end
%             % disp('半边仍然重复。');
%
%             nC = max(max(faces_Mutual));
%             rx = zeros(1, size(faces_Mutual,1));
%             for j = 1:nC
%                 [fi,~] = find(faces_Mutual == j);
%                 if length(fi)==2
%                     rx(fi) = 1;
%                 end
%             end
%             faces_Mutual = faces_Mutual(rx == 0,:);
%                     hedge = sparse([faces_Mutual(:,1);faces_Mutual(:,2);faces_Mutual(:,3)],...
%             [faces_Mutual(:,2);faces_Mutual(:,3);faces_Mutual(:,1)], 1);
%         end
%
%         if max(max(hedge)) == 2;disp('2度顶点没有去除干净。'); end
%     end
%
%     % iA中的iC点移除，更新
%     iA = find(stateP == 1);
% end
%
% % 修改面索引，去除孤立顶点
% idx_change = 1:np; % 用于修改最后的顶点索引
% removeP = find(stateP == 0);
% for j = removeP
%     idx_change(j+1:np) = idx_change(j+1:np) - 1;% adjust indexes
% end
%
% n_remove = length(removeP);
% n_remain = np - length(removeP);
%
% % adjust points
% vid = 1:np; vid(removeP) = 0;
% vertices_Mutual = [vertices(vid > 0, :); vertices_cand];
%
% % adjust faces
% idx_change2 = (np+1):(np + nCand);
% my_idx_change = [idx_change, idx_change2 - n_remove];
% faces_Mutual = my_idx_change(faces_Mutual);
% end
