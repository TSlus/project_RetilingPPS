function [vnew, fnew] = collapseShortEdge(vold, fold, targetL, CSP_idx)

vnew = vold; fnew = fold;
% targetL = targetL/sqrt(3);

fnew_turn = fnew(:,[2,3,1]);
hedges = [fnew(:), fnew_turn(:)];
eL = sqrt(sum((vnew(hedges(:,1),:) - vnew(hedges(:,2),:)).^2, 2));
shortEdges = hedges(eL<targetL,:);

while ~isempty(shortEdges)
    vnew_old = vnew; fnew_old = fnew;
    
    shortEdge = sort(shortEdges(1,:)); 
    id1 = shortEdge(1); id2 = shortEdge(2);
    member_idx = ismember(shortEdge, CSP_idx);
    
    % 点改变
    w1 = 0.5; w2 = 0.5;
    if sum(member_idx) == 1
        if member_idx(1) == 1
            w1 = 0.8; w2 = 0.2;
        else
            w1 = 0.2; w2 = 0.8;
        end
%         disp('get1')
    elseif sum(member_idx) == 2
        shortEdges = shortEdges(2:end,:);
%         disp('get2')
        continue;
    end
    meanv = w1*vnew(id1,:) + w2*vnew(id2,:);
    vnew(id1,:) = meanv; % vnew(id2,:) = meanv; 
        
    % 更新 CSP_idx 
    [isflag, ids] = ismember(id2, CSP_idx);
    if isflag %如果第二个所以是CSP_idx中元素，直接去掉就好
        CSP_idx = CSP_idx([1:(ids-1), (ids+1):end]);
    end
    np = size(vnew, 1);
    adj_idx = 1:np;
    adj_idx((id2+1):np) = adj_idx((id2+1):np)-1;
    adj_idx(id2) = id1;
    CSP_idx = adj_idx(CSP_idx);
    
    % 移除没有意义的顶点
    vnew = removerows(vnew, id2);
    
    fnew = adj_idx(fnew);
    % 移除没有意义的面
    diff1 = fnew(:,1) - fnew(:,2);
    diff2 = fnew(:,1) - fnew(:,3);
    diff3 = fnew(:,2) - fnew(:,3);
    diff = [diff1, diff2, diff3];
    
    [zx, ~] = find(diff == 0);
    zx = unique(zx);
    if length(zx)~=2
        shortEdges = shortEdges(2:end,:);
        vnew = vnew_old; fnew = fnew_old;
        continue;
    end
    fnew = removerows(fnew, zx);
        
    % 迭代条件
    fnew_turn = fnew(:,[2,3,1]);
    hedges = [fnew(:), fnew_turn(:)];
    eL = sqrt(sum((vnew(hedges(:,1),:) - vnew(hedges(:,2),:)).^2, 2)); 
    shortEdges = hedges(eL < targetL, :);
end

% 投影
vnew = Project(vnew, fnew,vold,fold);

end

