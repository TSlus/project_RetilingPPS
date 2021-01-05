function [faces_Mutual,hedge, P2] = removeValence2(faces_Mutual,hedge)

maxh = max(max(hedge));
while  maxh > 1
    P2 = fin(sum(hedge, 1) == 2); % 2度点
    rx = zeros(1, size(faces_Mutual,1));
    nP2 = length(P2);
    for j = 1:nP2
        [fi,~] = find(faces_Mutual == P2(j));
        rx(fi) = 1;
    end
    faces_Mutual = faces_Mutual(rx == 0,:);
    hedge = sparse([faces_Mutual(:,1);faces_Mutual(:,2);faces_Mutual(:,3)],...
        [faces_Mutual(:,2);faces_Mutual(:,3);faces_Mutual(:,1)], 1);
    maxh = max(max(hedge));
    
    if maxh > 2; warning('出现半边重复三次？');
end

end