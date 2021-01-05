function [re_flag, crirical_tri, crirical_tri3] = checkCriticalEdge(...
    nearP,faces_Mutual, nearPsp, xy, i)

n_near = length(nearP);
crirical_tri = zeros(n_near, 4); t_ce = 1;
crirical_tri3 = zeros(ceil(n_near^2/2)+1, 2); t_3 = 1;

re_flag = 1; loopf1 = 0;
for ip = nearP
    [rip, ~] = find(faces_Mutual == ip);
    triRing = faces_Mutual(rip,:);
    
    for j = 1:size(triRing, 1)
%         a = triRing(j, 1);
%         idx1 = find(nearP == a, 1); if isempty(idx1); continue; end
%         b = triRing(j, 2);
%         idx2 = find(nearP == b, 1); if isempty(idx2); continue; end
%         c = triRing(j, 3);
%         idx3 = find(nearP == c, 1); if isempty(idx3); continue; end
        
        test_arry = [-1,-1,-1];
        mid1 = 0; mid2 = 0; mid3 = 0;
        a = triRing(j, 1); b = triRing(j, 2); c = triRing(j, 3);
        idx1 = find(nearP == a, 1); 
        if isempty(idx1); test_arry = [a,b,c]; mid1 = 1; end
        idx2 = find(nearP == b, 1); 
        if isempty(idx2); test_arry = [b,c,a]; mid2 = 1; end
        idx3 = find(nearP == c, 1); 
        if isempty(idx3); test_arry = [c,a,b]; mid3 = 1; end
        if (mid1 + mid2 + mid3) > 1; continue; end
        
        % 半边信息，用于判定半边是否重复（1）
        if (mid1 + mid2 + mid3) == 1 && test_arry(1)~= i
        	ex_rep_3 = sum((crirical_tri3(1:t_3, :) - test_arry(2:3)).^2, 2);
            % min(ex_rep_3)
            if min(ex_rep_3) ~= 0
                crirical_tri3(t_3,:) =  test_arry(2:3); t_3 = t_3+1;
            end
        end
        
       % 保证crirical_tri不重复
        if (mid1 + mid2 + mid3) > 0; continue; end
        ex_rep1 = sum(([a,b,c] - crirical_tri(1:t_ce, 1:3)).^2, 2);
        ex_rep2 = sum(([b,c,a] - crirical_tri(1:t_ce, 1:3)).^2, 2);
        ex_rep3 = sum(([c,a,b] - crirical_tri(1:t_ce, 1:3)).^2, 2);
        if min([ex_rep1; ex_rep2; ex_rep3]) == 0; continue; end
        
        % 判定三角化是否合理
        aif = nearPsp(c, b); bif = nearPsp(a, c); cif = nearPsp(b, a);
        ceif = find(~[aif, bif, cif]);
        if length(ceif) ~= 1 % 如果三条边都不是critical edge
            continue;
        end
        
        % 半边信息，用于判定半边是否重复（2）
        abc_order = [a, b, c; b, c, a; c, a, b];
        myabc = abc_order(ceif, :);
        crirical_tri(t_ce, 1:3) = myabc;
        crirical_tri(t_ce, 4) = myabc(1); t_ce = t_ce+1;
        
        idxs = [idx3, idx2, idx1; idx1, idx3, idx2; idx2, idx1, idx3];
        idxs = idxs(ceif, :); % c, b, a 在邻域点钟的位置
        e1 = xy(idxs(2), :) - xy(idxs(3), :);
        e2 = xy(idxs(1), :) - xy(idxs(3), :);
        mcross = cross([e1, 0], [e2, 0]);
        if mcross(3) < 0 % 点1不可以去除，不做三角化
            re_flag = 0; loopf1 = 1;break;
        end
        break; % 每个邻域点只会有一种 critical顶点 情形。
    end
    if loopf1 ; break; end
end

crirical_tri = crirical_tri(1:(t_ce-1),:);
crirical_tri3 = crirical_tri3(1:(t_3-1),:);
end