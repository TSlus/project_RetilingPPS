% 去除重复点
function ubelong2 = DeleteRepPara( ubelong ,vertices, numF, vf_sparse)

v3d = vertices(ubelong(:,1),:).* ubelong(:,4) + vertices(ubelong(:,2),:).* ubelong(:,5) +...
        vertices(ubelong(:,3),:).* ubelong(:,6);
ubelong2 = ubelong;
   
face_para = zeros(size(ubelong,1),4);
face_para(:,1) = full(vf_sparse(sub2ind(size(vf_sparse), ubelong(:,1), ubelong(:,2))));
face_para(:,2:4) = v3d;
t = 1;
for i = 1:numF
    idx = face_para(:,1) == i;
    
    fpi = face_para(idx,2:4); % 其实每个点都是重复三次
    nfp = size(fpi,1);
    
    temp1 = repmat(fpi,nfp,1);
    arra = repmat(1:nfp,nfp,1);
    temp2 = fpi(arra(:),:);
    
    comp_diff = sum(abs(temp2 - temp1).^2,2);
    repPoint = find(comp_diff < 1e-4);
    inspect = ceil(repPoint/nfp);
    which = mod(repPoint,nfp);
    which(which ==0 ) = nfp;
    repP = [inspect, which];
    repP = repP(repP(:,1) ~= repP(:,2),:);
    repP = sort(repP,2);
    judge_arra = ones(1,nfp);
    judge_arra(repP(:,2)) = 0;
    
    ub = ubelong(idx,:);
    ubelong2(t:t + sum(judge_arra)-1,: ) = ub (judge_arra == 1,:);
    t = t + sum(judge_arra);
end
ubelong2 = ubelong2(1:t-1,:);

end
