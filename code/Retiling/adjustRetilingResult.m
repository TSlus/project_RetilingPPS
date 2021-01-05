function [vertices, faces, nfig] = adjustRetilingResult(vertices, faces, nfig, rp_flag)

[vertices, faces] = CleanPatch(vertices, faces);
% 1.判断网格是否为封闭网格
% 将网格调整为封闭网格
[isclose, vertices, faces] = isCloseMesh(vertices, faces);
if ~isclose
    warning('网格不封闭。');
    return;
end

% 移除三度顶点
vertices_final = vertices; faces_final = faces;
for kk = 1:6
np = size(vertices_final,1); nf = size(faces_final,1);
faces_turn = faces_final(:, [2,3,1]);
hedge = sparse(faces_final(:), faces_turn(:), [1:nf, 1:nf, 1:nf]');
%
[oneRing, val] = findNearPs(faces_final);
P3 = find(val == 3);
idxadj = 1:np; fadd = faces_final; ta = 1;
dete = zeros(1,nf);
for i = P3
    fadd(ta,:) = oneRing{i}; ta = ta+1;
    
    dete(hedge(oneRing{i}(1), oneRing{i}(2))) = 1;
    dete(hedge(oneRing{i}(2), oneRing{i}(3))) = 1;
    dete(hedge(oneRing{i}(3), oneRing{i}(1))) = 1;
    
    idxadj((i+1):np) = idxadj((i+1):np) - 1;
end
fadd = fadd(1:(ta-1), :);
faces_final = faces_final(dete == 0, :);
faces_final = idxadj([faces_final; fadd]);
%
vertices_final = removerows(vertices_final, P3);
end
% [oneRing, val] = findNearPs(faces_final); min(val)
vertices = vertices_final; faces = faces_final;


% 2.做remesh部分步骤
faces_turn = faces(:, [2,3,1]);
hedge = [faces(:), faces_turn(:)];
de = sum(abs(vertices(hedge(:, 1),:) - vertices(hedge(:, 2),:)).^2, 2).^(1/2);
L = mean(de);

vnew = vertices; fnew = faces;
% 重新计算曲率；
CSP_idx = [];
% if just do retiling CSP_idx = []
if rp_flag
FV.vertices = vertices;
FV.faces = faces;
[Cmean, Cgaussian, ~,~,~,~] = GetCurvature(FV, true);
abs_Cmean = abs(Cmean);
abs_Cgaussian = abs(Cgaussian);
np = size(vertices, 1);k_level = ceil(np*2/5);
[~, CSP_idx1] = maxk(abs_Cmean, ceil(k_level/2));
[~, CSP_idx2] = maxk(abs_Cgaussian, ceil(k_level/2));
CSP_idx = sort(unique([CSP_idx1; CSP_idx2]));
end

for i = 1:1
% [vnew, fnew, ~] = EdgeCollaps( vnew, fnew, 0.8*L ,vertices,faces);
[vertices, faces] = collapseShortEdge(vertices, faces, 0.8*L, CSP_idx);
% fnew = doFlip(fnew);
[vnew, fnew, ~] = RemoveBadTriangles( vnew, fnew, vertices, faces);
[vnew,fnew] = SubdivideLarge( vnew, fnew,1.5*L ,vertices, faces);
[vnew, fnew, ~] = RemoveBadTriangles( vnew, fnew,vertices, faces);
% fnew = doFlip(fnew);
end
vertices = vnew; faces = fnew;

% 3.对remesh后的网格再处理一次，包含移除二度，三度点
[isclose, vertices, faces] = isCloseMesh(vertices, faces);
if ~isclose
    warning('网格不封闭。');
    return;
end

% % 4.重心移动
% np = size(vertices,1); nf = size(faces, 1);
% [oneRing, val] = findNearPs(faces);
% disp(['Retiling结果最小度:', num2str(min(val))]);
% vertices2 = vertices;
% % 计算每个任意两个面之间的夹角，如果夹角度数在180度左右，就换成重心坐标。
% % 半边
% faces_turn = faces(:, [2,3,1]);
% hedge = sparse(faces(:), faces_turn(:), [1:nf, 1:nf, 1:nf]');
% % 每个面的法向量
% v1 = vertices(faces(:,2),:) - vertices(faces(:,1),:);
% v2 = vertices(faces(:,3),:) - vertices(faces(:,1),:);
% vc = cross(v1, v2, 2);
% vc_norm = sqrt(sum(vc.^2, 2));
% normF = vc ./ vc_norm;
% % 判断每个顶点邻域三角形夹角
% for i=1:np
%    nearP = oneRing{i}; n_near = length(nearP);
%    % 邻域半边
%    nearP2 = nearP([2:n_near, 1]);
%    ihedge = [nearP(:), nearP2(:)];
%    % 需要判断的三角形
%    temp = ones(n_near); temp = triu(temp, 1);
%    [xi, yi] = find(temp);
%    tri1 = ihedge(xi,:);
%    tri2 = ihedge(yi,:);
%    f1 = hedge(sub2ind(size(hedge), tri1(:,1), tri1(:,2)));
%    f2 = hedge(sub2ind(size(hedge), tri2(:,1), tri2(:,2)));
%    ni1 = normF(f1, :); ni2 = normF(f2, :);
%    % 夹角
%    theta = acos(dot(ni1, ni2, 2));
%    delta = abs(theta - pi)./pi;
% %    if min(delta) < 1e-1
%        vertices2(i, :) = sum(vertices(nearP,:),1) / n_near;
% %    end
% end
% vertices = vertices2;

% figure(nfig); nfig = nfig + 1;
% trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3));axis equal
% title('the adjust Retiling mesh')

end