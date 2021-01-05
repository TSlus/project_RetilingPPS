% compute the faces whic are less than dis to the face.
% the nearby 3 faces are choosed defaultly
function [FSD, num_FSD] = FacesShortDistance(v, f, dis, hedge_face)
if nargin < 3
    nf = size(f,1);
    v1 = v(f(:,2),:) - v(f(:,1),:);
    v2 = v(f(:,3),:) - v(f(:,1),:);
    norm_face = cross(v1, v2, 2);
    norm_face_normal = sqrt(sum((norm_face).^2, 2));
    
    si = norm_face_normal / 2; %每个小三角形面积
    sa = sum(si);              %网格表面积
    radius = 2 * sqrt(sa /(nf/4));
    dis = radius;
    
    f_turn = f(:,[2,3,1]);
    hedge_face = sparse(f(:), f_turn(:), [1:nf, 1:nf, 1:nf]');
end

nf = size(f, 1);
FSD{nf} = [];
num_FSD(nf) = 0;

% 计算每对面之间的距离
% 每个面的重心
vceng = (v(f(:,1),:) + v(f(:,2),:) + v(f(:,3),:)) / 3;

% 面与面之间的距离，只需要计算部分点
% 这部分先不改

for i = 1:nf
    a = f(i,1); b = f(i,2); c = f(i,3); 
    idxNot = [full(hedge_face(b,a)); full(hedge_face(c,b)); full(hedge_face(a,c)); i];
    dists = sum((vceng(i,:) - vceng).^2, 2);
    dists(idxNot) = 2 * dis^2;
    FSD{i} = find(dists < dis^2);
    num_FSD(i) = length(FSD{i});
end

% trimesh(f, v(:,1), v(:,2), v(:,3));
% axis equal; hold on;
% drawsphere(vceng(i,1),vceng(i,2),vceng(i,3),radius);
end