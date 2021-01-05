function [vrp, frp] = rusultRetilingPPS(vertices, faces, vertices_final, faces_final,vertices_ReT_old,...
    n_rem, nFout)

fL1 = sqrt(sum((vertices(faces(:,1),:) - vertices(faces(:,2),:)).^2, 2));
fL2 = sqrt(sum((vertices(faces(:,2),:) - vertices(faces(:,3),:)).^2, 2));
fL3 = sqrt(sum((vertices(faces(:,1),:) - vertices(faces(:,3),:)).^2, 2));

fL = min([fL1 , fL2 , fL3],[], 2);

mnp = size(vertices_ReT_old,1);
for i = (n_rem+1):mnp
    fi = nFout(i-n_rem);
    diff = norm(vertices_ReT_old(i,:) - vertices_final(i,:));
    if diff > fL(fi)/10
        vertices_final(i,:) = vertices_ReT_old(i,:);
    end
end

nfig = 100;
[vrp, frp] = adjustRetilingResult(...
    vertices_final, faces_final, nfig, 1);

end