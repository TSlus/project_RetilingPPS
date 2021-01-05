function [nearP, nearPsp, normal, delete_fM, n_rem] = localNearP(i,vertices, ...
    faces, vertices_Mutual, faces_Mutual)

n_rem = 0;
nearP = []; normal = [];
% Vertex local information.
vi = vertices(i, :);
[r_i, ~] = find(faces_Mutual == i);
delete_fM = zeros(1, size(faces_Mutual, 1));
delete_fM(r_i) = 1;
tri_1 = faces_Mutual(r_i,:);
tri_turn = tri_1(:, [2,3,1]);

% nearPsp :for critical edge.
nf_sp = size(tri_1, 1);
nearPsp = sparse(tri_1(:), tri_turn(:), [1:nf_sp, 1:nf_sp, 1:nf_sp]');
try
    nearP = findNearP_sp(nearPsp, i);
catch
    figure(100)
    trimesh(faces_Mutual, vertices_Mutual(:,1), vertices_Mutual(:,2), vertices_Mutual(:,3));
    axis equal; hold on
    plot3(vertices(i,1),vertices(i,2),vertices(i,3),'r*');
    
    figure(101)
    trimesh(faces, vertices(:,1), vertices(:,2), vertices(:,3));
    axis equal; hold on
    plot3(vertices(i,1),vertices(i,2),vertices(i,3),'r*');
    % warning('check valence 2.')
    % 查看奇异点半边情况
    n_rem = -1;
    % save odddata faces_Mutual vertices_Mutual
    return;
end

n_near = length(nearP);

% compute the normal of i vertex
v_nearp = vertices_Mutual(nearP, :);
normals_pi = cross(v_nearp - vi, v_nearp([2:n_near, 1], :) - vi, 2);
normals_pi_norm = sqrt(sum(normals_pi.^2, 2));
normal = sum(normals_pi ./ normals_pi_norm, 1);
normal = normal./norm(normal);

end