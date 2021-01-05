function [vertices, faces] = DoRemesh(vertices, faces)

faces_turn = faces(:, [2,3,1]);
hedge = [faces(:), faces_turn(:)];
de = sum(abs(vertices(hedge(:, 1),:) - vertices(hedge(:, 2),:)).^2, 2).^(1/2);
L = mean(de);

[vertices, faces, ~, ~] = Remesher(vertices, faces, 1.1*L, 1);

% delete the vertices of valence 2
[vertices, faces] = RemoveValence2(vertices, faces);

% findNearPs_remesh(faces);
disp('done Remesh.')
end

