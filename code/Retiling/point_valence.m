function vertex_valence = point_valence(faces, np)

x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1];
half_face = sparse(X, Y, 1);

% np = max(max(faces));
vertex_valence(np) = 0;
for P = 1:np
    neighbor_P = find(half_face(:,P));
    vertex_valence(P) = length(neighbor_P);
end

end