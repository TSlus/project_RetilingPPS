% p : n x 3
% triangle: n x 3
% (optional)normal_tri: n x 3
function ponit_on_plane = project_point_to_triangle3(p, vplane, normal_tri)
    d_plane = sum(vplane .* normal_tri, 2);
    d_p = sum(p .* normal_tri, 2);
    ponit_on_plane = p - (d_p - d_plane) .* normal_tri;
end

