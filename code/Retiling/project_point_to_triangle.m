% project_point_to_plane2D 
% 将空间点p投影到三角形所在平面
% 输入可以没有法向量
% p : n x 3
% triangle: 3 x 3
% (optional)normal_tri: 1 x 3
function ponit_on_plane = project_point_to_triangle(p, triangle, normal_tri)
    q1 = triangle(1,:);
    if(nargin < 3)
        q2 = triangle(2,:);
        q3 = triangle(3,:);
        cross_ = cross(q2 - q1, q3 -q1);
        normal_tri = cross_ / norm(cross_);
    end
    
    row_p = size(p, 1);
    pq = q1 - p; 
    pq_norm = sum(abs(pq).^2, 2).^(1/2);
    pq_normal = pq ./ pq_norm;    
    dot_ = dot(pq_normal', repmat(normal_tri, row_p, 1)');
    dot_(dot_ > 1) = 1; dot_(dot_ < -1) = -1;
    thetas = acos(dot_)';
    
    ponit_on_plane =  - (pq - pq_norm .* cos(thetas) .* normal_tri) + q1;
end

