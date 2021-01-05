% 判断 push 后的candidate point 是否在原三角形内部
function flag = isOnTriangle(fnum_cand, vcand, vertices, faces, norm_face)
    v1 = vertices(faces(fnum_cand, 1),:); 
    v2 = vertices(faces(fnum_cand, 2),:); 
    v3 = vertices(faces(fnum_cand, 3),:); 
    
    norm_fi = norm_face(fnum_cand, :);
    a1 = cross(v3 - v2, vcand - v2, 2);
    a2 = cross(vcand - v1, v3 - v1, 2);
    a3 = cross(v2 - v1, vcand - v1, 2);
    temp = norm_fi;
    
    dot1 = dot(a1', temp',1); dot2 = dot(a2', temp', 1); dot3 = dot(a3', temp',1);
    flag = (dot1 > 0) & (dot2 > 0) & (dot3 > 0);
end


