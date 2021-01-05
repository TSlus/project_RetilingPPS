function dis = distancePointPlane(q, p, norm_plane_p)
    row_q = size(q, 1); % 计算距离的个数
    qp = q - p;
    qp_norm = sum(abs(qp).^2,2).^(1/2);
    qp_normal = qp./ qp_norm; 
    coss = dot(qp_normal', repmat(norm_plane_p, row_q, 1)');
    dis = abs(qp_norm .* coss');
end

% distancePointPlane([1,1,0.1;2,1,0.005],[1,2,0],[0,0,1])
