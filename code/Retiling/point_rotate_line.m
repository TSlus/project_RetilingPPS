function rotated_position  = point_rotate_line(p3, line, v0, thetas, n1, n2, n3)
    % p3 -- the point to be rotated -- (n1+n2+n3) x 3
    % line -- Direction vector of rotation axis
    % v0 -- the point on the axis
    % thetas -- the angle for rotate, column vector
    
    if(nargin > 4)
        temp = [ones(1,n1), 2*ones(1,n2), 3*ones(1,n3)];
        line = line(temp, :);
        v0 = v0(temp, :);
        thetas = thetas(temp,:);
    end
    
    % 公式计算
    p3_bar = p3 - v0; % 将带旋转的点平移，再做，因为公式只对过原点的轴正确
    
    xrot =  (line(:,1).*line(:,1).*(1-cos(thetas)) + cos(thetas)) .* p3_bar(:,1) +...
            (line(:,1).*line(:,2).*(1-cos(thetas)) - line(:,3).*sin(thetas)) .* p3_bar(:,2) +...
            (line(:,1).*line(:,3).*(1-cos(thetas)) + line(:,2).*sin(thetas)) .* p3_bar(:,3);
        
    yrot =  (line(:,1).*line(:,2).*(1-cos(thetas)) + line(:,3).*sin(thetas)) .* p3_bar(:,1) +...
            (line(:,2).*line(:,2).*(1-cos(thetas)) + cos(thetas)) .* p3_bar(:,2) +...
            (line(:,2).*line(:,3).*(1-cos(thetas)) - line(:,1).*sin(thetas)) .* p3_bar(:,3);
        
    zrot =  (line(:,1).*line(:,3).*(1-cos(thetas)) - line(:,2).*sin(thetas)) .* p3_bar(:,1) +...
            (line(:,2).*line(:,3).*(1-cos(thetas)) + line(:,1).*sin(thetas)) .* p3_bar(:,2) +...
            (line(:,3).*line(:,3).*(1-cos(thetas)) + cos(thetas)) .* p3_bar(:,3);
    p3_bar_rot = [xrot, yrot, zrot];
    rotated_position = p3_bar_rot + v0; % p3旋转到 p 所在的平面
    
end