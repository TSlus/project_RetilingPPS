function [success,xy] = checkProject2DNearP(v_nearp, vi, normal)
n_near = size(v_nearp,1);
normal0 = normal;
n_test = 0; success = 0;
while n_test < 13
    switch n_test
        case 0
        case 1; normal = [.8507, .4472, .2764];
        case 2; normal = [-.8507, .4472, .2764];
        case 3; normal = [.8507, -.4472, .2764];
        case 4; normal = [.8507, .4472, -.2764];
        case 5; normal = [.5257, .4472, .7236];
        case 6; normal = [-.5257, .4472, .7236];
        case 7; normal = [.5257, -.4472, .7236];
        case 8; normal = [.5257, .4472, -.7236];
        case 9; normal = [.0, .4472, .8944];
        case 10; normal = [.0, -.4472, .8944];
        case 11; normal = [.0, .4472, -.8944];
        case 12; normal = [0, 1, 0];
        %case 13; normal = [1, 0,0];
        %case 14; normal = [0, 0, 1];
    end
    if sum(normal.*normal0) < 0; normal = -normal; end
    
    base1 = project_point_to_triangle3(v_nearp(1, :), vi, normal);
    base1 = base1/norm(base1);
    base2 = cross(normal, base1);
    
    % 先投影，再计算二维坐标对应
    proj_nears = project_point_to_triangle3(v_nearp, vi, normal);
    xy = zeros(n_near, 2);
    xy(:, 1) = sum(proj_nears.*base1, 2);
    xy(:, 2) = sum(proj_nears.*base2, 2);
    
    % 判断投影方向是否可用
    is_proj_flag = ProjectFormNormal(xy);
    if is_proj_flag
        success = 1;
        return;
    end
    n_test = n_test+1;
end
end