function drawsphere(a,b,c,R)
% 绘制球面
% 以(a,b,c)为球心，R为半径

    % 生成数据
    [x,y,z] = sphere(20);

    % 调整半径
    x = R*x; 
    y = R*y;
    z = R*z;

    % 调整球心
    x = x+a;
    y = y+b;
    z = z+c;

    % 使用mesh绘制
%     figure;
    axis equal;
    mesh(x,y,z);

%     % 使用surf绘制
%     figure;
%     axis equal;
%     surf(x,y,z);
end