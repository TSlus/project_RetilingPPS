% 每个面片上的控制顶点
function ctrPs = ContrlPoint_Omega2(vertices,vf_sparse,...
    v_valence, oneRingPs, loop_point, n_choose_k)

numP = size(vertices, 1);
ctrPs{numP}=[];
% eps = 1e-10;

% 根据度数，生成合适的权重矩阵 wight
max_vale = max(v_valence);
add_p_omega{max_vale} = [];
for mu = 3:max_vale
    L = abs(cos(pi/mu));
    len = 2 * ( mu + 1 );
    t = linspace(0,1,len); t = (t-0.5) * 2 * L;
    A = repmat(t,len,1);
    square = zeros(len * len, 2);
    square(:,2) = A(:);  %y坐标
    A = A'; square(:,1) = A(:);%x坐标
    % 将 square 中合适的点选出来
    
    %正多边形脚点
    pu = zeros(mu, 2);
    pu(:,1) = cos( (0:mu-1 )*( 2 * pi / mu));
    pu(:,2) = sin( (0:mu-1) * (2 * pi / mu) );
    
    xy_wight_fidx = zeros(len * len, 7); tw = 1;
    signs = cross([pu(1,:),0],[pu(2,:),0]); s = norm(signs); 
    for q = 1:size(square,1)
        x = square(q,1); y = square(q,2);
        theta = atan2(y,x);
        if theta < 0
            theta = theta + 2 * pi;
        end
        fidx = fix(theta/(2*pi/mu))+1;
        a1 = pu(fidx,:); fidxplus = fidx + 1;
        if fidx == mu; fidxplus = 1; end
        a2 = pu(fidxplus,:);
        
        xy = [x,y];
        signs1 = cross([a2-a1,0],[xy-a1,0]); c1 = sum(signs1.*signs, 2);
        signs2 = cross([xy, 0], [a2, 0]); c2 = sum(signs2.*signs, 2);
        signs3 = cross([a1, 0], [xy, 0]); c3 = sum(signs3.*signs, 2);
        % 
        if c1 > 0 && c2 > 0 && c3 > 0
            lambda = norm(signs1)/s; view = norm(signs2)/s;
            lambda(lambda > 1) = 1; lambda(lambda < 0) = 0;
            view(view > 1) = 1; view(view < 0) = 0;
            yita = 1 - lambda - view;
            yita(yita < 0) = 0;
            view = 1 - lambda - yita;
            
            xy_wight_fidx(tw,1) = x;
            xy_wight_fidx(tw,2) = y;
            xy_wight_fidx(tw,3) = lambda;
            xy_wight_fidx(tw,4) = view;
            xy_wight_fidx(tw,5) = yita;
            xy_wight_fidx(tw,6) = fidx;
            xy_wight_fidx(tw,7) = fidxplus;
            tw = tw+1;
        end
    end
    add_p_omega{mu} = xy_wight_fidx(1:tw-1, :);
end


% 计算每个点邻域上的控制点
for i = 1:numP
    mu = v_valence(i);
    m = mu+1; n = mu+1;
    xy_wight_fidx = add_p_omega{mu};
    
    n_omege = size(xy_wight_fidx, 1);
    if n_omege <(m+1) * (n+1)   %保证(M'*M)非奇异
        disp('第i和邻域Q的数量不够，增加点的选取。i=')
        disp(i);
        continue;
    end
    
    % 找到控制点 Q
    Q = zeros(n_omege, 3);
    for j = 1:n_omege
        p1_idx = oneRingPs{i}(xy_wight_fidx(j,6)); 
        p2_idx = oneRingPs{i}(xy_wight_fidx(j,7));
        
        %计算b(p)
        b1 = vertices(i,:); b2 = vertices(p1_idx,:);b3 = vertices(p2_idx,:);
        lambda = xy_wight_fidx(j,3);
        view = xy_wight_fidx(j,4);
        yita = xy_wight_fidx(j,5);
        qmesh = lambda * b1 + view * b2 + yita * b3;
        
        fi_u = vf_sparse(p1_idx, p2_idx);
        loop_point_fi = loop_point{fi_u};
        [~, pointidx] = min(sum(abs(loop_point_fi - qmesh).^2,2));
        Q(j,:) = loop_point_fi(pointidx, :);
    end
    
    cB = zeros(1, m+1);
    for ic=0:m; cB(ic+1) = n_choose_k(m+1, ic+1); end
    
    M = zeros(n_omege, (m+1)*(m+1)); % 计算最小二乘系数矩阵M
    Quv = xy_wight_fidx(:, 1:2);
    for im = 1:n_omege
        u = Quv(im,1:2); u=u'; u = u / (2 * L) + 0.5;%参数范围在(0,1)*(0,1)
        Bu = cB .*(u .^(0:m)).*((1-u).^(m-(0:m)));
        Buv = Bu(1,:)'.* Bu(2,:);
        M(im,:) = Buv(:)';
    end
    
    warning('off')
    ctrPs{i} = (M'*M)\(M'*Q);
    warning('on')
end

end



