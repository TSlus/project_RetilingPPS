% ubelong retiling
% 找到 candidate points 对应的面和重心坐标
% candidate point 坐标，faces_name
% vertices_cand(1:nCand,:), nameF_cand(1:nCand)
vcd = vertices_cand(1:nCand,:); nameFcd = nameF_cand(1:nCand);
ubelongcd = zeros(nCand, 15);
ubelongcd(:,15) = 1; % candidate points 不属于三个chart中任何一个
% 1.所属面
ubelongcd(:, 1:3) = faces(nameFcd, :);

% 计算每个 candidate_point 的权值
v1 = vertices(faces(nameFcd, 1),:);
v2 = vertices(faces(nameFcd, 2),:);
v3 = vertices(faces(nameFcd, 3),:);

msi = cross(v2 - v1, v3 - v1, 2);

cs1 = cross(v3 - v2, vcd - v2, 2);
cs2 = cross(vcd - v1, v3 - v1, 2);
cs3 = cross(v2 - v1, vcd - v1, 2);

% test
cs = [cs1(:,1), cs2(:,1), cs3(:,1)]./ msi(:,1);
cstest1 = [cs1(:,2), cs2(:,2), cs3(:,2)]./ msi(:, 2) - cs;
cstest2 = [cs1(:,3), cs2(:,3), cs3(:,3)]./ msi(:, 3) - cs;
centtest = max(cstest1, cstest2);

% 2.重心坐标
ss1 = sum(abs(cs1).^2, 2).^(1/2); ss2 = sum(abs(cs2).^2, 2).^(1/2);
ss3 = sum(abs(cs3).^2, 2).^(1/2); sss = sum(abs(msi).^2, 2).^(1/2);
ubelongcd(:, 4:6) = [ss1, ss2, ss3] ./ sss;

% 3.调整 candidate vertices 所属的 chart
for i = 1:nCand
    % 先判断是否属于Omegai,否则做一次调整
    mu = v_valence(ubelongcd(i, 1));
    uvpoly = [0, 0; 1, 0; cos(2*pi / mu), sin(2*pi / mu)];
    k = 0;
    while norm(ubelongcd(i,4:6)*uvpoly) > abs(cos(pi/mu))
        ubelongcd(i,1:3) = ubelongcd(i,[2,3,1]);
        ubelongcd(i,4:6) = ubelongcd(i,[5,6,4]);
        mu = v_valence(ubelongcd(i, 1));
        uvpoly = [0, 0; 1, 0; cos(2*pi / mu), sin(2*pi / mu)];
        k = k+1;
        if k == 3
            ubelongcd(i, 15) = 0;
            disp('存在不属于任何一个chart的 c-d.');
            break;
        end
    end
end

% 4.计算每个Cadidate point 在omegaj,k上的局部坐标
omegeCV = zeros(nCand, 6);
ifIN = ones(nCand, 3); % 判断是否在omegai,j,k上
for i = 1:nCand
    if ~ubelongcd(i,15) % 第十三个指标为0
        continue;
    end
    
    for j = 1:3 % 在omegai,j,k上的局部坐标
        pathi = oneRingPs{ubelongcd(i, j)};
        % 判断第二个点位于邻域点中第几个，
        jp = j+1; if jp > 3; jp = 1; end
        n_i2 = find(pathi == ubelongcd(i,jp),1);
        % 构成的三角形
        mu = v_valence(ubelongcd(i, j));
        uvpoly = [0, 0; cos(n_i2*2*pi/mu), sin(n_i2*2*pi/mu);
            cos((n_i2+1)*2*pi/mu), sin((n_i2+1)*2*pi/mu)];
        % omegeCV
        jpp = jp+1; if jpp > 3; jpp = 1; end
        uv_cd = ubelongcd(i,[j, jp, jpp]+3) * uvpoly;
        omegeCV(i,2*(j-1)+1:2*j) = uv_cd;
        % 判断是否在omegai,j,k上
        if norm(uv_cd) > abs(cos(pi/mu))
            ifIN(i, j) = 0;
        end
    end
end

% 5.给ubelongcd赋值
ubelongcd(:,[7:10,12,13]) = omegeCV;
ubelongcd(:,[11,14]) = ifIN(:,2:3);
