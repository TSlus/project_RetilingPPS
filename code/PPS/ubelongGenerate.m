nEveryF = 10;
nUV = ceil(sqrt(2 * nEveryF));  %ui,vi的个数
up = linspace(0.01,1-0.01,nUV);
upr = repmat(up,nUV,1);
uvPt = zeros(nUV*nUV, 2);
uvPt(:,2) = upr(:); upr=upr';
uvPt(:,1) = upr(:);
[iu, ~] = find(uvPt(:,1) + uvPt(:,2) < 1 & uvPt(:,1).* uvPt(:,2) > 0);
uvP = uvPt(iu,:); 
uvP(:,3) = 1 - uvP(:,1) - uvP(:,2);

%omega域上采点
ubelong = zeros(sum(v_valence) * size(uvP,1), 14); 
% 做点循环
ti=1;
for i=1:numP
    pathi = oneRingPs{i};
    mu = v_valence(i);
    k = length(pathi);  
    for j = 1:k       %邻域三角形循环
        uv1=[0,0];
        uv2=[cos((j-1) * 2*pi / mu), sin((j-1) * 2*pi / mu)];
        uv3=[cos(j * 2*pi / mu), sin(j * 2*pi / mu)];
        xy = uvP * [uv1;uv2;uv3];
        norm_row = sum(abs(xy).^2,2).^(1/2);
        valid_id = find(norm_row < abs(cos(pi/mu)));
        uvP = uvP(valid_id,:); % 排除了那些不在omega i 中的点
        xy = xy(valid_id,:);
        jplus = j + 1;
        jplus(jplus == k+1 ) = 1;        
        ube = zeros(length(valid_id), 8);        
        ube(:,1:3) = [i, pathi(j), pathi(jplus)] - ube(:,1:3);
        ube(:,4:6) = uvP; ube(:,7:8) = xy;
        
        %Ju
        %判断j点
        n_ji = find(oneRingPs{pathi(j)} == i);
        n_jjp = find(oneRingPs{pathi(j)} == pathi(jplus));
        mv = v_valence(pathi(j));
        u1 = [0, 0];
        u2 = [cos((n_ji-1) * 2*pi / mv), sin((n_ji-1) * 2*pi / mv)];
        u3 = [cos((n_jjp-1) * 2*pi / mv), sin((n_jjp-1) * 2*pi / mv)];
        ube(:,9:10) = uvP(:,[2,1,3]) * [u1; u2; u3];
        ube(sum(abs(ube(:,9:10)).^2,2).^(1/2) < abs(cos(pi/mv)),11) = 1;
        
        %判断j+1点jplus
        n_jpi = find(oneRingPs{pathi(jplus)} == i);
        n_jpj = find(oneRingPs{pathi(jplus)} == pathi(j));
        mv = v_valence(pathi(jplus));
        %u1 = [0, 0];
        u2 = [cos((n_jpi-1) * 2*pi / mv), sin((n_jpi-1) * 2*pi / mv)];
        u3 = [cos((n_jpj-1) * 2*pi / mv), sin((n_jpj-1) * 2*pi / mv)];
        ube(:,12:13) = uvP(:,[3,1,2]) * [u1; u2; u3];
        ube(sum(abs(ube(:,12:13)).^2,2).^(1/2) < abs(cos(pi/mv)),14) = 1;
               
        %ubelong = [ubelong; ube];
        ubelong(ti:ti + length(valid_id) - 1, :) =  ube;
        ti = ti+length(valid_id);
    end
end
ubelong = ubelong(1:ti-1, :); %预分配内存的强大！！！！！！！

% 去除重复点
% ubelong = DeleteRepPara( ubelong ,vertices, numF, vf_sparse);
