% 预计算
numP = size(vertices,1);numF = size(faces,1);
%% 邻域信息
[oneRingPs, v_valence] = findNearPs(faces); %利用元胞形式存邻域点
%% 整理三角面片，奇异点在第一个点
tri_point_valence = v_valence(faces);
[irregular_f,extra_v] = find(tri_point_valence ~= 6);
irregular_f2 = irregular_f(extra_v == 2);
irregular_f3 = irregular_f(extra_v == 3);
faces(irregular_f2, :) = faces(irregular_f2, [2,3,1]); 
faces(irregular_f3, :) = faces(irregular_f3, [3,1,2]); 

% %% 上述faces调整判断
% % 判断是否将奇异度调整到面的第一个索引
% flag1 = size(find(v_valence(faces(:,1)) ~= 6),2); %第一列度数不是6的面
% flag2 = find(v_valence(faces(:,2)) ~= 6 | v_valence(faces(:,3)) ~= 6);
% if (flag1~= length(irregular_f)) || size(flag2,2)
%     disp('三角形调整有误。')
% end

%% 非奇异面和奇异面
regular_f = setdiff(1:numF, irregular_f);      %不包含奇异点的三角形
faces_regular = faces(regular_f,:);

% 包含奇异点的面：1.度数含3的面。2.度数不含3的面。
tri_point_valence = v_valence(faces);
r3 = (tri_point_valence(irregular_f, 1) == 3);
irregular_f_3 = irregular_f(r3);
irregular_f = irregular_f(~r3);
faces_irregular = faces(irregular_f,:);

%% 构造面索引矩阵
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1]; R = [1:numF, 1:numF, 1:numF]';
vf_sparse = sparse(X, Y, R);

%% n_choose_k
max_val = max(v_valence);
n_choose_k = sparse(max_val + 2, max_val + 2);
for i = 0:max_val + 1
    for j = 0:max_val + 1
    if i >= j
        n_choose_k(i+1, j+1) = nchoosek(i,j);
    end
    end
end

%% 先计算 N=6的面的控制顶点
cr_idx = controlP_idx_regular_face(faces_regular, oneRingPs); % nf1 x 12
% 再计算N!=6的控制顶点
cir_idx = controlP_idx_irregular_face(faces_irregular, oneRingPs, v_valence);

%% 细分矩阵特征
myEigen = EIGENSTRUCT; % NMIN = 4; NMAX = 100;
