% 预计算
numP = size(vertices,1);numF = size(faces,1);
%% 邻域信息
[oneRingPs, v_valence] = findNearPs(faces); %利用元胞形式存邻域点

%% 构造面索引矩阵
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1]; R = [1:numF, 1:numF, 1:numF]';
vf_sparse = sparse(X, Y, R);

%% n_choose_k
max_val = max(v_valence);
n_choose_k = zeros(max_val + 2, max_val + 2);
for i = 0:max_val + 1
    for j = 0:max_val + 1
    if i >= j
        n_choose_k(i+1, j+1) = nchoosek(i,j);
    end
    end
end




