% Loop曲面采样
%%
% clear; clc;
% load('model_loop2_A.mat') %加载细分两次的结果

function loop_point = LoopSurfaceSample(vertices, faces)
% 输入：一般网格细分两次的网格

pre_compute_loopSurf;
nEveryF = 10; % 每个面上采15个点
%% 1、将三角片分组：包含奇异点的三角形、不包含奇异点的三角形
irregular_f_order = sort(irregular_f);       %包含奇异点的三角形
is_isolate = irregular_f_order(1:end - 1) - irregular_f_order(2 : end);
if find(is_isolate == 0)
    disp('奇异点不孤立');
end
n_irregular_f = length(irregular_f);
n_regular_f = length(regular_f);

%% 2、先来处理没有奇异点的三角形
%确定12个顶点过后，利用样条函数直接计算

%确定12个顶点的位置
cr_idx = cr_idx';
C0_regular = vertices(cr_idx(:), :)'; % 3 x 12*rf

% 最后的表示形式
b_vw_m; % 样条函数
temp_mat = reshape(1:12*n_regular_f, 12, n_regular_f)';
C0_regular = C0_regular(:, temp_mat(:));
regular_patch = reshape(C0_regular, 3*n_regular_f, 12) * b_vw'; % C0 * b

%% 3、处理包含奇异点的三角形
% 利用一个结构体，存细分矩阵 A 的特征结构
myEigen = EIGENSTRUCT; % NMIN = 4; NMAX = 100;

% 对于存在奇异点的三角形，找到周围的控制点C_0（(N+6) x 3）
% 控制点C_0。对于不同 N，找到(N+6)个控制顶点，返回值是一个元胞 1 x NMAX
cir_idx = controlP_idx_irregular_face(faces_irregular, oneRingPs, v_valence);

% 计算有奇异点的patch
irregular_surf = EvalSurf(vertices, faces, v_valence, irregular_f, cir_idx, myEigen, nEveryF);
loop_point_arran; % 将坐标整理成 n x 3

end