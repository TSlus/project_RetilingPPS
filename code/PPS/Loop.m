% 封闭网格的loop细分
%% 细分功能
% 这里的loop细分加入一些功能
% 1）让细分多次后的网格面，与原来网格面对应
% 2）每个原始面上有哪些细分点
%% 输入
clear;clc;
load('model.mat'); faces_old = faces; vertices_old = vertices;
load('loop_point_model.mat') %用于将loop曲面对应到原始网格三角形

% load('bronze.mat'); faces_old = faces; vertices_old = vertices;
% load('loop_point_bronze.mat')

% vertices : np x 3
% faces    : nf x 3
%% mesh
nf = size(faces, 1); np = size(vertices, 1);

%% face name
% 在每次loop结束后更新这个矩阵就可以达到目的1）
% 利用目的1）找到同一个原始面的loop点，就可以实现目的2）
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1]; F = [1:nf, 1:nf, 1:nf]';
hedge_face_old = sparse(X, Y, F);
%%
trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)); axis equal
for i = 1:2
    [vertices, faces, hedge_face] = LoopSubdivide(vertices, faces, hedge_face_old);
    hedge_face_old = hedge_face;
end
trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)); axis equal

%% 1.利用得到的hedge_face，查看原始面对应loop2上哪些新点
hold on;
[row, col] = find(hedge_face == 3);
points = unique([row; col]);
hold on
plot3(vertices(points,1),vertices(points,2),vertices(points,3),'*');

%% 2.将原始网格的面和细分曲面上参数点联系起来
% 一个原始面，对应细分两次后的16个面

% 细分后的半边--面
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1];
R = [1:size(faces, 1), 1:size(faces, 1), 1:size(faces, 1)]';
he_f = sparse(X, Y, R);

connect{nf} = []; % nf 是细分前faces个数
nEvery = 100; % loop细分曲面loop_point每个三角形采点最多100个
for i=1:nf
    [row, col] = find(hedge_face == i);
    fi = he_f(sub2ind(size(he_f), row, col));
    fi = unique(fi); n_fi = length(fi);
    cell_p = loop_point(fi);
    % 将这些元胞中的点取出来
    v_fi = zeros(n_fi*nEvery, 3);
    t = 1;
    for j = 1:n_fi
        k = size(cell_p{j}, 1);
        v_fi(t:t+k-1, :) = cell_p{j};
        t = t+k;
    end
     v_fi =  v_fi(1:t-1, :);
     connect{i} = v_fi;
end
loop_point = connect;

%% 画图验证
figure(2)
faces = faces_old; vertices = vertices_old;
trimesh(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3));axis equal
hold on;
v1 = vertices(faces(1,:),:);
plot3(v1(:,1), v1(:,2), v1(:,3),'r*');
v2 = loop_point{1};
plot3(v2(:,1), v2(:,2), v2(:,3),'b.');
