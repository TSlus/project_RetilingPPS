function [vertices, faces, hedge_face] = myLoop(vertices, faces, loopTimes)
% 封闭网格的loop细分
%% 细分功能
% 这里的loop细分加入一些功能
% 1）让细分多次后的网格面，与原来网格面对应（此处没有，可有2）得到）
% 2）每个原始面上有哪些细分点

% 输入
% vertices : np x 3
% faces    : nf x 3

% 输出：
% vertices : np x 3，细分结果的点
% faces    : nf x 3，细分结果的面
% hedge_face：loop细分后的（半边--面）数据结构，
%             可以找到每个原始面上有哪些loop点，见37-38行

%% mesh
nf = size(faces, 1);

%% face name
% 在每次loop结束后更新这个矩阵就可以达到目的1）
% 利用目的1）找到同一个原始面的loop点，就可以实现目的2）
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1]; F = [1:nf, 1:nf, 1:nf]';
hedge_face_old = sparse(X, Y, F);
%%
% trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)); axis equal
for i = 1:loopTimes
    [vertices, faces, hedge_face] = LoopSubdivide(vertices, faces, hedge_face_old);
    hedge_face_old = hedge_face;
end
% h = trimesh(faces, vertices(:,1), vertices(:,2), vertices(:,3));axis equal
% set(h,'EdgeColor',[1.0, 0.66, 0.49],'LineWid',0.0001); 
% %% 利用得到的hedge_face，查看原始面对应loop上哪些新点
% hold on;
% [row, col] = find(hedge_face == 3);
% points = unique([row; col]);
% hold on
% plot3(vertices(points,1),vertices(points,2),vertices(points,3),'r*');

end