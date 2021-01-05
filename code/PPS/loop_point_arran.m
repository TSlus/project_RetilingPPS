% 将坐标整理成 n x 3
% loop_point{n_regular_f} = [];
loop_point{numF} = [];
for i = 1:n_regular_f
    loop_point{regular_f(i)} = regular_patch(3*i - 2:3*i, :)';
end
for i = 1:n_irregular_f
    loop_point{irregular_f(i)} = irregular_surf{i};
end
%% 度数是3的点，在其所在面上直接取点
v = uvP(:,1); w = uvP(:,2); u = uvP(:,3);
faces_irregular_3 = faces(irregular_f_3,:);
for i = 1:length(irregular_f_3)
    loop_point{irregular_f_3(i)} = uvP * vertices(faces_irregular_3(i,:),:);
end

%% 作图
rows = 0;
L_N(numF) = 0;
for N=1:numF
    L_N(N) = size(loop_point{N},1);
    rows = rows + L_N(N);
end
vx = zeros(rows, 3);
t = 1;
for N=1:numF
    vx(t:t + L_N(N)-1,:) = loop_point{N};
    t = t + L_N(N);
end

% % 画图
% regular_patch = reshape(regular_patch(:), 3, size(regular_patch,1)*size(regular_patch,2)/3)';
% figure(1)
% % h = trimesh(faces, vertices(:,1), vertices(:,2), vertices(:,3));axis equal
% % set(h,'EdgeColor',[1.0, 0.66, 0.49],'LineWid',0.0001); hold on;
% plot3(regular_patch(:,1),regular_patch(:,2),regular_patch(:,3),'b.');hold on;
% title('regular patch')
% % legend('mesh - 3D','regular patch Loop-point')
% 
% figure(2)
% % h = trimesh(faces, vertices(:,1), vertices(:,2), vertices(:,3));axis equal
% % set(h,'EdgeColor',[1.0, 0.66, 0.49],'LineWid',0.0001); hold on;
% plot3(vx(:,1),vx(:,2),vx(:,3),'b.');
% title('Loop-point')
% % legend('mesh - 3D','Loop-point')
