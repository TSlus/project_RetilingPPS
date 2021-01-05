%function [CanF_Rotate, CDP_rotate] = CandidateAfterRotate3(vc, v, f, nameF, norm_face, hedge_face)
% 每个candidate point 旋转三次
function [CanF_Rotate, CDP_rotate] = CandidateAfterRotate3(vc, v, f, nameF, norm_face, hedge_face)

nc = size(vc, 1);

tempp = repmat(1:nc, 3, 1);
% 旋转点
points = vc(tempp(:), :);

faces_cand = f(nameF, :);
% 旋转轴
v0 = v(faces_cand(:,1),:);
v1 = v(faces_cand(:,2),:);
v2 = v(faces_cand(:,3),:);

line1 = v1 - v0;
line2 = v2 - v1;
line3 = v0 - v2;

lines = zeros(3*nc, 3);
temp = reshape(1:3*nc, 3, nc);
lines(temp(1,:), :) = line1;
lines(temp(2,:), :) = line2;
lines(temp(3,:), :) = line3;
clear line1 line2 line3;

% 夹角
norm0 = norm_face(nameF, :);
f1 = hedge_face(sub2ind(size(hedge_face), faces_cand(:,2), faces_cand(:,1)));
norm1 = norm_face(f1, :);
f2 = hedge_face(sub2ind(size(hedge_face), faces_cand(:,3), faces_cand(:,2)));
norm2 = norm_face(f2, :);
f3 = hedge_face(sub2ind(size(hedge_face), faces_cand(:,1), faces_cand(:,3)));
norm3 = norm_face(f3, :);

norm_4 = zeros(3*nc, 3);
norm_4(temp(1,:), :) = norm1;
norm_4(temp(2,:), :) = norm2;
norm_4(temp(3,:), :) = norm3;

theta1 = sum(norm0.*norm1, 2);
theta2 = sum(norm0.*norm2, 2);
theta3 = sum(norm0.*norm3, 2);

clear norm0 norm1 norm2 norm3;

thetas = zeros(3*nc, 1);
thetas(temp(1,:), 1) = theta1;
thetas(temp(2,:), 1) = theta2;
thetas(temp(3,:), 1) = theta3;

thetas(thetas > 1) = 1;
thetas(thetas < -1) = -1;
thetas = acos(thetas);

% 转轴上一点
faces_candT = faces_cand';
vline = v(faces_candT(:), :);

% 公式计算
p_rot1 = point_rotate_line(points, lines, vline, thetas);
p_rot2 = point_rotate_line(points, lines, vline, -thetas);

clear points lines thetas;

% 判断旋转点
dist1 = sum((p_rot1 - vline) .* norm_4, 2);
dist1 = abs(dist1);

dist2 = sum((p_rot2 - vline) .* norm_4, 2);
dist2 = abs(dist2);

diff_dis = dist1 - dist2;
diff_neg = diff_dis < 0;
diff_pos = ~diff_neg;

% 最终旋转点
CDP_rotate = zeros(3*nc, 3);
CDP_rotate(diff_neg, :) = p_rot1(diff_neg, :);
CDP_rotate(diff_pos, :) = p_rot2(diff_pos, :);

CDP_rotate = project_point_to_triangle3(CDP_rotate, vline, norm_4);

% 测试旋转过后，每个点到目标平面的距离
% dis = sum((CDP_rotate - vline).* norm_4, 2);
% max(abs(dis))

% 构造稀疏索引矩阵
X = repmat(1:nc, 3, 1); X = X(:);
Y = faces_candT(:);
CanF_Rotate = sparse(X, Y , 1:3*nc);

% trimesh(f, v(:,1), v(:,2), v(:,3));
% hold on
% plot3(CDP_rotate(:,1), CDP_rotate(:,2), CDP_rotate(:,3), 'r.')
% axis equal

end