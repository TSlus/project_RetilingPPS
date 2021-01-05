function is_proj_flag = ProjectFormNormal(xy)

% xy是一个环上的点
n = size(xy, 1);
xy1 = xy([2:n,1], :);
xy2 = xy([3:n,1,2], :);

% 判断相邻的边是否通向
vec1 = xy1 - xy; vec1 = vec1./sqrt(sum(vec1.^2,2));
vec2 = xy2 - xy1;vec2 = vec2./sqrt(sum(vec2.^2,2));
vec_k = vec2 ./ vec1 - ones(n,2);
k = find(sum(abs(vec_k),2) < 1e-10, 1);
if ~isempty(k)
    is_proj_flag = 0;
    return;
end

% 判断不相邻的边是否相交
temp = ones(n); temp(1,n) = 0;
[idx1, idx2] = find(triu(temp,2));
p1 = xy(idx1,:); p2 = xy1(idx1,:);
w1 = xy(idx2,:); w2 = xy1(idx2,:);

x1 = p1(:,1); y1 = p1(:,2); x2 = p2(:,1); y2 = p2(:,2);
u1 = w1(:,1); v1 = w1(:,2); u2 = w2(:,1); v2 = w2(:,2);
t1_up = -(u1 - x1).*(v2 - v1) + (v1 - y1).*(u2 - u1);
t2_up =  (x2 - x1).*(v1 - y1) - (y2 - y1).*(u1 - x1);
t_down = -(x2 - x1).*(v2 - v1) + (y2 - y1).*(u2 - u1);

try
    t1 = t1_up ./ t_down;
    t2 = t2_up ./ t_down;
catch
    1
    not_pall = abs(t_down) > 1e-12; % 两个斜率靠得太近无法判断
    
    t1_up = t1_up(not_pall);
    t2_up = t2_up(not_pall);
    t_down = t_down(not_pall);
    
    t1 = t1_up ./ t_down;
    t2 = t2_up ./ t_down;
end

is_proj_flag = ~any(t1>=0 & t1<=1 & t2>=0 & t2<=1);

end