function oneRingP = findNearP_sp(nearPsp, P)

% nearPsp: halfedge -- faces_name
half_face = nearPsp;

%初始化，第一个点
neighbor_P = find(half_face(:,P)); %返回非零元素位置，点P的邻域点（没有顺序）
a = neighbor_P(1);
% 存储邻域点
oneRingP(length(neighbor_P)) = 0;
% length(neighbor_P)
oneRingP(1) = a;
t = 2;
neighbor_res = neighbor_P(2:end);

%找下一个点
while  ~isempty(neighbor_res)
    a_nextP_idx = find(half_face(a, neighbor_res));
    k = 1; 
    
    b = neighbor_res(a_nextP_idx(k));
    while half_face(a,b) ~= half_face(P, a) 
        k = k + 1;
        b = neighbor_res(a_nextP_idx(k));    
    end
    oneRingP(t) = b;
    t = t + 1;
    a = b; neighbor_res = neighbor_res(neighbor_res ~= b);
end
end
