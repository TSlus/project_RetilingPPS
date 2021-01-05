addpath(genpath('Remesh')); 
addpath(genpath('Retiling')); addpath(genpath('PPS'))

%% the original mesh
nfig = 1;
% plot
figure(nfig); nfig = nfig + 1;
trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3));axis equal
title('the original mesh')

%% do Re_Tiling
disp('======== Do Re-Tiling ========')
[vertices_ReT, faces_ReT, nFout, n_rem, ubelong, nfig, xdelta, CSP_idx] = ...
    Re_tiling(vertices, faces, nCand, k_level, nfig, detail_plot,method,modelname);
% note：vertices_ReT中前n_rem个点是原始顶点
vertices_ReT_old = vertices_ReT; faces_ReT_old = faces_ReT;
% retiling 调整
[vertices_ReT2, faces_ReT2, nfig] = adjustRetilingResult(...
    vertices_ReT_old, faces_ReT_old, nfig, 0);

%% do PPS
% 1.loopSurface
disp('======== Do PPS ========')
disp('Constructing LoopSurface...')
loop_point = mesh_connect_LoopSurf(vertices, faces);

disp('Computing cadidate-PPS...')
vcdPPS = zeros(size(ubelong, 1), 3);

% generate ubelong to candidate points
idx_odd = ubelong(:, 15) == 0; % 有部分点可能不被Omega区域覆盖
a = ubelong(idx_odd, 1); vap = vertices(a,:);
b = ubelong(idx_odd, 2); vbp = vertices(b,:);
c = ubelong(idx_odd, 3); vcp = vertices(c,:);
wights = ubelong(idx_odd, 4:6);
vcdPPS(idx_odd, :) = wights(:,1).*vap + wights(:,2).*vbp + wights(:,3).*vcp;

idx_ord = ubelong(:, 15) == 1;
ube_part = ubelong(idx_ord,:);

% 2.do Retiling and PPS
vcdPPS(idx_ord, :) = surfaceConstruction(vertices, faces, loop_point, ube_part);
vertices_ReT(n_rem+1:end, :) = vcdPPS;
vertices_final = vertices_ReT; faces_final = faces_ReT;

%% the result of Retiling & PPS
[vrp, frp] = rusultRetilingPPS(vertices, faces,vertices_final,faces_final, vertices_ReT_old,n_rem, nFout);

figure(99); 
trimesh(frp, vrp(:,1),vrp(:,2),vrp(:,3));axis equal
title('Retiling & PPS');

