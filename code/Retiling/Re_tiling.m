function [vertices_ReT, faces_ReT, nFout, n_rem, ubelongcd, nfig, xdelta, CSP_idx] = ...
    Re_tiling(vertices, faces, nCand, k_level, nfig, detail_plot, method,modelname)
%%
pre_compute_ReT;

%% 1.generateVertices and constrain points.
disp('Generate random Candidate Vertices on surface...')
[nameF_cand, vertices_cand, CSP_idx] = generateVertices2(mymesh, nCand, k_level);
nCand = length(nameF_cand);

% CSP
vertices_CSP = vertices(CSP_idx, :); n_CSP = length(CSP_idx);

% plot v_cand and CSP
if detail_plot
figure(nfig); nfig = nfig + 1;
trimesh(faces, vertices(:,1), vertices(:,2), vertices(:,3));axis equal;
hold on
plot3(vertices_cand(:,1), vertices_cand(:,2), vertices_cand(:,3), 'b*');
plot3(vertices(CSP_idx,1), vertices(CSP_idx,2), vertices(CSP_idx,3), 'r*');axis equal
hold off; title('Candidate Points and CSP')
end

%% prapare for loop
vertices_cand = [vertices_cand; vertices_CSP];
nameF_cand = [nameF_cand, zeros(1, n_CSP)];

% compute the distance of faces to faces and faces to CSP
disp('computing FSD and VSD...')
% [FSD, num_FSD] = FacesShortDistance(vertices, faces, 1.2*radius, hedge_face);
[FSD, num_FSD] = FacesShortDistance(vertices, faces, 1.8*radius, hedge_face);
disp(['the minius number of num_FSD:', num2str(min(num_FSD))]);
[VSD, num_VSD] = VerticesShortDistance(vertices, faces, vertices_CSP, 1.8*radius);
disp(['the minius number of num_vSD:', num2str(min(num_VSD))]);

%% 2.loop k 
disp('Move candidate points on mesh, looping...')
for n_move = 1:30
    % (1).repulsive force
    forces = ComputeRepulsiveForce_adj(mymesh, nameF_cand, vertices_cand, FSD, VSD);
    nameF_cand3 = nameF_cand; % for forces test(test 1)

    % (2).move candidate points
    % we neednt consider the CSP when moving on mesh
    [vertices_cand2, nameF_cand2] = moveOnMesh_PLUS(mymesh,  vertices_cand(1:nCand, :), ...
        nameF_cand(1:nCand), forces(1:nCand, :));
    vertices_cand(1:nCand, :) = vertices_cand2;
    nameF_cand(1:nCand) = nameF_cand2;
end

% force_test; % my test 1, to check the force perpendicular to the face normal
aplane = sum(forces(1:nCand, :) .* norm_face(nameF_cand3(1:nCand), :),2);
xdelta = max(abs(aplane));

% plot, after moving candidate points
if detail_plot
figure(nfig); nfig = nfig + 1;
trimesh(faces, vertices(:,1), vertices(:,2), vertices(:,3));axis equal;
hold on
plot3(vertices_cand(:,1), vertices_cand(:,2), vertices_cand(:,3), 'b*');
plot3(vertices(CSP_idx,1), vertices(CSP_idx,2), vertices(CSP_idx,3), 'ro');axis equal
hold off
title('Candidate Points and CSP after moving')
end

%% UbelongCD
UbelongCD; % Re-Tiling to PPS

%% 3.doMutualTesselation
disp('Do MutualTesselation...')
% try
% save doMutualData mymesh nameF_cand nCand vertices_cand
faces_Mutual = doMutualTesselation2(mymesh, nameF_cand(1:nCand), vertices_cand(1:nCand,:));
% catch 
%     save doMutualData mymesh nameF_cand nCand vertices_cand
%     disp('数据已经保存。')
% end
nFout = nameF_cand(1:nCand);

vertices_Mutual = [vertices; vertices_cand(1:nCand,:)];
% save  vertices_Mutual.txt  vertices_Mutual -ascii 
% save  faces_Mutual.txt  faces_Mutual -ascii 
% save  CSP_idx.txt  CSP_idx -ascii 

mesh_test; % my test 2, check the candidate points in triangles or not.

% plot, MutualTesselation
if detail_plot
figure(nfig); nfig = nfig + 1;
trimesh(faces_Mutual, vertices_Mutual(:,1), vertices_Mutual(:,2), vertices_Mutual(:,3));
axis equal; title('Mutual Tesselation');
end

%% save data
% save dataBronze vertices_cand faces_Mutual CSP_idx nCand ubelongcd;
% if strcmp(modelname,'bronze.mat')
% mdata2 = load('dataBronze.mat');
% vertices_cand = mdata2.vertices_cand;
% faces_Mutual = mdata2.faces_Mutual;
% CSP_idx = mdata2.CSP_idx;
% nCand = mdata2.nCand;
% ubelongcd = mdata2.ubelongcd;
% 232
% end

%% 4.RemovingOldVertices
disp('Removing Old Vertices...')
% method 1
if method == 1 || method == 4
 tic
[vertices_ReT, faces_ReT, n_rem] = RemovingOldVertices_cpp2(...
    mymesh, vertices_cand(1:nCand, :), faces_Mutual, CSP_idx);
mt1 = toc;
disp(['移除顶点时，方法1用时：',num2str(mt1),'s'])
if n_rem<0; return; end
% plot 
if detail_plot
figure(nfig); nfig = nfig + 1;
trimesh(faces_ReT, vertices_ReT(:,1), vertices_ReT(:,2), vertices_ReT(:,3));axis equal;
title('mesh after kicking old points by method1');
end
end

% method 2
if (method == 2 || method == 4)
tic
[vertices_ReT, faces_ReT, n_rem] = removing_old_vertices2(...
    mymesh, vertices_cand(1:nCand, :), faces_Mutual, CSP_idx);
mt2 = toc;
disp(['移除顶点时，方法2用时：',num2str(mt2),'s'])
% plot 
if detail_plot
figure(nfig); nfig = nfig + 1;
trimesh(faces_ReT, vertices_ReT(:,1), vertices_ReT(:,2), vertices_ReT(:,3));axis equal;
title('mesh after kicking old points by method2');
end
end

% method 3
if (method == 3 || method == 4)
tic
try
[vertices_ReT, faces_ReT, n_rem, P2, CSP_idx] = RemoveOldVertices_new(...
    mymesh, vertices_cand(1:nCand, :), faces_Mutual, CSP_idx);
if ~isempty(P2)
    disp('P2非空！');
end

catch
    save dataBronze_new2 vertices_cand faces_Mutual CSP_idx nCand faces vertices;
    disp('移除顶点出错。')
end

mt3 = toc;
disp(['移除顶点时，方法3用时：',num2str(mt3),'s'])

% plot 
if detail_plot
figure(nfig); nfig = nfig + 1;
trimesh(faces_ReT, vertices_ReT(:,1), vertices_ReT(:,2), vertices_ReT(:,3));axis equal;
title('mesh after kicking old points by method3');
end
end

end