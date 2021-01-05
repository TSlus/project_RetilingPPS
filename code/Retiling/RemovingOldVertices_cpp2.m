function [vertices_Mutual, faces_Mutual, n_rem] = RemovingOldVertices_cpp2(...
    mymesh, vertices_cand, faces_Mutual, CSP_idx) %#ok<INUSL>
mesh_to_VFdetail;

lastwarn('');
nCand = size(vertices_cand, 1);
vertices_Mutual = [vertices; vertices_cand];
remove_idx = zeros(1, np); % equal 1 when v is removed.
idx_change = 1:np;

% faces_old = faces_faces_Mutual;
iA = 1:np; iA(CSP_idx) = 0; iA = iA(iA>0);
test_odd = 0;
for i = iA   
% for i = 1:np % be replaced.
    % if ismember(i, CSP_idx); continue; end % be replaced.
    vi = vertices(i, :);
        
    % local information.
    [nearP, nearPsp, normal, delete_fM, n_rem] = ...
        localNearP(i,vertices,faces, vertices_Mutual, faces_Mutual);
    if n_rem < 0
        return; 
    end    
    v_nearp = vertices_Mutual(nearP, :);

    % test1: check the project of 2D near p of i
    [success, xy] = checkProject2DNearP(v_nearp, vi, normal);
    if ~success; continue; end 
    
    % test2: check critical edges
    [re_flag, crirical_tri, crirical_tri3]=checkCriticalEdge(...
        nearP,faces_Mutual, nearPsp, xy, i);
    if ~re_flag; continue; end
    
    % add new faces.
    [faces_add2, re_flag, mcontinue, test_odd] = addFace(...
        xy, nearP, crirical_tri,crirical_tri3, test_odd);
    if (~re_flag || mcontinue); continue; end
    
    faces_add1 = faces_Mutual(delete_fM == 0, :);
    faces_Mutual = [faces_add1; faces_add2]; % all faces after remove i
    
    idx_change(i+1:np) = idx_change(i+1:np) - 1;% adjust indexes
    remove_idx(i) = 1;% removed vertex
end

% adjust points
remain_p = find(~remove_idx);
vertices_Mutual = [vertices(remain_p, :); vertices_cand];
n_rem = length(remain_p);

% adjust faces
idx_change2 = (np+1):(np + nCand); 
my_idx_change = [idx_change, idx_change2 - (np - n_rem)];
faces_Mutual = my_idx_change(faces_Mutual);
end

