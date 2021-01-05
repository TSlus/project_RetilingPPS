function cir_idx = controlP_idx_irregular_face(faces_irregular, oneRingPs, v_valence)

n_irf = size(faces_irregular, 1);
cir_idx{n_irf} = [];

for i = 1:n_irf
    faces_N = faces_irregular(i,:);
    N = v_valence(faces_N(1));
    
    cir_N = zeros(1, N + 6);
    cir_N([1, 2, N+1]) = faces_N;
    onering_tri = oneRingPs(faces_N);
    where = find(onering_tri{1} == faces_N(3));
    path1 = onering_tri{1}([where:end, 1:where-1]);
    where = find(onering_tri{2} == faces_N(1));
    path2 = onering_tri{2}([where:end, 1:where-1]);
    where = find(onering_tri{3} == faces_N(2));
    path3 = onering_tri{3}([where:end, 1:where-1]);
    cir_N([N:-1:3, N+4, N+3, N+2, N+5, N+6]) = [path1(2:end-1), path2(3:4), path3(2:4)];
    cir_idx{i} = cir_N;
    
end

end




