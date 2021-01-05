% 确定12个控制顶点的索引
function cr_idx = controlP_idx_regular_face(faces_regular, oneRingPs)

n_regular_f = size(faces_regular, 1);
cr_idx = zeros(n_regular_f,12);
cr_idx(:,[4,7,8]) = faces_regular;
onering_tri = oneRingPs(faces_regular);
for k = 1:n_regular_f
    where = find(onering_tri{k,1} == faces_regular(k,3));
    path1 = onering_tri{k,1}([where:end, 1:where-1]);
    where = find(onering_tri{k,2} == faces_regular(k,1));
    path2 = onering_tri{k,2}([where:end, 1:where-1]);
    where = find(onering_tri{k,3} == faces_regular(k,2));
    path3 = onering_tri{k,3}([where:end, 1:where-1]); 
    cr_idx(k, [5, 2, 1, 3, 6, 10, 11, 12, 9]) = [path1(2:4), path2(2:4), path3(2:4)];
end
end