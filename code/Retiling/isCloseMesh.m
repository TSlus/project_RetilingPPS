function [isclose, vertices, faces] = isCloseMesh(vertices, faces)

isclose = 1;

% 1.判断网格是否封闭，可以包含二度顶点
np = size(vertices, 1);
% 检查顶点
mtest1 = sort(unique(faces(:))) - (1:np)'; mtest1 = max(abs(mtest1));
% 检查半边
faces_turn = faces(:,[2,3,1]);
hedge = sparse(faces(:), faces_turn, 1);
mtest2 = find(hedge - hedge', 1);
if max(abs(mtest1)) == 0 && isempty(mtest2)
    % disp('网格封闭!');
else
    warning('网格不封闭.');
    isclose = 0;
    return;
end

if max(max(hedge)) > 1
    disp('存在半边重复，度二点，已经移除')
    [vertices, faces] = RemoveValence2(vertices, faces); % 移除二度点
end

try
    [oneRing, val] = findNearPs(faces);
catch
    warning('网格邻域信息出错.');
end

end