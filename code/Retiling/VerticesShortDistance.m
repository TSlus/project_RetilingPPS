% 每个面与哪些固定点“邻近”
function [VSD, num_VSD] = VerticesShortDistance(v, f, v_fix, dist)
nf = size(f, 1);
vceng = (v(f(:,1),:) + v(f(:,2),:) + v(f(:,3),:)) / 3;
VSD{nf} = [];
for i = 1:nf
    VSD{i} = find(sum((vceng(i,:) - v_fix).^2, 2) < dist.^2)';
    num_VSD = length(VSD{i});
end

end