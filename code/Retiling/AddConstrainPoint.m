function [CSP_idx,abs_Cmean, abs_Cgaussian] = AddConstrainPoint(mymesh, k_level)   %#ok<INUSL>
mesh_to_VFdetail; % 将mesh中的数据离散化

% 将曲率较大的点设为CSP
FV.vertices = vertices;
FV.faces = faces;
[Cmean, Cgaussian, ~,~,~,~] = GetCurvature(FV, true);
abs_Cmean = abs(Cmean);
abs_Cgaussian = abs(Cgaussian);

[~, CSP_idx1] = maxk(abs_Cmean, ceil(k_level/2));
[~, CSP_idx2] = maxk(abs_Cgaussian, ceil(k_level/2));

CSP_idx = sort(unique([CSP_idx1; CSP_idx2]));

% [~, CSP_idx] = maxk(abs_Cgaussian, ceil(k_level));
% CSP_idx = sort(CSP_idx);
% sizeCP = size(CSP_idx)

end