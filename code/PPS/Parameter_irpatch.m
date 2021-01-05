% 选参数，定n, k, v, w
function Para = Parameter_irpatch(nEveryF)
% nEveryF = 10; % 每个奇异面上采120个点
% nEveryF = 500; % 每个奇异面上采120个点
nUV = ceil(sqrt(2 * nEveryF));  %ui,vi的个数
up = linspace(0+0.01,1-0.01,nUV);
upr = repmat(up,nUV,1);
uvPt(:,2) = upr(:);upr = upr';
uvPt(:,1) = upr(:);
[iu, ~] = find(uvPt(:,1) + uvPt(:,2) <= 1);
uvP = uvPt(iu,:); uvP(:,3) = 1 - uvP(:,1) - uvP(:,2);
eps = 1e-20;
uvP(uvP < eps) = 0;
uxy = [0,0]; vxy = [1,0]; wxy = [0,1];

uvP2 = uvP * [uxy; vxy; wxy];
uvP2 = uvP2(uvP2(:,1) + uvP2(:,2) > 1/(2^8),:); % LOOP细分次数一般不大于10次假设
Para = zeros(size(uvP2,1), 4);
Para(:,[1,2]) = uvP2;
% plot(Para(:,1),Para(:,2),'*')

% n
Para(:,3) = floor(1-log2(Para(:,1) + Para(:,2)));
% k
Para(:,[1,2]) = Para(:,[1,2]) .* 2.^ (Para(:,3)-1);
k_1 = find(Para(:,1) > 0.5);
k_2 = find(Para(:,2) > 0.5);
k_3 = find(Para(:,1) <= 0.5 & Para(:,2) <= 0.5);

Para(k_1, [1,2,4]) = [2*Para(k_1,1)-1, 2*Para(k_1,2),ones(length(k_1),1)];
Para(k_2, [1,2,4]) = [2*Para(k_2,1), 2*Para(k_2,2)-1,ones(length(k_2),1)*3];
Para(k_3, [1,2,4]) = [1 - 2*Para(k_3,1), 1 - 2*Para(k_3,2),ones(length(k_3),1)*2];
% sort([k_1; k_2; k_3]

end


