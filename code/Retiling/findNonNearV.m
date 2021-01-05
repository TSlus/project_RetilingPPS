function iC = findNonNearV(iA, hedge)

% 与iA中第一个点，同一批的点 iC
iC = zeros(length(iA), 1); iC(1) = iA(1);
test_near = hedge(:, iA(1));
iB = find(~test_near); % iB之间的点也需要互不相邻
% 点集iB必须是iA的子集
% 计算 iA 交 iB
miA = max(iA);
maxm = max([miA, max(iB)]);
test1 = zeros(1, maxm); test2 = test1;
test1(iA) = 1; test2(iB) = 1;
test3 = test1 + test2;
test3(iC(1)) = 0;
iB = find(test3 == 2);

n_not_near = length(iB);
k = 1;
% while k < (ceil(n_not_near/4)+1) && ~isempty(iB) % ？？？？？？why not
% while k < (ceil(n_not_near/4)+2) && ~isempty(iB)
% while  k < (n_not_near+1) && ~isempty(iB)
while ~isempty(iB)
    iC(k+1) = iB(1);
    test_near = test_near + hedge(:,iB(1));
    iB = find(~test_near);
    
    maxm = max([miA , max(iB)]);
    test1 = zeros(1, maxm); test2 = test1;
    test1(iA) = 1; test2(iB) = 1;
    test3 = test1 + test2;
    test3(iC(1:(k+1))) = 0;
    iB = find(test3 == 2);
    
    k = k+1;
    if k > (n_not_near+1)
        disp('不能跳出循环')
    end
end

iC = iC(1:k);
% iC2 = unique(iC);
% if length(iC) ~= length(iC2)
%     disp('顶点重复移除')
% end
end