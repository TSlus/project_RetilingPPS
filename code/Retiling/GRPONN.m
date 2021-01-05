% GenerateRandomPointOnNNumber

% 输入是一个数组a

% 产生的随机点，与 a 中的点不重合
% a 的元素个数可能比较大，比如 4*10^4
% rand_num: 表示生成的随机数值（相对于a_sum）
% idxs: 表示随机点在哪个段上

% ex:
% 0..a(1)....a(1)+a(2)..x..a(1)+a(2)+a(3).....a(1)+..+a(n)
% [rand_num, idxs] = [x, 3];

function [rand_num, idxs, a_sum_i] = GRPONN(a, nCand) % 输入a, nCand
n = length(a);

a_sum = sum(a);
a_sum_i = zeros(n+1,1);
for i = 2:n+1
    a_sum_i(i) = a_sum_i(i-1) + a(i-1);
end

rand_num = zeros(nCand, 1);
idxs = zeros(nCand, 1);
% 还是一个点一个点来找
for i = 1:nCand
    while 1
        b = rand * a_sum; 
        diff = b - a_sum_i;
        midx = find((diff) <= 0,1);
        equal_zero = find(diff(midx) == 0, 1);
        if isempty(equal_zero) % 如果有随机点与 a_sum_i 中的点重合就重新选点
            rand_num(i) = b;
            idxs(i) = midx - 1;
            break;
        end
        disp('生成的随机点在三角形边上，重新生成...');
        continue;
    end
end 
end