function myEigen = EIGENSTRUCT()
%% 特征结构
% typedef
% struct{
% double L[K]; /* eigenvalues */
% double iV[K][K]; /* inverse of the eigenvectors */
% double Phi[12][K][3]; /* Coefficients of the eigenbasis */
% } EIGENSTRUCT;
% EIGENSTRUCT eigen[NMAX];
NMIN = 4; NMAX = 100;
myEigen(NMAX).L = [];
myEigen(NMAX).iV = [];
myEigen(NMAX).Phi = [];

%%
for N = NMIN:NMAX
%% 
K = N + 6; 
M = K + 6;
%%
%      
% A = ( S     0  )     % K = N + 6    
%     ( S11  S12 )     % A: K x K      S: (N+1) x (N+1)  S12: 5 x (N + 1)    
%                      % S12: 5 x 5
% 
% A_bar = ( S      0 )     % M = K + 6
%         ( S11  S12 )     % S21: 6 x (N+1)   S22: 6 x 5
%         ( S21  S22 )
% 
% 
% 
%% 矩阵块赋值
alpha_N = 5/8 - (3 + 2*cos(2*pi/N))^2 / 64;
a_N = 1 - alpha_N; b_N = alpha_N / N;
c = 3/8; d = 1/8;

% S
s1 = zeros(N); c1 = s1(1,:);
c1(1) = c; c1(2) = d; c1(N) = d;
for j = 1:N
s1(j,:) = circshift(c1,[0,j-1]);
end
S = zeros(N+1,N+1);
S(:,1) = c; S(1,1) = a_N;
S(:,2:end) = [ones(1,N)*b_N; s1];

% S12
S12 = 1/16 * [2,0,0,0,0;
              1,1,1,0,0;
              0,0,2,0,0;
              1,0,0,1,1;
              0,0,0,0,2];

% S11
S11 = zeros(5,N+1);
S11(:,1) = [2;1;2;1;2]; 
S11(:,2) = [6;10;6;1;0];
S11(:,3) = [0;1;6;0;0];
S11(:,N) = [0;0;0;1;6];
S11(:,N+1) = [6;1;0;10;6];
S11 = S11 / 16;

% W1
W1 = [0,-1,1,0,0;
      1,-1,1,0,1;
      1,0,0,0,0;
      0,0,1,1,0;
      0,1,0,0,0];

% S21
S21 = zeros(6, N+1);
S21(:,2) = [3;3;3;1;0;0];
S21(:,3) = [0;0;1;0;0;0];
S21(:,N) = [0;0;0;0;0;1];
S21(:,N+1) = [1;0;0;3;3;3];
S21 = S21 / 8;
  
%S22
S22 = [3,1,0,0,0;
       1,3,1,0,0;
       0,1,3,0,0;
       3,0,0,1,0;
       1,0,0,3,1;
       0,0,0,1,3] / 8;

%%
mu_3_end = 3/8 + 2*cos(2*pi*(1:N-1)/N) / 8;
mu = [1, 5/8 - alpha_N, mu_3_end]; % U_0特征值
delta =[1/8, 1/8, 1/8, 1/16, 1/16]; % S12特征值

%% 计算U_0
U_0 = zeros(N+1);
U_0(:,1:2) = 1;
U_0(1,2) = -3/8 * alpha_N;
% v_k, w_k
frequncy_coe = (1:N-1)'.*(1:N-1);
v_k0 = cos(2*pi*frequncy_coe / N);
w_k0 = sin(2*pi*frequncy_coe / N);
v_k = [zeros(1,N-1); ones(1,N-1); v_k0];
w_k = [zeros(2,N-1); w_k0]; 
vw_k = reshape([v_k; w_k], N+1, 2*(N-1)); % cos/sin交替，对应相同的特征值。
flag = mod(N,2); % N是奇数，为1

% 当N时奇数时
if flag
    temp = repmat(3:(N + 3)/2,2,1);
    sigma = [mu(1), mu(2), mu(temp(:))];
    % U_0(:,3:N+1),取v_k, w_k前(N-1)/2列
    U_0(:,3:N+1) = vw_k(:,1:N-1);
    U_1 = sylvester(-S12, diag(sigma), S11 * U_0);        % 求解U_1
    error = norm(-S12*U_1 + U_1*diag(sigma) - S11 * U_0);
end

% 当N是偶数时
if ~flag
    temp = repmat(3:(N + 2)/2,2,1);
    sigma = [mu(1), mu(2), mu(temp(:)), 1/8];
    % U_0(:,3:N+1),取v_k, w_k前(N-2)/2列，和新的一列
    U_0(:,3:N) = vw_k(:,1:N-2);
    U_0(:,N+1) = [0,1,cos((1:N-1)*pi)]';
    U_1 = sylvester(-S12, diag(sigma), S11 * U_0);        % 求解U_1
    U_1(:,N+1) = [0;8;0;-8;0];
    error = norm(-S12*U_1 + U_1*diag(sigma) - S11 * U_0);
end

%% 计算U_1
% U_1 * Sigma - S12 * U_1 = S11 * U_0
% !! matlab求解矩阵方程：AX + XB = C  -->  X = sylvester(A,B,C); 验证：norm(A*X + X*B - C )
if error > 1e-10
    disp('计算U_1矩阵方程，误差大于1e-10')
end

%% A特征值、特征向量 
lambda = [sigma, delta];
% Lambda = diag(lambda);
V = [U_0, zeros(N+1,5); U_1, W1];
% V的逆
iU_0 = inv(U_0); iW1 = inv(W1);
iV = [iU_0, zeros(N+1,5); -iW1*U_1*iU_0, iW1]; %#ok<MINV>

% 验证S的特征分解
% S = U_0*diag(sigma)/(U_0)

%% 
A = [S, zeros(N+1,5);S11, S12];
A_bar = [A; S21, S22];

%% P_k : 12 x M
% !! matlab中sub2ind给矩阵特定位置赋值。

P_k1 = zeros(12,M); P_k2 = P_k1; P_k3 = P_k1;
p1_col = [3, 1, N+4, 2, N+1, N+9, N+3, N+2, N+5, N+8, N+7, N+10];
P_k1(sub2ind(size(P_k1),1:12, p1_col)) = 1;
p2_col = [N+10, N+7, N+5. N+2, N+3, N+6, N+1, 2, N+4, N, 1, 3];
P_k2(sub2ind(size(P_k2),1:12, p2_col)) = 1;
p3_col = [1, N, 2, N+1, N+6, N+3, N+2, N+5, N+12, N+7, N+10, N+11];
P_k3(sub2ind(size(P_k3),1:12, p3_col)) = 1;
% sum(P_k1 + P_k2 + P_k3);

%% 构建结构体
myEigen(N).L = lambda;    % eigenvalues
myEigen(N).iV = iV;
myEigen(N).Phi = {P_k1 * A_bar * V, P_k2 * A_bar * V, P_k3 * A_bar * V};

end

end



