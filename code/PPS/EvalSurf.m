function irregular_surf = EvalSurf(vertices, faces, v_valence, irregular_f, cir_idx, myEigen, nEveryF)

%% Omega内参数信息
% 包含n, k, v, w
Para = Parameter_irpatch(nEveryF); % 与N无关，可以单独存起来

%% 不同 N surf 计算不同
n_ir = length(cir_idx);
irregular_surf{n_ir} = [];

% for N = [NMIN:5, 7:NMAX]
for i = 1:length(cir_idx)
    f_i = irregular_f(i);
    N = v_valence(faces(f_i,1));
    C0_irregular_N = vertices(cir_idx{i},:); % (N+6) x 3   
    % C_0_hat
    
    C_0_hat = myEigen(N).iV * C0_irregular_N; % C_0: (N+6) x 3
    
    %% 参数计算
    irregular_surf_N = zeros(3, size(Para,1));
    
    % irregular_surf_N = (C_0_hat)' * 
    
    for j = 1:size(Para,1)
        v = Para(j,1); w = Para(j,2); n = Para(j,3); k = Para(j,4);
        irregular_surf_N(:,j) = (C_0_hat)' * diag((myEigen(N).L).^(n - 1)) * ...
            myEigen(N).Phi{k}' * b_vw_div(v, w); % b_vw_div(v, w) : 12 * ni
        
    end
    
    irregular_surf{i} = irregular_surf_N'; % 每个数据为 ni x 3
end

end