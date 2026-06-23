function [H_star, H, W, S, rlabel, obj_history] = scrgnmf_engine(A, X, k, lambda1, lambda2, lambda3, k_dyn, idx_init)
    % scRG-NMTF Optimization Engine (with NaN-defense mechanism)
    K_views = size(A, 2); 
    n = size(A{1}, 1); 
    max_iter = 500; eta = 0.5; tol = 1e-5; my_eps = 1e-8;  
    
    % Initialization
    H_star = full(sparse(1:n, idx_init, 1, n, k)); H_star = H_star + 0.1;
    H = cell(1, K_views); W = cell(1, K_views); S = cell(1, K_views);  
    for v = 1:K_views
        d_v = size(X{v}, 2); H{v} = H_star; W{v} = rand(d_v, k); S{v} = eye(k); 
    end
    
    obj_history = zeros(max_iter, 1);
    
    for iter = 1:max_iter
        W_old = W; S_old = S; H_old = H; H_star_old = H_star;
        try
            Wnew = cell(1, K_views); Snew = cell(1, K_views); Hnew = cell(1, K_views);
            
            % Dynamic graph construction for sub-manifold alignment
            if lambda3 > 0
                S1_knn = my_construct_knn_local(H_star', k_dyn); 
                D_mat = diag(sum(S1_knn, 2)); 
            end
            
            % 1. Update W and S (Feature extractors)
            for v=1:K_views
                if lambda1 > 0
                    X_v = X{v}; 
                    term1_W = (X_v' * H_star) * S{v}'; mid = H_star' * X_v; term2_W = W{v} * S{v} * mid * W{v};
                    Wnew{v} = W{v} .* (term1_W ./ (term2_W + my_eps)).^eta; Wnew{v} = max(real(Wnew{v}), 1e-10);
                    
                    % Safe update for S (Avoids negative Laplacian values to prevent NaNs)
                    if lambda3 > 0
                        term1_S = lambda1 .* (Wnew{v}' * (X_v' * H_star)) + lambda3 .* (S1_knn * S{v});
                        term2_S = lambda1 .* S{v} + lambda3 .* (D_mat * S{v});
                    else
                        term1_S = lambda1 .* (Wnew{v}' * (X_v' * H_star));
                        term2_S = lambda1 .* S{v};
                    end
                    Snew{v} = S{v} .* (term1_S ./ (term2_S + my_eps)).^eta; Snew{v} = max(real(Snew{v}), 1e-10);
                else
                    Wnew{v} = W{v}; Snew{v} = S{v};
                end
            end
            
            % 2. Update H_v (View-specific latent features)
            for v=1:K_views
                term1_H = 2 * (A{v} * H{v}) + lambda2 .* H_star; 
                AtH = A{v} * H{v}; HtAtH = H{v}' * AtH; 
                term2_H = 2 * H{v} * HtAtH + lambda2 .* H{v} * (H_star' * H{v});
                Hnew{v} = H{v} .* (term1_H ./ (term2_H + my_eps)).^eta; Hnew{v} = max(real(Hnew{v}), 1e-10);
            end
            
            % 3. Update H_star (Consensus latent representation)
            num_Hs = zeros(size(H_star)); den_Hs = zeros(size(H_star));
            for v=1:K_views
                if lambda1 > 0
                    X_v = X{v}; 
                    num_Hs = num_Hs + (lambda1 * X_v * (Wnew{v} * Snew{v}) + lambda2 .* Hnew{v});
                    XtH = X_v' * H_star; WtXtH = Wnew{v}' * XtH; StWtXtH = Snew{v}' * WtXtH; 
                    den_Hs = den_Hs + (lambda1 * H_star * StWtXtH + lambda2 * H_star * (Hnew{v}' * H_star));
                else
                    num_Hs = num_Hs + lambda2 .* Hnew{v};
                    den_Hs = den_Hs + lambda2 * H_star * (Hnew{v}' * H_star);
                end
            end
            H_starnew = H_star .* (num_Hs ./ (den_Hs + my_eps)).^eta; H_starnew = max(real(H_starnew), 1e-10);
            
            % Check convergence
            for v=1:K_views, W{v}=Wnew{v}; S{v}=Snew{v}; H{v}=Hnew{v}; end; H_star=H_starnew;  
            
            curr_obj = norm(H_star - H_star_old, 'fro');
            if isnan(curr_obj) || isinf(curr_obj), W = W_old; S = S_old; H = H_old; H_star = H_star_old; break; end
            
            obj_history(iter) = curr_obj;
            
            if iter > 2 && curr_obj < tol, break; end
        catch
            W = W_old; S = S_old; H = H_old; H_star = H_star_old; break; 
        end
    end
    [~, rlabel] = max(H_star, [], 2);
end