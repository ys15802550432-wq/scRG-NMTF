function W = my_construct_knn_local(X, k)
    n = size(X, 1); 
    X_norm = sqrt(sum(X.^2, 2)); X(X_norm < 1e-8, 1) = 1e-8; 
    D = pdist2(X, X, 'cosine'); D(isnan(D)) = 1; 
    [~, idx] = sort(D, 2);
    W = sparse(n, n); for i = 1:n, id = idx(i, 2:min(k+1, n)); W(i, id) = 1; end
    W = full((W + W') / 2);
end