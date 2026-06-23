% =========================================================================
% Demo script for scRG-NMTF 
% Task: Integrative Single-Cell Clustering
% Dataset: PBMC-10k (10x Genomics, RNA + ATAC)
% =========================================================================

clear; clc; close all;
warning('off', 'all');

% 添加子文件夹路径（重要！这样才能调用核心算法和评估指标）
addpath('./src');
addpath('./utils');
addpath('./data'); % 假设你的数据放在 data 文件夹下

%% Configuration
n_runs = 5;                 % Number of independent runs
dataset_name = 'PBMC-10k';
k_nn_static = 20;           % K for static SNN graph
k_nn_dynamic = 1;           % K for dynamic cluster-manifold (Prevents over-smoothing)

% Optimal Hyper-parameters for PBMC-10k
lambda1 = 0.01;             % Feature reconstruction weight
lambda2 = 5.0;              % Multi-view consensus weight
lambda3 = 0.001;            % Dynamic graph regularization weight

fprintf('========================================================\n');
fprintf('  Running scRG-NMTF on %s \n', dataset_name);
fprintf('========================================================\n');

%% 1. Load Data
disp('[1/5] Loading multi-omics data...');
try
    data = load('PBMC-10k_GoldStandard.mat');
catch
    error('Dataset not found. Please ensure PBMC-10k_GoldStandard.mat is in the data/ path.');
end

if isfield(data, 'X_rna'), X_rna = double(data.X_rna); else, X_rna = double(data.X1); end
if isfield(data, 'X_atac'), X_atac = double(data.X_atac); else, X_atac = double(data.X2); end
if isfield(data, 'true_labels'), true_labels = double(data.true_labels); else, true_labels = double(data.Y); end
[~, ~, true_labels] = unique(true_labels);
n_samples = length(true_labels);
k_clusters = length(unique(true_labels));
fprintf('      Cells: %d | Modalities: RNA, ATAC | Cell Types: %d\n', n_samples, k_clusters);

%% 2. Preprocessing & ReLU-Mirrored Feature Expansion
disp('[2/5] Performing dimensionality reduction and ReLU expansion...');
if size(X_rna, 2) > 1000
    X_rna_norm = zscore(log(1 + X_rna ./ (sum(X_rna, 2) + eps) * 1e4), 0, 1); 
    X_rna_norm(isnan(X_rna_norm) | isinf(X_rna_norm)) = 0; 
    [~, X_rna_pca] = pca(full(X_rna_norm), 'NumComponents', 50);
    
    tf = log(1 + X_atac ./ (sum(X_atac, 2) + eps) * 1e4); 
    idf = log(1 + size(X_atac, 1) ./ (sum(X_atac > 0, 1) + eps));
    X_atac_tfidf = tf * spdiags(idf', 0, length(idf), length(idf));
    X_atac_tfidf(isnan(X_atac_tfidf) | isinf(X_atac_tfidf)) = 0; 
    [U, S_mat, ~] = svds(X_atac_tfidf, 50); 
    X_atac_lsi = U * S_mat; X_atac_lsi = X_atac_lsi(:, 2:end); 
else
    X_rna_pca = full(X_rna); X_atac_lsi = full(X_atac);
end

% NaN-defense mechanism
X_rna_pca = real(full(X_rna_pca)); 
X_atac_lsi = real(full(X_atac_lsi));
X_rna_pca(isnan(X_rna_pca) | isinf(X_rna_pca)) = 0;
X_atac_lsi(isnan(X_atac_lsi) | isinf(X_atac_lsi)) = 0;

% Zero-loss Non-negative Transformation (ReLU Mirroring)
X_base{1} = normalize([max(X_rna_pca, 0), max(-X_rna_pca, 0)], 2, 'norm');
X_base{2} = normalize([max(X_atac_lsi, 0), max(-X_atac_lsi, 0)], 2, 'norm');

%% 3. Static Graph Construction (Prior Manifold)
disp('[3/5] Building high-purity static SNN graphs...');
function A_snn = build_snn_graph(X, k)
    n = size(X, 1);
    [idx, ~] = knnsearch(X, X, 'K', k+1, 'Distance', 'cosine');
    row = repelem((1:n)', k); col = reshape(idx(:, 2:end)', [], 1);
    W = sparse(row, col, 1, n, n);
    SNN = W * W'; 
    A_snn = SNN ./ (2*k - SNN + eps);
    A_snn = A_snn - diag(diag(A_snn)); 
    A_snn(A_snn < 1/15) = 0; 
    A_snn = max(A_snn, A_snn'); 
end
A_SNN{1} = build_snn_graph(X_rna_pca, k_nn_static) * 5.0; 
A_SNN{2} = build_snn_graph(X_atac_lsi, k_nn_static) * 1.0;

%% 4. Spectral Initialization
disp('[4/5] Initializing latent factors via spectral fusion...');
W_fused = A_SNN{1} + A_SNN{2};
deg = sum(W_fused, 2); deg(deg == 0) = eps;
D_inv_sqrt = spdiags(deg.^(-0.5), 0, n_samples, n_samples);
L_sym = D_inv_sqrt * W_fused * D_inv_sqrt; L_sym = (L_sym + L_sym') / 2;
[V, ~] = svds(L_sym, k_clusters); V = normalize(V, 2, 'norm');
[idx_init, ~] = kmeans(V, k_clusters, 'MaxIter', 1000, 'Replicates', 15, 'Distance', 'cosine', 'Options', statset('Display','off'));

%% 5. Execute scRG-NMTF
disp('[5/5] Executing scRG-NMTF Optimization...');
fprintf('--------------------------------------------------------\n');
fprintf('| Run |   ACC    |   NMI    |   ARI    |   AMI    |\n'); 
fprintf('--------------------------------------------------------\n');

all_acc = zeros(n_runs, 1); all_nmi = zeros(n_runs, 1);
all_ari = zeros(n_runs, 1); all_ami = zeros(n_runs, 1);
best_ari = -1; best_results = struct();

for i = 1:n_runs
    % Core engine call
    [H_star, ~, ~, ~, pred_labels, obj_hist] = scrgnmf_engine(A_SNN, X_base, k_clusters, lambda1, lambda2, lambda3, k_nn_dynamic, idx_init);
    
    % Evaluation
    acc = get_acc(true_labels, pred_labels);
    nmi = get_nmi(true_labels, pred_labels);
    ari = get_ari(true_labels, pred_labels);
    ami = get_ami(true_labels, pred_labels); 
    
    all_acc(i) = acc; all_nmi(i) = nmi; all_ari(i) = ari; all_ami(i) = ami;
    fprintf('|  %d  |  %.4f  |  %.4f  |  %.4f  |  %.4f  |\n', i, acc, nmi, ari, ami);
    
    if ari > best_ari
        best_ari = ari;
        best_obj_hist = obj_hist; 
        
        best_results.H_star = H_star;
        best_results.pred_labels = pred_labels;
        best_results.acc = acc; best_results.nmi = nmi;
        best_results.ari = ari; best_results.ami = ami;
    end
end
fprintf('--------------------------------------------------------\n');

%% Summary & Output
fprintf('\n=== Final Performance (%d runs) ===\n', n_runs);
fprintf('ACC : %.4f ± %.4f\n', mean(all_acc), std(all_acc));
fprintf('NMI : %.4f ± %.4f\n', mean(all_nmi), std(all_nmi));
fprintf('ARI : %.4f ± %.4f\n', mean(all_ari), std(all_ari));
fprintf('AMI : %.4f ± %.4f\n', mean(all_ami), std(all_ami));

% Save final robust results
save('./data/scRGNMF_PBMC10k_Results.mat', 'best_results', 'true_labels', 'X_rna_pca', 'X_atac_lsi');
save('./data/Conv_PBMC10k.mat', 'best_obj_hist'); 

disp('✅ Optimization complete. Results saved to: ./data/scRGNMF_PBMC10k_Results.mat');