% =========================================================================
% Script to plot the convergence curve of scRG-NMTF
% =========================================================================

clear; clc; close all;

% Load the saved convergence history
try
    load('./data/Conv_PBMC10k.mat', 'best_obj_hist');
catch
    error('Convergence file not found. Please run demo_PBMC10k.m first.');
end

% Filter out zero entries if the algorithm converged early
valid_iters = best_obj_hist > 0;
iters = 1:sum(valid_iters);
obj_vals = best_obj_hist(valid_iters);

% Plotting
figure('Position', [100, 100, 600, 450], 'Color', 'w');
plot(iters, obj_vals, '-', 'LineWidth', 2.5, 'Color', [0, 0.4470, 0.7410]);
grid on;
set(gca, 'GridAlpha', 0.15, 'FontSize', 12, 'LineWidth', 1.2);
set(gca, 'YScale', 'log'); % Log scale for residuals

% Labels and Title
xlabel('Iterations', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Objective Function Residual (Frobenius Norm)', 'FontSize', 14, 'FontWeight', 'bold');
title('Convergence Curve of scRG-NMTF', 'FontSize', 16, 'FontWeight', 'bold');