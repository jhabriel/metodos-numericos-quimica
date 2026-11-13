% EJ_04_AX_B_JACOBI_Y_GAUSS_SEIDEL
% Metodos iterativos para un sistema diagonalmente dominante por filas.
%
% Requisito: haber ejecutado `iniciar` al abrir MATLAB/Octave (ver README).
% Ejecutar:  run('ejemplos/ej_04_ax_b_jacobi_y_gauss_seidel.m')

clear; clc; close all;

A = [10 -1  2  0;
     -1 11 -1  3;
      2 -1 10 -1;
      0  3 -1  8];
b = [6; 25; -11; 15];            % solucion exacta: [1; 2; -1; 1]

%% Verificar dominancia diagonal por filas: |a_ii| > sum_{j~=i} |a_ij|
d = abs(diag(A));
fuera = sum(abs(A), 2) - d;
disp('Fila   |a_ii|   sum|a_ij|   dominante?');
fprintf('%3d   %6.1f   %8.1f   %d\n', [(1:4)', d, fuera, d > fuera]');

%% Resolver
tol = 1e-8;
[xj, itj, hj] = jacobi(A, b, zeros(4, 1), tol);
[xg, itg, hg] = gauss_seidel(A, b, zeros(4, 1), tol);
xref = A \ b;

fprintf('\nJacobi       : %2d iteraciones, error = %.2e\n', itj, norm(xj - xref, inf));
fprintf('Gauss-Seidel : %2d iteraciones, error = %.2e\n', itg, norm(xg - xref, inf));
fprintf('x (Gauss-Seidel) = [%s]\n', num2str(xg', ' %.8f'));

%% Convergencia
figure;
semilogy(hj(:, 1), hj(:, 2), 'o-', hg(:, 1), hg(:, 2), 's-', 'LineWidth', 1.3);
xlabel('iteracion k'); ylabel('||x_{k+1} - x_k||_\infty'); grid on;
legend('Jacobi', 'Gauss-Seidel'); title('Metodos iterativos');
