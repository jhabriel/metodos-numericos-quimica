% EJ_01_ECUACION_NO_LINEAL
% Raiz de f(x) = x^3 - 2x - 5 = 0 (ecuacion clasica de Wallis/Newton)
% por biseccion y por Newton-Raphson.
%
% Ejecutar:  run('ejemplos/ej_01_ecuacion_no_lineal.m')

clear; clc; close all;
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src')));

f  = @(x) x.^3 - 2*x - 5;
df = @(x) 3*x.^2 - 2;
tol = 1e-10;

% f(2) = -1 < 0 y f(3) = 16 > 0  ->  hay una raiz en [2, 3]
[xb, itb, hb] = biseccion(f, 2, 3, tol);
[xn, itn, hn] = newton(f, df, 2, tol);
xref = fzero(f, 2);              % referencia (funcion nativa)

fprintf('Raiz de referencia (fzero) : %.12f\n\n', xref);
fprintf('Biseccion : x = %.12f  (%2d iteraciones, error %.2e)\n', xb, itb, abs(xb - xref));
fprintf('Newton    : x = %.12f  (%2d iteraciones, error %.2e)\n\n', xn, itn, abs(xn - xref));

fprintf('Newton, historia:\n   k        x_k            f(x_k)\n');
fprintf('%4d  %.12f  %12.3e\n', hn');

% --- Graficas ---
figure;
subplot(1, 2, 1);
xx = linspace(1.5, 3, 200);
plot(xx, f(xx), 'b-', xx, 0 * xx, 'k-', xn, 0, 'ro', 'LineWidth', 1.3);
xlabel('x'); ylabel('f(x)'); title('f(x) = x^3 - 2x - 5'); grid on;

subplot(1, 2, 2);
semilogy(hb(:, 1), max(abs(hb(:, 4) - xref), eps), 'o-', ...
         hn(:, 1), max(abs(hn(:, 2) - xref), eps), 's-', 'LineWidth', 1.3);
xlabel('iteracion'); ylabel('|x_k - x^*|'); grid on;
legend('Biseccion', 'Newton'); title('Convergencia');
