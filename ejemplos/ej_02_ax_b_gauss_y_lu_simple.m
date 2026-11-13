% EJ_02_AX_B_GAUSS_Y_LU_SIMPLE
% Resuelve A x = b (A y b fijos) con eliminacion de Gauss sin pivoteo y con
% factorizacion LU, mostrando los pasos intermedios.
%
% Ejecutar:  run('ejemplos/ej_02_ax_b_gauss_y_lu_simple.m')

clear; clc; close all;
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src')));

A = [ 4 -2  1;
     -2  4 -2;
      1 -2  4];
b = [11; -16; 17];

%% 1) Gauss naive
x_gauss = gauss_naive(A, b);

%% 2) Factorizacion LU:  P*A = L*U,  luego  L y = P b  (adelante),  U x = y (atras)
[L, U, P] = lu_factorizacion(A);
y = sustitucion_adelante(L, P * b);
x_lu = sustitucion_atras(U, y);

%% Resultados
disp('L ='), disp(L)
disp('U ='), disp(U)
disp('P ='), disp(P)
fprintf('||P*A - L*U|| = %.2e\n\n', norm(P * A - L * U));

fprintf('       Gauss naive        LU           A\\b (referencia)\n');
fprintf('x%d = %12.8f  %12.8f  %12.8f\n', [1:3; x_gauss'; x_lu'; (A \ b)']);

fprintf('\nResiduo ||A x - b||:  Gauss = %.2e,  LU = %.2e\n', ...
        norm(A * x_gauss - b), norm(A * x_lu - b));
fprintf('det(A) = prod(diag(U)) * (+-1) = %.4f\n', prod(diag(U)) * det(P));
