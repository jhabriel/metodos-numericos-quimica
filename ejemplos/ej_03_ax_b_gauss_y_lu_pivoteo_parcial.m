% EJ_03_AX_B_GAUSS_Y_LU_PIVOTEO_PARCIAL
% Sistemas donde Gauss sin pivoteo falla (pivote nulo) o es inexacto (pivote
% diminuto). El pivoteo parcial por filas lo resuelve en ambos casos.
%
% Ejecutar:  run('ejemplos/ej_03_ax_b_gauss_y_lu_pivoteo_parcial.m')

clear; clc; close all;
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src')));

%% Parte 1: pivote nulo (a11 = 0)
A = [0 2 1;
     1 1 1;
     2 1 3];
b = [3; 3; 7];

disp('--- Parte 1: A(1,1) = 0 ---');
try
    gauss_naive(A, b);
catch err
    fprintf('gauss_naive falla: %s\n', err.message);
end

x_piv = gauss_pivoteo(A, b);
[L, U, P] = lu_factorizacion(A);
x_lu  = sustitucion_atras(U, sustitucion_adelante(L, P * b));

disp('P ='), disp(P)
disp('L ='), disp(L)
disp('U ='), disp(U)
fprintf('||P*A - L*U|| = %.2e\n', norm(P * A - L * U));
fprintf('x (Gauss con pivoteo) = [%s]\n', num2str(x_piv', ' %.8f'));
fprintf('x (LU con pivoteo)    = [%s]\n', num2str(x_lu', ' %.8f'));
fprintf('x (A\\b, referencia)   = [%s]\n\n', num2str((A \ b)', ' %.8f'));

%% Parte 2: pivote diminuto (inestabilidad numerica)
epsi = 1e-20;
B = [epsi 1;
     1    1];
c = [1; 2];                      % solucion exacta ~ [1; 1]

disp('--- Parte 2: pivote 1e-20 ---');
x1 = gauss_naive(B, c);
x2 = gauss_pivoteo(B, c);
fprintf('Gauss naive    : x = [%.6f  %.6f]   (x1 es incorrecto)\n', x1);
fprintf('Gauss pivoteo  : x = [%.6f  %.6f]\n', x2);
fprintf('Exacta (aprox) : x = [1.000000  1.000000]\n');
