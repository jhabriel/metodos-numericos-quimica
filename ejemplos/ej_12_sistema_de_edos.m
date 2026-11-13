% EJ_12_SISTEMA_DE_EDOS
% Sistemas de EDOs con los tres metodos:
%   A) Reduccion de orden: oscilador amortiguado  x'' + 2*zeta*w*x' + w^2*x = 0
%      se escribe como sistema de primer orden con y1 = x, y2 = x'.
%   B) Modelo compartimental (farmacocinetica de dos compartimentos con
%      absorcion oral y eliminacion):
%        dA/dt  = -ka A                          (tracto gastrointestinal)
%        dC1/dt =  ka A - (k10 + k12) C1 + k21 C2  (sangre)
%        dC2/dt =  k12 C1 - k21 C2                 (tejido)
%        dE/dt  =  k10 C1                          (eliminado)
%      La masa total A + C1 + C2 + E se conserva.
%
% Ejecutar:  run('ejemplos/ej_12_sistema_de_edos.m')

clear; clc; close all;
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src')));

%% A) Oscilador amortiguado: reduccion de orden
w = 2 * pi;  zeta = 0.1;
% y = [x; v],   y' = [v; -2 zeta w v - w^2 x]
fA = @(t, y) [ y(2); -2*zeta*w*y(2) - w^2*y(1) ];
JA = @(t, y) [ 0, 1; -w^2, -2*zeta*w ];
yA0 = [1; 0];                                  % x(0) = 1, x'(0) = 0
tfA = 5; nA = 500;

wd = w * sqrt(1 - zeta^2);
xex = @(t) exp(-zeta*w*t) .* (cos(wd*t) + zeta*w/wd * sin(wd*t));

[tA, YE] = euler_explicito(fA, [0 tfA], yA0, nA);
[~,  YI] = euler_implicito(fA, [0 tfA], yA0, nA, JA);
[~,  YH] = heun(fA, [0 tfA], yA0, nA);

fprintf('A) Oscilador amortiguado, h = %.3f, error maximo |x - x_exacta|:\n', tfA / nA);
fprintf('   Euler explicito : %.3e  (la amplitud crece artificialmente)\n', max(abs(YE(:, 1) - xex(tA))));
fprintf('   Euler implicito : %.3e  (amortigua de mas)\n', max(abs(YI(:, 1) - xex(tA))));
fprintf('   Heun            : %.3e\n\n', max(abs(YH(:, 1) - xex(tA))));

figure;
plot(tA, xex(tA), 'k-', tA, YE(:, 1), 'r--', tA, YI(:, 1), 'g--', tA, YH(:, 1), 'b-.', 'LineWidth', 1.2);
xlabel('t'); ylabel('x(t)'); grid on;
legend('Exacta', 'E. explicito', 'E. implicito', 'Heun');
title('A) x'''' + 2\zeta\omega x'' + \omega^2 x = 0');

%% B) Modelo compartimental
ka = 1.2; k10 = 0.25; k12 = 0.4; k21 = 0.2;      % 1/h
K = [-ka,        0,          0,    0;
      ka, -(k10 + k12),     k21,   0;
       0,       k12,       -k21,   0;
       0,       k10,         0,    0];
fB = @(t, y) K * y;
JB = @(t, y) K;
yB0 = [100; 0; 0; 0];                            % dosis de 100 mg en A
tfB = 24; nB = 240;

[tB, ZE] = euler_explicito(fB, [0 tfB], yB0, nB);
[~,  ZI] = euler_implicito(fB, [0 tfB], yB0, nB, JB);
[~,  ZH] = heun(fB, [0 tfB], yB0, nB);

% Referencia: solucion exacta y(t) = expm(K t) y0
Zex = zeros(nB + 1, 4);
for i = 1:nB + 1
    Zex(i, :) = (expm(K * tB(i)) * yB0)';
end

fprintf('B) Modelo compartimental, h = %.2f h, error maximo en C1:\n', tfB / nB);
fprintf('   Euler explicito : %.3e\n', max(abs(ZE(:, 2) - Zex(:, 2))));
fprintf('   Euler implicito : %.3e\n', max(abs(ZI(:, 2) - Zex(:, 2))));
fprintf('   Heun            : %.3e\n', max(abs(ZH(:, 2) - Zex(:, 2))));
fprintf('   Masa total en t = %g h (Heun): %.6f mg (inicial 100)\n', tfB, sum(ZH(end, :)));

[Cmax, imax] = max(ZH(:, 2));
fprintf('   Pico en sangre (Heun): %.2f mg a las %.1f h\n', Cmax, tB(imax));

figure;
plot(tB, Zex(:, 2), 'k-', tB, ZH(:, 2), 'b--', tB, Zex(:, 3), 'k-', tB, ZH(:, 3), 'r--', 'LineWidth', 1.2);
xlabel('t (h)'); ylabel('cantidad (mg)'); grid on;
legend('Exacta', 'Heun C_1 (sangre)', '', 'Heun C_2 (tejido)');
title('B) Modelo de dos compartimentos');
