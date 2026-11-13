% EJEMPLO_CINETICA_ARRHENIUS
% 1) Ajuste de Arrhenius (regresion lineal generalizada) -> energia de activacion
% 2) Cinetica A -> B -> C resuelta con Euler explicito, Heun y Euler implicito
%
% Ejecutar desde la raiz del repositorio:  run('ejemplos/ejemplo_cinetica_arrhenius.m')

raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(raiz, 'src')));

%% 1) Arrhenius: ln k = ln A - (Ea/R) * (1/T)
R = 8.314;                                   % J/(mol K)
T = [298 308 318 328 338]';                  % K
k = [1.2e-4 3.9e-4 1.15e-3 3.1e-3 7.9e-3]';  % 1/s (datos de ejemplo)

[a, r2] = regresion_lineal_general(1 ./ T, log(k), {@(x) ones(size(x)), @(x) x});
Ea = -a(2) * R;
A  = exp(a(1));
fprintf('Ea = %.1f kJ/mol,  A = %.3e 1/s,  R^2 = %.5f\n', Ea / 1e3, A, r2);

%% 2) A -> B -> C
k1 = 1.0; k2 = 0.5;
f = @(t, c) [-k1 * c(1); k1 * c(1) - k2 * c(2); k2 * c(2)];
c0 = [1; 0; 0];

[t, Ce] = euler_explicito(f, [0 10], c0, 100);
[~, Ch] = heun(f, [0 10], c0, 100);
[~, Ci] = euler_implicito(f, [0 10], c0, 100);
cB = @(t) k1 / (k2 - k1) * (exp(-k1 * t) - exp(-k2 * t));   % solucion analitica

fprintf('\nError en [B](t=10), 100 pasos:\n');
fprintf('  Euler explicito : %.2e\n', abs(Ce(end, 2) - cB(10)));
fprintf('  Heun            : %.2e\n', abs(Ch(end, 2) - cB(10)));
fprintf('  Euler implicito : %.2e\n', abs(Ci(end, 2) - cB(10)));

figure; hold on;
plot(t, cB(t), 'k-', 'LineWidth', 1.5);
plot(t, Ce(:, 2), 'r--', t, Ch(:, 2), 'b-.', t, Ci(:, 2), 'g:', 'LineWidth', 1.2);
xlabel('t (s)'); ylabel('[B] (mol/L)');
legend('Analitica', 'Euler explicito', 'Heun', 'Euler implicito');
title('A \rightarrow B \rightarrow C');
grid on;
