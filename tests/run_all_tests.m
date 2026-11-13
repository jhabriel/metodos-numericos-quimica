% RUN_ALL_TESTS  Ejecuta todas las pruebas de tests/ (MATLAB y Octave).
%
%   Uso (desde la raiz del repositorio):  run('tests/run_all_tests.m')

raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(raiz, 'src')));
addpath(fullfile(raiz, 'tests'));

pruebas = {'test_no_lineales', 'test_sistemas_lineales', 'test_optimizacion', ...
           'test_interpolacion_regresion', 'test_edos'};

fallos = 0;
for i = 1:numel(pruebas)
    try
        feval(pruebas{i});
        fprintf('[ OK ] %s\n', pruebas{i});
    catch err
        fallos = fallos + 1;
        fprintf('[FAIL] %s: %s\n', pruebas{i}, err.message);
    end
end

if fallos == 0
    fprintf('\nTodas las pruebas pasaron (%d archivos).\n', numel(pruebas));
else
    error('%d archivo(s) de prueba fallaron.', fallos);
end
