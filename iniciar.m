% INICIAR  Prepara la sesion de MATLAB/Octave: agrega src/ al path.
%
%   EJECUTAR UNA VEZ cada vez que se abre MATLAB u Octave (el path se olvida
%   al cerrarlos), antes de usar las funciones o los ejemplos:
%
%       1. Ir a la carpeta raiz del repositorio (donde esta este archivo).
%       2. Escribir en la ventana de comandos:   iniciar
%
%   Despues, cualquier script (en ejemplos/, examen/ o donde sea) puede llamar a
%   las funciones de src/ sin agregar nada al path.
%
%   Si no esta en la carpeta raiz, ejecutar con la ruta completa:
%       run('/ruta/al/repositorio/iniciar.m')
%
%   Si aparecen errores como "'biseccion' undefined", falta ejecutar este archivo.

raiz_repo = fileparts(mfilename('fullpath'));
carpeta_src = fullfile(raiz_repo, 'src');
if ~exist(carpeta_src, 'dir')
    error('iniciar:sinSrc', 'No se encuentra la carpeta src/ junto a iniciar.m (%s).', raiz_repo);
end
addpath(genpath(carpeta_src));
fprintf('Listo: las funciones de src/ estan disponibles en esta sesion.\n');
clear raiz_repo carpeta_src
