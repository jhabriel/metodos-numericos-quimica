# Métodos numéricos para Ciencias Químicas

Implementaciones didácticas de métodos de cálculo numérico en **MATLAB/Octave**,
pensadas para estudiantes de la Facultad de Ciencias Químicas. El código prioriza
la claridad: cada método está escrito "a mano" (sin llamar a `\`, `lu`, `interp1`,
etc.) y se valida con pruebas automáticas.

## Estructura

| Carpeta | Contenido |
|---|---|
| `src/01_ecuaciones_no_lineales` | `biseccion`, `newton`, `newton_sistema` |
| `src/02_sistemas_lineales` | `sustitucion_adelante`, `sustitucion_atras`, `gauss_naive`, `gauss_pivoteo`, `lu_factorizacion`, `resolver_lu` |
| `src/03_optimizacion` | `seccion_aurea` (1D, sin restricciones) |
| `src/04_interpolacion_regresion` | `lagrange`, `regresion_lineal_general` |
| `src/05_diferenciacion` | `dif_primera` (adelante/atrás/central), `dif_segunda_central` |
| `src/06_integracion` | `trapecio`, `simpson13`, `simpson38`, `boole` |
| `src/07_edos` | `euler_explicito`, `euler_implicito`, `heun` (escalares y sistemas) |
| `tests/` | Pruebas automáticas |
| `ejemplos/` | Ejemplos con problemas de química |

## Uso

Desde la raíz del repositorio:

```matlab
addpath(genpath('src'));

f  = @(x) x.^2 - 2;
x  = biseccion(f, 0, 2, 1e-10);
x2 = newton(f, @(x) 2*x, 1);

A = [4 -2 1; -2 4 -2; 1 -2 4];  b = [11; -16; 17];
x = resolver_lu(A, b);

I = simpson13(@(x) sin(x), 0, pi, 20);
```

Cada función documenta sus argumentos con `help nombre_funcion`.

### Ejemplo de química

```matlab
run('ejemplos/ejemplo_cinetica_arrhenius.m')
```

Ajusta la ecuación de Arrhenius para obtener Ea y resuelve la cinética A → B → C
comparando los tres métodos de EDO con la solución analítica.

## Pruebas

```matlab
run('tests/run_all_tests.m')
```

o desde la terminal con Octave:

```bash
octave --no-gui --eval "run('tests/run_all_tests.m')"
```

## Convenciones

- Nombres de funciones en español, `snake_case`, un método por archivo.
- Compatible con MATLAB y Octave (sin sintaxis exclusiva de uno u otro).
- Sistemas de EDO: `f(t, y)` devuelve columna; la salida `Y` tiene una fila por instante.
- Restricciones de `n` en integración: Simpson 1/3 → par; Simpson 3/8 → múltiplo de 3; Boole → múltiplo de 4.

## Agregar un método nuevo

1. Crear `src/<tema>/nombre.m` con encabezado de ayuda (uso, argumentos, salidas).
2. Agregar pruebas en `tests/test_<tema>.m` (comparar contra solución exacta o contra el operador nativo).
3. Listarlo en esta tabla y ejecutar `run_all_tests`.
