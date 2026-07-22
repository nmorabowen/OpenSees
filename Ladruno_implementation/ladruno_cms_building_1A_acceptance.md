# Building 1A — protocolo de aceptación de LadrunoCMS

## Propósito

Este documento fija y registra el caso real de P2 del
[ADR 1000](1000_ladruno_cms_adr.md). La biblioteca ya fue implementada y compilada en
la rama `feature/ladruno-cms-adr1000`; el protocolo se ejecutó primero sin corrector
global, produjo un `FAIL` reproducible y, después de la remediación fundamentada en los
papers, cerró con veredicto final `PASS`. Se conservan ambos estados para no reescribir la
historia experimental.

El notebook [building_1A_manual_run.ipynb](building_1A_manual_run.ipynb) se conserva como
baseline y generador del modelo. El deck CMS se deriva mediante
`prepare_building_1A_cms.sh`; el notebook y el deck de referencia no se modifican.

## Baseline congelado

| Artefacto | SHA-256 |
|---|---|
| `building_1A_manual_run.ipynb` | `2c3ca6099574eb3554b2935aa5dc6ca143f8f57d48885cb2a75c5a30be9a8603` |
| `modal_parallel.tcl` | `74acc6ccded01e9b8dcf69a24c4b0028e9894c9c72947042daa4d7ef3aee2d22` |
| `eigenvalues.out` | `1d144d0f8c1b787974bbfab19c608ea3f622f5b19672f64726d6bf22fb6d94df` |
| `modal_results.h5` | `2837ee284cff25f55c216f076d5e0c69be3b27b0572256bb9f92a1bcb7db2287` |

Los tres últimos archivos están actualmente en:

```text
/home/pxpalacios/Dropbox/01. Brain/10. Ph.D U ANDES/10. obsidian PXP/SRF/1A/1A_manual_run/modal_parallel
```

El deck tiene 36 634 líneas y aproximadamente 2.1 MiB; el HDF5 de referencia ocupa
aproximadamente 18 MiB. Si un hash cambia, se crea una nueva versión del baseline y se
documenta la razón antes de comparar.

## Inventario del modelo

| Dato | Valor congelado |
|---|---:|
| nodos | 11 841 |
| elementos | 27 360 |
| particiones de autoría Gmsh | 4 |
| nodos de base fijos en 6 DOF | 1 333 |
| masa sísmica | 11 972.931 t |
| DOF por nodo | 6 |
| ecuaciones activas estimadas | 63 048 |
| coordenadas dinámicas estimadas | 31 524 |
| coordenadas sin masa estimadas | 31 524 |
| modos de aceptación | 8 |

La cuenta de ecuaciones es `6*(11841-1333)`; las cuentas `D/Z` son cada una
`3*(11841-1333)`. Deben compararse con el `Graph` y el pencil que recibe el SOE; si
difieren por constraints adicionales, los valores efectivos se convierten en los
reportados y se explica. Las líneas `mass` asignan masa a traslaciones y cero a rotaciones.
Éste es el caso que obliga a probar la condensación estática de coordenadas sin masa.

Las cuatro particiones de Gmsh no se entregan como partición CMS. El input modal es un
modelo plano replicado y la biblioteca aplica METIS al grafo de ecuaciones dentro del
SOE. Reutilizar los labels de Gmsh invalidaría esta aceptación porque probaría otro flujo
de información.

## Preparación reproducible

El script copia el deck congelado a un directorio de resultados CMS, conservando sin cambios
todos los comandos hasta `system UmfPack`. Sólo se reemplazará la sentencia de eigensolve
y se añadirá salida diagnóstica CMS; nodos, elementos, materiales, secciones, masses,
constraints, numberer y recorders permanecerán byte-equivalentes. El diff del deck será
un artefacto de la prueba.

La invocación inicial es:

```tcl
set _lam [eigen -ladrunoCMS \
    -hierarchy logical -level1 2 -level2 2 \
    -modesL2 24 -modesL1 48 \
    -tol 1.0e-8 -maxEnrich 4 -maxIter 500 \
    -refinement subspace -iterationVectors auto -maxRefineIter 160 \
    -partition metis -localEigen lanczos \
    -globalSolver dense -denseMax 2000 \
    -verifyAssembly signature -verbose \
    8]
```

Se ejecuta con el binario construido desde `feature/ladruno-cms-adr1000`:

```bash
mpiexec -np 4 /home/pxpalacios/builds/build-ladruno-cms-adr1000/OpenSeesMP modal_cms.tcl
```

El ejecutable verificado está en la ruta indicada. El script acepta overrides mediante
`CMS_MODES_LEVEL2`, `CMS_MODES_LEVEL1`, `CMS_TOLERANCE`, `CMS_MAX_ENRICH`,
`CMS_MAX_ITERATIONS`, `CMS_REFINEMENT`, `CMS_ITERATION_VECTORS`,
`CMS_MAX_REFINE_ITERATIONS` y `CMS_DENSE_MAX`, sin alterar el baseline congelado.

## Preflight obligatorio

Antes del solve denso, el log debe contener como mínimo:

- commit, opciones y tamaño MPI;
- `p=2`, `m=2` y composición de los tres comunicadores;
- número efectivo de ecuaciones y contribuciones `A/M`;
- conteos de coordenadas dinámicas/sin masa en cada eigensolve;
- tamaños interiores, de borde y retenidos para cada subdominio fino y grueso;
- dimensiones antes/después de `T_2`, `S_2`, `T_1` y `S_1`;
- `r_dup`, `r_raw`, `r_Z`, `r_D`, `sum(k1_g)`, `|union B_g|`, `denseMax` y memoria raíz
  estimada;
- memoria persistente final `8*n*numModes` por rank;
- residuales locales, reinicios, enriquecimientos y residual global por modo.

Si `r_D>2000`, el solver debe fallar colectivamente antes de asignar las matrices densas. Se
puede repetir con un `denseMax` mayor únicamente si la cota de memoria y el margen del nodo
han sido registrados. Si no cabe, no se reduce `modesL1/modesL2` hasta forzar un resultado
de residual insuficiente: el caso se convierte en puerta del backend distribuido de P3.

## Resultados de referencia

Los autovalores congelados, en orden ascendente, son:

```text
14.448045635233683
14.661827604161822
32.92074680793754
346.24936007060205
358.60147903440566
884.0394435484916
1102.3146722691856
1193.743437958608
```

Los periodos derivados son aproximadamente:

```text
1.65300937, 1.64091398, 1.09507794, 0.33766448,
0.33179805, 0.21132168, 0.18924616, 0.18185463
```

## Criterios PASS/FAIL

Todos son obligatorios:

1. cuatro ranks terminan con el mismo estado y sin deadlock;
2. el log demuestra ejecución no degenerada de `T_2`, `S_2`, `T_1`, `S_1` y del camino
   inverso;
3. se obtienen ocho autovalores finitos, positivos y ordenados;
4. error relativo por autovalor `<=1e-8` respecto del baseline;
5. residual del pencil original `<=1e-8` por modo, incluidas filas sin masa;
6. MAC `>=1-1e-8` para modos simples y MAC de subespacio para cualquier cluster;
7. salto de copias de interfaz dentro de la tolerancia de retro-sustitución;
8. `nodeEigenvector` tiene seis componentes por nodo y `modalProperties` termina;
9. dos corridas idénticas conservan valores y subespacios dentro de tolerancia;
10. una corrida con retención enriquecida no aumenta los Ritz fuera del margen de
    redondeo;
11. se archivan comando, hashes, log, autovalores, mode shapes, diagnósticos, tiempos y
    memoria máxima por rank.

El tiempo no decide P2. La comparación de rendimiento y la ablación controlada de niveles
pertenecen a P4; ningún resultado de un solo nivel puede sustituir esta aceptación.

## Ejecución inicial sin refinamiento global — evidencia histórica

La corrida activa usó cuatro ranks, `p=2`, `m=2`, `k2=24`, `k1=48`, ocho modos y
`denseMax=2000`. El preflight informó `n=63048`, `dim(Z)=31524`, `r_raw=r_D=486`; por
consiguiente, el límite denso no fue el causante del fallo. La cadena completa de dos
niveles y el camino inverso se ejecutaron sin deadlock.

Después de la reconstrucción estática distribuida de las coordenadas originales sin masa,
se obtuvieron:

```text
14.4550, 14.6689, 33.0778, 356.255,
367.264, 1019.29, 1200.53, 1377.41
```

Los errores relativos por modo fueron:

```text
0.0481 %, 0.0482 %, 0.4771 %, 2.8897 %,
2.4156 %, 15.2992 %, 8.9099 %, 15.3858 %
```

El máximo residual normalizado del pencil original fue `0.916338`; el tiempo de pared fue
`151.599 s`. En consecuencia, fallan los criterios 4 y 5 y no corresponde evaluar los
criterios de MAC como si la corrida hubiera sido aceptada. El log reproducible es
`building_1A_cms_run/modal_cms_global_z_cleanup.log`.

El segundo punto limpio duplicó la retención a `k2=48`, `k1=96`, sin enriquecimiento ni
correctores. Produjo `r_D=582`, `rho_max=0.902530` y `247.329 s`. Los autovalores
fueron `14.4549, 14.6679, 32.9464, 355.857, 366.949, 984.85, 1116.42, 1215.57`; sus
errores relativos aproximados fueron `0.0474 %, 0.0414 %, 0.0779 %, 2.7748 %, 2.3278 %,
11.4034 %, 1.2796 %, 1.8284 %`. La retención mayor mejora algunos modos, pero el residual
y el modo 6 permanecen lejos de la puerta. El log es
`building_1A_cms_run/modal_cms_active_48_96.log`.

Una corrida exploratoria con `k2=48`, `k1=96` y correcciones diagonal/Schwarz alcanzó
`r_D=582`, `rho_max=0.904959` y `249.696 s`. Esas correcciones se rechazaron y se retiraron
del camino activo: siete autovalores se aproximaron mejor, pero el modo 6 mantuvo un error
de `5.6353 %` y el residual siguió siendo inaceptable. El log se conserva únicamente como
evidencia negativa en `building_1A_cms_run/modal_cms_retention_48_96_corrected.log`.

| Criterio | Resultado |
|---|---|
| cuatro ranks y estado colectivo | PASS |
| `T_2`, `S_2`, `T_1`, `S_1` y reconstrucción | PASS |
| ocho Ritz positivos y ordenados antes de la puerta | PASS |
| error relativo por autovalor `<=1e-8` | **FAIL** |
| residual original por modo `<=1e-8` | **FAIL** |
| publicación ordinaria en OpenSees | no alcanzada: el solver falla de forma explícita |
| veredicto P2c | **FAIL** |

No debe aumentarse la tolerancia sólo para obtener un retorno exitoso. La próxima corrida
de aceptación requiere primero un estudio de convergencia de `k2/k1` y una estrategia de
refinamiento con respaldo matemático y pruebas unitarias. El criterio aproximado reportado
por el paper y la certificación estricta de esta ADR se documentarán por separado.

## Remediación ejecutada

La remediación conservó los dos niveles CMS y añadió una iteración de subespacios sobre
el pencil original. La jerarquía produce `q=max(p+8,2p)=16` vectores iniciales; cada
iteración resuelve `K Xbar=M X`, ortonormaliza respecto de `M`, ejecuta Rayleigh–Ritz y
evalúa el residual original. La rigidez se factoriza una vez. Cada acción inversa usa dos
correcciones de bloque `K Delta=B-KX`, necesarias porque el refinamiento iterativo interno
de MUMPS no opera con múltiples lados derechos.

El límite inicial de 20 iteraciones redujo el residual de `0.916338` a `0.00392268`, pero
no convergió. Con 160 iteraciones y sin correcciones de solve llegó al piso numérico
`1.63256e-8`. Las correcciones de bloque eliminaron ese piso sin cambiar la tolerancia.
Estos dos fallos se conservan en:

- `building_1A_cms_run/modal_cms_subspace_refinement_20iter.log`;
- `building_1A_cms_run/modal_cms_subspace_refinement_160iter.log`.

## Resultado final de aceptación

La configuración base fue `p=2`, `m=2`, `k2=24`, `k1=48`, `q=16`, `r_D=486`, ocho
modos y `tol=1e-8`. Dos corridas independientes terminaron sin deadlock y publicaron las
formas modales ordinarias de OpenSees:

| Métrica | Corrida 1 | Corrida 2 |
|---|---:|---:|
| iteraciones globales | 123 | 122 |
| residual original máximo | `9.09852e-9` | `9.81460e-9` |
| factorización original | `0.200568 s` | `0.197138 s` |
| solves de bloque totales | `58.0684 s` | `59.0481 s` |
| fracción atribuida a correcciones | `46.7490 s` | `47.2625 s` |
| Rayleigh–Ritz | `18.0845 s` | `17.6850 s` |

Los autovalores de la primera corrida fueron:

```text
14.448045635125181
14.661827604138683
32.920746807876327
346.24936007056868
358.60147903441879
884.03944354845726
1102.3146722691142
1193.7434379586321
```

El error relativo máximo frente al baseline congelado fue `7.510e-12`. La comparación de
los vectores publicados dio MAC diagonal mínimo `0.9999999999997498` frente al baseline.
Entre las dos corridas, la diferencia relativa máxima de autovalores fue `2.339e-12` y el
MAC diagonal mínimo `0.9999999999999278`.

La ablación de retención usó `k2=48`, `k1=96`, `q=16` y `r_D=582`. Convergió en 73
iteraciones con residual `9.68682e-9`; respecto de la corrida base, la diferencia relativa
máxima de autovalores fue `1.609e-12` y el MAC mínimo `0.9999999999997349`. La retención
mayor mejora el espacio inicial y reduce el número de iteraciones, pero no cambia el
resultado físico aceptado.

| Criterio | Resultado final |
|---|---|
| cuatro ranks y estado colectivo | PASS |
| `T_2`, `S_2`, `T_1`, `S_1` y reconstrucción | PASS |
| ocho autovalores finitos, positivos y ordenados | PASS |
| error relativo por autovalor `<=1e-8` | PASS |
| residual original por modo `<=1e-8` | PASS |
| MAC de modos publicados | PASS |
| dos corridas idénticas | PASS |
| retención enriquecida | PASS |
| publicación `nodeEigenvector`/recorders | PASS |
| veredicto P2c | **PASS** |

Los logs normativos son `modal_cms_subspace_refinement_run1.log`,
`modal_cms_subspace_refinement_run2.log` y `modal_cms_enriched_refinement.log`. Las salidas
de las dos corridas base están separadas en `accepted_run1/` y `accepted_run2/`. Durante
el ensamblaje replicado se observó por muestreo un máximo aproximado de 10.6 GB RSS
agregado; no es una medición instrumentada de pico y se conserva como riesgo para P3/P4.

## Extensión P3 — primer Building 1A físicamente distribuido

La rama `feature/ladruno-cms-physical-domain` ejecutó el mismo edificio mediante cuatro
fragmentos de input. Cada intérprete OpenSeesMP construyó únicamente su `Domain` local:

| Rank | nodos | elementos OpenSees |
|---:|---:|---:|
| 0 | 3,190 | 2,959 |
| 1 | 3,260 | 3,004 |
| 2 | 3,009 | 2,864 |
| 3 | 2,823 | 2,769 |

El universo emitido contiene 11,841 nodos y 11,596 elementos OpenSees. Ningún rank cargó
ese universo completo. Los 12,282 memberships nodales incluyen 441 copias de interfaz.
Las 27,360 entidades crudas de `FEMData` son un universo de generación distinto y se
registran separadamente en el manifest.

La corrida usó `-domainMode physical`, `p=2`, `m=2`, `k2=24`, `k1=48`, ocho modos,
`tol=1e-8` y `maxIter=500`. El problema tuvo `n=63048`, `r2=2940` y `r_D=564`.

| Criterio físico P3 | Resultado |
|---|---:|
| residual original máximo | `9.84323e-9` |
| error relativo máximo de autovalores | `4.474074e-12` |
| MAC mínimo | `0.9999999999997988` |
| salto máximo de interfaz | `0` |
| tiempo de pared | `319.410324 s` |
| primer particionado físico de cuatro ranks | **PASS** |

El caso reveló una meseta del solve local cuando se exigía `0.1*tol`: aumentar Lanczos de
500 a 1200 iteraciones no redujo un residual local de aproximadamente `3.7e-9`. La ruta
física usa ahora `min(1e-8,tol)` localmente y mantiene el refinamiento del pencil original
como criterio global estricto. El camino `replicatedReference` no cambió.

Los picos RSS observados fueron aproximadamente 4095, 4142, 2689 y 1272 MiB por rank. El
modelo no está replicado, pero la reducción Craig--Bampton actual crea matrices locales
densas y los vectores finales siguen siendo globales. Por ello este PASS demuestra
integración y corrección, no speedup ni escalabilidad de memoria.

La evidencia completa está fuera del repositorio fuente, en
`notebooks/building_1A_cms_physical_acceptance.md` y
`notebooks/building_1A_cms_physical_run/`. Permanecen abiertos la repetición, un segundo
particionado, el oráculo explícito `Kx/Mx`, fixtures negativos de manifest y la
instrumentación por fase.

## Baseline físico 2/4/6 y comparación de solver

La instrumentación posterior repitió el caso físico con ocho modos y snapshots frescos
para cada número de particiones. CMS-4 repitió el PASS con solve `311.606312 s`, residual
`9.84319e-9` y error relativo máximo `2.5883e-12`. CMS-6 pasó con solve `208.730637 s`,
residual `9.65914e-9` y error `2.80898e-12`. CMS-2 verificó dominios físicos pero abortó
en el ordenamiento MUMPS (`orderMinPriority`).

El desglose interno fue:

| Ranks | ensamblaje [s] | jerarquía [s] | refinamiento [s] | total [s] |
|---:|---:|---:|---:|---:|
| 4 | 0.21356 | 214.667 | 96.5068 | 311.389 |
| 6 | 0.146962 | 125.641 | 82.7916 | 208.581 |

En el mismo build, ARPACK secuencial resolvió en `30.126791 s`; FEAST certificado empleó
`376.666740`, `391.979597` y `448.826460 s` a 2, 4 y 6 ranks. CMS supera a FEAST en los
dos puntos no degenerados, pero no a ARPACK. Los picos agregados CMS de 10.29 GiB (4) y
6.72 GiB (6) impiden declarar escalabilidad de memoria, pese a la propiedad física
correcta.
