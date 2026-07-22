# CMS jerárquico — referencia matemática y validación independiente

## 1. Propósito

Esta carpeta contiene un oráculo numérico independiente para la biblioteca nueva
`OPS_LadrunoCMS` definida en [ADR 1000](../1000_ladruno_cms_adr.md). No importa OpenSees ni
reutiliza un eigensolver del repositorio. Su función es fijar las identidades matemáticas,
los flujos de compatibilidad y los criterios de aceptación antes de escribir C++ o MPI.

Existen dos niveles de verificación:

1. un oráculo algebraico sobre matrices globales `K`, `M` y labels de partición;
2. un oráculo basado en contribuciones de elemento que ejecuta de forma explícita las
   cuatro transformaciones de Yu et al.

El primero facilita casos densos, particiones no contiguas y patrones creados sólo por
masa. El segundo es la referencia normativa para la implementación jerárquica.

## 2. Contrato con OpenSees

La biblioteca C++ recibirá el pencil ya restringido y numerado por la ruta ordinaria:

```text
Domain / constraints / DOF numbering
              |
              v
setSize(Graph) -> zeroA/zeroM -> addA/addM
              |
              v
        LadrunoCMS::solve
              |
              v
getEigenvalue / getEigenvector
              |
              v
AnalysisModel -> Domain -> nodeEigenvector/modalProperties
```

Por tanto, el oráculo trabaja en espacio de ecuaciones. Las restricciones no se vuelven a
aplicar dentro de CMS. Los vectores reconstruidos tienen longitud `n`, donde `n` es el
número global de ecuaciones activas.

## 3. Craig–Bampton local

Para un subdominio `s`, las coordenadas se dividen en interiores `I` e interfaces `B`:

\[
K^{(s)}=
\begin{bmatrix}K_{II}&K_{IB}\\K_{BI}&K_{BB}\end{bmatrix},\qquad
M^{(s)}=
\begin{bmatrix}M_{II}&M_{IB}\\M_{BI}&M_{BB}\end{bmatrix}.
\]

Los modos de interfaz fija son los primeros `k_s` pares de

\[
K_{II}\Phi_s=M_{II}\Phi_s\Lambda_s.
\]

Los modos de restricción se obtienen de una sola factorización de `K_II` con múltiples
RHS:

\[
\Psi_s=-K_{II}^{-1}K_{IB}.
\]

La transformación local es

\[
T_s=
\begin{bmatrix}
\Phi_s & \Psi_s\\
0 & I_B
\end{bmatrix}.
\]

Las matrices reducidas son congruencias completas:

\[
\widetilde K_s=T_s^TK^{(s)}T_s,
\qquad
\widetilde M_s=T_s^TM^{(s)}T_s.
\]

Esto conserva todos los términos de masa consistente. En particular, el bloque reducido
de borde contiene

\[
M_{BB}+M_{BI}\Psi+\Psi^TM_{IB}+\Psi^TM_{II}\Psi.
\]

Omitir `M_IB` o sumar `M_BB` una vez por vecino produce un pencil distinto y queda cubierto
por los tests.

## 4. Las cuatro transformaciones del paper

### 4.1 Nivel 2: reducción local

Cada una de las `p*m` particiones finas posee contribuciones elementales disjuntas. Puede
duplicar coordenadas físicas de interfaz, pero nunca una contribución de elemento. Se
aplica `T_2^(s)` de manera independiente.

### 4.2 Nivel 2: compatibilidad

Dentro de cada partición de nivel 1, `S_2` mapea coordenadas únicas a sus copias locales:

\[
K_2=S_2^T\operatorname{blkdiag}(\widetilde K_2^{(s)})S_2,
\qquad
M_2=S_2^T\operatorname{blkdiag}(\widetilde M_2^{(s)})S_2.
\]

`S_2` es booleana y cada fila tiene exactamente un uno. Dos filas asociadas al mismo DOF
físico señalan la misma columna.

### 4.3 Nivel 1: segunda reducción

El pencil ensamblado de cada grupo se divide otra vez en coordenadas interiores y de
interfaz gruesa. Las coordenadas modales finas son interiores de nivel 1. Se aplica
`T_1^(g)` con su propia retención `k_1`.

### 4.4 Nivel 1: compatibilidad global

Los líderes unen las copias de interfaz gruesa con `S_1`:

\[
K_*=S_1^T\operatorname{blkdiag}(\widetilde K_1^{(g)})S_1,
\qquad
M_*=S_1^T\operatorname{blkdiag}(\widetilde M_1^{(g)})S_1.
\]

Una implementación que sólo calcula `T_fine @ T_coarse` sobre un pencil global verifica
una reducción algebraica de dos niveles, pero no verifica duplicación/compatibilidad. Por
eso `run_paper_hierarchical_cms` conserva explícitamente `T_2`, `S_2`, `T_1` y `S_1`.

## 5. Solución y retro-sustitución

El problema final es

\[
K_*q_j=\lambda_jM_*q_j.
\]

La reconstrucción invierte el flujo:

\[
q_j\to S_1q_j\to T_1S_1q_j\to S_2T_1S_1q_j
\to T_2S_2T_1S_1q_j=\phi_j.
\]

El oráculo compara las filas reconstruidas de todas las copias de un DOF antes de
promediarlas. La diferencia debe ser cero dentro de round-off; el promedio sólo traduce
las copias consistentes a un vector global.

El residual normalizado es

\[
\rho_j=
\frac{\|K\phi_j-\lambda_jM\phi_j\|_2}
{\|K\phi_j\|_2+|\lambda_j|\|M\phi_j\|_2}.
\]

Para autovalores simples se usa MAC ponderado por masa. Para autovalores repetidos o
agrupados se comparan subespacios mediante los cosenos principales, porque las bases
individuales pueden rotar sin cambiar el eigenspace.

## 6. Coordenadas sin masa

El edificio 1A tiene masa nodal traslacional y masa rotacional nula. Para no convertir
esta propiedad habitual de modelos frame/shell en un rechazo artificial, el oráculo
implementa la misma extensión acotada del ADR. Con coordenadas dinámicas `D` y sin masa
`Z`, sólo acepta

\[
M=\begin{bmatrix}M_{DD}&0\\0&0\end{bmatrix},\quad M_{DD}\succ0,
\quad K_{ZZ}\succ0.
\]

La transformación estática

\[
G=\begin{bmatrix}I\\-K_{ZZ}^{-1}K_{ZD}\end{bmatrix}
\]

produce `Khat=G^T K G` y `Mhat=M_DD`. Después del solve, `phi=G phihat` recupera también
las rotaciones sin masa y el residual se calcula con `K` y `M` originales. Las pruebas
cubren una y varias coordenadas `Z`, comparan los autovalores finitos del pencil singular
y rechazan tanto `K_ZZ` singular como una fila nula sólo en diagonal pero acoplada en
`M`. No se afirma soporte para un nullspace general.

## 7. Modelo determinista

Las pruebas principales usan una cadena fija-fija de `n` grados de libertad y `n+1`
elementos. Cada resorte tiene

\[
k_e=k\begin{bmatrix}1&-1\\-1&1\end{bmatrix}.
\]

Se prueban dos masas:

\[
m_e^{lumped}=\frac m2 I,
\qquad
m_e^{consistent}=\frac m6\begin{bmatrix}2&1\\1&2\end{bmatrix}.
\]

La masa consistente activa `M_IB != 0`. Para masa concentrada, la solución analítica es

\[
\lambda_j=\frac{4k}{m}\sin^2\left(\frac{j\pi}{2(n+1)}\right),
\qquad j=1,\ldots,n.
\]

Los elementos se reparten contiguamente y sin solape entre las particiones finas. Las
interfaces se descubren a partir de nodos usados por más de un dueño.

## 8. Batería actual

`test/test_cms_level2.py` contiene 30 tests:

| Grupo | Evidencia |
|---|---|
| referencia | cadena analítica fija-fija |
| extremos | partición sin interiores |
| ensamblaje | `K_BB/M_BB` una vez y `M_IB` consistente |
| exactitud | bases completas, masa lumped y consistent |
| aproximación | cota superior Ritz y enriquecimiento monótono |
| dos niveles algebraicos | segunda condensación real y caso completo exacto |
| patrones | pencil denso, interiores no contiguos y acoplamiento sólo por masa |
| propiedad efectiva | elemento que cruza grupos gruesos, dueño único y promoción a `S_1` |
| clusters | MAC de subespacio |
| masa semidefinida | condensación/reconstrucción `G` y residual original |
| guards | nullspace de masa no alineado, `K_ZZ` y `K_II` singulares |
| colectivas | unión estable, claves duplicadas, mapas y líder vacío |
| memoria | fórmula exacta de matrices densas y gather triangular |
| firmas | estructura exacta y métricas numéricas tolerantes |
| flujo Yu | cuatro transformaciones exactas con ambas masas |
| compatibilidad | una entrada por fila, duplicados realmente fusionados |
| truncación Yu | ambas reducciones activas y cotas Ritz |
| jerarquía | partición imposible y petición modal inválida rechazadas |
| refinamiento global | residual original decreciente y recuperación del pencil conocido |
| control negativo | iterar sólo sobre el pencil reducido no corrige el residual original |
| regla de Bathe | `q=max(p+8,2p)` supera el bloque mínimo `q=p` |
| guards globales | masa coordenadamente semidefinida, pérdida de rango `M` y rigidez singular |

Comando reproducible:

```bash
/home/pxpalacios/clark_kent/bin/python -m pytest -q test/test_cms_level2.py
```

Resultado del 2026-07-21:

```text
..............................                                           [100%]
30 passed in 2.62s
```

El informe se genera con:

```bash
/home/pxpalacios/clark_kent/bin/python test/cms_report.py
```

Para `n=120`, masa consistente y bases completas, el flujo de cuatro transformaciones
obtuvo:

| Métrica | Resultado |
|---|---:|
| error absoluto máximo de autovalor | `1.4211e-14` |
| residual máximo | `3.6553e-12` |
| salto máximo entre copias de interfaz | `0.0000e+00` |
| error de congruencia de `K` | `3.5527e-15` |
| error de congruencia de `M` | `3.5527e-15` |

Con truncación `k_L2=8`, `k_L1=10`, la dimensión final fue 43 frente a 120, y el error
relativo máximo de los primeros cinco autovalores fue `1.2057e-4`. El residual fue mucho
mayor (`7.1315e-2`), lo cual demuestra por qué el C++ no debe aceptar una retención sólo
por error aparente de autovalor: debe enriquecer usando el residual original.

## 9. Alcance de la evidencia

La batería demuestra:

- coherencia de las cuatro transformaciones;
- exactitud cuando las bases son completas;
- propiedad variacional bajo truncación;
- conservación de masa consistente;
- ensamblaje único y retro-sustitución coherente;
- comportamiento de guards matemáticos.

Por sí sola todavía no demuestra:

- parsing Tcl/Python;
- equivalencia del contribution store C++ con `addA/addM`;
- METIS/ParMETIS ni partición MPI real;
- Lanczos/MUMPS C++, LAPACK global o solve distribuido;
- ausencia de deadlocks;
- escalamiento o menor memoria en un modelo OpenSees real.

Las evidencias posteriores del 2026-07-21 completaron P1/P2: el espejo del repositorio
obtuvo 31 pruebas Python/C++ aprobadas, el refinador C++ pasó en `np=1,4`, el smoke Tcl
pasó en cuatro ranks y el edificio 1A cerró con residual `9.09852e-9`. Esos resultados no
provienen de este oráculo aislado; se documentan en el ADR y en el protocolo del edificio.
P3/P4 continúan siendo necesarios antes de afirmar escalabilidad o ventaja de rendimiento.
