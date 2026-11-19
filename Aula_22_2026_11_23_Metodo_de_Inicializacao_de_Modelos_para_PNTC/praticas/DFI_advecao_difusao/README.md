# DFI em um modelo 1D de adveccao-difusao (estilo HIRLAM) 

Equacao: `dc/dt + U dc/dx = nu d2c/dx2`, dominio periodico.
Versoes equivalentes em Python (`dfi_advdif.py`) e Fortran 90 (`dfi_advdif.f90`).

    make            # compila, roda as duas versoes e gera a figura
    # ou: python3 dfi_advdif.py ; gfortran -O2 -o dfi_advdif dfi_advdif.f90 && ./dfi_advdif

## O que o codigo faz
- **Modelo:** leapfrog centrado (Courant 0,5), difusÃ£o defasada, 1Âº passo por Runge-Kutta de ponto mÃ©dio.
  O leapfrog eh reversivel no tempo, entao aceita `dt < 0`.
- **Filtro:** Dolph-Chebyshev, pesos pela formula fechada (M = 10, N = 21, ts = 900 s, dt = 50 s,
  TS = 2MDt = 1000 s, r = 0,060). Pesos todos positivos, soma 1.
- **DFI em dois estagios (como no HIRLAM):**
  1. de t = 0 a -TS, para tras, **sem difusao**, acumulando soma h_k c(t_k) -> `cbar` valido em -TS/2;
  2. de -TS/2 a +TS/2, para frente, a partir de `cbar`, acumulando os mesmos pesos -> campo
     inicializado, valido em t = 0.
  A variante `diab` liga a difusao no estagio 2.
- **Teste:** gaussiana (sigma = 20 km) advectada a 10 m/s + ruido deterministico de ondas de 8dx a 3dx
  (fases por razao aurea, identico nas duas linguagens). A solucao analitica (com imagens
  periodicas) permite medir o erro.

## Resultados (valores identicos em Fortran e Python, diferenca maxima ~1e-15)
| RMSE vs analitica | t = 0 | t = 3 h |
|---|---|---|
| NIL (com ruido) | 0,0519 | 0,0127 |
| DFI  adiabatica | 0,0045 | 0,0046 |
| DFI, difusao no estagio 2 | 0,0046 | 0,0046 |

O NIL so melhora com o tempo porque a difusao amortece o ruido; a DFI o remove ja em t = 0.
O erro residual da DFI (0,0045) eh o alisamento do proprio sinal pelo filtro.

## Cuidados ao reutilizar
- **Separacao em frequencia.** O filtro age no tempo, ponto a ponto. Para um campo advectado,
  isso vira um filtro em numero de onda, com omega = U k. O sinal precisa estar bem abaixo do corte
  (aqui theta¸ <= 0,06 contra theta_s = 0,35). Com uma gaussiana de sigma = 8dx o filtro suavizava
  o proprio
  sinal, e foi preciso alarga-la.
- **Onda de 2dx.** No esquema centrado ela eh estacionaria (omega = 0) e **nao ehh removida** por filtro
  temporal. O ruio do teste usa ondas de 8dx a 3dx.
- **Courant baixo** empurra o ruido de grade para theta¸ pequeno; verifique tehta_s < arcsin(C sin kdx)`
  para o ruido que quer remover.
- Os parametros foram ajustados para este teste; para outro caso refaca a analise (figura, painel d).
