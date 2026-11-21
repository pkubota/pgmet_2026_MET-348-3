
# Modelos práticos de assimilação de dados (MET-348-3, PGMET)

Experimento de gêmeos idênticos com advecção-difusão 1D periódica:
`dc/dt + (U + α c) dc/dx = κ d²c/dx²`. Verdade: U=1.0. Modelo: U=0.8 (erro de modelo)
e background inicial deslocado/mais fraco. Observações: 1 a cada 5 pontos, a cada 5 passos, σ_o=0,02.
Todas as rodadas usam as MESMAS observações (semente fixa).

| Programa | Método | Módulo das notas | Verificação embutida |
|---|---|---|---|
| `oi_advec_diff.f90`   | OI ciclada, B gaussiana não diagonal | II/III | χ² = dᵀS⁻¹d/nobs ≈ 1 |
| `var3d_advec_diff.f90`| 3D-Var (variável de controle Lv, CG) | II | igual à OI (dif. ~1e-16) |
| `kf_advec_diff.f90`   | Filtro de Kalman linear (P evolui) | III | — |
| `ekf_advec_diff.f90`  | EKF, modelo não linear (α=0,3) | IV | Jacobiana x dif. finitas |
| `enkf_advec_diff.f90` | EnKF estocástico, inflação + localização GC | IV | espalhamento x erro |
| `var4d_advec_diff.f90`| 4D-Var incremental, TL + adjunto, CG | II/IV | teste do gradiente |

`da_mod.f90`: modelo, Jacobiana, B/Q, observações, Cholesky, CG, saída GrADS. `da.nml`: parâmetros.

Uso: `make run` (ou `make MODE=opt run`), depois `make plot`.
Saídas por método: `<m>.bin/.ctl` (GrADS: livre, prev, verd, anal), `<m>_rmse.dat`.

## Melhorias em relação ao código de partida (advec_diff_OI_B.f90)
1. `Lcorr=0,1` com `dx=1` tornava B praticamente diagonal (off-diagonal ~ e⁻⁵⁰). Agora Lcorr em pontos de grade (3).
2. A análise não realimentava o modelo (`Run_Model` usa estado interno): não havia ciclo. Agora é cíclico.
3. Observações sem ruído (`y = c_true + (c_model-c_true)/2`) e `random_seed()` a cada passo. Agora ruído gaussiano, semente fixa.
4. Verdade analítica sem difusão vs. modelo com difusão. Agora verdade = mesmo modelo com parâmetros "verdadeiros".
5. `update(i)=u(i)` com `i` indefinido no ramo `monotone=.false.`; o limitador monotônico distorcia a análise. Removido.
6. Gauss-Jordan a cada passo, com B e R constantes: agora Cholesky fatorado uma vez.
7. Makefile: `-ffast-math` removido (quebra reprodutibilidade bit a bit); dependência do .mod corrigida.

## Limitações a discutir em aula
- OI/3D-Var/EnKF: B estática ou amostral; 4D-Var strong-constraint ignora o erro de modelo (por isso o background de 20 passos é pior).
- No 4D-Var a "análise" dentro da janela usa observações futuras (suavizador): não é comparação estritamente igual à dos filtros.
- KF roda com α=0 (modelo linear); os demais com α=0,3.

## Parâmetros ajustados (varredura de sensibilidade)
- `nwin = 10` (4D-Var): RMSE 0,0174; com 20 passos sobe a 0,0225 e com 100 a 0,0817 (erro de modelo ignorado no strong-constraint).
- `nens = 20` (EnKF): platô a partir de 20 membros (0,0249; 80 membros: 0,0243). Localização GC (raio 8) é essencial com ≤10 membros.
- `infl` quase não influencia enquanto `sig_q > 0` (Q já mantém o espalhamento).

## DFI (inicializacao por filtro digital) como opcao
Portada de `dfi_advdif.f90` (Dolph-Chebyshev, 2 estagios, leapfrog reversivel) para `da_mod.f90`
(`dfi_init`). Ligue em `da.nml`: `use_dfi=.true.` (+ `dfi_M`, `dfi_tau`, `dfi_diab`). Aplicada ao estado
analisado em cada metodo (EnKF: a cada membro; 4D-Var: ao estado que inicia o proximo ciclo; KF/EKF: so o
estado, P inalterada).
- Validacao: `make dfi_demo.exe && ./dfi_demo.exe` reproduz seu teste (RMSE t=0: 0,0519 / 0,0045 / 0,0046).
- Nos experimentos de assimilacao a DFI NAO melhora (OI 0,0243 -> 0,0256 a 0,0328, conforme tau). O pulso
  (sigma=3 pontos) ocupa a mesma banda de frequencia do ruido da analise e o filtro suaviza o sinal.
  Com pulso largo (`pulse_w=8`) fica neutro. Ela so compensa quando ha ruido de alta frequencia/desbalanco
  bem separado do sinal (caso do dfi_demo). Padrao: desligada.
- `pulse_w` (namelist) controla a largura do pulso verdadeiro.
- `pulse_shape` (namelist): 0 = gaussiana (padrao); 1 = onda quadrada de mesma variancia (largura L = 2*sqrt(3)*pulse_w). Os executaveis aceitam o ar
quivo de parametros como argumento: `./oi.exe outro.nml`.
- Correcao: `build_cov` agora usa a gaussiana periodizada (soma de imagens), exatamente definida positiva;
  antes, com lcorr grande, o Cholesky podia falhar. Resultados com lcorr=3 inalterados.

## Cenario em que a DFI ajuda (`make dfi_cenario`  ou  `./rodar_cenario_dfi.sh`)
`da_dfi.nml`: nx=200, pulso sigma=20 pontos, ruido de alta frequencia (ondas 8dx-3dx, `noise_amp=0.02`)
no background inicial, modelo quase perfeito (u_model=0.98), DFI com tau=9, M=6 (escala do teste original:
theta_s ~ 0,35 e sinal com theta <= ~0,06). Roda OI, KF, EnKF e 4D-Var com e sem DFI e gera `cenario_dfi.png`.
RMSE da analise sem DFI -> com DFI:
| | passo 5 | passo 10 |
|---|---|---|
| OI | 0,0265 -> 0,0110 | 0,0159 -> 0,0108 |
| KF | 0,0264 -> 0,0091 | 0,0158 -> 0,0103 |
| EnKF | 0,0272 -> 0,0112 | 0,0170 -> 0,0114 |
| 4D-Var | 0,0251 -> 0,0064 | 0,0137 -> 0,0066 |
O ganho e no transiente: depois o ruido decai sozinho (difusao + assimilacao) e as curvas se aproximam.
O custo e um leve alisamento do pico (~0,97 contra 1,0). Regra pratica: o sinal precisa ter theta bem abaixo
de theta_s (o filtro de Dolph nao tem banda passante plana); por isso nx=50 com sigma=3..8 nao serve.
Novas opcoes no namelist: `noise_amp` (ruido no background inicial), `dfi_init0` (DFI tambem no passo 0).


