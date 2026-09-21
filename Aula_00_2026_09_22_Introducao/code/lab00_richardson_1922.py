"""
notebooks/lab00_richardson_1922.py -- MET-348-3, Aula prática 01
======================================================================
Recriando (em espírito) o experimento de L. F. Richardson (1922): a
primeira tentativa científica de integrar numericamente as equações da
atmosfera -- e por que ela falhou de forma catastrófica (tendência de
pressão prevista: 146 hPa em 6h, quando o valor real observado foi de
apenas ~1 hPa).

A causa raiz (explicada décadas depois por Peter Lynch, citado na aula
histórica): os dados de vento disponíveis a Richardson, estimados à mão a
partir de observações esparsas, tinham uma pequena componente divergente
espúria (erro de análise/diferenciação). Como a tendência de pressão vem
da DIVERGÊNCIA do vento (equação da continuidade), um erro pequeno no
vento produz um erro ENORME na tendência de pressão -- ondas de
gravidade/som não-filtradas "explodem" no cálculo.

Este notebook reproduz esse mecanismo, na mesma grade que Richardson usou
(~200 km), e mostra a solução: um "filtro" (aqui, suavização espacial do
campo de vento -- o mesmo espírito da inicialização por modos normais de
Lynch, que retira a componente de onda de gravidade não-balanceada antes
de integrar) traz a tendência de volta a valores realistas.

Rode como script ou copie para um notebook Jupyter/Colab.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

# ---- grade: a mesma resolução do experimento original de Richardson -----
nx, ny = 64, 64
dx = dy = 200e3  # 200 km, como no experimento de 1922
x = np.arange(nx) * dx
y = np.arange(ny) * dy
X, Y = np.meshgrid(x, y)


def divergence(u, v, dx, dy):
    return np.gradient(u, dx, axis=1) + np.gradient(v, dy, axis=0)


def pressure_tendency_hPa_per_6h(div, ps=1000.0):
    """Tendência de pressão de superfície, proxy simplificado da equação
    da continuidade integrada na coluna (dp/dt ~ -ps * divergência),
    escalada para hPa por 6 horas -- suficiente para reproduzir a ORDEM
    DE GRANDEZA do erro de Richardson, não os detalhes de sua análise
    multi-camada."""
    return -ps * div * 3600 * 6


# ---- escoamento sinótico "verdadeiro": suave e balanceado (não-divergente) ----
U0 = 8.0
k = 2 * np.pi / 3000e3   # comprimento de onda sinótico ~3000 km
l = 2 * np.pi / 2500e3
amp_psi = 3.0e6
psi_true = -U0 * Y + amp_psi * np.cos(k * X) * np.sin(l * Y)
u_true = -np.gradient(psi_true, dy, axis=0)
v_true = np.gradient(psi_true, dx, axis=1)

# ---- "dados de Richardson": o mesmo escoamento + um pequeno erro divergente ----
# (representa o erro de estimar/diferenciar à mão um campo de vento a
# partir de observações esparsas -- um erro médio de menos de 1 m/s basta)
rng = np.random.default_rng(0)
chi_noise = gaussian_filter(rng.normal(size=(ny, nx)), sigma=1.0)
chi = 3e5 * chi_noise
u_chi = np.gradient(chi, dx, axis=1)
v_chi = np.gradient(chi, dy, axis=0)
u_noisy, v_noisy = u_true + u_chi, v_true + v_chi
print(f"Erro médio no vento introduzido: {np.mean(np.sqrt(u_chi**2+v_chi**2)):.2f} m/s "
      f"(máximo: {np.max(np.sqrt(u_chi**2+v_chi**2)):.2f} m/s)")

# ---- "com filtro": suaviza o campo de vento antes de calcular a tendência ----
sigma_filter = 7
u_filtered = gaussian_filter(u_noisy, sigma=sigma_filter)
v_filtered = gaussian_filter(v_noisy, sigma=sigma_filter)

tend_true = pressure_tendency_hPa_per_6h(divergence(u_true, v_true, dx, dy))
tend_noisy = pressure_tendency_hPa_per_6h(divergence(u_noisy, v_noisy, dx, dy))
tend_filtered = pressure_tendency_hPa_per_6h(divergence(u_filtered, v_filtered, dx, dy))

print(f"\nTendência de pressão implícita (equação da continuidade):")
print(f"  verdade (balanceada):        max|.| = {np.abs(tend_true).max():7.2f} hPa/6h")
print(f"  SEM filtro (tipo Richardson): max|.| = {np.abs(tend_noisy).max():7.2f} hPa/6h")
print(f"  COM filtro (suavizado):       max|.| = {np.abs(tend_filtered).max():7.2f} hPa/6h")
print(f"\n  (para comparação: o valor real relatado por Richardson foi de 146 hPa/6h;")
print(f"   o valor observado de verdade naquele dia foi de apenas ~1 hPa/6h)")

# ---- Figura -----------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
vmax_wind = max(np.abs(u_noisy).max(), np.abs(u_true).max())
for ax, u, v, title in zip(axes[0], [u_true, u_noisy, u_filtered], [v_true, v_noisy, v_filtered],
                            ["Escoamento 'verdadeiro'\n(balanceado)",
                             "SEM filtro\n(dados 'à la Richardson')",
                             "COM filtro\n(suavizado, tipo Lynch)"]):
    skip = 4
    ax.quiver(X[::skip, ::skip] / 1e3, Y[::skip, ::skip] / 1e3,
              u[::skip, ::skip], v[::skip, ::skip], scale=200)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("x (km)")
axes[0, 0].set_ylabel("y (km)")

vmax_tend = max(np.abs(tend_noisy).max(), 1.0)
for ax, tend, title in zip(axes[1], [tend_true, tend_noisy, tend_filtered],
                            [f"máx |tendência| = {np.abs(tend_true).max():.1f} hPa/6h",
                             f"máx |tendência| = {np.abs(tend_noisy).max():.1f} hPa/6h  (Richardson: 146!)",
                             f"máx |tendência| = {np.abs(tend_filtered).max():.1f} hPa/6h"]):
    im = ax.pcolormesh(X / 1e3, Y / 1e3, tend, cmap="RdBu_r", vmin=-vmax_tend, vmax=vmax_tend, shading="auto")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("x (km)")
axes[1, 0].set_ylabel("y (km)")
fig.colorbar(im, ax=axes[1, :].tolist(), label="tendência de pressão (hPa/6h)", shrink=0.8)

fig.suptitle("Recriando o experimento de Richardson (1922): por que a primeira previsão numérica falhou")
fig.savefig("lab00_richardson_1922.png", dpi=120, bbox_inches="tight")
print("\nFigura salva: lab00_richardson_1922.png")
