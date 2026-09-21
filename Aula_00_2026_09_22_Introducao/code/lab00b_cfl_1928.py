"""
notebooks/lab00b_cfl_1928.py -- MET-348-3, Aula prática 01 (complemento)
=============================================================================
A condição de Courant-Friedrichs-Lewy (1928), citada na aula histórica como
um dos pilares teóricos que tornaram a PNT viável: um esquema numérico
explícito só é estável se a informação não "viajar" mais que um ponto de
grade por passo de tempo.

    C = U*dt/dx <= 1   (numero de Courant)

Demonstração mínima: a mesma equação de advecção 1D, resolvida com C<1
(estável) e C>1 (instável).
"""
import numpy as np
import matplotlib.pyplot as plt

nx = 100
dx = 1.0
U = 1.0
x = np.arange(nx) * dx
c0 = np.exp(-0.5 * ((x - 20) / 3) ** 2)  # pulso gaussiano inicial


def advect_upwind(c0, dt, n_steps):
    c = c0.copy()
    hist = [c.copy()]
    for _ in range(n_steps):
        c = c - (U * dt / dx) * (c - np.roll(c, 1))
        hist.append(c.copy())
    return np.array(hist)


fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for ax, dt, label in zip(axes, [0.8, 1.3], ["C = 0.8 (estável)", "C = 1.3 (instável)"]):
    n_steps = 60
    hist = advect_upwind(c0, dt, n_steps)
    for i in range(0, n_steps + 1, 15):
        ax.plot(x, hist[i], label=f"passo {i}")
    ax.set_title(f"{label}\nC = U·Δt/Δx = {U*dt/dx:.2f}")
    ax.set_xlabel("x"); ax.legend(fontsize=8)
    ax.set_ylim(-2, 2) if dt > 1.0 else None

fig.suptitle("Condição de CFL (Courant, Friedrichs e Lewy, 1928)")
fig.tight_layout()
fig.savefig("lab00b_cfl.png", dpi=120)
print("Figura salva: lab00b_cfl.png")

hist_stable = advect_upwind(c0, 0.8, 60)
hist_unstable = advect_upwind(c0, 1.3, 60)
print(f"\nC=0.8: amplitude máxima ao final = {np.abs(hist_stable[-1]).max():.2f} (esperado: ~1, preservada)")
print(f"C=1.3: amplitude máxima ao final = {np.abs(hist_unstable[-1]).max():.2e} (deveria estar explodindo)")
