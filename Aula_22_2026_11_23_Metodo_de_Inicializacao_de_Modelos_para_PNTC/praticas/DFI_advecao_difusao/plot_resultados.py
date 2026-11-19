#!/usr/bin/env python3
"""Compara saida_py x saida_f90 e gera figuras dos resultados."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dfi_advdif import P, resposta, grade

plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .3})
py = {n: np.loadtxt(f"saida_py/{n}.dat") for n in ("campos_t0", "rmse", "pesos", "campos_prev")}
f9 = {n: np.loadtxt(f"saida_f90/{n}.dat") for n in ("campos_t0", "rmse", "pesos", "campos_prev")}
print("Diferenca maxima |Fortran - Python|:")
for n in py:
    print(f"  {n:12s} {np.abs(py[n] - f9[n]).max():.3e}")

c = py["campos_t0"]; x = c[:, 0] / 1e3
fig, ax = plt.subplots(2, 2, figsize=(10, 7))
a = ax[0, 0]
a.plot(x, c[:, 2], c="0.65", lw=1, label="NIL (com ruído)")
a.plot(x, c[:, 1], "k--", lw=1, label="verdade")
a.plot(x, c[:, 3], c="C3", lw=1.4, label="DFI")
a.set_xlim(40, 200); a.set_xlabel("x (km)"); a.set_ylabel("c"); a.legend(fontsize=8)
a.set_title("(a) Campo em t = 0")

a = ax[0, 1]
r = py["rmse"]; t = r[:, 0] / 3600
a.semilogy(t, r[:, 1], c="0.5", label="NIL")
a.semilogy(t, r[:, 2], c="C3", label="DFI (adiabática)")
a.semilogy(t, r[:, 3], c="C0", ls="--", label="DFI (difusão no estágio 2)")
a.set_xlabel("tempo de previsão (h)"); a.set_ylabel("RMSE vs solução analítica"); a.legend(fontsize=8)
a.set_title("(b) Erro da previsão")

a = ax[1, 0]
M = P["M"]; k = py["pesos"][:, 0] * P["dt"] / 60
a.stem(k, py["pesos"][:, 1], basefmt=" ")
a.set_xlabel("tempo (min)"); a.set_ylabel("$h_k$"); a.set_title("(c) Pesos de Dolph ($M$=10)")

a = ax[1, 1]
th = np.linspace(0, np.pi, 2000); H = resposta(py["pesos"][:, 1], th)
a.plot(th, H ** 2, c="C3", label="$H^2$ (duas passagens)")
Cr = P["U"] * P["dt"] / (P["L"] / P["nx"])
ks = np.linspace(0, np.pi, 400)
th_k = np.arcsin(Cr * np.sin(ks))        # frequência numérica do leapfrog
sig_k = 2.5 / (P["sig0"] / (P["L"] / P["nx"]))
a.axvspan(0, np.arcsin(Cr * np.sin(sig_k)), color="C2", alpha=.25, label="banda do sinal")
nz = np.arcsin(Cr * np.sin(np.array([2 * np.pi / 8, 2 * np.pi / 3])))
a.axvspan(nz.min(), nz.max(), color="C1", alpha=.25, label="banda do ruído (8Δx–3Δx)")
a.set_xlim(0, 1.0); a.set_xlabel("$\\theta=\\omega\\Delta t$"); a.legend(fontsize=8)
a.set_title("(d) Resposta e bandas de frequência")
fig.tight_layout(); fig.savefig("resultados_dfi_advdif.png", dpi=150)
