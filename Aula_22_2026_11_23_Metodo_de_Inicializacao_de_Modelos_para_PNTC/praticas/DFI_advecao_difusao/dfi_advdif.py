#!/usr/bin/env python3
"""
Modelo 1D de advecção-difusão com Inicialização por Filtro Digital (DFI)
no estilo do HIRLAM (Lynch e Huang, 1992), usando o filtro de Dolph-Chebyshev.

    dc/dt + U dc/dx = nu d2c/dx2        (domínio periódico)

Esquema: leapfrog centrado no espaço; difusão defasada (lagged); primeiro
passo por Runge-Kutta de ponto médio. O leapfrog é reversível no tempo
(dt < 0), o que permite a integração para trás do estágio 1 da DFI.

DFI em dois estágios (como no HIRLAM):
  estágio 1: de t=0 a t=-TS, para trás, SEM difusão (processo irreversível),
             acumulando sum_k h_k c(t_k) -> campo cbar válido em t=-TS/2
  estágio 2: de t=-TS/2 a t=+TS/2, para frente, partindo de cbar,
             acumulando os mesmos pesos -> campo inicializado válido em t=0
Opção diab=True liga a difusão no estágio 2 (variante "diabática").

Uso:  python3 dfi_advdif.py [--outdir saida_py]
"""
import argparse, os
import numpy as np

# ----------------------------------------------------------------- parâmetros
P = dict(
    nx=400, L=400.0e3,        # pontos de grade, comprimento do domínio (m)
    U=10.0, nu=100.0,         # velocidade (m/s), difusividade (m2/s)
    dt=50.0,                  # passo de tempo (s) -> Courant = U dt/dx = 0.5
    tau_s=900.0,              # período de corte do filtro (s)
    M=10,                     # N = 2M+1 pesos; TS = 2 M dt = 1000 s
    x0c=120.0e3, sig0=20.0e3, # centro e largura da gaussiana inicial (m)
    nfc=216,                  # passos da previsão (216*50 s = 3 h)
    noise_amp=0.008,          # amplitude por modo do ruído
)

# ------------------------------------------------------ filtro de Dolph-Chebyshev
def cheb(n, x):
    """Polinômio de Chebyshev T_n(x), válido para |x|<=1 e |x|>1 (x>0)."""
    x = np.asarray(x, float)
    out = np.empty_like(x)
    m = np.abs(x) <= 1.0
    out[m] = np.cos(n * np.arccos(x[m]))
    out[~m] = np.cosh(n * np.arccosh(x[~m]))
    return out

def dolph(M, theta_s):
    """Pesos h_k, k=-M..M (retorna vetor de 2M+1, índice s=k+M), x0 e r."""
    N = 2 * M + 1
    x0 = 1.0 / np.cos(theta_s / 2.0)
    r = 1.0 / np.cosh(2 * M * np.arccosh(x0))
    m = np.arange(1, M + 1)
    th_m = 2.0 * np.pi * m / N
    k = np.arange(-M, M + 1)
    h = (1.0 + 2.0 * r * (cheb(2 * M, x0 * np.cos(th_m / 2.0))[None, :]
         * np.cos(np.outer(k, th_m))).sum(axis=1)) / N
    return h, x0, r

def resposta(h, theta):
    M = (len(h) - 1) // 2
    k = np.arange(-M, M + 1)
    return (h[None, :] * np.exp(-1j * np.outer(theta, k))).sum(axis=1).real

# ------------------------------------------------------------------ modelo
def grade(p):
    dx = p["L"] / p["nx"]
    x = (np.arange(p["nx"]) + 0.5) * dx
    return dx, x

def analitica(p, x, t):
    """Gaussiana advectada e difundida (imagens periódicas -1,0,+1)."""
    s2 = p["sig0"] ** 2 + 2.0 * p["nu"] * t
    amp = p["sig0"] / np.sqrt(s2)
    c = np.zeros_like(x)
    for j in (-1, 0, 1):
        xc = p["x0c"] + p["U"] * t + j * p["L"]
        c += amp * np.exp(-0.5 * (x - xc) ** 2 / s2)
    return c

def ruido(p, x):
    """Ruído determinístico (fases por razão áurea): ondas de 8dx a 3dx."""
    nx = p["nx"]
    ms = np.arange(nx // 8, nx // 3 + 1)
    phi = 2.0 * np.pi * np.mod(ms * 0.6180339887498949, 1.0)
    n = np.zeros_like(x)
    for m, ph in zip(ms, phi):
        n += p["noise_amp"] * np.cos(2.0 * np.pi * m * x / p["L"] + ph)
    return n

def tend(c, p, dx, adv=True, dif=True):
    cp, cm = np.roll(c, -1), np.roll(c, 1)
    f = np.zeros_like(c)
    if adv: f += -p["U"] * (cp - cm) / (2.0 * dx)
    if dif: f += p["nu"] * (cp - 2.0 * c + cm) / dx ** 2
    return f

def primeiro_passo(c, dt, p, dx, dif):
    """Runge-Kutta de ponto médio (mantém a amplitude, evita Euler instável)."""
    cs = c + 0.5 * dt * tend(c, p, dx, True, dif)
    return c + dt * tend(cs, p, dx, True, dif)

def integra(c0, dt, nsteps, p, dx, dif, hook=None):
    """Leapfrog: devolve a lista [c0, c1, ..., c_nsteps]. hook(s, c) é chamado
    em cada instante (s=0..nsteps)."""
    cold = c0.copy()
    if hook: hook(0, cold)
    if nsteps == 0: return cold
    cnow = primeiro_passo(cold, dt, p, dx, dif)
    if hook: hook(1, cnow)
    for s in range(2, nsteps + 1):
        cnew = cold + 2.0 * dt * tend(cnow, p, dx, True, False) \
                    + (2.0 * dt * tend(cold, p, dx, False, True) if dif else 0.0)
        cold, cnow = cnow, cnew
        if hook: hook(s, cnow)
    return cnow

# ------------------------------------------------------------------ DFI
def dfi(c0, p, dx, h, diab=False):
    M = p["M"]; dt = p["dt"]; ns = 2 * M
    # estágio 1: para trás (dt<0), sem difusão; pesos h[s] no passo s (s=0..2M)
    acc = np.zeros_like(c0)
    def hook1(s, c):
        nonlocal acc; acc += h[s] * c
    integra(c0, -dt, ns, p, dx, False, hook1)
    cbar = acc.copy()                       # válido em t = -TS/2
    # estágio 2: para frente, a partir de cbar
    acc = np.zeros_like(c0)
    def hook2(s, c):
        nonlocal acc; acc += h[s] * c
    integra(cbar, +dt, ns, p, dx, diab, hook2)
    return acc                               # válido em t = 0

def previsao(c0, p, dx, x):
    """Previsão com difusão; devolve RMSE(t) contra a analítica e a série N1."""
    n = p["nfc"]; dt = p["dt"]
    rmse = np.zeros(n + 1)
    def hook(s, c):
        rmse[s] = np.sqrt(np.mean((c - analitica(p, x, s * dt)) ** 2))
    cf = integra(c0, dt, n, p, dx, True, hook)
    return rmse, cf

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="saida_py")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    p = P
    dx, x = grade(p)
    theta_s = 2.0 * np.pi * p["dt"] / p["tau_s"]
    h, x0, r = dolph(p["M"], theta_s)
    print(f"Courant={p['U']*p['dt']/dx:.3f}  theta_s={theta_s:.5f}  M={p['M']}  "
          f"N={2*p['M']+1}  TS={2*p['M']*p['dt']:.0f} s")
    print(f"x0={x0:.6f}  r={r:.5f} ({20*np.log10(r):.2f} dB)  soma(h)={h.sum():.12f}  "
          f"h_min={h.min():.4e}")

    ctrue = analitica(p, x, 0.0)
    cnil = ctrue + ruido(p, x)                       # análise "não inicializada"
    cdfi = dfi(cnil, p, dx, h, diab=False)           # DFI adiabática (HIRLAM)
    cdfd = dfi(cnil, p, dx, h, diab=True)            # variante com difusão no estágio 2

    r_nil, f_nil = previsao(cnil, p, dx, x)
    r_dfi, f_dfi = previsao(cdfi, p, dx, x)
    r_dfd, f_dfd = previsao(cdfd, p, dx, x)
    rm = lambda c: np.sqrt(np.mean((c - ctrue) ** 2))
    print(f"RMSE em t=0:  NIL={rm(cnil):.5f}  DFI={rm(cdfi):.5f}  DFI-diab={rm(cdfd):.5f}")
    for k in (0, 36, 108, 216):
        print(f"RMSE t={k*p['dt']/3600:4.1f} h:  NIL={r_nil[k]:.5f}  DFI={r_dfi[k]:.5f}  "
              f"DFI-diab={r_dfd[k]:.5f}")

    np.savetxt(f"{a.outdir}/campos_t0.dat",
               np.column_stack([x, ctrue, cnil, cdfi, cdfd]),
               header="x ctrue cnil cdfi cdfi_diab", fmt="%.15e")
    np.savetxt(f"{a.outdir}/rmse.dat",
               np.column_stack([np.arange(p["nfc"] + 1) * p["dt"], r_nil, r_dfi, r_dfd]),
               header="t_s rmse_nil rmse_dfi rmse_dfi_diab", fmt="%.15e")
    np.savetxt(f"{a.outdir}/pesos.dat", np.column_stack([np.arange(-p["M"], p["M"] + 1), h]),
               header="k h_k", fmt="%.15e")
    np.savetxt(f"{a.outdir}/campos_prev.dat",
               np.column_stack([x, f_nil, f_dfi, f_dfd, analitica(p, x, p['nfc']*p['dt'])]),
               header="x c_nil c_dfi c_dfi_diab c_true (fim da previsao)", fmt="%.15e")

if __name__ == "__main__":
    main()
