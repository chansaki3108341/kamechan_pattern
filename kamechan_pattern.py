from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def make_kusagame_cmap(colors=("#121411", "#4A5237", "#F2E58C")):
    return mcolors.LinearSegmentedColormap.from_list("kusagame", list(colors))

@dataclass(frozen=True)
class GrayScottParams:
    Du: float = 0.16
    Dv: float = 0.08
    F: float = 0.045
    k: float = 0.060


@dataclass(frozen=True)
class InitParams:
    N: int = 100
    patches: int = 20
    patch_size: int = 5
    patch_U: float = 0.50
    patch_V: float = 0.25
    noise: float = 0.05
    U0: float = 1.0
    V0: float = 0.0

def laplacian_eigs_rfft(N: int) -> np.ndarray:
    kx = 2.0 * np.pi * np.fft.fftfreq(N)        
    ky = 2.0 * np.pi * np.fft.rfftfreq(N) 
    lam = (2.0 * np.cos(kx)[:, None] + 2.0 * np.cos(ky)[None, :] - 4.0)
    return lam 


def init_uv(rng: np.random.Generator, init: InitParams, dtype=np.float32):
    N = init.N
    U = np.full((N, N), init.U0, dtype=dtype)
    V = np.full((N, N), init.V0, dtype=dtype)

    s = init.patch_size
    for _ in range(init.patches):
        x, y = rng.integers(0, N - s, size=2)
        xs, ys = slice(x, x + s), slice(y, y + s)
        U[xs, ys] = init.patch_U
        V[xs, ys] = init.patch_V

    if init.noise:
        U += init.noise * rng.random((N, N), dtype=dtype)
        V += init.noise * rng.random((N, N), dtype=dtype)

    return U, V


def simulate_gray_scott_imex_fft(
    seed: int,
    *,
    total_time: float = 15000.0,
    dt: float = 5.0,  
    params: GrayScottParams = GrayScottParams(),
    init: InitParams = InitParams(),
    dtype=np.float32,  
):
    rng = np.random.default_rng(seed)
    U, V = init_uv(rng, init, dtype=dtype)

    N = init.N
    steps = int(total_time / dt)

    lam = laplacian_eigs_rfft(N).astype(np.float64, copy=False)  # 精度安定のため係数側はfloat64でもOK
    denom_u = 1.0 - dt * params.Du * lam
    denom_v = 1.0 - dt * params.Dv * lam

    F = params.F
    Fk = params.F + params.k

    for _ in range(steps):
        uvv = U * V * V
        Ru = -uvv + F * (1.0 - U)
        Rv =  uvv - Fk * V

        rhs_u = U + dt * Ru
        rhs_v = V + dt * Rv

        U_hat = np.fft.rfft2(rhs_u)
        V_hat = np.fft.rfft2(rhs_v)

        U = np.fft.irfft2(U_hat / denom_u, s=(N, N)).astype(dtype, copy=False)
        V = np.fft.irfft2(V_hat / denom_v, s=(N, N)).astype(dtype, copy=False)

        if not np.isfinite(U).all() or not np.isfinite(V).all():
            raise FloatingPointError("Numerical blow-up detected. Try smaller dt.")

    return U, V


def plot_gallery(
    seeds: np.ndarray,
    *,
    rows: int = 3,
    cols: int = 3,
    total_time: float = 15000.0,
    dt: float = 5.0,
):
    cmap = make_kusagame_cmap()
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
    axes = np.asarray(axes).reshape(rows, cols)

    for i, (ax, seed) in enumerate(zip(axes.flat, seeds), start=1):
        print(f"{i}枚目を計算中... (Seed: {seed})🐢")
        U, _ = simulate_gray_scott_imex_fft(seed, total_time=total_time, dt=dt)

        ax.imshow(U, cmap=cmap)
        ax.set_title(f"Seed: {seed}")
        ax.axis("off")

    print("完成しました！🐢")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    master_seed = 20260220  
    rng = np.random.default_rng(master_seed)
    seeds = rng.integers(0, 10000, size=9)

    plot_gallery(seeds, total_time=15000.0, dt=5.0)
