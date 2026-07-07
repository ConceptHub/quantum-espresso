# Data Flow and Computation Graph — `dv_vdW_DF.f90`

## Call Graph

```
dv_drho_vdwdf(rho, drho, nspin, q_point, dv_drho)
└── get_delta_v(rho, drho, nspin, q_point, delta_v)
    ├── fft_gradient_r2r(...)           builds gradient_rho(3,nnr)
    ├── fft_qgradient(...)              builds gradient_drho(3,nnr)
    ├── fill_q0_extended_on_grid()      builds q0, dq0_dq, d2q0_dq2,
    │                                   dq_dn_n, dn_dq_dn_n_n, dq_dgradn_n_gmod
    ├── initialize_spline_interpolation(q_mesh, d2y_dx2)   [once, saved]
    ├── [Loop 1]  get_abcdef + get_thetas_exentended × (nnr × Nqs)
    │   └── builds b1, b2, u, delta_u, dtheta_dgradn_save, h1part2_save
    ├── get_u_delta_u(u, delta_u, q_point)
    │   ├── fwfft × Nqs                 u, delta_u → G-space
    │   ├── interpolate_kernel × ngm    kernel_of_g(Nqs,Nqs) per G-vector
    │   └── invfft × Nqs                back to real space
    ├── [Accum]   delta_v += delta_u·b1 + u·b2
    ├── [Loop 2 — arithmetic only, no kernel calls after refactor]
    │   └── h1t, h2t from u·h1part2_save + delta_u·dtheta_dgradn_save
    └── [icar=1,2,3]  divergence step
        ├── fwfft(delta_h)
        ├── multiply by i(G+q)
        └── invfft → delta_v -= tpiba·delta_h_aux
```

---

## Array Inventory

| Array | Shape | Type | Produced by | Consumed by |
|-------|-------|------|-------------|-------------|
| `total_rho` | `(nnr)` | real | `get_delta_v` (copy of rho) | `fill_q0_extended_on_grid`, Loop 1 |
| `gradient_rho` | `(3,nnr)` | real | `fft_gradient_r2r` | Loop 1, fill_q0, Loop 2 divergence |
| `gradient_drho` | `(3,nnr)` | complex | `fft_qgradient` | `get_thetas_exentended`, Loop 2 divergence |
| `q0` | `(nnr)` | real | `fill_q0_extended_on_grid` | `get_abcdef` (binary search) |
| `q` | `(nnr)` | real | `fill_q0_extended_on_grid` | (diagnostic only) |
| `dq0_dq` | `(nnr)` | real | `fill_q0_extended_on_grid` | `get_thetas_exentended` |
| `d2q0_dq2` | `(nnr)` | real | `fill_q0_extended_on_grid` | `get_thetas_exentended` |
| `dq_dn_n` | `(nnr)` | real | `fill_q0_extended_on_grid` | `get_thetas_exentended` |
| `dn_dq_dn_n_n` | `(nnr)` | real | `fill_q0_extended_on_grid` | `get_thetas_exentended` |
| `dq_dgradn_n_gmod` | `(nnr)` | real | `fill_q0_extended_on_grid` | `get_thetas_exentended` |
| `d2y_dx2` | `(Nqs,Nqs)` | real | `initialize_spline_interpolation` | `get_thetas_exentended` (saved across calls) |
| `b1` | `(nnr,Nqs)` | real | Loop 1 (`dtheta_dn`) | Accum: `delta_v += delta_u·b1` |
| `b2` | `(nnr,Nqs)` | complex | Loop 1 (`d2theta_dn2…`) | Accum: `delta_v += u·b2` |
| `u` | `(nnr,Nqs)` | complex | Loop 1 (`theta`) | `get_u_delta_u` (in-place), Accum, Loop 2 |
| `delta_u` | `(nnr,Nqs)` | complex | Loop 1 (`dtheta_dn·drho + dtheta_dgradn·∇n·∇δn`) | `get_u_delta_u` (in-place), Accum, Loop 2 |
| `dtheta_dgradn_save` | `(nnr,Nqs)` | real | Loop 1 (saved from `get_thetas_exentended`) | Loop 2 |
| `h1part2_save` | `(nnr,Nqs)` | complex | Loop 1 (preformed `dn_dtheta_dgradn·(δn/n) + dgradn_dtheta_dgradn·(∇n·∇δn/n)`) | Loop 2 |
| `h1t` | `(nnr)` | complex | Loop 2 | Divergence step |
| `h2t` | `(nnr)` | complex | Loop 2 | Divergence step |
| `delta_h` | `(nnr)` | complex | Divergence step | Divergence step (FFT work buffer) |
| `delta_h_aux` | `(nnr)` | complex | Divergence step | subtracted from `delta_v` |
| `delta_v` | `(nnr)` | complex | Accum + Divergence step | returned to caller |

---

## Data Flow Diagram — `get_delta_v`

```
rho(:,1) ──────────────────────────────────────────── total_rho(nnr)
drho(:,1) ─────────────────────────────────────────── (used inline)
                                                              │
                          ┌───────────────────────────────────┘
                          │
          ┌───────────────▼───────────────┐
          │  fft_gradient_r2r             │
          │  fft_qgradient                │
          └───────┬───────────────┬───────┘
                  │               │
          gradient_rho(3,nnr)  gradient_drho(3,nnr)
                  │               │
          ┌───────▼───────────────┘
          │  fill_q0_extended_on_grid
          └───────┬──────────────────────────────────────────────────┐
                  │                                                  │
      q0, dq0_dq, d2q0_dq2,                                         │
      dq_dn_n, dn_dq_dn_n_n,                                        │
      dq_dgradn_n_gmod (all nnr)                                     │
                  │                                                  │
          ┌───────▼────────────────────────────────────────────────┐ │
          │         LOOP 1  — O(nnr × Nqs)                        │ │
          │                                                        │ │
          │  get_abcdef(q0, i_grid)  ─► q_hi, q_low, dq, a…f     │ │
          │       │                                                │ │
          │  get_thetas_exentended(…, P_i, i_grid)                │ │
          │       │  reads:  d2y_dx2(Nqs,Nqs)  ◄── (saved)       │ │
          │       │  reads:  gradient_rho, gradient_drho ─────────┼─┘
          │       │  reads:  dq0_dq, d2q0_dq2, dq_dn_n, …        │
          │       │                                                │
          │       ├─► theta         ──► u(i_grid, P_i)            │
          │       ├─► dtheta_dn     ──► b1(i_grid, P_i)           │
          │       │                     delta_u(i_grid, P_i)      │
          │       ├─► dtheta_dgradn ──► delta_u(i_grid, P_i)      │
          │       │                     dtheta_dgradn_save(i,P)   │
          │       ├─► d2theta_dn2   ──► b2(i_grid, P_i)           │
          │       ├─► dn_dtheta_dgradn ─► b2(i_grid, P_i)         │
          │       │                      h1part2_save(i,P)        │
          │       └─► dgradn_dtheta_dgradn ─► h1part2_save(i,P)   │
          └────────────────────────────────────────────────────────┘
                  │
                  ▼
          ┌───────────────────────────────────────────────────────┐
          │  get_u_delta_u                                        │
          │                                                       │
          │  fwfft(u, delta_u) × Nqs        (real → G-space)     │
          │  ∀ G:  interpolate_kernel(|G|)  → kernel_of_g        │
          │        interpolate_kernel(|G+q|) → kernel_of_gq      │
          │        u(G,α)      = Σ_β kernel_of_g(α,β)  u(G,β)   │
          │        δu(G,α)     = Σ_β kernel_of_gq(α,β) δu(G,β)  │
          │  invfft(u, delta_u) × Nqs       (G-space → real)     │
          └───────┬───────────────────────────────────────────────┘
                  │   u, delta_u now contain kernel-convoluted quantities
                  │
          ┌───────▼──────────────────────────────────────────┐
          │  ACCUMULATION  (pure arithmetic, nnr × Nqs)      │
          │                                                   │
          │  delta_v(i) += Σ_P  delta_u(i,P)·b1(i,P)        │
          │                    + u(i,P)·b2(i,P)              │
          └───────┬──────────────────────────────────────────┘
                  │
          ┌───────▼──────────────────────────────────────────┐
          │  LOOP 2  (arithmetic only, nnr × Nqs)            │
          │  [no get_abcdef / get_thetas_exentended calls]   │
          │                                                   │
          │  h1t(i) += Σ_P  delta_u(i,P)·dtheta_dgradn_save(i,P)   │
          │                + u(i,P)·h1part2_save(i,P)               │
          │  h2t(i) += Σ_P  u(i,P)·dtheta_dgradn_save(i,P)         │
          └───────┬──────────────────────────────────────────┘
                  │
          ┌───────▼────────────────────────────────────────────────────┐
          │  DIVERGENCE STEP  (icar = 1, 2, 3)                        │
          │                                                            │
          │  delta_h = h1t·∇n[icar] + h2t·∇δn[icar]                  │
          │  delta_h → G-space  (fwfft)                               │
          │  delta_h_aux(G) = i(G[icar]+q[icar]) · delta_h(G)        │
          │  delta_h_aux → real space  (invfft)                       │
          │  delta_v -= tpiba · delta_h_aux                           │
          └───────┬────────────────────────────────────────────────────┘
                  │
                  ▼
          delta_v(nnr)  [returned × e²]
```

---

## `fill_q0_extended_on_grid` — Internal Flow

```
total_rho(i), gradient_rho(:,i)
        │
        ├─► kF = (3π²ρ)^(1/3)
        ├─► r_s = (3/4πρ)^(1/3),  sqrt_r_s
        ├─► gc  = -Z_ab/(36 kF ρ²) |∇ρ|²
        ├─► gmod = |∇ρ|
        ├─► LDA_1, LDA_2  (Perdew-Wang LDA correlation pieces)
        │
        ├─► q(i) = kF + LDA_1·log(1+1/LDA_2) + gc    [DION eq. 11–12]
        │
        ├─► inner loop index=1..m_cut(=12):            [SOLER eq. 7]
        │       exponent  += (q/q_cut)^index / index
        │       dq0_dq    += (q/q_cut)^(index-1)
        │       d2q0_dq2  accumulates second-derivative terms
        │
        ├─► q0(i)     = q_cut (1 − exp(−exponent))
        ├─► dq0_dq(i) *= exp(−exponent)
        ├─► d2q0_dq2(i) assembled from expTemp2 and dq0_dq
        │
        ├─► dLDA_1_dn_n, dLDA_2_dn_n, d2LDA_1_dn2_n2, d2LDA_2_dn2_n2
        │
        ├─► dq_dn_n(i)          = ∂(ρq)/∂ρ / ρ   (chain rule)
        ├─► dn_dq_dn_n_n(i)     = ρ ∂²(ρq)/∂ρ² / ρ
        └─► dq_dgradn_n_gmod(i) = -Z_ab / (18 kF ρ)
```

---

## `get_thetas_exentended` — Spline Evaluation per (i_grid, P_i)

```
Inputs: q_hi, q_low, dq, a, b, c, d, e, f   [from get_abcdef]
        P_i                                   [q-mesh index]
        d2y_dx2(Nqs, Nqs)                     [spline 2nd derivatives, saved]
        gradient_rho, gradient_drho           [module arrays]
        dq0_dq, d2q0_dq2, dq_dn_n, …         [module arrays]

  y_qlow = δ(q_low, P_i),   y_qhi = δ(q_hi, P_i)   [unit-vector trick]

  d2P_dq02 = a·d2y_dx2(P_i,q_low) + b·d2y_dx2(P_i,q_hi)
  dP_dq0   = (y_qhi−y_qlow)/dq − e·d2y_dx2(P_i,q_low) + f·d2y_dx2(P_i,q_hi)
  P        = a·y_qlow + b·y_qhi + c·d2y_dx2(P_i,q_low) + d·d2y_dx2(P_i,q_hi)

  theta              = ρ · P
  dtheta_dn          = P + dP_dq0 · dq0_dq · dq_dn_n
  dtheta_dgradn      = dP_dq0 · dq0_dq · dq_dgradn_n_gmod
  d2theta_dn2        = f(dP_dq0, d2P_dq02, dq0_dq, dq_dn_n, d2q0_dq2, dn_dq_dn_n_n)
  dn_dtheta_dgradn   = f(d2P_dq02, dP_dq0, dq0_dq, dq_dn_n, dq_dgradn_n_gmod, d2q0_dq2)
  dgradn_dtheta_dgradn = f(d2P_dq02, dP_dq0, dq0_dq, dq_dgradn_n_gmod, d2q0_dq2)

  gmod              = |∇ρ|                       [P_i-independent]
  gradn_graddeltan  = ∇ρ · ∇δρ                   [P_i-independent]
```

---

## `get_u_delta_u` — Kernel Convolution in G-space

```
u(nnr, Nqs), delta_u(nnr, Nqs)  [real-space theta and delta-theta]

  fwfft(u[:,P])       for P = 1…Nqs
  fwfft(delta_u[:,P]) for P = 1…Nqs

  last_g = -1
  for G = 1…ngm:
      if igtongl(G) ≠ last_g:                      [new shell]
          gmod  = sqrt(gl(igtongl(G))) · tpiba
          interpolate_kernel(gmod, kernel_of_g)    [Nqs×Nqs spline]
          last_g = igtongl(G)
      gqmod = |G+q| · tpiba
      interpolate_kernel(gqmod, kernel_of_gq)      [always, no shell reuse]

      temp_u[nl(G), α]      += Σ_β kernel_of_g(α,β)  · u[nl(G), β]
      temp_delta_u[nl(G), α] += Σ_β kernel_of_gq(α,β) · delta_u[nl(G), β]

  invfft(temp_u[:,P])       for P = 1…Nqs
  invfft(temp_delta_u[:,P]) for P = 1…Nqs

  u ← temp_u,  delta_u ← temp_delta_u
```

**Note**: `kernel_of_g` is reused across G-vectors of the same shell magnitude.
`kernel_of_gq` is recomputed for every G-vector (|G+q| breaks shell symmetry).

---

## Key Dependencies Between Stages

```
Stage                       Reads                         Writes
─────────────────────────── ───────────────────────────── ──────────────────────────────
fft_gradient_r2r            total_rho                     gradient_rho
fft_qgradient               drho(:,1), q_point            gradient_drho
fill_q0_extended_on_grid    total_rho, gradient_rho       q0, dq*, d2q0_dq2
initialize_spline_interp    q_mesh                        d2y_dx2  (once, saved)
Loop 1                      q0, d2y_dx2, dq*, grad*       b1, b2, u, delta_u,
                            total_rho, drho                dtheta_dgradn_save, h1part2_save
get_u_delta_u               u, delta_u, q_point           u, delta_u  (in-place)
Accumulation                u, delta_u, b1, b2            delta_v
Loop 2                      u, delta_u,                   h1t, h2t
                            dtheta_dgradn_save,
                            h1part2_save
Divergence (icar=1,2,3)     h1t, h2t, gradient_rho,       delta_v
                            gradient_drho, g, q_point
```
