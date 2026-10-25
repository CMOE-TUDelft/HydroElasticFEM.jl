# Test helper — semi-analytical reference: 2D rectangular box (breadth 2b, draft d) floating
# in water of depth h, radiation by the matched eigenfunction expansion
# (Zheng, You & Shen 2004, Ocean Eng. 31; Black, Mei & Bray 1971).
# Time dependence e^{−iωt}; x from the box centre, z up, free surface z = 0.
# Modes (surge, heave, pitch) with the HydroElasticFEM 2D convention
# u = (ξ₁ + α r_z, ξ₃ − α r_x), r = x − x_r, reference point (0, z_r).
using LinearAlgebra

# Gauss–Legendre nodes and weights on [a, c] (Golub–Welsch).
function _box_gl(n, a, c)
  β = [k / sqrt(4k^2 - 1) for k in 1:n-1]
  E = eigen(SymTridiagonal(zeros(n), β))
  x, w = E.values, 2 .* E.vectors[1, :] .^ 2
  (c - a) / 2 .* x .+ (a + c) / 2, (c - a) / 2 .* w
end

# Exterior wavenumbers: k₀ tanh k₀h = ν, and kₙ tan kₙh = −ν (n ≥ 1).
function _box_wavenumbers(ν, h, N)
  f0(k) = k * tanh(k * h) - ν
  k0 = max(ν, sqrt(ν / h))
  for _ in 1:100
    k0 -= f0(k0) / (tanh(k0 * h) + k0 * h * sech(k0 * h)^2)
  end
  ks = Float64[]
  for n in 1:N
    lo, hi = (n - 0.5) * π / h + 1e-12, n * π / h - 1e-12
    f(k) = k * tan(k * h) + ν
    for _ in 1:200
      mid = (lo + hi) / 2
      f(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    push!(ks, (lo + hi) / 2)
  end
  k0, ks
end

"""
    box_radiation_reference(; b, d, h, ω, z_r = 0.0, ρ = 1025.0, g = 9.81, M = 40)
        -> (A, B)

Added mass and radiation damping (per unit width, dimensional; modes surge,
heave, pitch about (0, z_r)) of a 2D floating box by eigenfunction matching
with `M` gap modes.  Same definition as `radiation_coefficients`:
F_rad = (ω²A + iωB) ξ.
"""
function box_radiation_reference(; b, d, h, ω, z_r = 0.0, ρ = 1025.0, g = 9.81, M = 40)
  s = h - d
  ν = ω^2 / g
  N = ceil(Int, M * h / s)
  k0, ks = _box_wavenumbers(ν, h, N)
  # exterior modes Zₙ(z) and decay rates κₙ: φ_I = Σ Aₙ Zₙ e^{−κₙ(x−b)}
  Z = [z -> cosh(k0 * (z + h)) / cosh(k0 * h); [let k = k; z -> cos(k * (z + h)) / cos(k * h); end for k in ks]]
  κ = ComplexF64[-im * k0; ks]
  λ = [m * π / s for m in 0:M]
  cm(m, z) = cos(λ[m+1] * (z + h))
  nq = 600
  zg, wg = _box_gl(nq, -h, -d)        # gap
  zw, ww = _box_gl(nq, -d, 0.0)       # wall
  zf, wf = _box_gl(2nq, -h, 0.0)      # full depth
  xb, wb = _box_gl(nq, 0.0, b)        # body bottom (half)
  nN, nM = N + 1, M + 1
  Ezc = [sum(wg .* Z[n].(zg) .* cm.(m, zg)) for m in 0:M, n in 1:nN]  # ∫gap Zₙ cos_m
  Ezz = [sum(wf .* Z[n].(zf) .^ 2) for n in 1:nN]                     # ∫full Zₙ²
  Ecc = [m == 0 ? s : s / 2 for m in 0:M]
  F = zeros(ComplexF64, 3, 3)
  for mode in 1:3
    sym = mode == 2
    # interior φ_II = φ_p + Σ C_m Y_m(x) cos_m(z)
    if sym
      φp = (x, z) -> ((z + h)^2 - x^2) / (2s)
      dφp = z -> -b / s
      Y = (m, x) -> m == 0 ? 1.0 : cosh(λ[m+1] * x) / cosh(λ[m+1] * b)
      dY = m -> m == 0 ? 0.0 : λ[m+1] * tanh(λ[m+1] * b)
      uwall = z -> 0.0
    else
      if mode == 1
        φp = (x, z) -> 0.0; dφp = z -> 0.0; uwall = z -> 1.0
      else
        φp = (x, z) -> -x * (z + h)^2 / (2s) + x^3 / (6s)
        dφp = z -> -(z + h)^2 / (2s) + b^2 / (2s)
        uwall = z -> z - z_r
      end
      Y = (m, x) -> m == 0 ? x / b : sinh(λ[m+1] * x) / sinh(λ[m+1] * b)
      dY = m -> m == 0 ? 1 / b : λ[m+1] * coth(λ[m+1] * b)
    end
    # unknowns [A₀…A_N; C₀…C_M]
    K = zeros(ComplexF64, nN + nM, nN + nM)
    r = zeros(ComplexF64, nN + nM)
    for m in 0:M   # potential continuity on the gap, projected on cos_m
      K[m+1, 1:nN] = Ezc[m+1, :]
      K[m+1, nN+m+1] = -Ecc[m+1]
      r[m+1] = sum(wg .* φp.(b, zg) .* cm.(m, zg))
    end
    for n in 1:nN  # normal velocity over the full depth, projected on Zₙ
      row = nM + n
      K[row, n] = -κ[n] * Ezz[n]
      for m in 0:M
        K[row, nN+m+1] = -dY(m) * Ezc[m+1, n]
      end
      r[row] = sum(wg .* dφp.(zg) .* Z[n].(zg)) + sum(ww .* uwall.(zw) .* Z[n].(zw))
    end
    sol = K \ r
    An, Cm = sol[1:nN], sol[nN+1:end]
    φI(z) = sum(An[n] * Z[n](z) for n in 1:nN)
    φII(x, z) = φp(x, z) + sum(Cm[m+1] * Y(m, x) * cm(m, z) for m in 0:M)
    wall = sum(ww .* φI.(zw))                       # ∫ φ(b, z) dz on the wall
    wallr = sum(ww .* φI.(zw) .* (zw .- z_r))
    bot = sum(wb .* φII.(xb, -d))                   # ∫₀ᵇ φ(x, −d) dx
    botx = sum(wb .* xb .* φII.(xb, -d))
    # F_i = ∫ φ N_i over both halves (N from the fluid normal, into the body)
    F[1, mode] = sym ? 0.0 : -2wall
    F[2, mode] = sym ? 2bot : 0.0
    F[3, mode] = sym ? 0.0 : -2botx - 2wallr
  end
  # unit-velocity potentials: F_rad = ω²ρ ∫ φ N  ⇒  A = ρ Re, B = ρω Im
  return ρ .* real(F), ρ * ω .* imag(F)
end
