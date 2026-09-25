"""
    TimoshenkoBeam <: AbstractHydroelasticStructure

1-D Timoshenko beam for hydroelastic problems (two-field formulation).

Two displacement fields:
- `w`  — transverse deflection [m]
- `θ`  — cross-section rotation [rad]

The beam axis is defined by `tangent` (default `VectorValue(1.0, 0.0)` for
horizontal beams embedded in a 2-D fluid domain).  All computed quantities
are normalised by the fluid density `ρ_w`.

Use when the thickness-to-span ratio satisfies h_beam/L ≥ 1/20.  In the
thin limit (h_beam → 0 at fixed L) the model converges to the
Euler–Bernoulli beam; see `EulerBernoulliBeam`.

Mixed-order interpolation (`fe_w` at order p, `fe_θ` at order p-1) is
recommended to avoid shear locking.  The defaults (order 2 and 1) satisfy
this requirement.

Subtypes [`AbstractHydroelasticStructure`](@ref), but as a genuinely
two-field structure whose elastic operator couples `w` and `θ` together
(the shear term mixes `∂w/∂s` and `θ`), it cannot share the single-field
`mass` default: `mass` is implemented directly. It *does* share `stiffness`'s
hydrostatic term and all of `rhs`, since both only ever touch the
deflection field `w` (`variable_symbol(s) == s.symbol_w`); see
[`mass`](@ref)`(s::TimoshenkoBeam, ...)` and
[`stiffness_operator`](@ref)`(s::TimoshenkoBeam, ...)` for the parts that
are implemented directly. `f` passed to `rhs` must be a `FieldMap` with a
key equal to `s.symbol_w` (the `θ` component, if present, is ignored).
`has_damping_form(::TimoshenkoBeam) = false`, so `damping_parameter` is
not implemented.

# Fields
- `E::Float64`              — Young's modulus [Pa]
- `ν::Float64`              — Poisson's ratio [-]
- `h_beam::Float64`         — cross-section height (bending direction) [m]
- `b_beam::Float64`         — cross-section width [m]
- `ρ_s::Float64`            — structural density [kg/m³]
- `ρ_w::Float64`            — fluid density for normalisation [kg/m³];
                              default 1025.0
- `g::Float64`              — gravitational acceleration [m/s²];
                              default 9.81
- `κ::Float64`              — shear correction factor [-];
                              default 5/6 for rectangular cross-sections
- `tangent`                 — unit tangent along the beam axis; must match the
                              ambient dimension of the FE mesh (default
                              `VectorValue(1.0, 0.0)` for a horizontal beam in
                              a 2-D domain; use `VectorValue(1.0)` for 1-D
                              test meshes)
- `symbol_w::Symbol`        — deflection field symbol (default `:w`)
- `symbol_θ::Symbol`        — rotation field symbol (default `:θ`)
- `space_domain_symbol::Symbol` — key in `TankTriangulations` (default `:Γη`)
- `fe_w::FESpaceConfig`     — FE config for `w` (default `order = 2`)
- `fe_θ::FESpaceConfig`     — FE config for `θ` (default `order = 1`,
                              prevents shear locking)

# See also
[`EulerBernoulliBeam`](@ref), [`KirchhoffLovePlate`](@ref)

# References
- [A24] Agarwal, S., Colomes, O., & Metrikine, A. V. (2024).
  Dynamic analysis of viscoelastic floating membranes using monolithic
  finite element method. *Journal of Fluids and Structures*, 129, 104167.
  DOI: https://doi.org/10.1016/j.jfluidstructs.2024.104167
  Cited here only for related monolithic hydroelastic finite-element
  methodology; it is not a reference for the Timoshenko beam kinematics or
  constitutive model used by `TimoshenkoBeam`.
- [C23] Colomes, O., Verdugo, F., & Akkerman, I. (2023). A monolithic
  finite element formulation for the hydroelastic analysis of very large
  floating structures. *Int. J. Numer. Methods Eng.*, 124(3), 714-751.
  DOI: https://doi.org/10.1002/nme.7140
"""
@with_kw struct TimoshenkoBeam <: AbstractHydroelasticStructure
  E::Float64
  ν::Float64
  h_beam::Float64
  b_beam::Float64
  ρ_s::Float64
  ρ_w::Float64 = 1025.0
  g::Float64   = 9.81
  κ::Float64   = 5 / 6
  tangent::Any = VectorValue(1.0, 0.0)
  symbol_w::Symbol            = :w
  symbol_θ::Symbol            = :θ
  space_domain_symbol::Symbol = :Γη
  fe_w::FESpaceConfig = FESpaceConfig(order = 2)
  fe_θ::FESpaceConfig = FESpaceConfig(order = 1)
end

function print_parameters(s::TimoshenkoBeam)
  A_val = s.b_beam * s.h_beam
  I_val = s.b_beam * s.h_beam^3 / 12
  G     = s.E / (2 * (1 + s.ν))
  @printf("\n[MSG] TimoshenkoBeam Properties:\n")
  @printf("[VAL] E      = %.4e Pa\n",    s.E)
  @printf("[VAL] ν      = %.4f\n",       s.ν)
  @printf("[VAL] h_beam = %.4f m\n",     s.h_beam)
  @printf("[VAL] b_beam = %.4f m\n",     s.b_beam)
  @printf("[VAL] ρ_s    = %.4f kg/m³\n", s.ρ_s)
  @printf("[VAL] ρ_w    = %.4f kg/m³\n", s.ρ_w)
  @printf("[VAL] κ      = %.4f\n",        s.κ)
  @printf("[VAL] EI     = %.4e N·m²\n",  s.E * I_val)
  @printf("[VAL] κGA    = %.4e N\n",      s.κ * G * A_val)
  println()
end

variable_symbol(s::TimoshenkoBeam)  = s.symbol_w
variable_symbols(s::TimoshenkoBeam) = (s.symbol_w, s.symbol_θ)
field_fe_configs(s::TimoshenkoBeam) = (s.fe_w, s.fe_θ)

# Timoshenko beams do not have a standalone damping form (Rayleigh damping is not implemented).
has_damping_form(::TimoshenkoBeam) = false

# ── Two-field weak forms ────────────────────────────────────────────────────
#
# TimoshenkoBeam is a genuinely two-field structure (w, θ) with an
# elastic operator that *couples* both fields (the shear term mixes
# ∂w/∂s and θ), so it cannot share AbstractHydroelasticStructure's
# single-field `mass` default — `mass` is therefore implemented directly
# below, exactly as before this refactor.
#
# The hydrostatic restoring term `g·v_w·w`, however, only ever touches the
# deflection field `w` = variable_symbol(s), which is exactly what the
# shared `stiffness` default already assumes. So `stiffness_operator` below
# supplies only the coupled bending+shear part, and AbstractHydroelasticStructure's
# shared `stiffness` adds the `g·v_w·w` hydrostatic term around it — this
# reproduces the original combined bilinear form exactly (a single ∫(...)dΩ
# call vs. two summed ∫(...)dΩ calls over the same measure integrate to the
# same value). `rhs` (which also only ever touches `w`) is inherited
# unchanged from the shared default for the same reason.

"""
    mass(s::TimoshenkoBeam, dom, x_tt, y)

Timoshenko beam inertia form (normalised by ρ_w):

```math
m(\\ddot{u}, v) =
  \\frac{ρ_s A}{ρ_w}\\int_{Γ} v_w \\ddot{w}\\,dΓ
  + \\frac{ρ_s I}{ρ_w}\\int_{Γ} v_θ \\ddot{θ}\\,dΓ
```

where ``A = b \\cdot h``, ``I = b h^3/12``. Implemented directly (not via
[`AbstractHydroelasticStructure`](@ref)'s shared single-field `mass`
default) because the two fields `w` and `θ` carry different, independent
inertia coefficients (translational vs. rotational).
"""
function mass(s::TimoshenkoBeam, dom::IntegrationDomains, x_tt, y)
  w_tt = x_tt[s.symbol_w]
  θ_tt = x_tt[s.symbol_θ]
  v_w  = y[s.symbol_w]
  v_θ  = y[s.symbol_θ]

  A_val = s.b_beam * s.h_beam
  I_val = s.b_beam * s.h_beam^3 / 12
  m_A   = s.ρ_s * A_val / s.ρ_w   # translational mass per unit length / ρ_w
  m_I   = s.ρ_s * I_val / s.ρ_w   # rotational inertia per unit length / ρ_w
  dΩ    = _space_measure(dom, s)

  ∫(m_A * v_w * w_tt + m_I * v_θ * θ_tt)dΩ
end

"""
    stiffness_operator(s::TimoshenkoBeam, dom::IntegrationDomains, x, y)

Timoshenko beam elastic (bending + shear) stiffness operator, ``w``-``θ``
coupled, *excluding* the hydrostatic restoring term:

```math
k_{\\mathrm{op}}(u, v) =
  \\frac{EI}{ρ_w}\\int_{Γ} (∂_s θ)(∂_s v_θ)\\,dΓ
  + \\frac{κGA}{ρ_w}\\int_{Γ} (∂_s w - θ)(∂_s v_w - v_θ)\\,dΓ
```

where ``∂_s f = ∇(f) \\cdot t`` is the tangential derivative along the beam
axis defined by `s.tangent`, ``EI = E b h^3/12``, ``G = E/[2(1+ν)]``, and
``κGA = κ G b h``.

Combined with the shared hydrostatic term (in
[`stiffness`](@ref)`(s::AbstractHydroelasticStructure, ...)`, which uses
[`variable_symbol`](@ref)`(s) == s.symbol_w` to apply `g·v_w·w` — the
deflection field only) this reproduces the original beam bilinear form
exactly.
"""
function stiffness_operator(s::TimoshenkoBeam, dom::IntegrationDomains, x, y)
  w   = x[s.symbol_w]
  θ   = x[s.symbol_θ]
  v_w = y[s.symbol_w]
  v_θ = y[s.symbol_θ]

  A_val  = s.b_beam * s.h_beam
  I_val  = s.b_beam * s.h_beam^3 / 12
  G      = s.E / (2 * (1 + s.ν))
  EI_ρ   = s.E * I_val / s.ρ_w
  κGA_ρ  = s.κ * G * A_val / s.ρ_w
  t      = s.tangent
  dΩ     = _space_measure(dom, s)

  # Timoshenko two-field structural bilinear form (w, θ).
  # Includes bending and shear terms; thin limit recovers EB behavior.
  # Standard Timoshenko beam formulation; see [C23] for the
  # hydroelastic monolithic coupling context used in this model.
  ∫(
    EI_ρ  * (∇(θ)   ⋅ t) * (∇(v_θ) ⋅ t) +
    κGA_ρ * ((∇(w)  ⋅ t) - θ) * ((∇(v_w) ⋅ t) - v_θ)
  )dΩ
end
