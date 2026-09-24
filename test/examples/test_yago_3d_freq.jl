using Test
using Gridap

using HydroElasticFEM
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P

include(joinpath(@__DIR__, "..", "..", "examples", "YagoBenchmark3DFreq.jl"))
using .YagoBenchmark3DFreq

@testset "Yago 3D freq - graded mesh" begin
  domain = G.CartesianDomain(
    mins = (0.0, -540.0, 0.0),
    maxs = (3000.0, 540.0, 58.5),
    parts = (20, 36, 1),
    map = x -> G.map_fn(x, 58.5, 1; grading_base = 2.5),
  )

  Ω = G.triangulation(domain)
  zc = [sum(node[3] for node in cell) / length(cell)
        for cell in get_cell_coordinates(Ω)]
  @test maximum(zc) <= 58.5 + 1e-8
  @test minimum(zc) > 0.0
end

@testset "Yago 3D freq - plate mask" begin
  domain = G.CartesianDomain(
    mins = (0.0, -540.0, 0.0),
    maxs = (3000.0, 540.0, 58.5),
    parts = (20, 36, 1),
    map = x -> G.map_fn(x, 58.5, 1; grading_base = 2.5),
  )

  Γ = G.get_boundary(domain, "surface")
  L = 300.0
  B = 60.0
  Γb, Γf, _ = G.get_plate_triangulation(Γ, 4.5 * L, 5.5 * L, -B / 2, B / 2)

  xb = get_cell_coordinates(Γb)
  for cell in xb
    for node in cell
      @test 4.5 * L - 1e-10 <= node[1] <= 5.5 * L + 1e-10
      @test -B / 2 - 1e-10 <= node[2] <= B / 2 + 1e-10
    end
  end

  @test num_cells(Γb) + num_cells(Γf) == num_cells(Γ)
end

@testset "Yago 3D freq - C tensor symmetry" begin
  C = P.build_KL_tensor(11.9e9, 0.13, 2.0, 1025.0)
  for i in 1:3, j in 1:3, k in 1:3, l in 1:3
    @test C[i, j, k, l] ≈ C[k, l, i, j] atol=1e-10
  end
end

# Shared by the warm-up and the legacy-regression testsets (one solve only).
xs_warmup, η_warmup = nothing, nothing

@testset "Yago 3D freq - warm-up solve" begin
  global xs_warmup, η_warmup
  xs, η_rel = run_yago_3d_freq(
    nx = 2,
    ny = 2,
    nz = 1,
    order = 2,
    λfactor = 0.4,
    dfactor = 2.0,
    vtk_output = false,
    verbose = false,
  )

  xs_warmup, η_warmup = xs, η_rel
  @test length(xs) > 0
  @test all(0.0 .<= η_rel .<= 5.0)
  @test maximum(η_rel) > 1e-6
end

# -------------------------------------------------------------------------
# Regression: TankDomain geometry vs. legacy get_plate_triangulation path
#
# The example now builds all domains from a 3D TankDomain.  The legacy
# path below (plain CartesianDomain + get_plate_triangulation + hand-made
# measures) is kept only here, to check that the new geometry layer
# reproduces it exactly when fed to the same solver.
# -------------------------------------------------------------------------

function _yago_legacy_geo(c, nx, ny, nz, order)
  domain = G.CartesianDomain(
    mins = (0.0, -c.BΩ / 2, 0.0),
    maxs = (c.LΩ, c.BΩ / 2, c.H),
    parts = (c.nLΩ * nx, c.nBΩ * ny, nz),
    map = x -> G.map_fn(x, c.H, nz; grading_base = 2.5),
  )
  model = G.build_model(domain)
  trians = G.build_triangulations(domain, model)
  Ω, Γ, Γᵢₙ = trians[:Ω], trians[:Γfs], trians[:Γin]
  Γb, Γf, Λb = G.get_plate_triangulation(Γ, c.xb₀, c.xb₁, c.yb₀, c.yb₁)
  (Ω = Ω, Γf = Γf, Γb = Γb,
   dΩ = Measure(Ω, 2 * order), dΓᵢₙ = Measure(Γᵢₙ, 2 * order),
   dΓf = Measure(Γf, 2 * order), dΓb = Measure(Γb, 2 * order),
   dΛb = Measure(Λb, 2 * order), nΛb = get_normal_vector(Λb), Λb = Λb)
end

_sorted_centroids(trian) =
  sort([Tuple(round.(Tuple(sum(c) / length(c)); digits = 8)) for c in get_cell_coordinates(trian)])

@testset "Yago 3D freq - TankDomain geometry matches legacy plate mask" begin
  nx, ny, nz, order, dfactor = 2, 2, 1, 2, 2.0
  c = YagoBenchmark3DFreq._constants()
  legacy = _yago_legacy_geo(c, nx, ny, nz, order)

  tank = yago_tank(c, nx, ny, nz, dfactor)
  trians = G.build_triangulations(tank, G.build_model(tank))

  @test _sorted_centroids(trians[:Γη]) == _sorted_centroids(legacy.Γb)
  @test _sorted_centroids(trians[:Γκ]) == _sorted_centroids(legacy.Γf)
  @test num_cells(trians[:Λη]) == num_cells(legacy.Λb)
  # damping zones: 4 cells long (Ld = 600 m, Δx = 150 m) × full width
  @test num_cells(trians[:Γ_d_in]) == 4 * c.nBΩ * ny
  @test num_cells(trians[:Γ_d_out]) == 4 * c.nBΩ * ny
end

@testset "Yago 3D freq - TankDomain solve matches legacy solve" begin
  nx, ny, nz, order, λfactor, dfactor = 2, 2, 1, 2, 0.4, 2.0
  c = YagoBenchmark3DFreq._constants()
  wave = YagoBenchmark3DFreq._wave(c, λfactor)
  damping = YagoBenchmark3DFreq._damping(c, wave, dfactor)

  legacy = _yago_legacy_geo(c, nx, ny, nz, order)
  ξ_old, η_old = solve_yago(c, wave, damping, legacy;
    nx, order, λfactor, vtk_output = false, verbose = false)
  # Same parameters as the warm-up solve above (TankDomain geometry)
  @test xs_warmup == ξ_old
  @test η_warmup ≈ η_old rtol = 1e-8
end
