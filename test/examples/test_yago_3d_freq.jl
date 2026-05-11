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

@testset "Yago 3D freq - warm-up solve" begin
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

  @test length(xs) > 0
  @test all(0.0 .<= η_rel .<= 5.0)
  @test maximum(η_rel) > 1e-6
end
