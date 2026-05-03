using Test
using Gridap
import HydroElasticFEM.Geometry as G

# ─────────────────────────────────────────────────────────────
# Helper: locate the test fixture mesh
# ─────────────────────────────────────────────────────────────

const _FIXTURES_DIR = joinpath(
  @__DIR__, "..", "fixtures", "meshes",
)
const _TANK_MSH         = joinpath(_FIXTURES_DIR, "tank2d.msh")
const _MISSING_TAGS_MSH = joinpath(_FIXTURES_DIR, "tank2d_missing_tags.msh")

# ─────────────────────────────────────────────────────────────
# Platform guard
# ─────────────────────────────────────────────────────────────
#
# Julia 1.11 on macOS ARM64 (Apple Silicon) exhibits an LLVM JIT crash
# (SIGILL in buildEarlyOptimizerPipeline / _Znwm) when GridapGmsh methods
# are first compiled after ~14–46 M allocations in the same process.
# This is a confirmed upstream bug in the Julia 1.11 + GridapGmsh
# interaction on that platform.  Skip GmshDomain tests there to keep
# the test suite green, and verify separately on Linux/x86-64.
#
# Remove this guard once the upstream bug is fixed.
const _GMSH_SKIP = Sys.isapple() && Sys.ARCH == :aarch64

if _GMSH_SKIP
  @warn "GmshDomain tests skipped: Julia 1.11 JIT crash on macOS ARM64 " *
        "(upstream GridapGmsh incompatibility).  Run on Linux/x86-64 to verify."

  @testset "GmshDomain (skipped — macOS ARM64 JIT bug)" begin
    @test_skip "Skipped on macOS ARM64: Julia 1.11 + GridapGmsh SIGILL"
  end

else  # ── all other platforms ──────────────────────────────────────────────

  const _DOM = G.GmshDomain(_TANK_MSH)

  # ─────────────────────────────────────────────────────────────
  # GmshDomain — constructor
  # ─────────────────────────────────────────────────────────────

  @testset "GmshDomain — constructor loads mesh" begin
    @test _DOM isa G.GmshDomain
    @test _DOM isa G.AbstractDomain
  end

  @testset "GmshDomain — file not found raises error" begin
    @test_throws Exception G.GmshDomain("/nonexistent/path/mesh.msh")
  end

  @testset "GmshDomain — missing required tags raises error" begin
    @test_throws Exception G.GmshDomain(_MISSING_TAGS_MSH)
  end

  # ─────────────────────────────────────────────────────────────
  # GmshDomain — AbstractDomain interface
  # ─────────────────────────────────────────────────────────────

  @testset "GmshDomain — ambient_dimension / manifold_dimension" begin
    @test G.ambient_dimension(_DOM) == 2
    @test G.manifold_dimension(_DOM) == 2
  end

  @testset "GmshDomain — boundary_tags contains all standard tags" begin
    btags = G.boundary_tags(_DOM)
    @test btags isa Dict{String, String}
    for tag in G.STANDARD_TAGS
      @test haskey(btags, tag)
    end
  end

  @testset "GmshDomain — triangulation returns interior" begin
    Ω = G.triangulation(_DOM)
    @test num_cells(Ω) > 0
  end

  @testset "GmshDomain — get_boundary returns correct triangulations" begin
    for name in G.STANDARD_TAGS
      name == "fluid" && continue  # fluid is interior; not a boundary tag in Gmsh
      bnd = G.get_boundary(_DOM, name)
      @test num_cells(bnd) > 0
    end
  end

  @testset "GmshDomain — get_boundary invalid tag raises error" begin
    @test_throws Exception G.get_boundary(_DOM, "nonexistent_tag_xyz")
  end

  # ─────────────────────────────────────────────────────────────
  # validate_gmsh_tags
  # ─────────────────────────────────────────────────────────────

  @testset "validate_gmsh_tags — all tags present returns true" begin
    @test G.validate_gmsh_tags(_TANK_MSH) == true
  end

  @testset "validate_gmsh_tags — missing tags returns false" begin
    @test G.validate_gmsh_tags(_MISSING_TAGS_MSH) == false
  end

  @testset "validate_gmsh_tags — file not found returns false" begin
    @test G.validate_gmsh_tags("/nonexistent/path/mesh.msh") == false
  end

  @testset "validate_gmsh_tags — custom required list" begin
    result = G.validate_gmsh_tags(
      _MISSING_TAGS_MSH;
      required=["fluid", "seabed", "inlet", "outlet"],
    )
    @test result == true
  end

  # ─────────────────────────────────────────────────────────────
  # GmshDomain — Simulation pipeline bridge
  # ─────────────────────────────────────────────────────────────

  @testset "GmshDomain — build_model returns discrete model" begin
    model = G.build_model(_DOM)
    @test model === _DOM.model
  end

  @testset "GmshDomain — build_triangulations returns TankTriangulations" begin
    model  = G.build_model(_DOM)
    trians = G.build_triangulations(_DOM, model)
    @test trians isa G.TankTriangulations
    for key in [:Ω, :Γfs, :Γη, :Γbot, :Γin, :Γout,
                :Γ_structures, :Γ_dampings, :Λη, :Λ_joints, :joint_domains]
      @test haskey(trians, key)
    end
    @test num_cells(trians[:Ω]) > 0
    @test num_cells(trians[:Γfs]) > 0
    @test num_cells(trians[:Γη]) > 0
    @test length(trians[:Γ_structures]) == 1
    @test length(trians[:Γ_dampings])   == 0
  end

end  # _GMSH_SKIP
