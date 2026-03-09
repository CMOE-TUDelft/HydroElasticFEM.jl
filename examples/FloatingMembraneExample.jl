module FloatingMembraneExample
using Gridap
using HydroElasticFEM: PKG_ROOT
import HydroElasticFEM.Geometry as G

# ── Domain setup ────────────────────────────────────────────
# A 10 × 1 m tank with a 1 m floating membrane centred at x = 4.5–5.5 m
# on the free surface (y = 1.0), plus 0.5 m damping zones at each end
# to absorb outgoing waves.
structure_domain = G.StructureDomain1D(L=1.0, x₀=[4.5, 1.0])
damping_zone_left = G.DampingZone1D(L=0.5, x₀=[0.0, 1.0])
damping_zone_right = G.DampingZone1D(L=0.5, x₀=[9.5, 1.0])
tank_domain = G.TankDomain2D(L=10.0, H=1.0, nx=20, ny=2, 
  structure_domains=[structure_domain], 
  damping_zones=[damping_zone_left, damping_zone_right])

# ── Model & triangulations ──────────────────────────────────
# Build the Cartesian FE model, label boundaries, and partition
# the free surface into structure / damping / open-water zones.
model = G.build_model(tank_domain)
tank_trians = G.build_triangulations(tank_domain, model)

# ── VTK output ──────────────────────────────────────────────
# Write each triangulation to VTK for visualisation in ParaView.
filedir = joinpath(PKG_ROOT, "data", "VTK", "examples", "FloatingMembraneExample")

writevtk(model, joinpath(filedir, "floating_membrane_model"))            # Full model
writevtk(tank_trians.Ω, joinpath(filedir, "floating_membrane_interior")) # Fluid interior
writevtk(tank_trians.Γ, joinpath(filedir, "floating_membrane_surface"))  # Entire top surface
writevtk(tank_trians.Γ_structures[1], joinpath(filedir, "floating_membrane_membrane"))  # Membrane region
writevtk(tank_trians.Γ_dampings[1], joinpath(filedir, "floating_membrane_damping1"))    # Left damping zone
writevtk(tank_trians.Γ_dampings[2], joinpath(filedir, "floating_membrane_damping2"))    # Right damping zone
writevtk(tank_trians.Γin, joinpath(filedir, "floating_membrane_inlet"))   # Inlet wall
writevtk(tank_trians.Γout, joinpath(filedir, "floating_membrane_outlet")) # Outlet wall
writevtk(tank_trians.Γbot, joinpath(filedir, "floating_membrane_bottom")) # Bottom wall
writevtk(tank_trians.Γfs, joinpath(filedir, "floating_membrane_free_surface"))  # Free surface (no structure/damping)
writevtk(tank_trians.Γη, joinpath(filedir, "floating_membrane_eta"))      # All-structure surface (η)
writevtk(tank_trians.Γκ, joinpath(filedir, "floating_membrane_kappa"))    # Non-structure surface (κ)


end