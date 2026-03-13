module FloatingMembraneExample
using Gridap
using HydroElasticFEM: PKG_ROOT
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Simulation as S

# ── Domain setup ────────────────────────────────────────────
# A 10 × 1 m tank with a 1 m floating membrane centred at x = 4.5–5.5 m
# on the free surface (y = 1.0), plus 0.5 m damping zones at each end
# to absorb outgoing waves.
structure_domain = G.StructureDomain1D(L=1.0, x₀=[4.5, 1.0], domain_symbol=:Γs)
damping_zone_left = G.DampingZone1D(L=0.5, x₀=[0.0, 1.0], domain_symbol=:Γd_1)
damping_zone_right = G.DampingZone1D(L=0.5, x₀=[9.5, 1.0], domain_symbol=:Γd_2)
tank_domain = G.TankDomain2D(L=10.0, H=1.0, nx=20, ny=2, 
  structure_domains=[structure_domain], 
  damping_zones=[damping_zone_left, damping_zone_right])

# Physics properties
membrane = P.Membrane2D(L=1.0, mᵨ=0.9, Tᵨ=98.1, space_domain_symbol=:Γs)
freesurf = P.FreeSurface(ρw=1025.0, g=9.81, βₕ=0.5, space_domain_symbol=:Γκ)
potential = P.PotentialFlow(ρw=1025.0, g=9.81, space_domain_symbol=:Ω)

# Simulation configuration
sim_configuration = S.FreqDomainConfig(ω=2.0)

# build the problem (discrete model, triangulations, integration domains, etc.)
problem = S.build_problem(tank_domain, [membrane, freesurf, potential], sim_configuration)

# ── VTK output ──────────────────────────────────────────────
# Write each triangulation to VTK for visualisation in ParaView.
filedir = joinpath(PKG_ROOT, "data", "VTK", "examples", "FloatingMembraneExample")
model = S.get_model(problem)
tank_trians = S.get_triangulations(problem)

writevtk(model, joinpath(filedir, "floating_membrane_model"))            # Full model
writevtk(tank_trians[:Ω], joinpath(filedir, "floating_membrane_interior")) # Fluid interior
writevtk(tank_trians[:Γ], joinpath(filedir, "floating_membrane_surface"))  # Entire top surface
writevtk(tank_trians[:Γη], joinpath(filedir, "floating_membrane_membrane"))  # Membrane region
writevtk(tank_trians[:Γd_1], joinpath(filedir, "floating_membrane_damping1"))    # Left damping zone
writevtk(tank_trians[:Γd_2], joinpath(filedir, "floating_membrane_damping2"))    # Right damping zone
writevtk(tank_trians[:Γin], joinpath(filedir, "floating_membrane_inlet"))   # Inlet wall
writevtk(tank_trians[:Γout], joinpath(filedir, "floating_membrane_outlet")) # Outlet wall
writevtk(tank_trians[:Γbot], joinpath(filedir, "floating_membrane_bottom")) # Bottom wall
writevtk(tank_trians[:Γfs], joinpath(filedir, "floating_membrane_free_surface"))  # Free surface (no structure/damping)
writevtk(tank_trians[:Γη], joinpath(filedir, "floating_membrane_eta"))      # All-structure surface (η)
writevtk(tank_trians[:Γκ], joinpath(filedir, "floating_membrane_kappa"))    # Non-structure surface (κ)

sim_result = S.simulate(problem)

end