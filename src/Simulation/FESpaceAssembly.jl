"""
    module FESpaceAssembly

Automated FE space construction for HydroElasticFEM.

Provides `build_fe_spaces` which takes entity-triangulation pairs and
produces multi-field test/trial spaces plus the field-index mapping,
eliminating manual `TestFESpace`/`TrialFESpace`/`MultiFieldFESpace`
boilerplate.
"""
module FESpaceAssembly

using Gridap
using Gridap.FESpaces: ConstantFESpace, SingleFieldFESpace
using Gridap.ODEs: TransientTrialFESpace, TransientMultiFieldFESpace

import ...Physics as P
import ...Geometry as G
import ...ParameterHandler as PH

# ─────────────────────────────────────────────────────────────
# Per-entity FE space construction (dispatched on entity type)
# ─────────────────────────────────────────────────────────────

"""
    build_test_fe_space(entity::PhysicsParameters, trian, config::SimulationConfig)

Build a Gridap `TestFESpace` for `entity` on triangulation `trian`,
using the numerical parameters from `entity.fe`.
"""
function build_test_fe_space(entity::P.PhysicsParameters, trian, config::PH.SimulationConfig)
    fe = entity.fe
    reffe = ReferenceFE(fe.reffe_type, fe.space_type, fe.order)
    if config isa PH.FreqDomainConfig
        @assert fe.vector_type == Vector{ComplexF64} "Frequency-domain simulations require complex-valued FE spaces."    
    end
    if fe.dirichlet_tags !== nothing
        TestFESpace(trian, reffe;
            conformity   = fe.conformity,
            vector_type  = fe.vector_type,
            dirichlet_tags = fe.dirichlet_tags)
    else
        TestFESpace(trian, reffe;
            conformity  = fe.conformity,
            vector_type = fe.vector_type)
    end
end

"""
    build_test_fe_space(resn::Vector{P.ResonatorSingle}, trian, config::PH.SimulationConfig)

Build one `ConstantFESpace` per resonator on triangulation `trian`.
Returns a `Vector` of test spaces.
"""
function build_test_fe_space(resn::Vector{P.ResonatorSingle}, trian, config::PH.SimulationConfig)
     if config isa PH.FreqDomainConfig
        @assert all(r -> r.fe.vector_type == Vector{ComplexF64}, resn) "Frequency-domain simulations require complex-valued FE spaces."    
    end
    [ConstantFESpace(trian;
        vector_type = iresn.fe.vector_type,
        field_type  = iresn.fe.space_type) 
     for iresn in resn]
end

"""
    build_trial_fe_space(entity::P.PhysicsParameters, V_test, config::PH.SimulationConfig)

Build a `TrialFESpace` (or `TransientTrialFESpace` when `config` is `TimeDomainConfig`)
from the test space `V_test`, applying Dirichlet values from `entity.fe`
if present.
"""
function build_trial_fe_space(entity::P.PhysicsParameters, V_test, config::PH.SimulationConfig)
    fe = entity.fe
    if config isa PH.TimeDomainConfig
        if fe.dirichlet_value !== nothing
            TransientTrialFESpace(V_test, fe.dirichlet_value)
        else
            TransientTrialFESpace(V_test)
        end
    else
        if fe.dirichlet_value !== nothing
            TrialFESpace(V_test, fe.dirichlet_value)
        else
            TrialFESpace(V_test)
        end
    end
end

"""
    build_trial_fe_space(resn::Vector{P.ResonatorSingle}, Vs::Vector, config::PH.SimulationConfig)

Build trial spaces for each resonator's test space.
"""
function build_trial_fe_space(resn::Vector{P.ResonatorSingle}, Vs::Vector, config::PH.SimulationConfig)
    if config isa PH.TimeDomainConfig
        [TransientTrialFESpace(V) for V in Vs]
    else
        [TrialFESpace(V) for V in Vs]
    end
end

# ─────────────────────────────────────────────────────────────
# Main manager function
# ─────────────────────────────────────────────────────────────


"""
    build_fe_spaces(entities::Vector{P.PhysicsParameters}, trians::G.TankTriangulations, config::PH.SimulationConfig) -> (X, Y, fmap)

Build multi-field test and trial FE spaces from a list of entities and a TankTriangulations dictionary.
Each entity uses the triangulation from trians[entity.domain_symbol].

# Arguments
- `entities`: Vector of `P.PhysicsParameters` (or Vector{P.ResonatorSingle})
- `trians`: `G.TankTriangulations` (dictionary-like)
- `config`: `PH.SimulationConfig` (frequency or time-domain)

# Returns
- `X::MultiFieldFESpace` — trial space
- `Y::MultiFieldFESpace` — test space
- `fmap::Dict{Symbol,Int}` — field-symbol to positional-index map

# Example
```julia
X, Y, fmap = build_fe_spaces([fluid, fsurf, mem], trians)
# fmap == Dict(:ϕ => 1, :κ => 2, :η_m => 3)
```
"""
function build_fe_spaces(entities::Vector{<:P.PhysicsParameters}, 
                         trians::G.TankTriangulations, 
                         config::PH.SimulationConfig)
    test_spaces  = SingleFieldFESpace[]
    trial_spaces = SingleFieldFESpace[]
    fmap = Dict{Symbol, Int}()
    idx  = 0

    for entity in entities
        trian = trians[entity.space_domain_symbol]
        if entity isa Vector{P.ResonatorSingle}
            Vs = build_test_fe_space(entity, trian, config)
            Us = build_trial_fe_space(entity, Vs, config)
            for (i, (Vi, Ui)) in enumerate(zip(Vs, Us))
                idx += 1
                push!(test_spaces, Vi)
                push!(trial_spaces, Ui)
                fmap[Symbol("q_$i")] = idx
            end
        else
            V = build_test_fe_space(entity, trian, config)
            U = build_trial_fe_space(entity, V, config)
            idx += 1
            push!(test_spaces, V)
            push!(trial_spaces, U)
            fmap[P.variable_symbol(entity)] = idx
        end
    end

    Y = MultiFieldFESpace(test_spaces)
    X = MultiFieldFESpace(trial_spaces)
    return X, Y, fmap
end

# ─────────────────────────────────────────────────────────────
# Exports
# ─────────────────────────────────────────────────────────────

export build_fe_spaces, build_test_fe_space, build_trial_fe_space

end # module FESpaceAssembly
