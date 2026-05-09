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
# Private helpers: build from a FESpaceConfig directly
# ─────────────────────────────────────────────────────────────

function _build_fe_test_space(fe::PH.FESpaceConfig, trian, ::PH.FreqDomainConfig)
    reffe = ReferenceFE(fe.reffe_type, fe.space_type, fe.order)
    @assert fe.vector_type == Vector{ComplexF64} "Frequency-domain simulations require complex-valued FE spaces."
    if fe.dirichlet_tags !== nothing
        TestFESpace(trian, reffe;
            conformity     = fe.conformity,
            vector_type    = fe.vector_type,
            dirichlet_tags = fe.dirichlet_tags)
    else
        TestFESpace(trian, reffe;
            conformity  = fe.conformity,
            vector_type = fe.vector_type)
    end
end

function _build_fe_test_space(fe::PH.FESpaceConfig, trian, ::PH.TimeDomainConfig)
    reffe = ReferenceFE(fe.reffe_type, fe.space_type, fe.order)
    if fe.dirichlet_tags !== nothing
        TestFESpace(trian, reffe;
            conformity     = fe.conformity,
            vector_type    = fe.vector_type,
            dirichlet_tags = fe.dirichlet_tags)
    else
        TestFESpace(trian, reffe;
            conformity  = fe.conformity,
            vector_type = fe.vector_type)
    end
end

function _build_fe_trial_space(fe::PH.FESpaceConfig, V_test, ::PH.FreqDomainConfig)
    fe.dirichlet_value !== nothing ? TrialFESpace(V_test, fe.dirichlet_value) :
                                     TrialFESpace(V_test)
end

function _build_fe_trial_space(fe::PH.FESpaceConfig, V_test, ::PH.TimeDomainConfig)
    fe.dirichlet_value !== nothing ? TransientTrialFESpace(V_test, fe.dirichlet_value) :
                                     TransientTrialFESpace(V_test)
end

# ─────────────────────────────────────────────────────────────
# Per-entity FE space construction (dispatched on entity type)
# ─────────────────────────────────────────────────────────────

"""
    build_test_fe_space(entity::PhysicsParameters, trian, config::SimulationConfig)

Build a Gridap `TestFESpace` for a **single-field** entity on triangulation
`trian`, using the numerical parameters from `entity.fe`.

For multi-field entities (e.g. `TimoshenkoBeam`) use the higher-level
`build_fe_spaces`, which iterates over `variable_symbols`/`field_fe_configs`.
"""
function build_test_fe_space(entity::P.PhysicsParameters, trian, config::PH.FreqDomainConfig)
    _build_fe_test_space(entity.fe, trian, config)
end

function build_test_fe_space(entity::P.PhysicsParameters, trian, config::PH.TimeDomainConfig)
    _build_fe_test_space(entity.fe, trian, config)
end

"""
    build_test_fe_space(resn::Vector{P.ResonatorSingle}, trian, config::PH.SimulationConfig)

Build one `ConstantFESpace` per resonator on triangulation `trian`.
Returns a `Vector` of test spaces.
"""
function build_test_fe_space(resn::Vector{P.ResonatorSingle}, trian, ::PH.FreqDomainConfig)
    @assert all(r -> r.fe.vector_type == Vector{ComplexF64}, resn) "Frequency-domain simulations require complex-valued FE spaces."
    [ConstantFESpace(trian;
        vector_type = iresn.fe.vector_type,
        field_type  = iresn.fe.space_type)
     for iresn in resn]
end

function build_test_fe_space(resn::Vector{P.ResonatorSingle}, trian, ::PH.TimeDomainConfig)
    [ConstantFESpace(trian;
        vector_type = iresn.fe.vector_type,
        field_type  = iresn.fe.space_type) 
     for iresn in resn]
end

"""
    build_trial_fe_space(entity::P.PhysicsParameters, V_test, config::PH.SimulationConfig)

Build a `TrialFESpace` (or `TransientTrialFESpace` when `config` is
`TimeDomainConfig`) from the test space `V_test`, applying Dirichlet values
from `entity.fe` if present.

For multi-field entities use `build_fe_spaces`.
"""
function build_trial_fe_space(entity::P.PhysicsParameters, V_test, config::PH.FreqDomainConfig)
    _build_fe_trial_space(entity.fe, V_test, config)
end

function build_trial_fe_space(entity::P.PhysicsParameters, V_test, config::PH.TimeDomainConfig)
    _build_fe_trial_space(entity.fe, V_test, config)
end

"""
    build_trial_fe_space(resn::Vector{P.ResonatorSingle}, Vs::Vector, config::PH.SimulationConfig)

Build trial spaces for each resonator's test space.
"""
build_trial_fe_space(resn::Vector{P.ResonatorSingle}, Vs::Vector, ::PH.FreqDomainConfig) = [TrialFESpace(V) for V in Vs]
build_trial_fe_space(resn::Vector{P.ResonatorSingle}, Vs::Vector, ::PH.TimeDomainConfig) = [TransientTrialFESpace(V) for V in Vs]

# ─────────────────────────────────────────────────────────────
# Main manager function
# ─────────────────────────────────────────────────────────────


"""
    build_fe_spaces(entities, trians::G.TankTriangulations, config::PH.SimulationConfig) -> (X, Y, fmap)

Build multi-field test and trial FE spaces from a list of entities and a TankTriangulations dictionary.
Each entity uses the triangulation from trians[entity.domain_symbol].

# Arguments
- `entities`: iterable containing `P.PhysicsParameters` and/or `Vector{P.ResonatorSingle}`
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
function build_fe_spaces(entities, 
                         trians::G.TankTriangulations, 
                         config::PH.SimulationConfig)
    test_spaces  = SingleFieldFESpace[]
    trial_spaces = SingleFieldFESpace[]
    fmap = Dict{Symbol, Int}()
    idx  = 0

    for (ientity, entity) in enumerate(entities)
        if entity isa Vector{P.ResonatorSingle}
            isempty(entity) && throw(ArgumentError("Resonator array at entities[$ientity] must be non-empty."))
            domain_symbol = entity[1].space_domain_symbol
            all(r -> r.space_domain_symbol == domain_symbol, entity) ||
                throw(ArgumentError("Resonator array at entities[$ientity] has inconsistent `space_domain_symbol` values: $(unique(getfield.(entity, :space_domain_symbol)))"))
            trian = trians[domain_symbol]
            Vs = build_test_fe_space(entity, trian, config)
            Us = build_trial_fe_space(entity, Vs, config)
            for (i, (Vi, Ui)) in enumerate(zip(Vs, Us))
                idx += 1
                push!(test_spaces, Vi)
                push!(trial_spaces, Ui)
                fmap[Symbol("q_$i")] = idx
            end
        elseif entity isa P.PhysicsParameters
            trian = trians[entity.space_domain_symbol]
            for (sym, fe_cfg) in zip(P.variable_symbols(entity),
                                     P.field_fe_configs(entity))
                V = _build_fe_test_space(fe_cfg, trian, config)
                U = _build_fe_trial_space(fe_cfg, V, config)
                idx += 1
                push!(test_spaces, V)
                push!(trial_spaces, U)
                fmap[sym] = idx
            end
        else
            throw(ArgumentError("Unsupported entity type in `build_fe_spaces`: $(typeof(entity))"))
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
