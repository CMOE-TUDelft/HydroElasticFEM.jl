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

import ..Entities as E 

# ─────────────────────────────────────────────────────────────
# Per-entity FE space construction (dispatched on entity type)
# ─────────────────────────────────────────────────────────────

"""
    build_test_fe_space(entity::PhysicsParameters, trian)

Build a Gridap `TestFESpace` for `entity` on triangulation `trian`,
using the numerical parameters from `entity.fe`.
"""
function build_test_fe_space(entity::E.PhysicsParameters, trian)
    fe = entity.fe
    reffe = ReferenceFE(fe.reffe_type, fe.space_type, fe.order)
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
    build_test_fe_space(resn::Vector{E.ResonatorSingle}, trian)

Build one `ConstantFESpace` per resonator on triangulation `trian`.
Returns a `Vector` of test spaces.
"""
function build_test_fe_space(resn::Vector{E.ResonatorSingle}, trian)
    [ConstantFESpace(trian;
        vector_type = Vector{ComplexF64},
        field_type  = VectorValue{1, ComplexF64})
     for _ in resn]
end

"""
    build_trial_fe_space(entity::E.PhysicsParameters, V_test)

Build a `TrialFESpace` from the test space `V_test`, applying
Dirichlet values from `entity.fe` if present.
"""
function build_trial_fe_space(entity::E.PhysicsParameters, V_test)
    fe = entity.fe
    if fe.dirichlet_value !== nothing
        TrialFESpace(V_test, fe.dirichlet_value)
    else
        TrialFESpace(V_test)
    end
end

"""
    build_trial_fe_space(resn::Vector{E.ResonatorSingle}, Vs::Vector)

Build trial spaces for each resonator's test space.
"""
function build_trial_fe_space(resn::Vector{E.ResonatorSingle}, Vs::Vector)
    [TrialFESpace(V) for V in Vs]
end

# ─────────────────────────────────────────────────────────────
# Main manager function
# ─────────────────────────────────────────────────────────────

"""
    build_fe_spaces(pairs::Pair...) -> (X, Y, fmap)

Build multi-field test and trial FE spaces from entity => triangulation
pairs.  Automatically constructs per-entity spaces using each entity's
`FESpaceConfig`, bundles them into `MultiFieldFESpace`, and builds the
`FieldDict`-compatible field-index mapping.

# Arguments
- `pairs`: sequence of `entity => triangulation`.  Entity may be any
  `E.PhysicsParameters` subtype or a `Vector{E.ResonatorSingle}`.

# Returns
- `X::MultiFieldFESpace` — trial space
- `Y::MultiFieldFESpace` — test space
- `fmap::Dict{Symbol,Int}` — field-symbol to positional-index map

# Example
```julia
X, Y, fmap = build_fe_spaces(
    fluid => Ω,
    fsurf => Γκ,
    mem   => Γη,
)
# fmap == Dict(:ϕ => 1, :κ => 2, :η_m => 3)
```
"""
function build_fe_spaces(pairs::Pair...)
    test_spaces  = SingleFieldFESpace[]
    trial_spaces = SingleFieldFESpace[]
    fmap = Dict{Symbol, Int}()
    idx  = 0

    for (entity, trian) in pairs
        if entity isa Vector{E.ResonatorSingle}
            Vs = build_test_fe_space(entity, trian)
            Us = build_trial_fe_space(entity, Vs)
            for (i, (Vi, Ui)) in enumerate(zip(Vs, Us))
                idx += 1
                push!(test_spaces, Vi)
                push!(trial_spaces, Ui)
                fmap[Symbol("q_$i")] = idx
            end
        else
            V = build_test_fe_space(entity, trian)
            U = build_trial_fe_space(entity, V)
            idx += 1
            push!(test_spaces, V)
            push!(trial_spaces, U)
            fmap[E.variable_symbol(entity)] = idx
        end
    end

    Y = MultiFieldFESpace(test_spaces)
    X = MultiFieldFESpace(trial_spaces)
    return X, Y, fmap
end

# ─────────────────────────────────────────────────────────────
# Exports
# ─────────────────────────────────────────────────────────────

export FESpaceConfig
export build_fe_spaces, build_test_fe_space, build_trial_fe_space

end # module FESpaceAssembly
