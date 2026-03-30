"""
    WaveInput_FrequencyDomain

A module for defining various types of wave input functions 
in frequency domain formulation for hydro-elastic analysis.
Provides a consistent interface for different wave theories including Airy waves, 
Stokes waves, and other wave types. 
"""
module WaveInput_FrequencyDomain

using WaveSpec
# using .Constants
using Gridap

export AbstractWave, AiryWave, StokesWave
export surface_elevation, velocity_potential, potential_gradient

# ============================================================================
# Structs Definition
# ============================================================================

"""
    AbstractWave

Abstract base type for all wave input types.
"""
abstract type AbstractWave end


# ---------------------Start---------------------
"""
    AiryWave

Linear Airy wave theory implementation.

# Fields
- `η₀::Real`: Wave amplitude
- `T::Real`: Wave period
- `ω::Real`: Angular frequency  
- `λ::Real`: Linear Wavelength
- `k::Real`: Wave number
- `h::Real`: Still-water depth
- `α::Real`: Phase angle
"""
struct AiryWaveXZ <: AbstractWave
  η0::Real
  T::Real
  ω::Real
  λ::Real
  k::Real
  h::Real # Still-water depth
  α::Real

  # Derived quantities can be computed in the constructor
  kh::Real
  sinh_kh::Real
  cosh_kh::Real

  
  function AiryWaveXZ(h, ω, η0, α=0.0)        
    k = WaveSpec.AiryWaves.solve_wavenumber(ω, h)
    λ = 2π/k   
    T = 2π/ω
    kh = k * h
    sinh_kh = sinh(kh)
    cosh_kh = cosh(kh)
    new(
      η0, T, ω, λ, k, h, α,
      kh, sinh_kh, cosh_kh
    )
  end
end
# ----------------------End----------------------


# ---------------------Start---------------------
"""
    StokesWave

Second-order Stokes wave theory implementation.

# Fields
- `η₀::Real`: Wave amplitude
- `ω::Real`: Angular frequency
- `k::Real`: Wave number  
- `H₀::Real`: Water depth
- `α::Real`: Phase angle
- `g::Real`: Gravitational acceleration (default: 9.81)
"""
struct StokesWave <: AbstractWave
    η₀::Real
    ω::Real
    k::Real
    H₀::Real
    α::Real
    g::Real
    
    function StokesWave(η₀, ω, H₀, α=0.0; g=9.81)
        k = WaveSpec.AiryWaves.solve_wavenumber(ω, H₀)
        new(η₀, ω, k, H₀, α, g)
    end
end
# ----------------------End----------------------


# ============================================================================
# Functions Declarations
# ============================================================================

"""
    surface_elevation(wave::AbstractWave, x)

Compute the surface elevation η(x) for the given wave at position x.
"""
function surface_elevation end

"""
    velocity_potential(wave::AbstractWave, x)

Compute the velocity potential ϕ(x,z) for the given wave at position x=(x,z).
"""
function velocity_potential end

"""
    potential_gradient(wave::AbstractWave, x)

Compute the gradient of the velocity potential ∇ϕ(x,z) for the given wave.
"""
function potential_gradient end


# ============================================================================
# Airy Wave Implementations
# ============================================================================

function surface_elevation(wv::AiryWaveXZ, x::VectorValue{2, Float64})
  return wv.η0*exp( im*wv.k*x[1] + im*wv.α )
end

function velocity_potential(wv::AiryWaveXZ, x::VectorValue{2, Float64})
  return (-im)*( wv.η0 * wv.ω / wv.k )*
    (cosh(wv.k*(wv.h + x[2])) / wv.sinh_kh)*
    exp( im*wv.k*x[1] + im*wv.α )
end

function potential_gradient(wv::AiryWaveXZ, x::VectorValue{2, Float64})
    
  ∂ϕ∂x = (wv.η0 * wv.ω)*
    (cosh(wv.k*(wv.h + x[2])) / wv.sinh_kh)*
    exp( im*wv.k*x[1] + im*wv.α )
           
  ∂ϕ∂z = (-im)*(wv.η0 * wv.ω)*
    (sinh(wv.k*(wv.h + x[2])) / wv.sinh_kh)*
    exp( im*wv.k*x[1] + im*wv.α )
           
  return VectorValue(∂ϕ∂x, ∂ϕ∂z)
end

# ============================================================================
# Stokes Wave Implementations  
# NEED TO VERIFY THESE EXPRESSIONS
# ============================================================================

function surface_elevation(wave::StokesWave, x)
    # # First-order term
    # η₁ = wave.η₀ * exp(im * wave.k * x[1] + im * wave.α)
    
    # # Second-order correction
    # η₂ = (wave.η₀^2 * wave.k / 4) * 
    #      (2 + cosh(2 * wave.k * wave.H₀)) / (sinh(wave.k * wave.H₀))^2 *
    #      exp(2im * wave.k * x[1] + 2im * wave.α)
    
    # return η₁ + η₂
    return 0.0
end

function velocity_potential(wave::StokesWave, x)
    # # First-order potential
    # ϕ₁ = -im * (wave.η₀ * wave.ω / wave.k) * 
    #      (cosh(wave.k * (wave.H₀ + x[2])) / sinh(wave.k * wave.H₀)) *
    #      exp(im * wave.k * x[1] + im * wave.α)
    
    # # Second-order correction
    # ϕ₂ = -(wave.η₀^2 * wave.ω / 8) * 
    #      (3 * cosh(2 * wave.k * (wave.H₀ + x[2])) / sinh(2 * wave.k * wave.H₀) - 1) *
    #      exp(2im * wave.k * x[1] + 2im * wave.α)
    
    # return ϕ₁ + ϕ₂

    return 0.0
end


function potential_gradient(wave::StokesWave, x)
    # # First-order gradients
    # ∂ϕ₁∂x = (wave.η₀ * wave.ω) * 
    #         (cosh(wave.k * (wave.H₀ + x[2])) / sinh(wave.k * wave.H₀)) *
    #         exp(im * wave.k * x[1] + im * wave.α)
            
    # ∂ϕ₁∂z = -im * (wave.η₀ * wave.ω) * 
    #         (sinh(wave.k * (wave.H₀ + x[2])) / sinh(wave.k * wave.H₀)) *
    #         exp(im * wave.k * x[1] + im * wave.α)
    
    # # Second-order corrections
    # ∂ϕ₂∂x = -(wave.η₀^2 * wave.ω * wave.k / 4) * 
    #         (3 * cosh(2 * wave.k * (wave.H₀ + x[2])) / sinh(2 * wave.k * wave.H₀)) *
    #         exp(2im * wave.k * x[1] + 2im * wave.α)
            
    # ∂ϕ₂∂z = -im * (3 * wave.η₀^2 * wave.ω * wave.k / 4) *
    #         (sinh(2 * wave.k * (wave.H₀ + x[2])) / sinh(2 * wave.k * wave.H₀)) *
    #         exp(2im * wave.k * x[1] + 2im * wave.α)
    
    # return (∂ϕ₁∂x + ∂ϕ₂∂x, ∂ϕ₁∂z + ∂ϕ₂∂z)

    return VectorValue(0.0, 0.0)
end


# ============================================================================
# Utility Functions
# ============================================================================

"""
    wave_properties(wave::AbstractWave)

Return a named tuple with basic wave properties.
"""
function wave_properties(wv::AbstractWave)
  return (
    amplitude=wv.η0, frequency=wv.ω, period=wv.T, 
    wavelength=wv.λ, wavenumber=wv.k, depth=wv.h, 
    phase=wv.α, hbyL = wv.h/wv.λ
  )
end

"""
    is_deep_water(wave::AbstractWave; threshold=0.5)

Check if the wave is in deep water conditions (kh > threshold).
"""
function is_deep_water(wave::AbstractWave; threshold=0.5)
    return wave.kh > threshold
end

"""
    is_shallow_water(wave::AbstractWave; threshold=0.05)

Check if the wave is in shallow water conditions (kh < threshold).
"""
function is_shallow_water(wave::AbstractWave; threshold=0.05)
    return wave.kh < threshold
end

end # module WaveInputs
