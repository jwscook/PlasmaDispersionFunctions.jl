"""
Calculate the plasma dispersion function

```
principal(n, x) = x.^n .* exp(-x.^2) / √π
integrande(z, n, x) = P(x, n) / (x - z)
residue(z, n, x) = (imag(z) < 0 ? 2 : imag(z) == 0 ? 1 : 0) * im * π * principal(n, z)
Z(z, n) = ∫integrand(z, n, x)dx + residue(z, n, x)
```

# Example
```jldoctest
  plasmadispersionfunction(0.0) ≈ im * √π
```
"""
module PlasmaDispersionFunctions

using LinearAlgebra, QuadGK, SpecialFunctions

export plasmadispersionfunction
export generalisedplasmadispersionfunction

"""@docs
    plasmadispersionfunction(z::T,power::Unsigned=UInt64(0),Z₋₁=missing)

Return the value of the plasma dispersion function
This implementation includes the residue, which is easy to verify
because Z(0) = im sqrt(π).

...
# Arguments
- `z::T`: the argument to the plasma dispersion function
- `power::Unsigned=UInt64(0)`: the moment of the integral
- `Z₋₁=missing`: (optional) the value of the plasma dispersion function for
power-1. If already calculated, passing this in saves some time.
...

# Example
```julia
  plasmadispersionfunction(1.0)
```

# References
[1] D. Fried and S. D. Conte, "The Plasma Dispersion Relation", Elsevier, 1961
[2] S.D. Baalrud, Phys. Plasmas 20, 012118 (2013) and put ν = -Inf
[3] M. Sampoorna et al., Generalized Voigt functions and their derivatives,
  Journal of Quantitative Spectroscopy & Radiative Transfer (2006),
  doi:10.1016/j.jqsrt.2006.08.011
"""
function plasmadispersionfunction(z::T, power::Unsigned=UInt64(0), Z₋₁=missing
    ) where {T<:Number}
  function _const(i::Signed)
    i == 1 && return T(1) # make sure that 1 is handled quickly
    return iseven(i) ? T(0) : T(prod(1:2:(i - 2)) * T(2)^div(1 - i, 2))
  end
  if ismissing(Z₋₁)
    Z̄ = sqrt(π) * erfcx(Complex(imag(z), -real(z))) # hot loop; -im * z
    Z = Complex(-imag(Z̄), real(Z̄)) # hot loop; im * Z̄
    for i ∈ 1:power
      Z = z * Z + _const(Signed(i))
    end
    return Z
  else
    return z * Z₋₁ + _const(Signed(power))
  end
end

function plasmadispersionfunction(z::T, n::Integer, Z₋₁=missing
    ) where {T}
  n < 0 && throw(ArgumentError("n, $n, must be >= 0"))
  return plasmadispersionfunction(z::T, Unsigned(n), Z₋₁)
end

"""
Takes a function that when integrated between -Inf and +Inf returns value x,
and returns a new function that returns x when integrated between real(pole)
and +Inf.
"""
function foldnumeratoraboutpole(f::T, pole::Real) where {T<:Function}
  folded(v::Real) = (f(v + pole) - f(-v + pole)) / v
  return folded
end
function foldnumeratoraboutpole(f::T, pole::Number) where {T<:Function}
  r, i = reim(pole)
  function folded(v::Real)
    a, b, c = f(r + v), f(r - v), 1 / (v - Complex(0, i))
    return (a - b) * real(c) + (a + b) * Complex(0, imag(c))
  end
  return folded
end

"""
Transform the limits of an integrand
quadrature(foldnumeratoraboutpole(integrand, pole), limitsfolder(limits, pole)...)
"""
function limitsfolder(ls::AbstractVector{T}, pole) where {T}
  U = promote_type(T, typeof(pole))
  return vcat(abs.(U.(ls).- real(pole)), isreal(pole) ? eps(U) : zero(U))
end

struct WaveDirectionalityHandler{T}
  kz::T
end
function (wdh::WaveDirectionalityHandler)(x)
  # this way works with DualNumbers
  return real(x) + im * (real(wdh.kz) < 0 ? -imag(x) : imag(x))
end

residuesigma(pole::Number) = imag(pole) < 0 ? 2 : imag(pole) == 0 ? 1 : 0

function residue(numerator::F, pole::Number) where {F}
  principalpart = numerator(pole)
  σ = residuesigma(pole)
  output = im * (σ * π * principalpart)
  iszero(σ) && return zero(output) # defend against overflow
  return output
end

quadnorm(x) = maximum(norm.(x))

function frtol(z::T) where {T<:Number}
  U = float(real(T))
  return max(min(sqrt(eps(U)), angle(z)), eps(U))
end
frtol(z::T) where {T<:Real} = sqrt(eps(float(T)))
function generalisedplasmadispersionfunction(f::F, pole::Number, kz::Number=1;
    atol=eps(), rtol=frtol(pole), quadorder=8, quadnorm::F1=quadnorm
    ) where {F, F1}
  integrand = foldnumeratoraboutpole(f, pole) # also divide by (v - pole)
  principal = first(quadgk(integrand, 0, Inf, rtol=rtol, atol=atol,
                           order=quadorder, norm=quadnorm))
  wdh = WaveDirectionalityHandler(kz)
  residueatpole = wdh(residue(f, wdh(pole)))
  return Complex(principal) + residueatpole
end

end

