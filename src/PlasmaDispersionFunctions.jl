module PlasmaDispersionFunctions

using SpecialFunctions

export plasma_dispersion_function

"""
    plasma_dispersion_function(z::T,power::Unsigned=UInt64(0),Z₋₁=missing)

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
  plasma_dispersion_function(1.0)
```

# References
[1] S.D. Baalrud, Phys. Plasmas 20, 012118 (2013) and put ν = -Inf
[2] M. Sampoorna et al., Generalized Voigt functions and their derivatives,
  Journal of Quantitative Spectroscopy & Radiative Transfer (2006),
  doi:10.1016/j.jqsrt.2006.08.011
"""
function plasma_dispersion_function(z::T, power::Unsigned=UInt64(0), Z₋₁=missing
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

function plasma_dispersion_function(z::T, n::Integer, Z₋₁=missing
    ) where {T}
  n < 0 && throw(ArgumentError("n, $n, must be >= 0"))
  return plasma_dispersion_function(z::T, Unsigned(n), Z₋₁)
end

end
