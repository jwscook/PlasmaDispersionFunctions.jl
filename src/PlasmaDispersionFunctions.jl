module PlasmaDispersionFunctions

using SpecialFunctions

export plasma_dispersion_function

function constant(i::Unsigned, ::Type{T}) where {T}
  iseven(i) && return  T(0)
  isone(i) && return T(1)
  return T(prod(1:2:(i - 2)) / T(2)^div(i - 1, 2))
end
"""
Return the value of the plasma dispersion function
This implementation includes the residue, which is easy to verify
because Z(0) = im sqrt(π).

 - z is the argument to the plasma dispersion function
 - n is the moment of the integral

[1] S.D. Baalrud, Phys. Plasmas 20, 012118 (2013) and put ν = -Inf
[2] M. Sampoorna et al., Generalized Voigt functions and their derivatives,
  Journal of Quantitative Spectroscopy & Radiative Transfer (2006),
  doi:10.1016/j.jqsrt.2006.08.011
"""
function plasma_dispersion_function(z::T, n::Unsigned=UInt64(0)) where {T}
  iszero(n) && return im * (sqrt(π) * erfcx(-im * z))
  return z * plasma_dispersion_function(z, n - 1) + constant(n, T)
end
function plasma_dispersion_function(z::T, n::Integer) where {T}
  n < 0 && throw(ArgumentError("n, $n, must be >= 0"))
  return plasma_dispersion_function(z::T, Unsigned(n))
end

end
