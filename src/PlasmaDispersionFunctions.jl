module PlasmaDispersionFunctions

using SpecialFunctions

export plasma_dispersion_function

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
function plasma_dispersion_function(z::T, n::Int=0) where {T<:Number}
  n < 0 && throw(ArgumentError("n, $n, must be >= 0"))
  @assert isfinite(z) "z = $z"
  if n > 0
    Z_1 = plasma_dispersion_function(z, n - 1)
    numerator(i) = iseven(i) ? 0 : prod(1:2:(i - 2))
    denominator(i) = 2^((i - 1) / 2)
    return z * Z_1 + numerator(n) / denominator(n)
  end
  return im * (sqrt(π) * erfcx(-im * z)) # the parentheses have to be here!
end

end
