module PlasmaDispersionFunctions

using ContinuedFractions, SeriesAccelerators, SpecialFunctions
import SpecialFunctions.erfcx

export plasma_dispersion_function

"""
Return the value of the plasma dispersion function
This implementation includes the residue, which is easy to verify
because Z(0) = im sqrt(π).

 - z is the argument to the plasma disperion function
 - power is the moment of the integral

[1] S.D. Baalrud, Phys. Plasmas 20, 012118 (2013) and put ν = -Inf
[2] M. Sampoorna et al., Generalized Voigt functions and their derivatives,
  Journal of Quantitative Spectroscopy & Radiative Transfer (2006),
  doi:10.1016/j.jqsrt.2006.08.011
"""
function plasma_dispersion_function(z::T, power::Int=0) where {T<:Number}
  @assert power >= 0 "power, $power, must be >= 0"
  @assert isfinite(z) "z = $z"
  if power > 0
    Z_1 = plasma_dispersion_function(z, power - 1)
    numerator(i) = iseven(i) ? 0 : prod(1:2:(i - 2))
    denominator(i) = 2^((i - 1) / 2)
    return z * Z_1 + numerator(power) / denominator(power)
  end
  return im * (sqrt(π) * erfcx(-im * z)) # the parentheses have to be here!
end


logpochhammer(x, n) = log(gamma(x + n)) - log(gamma(x))

"""
erfcx for large real component e.g. Complex{BigFloat} arguments as described in
https://dlmf.nist.gov/7.12#E1
"""
erfcx(z::Number) = _erfcx(z) # need this to avoid weird dispatch method error
erfcx(z::Complex{BigFloat}) = _erfcx(z)
function _erfcx(z::Number)
  output, isconverged = abs(z) < 3 ? erfcx_smallarg(z) : erfcx_largearg(z)
  isconverged && return output
  output, isconverged = abs(z) < 3 ? erfcx_largearg(z) : erfcx_smallarg(z)
  isconverged && return output
  error("erfcx not converged for argument $z")
end
function erfcx_smallarg(z::T) where {T<:Number}
  y, isconverged = pFq(1, 1.5, z^2)
  return exp(z^2) - 2 * z / sqrt(π) * y, isconverged
end
function erfcx_largearg(z::T, rtol=eps(real(T))^(1/4)) where {T<:Number}
  summand(m, z) = (-1)^m * exp(logpochhammer(0.5, m) - (2m + 1) * log(z))
  if real(z) > 0
    output, isconverged = shanks(i->summand(i, z), 3, rtol=rtol)
    return output / sqrt(π), isconverged
  else
    output, isconverged = shanks(i->summand(i, -z), 3, rtol=rtol)
    return (2 * exp(z^2) - output / sqrt(π)), isconverged
  end
end

function pFq(a::S, b::U, z::V, rtol=sqrt(eps(real(promote_type(S, U, V))))
    ) where {S<:Number, U<:Number, V<:Number}
  numerator(i) = - z * prod(i .+ a) / (i + 1) / prod(i .+ b)
  denominator(i) = 1 - numerator(i)
  K, isconverged = continuedfraction(denominator, numerator, rtol=rtol)
  K -= denominator(0)
  return 1 + z * prod(a) / prod(b) / (1 + K), isconverged
end

end
