using PlasmaDispersionFunctions
using ArbNumerics, QuadGK, Random, SpecialFunctions, Test

import SpecialFunctions.erfcx
erfcx(x::Complex{BigFloat}) = erfcx(ArbComplex(x))

@testset "Tests" begin
Random.seed!(0)

@testset "Plasma dispersion function" begin
  @test plasma_dispersion_function(0.0) ≈ im*sqrt(pi) rtol=1e-5
  @test plasma_dispersion_function(im) ≈ im*0.757872156141312 rtol=1e-5
  @test plasma_dispersion_function(ComplexF64(-1.52, 0.47), 0) ≈
    ComplexF64(0.6088888957234254, 0.33494583882874024) rtol=1e-5
  @test plasma_dispersion_function(big(0.0+0im)) ≈
    ArbComplex(im*sqrt(pi)) rtol=1e-5
end

@testset "using Z₋₁" begin
  Z0 = plasma_dispersion_function(0.0, 0)
  Z1 = plasma_dispersion_function(0.0, 1)
  @test Z1 == plasma_dispersion_function(0.0, 1, Z0)
  Z2 = plasma_dispersion_function(0.0, 2)
  @test Z2 == plasma_dispersion_function(0.0, 2, Z1)
end

@testset "vs quadrature" begin
  function test_quad(Z, pow)
    for z in Z
      principal(x) = x.^pow .* exp(-x.^2)/sqrt(pi)
      integrand(x) = principal(x) / (x-z)
      principalpart = if iszero(imag(z))
        folded(v) = (principal(v + real(z)) - principal(-v + real(z))) / v
        QuadGK.quadgk(folded, nextfloat(0.0), Inf, rtol=eps())[1]
      else
        QuadGK.quadgk(integrand, -Inf, Inf, rtol=eps())[1]
      end
      σ = imag(z) < 0 ? 2 : imag(z) == 0 ? 1 : 0
      cauchyresidue = im * (σ * π * principal(z))
      b = plasma_dispersion_function(z, pow)
      a = principalpart + cauchyresidue
      @test real(a) ≈ real(b) rtol=sqrt(eps()) atol=0.0
      @test imag(a) ≈ imag(b) rtol=sqrt(eps()) atol=0.0
    end
  end
  Zs = []
  push!(Zs,  1.0)
  push!(Zs, -1.0)
  push!(Zs,  1.0 + im)
  push!(Zs, -1.0 + im)
  push!(Zs,  1.0 - im)
  push!(Zs, -1.0 - im)
  for j in 1:10
    push!(Zs, 4 * (rand() - 0.5) + im * 4 * (rand() - 0.5))
  end
  for i in Unsigned.(0:10)
    @testset "QuadGK of maxwellian with $(i)th moment" begin
      test_quad(Zs, i)
    end
  end
end

@testset "check errors are caught" begin
  @test_throws ArgumentError plasma_dispersion_function(1.0, -1)
end

end
