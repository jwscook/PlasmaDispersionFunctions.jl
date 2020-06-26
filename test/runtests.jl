using PlasmaDispersionFunctions
using SpecialFunctions, QuadGK, Test, Random

@testset "Tests" begin
  Random.seed!(3)

@testset "Plasma dispersion function" begin
  @test plasma_dispersion_function(0.0, 0) ≈ im*sqrt(pi) rtol=0.001
  @test plasma_dispersion_function(im, 0) ≈ im*0.757872156141312 rtol=0.001
  @test plasma_dispersion_function(ComplexF64(-1.52, 0.47), 0) ≈ ComplexF64(0.6088888957234254, 0.33494583882874024) rtol=0.001
end

@testset "erfcx" begin
  for A ∈ (1e0, 1e1, 1e2, 1e3), i in 1:10
    z = A * ((rand() - 0.5) + im * (rand() - 0.5))
    isfinite(erfcx(z)) || continue
    t = @elapsed erfcxbigz = PlasmaDispersionFunctions.erfcx(big(z))
    @test t < 3
    @test erfcxbigz ≈ erfcx(z) rtol=1e3eps()
  end
end

@testset "vs quadrature" begin
  function test_quad(Z, pow)
    for z in Z
      principal(x) = x.^pow .* exp(-x.^2)/sqrt(pi)
      integrand(x) = principal(x) / (x-z)
      principalpart = if iszero(imag(z))
        folded(v) = (principal(v + real(z)) - principal(-v + real(z))) / v
        QuadGK.quadgk(folded, eps(), 6 + 6*abs(z), rtol=eps())[1]
      else
        QuadGK.quadgk(integrand, -10*abs(z), 10*abs(z), rtol=eps())[1]
      end
      σ = imag(z) < 0 ? 2 : imag(z) == 0 ? 1 : 0
      cauchyresidue = im * (σ * π * principal(z))
      b = PlasmaDispersionFunctions.plasma_dispersion_function(z, pow)
      a = principalpart + cauchyresidue
      @test real(a) ≈ real(b) rtol=1.0e-8
      @test imag(a) ≈ imag(b) rtol=1.0e-8 atol=1.0e-8
    end
  end
  Zs = []
  push!(Zs,  1.0 + im*0)
  push!(Zs, -1.0 + im*0)
  push!(Zs,  1.0 + im/4)
  push!(Zs, -1.0 + im/4)
  push!(Zs,  1.0 - im/4)
  push!(Zs, -1.0 - im/4)
  @testset "QuadGK of maxwellian with 0th moment" begin
    test_quad(Zs, 0)
  end
  @testset "QuadGK of maxwellian with 1st moment" begin
    test_quad(Zs, 1)
  end
  @testset "QuadGK of maxwellian with 2nd moment" begin
    test_quad(Zs, 2)
  end
  @testset "QuadGK of maxwellian with 10th moment" begin
    test_quad(Zs, 10)
  end
end

end
