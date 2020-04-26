using PlasmaDispersionFunctions
using SpecialFunctions

@testset "Tests" begin

@testset "Plasma dispersion function" begin
  @test plasma_dispersion_function(0.0, 0) ≈ im*sqrt(pi) rtol=0.001
  @test plasma_dispersion_function(im, 0) ≈ im*0.757872156141312 rtol=0.001
  @test plasma_dispersion_function(ComplexF64(-1.52, 0.47), 0) ≈ ComplexF64(0.6088888957234254, 0.33494583882874024) rtol=0.001
end


@testset "erfcx" begin
  for A ∈ (1e0, 1e1, 1e2, 1e3)
    for i in 1:10
      z = A * ((rand() - 0.5) + im * (rand() - 0.5))
      isfinite(erfcx(z)) || continue
      t = @elapsed erfcxbigz = PlasmaDispersionFunctions.erfcx(big(z))
      @test t < 3
      @test erfcxbigz ≈ erfcx(z) rtol=1000eps()
    end
  end
end

end
