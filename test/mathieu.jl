###
#
# Test GSL Mathieu functions against Mathematica.
#

@testset "Mathieu characteristic values..." begin
  # Even Mathieu characteristic values with `q = -1/4`.
  mathA = [
    -0.031039395475617324
    0.7424288259866297
    4.025829084645603
    9.00366486704624
    16.002085290467196
    25.00130213222684
    36.000892873798435
    49.000651047848066
    64.00049603440671
    81.0003906262757
    100.00031565723008
    121.00026041703619
    144.0002185316853
    169.0001860120377
    196.0001602564949
    225.0001395089842
    256.0001225490572
    289.0001085069705
    324.00009674924445
    361.00008680556886
    400.0000783208117
  ]
  juliaA = [GSL.sf_mathieu_a(k, -1/4) for k in 0:(length(mathA)-1)]

  # Odd Mathieu characteristic values with `q = -1/4`.
  mathB = [
    1.2419411282429151
    3.994793078632119
    9.004152551546936
    16.002081901038174
    25.0013021454698
    36.00089287376532
    49.00065104784812
    64.00049603440671
    81.0003906262757
    100.00031565723008
    121.00026041703619
    144.0002185316853
    169.0001860120377
    196.0001602564949
    225.0001395089842
    256.0001225490572
    289.0001085069705
    324.00009674924445
    361.00008680556886
    400.0000783208117
  ]
  juliaB = [GSL.sf_mathieu_b(k, -1/4) for k in 1:length(mathB)]

  @test norm(mathA - juliaA, Inf) ≈ 0 atol=1e-12
  @test norm(mathB - juliaB, Inf) ≈ 0 atol=1e-12
end

@testset "Odd Mathieu characteristic functions: normalization" begin
  for k in 1:6
    func = x -> GSL.sf_mathieu_se(k, -1/4, x)^2
    val, err = hquadrature(func, -pi, pi)
    @test val ≈ pi
  end
end

@testset "Even Mathieu characteristic functions: normalization" begin
  for k in 0:5
    func = x -> GSL.sf_mathieu_ce(k, -1/4, x)^2
    val, err = hquadrature(func, -pi, pi)
    @test val ≈ pi
  end
end
