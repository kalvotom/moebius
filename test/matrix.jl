###
#
# Discretized matrix tests.
#

@testset "matrixM is available and gives sensible result..." begin
  a = 1.37
  R = 13.78
  indices = Moebius.basis_indices(R, 5)

  M = Moebius.M(indices, R, a)

  @test M ≈ M'
end