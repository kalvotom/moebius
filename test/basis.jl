##
#
# Tests orthonormality.
#

@testset "Checking orthonormality of χ basis..." begin
  N = 10

  matB = zeros(N, N)
  matI = diagm(0 => ones(N))

  for m in 1:N
    for n in 1:N
      func = u -> Moebius.χ(m)(u) * Moebius.χ(n)(u)
      (val, err) = hquadrature(func, -1, 1, abstol=1e-8)
      matB[m, n] = val
    end
  end

  @test norm(matB - matI, Inf) ≈ 0.0 atol=1e-10
end

@testset "Checking orthonormality of φ basis..." begin
  R = 1.3; N = 5

  matB = zeros(2N+1, 2N+1)
  matI = diagm(0 => ones(2N+1))

  for m in -N:N
    for n in -N:N
      func = u -> Moebius.φ(m, R)(u) * Moebius.φ(n, R)(u)
      (val, err) = hquadrature(func, -2*pi*R, 2*pi*R, abstol=1e-8)
      matB[m+N+1, n+N+1] = val
    end
  end

  @test norm(matB - matI, Inf) ≈ 0.0 atol=1e-10
end

@testset "Checking orthonormality of φ⊗χ basis..." begin
  R = 1.3; N = 10

  indices = Moebius.basis_indices(R, N)
  bN = size(indices)[1]

  matB = zeros(bN, bN)
  matI = diagm(0 => ones(bN))

  for m in 1:bN
    for n in 1:bN
      func = x -> Moebius.ψ(indices[m, :], R)(x) * Moebius.ψ(indices[n, :], R)(x)
      (val, err) = hcubature(func, [0.0, -1.0], [2.0*pi*R, 1.0], abstol=1e-10)
      matB[m, n] = val
    end
  end

  @test norm(matB - matI, Inf) ≈ 0.0 atol=1e-10
end

@testset "Checking orthonormality of ϕ1 ∪ ϕ2 basis..." begin
  R = 1.3; N = 10

  basis = vcat(
    [ Moebius.ϕ1((m, n), R) for m in 1:N for n in 1:N if isodd(m+n) ],
    [ Moebius.ϕ2((m, n), R) for m in 0:N for n in 1:N if isodd(m+n) ],
  )

  bN = length(basis)
  matB = zeros(bN, bN)
  matI = diagm(0 => ones(bN))

  for m in 1:bN
    for n in 1:bN
      func = x -> basis[m](x) * basis[n](x)
      (val, err) = hcubature(func, [0.0, -1.0], [2.0*pi*R, 1.0], abstol=1e-10)
      matB[m, n] = val
    end
  end

  @test norm(matB - matI, Inf) ≈ 0.0 atol=1e-10
end
