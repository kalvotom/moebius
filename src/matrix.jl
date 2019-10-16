###
#
# Discretized matrix construction.
#

using LinearAlgebra
using ProgressMeter, Cubature

"""
"""
function kronecker(mn, kl)
  if mn == kl
    return 1
  else
    return 0
  end
end

"""
    M(indices, R, a[, atol=1e-7])


"""
function M(indices, R, a; atol=1e-7)
  N = size(indices)[1]
  mat = zeros(N, N)

  # stuff
  f(s, t)  = sqrt((1 - t/R*cos(s/(2R)))^2 + (t/(2R))^2)
  fa(s, u) = f(s, a*u)

  fa1(s, u)  = a*u * (R - a*u*cos(s/(2R))) * sin(s/(2R)) / (2*R^3*fa(s, u))
  fa2(s, u)  = a * (-4R*cos(s/(2R)) + a*u*(3+2*cos(s/R))) / (4*R^2*fa(s, u))
  fa11(s, u) = -(a*u)/(32*R^6*fa(s, u)^3) * (
    -4R*(2*R^2 + 5*a^2*u^2)*cos(s/(2R)) +
    a*u*(
      12*R^2 + 3*a^2*u^2 + 6*(2*R^2 + a^2*u^2)*cos(s/R) -
      6*a*R*u*cos(3s/(2R)) + a^2*u^2*cos(2s/R)
    )
  )
  fa22(s, u) = a^2 / (4*R^2) / fa(s, u)^3

  # the potential
  V(s, u) = (
    -5/4 * (fa1(s, u)^2) / fa(s, u)^4 + 1/2 * fa11(s, u) / fa(s, u)^3 +
    -1/4 * (fa2(s, u)^2) / (a^2*fa(s, u)^2) + 1/2 * fa22(s, u) / (a^2 * fa(s, u))
  )

  pm = Progress(Int(N*(N+1)/2), 5)

  for m in 1:N
    for n in m:N
      mn = indices[m, :]; kl = indices[n, :]
      func1 = x -> 1 / fa(x[1], x[2])^2 * ψ([-mn[1], mn[2]], R)(x) * ψ([-kl[1], kl[2]], R)(x)
      funcV = x -> V(x[1], x[2]) * ψ(mn, R)(x) * ψ(kl, R)(x)

      # renormalized matrix: #1
      # mat[m, n] =
      #   mn[1] * kl[1] / (2R)^2 * hcubature(func1, [0.0, -1.0], [2.0*pi*R, 1.0], abstol=atol)[1] + 
      #   ((mn[2]^2 - 1) * pi^2 / (2a)^2) * kronecker(mn, kl) +
      #   hcubature(funcV, [0.0, -1.0], [2.0*pi*R, 1.0], abstol=atol)[1]

      # renormalized matrix: #2
      # mat[m, n] =
      #   a * mn[1] * kl[1] / (2R)^2 * hcubature(func1, [0.0, -1.0], [2.0*pi*R, 1.0], abstol=atol)[1] + 
      #   (mn[2]^2 - 1) * pi^2 / (4a) * kronecker(mn, kl) +
      #   a * hcubature(funcV, [0.0, -1.0], [2.0*pi*R, 1.0], abstol=atol)[1]

      # ordinary matrix
      mat[m, n] =
        mn[1] * kl[1] / (2R)^2 * hcubature(func1, [0.0, -1.0], [2.0*pi*R, 1.0], abstol=atol)[1] + 
        (1 / a^2 * (mn[2] * pi / 2)^2) * kronecker(mn, kl) +
        hcubature(funcV, [0.0, -1.0], [2.0*pi*R, 1.0], abstol=atol)[1]

      if m != n
        mat[n, m] = mat[m, n]
      end

      next!(pm)
    end
  end

  return mat
end

"""
"""
function getValue(eigenVector, indices, R, x)
  val = 0.0
  N = size(indices)[1]

  for m in 1:N
    val += eigenVector[m] * ψ(indices[m, :], R)(x)
  end

  return val
end

"""
"""
function getLValue(eigenVector, indices, R, a, x)
  val = 0.0
  N = size(indices)[1]

  # stuff
  f(s, t)  = sqrt((1 - t/R*cos(s/(2R)))^2 + (t/(2R))^2)
  fa(s, u) = f(s, a*u)

  fa1(s, u)  = a*u * (R - a*u*cos(s/(2R))) * sin(s/(2R)) / (2*R^3*fa(s, u))
  fa2(s, u)  = a * (-4R*cos(s/(2R)) + a*u*(3+2*cos(s/R))) / (4*R^2*fa(s, u))
  fa11(s, u) = -(a*u)/(32*R^6*fa(s, u)^3) * (
    -4R*(2*R^2 + 5*a^2*u^2)*cos(s/(2R)) +
    a*u*(
      12*R^2 + 3*a^2*u^2 + 6*(2*R^2 + a^2*u^2)*cos(s/R) -
      6*a*R*u*cos(3s/(2R)) + a^2*u^2*cos(2s/R)
    )
  )
  fa22(s, u) = a^2 / (4*R^2) / fa(s, u)^3

  # the potential
  V(s, u) = (
    -5/4 * (fa1(s, u)^2) / fa(s, u)^4 + 1/2 * fa11(s, u) / fa(s, u)^3 +
    -1/4 * (fa2(s, u)^2) / (a^2*fa(s, u)^2) + 1/2 * fa22(s, u) / (a^2 * fa(s, u))
  )

  # computation

  for m in 1:N
    mn = indices[m, :]
    val +=  eigenVector[m] * (
      2/fa(x[1], x[2])^3 * fa1(x[1], x[2]) * mn[1] / (2R) * ψ([-mn[1], mn[2]], R)(x) +
      1/fa(x[1], x[2])^2 * (mn[1]/(2R))^2 * ψ(mn, R)(x) +
      1/a^2 * (mn[2] * pi / 2)^2 * ψ(mn, R)(x) +
      V(x[1], x[2]) * ψ(mn, R)(x)
    )
  end

  return val
end

"""
"""
function getResidue(eigenVector, eigenValue, indices, R, a)
  func = x -> getValue(eigenVector, indices, R, x)^2
  normalization = sqrt(hcubature(func, [0.0, -1.0], [2.0*pi*R, 1.0], abstol=1e-10)[1])

  func = x -> (getLValue(eigenVector, indices, R, a, x)
    - eigenValue * getValue(eigenVector, indices, R, x))^2

  return sqrt(hcubature(func, [0.0, -1.0], [2.0*pi*R, 1.0], abstol=1e-6)[1]) / normalization
end

"""
  compareEigenvectors(f1, f2, R)

Compute ``L^2`` norm of `f1-f2`. `f1` is the result of our numerical
computation it is not expected to be normalized. `f2` is the normalized
eigenvector of the not-so-fake model.

This computation is tricky because of the possible change in sign. Also,
due to tight clustering of pairs of eigenvalues we have to be very careful
not to shuffle the eigenvectors. 
"""
function compareEigenvectors(f1, f2, R)
  # compute the normalization
  func = x -> f1(x)^2
  normalization = sqrt(hcubature(func, [0., -1.0], [2*pi*R, 1.0])[1])

  # attempt to heuristically align both functions
  tries   = 21
  aligned = 0
  sample  = () -> [(0.25 + rand()/2)*2*pi*R, -0.7+1.4*rand()]

  for j in 1:tries
    x = sample()
    if f1(x) * f2(x) > 0
      aligned += 1
    end
  end

  if aligned > 18 # sign is probably the same
    func = x -> (f1(x) / normalization - f2(x))^2
    return(sqrt(hcubature(func, [0., -1.0], [2*pi*R, 1.0])[1]))
  elseif aligned < 3 # sign is probably opposite
    func = x -> (f1(x) / normalization + f2(x))^2
    return(sqrt(hcubature(func, [0., -1.0], [2*pi*R, 1.0])[1]))
  else # not sure :-/, pick the smaller one...
    func = x -> (f1(x) / normalization - f2(x))^2
    a = sqrt(hcubature(func, [0., -1.0], [2*pi*R, 1.0])[1])
    func = x -> (f1(x) / normalization + f2(x))^2
    b = sqrt(hcubature(func, [0., -1.0], [2*pi*R, 1.0])[1])
    return min(a, b)
  end
end


"""
"""
function plotData(eigenVector, indices, R, a, width=120, height=80)
  du = 2 / (height - 1)
  ds = 2*pi*R / (width - 1)
  N  = width * height

  data = zeros(N, 4)
  n = 1

  for j in 1:width
    s = ds * (j - 1)

    for k in 1:height  
      u = -1 + du * (k - 1)

      val = getValue(eigenVector, indices, R, [s, u])

      if k == height
        data[n, :] = [s u val  1.0]
      else
        data[n, :] = [s u val -1.0]
      end

      n += 1
    end
  end

  return data
end