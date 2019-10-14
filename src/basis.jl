##
#
# Methods related to the construction of a particular
# orthonormal basis in ``L^2((0, 2πR)×(-1,1))``.
#

"""
    χ(n)

Element of orthonormal basis in ``L^2((-1, 1))``, ``n`` is
a positive integer index.
"""
function χ(n)
  if isodd(n)
    return u -> cos(n * pi * u / 2)
  else
    return u -> sin(n * pi * u / 2)
  end
end

"""
    φ(m, R)

Element of orthonormal basis in ``L^2((0, 2πR))``, ``m`` is
a integer index.
"""
function φ(m, R)
  if m == 0
    return s -> sqrt(1/(4*pi*R))
  elseif m > 0
    return s -> sqrt(1/(2*pi*R)) * cos(m * s / (2R))
  else
    return s -> sqrt(1/(2*pi*R)) * sin(m * s / (2R))
  end
end

"""
  ψ(mn, R)

Element of orthonormal basis in ``L^2((0, 2πR)×(-1,1))``,
in particular ``φ(m, R)⊗χ(n)`` where `(m, n) = mn`.
"""
function ψ(mn, R)
  m = mn[1]; n = mn[2]
  return x -> sqrt(2) * φ(m, R)(x[1]) * χ(n)(x[2])
end

"""
  basis_indices(R, maxE)

Find all indices ``(m, n)`` such that ``m + n`` is odd and
``(nπ/2)^2 + (m/2R)^2 < maxE``. The corresponding set of
tensor product functions ``φ(m, R)⊗χ(n)`` is then a orthonormal
set in ``L^2((0, 2πR)×(-1,1))``.
"""
function basis_indices(R, maxE)
  n       = 1
  indices = zeros(Int64, 0, 2)

  while (n*pi/2)^2 < maxE
    if isodd(0 + n)
      indices = vcat(indices, [0 n])
    end

    m = 1

    while (n*pi/2)^2 + (m/(2R))^2 < maxE
      if isodd(m + n)
        indices = vcat(indices, [ m n])
        indices = vcat(indices, [-m n])
      end

      m += 1
    end

    n += 1
  end

  @info "Using $(size(indices)[1]) basis elements."
  return indices
end

"""
  φ1(m, R)


"""
function φ1(m, R)
  return s -> GSL.sf_mathieu_se(m, -1/4, s / 2R) / sqrt(2*pi*R)
end

"""
  φ2(m, R)

 
"""
function φ2(m, R)
  return s -> GSL.sf_mathieu_ce(m, -1/4, s / 2R) / sqrt(2*pi*R)
end

"""
  ϕ1(mn, R)


"""
function ϕ1(mn, R)
  m = mn[1]; n = mn[2]
  return x -> sqrt(2) * φ1(m, R)(x[1]) * χ(n)(x[2])
end

"""
  ϕ2(mn, R)

 
"""
function ϕ2(mn, R)
  m = mn[1]; n = mn[2]
  return x -> sqrt(2) * φ2(m, R)(x[1]) * χ(n)(x[2])
end
