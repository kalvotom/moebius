###
#
# Eigenvalues of the "fake" and "not-so-fake" model.
#

using  MathieuFunctions

"""
"""
function fake_eigenvalues(R, a, maxE)
  vals = []
  n    = 1
  m    = 0

  val(m, n) = (m/2R)^2 + (n*pi/2a)^2

  while val(m, n) < maxE # loop over n
    while val(m, n) < maxE # loop over m
      if isodd(m + n)
        if m == 0
          append!(vals, val(m, n))
        else
          append!(vals, val( m, n)) # negative m too
          append!(vals, val(-m, n)) # (degeneracy)
        end
      end
      m += 1
    end

    n += 1
    m  = 0
  end

  return sort(vals)
end

"""
"""
function partA(R, a, maxE; kmax=20, kadd=5)
  # This is certainly sub-optimal.
  # Also, see https://github.com/BBN-Q/MathieuFunctions.jl/issues/5
  kcut   = kmax
  charAs = charA(-1/4, k=0:kcut)
  vals   = []

  m = 0
  n = 1

  valA(m, n) = (1/2R)^2 * charAs[m + 1]  + (n*pi/2a)^2

  while valA(m, n) < maxE # loop over m
    while valA(m, n) < maxE # loop over n
      if isodd(m + n)
        append!(vals, valA(m, n))
      end
      n += 1
    end

    n  = 1
    m += 1
    
    if m + 1 > length(charAs)
      kcut += kadd
      charAs = charA(-1/4, k=0:kcut)
    end
  end

  return vals
end

"""
"""
function partB(R, a, maxE; kmax=20, kadd=5)
  # This is certainly sub-optimal.
  # Also, see https://github.com/BBN-Q/MathieuFunctions.jl/issues/5
  kcut   = kmax
  charBs = charB(-1/4, k=1:kcut)
  vals   = []

  m = 1
  n = 1

  valB(m, n) = (1/2R)^2 * charBs[m]  + (n*pi/2a)^2

  while valB(m, n) < maxE # loop over m
    while valB(m, n) < maxE # loop over n
      if isodd(m + n)
        append!(vals, valB(m, n))
      end
      n += 1
    end

    n  = 1
    m += 1
    
    if m > length(charBs)
      kcut += kadd
      charBs = charB(-1/4, k=1:kcut)
    end
  end

  return vals
end

"""
"""
function not_so_fake_eigenvalues(R, a, maxE)
  vals = vcat(partA(R, a, maxE), partB(R, a, maxE))

  return sort(vals)
end
