###
#
# Eigenvalues of the "fake" and "not-so-fake" model.
#

using  MathieuFunctions

"""
  fake_eigenvalues(R, a, maxE)

Eigenvalues of the fake model, i.e. the flat Möbius strip of width ``2a``
and length ``R``, smaller then ``maxE``. The spectrum is explicitely given by

```math
  \\Big\\{ \\big(m/2R\\big)^2 + \\big(n\\pi/2a\\big)^2 \\,\\Big\\vert\\, m\\in\\mathbb{Z}, \\ n\\in\\mathbb{N}^*, \\ m+n \\ \\text{is odd} \\Big}.
```
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
  partA(R, a, maxE; kmax=20, kadd=5)

Compute the part of the not-so-fake spectrum originating from the
Mathieu characteristic value ``a_m``.
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
  partB(R, a, maxE; kmax=20, kadd=5)

Compute the part of the not-so-fake spectrum originating from the
Mathieu characteristic value ``b_m``.
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
  not_so_fake_eigenvalues(R, a, maxE)

Eigenvalues of the not-so-fake model, i.e. the flat Möbius strip
with effective potential ``V_{\\text{eff}}(s,t) = -\\frac{1}{8R^2}\\cos(s/R)``
of width ``2a`` and length ``R``, smaller then ``maxE``. The spectrum
is explicitely given by

```math
  \\Big\\{ (1/2R)^2 a_m(-1/4) + (n\\pi/2a)^2 \\,\\Big\\vert\\, m\\in\\mathbb{N}, \\ n\\in\\mathbb{N}^*, \\ m+n \\ \\text{is odd} \\Big}
  \\cup
  \\Big\\{ (1/2R)^2 b_m(-1/4) + (n\\pi/2a)^2 \\,\\Big\\vert\\, m\\in\\mathbb{N}^*, \\ n\\in\\mathbb{N}^*, \\ m+n \\ \\text{is odd} \\Big},
```

where ``a_m`` and ``b_m`` are the Mathieu characteristic values.
"""
function not_so_fake_eigenvalues(R, a, maxE)
  vals = vcat(partA(R, a, maxE), partB(R, a, maxE))

  return sort(vals)
end
