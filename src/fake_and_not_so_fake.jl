###
#
# Eigenvalues of the "fake" and "not-so-fake" model.
#

using  GSL

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
  partA(R, a, maxE[, eigenfunctions=false])

Compute the part of the not-so-fake spectrum originating from the
Mathieu characteristic value ``a_m``. If `eigenfunctions` is set to `true`
return corresponding eigenfunctions too.
"""
function partA(R, a, maxE; eigenfunctions=false)
  # This is certainly sub-optimal.
  vals   = []

  m = 0
  n = 1

  valA(m, n) = (1/2R)^2 * GSL.sf_mathieu_a(m, -1/4) + (n*pi/2a)^2

  while valA(m, n) < maxE # loop over m
    while valA(m, n) < maxE # loop over n
      if isodd(m + n)
        if eigenfunctions
          vals = vcat(vals, (valA(m, n), ϕ2((m, n), R)))
        else
          append!(vals, valA(m, n))
        end
      end
      n += 1
    end

    n  = 1
    m += 1
  end

  return vals
end

"""
  partB(R, a, maxE[, eigenfunctions=false])

Compute the part of the not-so-fake spectrum originating from the
Mathieu characteristic value ``b_m``. If `eigenfunctions` is set to `true`
return corresponding eigenfunctions too.
"""
function partB(R, a, maxE; eigenfunctions=false)
  # This is certainly sub-optimal.
  vals   = []

  m = 1
  n = 1

  valB(m, n) = (1/2R)^2 * GSL.sf_mathieu_b(m, -1/4)  + (n*pi/2a)^2

  while valB(m, n) < maxE # loop over m
    while valB(m, n) < maxE # loop over n
      if isodd(m + n)
        if eigenfunctions
          vals = vcat(vals, (valB(m, n), ϕ1((m, n), R)))
        else
          append!(vals, valB(m, n))
        end
      end
      n += 1
    end

    n  = 1
    m += 1
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

"""
"""
function not_so_fake_eigenpairs(R, a, maxE)
  vals = vcat(partA(R, a, maxE, eigenfunctions=true), partB(R, a, maxE, eigenfunctions=true))

  return sort(vals, lt=(x,y)->isless(x[1], y[1]))
end
