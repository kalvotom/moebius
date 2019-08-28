###
#
# Eigenvalues and eigenfunction of the Laplace-Beltrami
# operator on the Möbius strip.
#

using LinearAlgebra
using JLD, Printf

include("../lib/moebius.jl")

# Setup files
BUILD_DIR        = "build" # if you edit this you have to change Makefile too
MATRIX_CACHE     = joinpath(BUILD_DIR, "matrix.jld")
TABLE_OF_RESULTS = joinpath(BUILD_DIR, "results.tex")
PICTURE_2D       = j -> joinpath(BUILD_DIR, "eigenvector_2d_$j.tex")
PICTURE_3D       = j -> joinpath(BUILD_DIR, "eigenvector_3d_$j.tex")

if !isdir(BUILD_DIR)
  mkdir(BUILD_DIR)
end

# Setup parameters
SHOW_RESULTS = 20
a    = 1.3
L    = 18.0
R    = L / 2pi
maxE = 40

#
# The actual computation
#

indices = Moebius.basis_indices(R, maxE)
N = size(indices)[1]

if isfile(MATRIX_CACHE)
  @info "Loading previously computed matrix..."
  m = load(MATRIX_CACHE, "m")
else
  @info "Computing the matrix..."
  m = Moebius.M(indices, R, a, atol=1e-14)
  save(MATRIX_CACHE, "m", m)
end

@info "Computing its eigenvalues and eigenvectors..."
eigf = eigen(Symmetric(m))
eves = eigf.vectors
evas = eigf.values

@info "Writing out the table of results ($TABLE_OF_RESULTS)..."
if N < SHOW_RESULTS
  error("The cut off basis is too small to obtain required $SHOW_RESULTS eigenvalues!")
end

open(TABLE_OF_RESULTS, "w") do f
  for j in 1:SHOW_RESULTS
    print(f, j)
    print(f, " & ")
    print(f, evas[j])
    print(f, " & ")
    print(f, Moebius.getResidue(eves[:, j], evas[j], indices, R, a))
    print(f, "\\\\\n")
  end
end

@info "Generating plot data (2D)..."

for j in 1:SHOW_RESULTS
  Moebius.data_plot(
    Moebius.plotData(eves[:, j], indices, R, a),
    PICTURE_2D(j),
    title="\$\\tilde{\\lambda}^{(\\mathrm{true})}_{$(j)} = $(@sprintf "%.6f" evas[j])\$"
  )
end

@info "Generating plot data (3D)..."

for j in 1:SHOW_RESULTS
  Moebius.data_plot_3d(
    Moebius.plotData(eves[:, j], indices, R, a), R, a,
    PICTURE_3D(j),
    title="\$\\tilde{\\lambda}^{(\\mathrm{true})}_{$(j)} = $(@sprintf "%.6f" evas[j])\$"
  )
end
