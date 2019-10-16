###
#
# Example no. 1
#

@info "Running the computation for example 01."

# Setup files
BUILD_DIR        = "build" # if you edit this you have to change Makefile too
MATRIX_CACHE     = joinpath(BUILD_DIR, "ex01_matrix.jld")
TABLE_OF_RESULTS = joinpath(BUILD_DIR, "ex01_results.tex")
PICTURE_2D       = j -> joinpath(BUILD_DIR, "ex01_eigenvector_2d_$j.tex")
PICTURE_3D       = j -> joinpath(BUILD_DIR, "ex01_eigenvector_3d_$j.tex")

# Setup computation parameters
SHOW_RESULTS = 20
a    = 1.3
L    = 18.0
R    = L / 2pi
maxE = 40

include("eigenvalues_and_eigenvectors.jl")
