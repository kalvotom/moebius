###
#
# Example no. 2
#

@info "Running the computation for example 02."

# Setup files
BUILD_DIR        = "build" # if you edit this you have to change Makefile too
MATRIX_CACHE     = joinpath(BUILD_DIR, "ex02_matrix.jld")
TABLE_OF_RESULTS = joinpath(BUILD_DIR, "ex02_results.tex")
PICTURE_2D       = j -> joinpath(BUILD_DIR, "ex02_eigenvector_2d_$j.tex")
PICTURE_3D       = j -> joinpath(BUILD_DIR, "ex02_eigenvector_3d_$j.tex")

# Setup computation parameters
SHOW_RESULTS = 20
a    = 0.75
L    = 13.2
R    = L / 2pi
maxE = 45

include("eigenvalues_and_eigenvectors.jl")
