###
#
# Main file with the `Moebius` module.
#


"""
Module with functions related to the numerical computations
of eigenvalues and eigenfunctions on Möbius strip.
"""
module Moebius
	include("plot.jl")
	include("basis.jl")
	include("fake_and_not_so_fake.jl")
	include("matrix.jl")
end # module
