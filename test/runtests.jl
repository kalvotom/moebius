###
#
# Some very basic tests of the code.
#

using Test
using Cubature, LinearAlgebra
using GSL # for Mathieu functions

include("../src/moebius.jl")

include("basis.jl")
include("matrix.jl")
include("mathieu.jl")
