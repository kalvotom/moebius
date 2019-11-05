##
#
# Experiment with the width of the Möbius strip going to zero.
#

using LinearAlgebra, GenericSVD
using JLD, CSV, DataFrames, Cubature

include("../src/moebius.jl")

# Setup files
BUILD_DIR    = "build" # if you edit this you have to change Makefile too
MATRIX_CACHE = j -> joinpath(BUILD_DIR, "matrix_$j.jld")
RATIOS       = "ratios.csv"

if !isdir(BUILD_DIR)
  mkdir(BUILD_DIR)
end

# Setup computation parameters
# If you modify this you have to change `build/fig_*ratio*.tex`!
SHOW_RESULTS = 20
SHOW_EIGENVECTORS = 7

@info "Computing dependence on 'a'..."

# Möbius strip parameters
L    = 18.0
R    = L / (2pi)
maxE = 30

# Thickness
amax = 1.5
amin = 0.01
anum = 20
da   = (amax - amin) / (anum-1)
a    = amax

# Basis
indices = Moebius.basis_indices(R, maxE)
N       = size(indices)[1]

# Containers holding the results
as        = zeros(anum)
big_evals = big.(zeros(anum, N))
evs_conv  = zeros(anum, SHOW_EIGENVECTORS)

for j in 1:anum
  global a, as, big_evals

  println("Step $j of $anum.")
  
  file = MATRIX_CACHE(j)

  if isfile(file)
    @info "Loading previously computed data..."
    m = load(file, "m") 
  else
    @info "Running new computation..."
    m = Moebius.M(indices, R, a, atol=1e-10)
    save(file, "m", m)
  end

  as[j] = a
  big_evals[j, :] = sort(svdvals(big.(m)))
  # Note: We use `svdvals` functions from the GenericSVD package
  # because for small `a` the eigenvalues are pairwise close to each
  # other and the generic Julia eigenvalue finding method gives
  # incorrect results (probably due to the 64 bit floating point
  # arithmetic).

  # Convergence of true eigenvectors to the eigenvectors
  # of the not-so-fake model
  eigenpairs = Moebius.not_so_fake_eigenpairs(R, a, 1.2*big_evals[j, SHOW_EIGENVECTORS])
  svdfact = svd(big.(m))
  shuffle_pair = false # see below

  for k in 1:SHOW_EIGENVECTORS
    v = Float64.(svdfact.U[:, N+1-k])
    e = svdfact.S[N+1-k]
    func1 = x -> Moebius.getValue(v, indices, R, x)
    evs_conv[j, k] = Moebius.compareEigenvectors(func1, eigenpairs[k][2], R) / a^2

    if evs_conv[j, k] > 0.1
      println("Suspicious result for k=$(k)!")
      #
      # Disabled debugging attempt to take care of almost
      # degenerate eigenvalues.
      #
      if shuffle_pair
        # this is the second one, you should compare this with the previous one
        println("Second one of the suspicious pair ($k)!")
        println("true: $e")
        println("nsf:  $(eigenpairs[k][1])")
        evs_conv[j, k] = Moebius.compareEigenvectors(func1, eigenpairs[k-1][2], R) / a^2
        shuffle_pair = false
      else
        # this is the first one, you should compare this with the next one
        println("First one of the suspicious pair ($k)!")
        println("true: $e")
        println("nsf:  $(eigenpairs[k][1])")
        evs_conv[j, k] = Moebius.compareEigenvectors(func1, eigenpairs[k+1][2], R) / a^2
        shuffle_pair = true
      end
    end
  end

  a -= da
end

#
# Post computation jobs...
#

@info "Data export (raw)..."

evals = convert(Array{Float64,2}, big_evals)

CSV.write(joinpath(BUILD_DIR, "limit_as.csv"), DataFrame(a = as))
CSV.write(joinpath(BUILD_DIR, "limit_evals.csv"), DataFrame(evals))
# export concat version too

@info "Data export (eigenvectors convergence)..."

CSV.write(joinpath(BUILD_DIR, "limit_evectors.csv"), DataFrame(hcat(as, evs_conv)))

@info "Data export (ratios not-so-fake / true)..."

# Prepare the data frame
df = DataFrame()
df[!, :a] = as
for i in 1:SHOW_RESULTS
  df[!, Symbol("val$i")] = zero(as)
end

for j in 1:anum
  a = as[j]
  nsf = Moebius.not_so_fake_eigenvalues(R, a, 1.2 * evals[j, SHOW_RESULTS])[1:SHOW_RESULTS]

  for i in 1:SHOW_RESULTS
    df[j, Symbol("val$i")] = nsf[i] / evals[j, i]
  end
end

CSV.write(joinpath(BUILD_DIR, "ratios.csv"), df)

@info "Data export (ratios of differences |not-so-fake - true| / a)..."

# Prepare the data frame
df = DataFrame()
df[!, :a] = as
for i in 1:SHOW_RESULTS
  df[!, Symbol("val$i")] = zero(as)
end

for j in 1:anum
  a = as[j]
  nsf = Moebius.not_so_fake_eigenvalues(R, a, 1.2 * evals[j, SHOW_RESULTS])[1:SHOW_RESULTS]

  for i in 1:SHOW_RESULTS
    df[j, Symbol("val$i")] = abs(nsf[i] - evals[j, i]) / a^2
  end
end

CSV.write(joinpath(BUILD_DIR, "diff_ratios.csv"), df)
