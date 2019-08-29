###
#
# Makefile
#

BUILD=build

all: run pdfs

run:
	julia --project=@. computations/limit.jl
	julia --project=@. computations/eigenvalues_and_eigenvectors.jl

SRC=$(notdir $(wildcard $(BUILD)/eigen*.tex))

pdfs:
	$(foreach file, $(SRC), cd build; lualatex --interaction=batchmode $(file))
	cp figs/fig_*.tex $(BUILD)/
	cd build; lualatex --interaction=batchmode fig_ratio.tex
	cd build; lualatex --interaction=batchmode fig_diff_ratio.tex

tests:
	julia --project=@. tests/runtests.jl

