###
#
# Makefile
#

BUILD=build

all: run pdfs

run:
	julia --project=@. computations/limit.jl
	julia --project=@. computations/example01.jl
	julia --project=@. computations/example02.jl

SRC=$(notdir $(wildcard $(BUILD)/ex*eigen*.tex))

pdfs:
	$(foreach file, $(SRC), cd $(BUILD); lualatex --interaction=batchmode $(file))
	cp figs/fig_*.tex $(BUILD)/
	cd $(BUILD) && lualatex --interaction=batchmode fig_ratio.tex
	cd $(BUILD) && lualatex --interaction=batchmode fig_diff_ratio.tex
	cd $(BUILD) && lualatex --interaction=batchmode fig_diff_ratio_annotated.tex
	cd $(BUILD) && lualatex --interaction=batchmode fig_eigenvectors.tex

tests:
	julia --project=@. test/runtests.jl
