###
#
# Custom plotting functions based on the powerful
# PGFPlots LaTeX package.
#

"""
    data_plot(data, filename[, title="", density=true])

Custom matrix plot, produces LaTeX source file.
"""
function data_plot(data, filename; title="", density=true)

file = open(filename, "w")

write(file, """
\\documentclass{standalone}
\\usepackage[T1]{fontenc}
\\usepackage[utf8]{inputenc}
\\usepackage{pgfplots}
\\pgfplotsset{compat=1.16}
\\usepgfplotslibrary{patchplots}
\\begin{document}
\\begin{tikzpicture}[scale=0.8]
\\begin{axis}[""")

if title != ""
    write(file, "title={$(title)},")
end

write(file, """view={0}{90}, xlabel=\$s\$, ylabel=\$u\$]
\\addplot3[surf, shader=interp, patch type=bilinear]
coordinates {\n""")

for j in 1:size(data)[1]
    if density
        write(file, "($(data[j, 1]),$(data[j, 2]),$(data[j, 3]^2))\n")
    else
        write(file, "($(data[j, 1]),$(data[j, 2]),$(data[j, 3]))\n")
    end
    if data[j, 4] > 0
        write(file, "\n")
    end
end

write(file, """
};
\\end{axis}
\\end{tikzpicture}
\\end{document}
""")

close(file)
end # data_plot


"""
    data_plot_3d(data, filename[, title="", density=true])

Custom 3D visualization of eigenfunction on Möbius strip,
produces LaTeX source file.
"""
function data_plot_3d(data, R, a, filename; title="", density=true)

file = open(filename, "w")

write(file, """
\\documentclass{standalone}
\\usepackage[T1]{fontenc}
\\usepackage[utf8]{inputenc}
\\usepackage{pgfplots}
\\pgfplotsset{compat=1.16}
\\usepgfplotslibrary{patchplots}
\\begin{document}
\\begin{tikzpicture}[scale=0.8]
\\begin{axis}[""")

if title != ""
    write(file, "title={$(title)},")
end

write(file, """view={20}{50}, xlabel=\$x\$, ylabel=\$y\$, zlabel=\$z\$]
\\addplot3[surf, shader=interp, point meta=explicit]
coordinates {\n""")

# parametrization

x(s, u) = (R - a * u * cos(s/(2*R))) * cos(s/R)
y(s, u) = (R - a * u * cos(s/(2*R))) * sin(s/R)
z(s, u) = -a * u * sin(s/(2*R))

for j in 1:size(data)[1]
    if density
        write(file, "($(x(data[j, 1], data[j, 2])),$(y(data[j, 1], data[j, 2])),$(z(data[j, 1], data[j, 2]))) [$(data[j, 3]^2)]\n")
    else
        write(file, "($(x(data[j, 1], data[j, 2])),$(y(data[j, 1], data[j, 2])),$(z(data[j, 1], data[j, 2]))) [$(data[j, 3])]\n")    end
    if data[j, 4] > 0
        write(file, "\n")
    end
end

write(file, """
};
\\end{axis}
\\end{tikzpicture}
\\end{document}
""")

close(file)
end # data_plot_3d
