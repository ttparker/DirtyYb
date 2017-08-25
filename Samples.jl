const titleN = "N = 2116"
const dir = "Delta0p4/L46/"
const files = ["Real1", "Real2", "Real3", "Real4", "Real5", "Real6"]
const TsToSample = [1.0, 1.2, 1.4, 1.5, 1.53, 1.56, 1.58, 1.6, 1.7, 1.8, 2.0]
const nSteps = 500

import JLD
include("CustomTypes.jl")
using Plots
pyplot()
Plots.scalefontsizes(1.3)

for file in files
  runs = JLD.load(dir * file * ".jld", "runs")::Vector{RawRunData}
  mkpath("Samples/" * dir * file * "/")
  for T in TsToSample
    psis = runs[findfirst(x -> x.T == T, runs)].psiList
    stepSize = div(length(psis), nSteps)
    selectedPsis = psis[stepSize:stepSize:length(psis)]
    plot(Plots.partialcircle(0, 2 * pi, 100, 1.), xlims = (-1, 1), ylims = (-1, 1), label = "", seriescolor = :grey, aspect_ratio = :equal, xlabel = "Re \$\\psi\$", ylabel = "Im \$\\psi\$", title = titleN * ", T = " * string(T))
    plot!(real(selectedPsis), imag(selectedPsis), seriestype = :scatter, label = "")
    savefig("Samples/" * dir * file * "/T" * string(T) * ".png")
    println("T = " * string(T) * " analyzed")
  end
  println("Realization ", file, " analyzed")
end
