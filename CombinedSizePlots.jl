const dir = "JNeg/Delta0p4/"
const Delta = 0.4
const sizes = ["L46/","L64/","L90/","L128/"]
const filename = "Averaged.jld"

import JLD
include("CustomTypes.jl")
using Plots

function combineSizes(plotylabel::AbstractString, f::Function, plotFilename::String; plotylims = (-Inf, Inf))
  plot(xlabel = "T", ylabel = plotylabel, title = "In-plane phase, \$\\Delta = " * string(Delta) * "\$", ylims = plotylims)
  for i in 1:length(systems)
    data = f(systems[i])
    plot!(systems[i].Ts, data, seriestype = :line, color = i, label = "")
    plot!(systems[i].Ts, data, seriestype = :scatter, color = i, label = "N = " * string(systems[i].params.L^2))
  end
  savefig("Analyzed/" * dir * "Plots/" * plotFilename)
end

systems = [JLD.load("Analyzed/" * dir * size * filename, "system")::SystemSummary for size in sizes]
pyplot(markersize = 6)
Plots.scalefontsizes(1.5)
mkpath("Analyzed/" * dir * "Plots/")

combineSizes("E / N", x -> x.avgEnergies, "Energy.png")
combineSizes("C", x -> x.energyVariances .* x.params.L.^2 ./ x.Ts.^2, "HeatCapacity.png", plotylims = (-Inf, 2.0))
combineSizes("\$\\left| \\left \\langle \\psi(x) \\right \\rangle \\right|\$", x -> x.absAvgPsis, "AbsAvgPsis.png")
combineSizes("\$\\left \\langle \\left| \\psi(x) \\right| \\right \\rangle \$", x -> x.avgAbsPsis, "AvgAbsPsis.png")
combineSizes("\$\\left| \\left \\langle \\psi(x)^2 \\right \\rangle \\right|\$", x -> x.absAvgPsis, "AbsAvgPsi2s.png")
combineSizes("\$\\left \\langle \\left| \\psi^2(x) \\right| \\right \\rangle \$", x -> x.avgAbsPsis, "AvgAbsPsi2s.png")
combineSizes("\$\\left| \\left \\langle \\psi(x)^3 \\right \\rangle \\right|\$", x -> x.absAvgPsis, "AbsAvgPsi3s.png")
combineSizes("\$\\left \\langle \\left| \\psi(x)^3 \\right| \\right \\rangle \$", x -> x.avgAbsPsis, "AvgAbsPsi3s.png")
